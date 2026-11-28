#define _POSIX_C_SOURCE 200809L
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <time.h>
#endif
#include "zstd.h"

typedef struct { unsigned char *data; size_t len; } Buffer;
typedef struct { Buffer compressed, plain; } Frame;

static void fail(const char *reason) {
    fprintf(stderr, "%s\n", reason);
    exit(1);
}
static unsigned char *allocate(size_t len) {
    unsigned char *bytes = malloc(len == 0 ? 1 : len);
    if (!bytes)
        fail("allocation failed");
    return bytes;
}
static Buffer read_file(const char *path) {
    FILE *f = fopen(path, "rb");
    if (!f || fseek(f, 0, SEEK_END))
        fail("unable to open fixture");
    long len = ftell(f);
    if (len < 0 || fseek(f, 0, SEEK_SET))
        fail("unable to size fixture");
    Buffer out = {allocate((size_t)len), (size_t)len};
    if (fread(out.data, 1, out.len, f) != out.len || fclose(f))
        fail("unable to read fixture");
    return out;
}
static size_t extract_frames(Frame *frames, size_t capacity, Buffer file) {
    static const unsigned char magic[] = {0x28, 0xb5, 0x2f, 0xfd};
    size_t n = 0;
    for (size_t off = 0; off + sizeof magic <= file.len; ++off) {
        if (memcmp(file.data + off, magic, sizeof magic))
            continue;
        size_t len = ZSTD_findFrameCompressedSize(file.data + off, file.len - off);
        if (ZSTD_isError(len))
            continue;
        unsigned long long plain = ZSTD_getFrameContentSize(file.data + off, len);
        if (n == capacity || plain > SIZE_MAX || plain == ZSTD_CONTENTSIZE_ERROR ||
            plain == ZSTD_CONTENTSIZE_UNKNOWN)
            fail("unexpected fixture frame");
        frames[n++] = (Frame){{file.data + off, len}, {allocate((size_t)plain), (size_t)plain}};
        off += len - 1;
    }
    return n;
}
static uint64_t elapsed(void) {
#ifdef _WIN32
    FILETIME created, finished, kernel, user;
    if (!GetThreadTimes(GetCurrentThread(), &created, &finished, &kernel, &user))
        fail("unable to read thread time");
    const uint64_t kernel_ticks = ((uint64_t)kernel.dwHighDateTime << 32) | kernel.dwLowDateTime;
    const uint64_t user_ticks = ((uint64_t)user.dwHighDateTime << 32) | user.dwLowDateTime;
    return (kernel_ticks + user_ticks) * 100;
#else
    struct timespec ts;
    if (clock_gettime(CLOCK_THREAD_CPUTIME_ID, &ts))
        fail("unable to read thread time");
    return (uint64_t)ts.tv_sec * UINT64_C(1000000000) + (uint64_t)ts.tv_nsec;
#endif
}
static void decode(ZSTD_DCtx *ctx, const Frame *frames, size_t count) {
    for (size_t i = 0; i < count; ++i) {
        const Frame *f = &frames[i];
        if (ZSTD_decompressDCtx(ctx, f->plain.data, f->plain.len, f->compressed.data,
                                f->compressed.len) != f->plain.len)
            fail("fixture did not decode");
    }
}
static size_t compress_input(ZSTD_CCtx *ctx, Buffer source, Buffer dict,
                             Buffer output, int min_match_four) {
    if (!min_match_four)
        return ZSTD_compress_usingDict(ctx, output.data, output.len,
                                       source.data, source.len, dict.data, dict.len, 3);
    size_t status = ZSTD_CCtx_refPrefix(ctx, dict.data, dict.len);
    if (ZSTD_isError(status))
        return status;
    return ZSTD_compress2(ctx, output.data, output.len, source.data, source.len);
}

static int compress_fixture(int argc, char **argv) {
    if (argc < 4 || argc > 5)
        fail("usage: runner --compress|--compress-min4 rounds source [raw-dict]");
    int min_match_four = !strcmp(argv[1], "--compress-min4");
    char *tail;
    unsigned long rounds = strtoul(argv[2], &tail, 10);
    if (*tail)
        fail("invalid rounds");
    Buffer source = read_file(argv[3]);
    Buffer dict = argc == 5 ? read_file(argv[4]) : (Buffer){NULL, 0};
    size_t capacity = ZSTD_compressBound(source.len);
    unsigned char *encoded = allocate(capacity);
    unsigned char *decoded = allocate(source.len);
    ZSTD_CCtx *compressor = ZSTD_createCCtx();
    ZSTD_DCtx *decoder = ZSTD_createDCtx();
    if (!compressor || !decoder)
        fail("unable to create codec contexts");
    if (min_match_four &&
        (ZSTD_isError(ZSTD_CCtx_setParameter(compressor, ZSTD_c_compressionLevel, 3)) ||
         ZSTD_isError(ZSTD_CCtx_setParameter(compressor, ZSTD_c_strategy, ZSTD_dfast)) ||
         ZSTD_isError(ZSTD_CCtx_setParameter(compressor, ZSTD_c_minMatch, 4))))
        fail("unable to configure doubleFast with minMatch=4");
    Buffer output = {encoded, capacity};
    size_t len = compress_input(compressor, source, dict, output, min_match_four);
    if (ZSTD_isError(len))
        fail(ZSTD_getErrorName(len));
    const uint64_t start = elapsed();
    for (unsigned long i = 0; i < rounds; ++i) {
        len = compress_input(compressor, source, dict, output, min_match_four);
        if (ZSTD_isError(len))
            fail(ZSTD_getErrorName(len));
    }
    const uint64_t time = elapsed() - start;
    if (ZSTD_decompress_usingDict(decoder, decoded, source.len, encoded, len,
                                  dict.data, dict.len) != source.len ||
        memcmp(source.data, decoded, source.len))
        fail("fixture roundtrip differed");
    uint64_t hash = UINT64_C(14695981039346656037);
    for (size_t i = 0; i < len; ++i)
        hash = (hash ^ encoded[i]) * UINT64_C(1099511628211);
    printf("%" PRIu64 " %016" PRIx64 "\n", time, hash);

    ZSTD_freeCCtx(compressor);
    ZSTD_freeDCtx(decoder);
    free(encoded);
    free(decoded);
    free(dict.data);
    free(source.data);
    return 0;
}

static int decode_cold_dict(int argc, char **argv) {
    if (argc != 5)
        fail("usage: runner --decode-cold-dict rounds source raw-dict");
    char *tail;
    unsigned long rounds = strtoul(argv[2], &tail, 10);
    if (*tail)
        fail("invalid rounds");
    Buffer source = read_file(argv[3]);
    Buffer dict = read_file(argv[4]);
    Buffer encoded = {NULL, ZSTD_compressBound(source.len)};
    encoded.data = allocate(encoded.len);
    unsigned char *decoded = allocate(source.len);
    ZSTD_CCtx *compressor = ZSTD_createCCtx();
    ZSTD_DCtx *decoder = ZSTD_createDCtx();
    ZSTD_DDict *dicts[2] = {ZSTD_createDDict(dict.data, dict.len),
                           ZSTD_createDDict(dict.data, dict.len)};
    if (!compressor || !decoder || !dicts[0] || !dicts[1])
        fail("unable to create codec contexts");
    size_t len = compress_input(compressor, source, dict, encoded, 0);
    if (ZSTD_isError(len))
        fail(ZSTD_getErrorName(len));
    size_t plain = ZSTD_decompress_usingDDict(decoder, decoded, source.len,
                                             encoded.data, len, dicts[1]);
    if (plain != source.len || memcmp(decoded, source.data, source.len))
        fail("fixture warmup differed");
    // Different DDict backing addresses keep the decoder's cold-dictionary
    // path active even though the two dictionaries have identical contents.
    const uint64_t start = elapsed();
    for (unsigned long i = 0; i < rounds; ++i) {
        plain = ZSTD_decompress_usingDDict(decoder, decoded, source.len,
                                           encoded.data, len, dicts[i % 2]);
        if (plain != source.len)
            fail("fixture did not decode");
    }
    const uint64_t time = elapsed() - start;
    if (memcmp(decoded, source.data, source.len))
        fail("fixture roundtrip differed");
    uint64_t hash = UINT64_C(14695981039346656037);
    for (size_t i = 0; i < source.len; ++i)
        hash = (hash ^ decoded[i]) * UINT64_C(1099511628211);
    printf("%" PRIu64 " %016" PRIx64 "\n", time, hash);
    ZSTD_freeCCtx(compressor);
    ZSTD_freeDCtx(decoder);
    ZSTD_freeDDict(dicts[0]);
    ZSTD_freeDDict(dicts[1]);
    free(encoded.data);
    free(decoded);
    free(dict.data);
    free(source.data);
    return 0;
}

int main(int argc, char **argv) {
    if (argc >= 2 && !strcmp(argv[1], "--decode-cold-dict"))
        return decode_cold_dict(argc, argv);
    if (argc >= 2 && (!strcmp(argv[1], "--compress") ||
                      !strcmp(argv[1], "--compress-min4")))
        return compress_fixture(argc, argv);
    if (argc < 4 || argc > 5)
        fail("usage: runner rounds parquet hdf5 [bulk-file]");
    char *tail;
    unsigned long rounds = strtoul(argv[1], &tail, 10);
    if (*tail)
        fail("invalid rounds");
    Frame frames[16];
    Buffer parquet = read_file(argv[2]), hdf5 = read_file(argv[3]);
    size_t count = extract_frames(frames, 16, parquet);
    count += extract_frames(frames + count, 16 - count, hdf5);
    if (count != 8)
        fail("fixture frame count changed");
    Buffer bulk = {NULL, 0};
    if (argc == 5) {
        bulk = read_file(argv[4]);
        size_t capacity = ZSTD_compressBound(bulk.len);
        unsigned char *temp = allocate(capacity);
        size_t len = ZSTD_compress(temp, capacity, bulk.data, bulk.len, 3);
        if (ZSTD_isError(len))
            fail("unable to encode bulk corpus");
        frames[count++] = (Frame){{temp, len}, {allocate(bulk.len), bulk.len}};
    }
    ZSTD_DCtx *ctx = ZSTD_createDCtx();
    if (!ctx)
        fail("unable to create decompression context");
    decode(ctx, frames, count);
    const uint64_t start = elapsed();
    for (unsigned long i = 0; i < rounds; ++i)
        decode(ctx, frames, count);
    const uint64_t time = elapsed() - start;
    uint64_t hash = UINT64_C(14695981039346656037);
    for (size_t i = 0; i < count; ++i)
        for (size_t j = 0; j < frames[i].plain.len; ++j)
            hash = (hash ^ frames[i].plain.data[j]) * UINT64_C(1099511628211);
    if (bulk.data) {
        if (memcmp(bulk.data, frames[count - 1].plain.data, bulk.len))
            fail("bulk roundtrip differed");
        free(frames[count - 1].compressed.data);
        free(bulk.data);
    }
    for (size_t i = 0; i < count; ++i)
        free(frames[i].plain.data);
    free(parquet.data);
    free(hdf5.data);
    ZSTD_freeDCtx(ctx);
    printf("%" PRIu64 " %016" PRIx64 "\n", time, hash);
}
