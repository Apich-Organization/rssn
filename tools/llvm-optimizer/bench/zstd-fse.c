#define _POSIX_C_SOURCE 200809L
#define FSE_STATIC_LINKING_ONLY
#include "fse.h"
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* Replay the compressed Huffman weights logged by analysis/zstd_entropy.gdb. */
typedef struct {
    unsigned char encoded[256], decoded[256];
    size_t compressed_size, decoded_size;
    unsigned capacity;
} Header;

static uint64_t now(void) {
    struct timespec ts;
    if (clock_gettime(CLOCK_THREAD_CPUTIME_ID, &ts))
        abort();
    return (uint64_t)ts.tv_sec * 1000000000 + ts.tv_nsec;
}

static size_t replay(const Header *header, unsigned char *output) {
    unsigned work[FSE_DECOMPRESS_WKSP_SIZE_U32(6, 12)];
    return FSE_decompress_wksp_bmi2(output, header->capacity, header->encoded,
                                   header->compressed_size, 6, work, sizeof work, 1);
}

int main(int argc, char **argv) {
    if (argc != 3) {
        fputs("usage: replay headers.txt rounds\n", stderr);
        return 1;
    }
    FILE *file = fopen(argv[1], "r");
    if (!file)
        abort();
    Header headers[128];
    unsigned count = 0, capacity, log;
    char hex[513];
    int fields;
    while ((fields = fscanf(file, "FSE %u %u %512s\n", &capacity, &log, hex)) == 3) {
        size_t len = strlen(hex);
        if (!len || len % 2 || capacity > sizeof headers[0].decoded ||
            log != 6 || count == sizeof headers / sizeof headers[0])
            abort();
        Header *header = &headers[count++];
        header->capacity = capacity;
        header->compressed_size = len / 2;
        for (size_t i = 0; i < header->compressed_size; ++i) {
            unsigned byte;
            if (sscanf(hex + 2 * i, "%2x", &byte) != 1)
                abort();
            header->encoded[i] = byte;
        }
        header->decoded_size = replay(header, header->decoded);
        if (header->decoded_size > capacity)
            abort();
    }
    if (fields != EOF || ferror(file) || fclose(file) || !count)
        abort();
    char *tail;
    unsigned long rounds = strtoul(argv[2], &tail, 10);
    if (*tail || !rounds)
        abort();

    unsigned char output[256];
    const uint64_t start = now();
    for (unsigned long round = 0; round < rounds; ++round)
        for (unsigned i = 0; i < count; ++i)
            if (replay(&headers[i], output) != headers[i].decoded_size)
                abort();
    const uint64_t elapsed = now() - start;

    uint64_t hash = UINT64_C(14695981039346656037);
    for (unsigned i = 0; i < count; ++i)
        for (size_t j = 0; j < headers[i].decoded_size; ++j)
            hash = (hash ^ headers[i].decoded[j]) * UINT64_C(1099511628211);
    printf("%" PRIu64 " %016" PRIx64 " %u\n", elapsed, hash, count);
}
