/*
 * rssn_r.c -- the C bridge between R and the rssn C ABI (`rssn.h`).
 *
 * R can only call a C function through `.Call`, and `.Call` speaks SEXP, while
 * the rssn ABI speaks `void *` handles, out-parameters and C strings. This file
 * is that translation, and nothing else: it holds no policy about what a
 * session or a term means. `R/session.R`, `R/term.R` and friends build the
 * ergonomic layer on top.
 *
 * Three properties are worth knowing before reading on.
 *
 * 1. The library is loaded at run time, not linked. `rssn_load()` dlopen()s
 *    `rssn.dll` / `librssn.so` / `librssn.dylib` from `RSSN_LIB`, from the
 *    package's `libs/`, or from a source checkout's build directories, and
 *    resolves every entry point once into `api`. The package therefore has no
 *    link-time dependency on rssn: any build of the shared library works with
 *    any build of this package, exactly as the Python and Julia wrappers work.
 *
 * 2. Failures become R conditions, not C errors. A status code is turned into
 *    a condition of class `c("rssn_error", "error", "condition")` carrying
 *    `status` (the integer code) and `status_name`, so `tryCatch(expr,
 *    rssn_error = function(e) e$status)` works and the message is the one rssn
 *    reported.
 *
 * 3. A host callback is an R closure. The engine calls a C trampoline, which
 *    looks the closure up in a registry, hands it the argument cell ids, and
 *    takes back a cell id. The R-level conversion (cells -> R objects and back)
 *    lives in `R/callback.R` where it is easy to read and to test; the
 *    trampoline only marshals integers. R errors are caught by `R_tryEval` so
 *    they never unwind through the engine's stack frames.
 */

#define R_NO_REMAP

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#ifdef _WIN32
#  define WIN32_LEAN_AND_MEAN
#  include <windows.h>
typedef HMODULE rssn_lib_handle;
static rssn_lib_handle lib_open(const char *path) { return LoadLibraryA(path); }
static void *lib_sym(rssn_lib_handle handle, const char *name) {
    return (void *)GetProcAddress(handle, name);
}
static void lib_close(rssn_lib_handle handle) {
    if (handle != NULL) FreeLibrary(handle);
}
#else
#  include <dlfcn.h>
typedef void *rssn_lib_handle;
static rssn_lib_handle lib_open(const char *path) { return dlopen(path, RTLD_LAZY | RTLD_GLOBAL); }
static void *lib_sym(rssn_lib_handle handle, const char *name) { return dlsym(handle, name); }
static void lib_close(rssn_lib_handle handle) {
    if (handle != NULL) dlclose(handle);
}
#endif

/* ------------------------------------------------------------------------- */
/* the ABI, as function pointers                                              */
/* ------------------------------------------------------------------------- */

#define RSSN_OK 0
#define RSSN_NO_TERM 0xFFFFFFFFu

typedef struct RssnAnswer {
    uint32_t term;
    bool has_value;
    double value;
    double error;
    bool reduced;
    size_t iterations;
} RssnAnswer;

typedef struct RssnAnswerMetadata {
    char *dimension;
    char *diagnostic_codes;
    char *message;
} RssnAnswerMetadata;

typedef struct RssnTensor {
    size_t rank;
    size_t *shape;
    size_t len;
    double *data;
} RssnTensor;

typedef int (*rssn_callback_fn)(void *user, const void *session, const uint32_t *args,
                                size_t count, uint32_t *out);

typedef struct {
    int loaded;
    char path[4096];

    const char *(*version)(void);
    const char *(*last_error)(void);
    char *(*manifest_json)(void);
    int (*op_info_json)(const char *, char **);
    uint64_t (*manifest_hash)(void);
    void (*string_free)(char *);
    void (*bytes_free)(uint8_t *, size_t);

    void *(*session_new)(void);
    void (*session_free)(void *);

    void *(*config_new)(void);
    void (*config_free)(void *);
    int (*config_symbolic)(void *);
    int (*config_numeric)(void *, double);
    int (*config_bind)(void *, const char *, double);
    int (*config_assume)(void *, const char *, uint8_t);
    int (*config_budget)(void *, size_t, size_t);
    int (*config_dimension_mode)(void *, uint8_t);
    int (*config_bypass_dimensions)(void *);
    int (*config_dimension)(void *, const char *, const char *);
    int (*config_define_base_dimension)(void *, const char *);

    uint32_t (*term_parse)(void *, const char *);
    uint32_t (*term_sym)(void *, const char *);
    uint32_t (*term_int)(void *, int64_t);
    uint32_t (*term_float)(void *, double);
    uint32_t (*term_rational)(void *, int64_t, int64_t);
    uint32_t (*term_apply)(void *, const char *, const uint32_t *, size_t);
    int (*term_eval)(void *, uint32_t, const char *const *, const double *, size_t, double *);
    int (*term_to_f64)(void *, uint32_t, double *);
    int (*term_tensor_data)(void *, uint32_t, RssnTensor *);
    void (*tensor_free)(RssnTensor *);
    char *(*term_to_string)(void *, uint32_t);
    char *(*term_to_latex)(void *, uint32_t);
    char *(*term_to_typst)(void *, uint32_t);
    char *(*term_to_pretty)(void *, uint32_t);

    int (*compute)(void *, uint32_t, void *, RssnAnswer *);
    int (*compute_all)(void *, const uint32_t *, size_t, void *, RssnAnswer *, int *);
    int (*compute_all_with_metadata)(void *, const uint32_t *, size_t, void *,
                                     RssnAnswer *, int *, RssnAnswerMetadata *);
    int (*compute_json)(void *, uint32_t, void *, char **);
    int (*compute_bincode)(void *, uint32_t, void *, uint8_t **, size_t *);
    int (*compute_all_bincode)(void *, const uint32_t *, size_t, void *, uint8_t **, size_t *);
    uint32_t (*simplify)(void *, const char *);
    int (*dimension)(void *, uint32_t, void *, char **);
    size_t (*dimension_count)(void);
    const char *(*dimension_code)(size_t);

    int (*value_with_coherent_dimension)(void *, uint32_t, const uint8_t *, size_t, uint32_t *);
    int (*value_new_int)(void *, int64_t, uint32_t *);
    int (*value_new_int_decimal)(void *, const char *, uint32_t *);
    int (*value_new_float)(void *, double, uint32_t *);
    int (*value_new_text)(void *, const char *, uint32_t *);
    int (*value_new_symbol)(void *, const char *, uint32_t *);
    int (*value_new_term)(void *, uint32_t, uint32_t *);
    int (*value_new_list)(void *, const uint32_t *, size_t, uint32_t *);
    int (*value_new_tensor)(void *, const double *, size_t, const size_t *, size_t, uint32_t *);
    int (*value_new_tensor_cells)(void *, const uint32_t *, size_t, const size_t *, size_t,
                                  uint32_t, uint32_t *);
    int (*value_new_bool)(void *, int, uint32_t *);
    int (*value_new_complex)(void *, double, double, uint32_t *);
    int (*value_new_record)(void *, const char *const *, const uint32_t *, size_t, uint32_t *);
    int (*value_call)(void *, const char *, const uint32_t *, size_t, uint32_t *);
    int (*value_call_key)(void *, const char *, const char *, const uint8_t *,
                          const uint32_t *, size_t, uint32_t *);
    int (*value_call_json)(void *, const char *, const uint32_t *, size_t, char **);
    int (*value_kind)(void *, uint32_t, uint8_t *);
    int (*value_to_f64)(void *, uint32_t, double *);
    int (*value_to_text)(void *, uint32_t, char **);
    int (*value_to_term)(void *, uint32_t, uint32_t *);
    int (*value_to_tensor)(void *, uint32_t, RssnTensor *);
    int (*value_len)(void *, uint32_t, size_t *);
    int (*value_at)(void *, uint32_t, size_t, uint32_t *);
    int (*value_field_name)(void *, uint32_t, size_t, char **);
    int (*value_tag)(void *, uint32_t, char **);
    int (*value_json)(void *, uint32_t, char **);
    int (*value_clone)(void *, uint32_t, uint32_t *);
    void (*value_free)(void *, uint32_t);
    size_t (*value_count)(void *);

    int (*callback_new)(void *, const char *, rssn_callback_fn, void *, uint32_t *);
    int (*callback_new_for)(void *, const char *, size_t, rssn_callback_fn, void *, uint32_t *);
    int (*callback_signatures)(void *, const char *, size_t, char **);
    void (*callback_free)(void *, uint32_t);
    size_t (*callback_count)(void *);
    void (*callback_fail)(const char *);

    int (*sim_run)(const char *, const char *, char **);
} rssn_api;

static rssn_api api;
static rssn_lib_handle api_handle = NULL;

/* ------------------------------------------------------------------------- */
/* status names and conditions                                               */
/* ------------------------------------------------------------------------- */

static const char *const STATUS_NAMES[] = {
    "ok", "null_argument", "invalid_utf8", "parse", "unknown_operator", "not_numeric",
    "inconsistent", "rule", "no_such_term", "no_value", "simulation", "panic",
    "dimension", "invalid_argument", "cancelled",
};
#define STATUS_COUNT ((int)(sizeof(STATUS_NAMES) / sizeof(STATUS_NAMES[0])))

static const char *status_name(int code) {
    return (code >= 0 && code < STATUS_COUNT) ? STATUS_NAMES[code] : "unknown";
}

/* The condition `stop()` is given: a list with `message`, `call`, `status` and
 * `status_name`, classed `rssn_error`. The status fields are elements rather
 * than attributes so that `e$status` and `e$status_name` read naturally. */
static void raise_status(int code, const char *fallback) {
    const char *detail = api.loaded && api.last_error != NULL ? api.last_error() : "";
    if (detail == NULL || detail[0] == '\0') detail = fallback;

    SEXP condition = PROTECT(Rf_allocVector(VECSXP, 4));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 4));
    SET_STRING_ELT(names, 0, Rf_mkChar("message"));
    SET_STRING_ELT(names, 1, Rf_mkChar("call"));
    SET_STRING_ELT(names, 2, Rf_mkChar("status"));
    SET_STRING_ELT(names, 3, Rf_mkChar("status_name"));
    SET_VECTOR_ELT(condition, 0, Rf_mkString(detail));
    SET_VECTOR_ELT(condition, 1, R_NilValue);
    SET_VECTOR_ELT(condition, 2, Rf_ScalarInteger(code));
    SET_VECTOR_ELT(condition, 3, Rf_mkString(status_name(code)));
    Rf_setAttrib(condition, R_NamesSymbol, names);

    SEXP classes = PROTECT(Rf_allocVector(STRSXP, 3));
    SET_STRING_ELT(classes, 0, Rf_mkChar("rssn_error"));
    SET_STRING_ELT(classes, 1, Rf_mkChar("error"));
    SET_STRING_ELT(classes, 2, Rf_mkChar("condition"));
    Rf_setAttrib(condition, R_ClassSymbol, classes);

    SEXP call = PROTECT(Rf_lang2(Rf_install("stop"), condition));
    Rf_eval(call, R_BaseEnv);  /* does not return */
    UNPROTECT(4);
}

static void check_status(int code) {
    if (code != RSSN_OK) raise_status(code, "rssn call failed");
}

/* ------------------------------------------------------------------------- */
/* small helpers                                                             */
/* ------------------------------------------------------------------------- */

static const char *as_cstring(SEXP x) {
    if (x == NA_STRING) Rf_error("a string argument is NA");
    return Rf_translateCharUTF8(Rf_asChar(x));
}

static void *as_handle(SEXP x, const char *what) {
    void *pointer = R_ExternalPtrAddr(x);
    if (pointer == NULL) Rf_error("this %s has already been released", what);
    return pointer;
}

static void *optional_handle(SEXP x) { return x == R_NilValue ? NULL : as_handle(x, "object"); }

/* A caller-owned `char *` is copied into an R string and then released. */
static SEXP take_string(char *text) {
    if (text == NULL) return R_NilValue;
    SEXP out = PROTECT(Rf_mkString(text));
    api.string_free(text);
    UNPROTECT(1);
    return out;
}

static SEXP take_bytes(uint8_t *data, size_t length) {
    if (data == NULL) return Rf_allocVector(RAWSXP, 0);
    SEXP out = PROTECT(Rf_allocVector(RAWSXP, (R_xlen_t)length));
    if (length > 0) memcpy(RAW(out), data, length);
    api.bytes_free(data, length);
    UNPROTECT(1);
    return out;
}

/* A `size_t` vector as an R integer vector (the extents of a tensor). */
static SEXP shape_vector(const RssnTensor *tensor) {
    SEXP out = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t)tensor->rank));
    for (size_t index = 0; index < tensor->rank; index++) {
        INTEGER(out)[index] = (int)tensor->shape[index];
    }
    UNPROTECT(1);
    return out;
}

/* Copies a tensor out and frees the engine's buffers: `list(data, dim)`. */
static SEXP take_tensor(RssnTensor *tensor) {
    SEXP data = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)tensor->len));
    if (tensor->len > 0) memcpy(REAL(data), tensor->data, tensor->len * sizeof(double));
    SEXP dim = PROTECT(shape_vector(tensor));
    api.tensor_free(tensor);

    SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(out, 0, data);
    SET_VECTOR_ELT(out, 1, dim);
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 2));
    SET_STRING_ELT(names, 0, Rf_mkChar("data"));
    SET_STRING_ELT(names, 1, Rf_mkChar("dim"));
    Rf_setAttrib(out, R_NamesSymbol, names);
    UNPROTECT(4);
    return out;
}

/* A `uint32_t *` view of an R integer vector of handles. */
static uint32_t *handle_buffer(SEXP x, size_t *count) {
    if (!Rf_isInteger(x)) Rf_error("expected an integer vector of handles");
    *count = (size_t)Rf_xlength(x);
    return (uint32_t *)INTEGER(x);
}

/* ------------------------------------------------------------------------- */
/* loading                                                                   */
/* ------------------------------------------------------------------------- */

static const char *LIB_NAMES[] = {
#ifdef _WIN32
    "rssn.dll", "librssn.dll",
#elif defined(__APPLE__)
    "librssn.dylib", "rssn.dylib",
#else
    "librssn.so", "librssn.so.0", "rssn.so",
#endif
};
#define LIB_NAME_COUNT ((size_t)(sizeof(LIB_NAMES) / sizeof(LIB_NAMES[0])))

static int file_exists(const char *path) {
#ifdef _WIN32
    DWORD attributes = GetFileAttributesA(path);
    return attributes != INVALID_FILE_ATTRIBUTES && !(attributes & FILE_ATTRIBUTE_DIRECTORY);
#else
    struct stat metadata;
    return stat(path, &metadata) == 0 && S_ISREG(metadata.st_mode);
#endif
}

static void join_path(char *out, size_t size, const char *dir, const char *name) {
    size_t length = strlen(dir);
    int needs_separator = length > 0 && dir[length - 1] != '/' && dir[length - 1] != '\\';
    snprintf(out, size, "%s%s%s", dir, needs_separator ? "/" : "", name);
}

/*
 * One candidate: a file, or a directory that may hold the platform library
 * name. Writes the resolved path to `out` and answers whether it is a file.
 */
static int candidate_exists(const char *candidate, char *out, size_t size) {
    if (candidate[0] == '\0') return 0;
    if (file_exists(candidate)) {
        // dlopen searches its loader path for bare names. This branch has
        // already found a local file, so resolve that file explicitly.
        const char *prefix = strpbrk(candidate, "/\\:") == NULL ? "./" : "";
        int length = snprintf(out, size, "%s%s", prefix, candidate);
        return length >= 0 && (size_t)length < size;
    }
    for (size_t index = 0; index < LIB_NAME_COUNT; index++) {
        join_path(out, size, candidate, LIB_NAMES[index]);
        if (file_exists(out)) return 1;
    }
    return 0;
}

#define RESOLVE(field, symbol)                                                      \
    do {                                                                            \
        next.field = (void *)lib_sym(handle, symbol);                               \
        if (next.field == NULL) {                                                    \
            lib_close(handle);                                                      \
            Rf_error("the rssn library at '%s' has no symbol %s", path, symbol);    \
        }                                                                           \
    } while (0)

/*
 * `candidates` is the search order `R/zzz.R` built: `RSSN_LIB`, the package's
 * own `libs/`, the build directories of a source checkout, and finally the bare
 * library name, which lets the platform loader search its own path.
 */
SEXP C_rssn_load(SEXP candidates) {
    char path[4096];
    rssn_lib_handle handle = NULL;
    for (R_xlen_t index = 0; index < Rf_xlength(candidates); index++) {
        const char *candidate = Rf_translateCharUTF8(STRING_ELT(candidates, index));
        if (candidate_exists(candidate, path, sizeof(path))) {
            handle = lib_open(path);
            if (handle == NULL) Rf_error("cannot load the rssn shared library at '%s'", path);
            break;
        }
        if (candidate[0] != '\0' && strpbrk(candidate, "/\\:") == NULL) {
            int length = snprintf(path, sizeof(path), "%s", candidate);
            if (length >= 0 && (size_t)length < sizeof(path)) {
                handle = lib_open(path);
                if (handle != NULL) break;
            }
        }
    }
    if (handle == NULL) {
        Rf_error("no rssn shared library found; set RSSN_LIB to the library file or its directory, "
                 "or build one with `bash scripts/build_ffi.sh`");
    }

    rssn_api next = {0};
    RESOLVE(version, "rssn_version");
    RESOLVE(last_error, "rssn_last_error");
    RESOLVE(manifest_json, "rssn_manifest_json");
    next.op_info_json = (void *)lib_sym(handle, "rssn_op_info_json");
    RESOLVE(manifest_hash, "rssn_manifest_hash");
    RESOLVE(string_free, "rssn_string_free");
    RESOLVE(bytes_free, "rssn_bytes_free");

    RESOLVE(session_new, "rssn_session_new");
    RESOLVE(session_free, "rssn_session_free");

    RESOLVE(config_new, "rssn_config_new");
    RESOLVE(config_free, "rssn_config_free");
    RESOLVE(config_symbolic, "rssn_config_symbolic");
    RESOLVE(config_numeric, "rssn_config_numeric");
    RESOLVE(config_bind, "rssn_config_bind");
    RESOLVE(config_assume, "rssn_config_assume");
    RESOLVE(config_budget, "rssn_config_budget");
    RESOLVE(config_dimension_mode, "rssn_config_dimension_mode");
    RESOLVE(config_bypass_dimensions, "rssn_config_bypass_dimensions");
    RESOLVE(config_dimension, "rssn_config_dimension");
    RESOLVE(config_define_base_dimension, "rssn_config_define_base_dimension");

    RESOLVE(term_parse, "rssn_term_parse");
    RESOLVE(term_sym, "rssn_term_sym");
    RESOLVE(term_int, "rssn_term_int");
    RESOLVE(term_float, "rssn_term_float");
    RESOLVE(term_rational, "rssn_term_rational");
    RESOLVE(term_apply, "rssn_term_apply");
    RESOLVE(term_eval, "rssn_term_eval");
    RESOLVE(term_to_f64, "rssn_term_to_f64");
    RESOLVE(term_tensor_data, "rssn_term_tensor_data");
    RESOLVE(tensor_free, "rssn_tensor_free");
    RESOLVE(term_to_string, "rssn_term_to_string");
    RESOLVE(term_to_latex, "rssn_term_to_latex");
    RESOLVE(term_to_typst, "rssn_term_to_typst");
    RESOLVE(term_to_pretty, "rssn_term_to_pretty");

    RESOLVE(compute, "rssn_compute");
    RESOLVE(compute_all, "rssn_compute_all");
    next.compute_all_with_metadata = (void *)lib_sym(handle, "rssn_compute_all_with_metadata");
    RESOLVE(compute_json, "rssn_compute_json");
    RESOLVE(compute_bincode, "rssn_compute_bincode");
    RESOLVE(compute_all_bincode, "rssn_compute_all_bincode");
    RESOLVE(simplify, "rssn_simplify");
    RESOLVE(dimension, "rssn_dimension");
    RESOLVE(dimension_count, "rssn_dimension_count");
    RESOLVE(dimension_code, "rssn_dimension_code");

    next.value_with_coherent_dimension = (void *)lib_sym(handle, "rssn_value_with_coherent_dimension");
    RESOLVE(value_new_int, "rssn_value_new_int");
    next.value_new_int_decimal = (void *)lib_sym(handle, "rssn_value_new_int_decimal");
    RESOLVE(value_new_float, "rssn_value_new_float");
    RESOLVE(value_new_text, "rssn_value_new_text");
    RESOLVE(value_new_symbol, "rssn_value_new_symbol");
    RESOLVE(value_new_term, "rssn_value_new_term");
    RESOLVE(value_new_list, "rssn_value_new_list");
    RESOLVE(value_new_tensor, "rssn_value_new_tensor");
    next.value_new_tensor_cells = (void *)lib_sym(handle, "rssn_value_new_tensor_cells");
    RESOLVE(value_new_bool, "rssn_value_new_bool");
    RESOLVE(value_new_complex, "rssn_value_new_complex");
    RESOLVE(value_new_record, "rssn_value_new_record");
    RESOLVE(value_call, "rssn_value_call");
    next.value_call_key = (void *)lib_sym(handle, "rssn_value_call_key");
    RESOLVE(value_call_json, "rssn_value_call_json");
    RESOLVE(value_kind, "rssn_value_kind");
    RESOLVE(value_to_f64, "rssn_value_to_f64");
    RESOLVE(value_to_text, "rssn_value_to_text");
    RESOLVE(value_to_term, "rssn_value_to_term");
    RESOLVE(value_to_tensor, "rssn_value_to_tensor");
    RESOLVE(value_len, "rssn_value_len");
    RESOLVE(value_at, "rssn_value_at");
    RESOLVE(value_field_name, "rssn_value_field_name");
    RESOLVE(value_tag, "rssn_value_tag");
    RESOLVE(value_json, "rssn_value_json");
    next.value_clone = (void *)lib_sym(handle, "rssn_value_clone");
    RESOLVE(value_free, "rssn_value_free");
    RESOLVE(value_count, "rssn_value_count");

    RESOLVE(callback_new, "rssn_callback_new");
    RESOLVE(callback_new_for, "rssn_callback_new_for");
    RESOLVE(callback_signatures, "rssn_callback_signatures");
    RESOLVE(callback_free, "rssn_callback_free");
    RESOLVE(callback_count, "rssn_callback_count");
    RESOLVE(callback_fail, "rssn_callback_fail");

    RESOLVE(sim_run, "rssn_sim_run");

    if (api_handle != NULL) lib_close(api_handle);
    api_handle = handle;
    next.loaded = 1;
    snprintf(next.path, sizeof(next.path), "%s", path);
    api = next;
    return Rf_mkString(api.path);
}

SEXP C_rssn_libpath(void) {
    if (!api.loaded) return Rf_mkString("");
    return Rf_mkString(api.path);
}

SEXP C_rssn_version(void) {
    if (!api.loaded) Rf_error("the rssn library is not loaded");
    const char *text = api.version();
    return Rf_mkString(text == NULL ? "" : text);
}

SEXP C_rssn_last_error(void) {
    if (!api.loaded) return Rf_mkString("");
    const char *text = api.last_error();
    return Rf_mkString(text == NULL ? "" : text);
}

/* The hash is a `uint64_t`; R's double cannot hold it exactly, so it is text. */
SEXP C_rssn_manifest_hash(void) {
    if (!api.loaded) Rf_error("the rssn library is not loaded");
    char buffer[32];
    snprintf(buffer, sizeof(buffer), "%llu", (unsigned long long)api.manifest_hash());
    return Rf_mkString(buffer);
}

SEXP C_rssn_manifest_json(void) {
    if (!api.loaded) Rf_error("the rssn library is not loaded");
    return take_string(api.manifest_json());
}

SEXP C_rssn_op_info_json(SEXP name) {
    if (api.op_info_json == NULL) Rf_error("this RSSN library lacks operator contract lookup");
    char *json = NULL;
    check_status(api.op_info_json(as_cstring(name), &json));
    return take_string(json);
}

/* ------------------------------------------------------------------------- */
/* sessions and configurations                                               */
/* ------------------------------------------------------------------------- */

static void session_finalizer(SEXP pointer) {
    void *session = R_ExternalPtrAddr(pointer);
    if (session != NULL) {
        api.session_free(session);
        R_ClearExternalPtr(pointer);
    }
}

static void config_finalizer(SEXP pointer) {
    void *config = R_ExternalPtrAddr(pointer);
    if (config != NULL) {
        api.config_free(config);
        R_ClearExternalPtr(pointer);
    }
}

SEXP C_rssn_session_new(void) {
    if (!api.loaded) Rf_error("the rssn library is not loaded");
    void *session = api.session_new();
    if (session == NULL) raise_status(7, "rssn_session_new failed");
    SEXP pointer = PROTECT(R_MakeExternalPtr(session, Rf_install("rssn_session"), R_NilValue));
    R_RegisterCFinalizerEx(pointer, session_finalizer, TRUE);
    UNPROTECT(1);
    return pointer;
}

SEXP C_rssn_session_free(SEXP pointer) {
    session_finalizer(pointer);
    return R_NilValue;
}

SEXP C_rssn_config_new(void) {
    if (!api.loaded) Rf_error("the rssn library is not loaded");
    void *config = api.config_new();
    if (config == NULL) raise_status(7, "rssn_config_new failed");
    SEXP pointer = PROTECT(R_MakeExternalPtr(config, Rf_install("rssn_config"), R_NilValue));
    R_RegisterCFinalizerEx(pointer, config_finalizer, TRUE);
    UNPROTECT(1);
    return pointer;
}

SEXP C_rssn_config_free(SEXP pointer) {
    config_finalizer(pointer);
    return R_NilValue;
}

#define CONFIG_SETTER(name, call)                                                   \
    SEXP name(SEXP config) {                                                        \
        check_status(call);                                                         \
        return config;                                                              \
    }

CONFIG_SETTER(C_rssn_config_symbolic, api.config_symbolic(as_handle(config, "config")))
CONFIG_SETTER(C_rssn_config_bypass_dimensions,
              api.config_bypass_dimensions(as_handle(config, "config")))

SEXP C_rssn_config_numeric(SEXP config, SEXP tolerance) {
    check_status(api.config_numeric(as_handle(config, "config"), Rf_asReal(tolerance)));
    return config;
}

SEXP C_rssn_config_bind(SEXP config, SEXP name, SEXP value) {
    check_status(api.config_bind(as_handle(config, "config"), as_cstring(name),
                                 Rf_asReal(value)));
    return config;
}

SEXP C_rssn_config_assume(SEXP config, SEXP name, SEXP facts) {
    check_status(api.config_assume(as_handle(config, "config"), as_cstring(name),
                                   (uint8_t)Rf_asInteger(facts)));
    return config;
}

SEXP C_rssn_config_budget(SEXP config, SEXP iterations, SEXP nodes) {
    check_status(api.config_budget(as_handle(config, "config"),
                                   (size_t)Rf_asReal(iterations), (size_t)Rf_asReal(nodes)));
    return config;
}

SEXP C_rssn_config_dimension_mode(SEXP config, SEXP mode) {
    check_status(api.config_dimension_mode(as_handle(config, "config"),
                                           (uint8_t)Rf_asInteger(mode)));
    return config;
}

SEXP C_rssn_config_dimension(SEXP config, SEXP symbol, SEXP text) {
    check_status(api.config_dimension(as_handle(config, "config"), as_cstring(symbol),
                                      as_cstring(text)));
    return config;
}

SEXP C_rssn_config_define_base_dimension(SEXP config, SEXP name) {
    check_status(api.config_define_base_dimension(as_handle(config, "config"),
                                                  as_cstring(name)));
    return config;
}

/* ------------------------------------------------------------------------- */
/* terms                                                                     */
/* ------------------------------------------------------------------------- */

static SEXP term_handle(uint32_t handle, int hint) {
    if (handle == RSSN_NO_TERM) raise_status(hint, "the term could not be built");
    return Rf_ScalarInteger((int)handle);
}

SEXP C_rssn_term_parse(SEXP session, SEXP source) {
    return term_handle(api.term_parse(as_handle(session, "session"), as_cstring(source)), 3);
}

SEXP C_rssn_term_sym(SEXP session, SEXP name) {
    return term_handle(api.term_sym(as_handle(session, "session"), as_cstring(name)), 4);
}

static int64_t as_i64(SEXP value) {
    double number = Rf_asReal(value);
    const double signed_limit = 9223372036854775808.0;
    if (!R_FINITE(number) || number < -signed_limit || number >= signed_limit ||
        number != trunc(number)) {
        Rf_error("expected a whole number in the signed 64-bit range; use value_integer_decimal for larger values");
    }
    return (int64_t)number;
}

SEXP C_rssn_term_int(SEXP session, SEXP value) {
    return term_handle(api.term_int(as_handle(session, "session"), as_i64(value)), 13);
}

SEXP C_rssn_term_float(SEXP session, SEXP value) {
    return term_handle(api.term_float(as_handle(session, "session"), Rf_asReal(value)), 13);
}

SEXP C_rssn_term_rational(SEXP session, SEXP numerator, SEXP denominator) {
    return term_handle(api.term_rational(as_handle(session, "session"),
                                         as_i64(numerator),
                                         as_i64(denominator)),
                       13);
}

SEXP C_rssn_term_apply(SEXP session, SEXP op, SEXP handles) {
    size_t count = 0;
    uint32_t *buffer = handle_buffer(handles, &count);
    return term_handle(api.term_apply(as_handle(session, "session"), as_cstring(op), buffer,
                                      count),
                       4);
}

SEXP C_rssn_term_eval(SEXP session, SEXP term, SEXP names, SEXP values) {
    size_t count = (size_t)Rf_xlength(names);
    const char **name_buffer = (const char **)R_alloc(count > 0 ? count : 1, sizeof(char *));
    for (size_t index = 0; index < count; index++) {
        name_buffer[index] = Rf_translateCharUTF8(STRING_ELT(names, (R_xlen_t)index));
    }
    double out = 0.0;
    check_status(api.term_eval(as_handle(session, "session"), (uint32_t)Rf_asInteger(term),
                               name_buffer, REAL(values), count, &out));
    return Rf_ScalarReal(out);
}

/* The value of a literal, or `NA` when the term is symbolic. */
SEXP C_rssn_term_to_f64(SEXP session, SEXP term) {
    double out = 0.0;
    int code = api.term_to_f64(as_handle(session, "session"), (uint32_t)Rf_asInteger(term), &out);
    if (code == 9) return Rf_ScalarReal(NA_REAL);  /* no_value */
    check_status(code);
    return Rf_ScalarReal(out);
}

SEXP C_rssn_term_tensor(SEXP session, SEXP term) {
    RssnTensor tensor;
    memset(&tensor, 0, sizeof(tensor));
    check_status(api.term_tensor_data(as_handle(session, "session"),
                                      (uint32_t)Rf_asInteger(term), &tensor));
    return take_tensor(&tensor);
}

#define RENDERER(name, field)                                                       \
    SEXP name(SEXP session, SEXP term) {                                            \
        return take_string(api.field(as_handle(session, "session"),                 \
                                     (uint32_t)Rf_asInteger(term)));                \
    }

RENDERER(C_rssn_term_to_string, term_to_string)
RENDERER(C_rssn_term_to_latex, term_to_latex)
RENDERER(C_rssn_term_to_typst, term_to_typst)
RENDERER(C_rssn_term_to_pretty, term_to_pretty)

/* ------------------------------------------------------------------------- */
/* computation                                                               */
/* ------------------------------------------------------------------------- */

/* Columnar, so `compute_all` does not build one list per answer. */
SEXP C_rssn_compute(SEXP session, SEXP term, SEXP config) {
    RssnAnswer answer;
    memset(&answer, 0, sizeof(answer));
    check_status(api.compute(as_handle(session, "session"), (uint32_t)Rf_asInteger(term),
                             optional_handle(config), &answer));

    SEXP out = PROTECT(Rf_allocVector(VECSXP, 6));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 6));
    const char *labels[] = {"term", "has_value", "value", "error", "reduced", "iterations"};
    for (int index = 0; index < 6; index++) SET_STRING_ELT(names, index, Rf_mkChar(labels[index]));
    SET_VECTOR_ELT(out, 0, Rf_ScalarInteger((int)answer.term));
    SET_VECTOR_ELT(out, 1, Rf_ScalarLogical(answer.has_value ? TRUE : FALSE));
    SET_VECTOR_ELT(out, 2, Rf_ScalarReal(answer.has_value ? answer.value : NA_REAL));
    SET_VECTOR_ELT(out, 3, Rf_ScalarReal(answer.has_value ? answer.error : NA_REAL));
    SET_VECTOR_ELT(out, 4, Rf_ScalarLogical(answer.reduced ? TRUE : FALSE));
    SET_VECTOR_ELT(out, 5, Rf_ScalarReal((double)answer.iterations));
    Rf_setAttrib(out, R_NamesSymbol, names);
    UNPROTECT(2);
    return out;
}

typedef struct BatchMetadata {
    RssnAnswerMetadata *rows;
    size_t count;
    SEXP dimensions;
    SEXP findings;
    SEXP messages;
} BatchMetadata;

static SEXP copy_batch_metadata(void *data) {
    BatchMetadata *batch = data;
    for (size_t index = 0; index < batch->count; index++) {
        RssnAnswerMetadata *row = &batch->rows[index];
        SET_STRING_ELT(batch->dimensions, (R_xlen_t)index,
                       Rf_mkCharCE(row->dimension ? row->dimension : "", CE_UTF8));
        SET_STRING_ELT(batch->findings, (R_xlen_t)index,
                       Rf_mkCharCE(row->diagnostic_codes ? row->diagnostic_codes : "", CE_UTF8));
        SET_STRING_ELT(batch->messages, (R_xlen_t)index,
                       Rf_mkCharCE(row->message ? row->message : "", CE_UTF8));
    }
    return R_NilValue;
}

static void release_batch_metadata(void *data) {
    BatchMetadata *batch = data;
    for (size_t index = 0; index < batch->count; index++) {
        api.string_free(batch->rows[index].dimension);
        api.string_free(batch->rows[index].diagnostic_codes);
        api.string_free(batch->rows[index].message);
    }
}

SEXP C_rssn_compute_all(SEXP session, SEXP handles, SEXP config) {
    size_t count = 0;
    uint32_t *buffer = handle_buffer(handles, &count);
    RssnAnswer *answers = (RssnAnswer *)R_alloc(count > 0 ? count : 1, sizeof(RssnAnswer));
    int *statuses = (int *)R_alloc(count > 0 ? count : 1, sizeof(int));
    memset(answers, 0, (count > 0 ? count : 1) * sizeof(RssnAnswer));
    void *session_handle = as_handle(session, "session");
    void *config_handle = optional_handle(config);
    const bool with_metadata = api.compute_all_with_metadata != NULL;

    SEXP terms = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t)count));
    SEXP has_value = PROTECT(Rf_allocVector(LGLSXP, (R_xlen_t)count));
    SEXP values = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)count));
    SEXP errors = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)count));
    SEXP reduced = PROTECT(Rf_allocVector(LGLSXP, (R_xlen_t)count));
    SEXP iterations = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)count));
    SEXP codes = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t)count));
    SEXP dimensions = PROTECT(with_metadata ? Rf_allocVector(STRSXP, (R_xlen_t)count) : R_NilValue);
    SEXP findings = PROTECT(with_metadata ? Rf_allocVector(STRSXP, (R_xlen_t)count) : R_NilValue);
    SEXP messages = PROTECT(with_metadata ? Rf_allocVector(STRSXP, (R_xlen_t)count) : R_NilValue);
    if (with_metadata) {
        RssnAnswerMetadata *rows = (RssnAnswerMetadata *)R_alloc(
            count > 0 ? count : 1, sizeof(RssnAnswerMetadata));
        check_status(api.compute_all_with_metadata(session_handle, buffer, count, config_handle,
                                                   answers, statuses, rows));
        BatchMetadata batch = {rows, count, dimensions, findings, messages};
        R_ExecWithCleanup(copy_batch_metadata, &batch, release_batch_metadata, &batch);
    } else {
        check_status(api.compute_all(session_handle, buffer, count, config_handle, answers, statuses));
    }
    for (size_t index = 0; index < count; index++) {
        INTEGER(terms)[index] = (int)answers[index].term;
        LOGICAL(has_value)[index] = answers[index].has_value ? TRUE : FALSE;
        REAL(values)[index] = answers[index].has_value ? answers[index].value : NA_REAL;
        REAL(errors)[index] = answers[index].has_value ? answers[index].error : NA_REAL;
        LOGICAL(reduced)[index] = answers[index].reduced ? TRUE : FALSE;
        REAL(iterations)[index] = (double)answers[index].iterations;
        INTEGER(codes)[index] = statuses[index];
    }

    SEXP out = PROTECT(Rf_allocVector(VECSXP, 10));
    SET_VECTOR_ELT(out, 0, terms);
    SET_VECTOR_ELT(out, 1, has_value);
    SET_VECTOR_ELT(out, 2, values);
    SET_VECTOR_ELT(out, 3, errors);
    SET_VECTOR_ELT(out, 4, reduced);
    SET_VECTOR_ELT(out, 5, iterations);
    SET_VECTOR_ELT(out, 6, codes);
    SET_VECTOR_ELT(out, 7, dimensions);
    SET_VECTOR_ELT(out, 8, findings);
    SET_VECTOR_ELT(out, 9, messages);
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 10));
    const char *labels[] = {"term", "has_value", "value", "error", "reduced", "iterations",
                            "status", "dimension", "diagnostic_codes", "message"};
    for (int index = 0; index < 10; index++) SET_STRING_ELT(names, index, Rf_mkChar(labels[index]));
    Rf_setAttrib(out, R_NamesSymbol, names);
    UNPROTECT(12);
    return out;
}

SEXP C_rssn_compute_json(SEXP session, SEXP term, SEXP config) {
    char *text = NULL;
    check_status(api.compute_json(as_handle(session, "session"), (uint32_t)Rf_asInteger(term),
                                  optional_handle(config), &text));
    return take_string(text);
}

SEXP C_rssn_compute_bincode(SEXP session, SEXP term, SEXP config) {
    uint8_t *data = NULL;
    size_t length = 0;
    check_status(api.compute_bincode(as_handle(session, "session"), (uint32_t)Rf_asInteger(term),
                                     optional_handle(config), &data, &length));
    return take_bytes(data, length);
}

SEXP C_rssn_compute_all_bincode(SEXP session, SEXP handles, SEXP config) {
    size_t count = 0;
    uint32_t *buffer = handle_buffer(handles, &count);
    uint8_t *data = NULL;
    size_t length = 0;
    check_status(api.compute_all_bincode(as_handle(session, "session"), buffer, count,
                                         optional_handle(config), &data, &length));
    return take_bytes(data, length);
}

SEXP C_rssn_simplify(SEXP session, SEXP source) {
    return term_handle(api.simplify(as_handle(session, "session"), as_cstring(source)), 3);
}

SEXP C_rssn_dimension(SEXP session, SEXP term, SEXP config) {
    char *text = NULL;
    check_status(api.dimension(as_handle(session, "session"), (uint32_t)Rf_asInteger(term),
                               optional_handle(config), &text));
    SEXP out = take_string(text);
    return out == R_NilValue ? Rf_mkString("") : out;
}

SEXP C_rssn_dimension_count(void) { return Rf_ScalarInteger((int)api.dimension_count()); }

/* `index` is zero-based, matching the C report. */
SEXP C_rssn_dimension_code(SEXP index) {
    const char *code = api.dimension_code((size_t)Rf_asInteger(index));
    return Rf_mkString(code == NULL ? "" : code);
}

/* ------------------------------------------------------------------------- */
/* values (the cell ABI)                                                     */
/* ------------------------------------------------------------------------- */

SEXP C_rssn_value_with_coherent_dimension(SEXP session, SEXP source, SEXP json) {
    if (api.value_with_coherent_dimension == NULL) {
        Rf_error("this RSSN library lacks coherent dimension cells");
    }
    if (TYPEOF(json) != RAWSXP) Rf_error("dimension JSON must be raw UTF-8 bytes");
    uint32_t cell = 0;
    check_status(api.value_with_coherent_dimension(as_handle(session, "session"),
        (uint32_t)Rf_asInteger(source), RAW(json), (size_t)XLENGTH(json), &cell));
    return Rf_ScalarInteger((int)cell);
}

SEXP C_rssn_value_new_int(SEXP session, SEXP value) {
    uint32_t cell = 0;
    check_status(api.value_new_int(as_handle(session, "session"), as_i64(value), &cell));
    return Rf_ScalarInteger((int)cell);
}

SEXP C_rssn_value_new_int_decimal(SEXP session, SEXP digits) {
    if (api.value_new_int_decimal == NULL) {
        Rf_error("this RSSN library only supports signed 64-bit integer cells");
    }
    uint32_t cell = 0;
    check_status(api.value_new_int_decimal(as_handle(session, "session"), as_cstring(digits), &cell));
    return Rf_ScalarInteger((int)cell);
}

SEXP C_rssn_value_new_float(SEXP session, SEXP value) {
    uint32_t cell = 0;
    check_status(api.value_new_float(as_handle(session, "session"), Rf_asReal(value), &cell));
    return Rf_ScalarInteger((int)cell);
}

SEXP C_rssn_value_new_text(SEXP session, SEXP text) {
    uint32_t cell = 0;
    check_status(api.value_new_text(as_handle(session, "session"), as_cstring(text), &cell));
    return Rf_ScalarInteger((int)cell);
}

SEXP C_rssn_value_new_symbol(SEXP session, SEXP name) {
    uint32_t cell = 0;
    check_status(api.value_new_symbol(as_handle(session, "session"), as_cstring(name), &cell));
    return Rf_ScalarInteger((int)cell);
}

SEXP C_rssn_value_new_bool(SEXP session, SEXP value) {
    uint32_t cell = 0;
    check_status(api.value_new_bool(as_handle(session, "session"), Rf_asLogical(value), &cell));
    return Rf_ScalarInteger((int)cell);
}

SEXP C_rssn_value_new_complex(SEXP session, SEXP real_part, SEXP imaginary_part) {
    uint32_t cell = 0;
    check_status(api.value_new_complex(as_handle(session, "session"), Rf_asReal(real_part),
                                       Rf_asReal(imaginary_part), &cell));
    return Rf_ScalarInteger((int)cell);
}

SEXP C_rssn_value_new_term(SEXP session, SEXP term) {
    uint32_t cell = 0;
    check_status(api.value_new_term(as_handle(session, "session"), (uint32_t)Rf_asInteger(term),
                                    &cell));
    return Rf_ScalarInteger((int)cell);
}

SEXP C_rssn_value_new_list(SEXP session, SEXP handles) {
    size_t count = 0;
    uint32_t *buffer = handle_buffer(handles, &count);
    uint32_t cell = 0;
    check_status(api.value_new_list(as_handle(session, "session"), buffer, count, &cell));
    return Rf_ScalarInteger((int)cell);
}

/* `data` is row-major; `dim` holds the extents in reading order. */
SEXP C_rssn_value_new_tensor(SEXP session, SEXP data, SEXP dim) {
    size_t rank = (size_t)Rf_xlength(dim);
    size_t *shape = (size_t *)R_alloc(rank > 0 ? rank : 1, sizeof(size_t));
    size_t length = 1;
    for (size_t index = 0; index < rank; index++) {
        shape[index] = (size_t)INTEGER(dim)[index];
        length *= shape[index];
    }
    if ((size_t)Rf_xlength(data) != length) {
        Rf_error("the data length does not match the extents");
    }
    uint32_t cell = 0;
    check_status(api.value_new_tensor(as_handle(session, "session"), REAL(data), length, shape,
                                      rank, &cell));
    return Rf_ScalarInteger((int)cell);
}

SEXP C_rssn_value_new_tensor_cells(SEXP session, SEXP handles, SEXP dim, SEXP exemplar) {
    if (api.value_new_tensor_cells == NULL) {
        Rf_error("this RSSN library lacks tensor cell construction");
    }
    if (TYPEOF(dim) != INTSXP) Rf_error("tensor extents must be integers");
    size_t rank = (size_t)XLENGTH(dim);
    size_t *shape = (size_t *)R_alloc(rank > 0 ? rank : 1, sizeof(size_t));
    for (size_t axis = 0; axis < rank; axis++) {
        int extent = INTEGER(dim)[axis];
        if (extent == NA_INTEGER || extent < 0) Rf_error("invalid tensor extent");
        shape[axis] = (size_t)extent;
    }
    size_t count = 0;
    uint32_t *cells = handle_buffer(handles, &count);
    uint32_t output = 0;
    check_status(api.value_new_tensor_cells(as_handle(session, "session"), cells, count, shape,
                                            rank, (uint32_t)Rf_asInteger(exemplar), &output));
    return Rf_ScalarInteger((int)output);
}

SEXP C_rssn_value_new_record(SEXP session, SEXP names, SEXP handles) {
    size_t count = 0;
    uint32_t *cells = handle_buffer(handles, &count);
    const char **fields = (const char **)R_alloc(count > 0 ? count : 1, sizeof(char *));
    for (size_t index = 0; index < count; index++) {
        fields[index] = Rf_translateCharUTF8(STRING_ELT(names, (R_xlen_t)index));
    }
    uint32_t cell = 0;
    check_status(api.value_new_record(as_handle(session, "session"), fields, cells, count, &cell));
    return Rf_ScalarInteger((int)cell);
}

SEXP C_rssn_value_call(SEXP session, SEXP op, SEXP handles) {
    size_t count = 0;
    uint32_t *buffer = handle_buffer(handles, &count);
    uint32_t cell = 0;
    check_status(api.value_call(as_handle(session, "session"), as_cstring(op), buffer, count,
                                &cell));
    return Rf_ScalarInteger((int)cell);
}

SEXP C_rssn_value_call_key(SEXP session, SEXP rule_set, SEXP name, SEXP hash, SEXP handles) {
    if (api.value_call_key == NULL) Rf_error("this RSSN library lacks keyed calls");
    if (TYPEOF(hash) != RAWSXP || XLENGTH(hash) != 16) {
        Rf_error("operator key requires 16 signature hash bytes");
    }
    size_t count = 0;
    uint32_t *buffer = handle_buffer(handles, &count);
    uint32_t result = 0;
    check_status(api.value_call_key(as_handle(session, "session"), as_cstring(rule_set),
                                    as_cstring(name), RAW(hash), buffer, count, &result));
    return Rf_ScalarInteger((int)result);
}

SEXP C_rssn_value_call_json(SEXP session, SEXP op, SEXP handles) {
    size_t count = 0;
    uint32_t *buffer = handle_buffer(handles, &count);
    char *text = NULL;
    check_status(api.value_call_json(as_handle(session, "session"), as_cstring(op), buffer, count,
                                     &text));
    return take_string(text);
}

SEXP C_rssn_value_kind(SEXP session, SEXP cell) {
    uint8_t kind = 0;
    check_status(api.value_kind(as_handle(session, "session"), (uint32_t)Rf_asInteger(cell), &kind));
    return Rf_ScalarInteger((int)kind);
}

SEXP C_rssn_value_to_f64(SEXP session, SEXP cell) {
    double out = 0.0;
    check_status(api.value_to_f64(as_handle(session, "session"), (uint32_t)Rf_asInteger(cell),
                                  &out));
    return Rf_ScalarReal(out);
}

SEXP C_rssn_value_to_text(SEXP session, SEXP cell) {
    char *text = NULL;
    check_status(api.value_to_text(as_handle(session, "session"), (uint32_t)Rf_asInteger(cell),
                                   &text));
    SEXP out = take_string(text);
    return out == R_NilValue ? Rf_mkString("") : out;
}

SEXP C_rssn_value_to_term(SEXP session, SEXP cell) {
    uint32_t term = 0;
    check_status(api.value_to_term(as_handle(session, "session"), (uint32_t)Rf_asInteger(cell),
                                   &term));
    return Rf_ScalarInteger((int)term);
}

SEXP C_rssn_value_to_tensor(SEXP session, SEXP cell) {
    RssnTensor tensor;
    memset(&tensor, 0, sizeof(tensor));
    check_status(api.value_to_tensor(as_handle(session, "session"), (uint32_t)Rf_asInteger(cell),
                                     &tensor));
    return take_tensor(&tensor);
}

SEXP C_rssn_value_len(SEXP session, SEXP cell) {
    size_t length = 0;
    check_status(api.value_len(as_handle(session, "session"), (uint32_t)Rf_asInteger(cell),
                               &length));
    return Rf_ScalarInteger((int)length);
}

SEXP C_rssn_value_at(SEXP session, SEXP cell, SEXP index) {
    uint32_t out = 0;
    check_status(api.value_at(as_handle(session, "session"), (uint32_t)Rf_asInteger(cell),
                              (size_t)Rf_asInteger(index), &out));
    return Rf_ScalarInteger((int)out);
}

SEXP C_rssn_value_field_name(SEXP session, SEXP cell, SEXP index) {
    char *name = NULL;
    check_status(api.value_field_name(as_handle(session, "session"),
                                      (uint32_t)Rf_asInteger(cell),
                                      (size_t)Rf_asInteger(index), &name));
    return take_string(name);
}

SEXP C_rssn_value_tag(SEXP session, SEXP cell) {
    char *tag = NULL;
    check_status(api.value_tag(as_handle(session, "session"), (uint32_t)Rf_asInteger(cell), &tag));
    return take_string(tag);
}

SEXP C_rssn_value_json(SEXP session, SEXP cell) {
    char *text = NULL;
    check_status(api.value_json(as_handle(session, "session"), (uint32_t)Rf_asInteger(cell),
                                &text));
    return take_string(text);
}

SEXP C_rssn_value_clone(SEXP session, SEXP cell) {
    void *owner = as_handle(session, "session");
    uint32_t source = (uint32_t)Rf_asInteger(cell);
    uint32_t copy = 0;
    if (api.value_clone != NULL) {
        check_status(api.value_clone(owner, source, &copy));
    } else {
        uint32_t holder = 0;
        check_status(api.value_new_list(owner, &source, 1, &holder));
        int status = api.value_at(owner, holder, 0, &copy);
        api.value_free(owner, holder);
        check_status(status);
    }
    return Rf_ScalarInteger((int)copy);
}

SEXP C_rssn_value_free(SEXP session, SEXP cell) {
    api.value_free(as_handle(session, "session"), (uint32_t)Rf_asInteger(cell));
    return R_NilValue;
}

SEXP C_rssn_value_count(SEXP session) {
    return Rf_ScalarInteger((int)api.value_count(as_handle(session, "session")));
}

/* ------------------------------------------------------------------------- */
/* host callbacks                                                            */
/* ------------------------------------------------------------------------- */

/*
 * The registry maps the integer the engine carries as `user` to the R closure
 * that answers the call. Slots are reused, so a long run of one-shot callbacks
 * does not grow the table.
 */
static SEXP *registry_closures = NULL;
static size_t registry_capacity = 0;

static uintptr_t registry_add(SEXP closure) {
    for (size_t index = 0; index < registry_capacity; index++) {
        if (registry_closures[index] == NULL) {
            registry_closures[index] = closure;
            R_PreserveObject(closure);
            return (uintptr_t)(index + 1);
        }
    }
    size_t next = registry_capacity + 16;
    SEXP *grown = (SEXP *)realloc(registry_closures, next * sizeof(SEXP));
    if (grown == NULL) Rf_error("cannot grow the callback registry");
    for (size_t index = registry_capacity; index < next; index++) grown[index] = NULL;
    registry_closures = grown;
    registry_capacity = next;
    registry_closures[0] = closure;
    R_PreserveObject(closure);
    return (uintptr_t)1;
}

static SEXP registry_get(uintptr_t id) {
    if (id == 0 || id > registry_capacity) return NULL;
    return registry_closures[id - 1];
}

static void registry_remove(uintptr_t id) {
    if (id == 0 || id > registry_capacity) return;
    SEXP closure = registry_closures[id - 1];
    if (closure != NULL) {
        registry_closures[id - 1] = NULL;
        R_ReleaseObject(closure);
    }
}

/*
 * The trampoline. It is a plain C function, so it captures nothing: the R
 * closure comes from the registry through `user`. The closure receives the
 * argument cell ids as an integer vector and answers one integer: a cell id,
 * or a negative status after it has reported the message itself.
 *
 * `R_tryEval` catches an R error here rather than letting it unwind through the
 * engine's Rust frames, which would be undefined behaviour.
 */
static int r_callback_thunk(void *user, const void *session, const uint32_t *args, size_t count,
                            uint32_t *out) {
    SEXP closure = registry_get((uintptr_t)user);
    if (closure == NULL) return 13;

    SEXP cells = PROTECT(Rf_allocVector(INTSXP, (R_xlen_t)count));
    for (size_t index = 0; index < count; index++) INTEGER(cells)[index] = (int)args[index];
    SEXP call = PROTECT(Rf_lang2(closure, cells));
    int failed = 0;
    SEXP result = R_tryEval(call, R_GlobalEnv, &failed);
    if (failed || result == R_NilValue || !Rf_isInteger(result) || Rf_xlength(result) < 1) {
        UNPROTECT(2);
        api.callback_fail("the R callback did not answer a cell id");
        return 13;
    }
    PROTECT(result);
    int code = INTEGER(result)[0];
    UNPROTECT(3);
    if (code < 0) return 13;
    *out = (uint32_t)code;
    return 0;
}

/*
 * Registers `closure` as a host callback. `signature` is a JSON type expression
 * or `R_NilValue`; `op`/`index` name the argument of an operator instead.
 *
 * The answer is `c(cell, registry_id)`: the cell is what the engine takes as an
 * argument, and the id is what releases the closure once the callback is
 * retired.
 */
static SEXP callback_register(SEXP session, SEXP signature, SEXP op, SEXP index, SEXP closure) {
    if (!Rf_isFunction(closure)) Rf_error("a callback must be a function");
    uintptr_t id = registry_add(closure);
    uint32_t cell = 0;
    int code;
    if (signature != R_NilValue) {
        code = api.callback_new(as_handle(session, "session"), as_cstring(signature),
                                r_callback_thunk, (void *)id, &cell);
    } else if (op != R_NilValue) {
        code = api.callback_new_for(as_handle(session, "session"), as_cstring(op),
                                    (size_t)Rf_asInteger(index), r_callback_thunk, (void *)id,
                                    &cell);
    } else {
        registry_remove(id);
        Rf_error("a callback needs a signature, or an operator and an argument position");
    }
    if (code != RSSN_OK) {
        registry_remove(id);
        raise_status(code, "the callback could not be registered");
    }
    SEXP out = PROTECT(Rf_allocVector(INTSXP, 2));
    INTEGER(out)[0] = (int)cell;
    INTEGER(out)[1] = (int)id;
    UNPROTECT(1);
    return out;
}

SEXP C_rssn_callback_new(SEXP session, SEXP signature, SEXP closure) {
    return callback_register(session, signature, R_NilValue, R_NilValue, closure);
}

SEXP C_rssn_callback_new_for(SEXP session, SEXP op, SEXP index, SEXP closure) {
    return callback_register(session, R_NilValue, op, index, closure);
}

SEXP C_rssn_callback_signatures(SEXP session, SEXP op, SEXP index) {
    char *text = NULL;
    check_status(api.callback_signatures(as_handle(session, "session"), as_cstring(op),
                                         (size_t)Rf_asInteger(index), &text));
    SEXP out = take_string(text);
    return out == R_NilValue ? Rf_mkString("[]") : out;
}

SEXP C_rssn_callback_free(SEXP session, SEXP cell) {
    api.callback_free(as_handle(session, "session"), (uint32_t)Rf_asInteger(cell));
    return R_NilValue;
}

SEXP C_rssn_callback_count(SEXP session) {
    return Rf_ScalarInteger((int)api.callback_count(as_handle(session, "session")));
}

SEXP C_rssn_callback_fail(SEXP message) {
    api.callback_fail(as_cstring(message));
    return R_NilValue;
}

/* Releases a registry slot when a callback object is collected. */
SEXP C_rssn_registry_release(SEXP id) {
    registry_remove((uintptr_t)Rf_asReal(id));
    return R_NilValue;
}

/* ------------------------------------------------------------------------- */
/* simulations                                                               */
/* ------------------------------------------------------------------------- */

SEXP C_rssn_sim_run(SEXP name, SEXP parameters) {
    char *text = NULL;
    check_status(api.sim_run(as_cstring(name), as_cstring(parameters), &text));
    return take_string(text);
}

/* ------------------------------------------------------------------------- */
/* registration                                                              */
/* ------------------------------------------------------------------------- */

static const R_CallMethodDef call_methods[] = {
    {"C_rssn_load", (DL_FUNC)&C_rssn_load, 1},
    {"C_rssn_libpath", (DL_FUNC)&C_rssn_libpath, 0},
    {"C_rssn_version", (DL_FUNC)&C_rssn_version, 0},
    {"C_rssn_last_error", (DL_FUNC)&C_rssn_last_error, 0},
    {"C_rssn_manifest_hash", (DL_FUNC)&C_rssn_manifest_hash, 0},
    {"C_rssn_manifest_json", (DL_FUNC)&C_rssn_manifest_json, 0},
    {"C_rssn_op_info_json", (DL_FUNC)&C_rssn_op_info_json, 1},

    {"C_rssn_session_new", (DL_FUNC)&C_rssn_session_new, 0},
    {"C_rssn_session_free", (DL_FUNC)&C_rssn_session_free, 1},
    {"C_rssn_config_new", (DL_FUNC)&C_rssn_config_new, 0},
    {"C_rssn_config_free", (DL_FUNC)&C_rssn_config_free, 1},
    {"C_rssn_config_symbolic", (DL_FUNC)&C_rssn_config_symbolic, 1},
    {"C_rssn_config_numeric", (DL_FUNC)&C_rssn_config_numeric, 2},
    {"C_rssn_config_bind", (DL_FUNC)&C_rssn_config_bind, 3},
    {"C_rssn_config_assume", (DL_FUNC)&C_rssn_config_assume, 3},
    {"C_rssn_config_budget", (DL_FUNC)&C_rssn_config_budget, 3},
    {"C_rssn_config_dimension_mode", (DL_FUNC)&C_rssn_config_dimension_mode, 2},
    {"C_rssn_config_bypass_dimensions", (DL_FUNC)&C_rssn_config_bypass_dimensions, 1},
    {"C_rssn_config_dimension", (DL_FUNC)&C_rssn_config_dimension, 3},
    {"C_rssn_config_define_base_dimension", (DL_FUNC)&C_rssn_config_define_base_dimension, 2},

    {"C_rssn_term_parse", (DL_FUNC)&C_rssn_term_parse, 2},
    {"C_rssn_term_sym", (DL_FUNC)&C_rssn_term_sym, 2},
    {"C_rssn_term_int", (DL_FUNC)&C_rssn_term_int, 2},
    {"C_rssn_term_float", (DL_FUNC)&C_rssn_term_float, 2},
    {"C_rssn_term_rational", (DL_FUNC)&C_rssn_term_rational, 3},
    {"C_rssn_term_apply", (DL_FUNC)&C_rssn_term_apply, 3},
    {"C_rssn_term_eval", (DL_FUNC)&C_rssn_term_eval, 4},
    {"C_rssn_term_to_f64", (DL_FUNC)&C_rssn_term_to_f64, 2},
    {"C_rssn_term_tensor", (DL_FUNC)&C_rssn_term_tensor, 2},
    {"C_rssn_term_to_string", (DL_FUNC)&C_rssn_term_to_string, 2},
    {"C_rssn_term_to_latex", (DL_FUNC)&C_rssn_term_to_latex, 2},
    {"C_rssn_term_to_typst", (DL_FUNC)&C_rssn_term_to_typst, 2},
    {"C_rssn_term_to_pretty", (DL_FUNC)&C_rssn_term_to_pretty, 2},

    {"C_rssn_compute", (DL_FUNC)&C_rssn_compute, 3},
    {"C_rssn_compute_all", (DL_FUNC)&C_rssn_compute_all, 3},
    {"C_rssn_compute_json", (DL_FUNC)&C_rssn_compute_json, 3},
    {"C_rssn_compute_bincode", (DL_FUNC)&C_rssn_compute_bincode, 3},
    {"C_rssn_compute_all_bincode", (DL_FUNC)&C_rssn_compute_all_bincode, 3},
    {"C_rssn_simplify", (DL_FUNC)&C_rssn_simplify, 2},
    {"C_rssn_dimension", (DL_FUNC)&C_rssn_dimension, 3},
    {"C_rssn_dimension_count", (DL_FUNC)&C_rssn_dimension_count, 0},
    {"C_rssn_dimension_code", (DL_FUNC)&C_rssn_dimension_code, 1},

    {"C_rssn_value_with_coherent_dimension", (DL_FUNC)&C_rssn_value_with_coherent_dimension, 3},
    {"C_rssn_value_new_int", (DL_FUNC)&C_rssn_value_new_int, 2},
    {"C_rssn_value_new_int_decimal", (DL_FUNC)&C_rssn_value_new_int_decimal, 2},
    {"C_rssn_value_new_float", (DL_FUNC)&C_rssn_value_new_float, 2},
    {"C_rssn_value_new_text", (DL_FUNC)&C_rssn_value_new_text, 2},
    {"C_rssn_value_new_symbol", (DL_FUNC)&C_rssn_value_new_symbol, 2},
    {"C_rssn_value_new_bool", (DL_FUNC)&C_rssn_value_new_bool, 2},
    {"C_rssn_value_new_complex", (DL_FUNC)&C_rssn_value_new_complex, 3},
    {"C_rssn_value_new_term", (DL_FUNC)&C_rssn_value_new_term, 2},
    {"C_rssn_value_new_list", (DL_FUNC)&C_rssn_value_new_list, 2},
    {"C_rssn_value_new_tensor", (DL_FUNC)&C_rssn_value_new_tensor, 3},
    {"C_rssn_value_new_tensor_cells", (DL_FUNC)&C_rssn_value_new_tensor_cells, 4},
    {"C_rssn_value_new_record", (DL_FUNC)&C_rssn_value_new_record, 3},
    {"C_rssn_value_call", (DL_FUNC)&C_rssn_value_call, 3},
    {"C_rssn_value_call_key", (DL_FUNC)&C_rssn_value_call_key, 5},
    {"C_rssn_value_call_json", (DL_FUNC)&C_rssn_value_call_json, 3},
    {"C_rssn_value_kind", (DL_FUNC)&C_rssn_value_kind, 2},
    {"C_rssn_value_to_f64", (DL_FUNC)&C_rssn_value_to_f64, 2},
    {"C_rssn_value_to_text", (DL_FUNC)&C_rssn_value_to_text, 2},
    {"C_rssn_value_to_term", (DL_FUNC)&C_rssn_value_to_term, 2},
    {"C_rssn_value_to_tensor", (DL_FUNC)&C_rssn_value_to_tensor, 2},
    {"C_rssn_value_len", (DL_FUNC)&C_rssn_value_len, 2},
    {"C_rssn_value_at", (DL_FUNC)&C_rssn_value_at, 3},
    {"C_rssn_value_field_name", (DL_FUNC)&C_rssn_value_field_name, 3},
    {"C_rssn_value_tag", (DL_FUNC)&C_rssn_value_tag, 2},
    {"C_rssn_value_json", (DL_FUNC)&C_rssn_value_json, 2},
    {"C_rssn_value_clone", (DL_FUNC)&C_rssn_value_clone, 2},
    {"C_rssn_value_free", (DL_FUNC)&C_rssn_value_free, 2},
    {"C_rssn_value_count", (DL_FUNC)&C_rssn_value_count, 1},

    {"C_rssn_callback_new", (DL_FUNC)&C_rssn_callback_new, 3},
    {"C_rssn_callback_new_for", (DL_FUNC)&C_rssn_callback_new_for, 4},
    {"C_rssn_callback_signatures", (DL_FUNC)&C_rssn_callback_signatures, 3},
    {"C_rssn_callback_free", (DL_FUNC)&C_rssn_callback_free, 2},
    {"C_rssn_callback_count", (DL_FUNC)&C_rssn_callback_count, 1},
    {"C_rssn_callback_fail", (DL_FUNC)&C_rssn_callback_fail, 1},
    {"C_rssn_registry_release", (DL_FUNC)&C_rssn_registry_release, 1},

    {"C_rssn_sim_run", (DL_FUNC)&C_rssn_sim_run, 2},
    {NULL, NULL, 0},
};

void R_init_rssn(DllInfo *info) {
    R_registerRoutines(info, NULL, call_methods, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
    R_forceSymbols(info, TRUE);
}
