/* Minimal PE runner for Appell and Lerch kernels extracted from rustc's fat LTO IR. */
typedef unsigned long long U64;
typedef unsigned long Dword;
typedef struct {
    Dword low;
    Dword high;
} FileTime;
#define IMPORT __declspec(dllimport)

IMPORT int QueryPerformanceCounter(long long *);
IMPORT int QueryPerformanceFrequency(long long *);
IMPORT void *GetCurrentThread(void);
IMPORT int GetThreadTimes(void *, FileTime *, FileTime *, FileTime *, FileTime *);
IMPORT void *GetStdHandle(Dword);
IMPORT int WriteFile(void *, const void *, Dword, Dword *, void *);
IMPORT __declspec(noreturn) void ExitProcess(Dword);
IMPORT Dword GetEnvironmentVariableA(const char *, char *, Dword);
IMPORT void *LoadLibraryA(const char *);
IMPORT void (*GetProcAddress(void *, const char *))(void);
int _fltused;

/* These are the same integration cases as the Rust benchmark. */
#if defined(RSSN_WIN_LERCH_SMALL) || defined(RSSN_WIN_LERCH_POSITIVE) || defined(RSSN_WIN_LERCH_NEGATIVE)
typedef double (*Kernel)(double, double, double);
#ifdef RSSN_WIN_LERCH_POSITIVE
static const char symbol[] = "rssn_de_positive";
static const double cases[][3] = {
    {-2.934, 0.805, 1.0},
    {0.999, 0.805, 1.0},
    {-6.5, 0.5, 1.7},
    {0.999, 0.15, 0.001},
};
#elif defined(RSSN_WIN_LERCH_NEGATIVE)
static const char symbol[] = "rssn_de_negative";
static const double cases[][3] = {
    {-1.0, -0.25, 1.0},
    {0.95, -0.5, 1.2},
    {-4.0, -3.5, 0.7},
    {0.99, -2.25, 1.7},
};
#else
static const char symbol[] = "rssn_de_small";
static const double cases[][3] = {
    {-2.934, 0.09, 1.0},
    {0.999, 0.001, 0.001},
    {-1234.0, 0.1, 1.2},
    {0.9999999, 0.12, 1.0},
};
#endif
#else
typedef double (*Kernel)(double, double, double, double, double, double);
static const char symbol[] = "rssn_de_appell";
static const double cases[][6] = {
    {0.8, 1.25, 0.5, 2.1, 0.96, 0.3},
    {1.2, 0.25, 2.0, 2.5, 0.4, 0.99},
    {0.05, 2.0, 0.25, 0.9, -1.2, 0.97},
    {1.55, 0.05, 0.25, 2.1, -2.1, -3.1},
};
#endif

static U64 evaluate(Kernel kernel, const double *x) {
    union {
        double number;
        U64 bits;
    } value;
#if defined(RSSN_WIN_LERCH_SMALL) || defined(RSSN_WIN_LERCH_POSITIVE) || defined(RSSN_WIN_LERCH_NEGATIVE)
    value.number = kernel(x[0], x[1], x[2]);
#else
    value.number = kernel(x[0], x[1], x[2], x[3], x[4], x[5]);
#endif
    return value.bits;
}

static unsigned get_number(const char *name, unsigned fallback) {
    char argument[12];
    Dword length = GetEnvironmentVariableA(name, argument, sizeof argument);
    if (length == 0)
        return fallback;
    if (length > 5)
        ExitProcess(2);
    unsigned number = 0;
    for (unsigned i = 0; i < length; ++i) {
        if (argument[i] < '0' || argument[i] > '9')
            ExitProcess(2);
        number = 10 * number + (argument[i] - '0');
    }
    return number;
}

static unsigned append_decimal(char *text, unsigned used, U64 number) {
    U64 divisor = 1;
    while (number / divisor >= 10)
        divisor *= 10;
    for (; divisor; divisor /= 10)
        text[used++] = '0' + (char)(number / divisor % 10);
    return used;
}

static U64 ticks(FileTime time) {
    return ((U64)time.high << 32) | time.low;
}

static U64 thread_cpu_ns(void) {
    FileTime creation, exit, kernel, user;
    if (!GetThreadTimes(GetCurrentThread(), &creation, &exit, &kernel, &user))
        ExitProcess(4);
    return 100 * (ticks(kernel) + ticks(user));
}

static void display(const char *mode, U64 ns, U64 cpu_ns, U64 hash) {
    char text[68];
    unsigned used = 0;
    for (; *mode; ++mode)
        text[used++] = *mode;
    text[used++] = ' ';
    used = append_decimal(text, used, ns);
    text[used++] = ' ';
    used = append_decimal(text, used, cpu_ns);
    text[used++] = ' ';
    for (unsigned i = 0; i < 16; ++i)
        text[used++] = "0123456789abcdef"[(hash >> (60 - 4 * i)) & 15];
    text[used++] = '\n';
    Dword written = 0;
    if (!WriteFile(GetStdHandle((Dword)-11), text, used, &written, 0) || written != used)
        ExitProcess(5);
}

void mainCRTStartup(void) {
    const unsigned rounds = get_number("RSSN_WIN_ROUNDS", 55);
    const unsigned pairs = get_number("RSSN_WIN_PAIRS", 24);
    const unsigned variants = get_number("RSSN_WIN_VARIANTS", 4);
    const int compression = get_number("RSSN_WIN_EXIT_COMPRESSION", 0) != 0;
    const int endpoint = get_number("RSSN_WIN_FP_ENDPOINT", 0) != 0;
    const int bank = get_number("RSSN_WIN_FP_BANK", 0) != 0;
    const int unit = get_number("RSSN_WIN_POW_UNIT_BASE", 0) != 0;
    const int log_cache = get_number("RSSN_WIN_LIBM_UNIT_LOG", 0) != 0;
    const int tail = get_number("RSSN_WIN_LIBM_TAIL_CUTOFF", 0) != 0;
    if (compression + endpoint + bank + unit + log_cache + tail > 1 ||
        (variants != 2 && variants != 4))
        ExitProcess(2);
    const char *modes[] = {"base", "pack", "reuse", "both"};
    const char *libraries[] = {"base.dll", "pack.dll", "reuse.dll", "both.dll"};
    if (compression) {
        modes[1] = "compression";
        libraries[1] = "compression.dll";
    } else if (endpoint) {
        modes[1] = "endpoint";
        libraries[1] = "endpoint.dll";
    } else if (unit) {
        modes[1] = "unit";
        libraries[1] = "unit.dll";
    } else if (log_cache) {
        modes[1] = "unit-log";
        libraries[1] = "unit-log.dll";
        modes[2] = "endpoint";
        libraries[2] = "endpoint.dll";
    } else if (tail) {
        modes[1] = "tail";
        libraries[1] = "tail.dll";
        modes[2] = "unit-log";
        libraries[2] = "unit-log.dll";
    } else if (bank) {
        modes[1] = "bank";
        libraries[1] = "bank.dll";
#ifdef RSSN_WIN_LERCH_SMALL
        modes[2] = "endpoint";
        libraries[2] = "endpoint.dll";
#endif
    }
    Kernel kernels[4];
    for (unsigned mode = 0; mode < variants; ++mode) {
        void *module = LoadLibraryA(libraries[mode]);
        if (!module)
            ExitProcess(6);
        kernels[mode] = (Kernel)GetProcAddress(module, symbol);
        if (!kernels[mode])
            ExitProcess(6);
    }
    for (unsigned i = 0; i < sizeof cases / sizeof cases[0]; ++i) {
        const U64 reference = evaluate(kernels[0], cases[i]);
        if (((reference >> 52) & 2047) == 2047)
            ExitProcess(3);
        for (unsigned mode = 1; mode < variants; ++mode)
            if (evaluate(kernels[mode], cases[i]) != reference)
                ExitProcess(3);
    }
    long long frequency;
    if (!QueryPerformanceFrequency(&frequency) || frequency <= 0)
        ExitProcess(4);
    const U64 hz = (U64)frequency;
    U64 reference = 0;
    for (unsigned pair = 0; pair < pairs; ++pair) {
        for (unsigned i = 0; i < variants; ++i) {
            const unsigned direction = (pair & variants) ? variants - 1 - i : i;
            const unsigned mode = (pair + direction) % variants;
            long long start, stop;
            const U64 cpu_start = thread_cpu_ns();
            if (!QueryPerformanceCounter(&start))
                ExitProcess(4);
            U64 hash = 0;
            for (unsigned n = 0; n < rounds; ++n)
                for (unsigned j = 0; j < sizeof cases / sizeof cases[0]; ++j) {
                    hash = (((hash << 7) | (hash >> 57)) ^ evaluate(kernels[mode], cases[j]))
                        * 0x9e3779b97f4a7c15ULL;
                }
            if (!QueryPerformanceCounter(&stop) || stop < start)
                ExitProcess(4);
            const U64 cpu_stop = thread_cpu_ns();
            if (cpu_stop < cpu_start)
                ExitProcess(4);
            if (pair == 0 && i == 0)
                reference = hash;
            if (hash != reference)
                ExitProcess(3);
            const U64 ticks = (U64)(stop - start);
            const U64 ns = (ticks / hz) * 1000000000ULL + (ticks % hz) * 1000000000ULL / hz;
            display(modes[mode], ns, cpu_stop - cpu_start, hash);
        }
    }
    ExitProcess(0);
}
