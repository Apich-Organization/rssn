#define _GNU_SOURCE
#include <dlfcn.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

// Run the single-threaded benchmark with this probe separately from timing.
#ifdef RSSN_COUNT_POW
static double (*real_pow)(double, double);
static uint64_t calls, unit_bases;

double pow(double x, double y) {
    if (real_pow == NULL) {
        real_pow = dlsym(RTLD_NEXT, "pow");
        if (real_pow == NULL) {
            fputs("unable to resolve the next pow\n", stderr);
            abort();
        }
    }
    ++calls;
    unit_bases += x == 1.0;
    return real_pow(x, y);
}

__attribute__((destructor))
static void report(void) {
    fprintf(stderr, "%llu %llu\n", (unsigned long long)calls,
            (unsigned long long)unit_bases);
}
#else
static double (*real_exp)(double);
static uint64_t calls;

double exp(double x) {
    if (real_exp == NULL) {
        real_exp = dlsym(RTLD_NEXT, "exp");
        if (real_exp == NULL) {
            fputs("unable to resolve the next exp\n", stderr);
            abort();
        }
    }
    ++calls;
    return real_exp(x);
}

#ifdef RSSN_COUNT_LOG
static double (*real_log)(double);
static uint64_t log_calls, log_ones;

double log(double x) {
    if (real_log == NULL) {
        real_log = dlsym(RTLD_NEXT, "log");
        if (real_log == NULL) {
            fputs("unable to resolve the next log\n", stderr);
            abort();
        }
    }
    ++log_calls;
    log_ones += x == 1.0;
    return real_log(x);
}

__attribute__((destructor))
static void report(void) {
    fprintf(stderr, "%llu %llu %llu\n", (unsigned long long)calls,
            (unsigned long long)log_calls, (unsigned long long)log_ones);
}
#else
__attribute__((destructor))
static void report(void) {
    fprintf(stderr, "%llu\n", (unsigned long long)calls);
}
#endif
#endif
