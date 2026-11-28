#define _POSIX_C_SOURCE 200809L
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

enum { ELEMENTS = 48, ROWS = 256 };

extern void rssn_callback_kernel(const double *, double *, uint64_t);

double rssn_gain(uint64_t key) {
  key ^= key >> 17;
  key *= UINT64_C(0x9e3779b97f4a7c15);
  return 0.5 + (double)(key >> 52) * (1.0 / 4096.0);
}

static uint64_t now_ns(void) {
  struct timespec ts;
  if (clock_gettime(CLOCK_MONOTONIC, &ts) != 0)
    exit(EXIT_FAILURE);
  return (uint64_t)ts.tv_sec * UINT64_C(1000000000) + ts.tv_nsec;
}

int main(int argc, char **argv) {
  if (argc != 2)
    return EXIT_FAILURE;
  const uint64_t rounds = strtoull(argv[1], NULL, 10);
  static double input[ROWS][ELEMENTS], output[ROWS][ELEMENTS];
  for (unsigned row = 0; row < ROWS; ++row) {
    for (unsigned col = 0; col < ELEMENTS; ++col)
      input[row][col] = (double)((row * 13 + col) % 97) / 97.0;
    rssn_callback_kernel(input[row], output[row], row);
  }

  const uint64_t start = now_ns();
  for (uint64_t i = 0; i < rounds; ++i) {
    const unsigned row = i % ROWS;
    rssn_callback_kernel(input[row], output[row], i);
  }
  const uint64_t elapsed = now_ns() - start;

  uint64_t hash = 0;
  for (unsigned row = 0; row < ROWS; ++row)
    for (unsigned col = 0; col < ELEMENTS; ++col) {
      uint64_t bits;
      memcpy(&bits, &output[row][col], sizeof(bits));
      hash = (hash ^ bits) * UINT64_C(0x100000001b3);
    }
  printf("%" PRIu64 " %016" PRIx64 "\n", elapsed, hash);
}
