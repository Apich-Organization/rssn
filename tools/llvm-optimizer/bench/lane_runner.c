#define _POSIX_C_SOURCE 200809L
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

enum { ROWS = 256 };
extern double pack_live_scalars(double *);

double opaque_scale(double *input) {
  double scale = 0.5 + input[1] * 0.25;
  input[0] += 0.125;
  input[1] += 0.0625;
  return scale;
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
  double result[ROWS] = {0};
  double input[ROWS][2];
  for (unsigned row = 0; row < ROWS; ++row) {
    input[row][0] = row * (1.0 / 256.0);
    input[row][1] = 1.0 - row * (1.0 / 512.0);
    result[row] = pack_live_scalars(input[row]);
  }

  const uint64_t start = now_ns();
  for (uint64_t i = 0; i < rounds; ++i) {
    const unsigned row = i % ROWS;
    result[row] = pack_live_scalars(input[row]);
  }
  const uint64_t elapsed = now_ns() - start;

  uint64_t hash = 0;
  for (unsigned row = 0; row < ROWS; ++row) {
    uint64_t bits;
    memcpy(&bits, &result[row], sizeof(bits));
    hash = (hash ^ bits) * UINT64_C(0x100000001b3);
    memcpy(&bits, &input[row][0], sizeof(bits));
    hash = (hash ^ bits) * UINT64_C(0x100000001b3);
    memcpy(&bits, &input[row][1], sizeof(bits));
    hash = (hash ^ bits) * UINT64_C(0x100000001b3);
  }
  printf("%" PRIu64 " %016" PRIx64 "\n", elapsed, hash);
}
