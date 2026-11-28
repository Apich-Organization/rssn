// rssn kernel: 2 inputs, 1 params, 1 outputs, 16 ops, fold (float)
#include <math.h>
__device__ int rssn_toidx(float x) {
    if (x != x) { return 0; }
    if (x >= 2147483500.0f) { return 2147483647; }
    if (x <= -2147483500.0f) { return -2147483647; }
    return ((int)(x));
}

__device__ int rssn_idiv(int a, int b) {
    if (b == 0) { return 0; }
    return a / b;
}

__device__ int rssn_imod(int a, int b) {
    if (b == 0) { return 0; }
    return a % b;
}

#define RSSN_LOAD(col, i) (((i) >= 0 && ((unsigned int)(i)) < rssn_cols[2u * (col) + 1u]) ? rssn_in[rssn_cols[2u * (col)] + ((unsigned int)(i))] : ((float)(NAN)))

extern "C" __global__ void rssn_main(const float* rssn_in, const float* rssn_params, float* rssn_out, const unsigned int* rssn_cols, unsigned int rssn_rows, unsigned int rssn_inputs) {
    const unsigned int row = blockIdx.x * blockDim.x + threadIdx.x;
    if (row >= rssn_rows) { return; }
    const int irow = (int)row;
    const float v0 = rssn_params[0u];
    const int v1 = rssn_toidx(v0);
    const int v2 = ((int)(rssn_cols[3u]));
    const int v3 = rssn_idiv(v2, v1);
    const int v4 = irow;
    const int v5 = rssn_idiv(v4, v1);
    const int v6 = rssn_imod(v4, v1);
    const int v7 = v5 * v3;
    float a0 = 0.0f;
    for (int j = 0; j < v3; j = j + 1) {
        const int v8 = j;
        const int v9 = v7 + v8;
        const int v10 = v8 * v1;
        const int v11 = v10 + v6;
        const float v12 = RSSN_LOAD(0u, v9);
        const float v13 = RSSN_LOAD(1u, v11);
        const float v14 = v12 * v13;
        a0 = a0 + v14;
    }
    const float v15 = a0;
    rssn_out[0u * rssn_rows + row] = v15;
}
