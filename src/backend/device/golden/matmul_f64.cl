// rssn kernel: 2 inputs, 1 params, 1 outputs, 16 ops, fold (double)
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
static inline int rssn_toidx(double x) {
    if (x != x) { return 0; }
    if (x >= 2147483520.0) { return 2147483647; }
    if (x <= -2147483520.0) { return -2147483647; }
    return ((int)(x));
}

static inline int rssn_idiv(int a, int b) {
    if (b == 0) { return 0; }
    return a / b;
}

static inline int rssn_imod(int a, int b) {
    if (b == 0) { return 0; }
    return a % b;
}

#define RSSN_LOAD(col, i) (((i) >= 0 && ((uint)(i)) < rssn_cols[2u * (col) + 1u]) ? rssn_in[rssn_cols[2u * (col)] + ((uint)(i))] : ((double)(NAN)))

__kernel void rssn_main(__global const double* rssn_in, __global const double* rssn_params, __global double* rssn_out, __global const uint* rssn_cols, const uint rssn_rows, const uint rssn_inputs) {
    const uint row = (uint)get_global_id(0);
    if (row >= rssn_rows) { return; }
    const int irow = (int)row;
    const double v0 = rssn_params[0u];
    const int v1 = rssn_toidx(v0);
    const int v2 = ((int)(rssn_cols[3u]));
    const int v3 = rssn_idiv(v2, v1);
    const int v4 = irow;
    const int v5 = rssn_idiv(v4, v1);
    const int v6 = rssn_imod(v4, v1);
    const int v7 = v5 * v3;
    double a0 = 0.0;
    for (int j = 0; j < v3; j = j + 1) {
        const int v8 = j;
        const int v9 = v7 + v8;
        const int v10 = v8 * v1;
        const int v11 = v10 + v6;
        const double v12 = RSSN_LOAD(0u, v9);
        const double v13 = RSSN_LOAD(1u, v11);
        const double v14 = v12 * v13;
        a0 = a0 + v14;
    }
    const double v15 = a0;
    rssn_out[0u * rssn_rows + row] = v15;
}
