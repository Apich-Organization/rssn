// rssn kernel: 1 inputs, 0 params, 2 outputs, 15 ops (float)
static inline float rssn_powi(float x, int n) {
    uint e = ((n < 0) ? ((uint)(-n)) : ((uint)(n)));
    float b = x;
    float acc = 1.0f;
    while (e != 0u) {
        if ((e & 1u) == 1u) { acc = acc * b; }
        e = e >> 1u;
        if (e != 0u) { b = b * b; }
    }
    return ((n < 0) ? 1.0f / acc : acc);
}

#define RSSN_LOAD(col, i) (((i) >= 0 && ((uint)(i)) < rssn_cols[2u * (col) + 1u]) ? rssn_in[rssn_cols[2u * (col)] + ((uint)(i))] : ((float)(NAN)))
#define RSSN_GATHER(base, off, i) ((((off) >= 0.0f) && ((off) < ((float)(rssn_inputs - (base))))) ? RSSN_LOAD((base) + ((uint)(off)), i) : ((float)(NAN)))

__kernel void rssn_main(__global const float* rssn_in, __global const float* rssn_params, __global float* rssn_out, __global const uint* rssn_cols, const uint rssn_rows, const uint rssn_inputs) {
    const uint row = (uint)get_global_id(0);
    if (row >= rssn_rows) { return; }
    const int irow = (int)row;
    const float v0 = RSSN_LOAD(0u, irow);
    const float v1 = ((float)(NAN));
    const float v2 = -2.5e-8f;
    const float v3 = ((float)(INFINITY));
    const float v4 = rssn_powi(v0, -3);
    const float v5 = pow(v0, v2);
    const float v6 = cbrt(v0);
    const float v7 = log(v0);
    const float v8 = ((v0 >= v2) ? 1.0f : 0.0f);
    const float v9 = (((v8 != 0.0f) || (v8 != v8)) ? v4 : v5);
    const float v10 = RSSN_GATHER(0u, v3, irow);
    const float v11 = atan2(v6, v7);
    const float v12 = v9 + v10;
    const float v13 = v12 + v11;
    const float v14 = v13 + v1;
    rssn_out[0u * rssn_rows + row] = v14;
    rssn_out[1u * rssn_rows + row] = v6;
}
