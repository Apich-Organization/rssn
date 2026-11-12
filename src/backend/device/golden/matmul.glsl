// rssn kernel: 2 inputs, 1 params, 1 outputs, 16 ops, fold (float)
#version 450
layout(local_size_x = 64) in;
layout(std430, set = 0, binding = 0) readonly buffer RssnIn { float rssn_in[]; };
layout(std430, set = 0, binding = 1) readonly buffer RssnParams { float rssn_params[]; };
layout(std430, set = 0, binding = 2) buffer RssnOut { float rssn_out[]; };
layout(std140, set = 0, binding = 3) uniform RssnDims { uint rssn_rows; uint rssn_inputs; uint rssn_outputs; uint rssn_pad; };
layout(std430, set = 0, binding = 4) readonly buffer RssnCols { uint rssn_cols[]; };

int rssn_toidx(float x) {
    if (x != x) { return 0; }
    if (x >= 2147483500.0) { return 2147483647; }
    if (x <= -2147483500.0) { return -2147483647; }
    return int(x);
}

int rssn_idiv(int a, int b) {
    if (b == 0) { return 0; }
    return a / b;
}

int rssn_imod(int a, int b) {
    if (b == 0) { return 0; }
    return a % b;
}

float rssn_load(uint col, int i) {
    if (i < 0 || uint(i) >= rssn_cols[2u * col + 1u]) { return float(uintBitsToFloat(0x7fc00000u)); }
    return rssn_in[rssn_cols[2u * col] + uint(i)];
}

void main() {
    uint row = gl_GlobalInvocationID.y * (gl_NumWorkGroups.x * 64u) + gl_GlobalInvocationID.x;
    if (row >= rssn_rows) { return; }
    int irow = int(row);
    float v0 = rssn_params[0u];
    int v1 = rssn_toidx(v0);
    int v2 = int(rssn_cols[3u]);
    int v3 = rssn_idiv(v2, v1);
    int v4 = irow;
    int v5 = rssn_idiv(v4, v1);
    int v6 = rssn_imod(v4, v1);
    int v7 = v5 * v3;
    float a0 = 0.0;
    for (int j = 0; j < v3; j = j + 1) {
        int v8 = j;
        int v9 = v7 + v8;
        int v10 = v8 * v1;
        int v11 = v10 + v6;
        float v12 = rssn_load(0u, v9);
        float v13 = rssn_load(1u, v11);
        float v14 = v12 * v13;
        a0 = a0 + v14;
    }
    float v15 = a0;
    rssn_out[0u * rssn_rows + row] = v15;
}
