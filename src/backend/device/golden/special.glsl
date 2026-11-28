// rssn kernel: 1 inputs, 0 params, 2 outputs, 15 ops (float)
#version 450
layout(local_size_x = 64) in;
layout(std430, set = 0, binding = 0) readonly buffer RssnIn { float rssn_in[]; };
layout(std430, set = 0, binding = 1) readonly buffer RssnParams { float rssn_params[]; };
layout(std430, set = 0, binding = 2) buffer RssnOut { float rssn_out[]; };
layout(std140, set = 0, binding = 3) uniform RssnDims { uint rssn_rows; uint rssn_inputs; uint rssn_outputs; uint rssn_pad; };
layout(std430, set = 0, binding = 4) readonly buffer RssnCols { uint rssn_cols[]; };

float rssn_powi(float x, int n) {
    uint e = ((n < 0) ? uint(-n) : uint(n));
    float b = x;
    float acc = 1.0;
    while (e != 0u) {
        if ((e & 1u) == 1u) { acc = acc * b; }
        e = e >> 1u;
        if (e != 0u) { b = b * b; }
    }
    return ((n < 0) ? 1.0 / acc : acc);
}

float rssn_powf(float x, float y) {
    if (y == 0.0) { return 1.0; }
    if (x > 0.0) { return pow(x, y); }
    if (x == 0.0) {
        if (y > 0.0) { return 0.0; }
        return float(uintBitsToFloat(0x7f800000u));
    }
    if (y == floor(y)) {
        float m = pow(-x, y);
        return ((mod(y, 2.0) != 0.0) ? -m : m);
    }
    return float(uintBitsToFloat(0x7fc00000u));
}

float rssn_cbrt(float x) {
    if (x == 0.0) { return x; }
    float m = pow(abs(x), 0.33333334);
    return ((x < 0.0) ? -m : m);
}

float rssn_load(uint col, int i) {
    if (i < 0 || uint(i) >= rssn_cols[2u * col + 1u]) { return float(uintBitsToFloat(0x7fc00000u)); }
    return rssn_in[rssn_cols[2u * col] + uint(i)];
}

float rssn_gather(uint base, float off, int i) {
    if (!(off >= 0.0) || !(off < float(rssn_inputs - base))) { return float(uintBitsToFloat(0x7fc00000u)); }
    return rssn_load(base + uint(off), i);
}

void main() {
    uint row = gl_GlobalInvocationID.y * (gl_NumWorkGroups.x * 64u) + gl_GlobalInvocationID.x;
    if (row >= rssn_rows) { return; }
    int irow = int(row);
    float v0 = rssn_load(0u, irow);
    float v1 = float(uintBitsToFloat(0x7fc00000u));
    float v2 = -2.5e-8;
    float v3 = float(uintBitsToFloat(0x7f800000u));
    float v4 = rssn_powi(v0, -3);
    float v5 = rssn_powf(v0, v2);
    float v6 = rssn_cbrt(v0);
    float v7 = log(v0);
    float v8 = ((v0 >= v2) ? 1.0 : 0.0);
    float v9 = (((v8 != 0.0) || (v8 != v8)) ? v4 : v5);
    float v10 = rssn_gather(0u, v3, irow);
    float v11 = atan(v6, v7);
    float v12 = v9 + v10;
    float v13 = v12 + v11;
    float v14 = v13 + v1;
    rssn_out[0u * rssn_rows + row] = v14;
    rssn_out[1u * rssn_rows + row] = v6;
}
