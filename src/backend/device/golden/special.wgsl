// rssn kernel: 1 inputs, 0 params, 2 outputs, 15 ops (f32)
struct Dims { rows: u32, inputs: u32, outputs: u32, pad: u32 }
@group(0) @binding(0) var<storage, read> rssn_in: array<f32>;
@group(0) @binding(1) var<storage, read> rssn_params: array<f32>;
@group(0) @binding(2) var<storage, read_write> rssn_out: array<f32>;
@group(0) @binding(3) var<uniform> rssn_dims: Dims;
@group(0) @binding(4) var<storage, read> rssn_cols: array<u32>;

fn rssn_powi(x: f32, n: i32) -> f32 {
    var e: u32 = select(u32(n), u32(-n), n < 0);
    var b: f32 = x;
    var acc: f32 = 1.0;
    loop {
        if (e == 0u) { break; }
        if ((e & 1u) == 1u) { acc = acc * b; }
        e = e >> 1u;
        if (e > 0u) { b = b * b; }
    }
    return select(acc, 1.0 / acc, n < 0);
}

fn rssn_powf(x: f32, y: f32) -> f32 {
    if (y == 0.0) { return 1.0; }
    if (x > 0.0) { return pow(x, y); }
    if (x == 0.0) {
        if (y > 0.0) { return 0.0; }
        return bitcast<f32>(0x7f800000u);
    }
    if (y == floor(y)) {
        let m: f32 = pow(-x, y);
        return select(m, -m, (y % 2.0) != 0.0);
    }
    return bitcast<f32>(0x7fc00000u);
}

fn rssn_cbrt(x: f32) -> f32 {
    if (x == 0.0) { return x; }
    let m: f32 = pow(abs(x), 0.33333334);
    return select(m, -m, x < 0.0);
}

fn rssn_load(col: u32, i: i32) -> f32 {
    let len = rssn_cols[2u * col + 1u];
    if (i < 0 || u32(i) >= len) { return bitcast<f32>(0x7fc00000u); }
    return rssn_in[rssn_cols[2u * col] + u32(i)];
}

fn rssn_gather(base: u32, off: f32, i: i32) -> f32 {
    if (!(off >= 0.0) || !(off < f32(rssn_dims.inputs - base))) { return bitcast<f32>(0x7fc00000u); }
    return rssn_load(base + u32(off), i);
}

@compute @workgroup_size(128)
fn kernel_main(@builtin(global_invocation_id) gid: vec3<u32>, @builtin(num_workgroups) nwg: vec3<u32>) {
    let row = gid.y * (nwg.x * 128u) + gid.x;
    let rssn_rows = rssn_dims.rows;
    if (row >= rssn_rows) { return; }
    let irow = i32(row);
    let v0: f32 = rssn_load(0u, irow);
    let v1: f32 = bitcast<f32>(0x7fc00000u);
    let v2: f32 = -2.5e-8;
    let v3: f32 = bitcast<f32>(0x7f800000u);
    let v4: f32 = rssn_powi(v0, -3i);
    let v5: f32 = rssn_powf(v0, v2);
    let v6: f32 = rssn_cbrt(v0);
    let v7: f32 = log(v0);
    let v8: f32 = select(0.0, 1.0, v0 >= v2);
    let v9: f32 = select(v5, v4, (v8 != 0.0) || (v8 != v8));
    let v10: f32 = rssn_gather(0u, v3, irow);
    let v11: f32 = atan2(v6, v7);
    let v12: f32 = v9 + v10;
    let v13: f32 = v12 + v11;
    let v14: f32 = v13 + v1;
    rssn_out[0u * rssn_rows + row] = v14;
    rssn_out[1u * rssn_rows + row] = v6;
}
