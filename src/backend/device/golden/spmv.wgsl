// rssn kernel: 4 inputs, 0 params, 1 outputs, 16 ops, fold (f32)
struct Dims { rows: u32, inputs: u32, outputs: u32, pad: u32 }
@group(0) @binding(0) var<storage, read> rssn_in: array<f32>;
@group(0) @binding(1) var<storage, read> rssn_params: array<f32>;
@group(0) @binding(2) var<storage, read_write> rssn_out: array<f32>;
@group(0) @binding(3) var<uniform> rssn_dims: Dims;
@group(0) @binding(4) var<storage, read> rssn_cols: array<u32>;

fn rssn_toidx(x: f32) -> i32 {
    if (x != x) { return 0; }
    if (x >= 2147483500.0) { return 2147483647; }
    if (x <= -2147483500.0) { return -2147483647; }
    return i32(x);
}

fn rssn_load(col: u32, i: i32) -> f32 {
    let len = rssn_cols[2u * col + 1u];
    if (i < 0 || u32(i) >= len) { return bitcast<f32>(0x7fc00000u); }
    return rssn_in[rssn_cols[2u * col] + u32(i)];
}

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>, @builtin(num_workgroups) nwg: vec3<u32>) {
    let row = gid.y * (nwg.x * 64u) + gid.x;
    let rssn_rows = rssn_dims.rows;
    if (row >= rssn_rows) { return; }
    let irow = i32(row);
    let v0: i32 = irow;
    let v1: i32 = 1i;
    let v2: i32 = v0 + v1;
    let v3: f32 = rssn_load(0u, v0);
    let v4: f32 = rssn_load(0u, v2);
    let v5: i32 = rssn_toidx(v3);
    let v6: i32 = rssn_toidx(v4);
    let v7: i32 = v6 - v5;
    var a0: f32 = 0.0;
    for (var j: i32 = 0; j < v7; j = j + 1) {
        let v8: i32 = j;
        let v9: i32 = v5 + v8;
        let v10: f32 = rssn_load(2u, v9);
        let v11: f32 = rssn_load(1u, v9);
        let v12: i32 = rssn_toidx(v11);
        let v13: f32 = rssn_load(3u, v12);
        let v14: f32 = v10 * v13;
        a0 = a0 + v14;
    }
    let v15: f32 = a0;
    rssn_out[0u * rssn_rows + row] = v15;
}
