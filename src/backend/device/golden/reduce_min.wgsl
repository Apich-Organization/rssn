// rssn kernel: 1 inputs, 0 params, 1 outputs, 11 ops, fold (f32)
struct Dims { rows: u32, inputs: u32, outputs: u32, pad: u32 }
@group(0) @binding(0) var<storage, read> rssn_in: array<f32>;
@group(0) @binding(1) var<storage, read> rssn_params: array<f32>;
@group(0) @binding(2) var<storage, read_write> rssn_out: array<f32>;
@group(0) @binding(3) var<uniform> rssn_dims: Dims;
@group(0) @binding(4) var<storage, read> rssn_cols: array<u32>;

fn rssn_min(a: f32, b: f32) -> f32 {
    if (a != a || b != b) { return bitcast<f32>(0x7fc00000u); }
    return select(b, a, a < b);
}

fn rssn_max(a: f32, b: f32) -> f32 {
    if (a != a || b != b) { return bitcast<f32>(0x7fc00000u); }
    return select(b, a, a > b);
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
    let v0: i32 = 256i;
    let v1: i32 = irow;
    let v2: i32 = v1 * v0;
    let v5: i32 = i32(rssn_cols[1u]);
    let v8: f32 = bitcast<f32>(0x7f800000u);
    var a0: f32 = bitcast<f32>(0x7f800000u);
    for (var j: i32 = 0; j < v0; j = j + 1) {
        let v3: i32 = j;
        let v4: i32 = v2 + v3;
        let v6: f32 = select(0.0, 1.0, v4 < v5);
        let v7: f32 = rssn_load(0u, v4);
        let v9: f32 = select(v8, v7, (v6 != 0.0) || (v6 != v6));
        a0 = rssn_min(a0, v9);
    }
    let v10: f32 = a0;
    rssn_out[0u * rssn_rows + row] = v10;
}
