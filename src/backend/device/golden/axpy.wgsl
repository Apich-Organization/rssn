// rssn kernel: 2 inputs, 1 params, 1 outputs, 5 ops (f32)
struct Dims { rows: u32, inputs: u32, outputs: u32, pad: u32 }
@group(0) @binding(0) var<storage, read> rssn_in: array<f32>;
@group(0) @binding(1) var<storage, read> rssn_params: array<f32>;
@group(0) @binding(2) var<storage, read_write> rssn_out: array<f32>;
@group(0) @binding(3) var<uniform> rssn_dims: Dims;
@group(0) @binding(4) var<storage, read> rssn_cols: array<u32>;

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
    let v0: f32 = rssn_load(0u, irow);
    let v1: f32 = rssn_load(1u, irow);
    let v2: f32 = rssn_params[0u];
    let v3: f32 = v2 * v0;
    let v4: f32 = v3 + v1;
    rssn_out[0u * rssn_rows + row] = v4;
}
