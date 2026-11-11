// rssn kernel: 2 inputs, 1 params, 1 outputs, 5 ops (f32)
struct Dims { rows: u32, inputs: u32, outputs: u32, pad: u32 }
@group(0) @binding(0) var<storage, read> rssn_in: array<f32>;
@group(0) @binding(1) var<storage, read> rssn_params: array<f32>;
@group(0) @binding(2) var<storage, read_write> rssn_out: array<f32>;
@group(0) @binding(3) var<uniform> rssn_dims: Dims;

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>, @builtin(num_workgroups) nwg: vec3<u32>) {
    let row = gid.y * (nwg.x * 64u) + gid.x;
    if (row >= rssn_dims.rows) { return; }
    let v0: f32 = rssn_in[0u * rssn_dims.rows + row];
    let v1: f32 = rssn_in[1u * rssn_dims.rows + row];
    let v2: f32 = rssn_params[0u];
    let v3: f32 = v2 * v0;
    let v4: f32 = v3 + v1;
    rssn_out[0u * rssn_dims.rows + row] = v4;
}
