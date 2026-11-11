//! WGSL source from a [`KernelIr`]: pure string generation, no GPU needed.
//!
//! The emitted module has this interface (the host side, `gpu.rs` or any
//! other WebGPU embedder, binds exactly this):
//!
//! | binding | WGSL | contents |
//! |---|---|---|
//! | `@group(0) @binding(0)` | `var<storage, read> rssn_in: array<f32>` | input columns, column-major: `rssn_in[i * rows + row]` |
//! | `@group(0) @binding(1)` | `var<storage, read> rssn_params: array<f32>` | the kernel's parameters |
//! | `@group(0) @binding(2)` | `var<storage, read_write> rssn_out: array<f32>` | output columns, column-major |
//! | `@group(0) @binding(3)` | `var<uniform> rssn_dims: Dims` | `rows`, `inputs`, `outputs`, padding |
//!
//! One invocation computes one row. The dispatch is two-dimensional so that
//! more than `65535 * workgroup_size` rows fit:
//! `row = gid.y * (num_workgroups.x * workgroup_size) + gid.x`; the host
//! dispatches `ceil(rows / workgroup_size)` workgroups spread over x and y.
//!
//! ## Precision and semantics (deliberate differences from the CPU)
//!
//! * WGSL has no `f64`; everything is `f32`. Constants are rounded to `f32`.
//!   A device using this emitter reports `Dtype::F32` only.
//! * WGSL does not promise IEEE NaN/infinity propagation. Constants that are
//!   NaN or infinite are materialised with `bitcast`, and the helpers return
//!   NaN where the CPU does, but an implementation is allowed to optimise
//!   assuming no NaNs. Do not rely on NaN-valued results from a GPU.
//! * `pow` of a negative base is undefined in WGSL; `rssn_powf` handles
//!   integral exponents and returns NaN otherwise, as `f64::powf` does.
//! * Transcendental functions are implementation-defined to a few ulp of
//!   `f32`; they are not the correctly rounded `libm` results.

use std::fmt::Write as _;

use super::DeviceError;
use super::KernelIr;
use super::KernelOp;
use crate::backend::Intrinsic;
use crate::backend::tape::Cmp;

/// A source emitter for one target language. Implement it for PTX, SPIR-V
/// assembly, MSL, OpenCL C, ... (`docs/BACKENDS.md`).
pub trait KernelEmitter {
    /// Target name (`wgsl`, `ptx`, ...).
    fn target(&self) -> &'static str;

    /// Emits the program.
    ///
    /// # Errors
    /// [`DeviceError::InvalidKernel`] for a malformed program,
    /// [`DeviceError::Unsupported`] for something the target cannot express.
    fn emit(
        &self,
        ir: &KernelIr,
    ) -> Result<String, DeviceError>;
}

/// Options of the WGSL emitter.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct WgslOptions {
    /// `@workgroup_size(n)`.
    pub workgroup_size: u32,
    /// Entry point name.
    pub entry: String,
}

impl Default for WgslOptions {
    fn default() -> Self {
        Self {
            workgroup_size: 64,
            entry: "main".into(),
        }
    }
}

/// The WGSL emitter.
#[derive(Clone, Debug, Default)]
pub struct WgslEmitter {
    /// Options.
    pub options: WgslOptions,
}

impl KernelEmitter for WgslEmitter {
    fn target(&self) -> &'static str {
        "wgsl"
    }

    fn emit(
        &self,
        ir: &KernelIr,
    ) -> Result<String, DeviceError> {
        emit_wgsl(ir, &self.options)
    }
}

const NAN_EXPR: &str = "bitcast<f32>(0x7fc00000u)";
const POS_INF_EXPR: &str = "bitcast<f32>(0x7f800000u)";
const NEG_INF_EXPR: &str = "bitcast<f32>(0xff800000u)";

/// A WGSL `f32` expression for a constant.
#[must_use]
pub fn f32_literal(c: f64) -> String {
    #[allow(clippy::cast_possible_truncation)]
    let x = c as f32;
    if x.is_nan() {
        NAN_EXPR.into()
    } else if x == f32::INFINITY {
        POS_INF_EXPR.into()
    } else if x == f32::NEG_INFINITY {
        NEG_INF_EXPR.into()
    } else {
        format!("{x:?}")
    }
}

const fn intrinsic_name(f: Intrinsic) -> &'static str {
    match f {
        | Intrinsic::Exp => "exp",
        | Intrinsic::Ln => "log",
        | Intrinsic::Sin => "sin",
        | Intrinsic::Cos => "cos",
        | Intrinsic::Tan => "tan",
        | Intrinsic::Asin => "asin",
        | Intrinsic::Acos => "acos",
        | Intrinsic::Atan => "atan",
        | Intrinsic::Sinh => "sinh",
        | Intrinsic::Cosh => "cosh",
        | Intrinsic::Tanh => "tanh",
        | Intrinsic::Asinh => "asinh",
        | Intrinsic::Acosh => "acosh",
        | Intrinsic::Atanh => "atanh",
        | Intrinsic::Sqrt => "sqrt",
        | Intrinsic::Abs => "abs",
        | Intrinsic::Floor => "floor",
        | Intrinsic::Ceil => "ceil",
        | Intrinsic::Cbrt => "rssn_cbrt",
    }
}

const fn cmp_symbol(c: Cmp) -> &'static str {
    match c {
        | Cmp::Lt => "<",
        | Cmp::Le => "<=",
        | Cmp::Gt => ">",
        | Cmp::Ge => ">=",
        | Cmp::Eq => "==",
        | Cmp::Ne => "!=",
    }
}

const HELPER_POWI: &str = "\
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
";

const HELPER_POWF: &str = "\
fn rssn_powf(x: f32, y: f32) -> f32 {
    if (y == 0.0) { return 1.0; }
    if (x > 0.0) { return pow(x, y); }
    if (x == 0.0) {
        if (y > 0.0) { return 0.0; }
        return bitcast<f32>(0x7f800000u);
    }
    if (y == floor(y)) {
        let m = pow(-x, y);
        return select(m, -m, (y % 2.0) != 0.0);
    }
    return bitcast<f32>(0x7fc00000u);
}
";

const HELPER_CBRT: &str = "\
fn rssn_cbrt(x: f32) -> f32 {
    if (x == 0.0) { return x; }
    return sign(x) * pow(abs(x), 0.33333334);
}
";

const HELPER_GATHER: &str = "\
fn rssn_gather(base: u32, off: f32, row: u32) -> f32 {
    if (!(off >= 0.0) || !(off < f32(rssn_dims.inputs - base))) { return bitcast<f32>(0x7fc00000u); }
    return rssn_in[(base + u32(off)) * rssn_dims.rows + row];
}
";

/// Emits the WGSL module.
///
/// # Errors
/// [`DeviceError::InvalidKernel`] when `ir` does not validate or the
/// workgroup size is zero.
pub fn emit_wgsl(
    ir: &KernelIr,
    options: &WgslOptions,
) -> Result<String, DeviceError> {
    ir.validate()?;
    if options.workgroup_size == 0 {
        return Err(DeviceError::InvalidKernel(
            "workgroup size must be positive".into(),
        ));
    }
    let uses = |pred: &dyn Fn(&KernelOp) -> bool| ir.ops.iter().any(pred);
    let mut s = String::new();
    let _ = writeln!(
        s,
        "// rssn kernel: {} inputs, {} params, {} outputs, {} ops (f32)",
        ir.inputs,
        ir.params,
        ir.outputs.len(),
        ir.ops.len()
    );
    s.push_str("struct Dims { rows: u32, inputs: u32, outputs: u32, pad: u32 }\n");
    s.push_str("@group(0) @binding(0) var<storage, read> rssn_in: array<f32>;\n");
    s.push_str("@group(0) @binding(1) var<storage, read> rssn_params: array<f32>;\n");
    s.push_str("@group(0) @binding(2) var<storage, read_write> rssn_out: array<f32>;\n");
    s.push_str("@group(0) @binding(3) var<uniform> rssn_dims: Dims;\n\n");
    if uses(&|o| matches!(o, KernelOp::PowI(..))) {
        s.push_str(HELPER_POWI);
        s.push('\n');
    }
    if uses(&|o| matches!(o, KernelOp::PowF(..))) {
        s.push_str(HELPER_POWF);
        s.push('\n');
    }
    if uses(&|o| matches!(o, KernelOp::Unary(Intrinsic::Cbrt, _))) {
        s.push_str(HELPER_CBRT);
        s.push('\n');
    }
    if uses(&|o| matches!(o, KernelOp::Gather(..))) {
        s.push_str(HELPER_GATHER);
        s.push('\n');
    }
    let w = options.workgroup_size;
    let _ = writeln!(s, "@compute @workgroup_size({w})");
    let _ = writeln!(
        s,
        "fn {}(@builtin(global_invocation_id) gid: vec3<u32>, @builtin(num_workgroups) nwg: \
         vec3<u32>) {{",
        options.entry
    );
    let _ = writeln!(s, "    let row = gid.y * (nwg.x * {w}u) + gid.x;");
    s.push_str("    if (row >= rssn_dims.rows) { return; }\n");
    for (d, op) in ir.ops.iter().enumerate() {
        let expr = match *op {
            | KernelOp::Const(c) => f32_literal(c),
            | KernelOp::Input(i) => format!("rssn_in[{i}u * rssn_dims.rows + row]"),
            | KernelOp::Param(i) => format!("rssn_params[{i}u]"),
            | KernelOp::Add(a, b) => format!("v{a} + v{b}"),
            | KernelOp::Sub(a, b) => format!("v{a} - v{b}"),
            | KernelOp::Mul(a, b) => format!("v{a} * v{b}"),
            | KernelOp::Div(a, b) => format!("v{a} / v{b}"),
            | KernelOp::Neg(a) => format!("-v{a}"),
            | KernelOp::PowI(a, n) => format!("rssn_powi(v{a}, {n}i)"),
            | KernelOp::PowF(a, b) => format!("rssn_powf(v{a}, v{b})"),
            | KernelOp::Unary(f, a) => format!("{}(v{a})", intrinsic_name(f)),
            | KernelOp::Atan2(a, b) => format!("atan2(v{a}, v{b})"),
            | KernelOp::Cmp(c, a, b) => {
                format!("select(0.0, 1.0, v{a} {} v{b})", cmp_symbol(c))
            },
            | KernelOp::Select(c, a, b) => {
                format!("select(v{b}, v{a}, (v{c} != 0.0) || (v{c} != v{c}))")
            },
            | KernelOp::Gather(base, i) => format!("rssn_gather({base}u, v{i}, row)"),
        };
        let _ = writeln!(s, "    let v{d}: f32 = {expr};");
    }
    for (k, &r) in ir.outputs.iter().enumerate() {
        let _ = writeln!(s, "    rssn_out[{k}u * rssn_dims.rows + row] = v{r};");
    }
    s.push_str("}\n");
    Ok(s)
}
