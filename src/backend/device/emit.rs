//! Source emitters: one walk of a [`KernelIr`], five target languages.
//!
//! The first emitter was a WGSL printer. Adding a second target by copying it
//! would have duplicated the part that is hard to get right (phases, the
//! fold, bounds-checked loads, NaN conventions) and kept only the part that
//! is trivial (spelling). This module inverts that: [`emit`] walks the
//! program once and asks a [`Target`] how to *spell* each construct.
//!
//! | target | language | precisions | status |
//! |---|---|---|---|
//! | [`Target::Wgsl`] | WebGPU Shading Language | `f32` | executed by `gpu-wgpu`, validated by naga |
//! | [`Target::Glsl`] | GLSL 4.50 compute (Vulkan / OpenGL) | `f32`, `f64` | emitter, golden-tested |
//! | [`Target::OpenCl`] | OpenCL C 1.2 | `f32`, `f64` | emitter, golden-tested |
//! | [`Target::Cuda`] | CUDA C++ | `f32`, `f64` | emitter, golden-tested |
//! | [`Target::Metal`] | Metal Shading Language | `f32` | emitter, golden-tested |
//!
//! "Emitter" means exactly that: the text is generated and pinned by golden
//! files, but this crate ships no runtime that compiles and launches it. A
//! device for one of those APIs supplies the launch half and reuses this
//! half unchanged (`docs/BACKENDS.md`).
//!
//! # The interface every target shares
//!
//! | name | contents |
//! |---|---|
//! | `rssn_in` | every input column, concatenated |
//! | `rssn_params` | the kernel's parameters |
//! | `rssn_out` | output columns, column-major: `rssn_out[k * rows + row]` |
//! | `rssn_cols` | two unsigned integers per input column: its offset in `rssn_in` and its length |
//! | `rssn_rows`, `rssn_inputs` | the number of rows of the launch and of input columns |
//!
//! WGSL and GLSL bind them as `@group(0)`/`set = 0` bindings 0 (`rssn_in`),
//! 1 (`rssn_params`), 2 (`rssn_out`), 3 (a uniform `rows, inputs, outputs,
//! pad`) and 4 (`rssn_cols`); OpenCL and CUDA take them as kernel arguments in
//! the order of the table; Metal uses buffer indices 0 to 4 like the
//! bindings. One invocation computes one row.
//!
//! # Precision and semantics (deliberate differences from the CPU)
//!
//! * `f32` targets round constants to `f32`; a device using them reports
//!   [`Dtype::F32`](super::Dtype::F32) only. Indices are 32-bit.
//! * Shading languages do not promise IEEE NaN/infinity propagation. NaN and
//!   infinite constants are materialised from bit patterns and the helpers
//!   return NaN where the CPU does, but an implementation may optimise
//!   assuming no NaNs. Do not rely on NaN-valued results from a GPU.
//! * `pow` of a negative base is undefined in WGSL and GLSL; `rssn_powf`
//!   handles integral exponents and returns NaN otherwise, as `f64::powf`.
//! * Transcendental functions are implementation-defined to a few ulp.
//! * Index data that a kernel reads from a column (`ToIdx(Load(..))`, as in
//!   a CSR structure) is stored as a float: it is exact below `2^24` on an
//!   `f32` target and below `2^53` on an `f64` one.

use std::fmt::Write as _;

use serde::Deserialize;
use serde::Serialize;

use super::DeviceError;
use super::KernelIr;
use super::KernelOp;
use super::ReduceOp;
use super::kernel::Analysis;
use super::kernel::Ty;
use crate::backend::Intrinsic;
use crate::backend::tape::Cmp;
use crate::backend::tape::Reg;

/// A target language.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Target {
    /// WebGPU Shading Language.
    Wgsl,
    /// GLSL 4.50 compute shader.
    Glsl,
    /// OpenCL C.
    OpenCl,
    /// CUDA C++.
    Cuda,
    /// Metal Shading Language.
    Metal,
}

impl Target {
    /// Every target.
    pub const ALL: [Self; 5] = [
        Self::Wgsl,
        Self::Glsl,
        Self::OpenCl,
        Self::Cuda,
        Self::Metal,
    ];

    /// A short name, also the extension of the golden files.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            | Self::Wgsl => "wgsl",
            | Self::Glsl => "glsl",
            | Self::OpenCl => "cl",
            | Self::Cuda => "cu",
            | Self::Metal => "metal",
        }
    }

    /// Whether the language has a 64-bit float.
    #[must_use]
    pub const fn has_f64(self) -> bool {
        matches!(self, Self::Glsl | Self::OpenCl | Self::Cuda)
    }
}

/// The floating-point type of a kernel's numbers.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Precision {
    /// 32-bit.
    F32,
    /// 64-bit.
    F64,
}

/// Options of [`emit`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct EmitOptions {
    /// Target language.
    pub target: Target,
    /// Float type.
    pub precision: Precision,
    /// Threads per workgroup (WGSL, GLSL); ignored elsewhere.
    pub workgroup_size: u32,
    /// Entry point name.
    pub entry: String,
}

impl EmitOptions {
    /// The defaults for `target`: `f32`, workgroups of 64, entry `main`
    /// (`rssn_main` where `main` is reserved).
    #[must_use]
    pub fn new(target: Target) -> Self {
        Self {
            target,
            precision: Precision::F32,
            workgroup_size: 64,
            entry: match target {
                | Target::Wgsl | Target::Glsl => "main".into(),
                | _ => "rssn_main".into(),
            },
        }
    }

    /// The same options with another precision.
    #[must_use]
    pub const fn precision(
        mut self,
        precision: Precision,
    ) -> Self {
        self.precision = precision;
        self
    }
}

/// How one target spells things.
#[derive(Clone, Copy)]
struct Lang {
    target: Target,
    precision: Precision,
}

impl Lang {
    const fn wgsl(self) -> bool {
        matches!(self.target, Target::Wgsl)
    }

    const fn glsl(self) -> bool {
        matches!(self.target, Target::Glsl)
    }

    const fn real(self) -> &'static str {
        match (self.target, self.precision) {
            | (Target::Wgsl, _) => "f32",
            | (_, Precision::F32) => "float",
            | (_, Precision::F64) => "double",
        }
    }

    const fn int(self) -> &'static str {
        if self.wgsl() { "i32" } else { "int" }
    }

    const fn uint(self) -> &'static str {
        match self.target {
            | Target::Wgsl => "u32",
            | Target::Glsl | Target::OpenCl | Target::Metal => "uint",
            | Target::Cuda => "unsigned int",
        }
    }

    /// A cast of `expr` to the type named `ty`.
    fn cast(
        self,
        ty: &str,
        expr: &str,
    ) -> String {
        if self.wgsl() || self.glsl() {
            format!("{ty}({expr})")
        } else {
            format!("(({ty})({expr}))")
        }
    }

    fn nan(self) -> String {
        match self.target {
            | Target::Wgsl => "bitcast<f32>(0x7fc00000u)".into(),
            | Target::Glsl => self.cast(self.real(), "uintBitsToFloat(0x7fc00000u)"),
            | _ => self.cast(self.real(), "NAN"),
        }
    }

    fn inf(
        self,
        negative: bool,
    ) -> String {
        match (self.target, negative) {
            | (Target::Wgsl, false) => "bitcast<f32>(0x7f800000u)".into(),
            | (Target::Wgsl, true) => "bitcast<f32>(0xff800000u)".into(),
            | (Target::Glsl, false) => self.cast(self.real(), "uintBitsToFloat(0x7f800000u)"),
            | (Target::Glsl, true) => self.cast(self.real(), "uintBitsToFloat(0xff800000u)"),
            | (_, false) => self.cast(self.real(), "INFINITY"),
            | (_, true) => self.cast(self.real(), "-INFINITY"),
        }
    }

    /// A literal of the float type.
    fn lit(
        self,
        c: f64,
    ) -> String {
        #[allow(clippy::cast_possible_truncation)]
        let rounded = match self.precision {
            | Precision::F32 => f64::from(c as f32),
            | Precision::F64 => c,
        };
        if rounded.is_nan() {
            return self.nan();
        }
        if rounded.is_infinite() {
            return self.inf(rounded < 0.0);
        }
        #[allow(clippy::cast_possible_truncation)]
        let text = match self.precision {
            | Precision::F32 => format!("{:?}", c as f32),
            | Precision::F64 => format!("{c:?}"),
        };
        match (self.target, self.precision) {
            | (Target::Wgsl, _) | (Target::Glsl, Precision::F32) => text,
            | (Target::Glsl, Precision::F64) => format!("{text}lf"),
            | (_, Precision::F32) => format!("{text}f"),
            | (_, Precision::F64) => text,
        }
    }

    /// `if cond { a } else { b }` for a boolean expression `cond`.
    fn pick(
        self,
        cond: &str,
        a: &str,
        b: &str,
    ) -> String {
        if self.wgsl() {
            format!("select({b}, {a}, {cond})")
        } else {
            format!("(({cond}) ? {a} : {b})")
        }
    }

    /// A declaration of an immutable local.
    fn bind(
        self,
        ty: &str,
        name: &str,
        expr: &str,
    ) -> String {
        match self.target {
            | Target::Wgsl => format!("let {name}: {ty} = {expr};"),
            // GLSL `const` requires a constant expression.
            | Target::Glsl => format!("{ty} {name} = {expr};"),
            | _ => format!("const {ty} {name} = {expr};"),
        }
    }

    /// A declaration of a mutable local.
    fn var(
        self,
        ty: &str,
        name: &str,
        expr: &str,
    ) -> String {
        if self.wgsl() {
            format!("var {name}: {ty} = {expr};")
        } else {
            format!("{ty} {name} = {expr};")
        }
    }

    fn intrinsic(
        self,
        f: Intrinsic,
    ) -> &'static str {
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
            | Intrinsic::Abs => {
                if self.wgsl() || self.glsl() {
                    "abs"
                } else {
                    "fabs"
                }
            },
            | Intrinsic::Floor => "floor",
            | Intrinsic::Ceil => "ceil",
            | Intrinsic::Cbrt => {
                match self.target {
                    | Target::OpenCl | Target::Cuda => "cbrt",
                    | _ => "rssn_cbrt",
                }
            },
        }
    }

    /// Whether `pow` needs the `rssn_powf` wrapper (undefined for a negative
    /// base).
    const fn wraps_pow(self) -> bool {
        matches!(self.target, Target::Wgsl | Target::Glsl)
    }

    /// The qualifier in front of a helper function.
    const fn helper_prefix(self) -> &'static str {
        match self.target {
            | Target::Cuda => "__device__ ",
            | Target::Metal | Target::OpenCl => "static inline ",
            | Target::Wgsl | Target::Glsl => "",
        }
    }

    /// A helper function definition: `name(params) -> ret { body }`.
    fn function(
        self,
        ret: &str,
        name: &str,
        params: &[(&str, &str)],
        body: &str,
    ) -> String {
        let mut s = String::new();
        if self.wgsl() {
            let args: Vec<String> = params.iter().map(|(n, t)| format!("{n}: {t}")).collect();
            let _ = writeln!(s, "fn {name}({}) -> {ret} {{", args.join(", "));
        } else {
            let args: Vec<String> = params.iter().map(|(n, t)| format!("{t} {n}")).collect();
            let _ = writeln!(
                s,
                "{}{ret} {name}({}) {{",
                self.helper_prefix(),
                args.join(", ")
            );
        }
        s.push_str(body);
        s.push_str("}\n");
        s
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

/// What a program uses, which decides the helpers that are emitted.
#[derive(Default)]
struct Uses {
    powi: bool,
    powf: bool,
    cbrt: bool,
    load: bool,
    gather: bool,
    toidx: bool,
    idiv: bool,
    imod: bool,
    minmax: bool,
}

impl Uses {
    fn of(ir: &KernelIr) -> Self {
        let mut u = Self::default();
        for op in &ir.ops {
            match op {
                | KernelOp::PowI(..) => u.powi = true,
                | KernelOp::PowF(..) => u.powf = true,
                | KernelOp::Unary(Intrinsic::Cbrt, _) => u.cbrt = true,
                | KernelOp::Input(_) | KernelOp::Load(..) => u.load = true,
                | KernelOp::Gather(..) => {
                    u.gather = true;
                    u.load = true;
                },
                | KernelOp::ToIdx(_) => u.toidx = true,
                | KernelOp::IDiv(..) => u.idiv = true,
                | KernelOp::IMod(..) => u.imod = true,
                | _ => {},
            }
        }
        if let Some(fold) = &ir.fold {
            u.minmax = fold
                .accs
                .iter()
                .any(|(op, _)| matches!(op, ReduceOp::Min | ReduceOp::Max));
        }
        u
    }
}

/// Emits `ir` for the target in `options`.
///
/// # Errors
/// [`DeviceError::InvalidKernel`] when `ir` does not validate or the
/// workgroup size is zero; [`DeviceError::Unsupported`] when the target has
/// no type of the requested precision.
pub fn emit(
    ir: &KernelIr,
    options: &EmitOptions,
) -> Result<String, DeviceError> {
    let analysis = ir.analyse()?;
    if options.workgroup_size == 0 {
        return Err(DeviceError::InvalidKernel(
            "workgroup size must be positive".into(),
        ));
    }
    if options.precision == Precision::F64 && !options.target.has_f64() {
        return Err(DeviceError::Unsupported(format!(
            "{} has no 64-bit float",
            options.target.name()
        )));
    }
    let lang = Lang {
        target: options.target,
        precision: options.precision,
    };
    let uses = Uses::of(ir);
    let mut s = String::new();
    let _ = writeln!(
        s,
        "// rssn kernel: {} inputs, {} params, {} outputs, {} ops{} ({})",
        ir.inputs,
        ir.params,
        ir.outputs.len(),
        ir.ops.len(),
        if ir.fold.is_some() {
            ", fold"
        } else {
            ""
        },
        lang.real()
    );
    prelude(&mut s, lang, options);
    helpers(&mut s, lang, &uses);
    entry_open(&mut s, lang, options);
    body(&mut s, lang, ir, &analysis);
    s.push_str("}\n");
    Ok(s)
}

/// Declarations before the helpers: bindings (WGSL, GLSL), pragmas, macros.
fn prelude(
    s: &mut String,
    lang: Lang,
    options: &EmitOptions,
) {
    let real = lang.real();
    match lang.target {
        | Target::Wgsl => {
            s.push_str("struct Dims { rows: u32, inputs: u32, outputs: u32, pad: u32 }\n");
            s.push_str("@group(0) @binding(0) var<storage, read> rssn_in: array<f32>;\n");
            s.push_str("@group(0) @binding(1) var<storage, read> rssn_params: array<f32>;\n");
            s.push_str("@group(0) @binding(2) var<storage, read_write> rssn_out: array<f32>;\n");
            s.push_str("@group(0) @binding(3) var<uniform> rssn_dims: Dims;\n");
            s.push_str("@group(0) @binding(4) var<storage, read> rssn_cols: array<u32>;\n\n");
        },
        | Target::Glsl => {
            s.push_str("#version 450\n");
            let _ = writeln!(s, "layout(local_size_x = {}) in;", options.workgroup_size);
            let _ = writeln!(
                s,
                "layout(std430, set = 0, binding = 0) readonly buffer RssnIn {{ {real} rssn_in[]; \
                 }};"
            );
            let _ = writeln!(
                s,
                "layout(std430, set = 0, binding = 1) readonly buffer RssnParams {{ {real} \
                 rssn_params[]; }};"
            );
            let _ = writeln!(
                s,
                "layout(std430, set = 0, binding = 2) buffer RssnOut {{ {real} rssn_out[]; }};"
            );
            s.push_str(
                "layout(std140, set = 0, binding = 3) uniform RssnDims { uint rssn_rows; uint \
                 rssn_inputs; uint rssn_outputs; uint rssn_pad; };\n",
            );
            s.push_str(
                "layout(std430, set = 0, binding = 4) readonly buffer RssnCols { uint \
                 rssn_cols[]; };\n\n",
            );
        },
        | Target::OpenCl => {
            if lang.precision == Precision::F64 {
                s.push_str("#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n");
            }
        },
        | Target::Cuda => s.push_str("#include <math.h>\n"),
        | Target::Metal => s.push_str("#include <metal_stdlib>\nusing namespace metal;\n"),
    }
}

/// Helper functions and macros, each emitted only when the program uses it.
#[allow(clippy::too_many_lines)]
fn helpers(
    s: &mut String,
    lang: Lang,
    uses: &Uses,
) {
    let real = lang.real();
    let int = lang.int();
    let uint = lang.uint();
    let nan = lang.nan();
    let one = lang.lit(1.0);
    let zero = lang.lit(0.0);
    if uses.powi {
        let e0 = lang.pick("n < 0", &lang.cast(uint, "-n"), &lang.cast(uint, "n"));
        let mut b = String::new();
        let _ = writeln!(b, "    {}", lang.var(uint, "e", &e0));
        let _ = writeln!(b, "    {}", lang.var(real, "b", "x"));
        let _ = writeln!(b, "    {}", lang.var(real, "acc", &one));
        if lang.wgsl() {
            b.push_str(
                "    loop {\n        if (e == 0u) { break; }\n        if ((e & 1u) == 1u) { acc = \
                 acc * b; }\n        e = e >> 1u;\n        if (e > 0u) { b = b * b; }\n    }\n",
            );
        } else {
            b.push_str(
                "    while (e != 0u) {\n        if ((e & 1u) == 1u) { acc = acc * b; }\n        e \
                 = e >> 1u;\n        if (e != 0u) { b = b * b; }\n    }\n",
            );
        }
        let _ = writeln!(
            b,
            "    return {};",
            lang.pick("n < 0", &format!("{one} / acc"), "acc")
        );
        s.push_str(&lang.function(real, "rssn_powi", &[("x", real), ("n", int)], &b));
        s.push('\n');
    }
    if uses.powf && lang.wraps_pow() {
        let mut b = String::new();
        let _ = writeln!(b, "    if (y == {zero}) {{ return {one}; }}");
        let _ = writeln!(b, "    if (x > {zero}) {{ return pow(x, y); }}");
        let _ = writeln!(b, "    if (x == {zero}) {{");
        let _ = writeln!(b, "        if (y > {zero}) {{ return {zero}; }}");
        let _ = writeln!(b, "        return {};", lang.inf(false));
        b.push_str("    }\n");
        b.push_str("    if (y == floor(y)) {\n");
        let _ = writeln!(b, "        {}", lang.bind(real, "m", "pow(-x, y)"));
        let two = lang.lit(2.0);
        let odd = if lang.wgsl() {
            format!("(y % {two}) != {zero}")
        } else {
            format!("mod(y, {two}) != {zero}")
        };
        let _ = writeln!(b, "        return {};", lang.pick(&odd, "-m", "m"));
        b.push_str("    }\n");
        let _ = writeln!(b, "    return {nan};");
        s.push_str(&lang.function(real, "rssn_powf", &[("x", real), ("y", real)], &b));
        s.push('\n');
    }
    if uses.cbrt && lang.intrinsic(Intrinsic::Cbrt) == "rssn_cbrt" {
        let third = lang.lit(1.0 / 3.0);
        let abs = lang.intrinsic(Intrinsic::Abs);
        let mut b = String::new();
        let _ = writeln!(b, "    if (x == {zero}) {{ return x; }}");
        let _ = writeln!(
            b,
            "    {}",
            lang.bind(real, "m", &format!("pow({abs}(x), {third})"))
        );
        let _ = writeln!(
            b,
            "    return {};",
            lang.pick(&format!("x < {zero}"), "-m", "m")
        );
        s.push_str(&lang.function(real, "rssn_cbrt", &[("x", real)], &b));
        s.push('\n');
    }
    if uses.minmax {
        for (name, cmp) in [("rssn_min", "<"), ("rssn_max", ">")] {
            let mut b = String::new();
            let _ = writeln!(b, "    if (a != a || b != b) {{ return {nan}; }}");
            let _ = writeln!(
                b,
                "    return {};",
                lang.pick(&format!("a {cmp} b"), "a", "b")
            );
            s.push_str(&lang.function(real, name, &[("a", real), ("b", real)], &b));
            s.push('\n');
        }
    }
    if uses.toidx {
        let mut b = String::new();
        b.push_str("    if (x != x) { return 0; }\n");
        let hi = lang.lit(2_147_483_520.0);
        let lo = lang.lit(-2_147_483_520.0);
        let _ = writeln!(b, "    if (x >= {hi}) {{ return 2147483647; }}");
        let _ = writeln!(b, "    if (x <= {lo}) {{ return -2147483647; }}");
        let _ = writeln!(b, "    return {};", lang.cast(int, "x"));
        s.push_str(&lang.function(int, "rssn_toidx", &[("x", real)], &b));
        s.push('\n');
    }
    for (used, name, symbol) in [(uses.idiv, "rssn_idiv", "/"), (uses.imod, "rssn_imod", "%")] {
        if used {
            let b = format!("    if (b == 0) {{ return 0; }}\n    return a {symbol} b;\n");
            s.push_str(&lang.function(int, name, &[("a", int), ("b", int)], &b));
            s.push('\n');
        }
    }
    // Buffer accessors. WGSL and GLSL see the buffers as globals, so these
    // are functions; the C-like targets receive them as kernel arguments, so
    // they are macros expanded inside the entry point.
    if uses.load {
        if lang.wgsl() {
            s.push_str(
                "fn rssn_load(col: u32, i: i32) -> f32 {\n    let len = rssn_cols[2u * col + \
                 1u];\n    if (i < 0 || u32(i) >= len) { return bitcast<f32>(0x7fc00000u); }\n    \
                 return rssn_in[rssn_cols[2u * col] + u32(i)];\n}\n\n",
            );
        } else if lang.glsl() {
            let _ = writeln!(
                s,
                "{real} rssn_load(uint col, int i) {{\n    if (i < 0 || uint(i) >= rssn_cols[2u * \
                 col + 1u]) {{ return {nan}; }}\n    return rssn_in[rssn_cols[2u * col] + \
                 uint(i)];\n}}\n"
            );
        } else {
            let u = |e: &str| lang.cast(uint, e);
            let _ = writeln!(
                s,
                "#define RSSN_LOAD(col, i) (((i) >= 0 && {} < rssn_cols[2u * (col) + 1u]) ? \
                 rssn_in[rssn_cols[2u * (col)] + {}] : {nan})",
                u("i"),
                u("i")
            );
        }
    }
    if uses.gather {
        if lang.wgsl() {
            s.push_str(
                "fn rssn_gather(base: u32, off: f32, i: i32) -> f32 {\n    if (!(off >= 0.0) || \
                 !(off < f32(rssn_dims.inputs - base))) { return bitcast<f32>(0x7fc00000u); }\n    \
                 return rssn_load(base + u32(off), i);\n}\n\n",
            );
        } else if lang.glsl() {
            let _ = writeln!(
                s,
                "{real} rssn_gather(uint base, {real} off, int i) {{\n    if (!(off >= {zero}) || \
                 !(off < {real}(rssn_inputs - base))) {{ return {nan}; }}\n    return \
                 rssn_load(base + uint(off), i);\n}}\n"
            );
        } else {
            let _ = writeln!(
                s,
                "#define RSSN_GATHER(base, off, i) ((((off) >= {zero}) && ((off) < {})) ? \
                 RSSN_LOAD((base) + {}, i) : {nan})",
                lang.cast(real, "rssn_inputs - (base)"),
                lang.cast(uint, "off")
            );
        }
    }
    if !(lang.wgsl() || lang.glsl()) && (uses.load || uses.gather) {
        s.push('\n');
    }
}

/// The entry point up to and including the row computation and the range
/// guard.
fn entry_open(
    s: &mut String,
    lang: Lang,
    options: &EmitOptions,
) {
    let real = lang.real();
    let entry = &options.entry;
    let w = options.workgroup_size;
    match lang.target {
        | Target::Wgsl => {
            let _ = writeln!(s, "@compute @workgroup_size({w})");
            let _ = writeln!(
                s,
                "fn {entry}(@builtin(global_invocation_id) gid: vec3<u32>, \
                 @builtin(num_workgroups) nwg: vec3<u32>) {{"
            );
            let _ = writeln!(s, "    let row = gid.y * (nwg.x * {w}u) + gid.x;");
            s.push_str("    let rssn_rows = rssn_dims.rows;\n");
            s.push_str("    if (row >= rssn_rows) { return; }\n");
            s.push_str("    let irow = i32(row);\n");
        },
        | Target::Glsl => {
            let _ = writeln!(s, "void {entry}() {{");
            let _ = writeln!(
                s,
                "    uint row = gl_GlobalInvocationID.y * (gl_NumWorkGroups.x * {w}u) + \
                 gl_GlobalInvocationID.x;"
            );
            s.push_str("    if (row >= rssn_rows) { return; }\n");
            s.push_str("    int irow = int(row);\n");
        },
        | Target::OpenCl => {
            let _ = writeln!(
                s,
                "__kernel void {entry}(__global const {real}* rssn_in, __global const {real}* \
                 rssn_params, __global {real}* rssn_out, __global const uint* rssn_cols, const \
                 uint rssn_rows, const uint rssn_inputs) {{"
            );
            s.push_str("    const uint row = (uint)get_global_id(0);\n");
            s.push_str("    if (row >= rssn_rows) { return; }\n");
            s.push_str("    const int irow = (int)row;\n");
        },
        | Target::Cuda => {
            let _ = writeln!(
                s,
                "extern \"C\" __global__ void {entry}(const {real}* rssn_in, const {real}* \
                 rssn_params, {real}* rssn_out, const unsigned int* rssn_cols, unsigned int \
                 rssn_rows, unsigned int rssn_inputs) {{"
            );
            s.push_str("    const unsigned int row = blockIdx.x * blockDim.x + threadIdx.x;\n");
            s.push_str("    if (row >= rssn_rows) { return; }\n");
            s.push_str("    const int irow = (int)row;\n");
        },
        | Target::Metal => {
            let _ = writeln!(
                s,
                "kernel void {entry}(device const float* rssn_in [[buffer(0)]], device const \
                 float* rssn_params [[buffer(1)]], device float* rssn_out [[buffer(2)]], constant \
                 uint* rssn_dims [[buffer(3)]], device const uint* rssn_cols [[buffer(4)]], uint \
                 row [[thread_position_in_grid]]) {{"
            );
            s.push_str("    const uint rssn_rows = rssn_dims[0];\n");
            s.push_str("    const uint rssn_inputs = rssn_dims[1];\n");
            s.push_str("    if (row >= rssn_rows) { return; }\n");
            s.push_str("    const int irow = (int)row;\n");
        },
    }
}

/// The three phases and the output stores.
fn body(
    s: &mut String,
    lang: Lang,
    ir: &KernelIr,
    analysis: &Analysis,
) {
    let line = |s: &mut String, indent: &str, d: usize| {
        let Some(op) = ir.ops.get(d) else {
            return;
        };
        let ty = match analysis.ty.get(d) {
            | Some(Ty::Idx) => lang.int(),
            | _ => lang.real(),
        };
        let expr = expression(lang, ir, analysis, op);
        let _ = writeln!(s, "{indent}{}", lang.bind(ty, &format!("v{d}"), &expr));
    };
    for &d in &analysis.pre {
        line(s, "    ", d);
    }
    if let Some(fold) = &ir.fold {
        for (k, (op, _)) in fold.accs.iter().enumerate() {
            let _ = writeln!(
                s,
                "    {}",
                lang.var(lang.real(), &format!("a{k}"), &lang.lit(op.identity()))
            );
        }
        if lang.wgsl() {
            let _ = writeln!(
                s,
                "    for (var j: i32 = 0; j < v{}; j = j + 1) {{",
                fold.count
            );
        } else {
            let _ = writeln!(s, "    for (int j = 0; j < v{}; j = j + 1) {{", fold.count);
        }
        for &d in &analysis.body {
            line(s, "        ", d);
        }
        for (k, &(op, r)) in fold.accs.iter().enumerate() {
            let combined = match op {
                | ReduceOp::Sum => format!("a{k} + v{r}"),
                | ReduceOp::Prod => format!("a{k} * v{r}"),
                | ReduceOp::Min => format!("rssn_min(a{k}, v{r})"),
                | ReduceOp::Max => format!("rssn_max(a{k}, v{r})"),
            };
            let _ = writeln!(s, "        a{k} = {combined};");
        }
        s.push_str("    }\n");
        for &d in &analysis.post {
            line(s, "    ", d);
        }
    }
    for (k, &r) in ir.outputs.iter().enumerate() {
        let _ = writeln!(s, "    rssn_out[{k}u * rssn_rows + row] = v{r};");
    }
}

/// The right-hand side of one instruction.
#[allow(clippy::too_many_lines)]
fn expression(
    lang: Lang,
    ir: &KernelIr,
    analysis: &Analysis,
    op: &KernelOp,
) -> String {
    let real = lang.real();
    let int = lang.int();
    let one = lang.lit(1.0);
    let zero = lang.lit(0.0);
    let load = |col: &str, at: &str| {
        if lang.wgsl() || lang.glsl() {
            format!("rssn_load({col}, {at})")
        } else {
            format!("RSSN_LOAD({col}, {at})")
        }
    };
    let truthy = |c: Reg| format!("(v{c} != {zero}) || (v{c} != v{c})");
    match *op {
        | KernelOp::Const(c) => lang.lit(c),
        | KernelOp::Input(i) => load(&format!("{i}u"), "irow"),
        | KernelOp::Param(i) => format!("rssn_params[{i}u]"),
        | KernelOp::Add(a, b) => format!("v{a} + v{b}"),
        | KernelOp::Sub(a, b) => format!("v{a} - v{b}"),
        | KernelOp::Mul(a, b) => format!("v{a} * v{b}"),
        | KernelOp::Div(a, b) => format!("v{a} / v{b}"),
        | KernelOp::Neg(a) => format!("-v{a}"),
        | KernelOp::PowI(a, n) => {
            if lang.wgsl() {
                format!("rssn_powi(v{a}, {n}i)")
            } else {
                format!("rssn_powi(v{a}, {n})")
            }
        },
        | KernelOp::PowF(a, b) => {
            if lang.wraps_pow() {
                format!("rssn_powf(v{a}, v{b})")
            } else {
                format!("pow(v{a}, v{b})")
            }
        },
        | KernelOp::Unary(f, a) => format!("{}(v{a})", lang.intrinsic(f)),
        | KernelOp::Atan2(a, b) => {
            if lang.glsl() {
                format!("atan(v{a}, v{b})")
            } else {
                format!("atan2(v{a}, v{b})")
            }
        },
        | KernelOp::Cmp(c, a, b) | KernelOp::ICmp(c, a, b) => {
            lang.pick(&format!("v{a} {} v{b}", cmp_symbol(c)), &one, &zero)
        },
        | KernelOp::Select(c, a, b) => lang.pick(&truthy(c), &format!("v{a}"), &format!("v{b}")),
        | KernelOp::Gather(base, off) => {
            if lang.wgsl() || lang.glsl() {
                format!("rssn_gather({base}u, v{off}, irow)")
            } else {
                format!("RSSN_GATHER({base}u, v{off}, irow)")
            }
        },
        | KernelOp::Row => "irow".into(),
        | KernelOp::Rows => lang.cast(int, "rssn_rows"),
        | KernelOp::Len(i) => lang.cast(int, &format!("rssn_cols[{}u]", 2 * u64::from(i) + 1)),
        | KernelOp::Inner => "j".into(),
        | KernelOp::IConst(c) => {
            if lang.wgsl() {
                format!("{c}i")
            } else {
                format!("{c}")
            }
        },
        | KernelOp::IAdd(a, b) => format!("v{a} + v{b}"),
        | KernelOp::ISub(a, b) => format!("v{a} - v{b}"),
        | KernelOp::IMul(a, b) => format!("v{a} * v{b}"),
        | KernelOp::IDiv(a, b) => format!("rssn_idiv(v{a}, v{b})"),
        | KernelOp::IMod(a, b) => format!("rssn_imod(v{a}, v{b})"),
        | KernelOp::IMin(a, b) => {
            lang.pick(&format!("v{a} < v{b}"), &format!("v{a}"), &format!("v{b}"))
        },
        | KernelOp::IMax(a, b) => {
            lang.pick(&format!("v{a} > v{b}"), &format!("v{a}"), &format!("v{b}"))
        },
        | KernelOp::ToNum(a) => lang.cast(real, &format!("v{a}")),
        | KernelOp::ToIdx(a) => format!("rssn_toidx(v{a})"),
        | KernelOp::Load(col, at) => load(&format!("{col}u"), &format!("v{at}")),
        | KernelOp::Acc(k) => {
            // Reading an accumulator of a program without one is rejected by
            // `analyse`; the bound keeps the emitted name honest regardless.
            let _ = (ir, analysis);
            format!("a{k}")
        },
    }
}
