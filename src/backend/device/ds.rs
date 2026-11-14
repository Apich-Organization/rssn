//! Double-single arithmetic for WGSL: `f64`-class kernels on devices that
//! only have `f32`.
//!
//! A number is a pair `vec2<f32>(hi, lo)` with `hi + lo` the value and
//! `|lo| <= ulp(hi) / 2`, which carries about 48 significant bits (relative
//! rounding error near `2^-46`, against `2^-53` for an IEEE double and `2^-24`
//! for a float). [`emit_wgsl_ds`] walks a [`KernelIr`] exactly like
//! [`super::emit::emit`] does and spells every operation as a call of the
//! library below, so the interface of the generated module is the one of
//! [`super::emit`] with every element of `rssn_in`, `rssn_params` and
//! `rssn_out` a pair: element `i` of a column occupies the 8 bytes of
//! `array<vec2<f32>>[i]`. The host converts on upload and download
//! ([`split`], [`join`]).
//!
//! # Accuracy and limits
//!
//! * `+ - * /` and `sqrt` are correct to a few units of `2^-46`; the splitting
//!   used for products does not rely on a fused multiply-add, so it holds
//!   whether or not the compiler contracts `a * b + c`.
//! * `exp`, `ln`, `sin`, `cos`, `atan` and everything derived from them are
//!   evaluated in double-single (range reduction plus a Taylor series, or a
//!   Newton step from the `f32` result), about `1e-13` relative for moderate
//!   arguments. `sin`/`cos`/`tan` reduce by `pi/2` held to double-single, so
//!   the absolute error grows like `|x| * 1e-14`; beyond `|x| = 1e6` they
//!   return NaN.
//! * The exponent range is the one of `f32`: magnitudes above `~1e38`
//!   overflow and below `~1e-38` flush, and products of operands above
//!   `~1e34` lose their low half (the splitting constant overflows). Do not
//!   use it for data outside `[1e-30, 1e30]` in magnitude.
//! * As on every shading-language target, NaN/infinity propagation is only
//!   best effort (see [`super::emit`]).

use std::fmt::Write as _;

use super::DeviceError;
use super::KernelIr;
use super::KernelOp;
use super::ReduceOp;
use super::emit::EmitOptions;
use super::kernel::Analysis;
use super::kernel::Ty;
use crate::backend::Intrinsic;
use crate::backend::tape::Cmp;
use crate::backend::tape::Reg;

/// The pair `(hi, lo)` of an `f64`: `hi` is the nearest `f32`, `lo` the
/// nearest `f32` to the remainder.
#[must_use]
#[allow(clippy::cast_possible_truncation)]
pub fn split(x: f64) -> [f32; 2] {
    let hi = x as f32;
    if !hi.is_finite() {
        return [hi, 0.0];
    }
    [hi, (x - f64::from(hi)) as f32]
}

/// The `f64` a pair stands for.
#[must_use]
pub fn join(pair: [f32; 2]) -> f64 {
    f64::from(pair[0]) + f64::from(pair[1])
}

/// A WGSL expression for a constant.
fn lit(c: f64) -> String {
    if c.is_nan() {
        return "vec2<f32>(bitcast<f32>(0x7fc00000u), bitcast<f32>(0x7fc00000u))".into();
    }
    let [hi, lo] = split(c);
    if hi == f32::INFINITY {
        return "vec2<f32>(bitcast<f32>(0x7f800000u), 0.0)".into();
    }
    if hi == f32::NEG_INFINITY {
        return "vec2<f32>(bitcast<f32>(0xff800000u), 0.0)".into();
    }
    format!("vec2<f32>({hi:?}, {lo:?})")
}

/// `c0 + x (c1 + x (c2 + ...))` in double-single, `coefficients` ascending.
fn horner(
    x: &str,
    coefficients: &[f64],
) -> String {
    let mut acc = lit(*coefficients.last().unwrap_or(&0.0));
    for &c in coefficients.iter().rev().skip(1) {
        acc = format!("ds_add(ds_mul({acc}, {x}), {})", lit(c));
    }
    acc
}

fn factorial(n: usize) -> f64 {
    (1..=n).map(|k| k as f64).product()
}

/// The double-single library: arithmetic, then the elementary functions.
/// Constants are filled in from `f64` values so the table is exact to the
/// pair.
#[allow(clippy::too_many_lines)]
fn library() -> String {
    let mut s = String::new();
    s.push_str(
        // The error-free transformations rely on every `+`, `-` and `*` being
        // rounded to `f32` exactly where it is written. Shader compilers may
        // reassociate (`x - (x - y)` becomes `y`) or contract `a * b + c` into
        // a fused multiply-add, and the NVIDIA Vulkan driver does both, which
        // silently reduces the pair to one float. Passing a rounded value
        // through an integer OR with a zero that is only known at run time
        // (it comes from the uniform block) forbids both rewrites, because
        // the optimiser cannot see through integer operations.
        "fn ds_hide(x: f32) -> f32 { return bitcast<f32>(bitcast<u32>(x) | rssn_dims.pad); }\n\
         fn ds_two_sum(a: f32, b: f32) -> vec2<f32> {\n    let s = ds_hide(a + b);\n    let bb = \
         ds_hide(s - a);\n    let u = ds_hide(s - bb);\n    let e = ds_hide(a - u) + ds_hide(b - bb);\n    return vec2<f32>(s, e);\n}\n\
         fn ds_quick(a: f32, b: f32) -> vec2<f32> {\n    let s = ds_hide(a + b);\n    return \
         vec2<f32>(s, b - ds_hide(s - a));\n}\n\
         fn ds_split(a: f32) -> vec2<f32> {\n    let t = ds_hide(4097.0 * a);\n    let d = \
         ds_hide(t - a);\n    let hi = t - d;\n    return vec2<f32>(hi, a - hi);\n}\n\
         fn ds_two_prod(a: f32, b: f32) -> vec2<f32> {\n    let p = ds_hide(a * b);\n    let sa = \
         ds_split(a);\n    let sb = ds_split(b);\n    let m1 = ds_hide(sa.x * sb.x);\n    let m2 = \
         ds_hide(sa.x * sb.y);\n    let m3 = ds_hide(sa.y * sb.x);\n    let m4 = ds_hide(sa.y * \
         sb.y);\n    let r = ds_hide(m1 - p);\n    let e = ((r + m2) + m3) + m4;\n    return vec2<f32>(p, e);\n}\n\
         fn ds_neg(a: vec2<f32>) -> vec2<f32> { return vec2<f32>(-a.x, -a.y); }\n\
         fn ds_add(a: vec2<f32>, b: vec2<f32>) -> vec2<f32> {\n    let s = ds_two_sum(a.x, b.x);\n    \
         let t = ds_two_sum(a.y, b.y);\n    let c = s.y + t.x;\n    let v = ds_quick(s.x, c);\n    \
         let w = t.y + v.y;\n    return ds_quick(v.x, w);\n}\n\
         fn ds_sub(a: vec2<f32>, b: vec2<f32>) -> vec2<f32> { return ds_add(a, ds_neg(b)); }\n\
         fn ds_mul(a: vec2<f32>, b: vec2<f32>) -> vec2<f32> {\n    let p = ds_two_prod(a.x, b.x);\n    \
         let q = p.y + (ds_hide(a.x * b.y) + ds_hide(a.y * b.x));\n    return ds_quick(p.x, q);\n}\n\
         fn ds_mulf(a: vec2<f32>, b: f32) -> vec2<f32> {\n    let p = ds_two_prod(a.x, b);\n    \
         let q = p.y + ds_hide(a.y * b);\n    return ds_quick(p.x, q);\n}\n\
         fn ds_div(a: vec2<f32>, b: vec2<f32>) -> vec2<f32> {\n    let q1 = a.x / b.x;\n    \
         if (q1 != q1 || abs(q1) > 1.0e38 || b.x == 0.0) { return vec2<f32>(q1, 0.0); }\n    \
         let r1 = ds_sub(a, ds_mulf(b, q1));\n    let q2 = r1.x / b.x;\n    let r2 = ds_sub(r1, \
         ds_mulf(b, q2));\n    let q3 = r2.x / b.x;\n    return ds_add(ds_quick(q1, q2), \
         vec2<f32>(q3, 0.0));\n}\n\
         fn ds_sqrt(a: vec2<f32>) -> vec2<f32> {\n    if (a.x != a.x || a.x < 0.0) { return \
         vec2<f32>(bitcast<f32>(0x7fc00000u), bitcast<f32>(0x7fc00000u)); }\n    if (a.x == 0.0) { \
         return vec2<f32>(0.0, 0.0); }\n    let s = sqrt(a.x);\n    if (s > 1.0e18) { return \
         vec2<f32>(s, 0.0); }\n    let r = ds_sub(a, ds_two_prod(s, s));\n    return ds_quick(s, r.x \
         / (2.0 * s));\n}\n\
         fn ds_lt(a: vec2<f32>, b: vec2<f32>) -> bool { return a.x < b.x || (a.x == b.x && a.y < \
         b.y); }\n\
         fn ds_le(a: vec2<f32>, b: vec2<f32>) -> bool { return a.x < b.x || (a.x == b.x && a.y <= \
         b.y); }\n\
         fn ds_eq(a: vec2<f32>, b: vec2<f32>) -> bool { return a.x == b.x && a.y == b.y; }\n\
         fn ds_nan(a: vec2<f32>) -> bool { return a.x != a.x; }\n\
         fn ds_abs(a: vec2<f32>) -> vec2<f32> { if (a.x < 0.0 || (a.x == 0.0 && a.y < 0.0)) { \
         return ds_neg(a); } return a; }\n\
         fn ds_floor(a: vec2<f32>) -> vec2<f32> {\n    let f = floor(a.x);\n    if (f == a.x) { \
         return ds_quick(f, floor(a.y)); }\n    return vec2<f32>(f, 0.0);\n}\n\
         fn ds_ceil(a: vec2<f32>) -> vec2<f32> { return ds_neg(ds_floor(ds_neg(a))); }\n\
         fn ds_pow2(k: f32) -> f32 { return bitcast<f32>(u32(i32(k) + 127) << 23u); }\n",
    );
    let ln2 = std::f64::consts::LN_2;
    let pi = std::f64::consts::PI;
    // exp: exact range reduction by `k ln 2`, then a degree-13 Taylor series on
    // |r| <= ln 2 / 2.
    let exp_coefficients: Vec<f64> = (0..=13).map(|n| 1.0 / factorial(n)).collect();
    let _ = write!(
        s,
        "fn ds_exp(x: vec2<f32>) -> vec2<f32> {{\n    if (x.x != x.x) {{ return x; }}\n    if (x.x > \
         88.0) {{ return vec2<f32>(bitcast<f32>(0x7f800000u), 0.0); }}\n    if (x.x < -87.0) {{ \
         return vec2<f32>(0.0, 0.0); }}\n    let k = floor(x.x * 1.4426950408889634 + 0.5);\n    \
         let r = ds_sub(x, ds_mulf({ln2}, k));\n    let p = {poly};\n    let sc = ds_pow2(k);\n    \
         return vec2<f32>(p.x * sc, p.y * sc);\n}}\n",
        ln2 = lit(ln2),
        poly = horner("r", &exp_coefficients),
    );
    let _ = write!(
        s,
        "fn ds_ln(x: vec2<f32>) -> vec2<f32> {{\n    if (x.x != x.x || x.x < 0.0) {{ return \
         vec2<f32>(bitcast<f32>(0x7fc00000u), bitcast<f32>(0x7fc00000u)); }}\n    if (x.x == 0.0) {{ \
         return vec2<f32>(bitcast<f32>(0xff800000u), 0.0); }}\n    if (x.x > 1.0e38) {{ return \
         vec2<f32>(log(x.x), 0.0); }}\n    var xs = x;\n    var shift = 0.0;\n    if (x.x < 1.0e-30) \
         {{ xs = vec2<f32>(x.x * 1.8446744073709552e19, x.y * 1.8446744073709552e19); shift = \
         64.0; }}\n    if (x.x > 1.0e30) {{ xs = vec2<f32>(x.x * 5.421010862427522e-20, x.y * \
         5.421010862427522e-20); shift = -64.0; }}\n    let y0 = log(xs.x);\n    let e = \
         ds_exp(vec2<f32>(y0, 0.0));\n    let d = ds_div(ds_sub(xs, e), ds_add(xs, e));\n    let y = \
         ds_add(vec2<f32>(y0, 0.0), ds_mulf(d, 2.0));\n    return ds_add(y, ds_mulf({ln2}, \
         shift));\n}}\n",
        ln2 = lit(ln2),
    );
    // sin and cos on |r| <= pi / 4.
    let sin_coefficients: Vec<f64> = (0..=8)
        .map(|n| {
            let sign = if n % 2 == 0 { 1.0 } else { -1.0 };
            sign / factorial(2 * n + 1)
        })
        .collect();
    let cos_coefficients: Vec<f64> = (0..=9)
        .map(|n| {
            let sign = if n % 2 == 0 { 1.0 } else { -1.0 };
            sign / factorial(2 * n)
        })
        .collect();
    let _ = write!(
        s,
        "fn ds_sincos(x: vec2<f32>) -> vec4<f32> {{\n    let nan = bitcast<f32>(0x7fc00000u);\n    if \
         (!(abs(x.x) < 1.0e6)) {{ return vec4<f32>(nan, nan, nan, nan); }}\n    let k = floor(x.x * \
         0.6366197723675814 + 0.5);\n    let r = ds_sub(x, ds_mulf({pio2}, k));\n    let r2 = \
         ds_mul(r, r);\n    let s = ds_mul(r, {sin});\n    let c = {cos};\n    let q = ((i32(k) % \
         4) + 4) % 4;\n    if (q == 0) {{ return vec4<f32>(s.x, s.y, c.x, c.y); }}\n    if (q == 1) \
         {{ return vec4<f32>(c.x, c.y, -s.x, -s.y); }}\n    if (q == 2) {{ return vec4<f32>(-s.x, \
         -s.y, -c.x, -c.y); }}\n    return vec4<f32>(-c.x, -c.y, s.x, s.y);\n}}\n",
        pio2 = lit(pi / 2.0),
        sin = horner("r2", &sin_coefficients),
        cos = horner("r2", &cos_coefficients),
    );
    let _ = write!(
        s,
        "fn ds_sin(x: vec2<f32>) -> vec2<f32> {{ let t = ds_sincos(x); return vec2<f32>(t.x, t.y); \
         }}\nfn ds_cos(x: vec2<f32>) -> vec2<f32> {{ let t = ds_sincos(x); return vec2<f32>(t.z, \
         t.w); }}\nfn ds_tan(x: vec2<f32>) -> vec2<f32> {{ let t = ds_sincos(x); return \
         ds_div(vec2<f32>(t.x, t.y), vec2<f32>(t.z, t.w)); }}\n\
         fn ds_atan(x: vec2<f32>) -> vec2<f32> {{\n    if (x.x != x.x) {{ return x; }}\n    let sg = \
         select(-1.0, 1.0, x.x >= 0.0);\n    if (abs(x.x) > 1.0e18) {{ return ds_mulf({pio2}, sg); \
         }}\n    let y0 = atan(x.x);\n    let t = ds_sincos(vec2<f32>(y0, 0.0));\n    let s = \
         vec2<f32>(t.x, t.y);\n    let c = vec2<f32>(t.z, t.w);\n    let d = ds_div(ds_sub(ds_mul(x, \
         c), s), ds_add(c, ds_mul(x, s)));\n    return ds_add(vec2<f32>(y0, 0.0), d);\n}}\n\
         fn ds_atan2(y: vec2<f32>, x: vec2<f32>) -> vec2<f32> {{\n    if (x.x == 0.0 && y.x == \
         0.0) {{ return vec2<f32>(0.0, 0.0); }}\n    if (x.x > 0.0) {{ return ds_atan(ds_div(y, \
         x)); }}\n    if (x.x < 0.0) {{\n        let a = ds_atan(ds_div(y, x));\n        if (y.x >= \
         0.0) {{ return ds_add(a, {pi}); }}\n        return ds_sub(a, {pi});\n    }}\n    if (y.x > \
         0.0) {{ return {pio2}; }}\n    return ds_neg({pio2});\n}}\n\
         fn ds_asin(x: vec2<f32>) -> vec2<f32> {{\n    if (x.x != x.x || abs(x.x) > 1.0) {{ return \
         vec2<f32>(bitcast<f32>(0x7fc00000u), bitcast<f32>(0x7fc00000u)); }}\n    if (abs(x.x) == \
         1.0) {{ return ds_mulf({pio2}, x.x); }}\n    let one = vec2<f32>(1.0, 0.0);\n    return \
         ds_atan(ds_div(x, ds_sqrt(ds_mul(ds_sub(one, x), ds_add(one, x)))));\n}}\n\
         fn ds_acos(x: vec2<f32>) -> vec2<f32> {{ return ds_sub({pio2}, ds_asin(x)); }}\n",
        pio2 = lit(pi / 2.0),
        pi = lit(pi),
    );
    // Hyperbolic functions: a short series for small arguments, `exp` above.
    let sinh_coefficients: Vec<f64> = (0..=6).map(|n| 1.0 / factorial(2 * n + 1)).collect();
    let cosh_coefficients: Vec<f64> = (0..=6).map(|n| 1.0 / factorial(2 * n)).collect();
    let _ = write!(
        s,
        "fn ds_sinh(x: vec2<f32>) -> vec2<f32> {{\n    if (abs(x.x) < 0.1) {{ let x2 = ds_mul(x, \
         x); return ds_mul(x, {sinh}); }}\n    let e = ds_exp(x);\n    let ei = ds_div(vec2<f32>(1.0, \
         0.0), e);\n    return ds_mulf(ds_sub(e, ei), 0.5);\n}}\n\
         fn ds_cosh(x: vec2<f32>) -> vec2<f32> {{\n    if (abs(x.x) < 0.1) {{ let x2 = ds_mul(x, \
         x); return {cosh}; }}\n    let e = ds_exp(x);\n    let ei = ds_div(vec2<f32>(1.0, 0.0), \
         e);\n    return ds_mulf(ds_add(e, ei), 0.5);\n}}\n\
         fn ds_tanh(x: vec2<f32>) -> vec2<f32> {{\n    if (x.x > 20.0) {{ return vec2<f32>(1.0, \
         0.0); }}\n    if (x.x < -20.0) {{ return vec2<f32>(-1.0, 0.0); }}\n    if (abs(x.x) < 0.1) \
         {{ return ds_div(ds_sinh(x), ds_cosh(x)); }}\n    let e = ds_exp(ds_mulf(x, 2.0));\n    let \
         one = vec2<f32>(1.0, 0.0);\n    return ds_div(ds_sub(e, one), ds_add(e, one));\n}}\n",
        sinh = horner("x2", &sinh_coefficients),
        cosh = horner("x2", &cosh_coefficients),
    );
    let _ = write!(
        s,
        "fn ds_asinh(x: vec2<f32>) -> vec2<f32> {{\n    if (x.x != x.x) {{ return x; }}\n    let sg = \
         select(-1.0, 1.0, x.x >= 0.0);\n    let ax = ds_abs(x);\n    if (ax.x > 1.0e18) {{ return \
         ds_mulf(ds_ln(ds_mulf(ax, 2.0)), sg); }}\n    var y = vec2<f32>(select(log(ax.x + \
         sqrt(ax.x * ax.x + 1.0)), ax.x, ax.x < 1.0e-3), 0.0);\n    for (var i = 0; i < 2; i = i + \
         1) {{\n        y = ds_sub(y, ds_div(ds_sub(ds_sinh(y), ax), ds_cosh(y)));\n    }}\n    \
         return ds_mulf(y, sg);\n}}\n\
         fn ds_acosh(x: vec2<f32>) -> vec2<f32> {{\n    if (x.x != x.x || x.x < 1.0) {{ return \
         vec2<f32>(bitcast<f32>(0x7fc00000u), bitcast<f32>(0x7fc00000u)); }}\n    if (x.x == 1.0 \
         && x.y <= 0.0) {{ return vec2<f32>(0.0, 0.0); }}\n    if (x.x > 1.0e18) {{ return \
         ds_ln(ds_mulf(x, 2.0)); }}\n    var y = vec2<f32>(select(sqrt(max(2.0 * (x.x - 1.0), 1.0e-12)), log(x.x + \
         sqrt((x.x - 1.0) * (x.x + 1.0))), x.x > 1.1), 0.0);\n    for (var i = 0; i < 3; i = i + \
         1) {{\n        let sh = ds_sinh(y);\n        if (sh.x == 0.0) {{ break; }}\n        y = \
         ds_sub(y, ds_div(ds_sub(ds_cosh(y), x), sh));\n    }}\n    return y;\n}}\n\
         fn ds_atanh(x: vec2<f32>) -> vec2<f32> {{\n    if (x.x != x.x || abs(x.x) > 1.0) {{ return \
         vec2<f32>(bitcast<f32>(0x7fc00000u), bitcast<f32>(0x7fc00000u)); }}\n    if (abs(x.x) == \
         1.0) {{ return vec2<f32>(x.x * bitcast<f32>(0x7f800000u), 0.0); }}\n    var y = \
         vec2<f32>(select(0.5 * log((1.0 + x.x) / (1.0 - x.x)), x.x, abs(x.x) < 1.0e-4), 0.0);\n    \
         for (var i = 0; i < 2; i = i + 1) {{\n        let t = ds_tanh(y);\n        let one = \
         vec2<f32>(1.0, 0.0);\n        y = ds_sub(y, ds_div(ds_sub(t, x), ds_sub(one, \
         ds_mul(t, t))));\n    }}\n    return y;\n}}\n\
         fn ds_cbrt(x: vec2<f32>) -> vec2<f32> {{\n    if (x.x == 0.0 || x.x != x.x) {{ return x; \
         }}\n    let sg = select(-1.0, 1.0, x.x > 0.0);\n    let ax = ds_abs(x);\n    var y = \
         vec2<f32>(pow(ax.x, 0.3333333333333333), 0.0);\n    let three = vec2<f32>(3.0, 0.0);\n    \
         for (var i = 0; i < 2; i = i + 1) {{\n        y = ds_div(ds_add(ds_mulf(y, 2.0), \
         ds_div(ax, ds_mul(y, y))), three);\n    }}\n    return ds_mulf(y, sg);\n}}\n\
         fn ds_min(a: vec2<f32>, b: vec2<f32>) -> vec2<f32> {{\n    if (ds_nan(a) || ds_nan(b)) {{ \
         return vec2<f32>(bitcast<f32>(0x7fc00000u), bitcast<f32>(0x7fc00000u)); }}\n    if \
         (ds_lt(a, b)) {{ return a; }}\n    return b;\n}}\n\
         fn ds_max(a: vec2<f32>, b: vec2<f32>) -> vec2<f32> {{\n    if (ds_nan(a) || ds_nan(b)) {{ \
         return vec2<f32>(bitcast<f32>(0x7fc00000u), bitcast<f32>(0x7fc00000u)); }}\n    if \
         (ds_lt(b, a)) {{ return a; }}\n    return b;\n}}\n",
    );
    s
}

/// Helpers that depend on what the program uses.
fn helpers(
    s: &mut String,
    ir: &KernelIr,
) {
    let uses = |f: &dyn Fn(&KernelOp) -> bool| ir.ops.iter().any(f);
    let nan = lit(f64::NAN);
    if uses(&|op| matches!(op, KernelOp::PowI(..))) {
        s.push_str(
            "fn rssn_powi(x: vec2<f32>, n: i32) -> vec2<f32> {\n    var e = select(u32(n), u32(-n), n \
             < 0);\n    var b = x;\n    var acc = vec2<f32>(1.0, 0.0);\n    loop {\n        if (e == \
             0u) { break; }\n        if ((e & 1u) == 1u) { acc = ds_mul(acc, b); }\n        e = e >> \
             1u;\n        if (e > 0u) { b = ds_mul(b, b); }\n    }\n    if (n < 0) { return \
             ds_div(vec2<f32>(1.0, 0.0), acc); }\n    return acc;\n}\n",
        );
    }
    if uses(&|op| matches!(op, KernelOp::PowF(..))) {
        let _ = write!(
            s,
            "fn rssn_powf(x: vec2<f32>, y: vec2<f32>) -> vec2<f32> {{\n    if (y.x == 0.0 && y.y == \
             0.0) {{ return vec2<f32>(1.0, 0.0); }}\n    if (x.x > 0.0) {{ return \
             ds_exp(ds_mul(y, ds_ln(x))); }}\n    if (x.x == 0.0) {{\n        if (y.x > 0.0) {{ \
             return vec2<f32>(0.0, 0.0); }}\n        return vec2<f32>(bitcast<f32>(0x7f800000u), \
             0.0);\n    }}\n    if (y.x == floor(y.x) && y.y == 0.0) {{\n        let m = \
             ds_exp(ds_mul(y, ds_ln(ds_neg(x))));\n        if (fract(y.x * 0.5) != 0.0) {{ return \
             ds_neg(m); }}\n        return m;\n    }}\n    return {nan};\n}}\n"
        );
    }
    if uses(&|op| matches!(op, KernelOp::ToIdx(_))) {
        s.push_str(
            "fn rssn_toidx(x: vec2<f32>) -> i32 {\n    if (x.x != x.x) { return 0; }\n    if (x.x >= \
             2147483520.0) { return 2147483647; }\n    if (x.x <= -2147483520.0) { return \
             -2147483647; }\n    let t = i32(x.x);\n    let r = (x.x - f32(t)) + x.y;\n    if (x.x \
             >= 0.0) { return t + i32(floor(r)); }\n    return t + i32(ceil(r));\n}\n\
             fn rssn_tonum(v: i32) -> vec2<f32> {\n    let h = f32(v);\n    return vec2<f32>(h, \
             f32(v - i32(h)));\n}\n",
        );
    }
    if uses(&|op| matches!(op, KernelOp::ToNum(_))) && !uses(&|op| matches!(op, KernelOp::ToIdx(_)))
    {
        s.push_str(
            "fn rssn_tonum(v: i32) -> vec2<f32> {\n    let h = f32(v);\n    return vec2<f32>(h, \
             f32(v - i32(h)));\n}\n",
        );
    }
    if uses(&|op| matches!(op, KernelOp::IDiv(..))) {
        s.push_str(
            "fn rssn_idiv(a: i32, b: i32) -> i32 {\n    if (b == 0) { return 0; }\n    return a / \
             b;\n}\n",
        );
    }
    if uses(&|op| matches!(op, KernelOp::IMod(..))) {
        s.push_str(
            "fn rssn_imod(a: i32, b: i32) -> i32 {\n    if (b == 0) { return 0; }\n    return a % \
             b;\n}\n",
        );
    }
    if uses(&|op| {
        matches!(
            op,
            KernelOp::Input(_) | KernelOp::Load(..) | KernelOp::Gather(..)
        )
    }) {
        let _ = write!(
            s,
            "fn rssn_load(col: u32, i: i32) -> vec2<f32> {{\n    let len = rssn_cols[2u * col + \
             1u];\n    if (i < 0 || u32(i) >= len) {{ return {nan}; }}\n    return \
             rssn_in[rssn_cols[2u * col] + u32(i)];\n}}\n"
        );
    }
    if uses(&|op| matches!(op, KernelOp::Gather(..))) {
        let _ = write!(
            s,
            "fn rssn_gather(base: u32, off: vec2<f32>, i: i32) -> vec2<f32> {{\n    if (!(off.x >= \
             0.0) || !(off.x < f32(rssn_dims.inputs - base))) {{ return {nan}; }}\n    return \
             rssn_load(base + u32(off.x), i);\n}}\n"
        );
    }
}

const fn intrinsic(f: Intrinsic) -> &'static str {
    match f {
        | Intrinsic::Exp => "ds_exp",
        | Intrinsic::Ln => "ds_ln",
        | Intrinsic::Sin => "ds_sin",
        | Intrinsic::Cos => "ds_cos",
        | Intrinsic::Tan => "ds_tan",
        | Intrinsic::Asin => "ds_asin",
        | Intrinsic::Acos => "ds_acos",
        | Intrinsic::Atan => "ds_atan",
        | Intrinsic::Sinh => "ds_sinh",
        | Intrinsic::Cosh => "ds_cosh",
        | Intrinsic::Tanh => "ds_tanh",
        | Intrinsic::Asinh => "ds_asinh",
        | Intrinsic::Acosh => "ds_acosh",
        | Intrinsic::Atanh => "ds_atanh",
        | Intrinsic::Sqrt => "ds_sqrt",
        | Intrinsic::Abs => "ds_abs",
        | Intrinsic::Floor => "ds_floor",
        | Intrinsic::Ceil => "ds_ceil",
        | Intrinsic::Cbrt => "ds_cbrt",
    }
}

fn compare(
    c: Cmp,
    a: Reg,
    b: Reg,
) -> String {
    match c {
        | Cmp::Lt => format!("ds_lt(v{a}, v{b})"),
        | Cmp::Le => format!("ds_le(v{a}, v{b})"),
        | Cmp::Gt => format!("ds_lt(v{b}, v{a})"),
        | Cmp::Ge => format!("ds_le(v{b}, v{a})"),
        | Cmp::Eq => format!("ds_eq(v{a}, v{b})"),
        | Cmp::Ne => format!("(!ds_eq(v{a}, v{b}))"),
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

/// The right-hand side of one instruction.
#[allow(clippy::too_many_lines)]
fn expression(op: &KernelOp) -> String {
    let one = lit(1.0);
    let zero = lit(0.0);
    let truthy = |c: Reg| format!("(v{c}.x != 0.0 || v{c}.y != 0.0 || v{c}.x != v{c}.x)");
    match *op {
        | KernelOp::Const(c) => lit(c),
        | KernelOp::Input(i) => format!("rssn_load({i}u, irow)"),
        | KernelOp::Param(i) => format!("rssn_params[{i}u]"),
        | KernelOp::Add(a, b) => format!("ds_add(v{a}, v{b})"),
        | KernelOp::Sub(a, b) => format!("ds_sub(v{a}, v{b})"),
        | KernelOp::Mul(a, b) => format!("ds_mul(v{a}, v{b})"),
        | KernelOp::Div(a, b) => format!("ds_div(v{a}, v{b})"),
        | KernelOp::Neg(a) => format!("ds_neg(v{a})"),
        | KernelOp::PowI(a, n) => format!("rssn_powi(v{a}, {n}i)"),
        | KernelOp::PowF(a, b) => format!("rssn_powf(v{a}, v{b})"),
        | KernelOp::Unary(f, a) => format!("{}(v{a})", intrinsic(f)),
        | KernelOp::Atan2(a, b) => format!("ds_atan2(v{a}, v{b})"),
        | KernelOp::Cmp(c, a, b) => format!("select({zero}, {one}, {})", compare(c, a, b)),
        | KernelOp::ICmp(c, a, b) => {
            format!("select({zero}, {one}, v{a} {} v{b})", cmp_symbol(c))
        },
        | KernelOp::Select(c, a, b) => format!("select(v{b}, v{a}, {})", truthy(c)),
        | KernelOp::Gather(base, off) => format!("rssn_gather({base}u, v{off}, irow)"),
        | KernelOp::Row => "irow".into(),
        | KernelOp::Rows => "i32(rssn_rows)".into(),
        | KernelOp::Len(i) => format!("i32(rssn_cols[{}u])", 2 * u64::from(i) + 1),
        | KernelOp::Inner => "j".into(),
        | KernelOp::IConst(c) => format!("{c}i"),
        | KernelOp::IAdd(a, b) => format!("v{a} + v{b}"),
        | KernelOp::ISub(a, b) => format!("v{a} - v{b}"),
        | KernelOp::IMul(a, b) => format!("v{a} * v{b}"),
        | KernelOp::IDiv(a, b) => format!("rssn_idiv(v{a}, v{b})"),
        | KernelOp::IMod(a, b) => format!("rssn_imod(v{a}, v{b})"),
        | KernelOp::IMin(a, b) => format!("select(v{b}, v{a}, v{a} < v{b})"),
        | KernelOp::IMax(a, b) => format!("select(v{b}, v{a}, v{a} > v{b})"),
        | KernelOp::ToNum(a) => format!("rssn_tonum(v{a})"),
        | KernelOp::ToIdx(a) => format!("rssn_toidx(v{a})"),
        | KernelOp::Load(col, at) => format!("rssn_load({col}u, v{at})"),
        | KernelOp::Acc(k) => format!("a{k}"),
    }
}

/// Emits the double-single WGSL module for `ir`.
///
/// The bindings are those of [`super::emit`] with `array<vec2<f32>>` for
/// `rssn_in`, `rssn_params` and `rssn_out`.
///
/// # Errors
/// [`DeviceError::InvalidKernel`] when `ir` does not validate or the
/// workgroup size is zero.
pub fn emit_wgsl_ds(
    ir: &KernelIr,
    options: &EmitOptions,
) -> Result<String, DeviceError> {
    let analysis: Analysis = ir.analyse()?;
    if options.workgroup_size == 0 {
        return Err(DeviceError::InvalidKernel(
            "workgroup size must be positive".into(),
        ));
    }
    let w = options.workgroup_size;
    let entry = &options.entry;
    let mut s = String::new();
    let _ = writeln!(
        s,
        "// rssn kernel: {} inputs, {} params, {} outputs, {} ops{} (double-single f32 pairs)",
        ir.inputs,
        ir.params,
        ir.outputs.len(),
        ir.ops.len(),
        if ir.fold.is_some() {
            ", fold"
        } else {
            ""
        },
    );
    s.push_str("struct Dims { rows: u32, inputs: u32, outputs: u32, pad: u32 }\n");
    s.push_str("@group(0) @binding(0) var<storage, read> rssn_in: array<vec2<f32>>;\n");
    s.push_str("@group(0) @binding(1) var<storage, read> rssn_params: array<vec2<f32>>;\n");
    s.push_str("@group(0) @binding(2) var<storage, read_write> rssn_out: array<vec2<f32>>;\n");
    s.push_str("@group(0) @binding(3) var<uniform> rssn_dims: Dims;\n");
    s.push_str("@group(0) @binding(4) var<storage, read> rssn_cols: array<u32>;\n\n");
    s.push_str(&library());
    helpers(&mut s, ir);
    let _ = writeln!(s, "@compute @workgroup_size({w})");
    let _ = writeln!(
        s,
        "fn {entry}(@builtin(global_invocation_id) gid: vec3<u32>, @builtin(num_workgroups) nwg: \
         vec3<u32>) {{"
    );
    let _ = writeln!(s, "    let row = gid.y * (nwg.x * {w}u) + gid.x;");
    s.push_str("    let rssn_rows = rssn_dims.rows;\n");
    s.push_str("    if (row >= rssn_rows) { return; }\n");
    s.push_str("    let irow = i32(row);\n");
    let line = |s: &mut String, indent: &str, d: usize| {
        let Some(op) = ir.ops.get(d) else {
            return;
        };
        let ty = match analysis.ty.get(d) {
            | Some(Ty::Idx) => "i32",
            | _ => "vec2<f32>",
        };
        let _ = writeln!(s, "{indent}let v{d}: {ty} = {};", expression(op));
    };
    for &d in &analysis.pre {
        line(&mut s, "    ", d);
    }
    if let Some(fold) = &ir.fold {
        for (k, (op, _)) in fold.accs.iter().enumerate() {
            let _ = writeln!(s, "    var a{k}: vec2<f32> = {};", lit(op.identity()));
        }
        let _ = writeln!(
            s,
            "    for (var j: i32 = 0; j < v{}; j = j + 1) {{",
            fold.count
        );
        for &d in &analysis.body {
            line(&mut s, "        ", d);
        }
        for (k, &(op, r)) in fold.accs.iter().enumerate() {
            let combined = match op {
                | ReduceOp::Sum => format!("ds_add(a{k}, v{r})"),
                | ReduceOp::Prod => format!("ds_mul(a{k}, v{r})"),
                | ReduceOp::Min => format!("ds_min(a{k}, v{r})"),
                | ReduceOp::Max => format!("ds_max(a{k}, v{r})"),
            };
            let _ = writeln!(s, "        a{k} = {combined};");
        }
        s.push_str("    }\n");
        for &d in &analysis.post {
            line(&mut s, "    ", d);
        }
    }
    for (k, &r) in ir.outputs.iter().enumerate() {
        let _ = writeln!(s, "    rssn_out[{k}u * rssn_rows + row] = v{r};");
    }
    s.push_str("}\n");
    Ok(s)
}

#[cfg(test)]
mod tests {
    use super::join;
    use super::lit;
    use super::split;

    #[test]
    fn a_pair_carries_about_48_bits() {
        for x in [
            1.0 / 3.0,
            std::f64::consts::PI,
            1.0e10 + 0.123_456_789,
            -std::f64::consts::E,
            1.0e-20,
        ] {
            let back = join(split(x));
            assert!(((back - x) / x).abs() < 1.0e-14, "{x}: {back}");
        }
        // Whole numbers below 2^48 are exact: index data survives.
        for k in [0.0, 1.0, 16_777_217.0, 281_474_976_710_655.0] {
            assert_eq!(join(split(k)), k);
        }
        assert_eq!(split(f64::INFINITY), [f32::INFINITY, 0.0]);
        assert!(join(split(f64::NAN)).is_nan());
    }

    #[test]
    fn constants_are_spelled_as_pairs() {
        assert_eq!(lit(1.0), "vec2<f32>(1.0, 0.0)");
        assert!(lit(std::f64::consts::PI).starts_with("vec2<f32>(3.1415927, "));
        assert!(lit(f64::INFINITY).contains("0x7f800000u"));
        assert!(lit(f64::NAN).contains("0x7fc00000u"));
    }
}
