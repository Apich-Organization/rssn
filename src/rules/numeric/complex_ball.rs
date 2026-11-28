//! Complex ball (`list(re, im, rad)`) and rectangular complex interval
//! (`list(re_lo, re_hi, im_lo, im_hi)`) operators.

use num_complex::Complex64;

use super::bridge::Data;
use super::bridge::bad;
use super::bridge::vector;
use crate::kernel::diag::DiagKind;
use crate::kernels::ball::Cball;
use crate::kernels::complex_ball::ComplexBall;

reviewed_entries! {
    (
        "cball_add", "vv", cball_add, complex_ball::ComplexBall::add_ball,
        "`cball_add(a, b)`: sum of two complex balls `list(re, im, rad)`."
    ),
    (
        "cball_sub", "vv", cball_sub, complex_ball::ComplexBall::sub_ball,
        "`cball_sub(a, b)`: difference of two complex balls `list(re, im, rad)`."
    ),
    (
        "cball_mul", "vv", cball_mul, complex_ball::ComplexBall::mul_ball,
        "`cball_mul(a, b)`: product of two complex balls `list(re, im, rad)`."
    ),
    (
        "cball_div", "vv", cball_div, complex_ball::ComplexBall::div_ball,
        "`cball_div(a, b)`: quotient of two complex balls `list(re, im, rad)` when `b` excludes zero."
    ),
    (
        "cball_pow", "vv", cball_pow, complex_ball::ComplexBall::powc,
        "`cball_pow(a, b)`: principal complex power `a^b` over complex balls `list(re, im, rad)`."
    ),
    (
        "cball_exp", "v", cball_exp, complex_ball::ComplexBall::exp,
        "`cball_exp(z)`: complex exponential enclosure of `list(re, im, rad)`."
    ),
    (
        "cball_ln", "v", cball_ln, complex_ball::ComplexBall::ln,
        "`cball_ln(z)`: principal complex logarithm enclosure of `list(re, im, rad)`."
    ),
    (
        "cball_sqrt", "v", cball_sqrt, complex_ball::ComplexBall::sqrt,
        "`cball_sqrt(z)`: principal complex square root enclosure of `list(re, im, rad)`."
    ),
    (
        "cball_sin", "v", cball_sin, complex_ball::ComplexBall::sin,
        "`cball_sin(z)`: complex sine enclosure of `list(re, im, rad)`."
    ),
    (
        "cball_cos", "v", cball_cos, complex_ball::ComplexBall::cos,
        "`cball_cos(z)`: complex cosine enclosure of `list(re, im, rad)`."
    ),
    (
        "cball_tan", "v", cball_tan, complex_ball::ComplexBall::tan,
        "`cball_tan(z)`: complex tangent enclosure of `list(re, im, rad)`."
    ),
    (
        "cball_sinh", "v", cball_sinh, complex_ball::ComplexBall::sinh,
        "`cball_sinh(z)`: complex hyperbolic sine enclosure of `list(re, im, rad)`."
    ),
    (
        "cball_cosh", "v", cball_cosh, complex_ball::ComplexBall::cosh,
        "`cball_cosh(z)`: complex hyperbolic cosine enclosure of `list(re, im, rad)`."
    ),
    (
        "cball_tanh", "v", cball_tanh, complex_ball::ComplexBall::tanh,
        "`cball_tanh(z)`: complex hyperbolic tangent enclosure of `list(re, im, rad)`."
    ),
    (
        "cball_asin", "v", cball_asin, complex_ball::ComplexBall::asin,
        "`cball_asin(z)`: principal inverse sine enclosure of `list(re, im, rad)`."
    ),
    (
        "cball_acos", "v", cball_acos, complex_ball::ComplexBall::acos,
        "`cball_acos(z)`: principal inverse cosine enclosure of `list(re, im, rad)`."
    ),
    (
        "cball_atan", "v", cball_atan, complex_ball::ComplexBall::atan,
        "`cball_atan(z)`: principal inverse tangent enclosure of `list(re, im, rad)`."
    ),
    (
        "cball_asinh", "v", cball_asinh, complex_ball::ComplexBall::asinh,
        "`cball_asinh(z)`: principal inverse hyperbolic sine enclosure of `list(re, im, rad)`."
    ),
    (
        "cball_acosh", "v", cball_acosh, complex_ball::ComplexBall::acosh,
        "`cball_acosh(z)`: principal inverse hyperbolic cosine enclosure of `list(re, im, rad)`."
    ),
    (
        "cball_atanh", "v", cball_atanh, complex_ball::ComplexBall::atanh,
        "`cball_atanh(z)`: principal inverse hyperbolic tangent enclosure of `list(re, im, rad)`."
    ),
    (
        "cball_conj", "v", cball_conj, complex_ball::ComplexBall::conj,
        "`cball_conj(z)`: complex conjugate enclosure of `list(re, im, rad)`."
    ),
    (
        "cball_powi", "vs", cball_powi, complex_ball::ComplexBall::powi,
        "`cball_powi(z, n)`: integer power enclosure `z^n` of `list(re, im, rad)`."
    ),
    (
        "cball_abs", "v", cball_abs, complex_ball::ComplexBall::abs_ball,
        "`cball_abs(z)`: modulus enclosure `|z|` as a complex ball `list(re, 0, rad)`."
    ),
    (
        "cball_arg", "v", cball_arg, complex_ball::ComplexBall::arg_ball,
        "`cball_arg(z)`: principal argument enclosure `arg(z)` as a complex ball `list(re, 0, rad)`."
    ),
    (
        "cball_contains", "vv", cball_contains,
        "`cball_contains(b, pt)`: returns `1` if `pt` (`list(re, im, 0)` or `list(re, im, rad)`) lies in `b`, else `0`."
    ),
    (
        "cball_overlaps", "vv", cball_overlaps, complex_ball::ComplexBall::overlaps,
        "`cball_overlaps(a, b)`: returns `1` if complex balls `a` and `b` overlap, else `0`."
    ),
    (
        "cball_to_cinterval", "v", cball_to_cinterval, complex_ball::ComplexBall::to_cball,
        "`cball_to_cinterval(b)`: converts a complex ball `list(re, im, rad)` to a rectangular interval `list(re_lo, re_hi, im_lo, im_hi)`."
    ),
    (
        "cinterval_to_cball", "v", cinterval_to_cball, complex_ball::ComplexBall::from_cball,
        "`cinterval_to_cball(rect)`: converts a rectangular complex interval `list(re_lo, re_hi, im_lo, im_hi)` to a complex ball `list(re, im, rad)`."
    ),
    (
        "cinterval_add", "vv", cinterval_add,
        "`cinterval_add(a, b)`: sum of two complex intervals `list(re_lo, re_hi, im_lo, im_hi)`."
    ),
    (
        "cinterval_sub", "vv", cinterval_sub,
        "`cinterval_sub(a, b)`: difference of two complex intervals `list(re_lo, re_hi, im_lo, im_hi)`."
    ),
    (
        "cinterval_mul", "vv", cinterval_mul,
        "`cinterval_mul(a, b)`: product of two complex intervals `list(re_lo, re_hi, im_lo, im_hi)`."
    ),
    (
        "cinterval_div", "vv", cinterval_div, ball::Cball::div,
        "`cinterval_div(a, b)`: quotient of two complex intervals `list(re_lo, re_hi, im_lo, im_hi)`."
    ),
    (
        "cinterval_exp", "v", cinterval_exp, ball::Cball::exp,
        "`cinterval_exp(z)`: complex exponential enclosure of `list(re_lo, re_hi, im_lo, im_hi)`."
    ),
    (
        "cinterval_ln", "v", cinterval_ln, ball::Cball::ln,
        "`cinterval_ln(z)`: principal complex logarithm enclosure of `list(re_lo, re_hi, im_lo, im_hi)`."
    ),
    (
        "cinterval_sqrt", "v", cinterval_sqrt, ball::Cball::sqrt,
        "`cinterval_sqrt(z)`: principal complex square root enclosure of `list(re_lo, re_hi, im_lo, im_hi)`."
    ),
    (
        "cinterval_sin", "v", cinterval_sin, ball::Cball::sin,
        "`cinterval_sin(z)`: complex sine enclosure of `list(re_lo, re_hi, im_lo, im_hi)`."
    ),
    (
        "cinterval_cos", "v", cinterval_cos, ball::Cball::cos,
        "`cinterval_cos(z)`: complex cosine enclosure of `list(re_lo, re_hi, im_lo, im_hi)`."
    ),
    (
        "cinterval_tan", "v", cinterval_tan, ball::Cball::tan,
        "`cinterval_tan(z)`: complex tangent enclosure of `list(re_lo, re_hi, im_lo, im_hi)`."
    ),
    (
        "cinterval_sinh", "v", cinterval_sinh, ball::Cball::sinh,
        "`cinterval_sinh(z)`: complex hyperbolic sine enclosure of `list(re_lo, re_hi, im_lo, im_hi)`."
    ),
    (
        "cinterval_cosh", "v", cinterval_cosh, ball::Cball::cosh,
        "`cinterval_cosh(z)`: complex hyperbolic cosine enclosure of `list(re_lo, re_hi, im_lo, im_hi)`."
    ),
    (
        "cinterval_tanh", "v", cinterval_tanh, ball::Cball::tanh,
        "`cinterval_tanh(z)`: complex hyperbolic tangent enclosure of `list(re_lo, re_hi, im_lo, im_hi)`."
    ),
    (
        "cinterval_pow", "vv", cinterval_pow, ball::Cball::powc,
        "`cinterval_pow(a, b)`: principal complex power enclosure `a^b`."
    ),
    (
        "cinterval_asin", "v", cinterval_asin, ball::Cball::asin,
        "`cinterval_asin(z)`: principal inverse sine enclosure of `list(re_lo, re_hi, im_lo, im_hi)`."
    ),
    (
        "cinterval_acos", "v", cinterval_acos, ball::Cball::acos,
        "`cinterval_acos(z)`: principal inverse cosine enclosure of `list(re_lo, re_hi, im_lo, im_hi)`."
    ),
    (
        "cinterval_atan", "v", cinterval_atan, ball::Cball::atan,
        "`cinterval_atan(z)`: principal inverse tangent enclosure of `list(re_lo, re_hi, im_lo, im_hi)`."
    ),
    (
        "cinterval_asinh", "v", cinterval_asinh, ball::Cball::asinh,
        "`cinterval_asinh(z)`: principal inverse hyperbolic sine enclosure of dimensionless bounds `list(re_lo, re_hi, im_lo, im_hi)`. Symbolic calls stay unevaluated when a bound cannot be formed. Checked calls report a domain diagnostic."
    ),
    (
        "cinterval_acosh", "v", cinterval_acosh, ball::Cball::acosh,
        "`cinterval_acosh(z)`: principal inverse hyperbolic cosine enclosure of dimensionless bounds `list(re_lo, re_hi, im_lo, im_hi)`. Symbolic calls stay unevaluated when a bound cannot be formed. Checked calls report a domain diagnostic."
    ),
    (
        "cinterval_atanh", "v", cinterval_atanh, ball::Cball::atanh,
        "`cinterval_atanh(z)`: principal inverse hyperbolic tangent enclosure of dimensionless bounds `list(re_lo, re_hi, im_lo, im_hi)`. Symbolic calls stay unevaluated when a bound cannot be formed. Checked calls report a domain diagnostic."
    ),
    (
        "cinterval_intersect", "vv", cinterval_intersect, ball::Cball::intersect,
        "`cinterval_intersect(a, b)`: intersection of two rectangular complex intervals."
    ),
}

fn decode_ball(
    args: &[Data],
    at: usize,
) -> Option<ComplexBall> {
    let v = vector(args, at)?;
    if v.len() != 3 {
        return bad(
            at,
            DiagKind::Shape,
            format_args!("expected list(re, im, rad) of length 3, got {}", v.len()),
        );
    }
    let (re, im, rad) = (*v.first()?, *v.get(1)?, *v.get(2)?);
    if !re.is_finite() || !im.is_finite() || !rad.is_finite() || rad < 0.0 {
        return bad(
            at,
            DiagKind::Domain,
            format_args!("complex ball components must be finite with rad >= 0"),
        );
    }
    ComplexBall::new(Complex64::new(re, im), rad)
}

fn encode_ball(b: ComplexBall) -> Option<Data> {
    if !b.is_finite() {
        return bad(
            0,
            DiagKind::Domain,
            format_args!("complex ball result overflowed to non-finite bounds"),
        );
    }
    Some(Data::vector(vec![b.mid.re, b.mid.im, b.rad]))
}

fn decode_cball(
    args: &[Data],
    at: usize,
) -> Option<Cball> {
    let v = vector(args, at)?;
    if v.len() != 4 {
        return bad(
            at,
            DiagKind::Shape,
            format_args!(
                "expected list(re_lo, re_hi, im_lo, im_hi) of length 4, got {}",
                v.len()
            ),
        );
    }
    let (rlo, rhi, ilo, ihi) = (*v.first()?, *v.get(1)?, *v.get(2)?, *v.get(3)?);
    if !rlo.is_finite()
        || !rhi.is_finite()
        || !ilo.is_finite()
        || !ihi.is_finite()
        || rlo > rhi
        || ilo > ihi
    {
        return bad(
            at,
            DiagKind::Domain,
            format_args!("complex interval bounds must be finite with lo <= hi"),
        );
    }
    Some(Cball::rectangle((rlo, rhi), (ilo, ihi)))
}

fn encode_cball(b: Cball) -> Option<Data> {
    if !b.is_finite() {
        return bad(
            0,
            DiagKind::Domain,
            format_args!("complex interval result overflowed to non-finite bounds"),
        );
    }
    Some(Data::vector(vec![b.re.lo, b.re.hi, b.im.lo, b.im.hi]))
}

fn cball_add(args: &[Data]) -> Option<Data> {
    encode_ball(decode_ball(args, 0)?.add_ball(&decode_ball(args, 1)?))
}

fn cball_sub(args: &[Data]) -> Option<Data> {
    encode_ball(decode_ball(args, 0)?.sub_ball(&decode_ball(args, 1)?))
}

fn cball_mul(args: &[Data]) -> Option<Data> {
    encode_ball(decode_ball(args, 0)?.mul_ball(&decode_ball(args, 1)?))
}

fn cball_div(args: &[Data]) -> Option<Data> {
    let a = decode_ball(args, 0)?;
    let b = decode_ball(args, 1)?;
    let Some(res) = a.div_ball(&b) else {
        return bad(
            1,
            DiagKind::Domain,
            format_args!("complex ball divisor cannot be separated from zero"),
        );
    };
    encode_ball(res)
}

fn cball_pow(args: &[Data]) -> Option<Data> {
    let a = decode_ball(args, 0)?;
    let b = decode_ball(args, 1)?;
    let Some(res) = a.powc(&b) else {
        return bad(
            0,
            DiagKind::Domain,
            format_args!("complex ball base meets the logarithm branch cut"),
        );
    };
    encode_ball(res)
}

fn cball_exp(args: &[Data]) -> Option<Data> {
    encode_ball(decode_ball(args, 0)?.exp())
}

fn cball_ln(args: &[Data]) -> Option<Data> {
    let a = decode_ball(args, 0)?;
    let Some(res) = a.ln() else {
        return bad(
            0,
            DiagKind::Domain,
            format_args!("complex ball meets the logarithm branch cut (-inf, 0]"),
        );
    };
    encode_ball(res)
}

fn cball_sqrt(args: &[Data]) -> Option<Data> {
    let a = decode_ball(args, 0)?;
    let Some(res) = a.sqrt() else {
        return bad(
            0,
            DiagKind::Domain,
            format_args!("complex ball crosses the square-root branch cut (-inf, 0)"),
        );
    };
    encode_ball(res)
}

fn cball_sin(args: &[Data]) -> Option<Data> {
    encode_ball(decode_ball(args, 0)?.sin())
}

fn cball_cos(args: &[Data]) -> Option<Data> {
    encode_ball(decode_ball(args, 0)?.cos())
}

fn cball_tan(args: &[Data]) -> Option<Data> {
    let a = decode_ball(args, 0)?;
    let Some(res) = a.tan() else {
        return bad(
            0,
            DiagKind::Domain,
            format_args!("complex ball tangent denominator may contain zero"),
        );
    };
    encode_ball(res)
}

fn cball_sinh(args: &[Data]) -> Option<Data> {
    encode_ball(decode_ball(args, 0)?.sinh())
}

fn cball_cosh(args: &[Data]) -> Option<Data> {
    encode_ball(decode_ball(args, 0)?.cosh())
}

fn cball_tanh(args: &[Data]) -> Option<Data> {
    let a = decode_ball(args, 0)?;
    let Some(res) = a.tanh() else {
        return bad(
            0,
            DiagKind::Domain,
            format_args!("complex ball tanh denominator may contain zero"),
        );
    };
    encode_ball(res)
}

fn cinterval_add(args: &[Data]) -> Option<Data> {
    encode_cball(decode_cball(args, 0)? + decode_cball(args, 1)?)
}

fn cinterval_sub(args: &[Data]) -> Option<Data> {
    encode_cball(decode_cball(args, 0)? - decode_cball(args, 1)?)
}

fn cinterval_mul(args: &[Data]) -> Option<Data> {
    encode_cball(decode_cball(args, 0)? * decode_cball(args, 1)?)
}

fn cinterval_div(args: &[Data]) -> Option<Data> {
    let a = decode_cball(args, 0)?;
    let b = decode_cball(args, 1)?;
    let Some(res) = a.div(&b) else {
        return bad(
            1,
            DiagKind::Domain,
            format_args!("complex interval divisor cannot be separated from zero"),
        );
    };
    encode_cball(res)
}

fn cinterval_exp(args: &[Data]) -> Option<Data> {
    encode_cball(decode_cball(args, 0)?.exp())
}

fn cinterval_ln(args: &[Data]) -> Option<Data> {
    let a = decode_cball(args, 0)?;
    let Some(res) = a.ln() else {
        return bad(
            0,
            DiagKind::Domain,
            format_args!("complex interval meets the logarithm branch cut (-inf, 0]"),
        );
    };
    encode_cball(res)
}

fn cinterval_sqrt(args: &[Data]) -> Option<Data> {
    let a = decode_cball(args, 0)?;
    let Some(res) = a.sqrt() else {
        return bad(
            0,
            DiagKind::Domain,
            format_args!("complex interval crosses the square-root branch cut (-inf, 0)"),
        );
    };
    encode_cball(res)
}

fn cinterval_sin(args: &[Data]) -> Option<Data> {
    encode_cball(decode_cball(args, 0)?.sin())
}

fn cinterval_cos(args: &[Data]) -> Option<Data> {
    encode_cball(decode_cball(args, 0)?.cos())
}

fn cball_asin(args: &[Data]) -> Option<Data> {
    let a = decode_ball(args, 0)?;
    let Some(res) = a.asin() else {
        return bad(
            0,
            DiagKind::Domain,
            format_args!("complex ball meets asin branch cut"),
        );
    };
    encode_ball(res)
}

fn cball_acos(args: &[Data]) -> Option<Data> {
    let a = decode_ball(args, 0)?;
    let Some(res) = a.acos() else {
        return bad(
            0,
            DiagKind::Domain,
            format_args!("complex ball meets acos branch cut"),
        );
    };
    encode_ball(res)
}

fn cball_atan(args: &[Data]) -> Option<Data> {
    let a = decode_ball(args, 0)?;
    let Some(res) = a.atan() else {
        return bad(
            0,
            DiagKind::Domain,
            format_args!("complex ball meets atan branch cut"),
        );
    };
    encode_ball(res)
}

fn cball_asinh(args: &[Data]) -> Option<Data> {
    let a = decode_ball(args, 0)?;
    let Some(res) = a.asinh() else {
        return bad(
            0,
            DiagKind::Domain,
            format_args!("complex ball meets asinh branch cut"),
        );
    };
    encode_ball(res)
}

fn cball_acosh(args: &[Data]) -> Option<Data> {
    let a = decode_ball(args, 0)?;
    let Some(res) = a.acosh() else {
        return bad(
            0,
            DiagKind::Domain,
            format_args!("complex ball meets acosh branch cut"),
        );
    };
    encode_ball(res)
}

fn cball_atanh(args: &[Data]) -> Option<Data> {
    let a = decode_ball(args, 0)?;
    let Some(res) = a.atanh() else {
        return bad(
            0,
            DiagKind::Domain,
            format_args!("complex ball meets atanh branch cut"),
        );
    };
    encode_ball(res)
}

fn cball_conj(args: &[Data]) -> Option<Data> {
    encode_ball(decode_ball(args, 0)?.conj())
}

fn cball_powi(args: &[Data]) -> Option<Data> {
    let a = decode_ball(args, 0)?;
    let n_f = super::bridge::scalar(args, 1)?;
    if !n_f.is_finite()
        || n_f.fract() != 0.0
        || n_f < f64::from(i32::MIN)
        || n_f > f64::from(i32::MAX)
    {
        return bad(
            1,
            DiagKind::Domain,
            format_args!("cball_powi requires an integer exponent"),
        );
    }
    let Some(res) = a.powi(n_f as i32) else {
        return bad(
            0,
            DiagKind::Domain,
            format_args!("complex ball powi failed (zero base with negative exponent)"),
        );
    };
    encode_ball(res)
}

fn cball_abs(args: &[Data]) -> Option<Data> {
    let a = decode_ball(args, 0)?;
    let Some(res) = a.abs_ball() else {
        return bad(0, DiagKind::Domain, format_args!("non-finite complex ball"));
    };
    encode_ball(res)
}

fn cball_arg(args: &[Data]) -> Option<Data> {
    let a = decode_ball(args, 0)?;
    let Some(res) = a.arg_ball() else {
        return bad(
            0,
            DiagKind::Domain,
            format_args!("complex ball meets arg branch cut (-inf, 0]"),
        );
    };
    encode_ball(res)
}

fn cball_contains(args: &[Data]) -> Option<Data> {
    let a = decode_ball(args, 0)?;
    let b = decode_ball(args, 1)?;
    Some(Data::Scalar(
        if b.inside_strictly(&a) || (b.rad == 0.0 && a.contains(b.mid)) {
            1.0
        } else {
            0.0
        },
    ))
}

fn cball_overlaps(args: &[Data]) -> Option<Data> {
    let a = decode_ball(args, 0)?;
    let b = decode_ball(args, 1)?;
    Some(Data::Scalar(if a.overlaps(&b) { 1.0 } else { 0.0 }))
}

fn cball_to_cinterval(args: &[Data]) -> Option<Data> {
    encode_cball(decode_ball(args, 0)?.to_cball())
}

fn cinterval_to_cball(args: &[Data]) -> Option<Data> {
    encode_ball(ComplexBall::from_cball(decode_cball(args, 0)?))
}

fn cinterval_tan(args: &[Data]) -> Option<Data> {
    let a = decode_cball(args, 0)?;
    let Some(res) = a.tan() else {
        return bad(
            0,
            DiagKind::Domain,
            format_args!("complex interval tangent denominator may contain zero"),
        );
    };
    encode_cball(res)
}

fn cinterval_sinh(args: &[Data]) -> Option<Data> {
    encode_cball(decode_cball(args, 0)?.sinh())
}

fn cinterval_cosh(args: &[Data]) -> Option<Data> {
    encode_cball(decode_cball(args, 0)?.cosh())
}

fn cinterval_tanh(args: &[Data]) -> Option<Data> {
    let a = decode_cball(args, 0)?;
    let Some(res) = a.tanh() else {
        return bad(
            0,
            DiagKind::Domain,
            format_args!("complex interval tanh denominator may contain zero"),
        );
    };
    encode_cball(res)
}

fn cinterval_pow(args: &[Data]) -> Option<Data> {
    let a = decode_cball(args, 0)?;
    let b = decode_cball(args, 1)?;
    let Some(res) = a.powc(&b) else {
        return bad(
            0,
            DiagKind::Domain,
            format_args!("complex interval base meets logarithm branch cut"),
        );
    };
    encode_cball(res)
}

fn cinterval_asin(args: &[Data]) -> Option<Data> {
    let a = decode_cball(args, 0)?;
    let Some(res) = a.asin() else {
        return bad(
            0,
            DiagKind::Domain,
            format_args!("complex interval meets asin branch cut"),
        );
    };
    encode_cball(res)
}

fn cinterval_acos(args: &[Data]) -> Option<Data> {
    let a = decode_cball(args, 0)?;
    let Some(res) = a.acos() else {
        return bad(
            0,
            DiagKind::Domain,
            format_args!("complex interval meets acos branch cut"),
        );
    };
    encode_cball(res)
}

fn cinterval_atan(args: &[Data]) -> Option<Data> {
    let a = decode_cball(args, 0)?;
    let Some(res) = a.atan() else {
        return bad(
            0,
            DiagKind::Domain,
            format_args!("complex interval meets atan branch cut"),
        );
    };
    encode_cball(res)
}

fn inverse_hyperbolic_interval(
    args: &[Data],
    name: &str,
    enclose: fn(&Cball) -> Option<Cball>,
) -> Option<Data> {
    let interval = decode_cball(args, 0)?;
    let Some(result) = enclose(&interval) else {
        return bad(
            0,
            DiagKind::Domain,
            format_args!("cannot enclose principal {name} on this complex interval"),
        );
    };
    encode_cball(result)
}

fn cinterval_asinh(args: &[Data]) -> Option<Data> {
    inverse_hyperbolic_interval(args, "asinh", Cball::asinh)
}

fn cinterval_acosh(args: &[Data]) -> Option<Data> {
    inverse_hyperbolic_interval(args, "acosh", Cball::acosh)
}

fn cinterval_atanh(args: &[Data]) -> Option<Data> {
    inverse_hyperbolic_interval(args, "atanh", Cball::atanh)
}

fn cinterval_intersect(args: &[Data]) -> Option<Data> {
    let a = decode_cball(args, 0)?;
    let b = decode_cball(args, 1)?;
    let Some(res) = a.intersect(&b) else {
        return bad(
            0,
            DiagKind::Domain,
            format_args!("complex intervals are disjoint"),
        );
    };
    encode_cball(res)
}
