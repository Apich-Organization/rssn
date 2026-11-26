//! Complex ball (`list(re, im, rad)`) and rectangular complex interval
//! (`list(re_lo, re_hi, im_lo, im_hi)`) operators.

use num_complex::Complex64;

use super::bridge::Data;
use super::bridge::Entry;
use super::bridge::bad;
use super::bridge::vector;
use crate::kernel::diag::DiagKind;
use crate::kernels::ball::Cball;
use crate::kernels::complex_ball::ComplexBall;

/// Every operator of this family.
pub(super) const ENTRIES: &[Entry] = &[
    Entry {
        name: "cball_add",
        params: "vv",
        eval: cball_add,
        doc: "`cball_add(a, b)`: sum of two complex balls `list(re, im, rad)`.",
    },
    Entry {
        name: "cball_sub",
        params: "vv",
        eval: cball_sub,
        doc: "`cball_sub(a, b)`: difference of two complex balls `list(re, im, rad)`.",
    },
    Entry {
        name: "cball_mul",
        params: "vv",
        eval: cball_mul,
        doc: "`cball_mul(a, b)`: product of two complex balls `list(re, im, rad)`.",
    },
    Entry {
        name: "cball_div",
        params: "vv",
        eval: cball_div,
        doc: "`cball_div(a, b)`: quotient of two complex balls `list(re, im, rad)` when `b` excludes zero.",
    },
    Entry {
        name: "cball_pow",
        params: "vv",
        eval: cball_pow,
        doc: "`cball_pow(a, b)`: principal complex power `a^b` over complex balls `list(re, im, rad)`.",
    },
    Entry {
        name: "cball_exp",
        params: "v",
        eval: cball_exp,
        doc: "`cball_exp(z)`: complex exponential enclosure of `list(re, im, rad)`.",
    },
    Entry {
        name: "cball_ln",
        params: "v",
        eval: cball_ln,
        doc: "`cball_ln(z)`: principal complex logarithm enclosure of `list(re, im, rad)`.",
    },
    Entry {
        name: "cball_sqrt",
        params: "v",
        eval: cball_sqrt,
        doc: "`cball_sqrt(z)`: principal complex square root enclosure of `list(re, im, rad)`.",
    },
    Entry {
        name: "cball_sin",
        params: "v",
        eval: cball_sin,
        doc: "`cball_sin(z)`: complex sine enclosure of `list(re, im, rad)`.",
    },
    Entry {
        name: "cball_cos",
        params: "v",
        eval: cball_cos,
        doc: "`cball_cos(z)`: complex cosine enclosure of `list(re, im, rad)`.",
    },
    Entry {
        name: "cball_tan",
        params: "v",
        eval: cball_tan,
        doc: "`cball_tan(z)`: complex tangent enclosure of `list(re, im, rad)`.",
    },
    Entry {
        name: "cball_sinh",
        params: "v",
        eval: cball_sinh,
        doc: "`cball_sinh(z)`: complex hyperbolic sine enclosure of `list(re, im, rad)`.",
    },
    Entry {
        name: "cball_cosh",
        params: "v",
        eval: cball_cosh,
        doc: "`cball_cosh(z)`: complex hyperbolic cosine enclosure of `list(re, im, rad)`.",
    },
    Entry {
        name: "cball_tanh",
        params: "v",
        eval: cball_tanh,
        doc: "`cball_tanh(z)`: complex hyperbolic tangent enclosure of `list(re, im, rad)`.",
    },
    Entry {
        name: "cinterval_add",
        params: "vv",
        eval: cinterval_add,
        doc: "`cinterval_add(a, b)`: sum of two complex intervals `list(re_lo, re_hi, im_lo, im_hi)`.",
    },
    Entry {
        name: "cinterval_sub",
        params: "vv",
        eval: cinterval_sub,
        doc: "`cinterval_sub(a, b)`: difference of two complex intervals `list(re_lo, re_hi, im_lo, im_hi)`.",
    },
    Entry {
        name: "cinterval_mul",
        params: "vv",
        eval: cinterval_mul,
        doc: "`cinterval_mul(a, b)`: product of two complex intervals `list(re_lo, re_hi, im_lo, im_hi)`.",
    },
    Entry {
        name: "cinterval_div",
        params: "vv",
        eval: cinterval_div,
        doc: "`cinterval_div(a, b)`: quotient of two complex intervals `list(re_lo, re_hi, im_lo, im_hi)`.",
    },
    Entry {
        name: "cinterval_exp",
        params: "v",
        eval: cinterval_exp,
        doc: "`cinterval_exp(z)`: complex exponential enclosure of `list(re_lo, re_hi, im_lo, im_hi)`.",
    },
    Entry {
        name: "cinterval_ln",
        params: "v",
        eval: cinterval_ln,
        doc: "`cinterval_ln(z)`: principal complex logarithm enclosure of `list(re_lo, re_hi, im_lo, im_hi)`.",
    },
    Entry {
        name: "cinterval_sqrt",
        params: "v",
        eval: cinterval_sqrt,
        doc: "`cinterval_sqrt(z)`: principal complex square root enclosure of `list(re_lo, re_hi, im_lo, im_hi)`.",
    },
    Entry {
        name: "cinterval_sin",
        params: "v",
        eval: cinterval_sin,
        doc: "`cinterval_sin(z)`: complex sine enclosure of `list(re_lo, re_hi, im_lo, im_hi)`.",
    },
    Entry {
        name: "cinterval_cos",
        params: "v",
        eval: cinterval_cos,
        doc: "`cinterval_cos(z)`: complex cosine enclosure of `list(re_lo, re_hi, im_lo, im_hi)`.",
    },
];

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
