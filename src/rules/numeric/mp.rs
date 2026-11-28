//! Arbitrary-precision operators.
//!
//! | operator | value |
//! |---|---|
//! | `mp_eval(expr, digits)` | the value of `expr` to `digits` significant decimal digits, as text |
//! | `mp_const(name, digits)` | a constant (`pi`, `e`, `tau`, `ln2`, `sqrt2`, `phi`) to `digits` digits |
//! | `mp_round(x, digits)` | `x` rounded to `digits` significant digits, nearest even |
//! | `mp_round_mode(x, digits, mode)` | the same with `nearest`, `zero`, `away`, `up` or `down` |
//! | `mp_add(a, b, digits)` | `a + b` (likewise `mp_sub`, `mp_mul`, `mp_div`) |
//! | `mp_compare(a, b)` | `-1`, `0` or `1` as `a` is below, equal to or above `b`, exactly |
//! | `mp_to_float(x)` | `x` narrowed to a double, correctly rounded |
//! | `mp_from_float(x, digits)` | the exact decimal expansion of the double `x` |
//! | `mp_ball(expr, digits)` | `list(low, high)`: a rigorous enclosure of `expr` |
//! | `mp_digits(x)` | how many significant decimal digits `x` carries |
//!
//! A number crosses this boundary as **decimal text**, which is the design
//! point: the digits are part of the value, so precision is hash-consed,
//! cached and shared like any other literal, and there is no session-wide
//! rounding mode to reason about. `mp_eval` accepts the small language of
//! `crate::kernels::mp::parse` (`+ - * / ^`, parentheses, the constants and
//! `exp ln log log2 log10 sqrt cbrt abs sin cos tan asin acos atan atan2 sinh
//! cosh tanh exp2 hypot pow`).
//!
//! A result is returned only when the enclosure pins the requested digits, so
//! `mp_eval("pi", 50)` is 50 correct digits rather than 50 plausible ones.
//! `mp_ball` is the same computation without that requirement: it hands back
//! the enclosure, for a host that wants to see the uncertainty itself.

use super::bridge::Data;
use super::bridge::bad;
use super::bridge::count;
use super::bridge::name;
use super::bridge::num;
use super::bridge::text;
use super::bridge::text_opt;
use crate::kernel::diag::DiagKind;
use crate::kernels::mp;

reviewed_entries! {
    (
        "mp_eval", "td", mp_eval, mp::eval_text,
        "`mp_eval(expr, digits)`: the value of the expression `expr` to `digits` significant decimal digits, as text, correctly rounded."
    ),
    (
        "mp_const", "nd", mp_const, mp::const_text,
        "`mp_const(name, digits)`: the constant `pi`, `e`, `tau`, `ln2`, `sqrt2` or `phi` to `digits` significant decimal digits."
    ),
    (
        "mp_round", "td", mp_round, mp::round_text,
        "`mp_round(x, digits)`: the decimal `x` rounded to `digits` significant digits, nearest even."
    ),
    (
        "mp_round_mode", "tdn", mp_round_mode, mp::round_text,
        "`mp_round_mode(x, digits, mode)`: `x` rounded to `digits` significant digits in `mode` (`nearest`, `zero`, `away`, `up`, `down`)."
    ),
    (
        "mp_add", "ttd", mp_add, mp::arithmetic_text,
        "`mp_add(a, b, digits)`: the sum of two decimal numbers, to `digits` significant digits."
    ),
    (
        "mp_sub", "ttd", mp_sub, mp::arithmetic_text,
        "`mp_sub(a, b, digits)`: the difference of two decimal numbers, to `digits` significant digits."
    ),
    (
        "mp_mul", "ttd", mp_mul, mp::arithmetic_text,
        "`mp_mul(a, b, digits)`: the product of two decimal numbers, to `digits` significant digits."
    ),
    (
        "mp_div", "ttd", mp_div, mp::arithmetic_text,
        "`mp_div(a, b, digits)`: the quotient of two decimal numbers, to `digits` significant digits."
    ),
    (
        "mp_compare", "tt", mp_compare, mp::compare_text,
        "`mp_compare(a, b)`: `-1`, `0` or `1` as the decimal `a` is below, equal to or above `b`, compared exactly."
    ),
    (
        "mp_to_float", "t", mp_to_float, mp::to_float_text,
        "`mp_to_float(x)`: the decimal `x` narrowed to a double, correctly rounded (no double rounding through an intermediate)."
    ),
    (
        "mp_from_float", "sd", mp_from_float, mp::from_float_text,
        "`mp_from_float(x, digits)`: the decimal expansion of the double `x` to `digits` significant digits; the double is exact, so the digits are its true value."
    ),
    (
        "mp_ball", "td", mp_ball, mp::enclosure_text,
        "`mp_ball(expr, digits)`: `list(low, high)`, a rigorous enclosure of `expr` at `digits` significant digits."
    ),
    (
        "mp_complex_ball", "td", mp_complex_ball, mp::enclosure_complex_text,
        "`mp_complex_ball(expr, digits)`: `list(list(re_lo, re_hi), list(im_lo, im_hi))`, a rigorous complex enclosure of `expr` at `digits` significant digits."
    ),
    (
        "mp_complex_eval", "td", mp_complex_eval, mp::eval_complex_text,
        "`mp_complex_eval(expr, digits)`: `list(re, im)` of the complex expression `expr` to `digits` significant decimal digits."
    ),
    (
        "mp_digits", "t", mp_digits, mp::digits_of,
        "`mp_digits(x)`: the number of significant decimal digits the literal `x` carries."
    ),
}

/// Maps a kernel failure onto a recorded reason.
fn reason<T>(
    at: usize,
    result: Result<T, String>,
) -> Option<T> {
    match result {
        | Ok(value) => Some(value),
        | Err(message) => bad(at, DiagKind::Domain, format_args!("{message}")),
    }
}

/// The digit count argument.
fn digits(
    args: &[Data],
    at: usize,
) -> Option<usize> {
    let value = count(args, at, mp::MAX_DIGITS)?;
    if value == 0 {
        return bad(
            at,
            DiagKind::Domain,
            format_args!("at least one significant digit is required"),
        );
    }
    Some(value)
}

fn mp_eval(args: &[Data]) -> Option<Data> {
    let expression = text(args, 0)?;
    let places = digits(args, 1)?;
    let rendered = reason(0, mp::eval_text(expression, places))?;
    super::bridge::text_out(rendered)
}

fn name_or_text(
    args: &[Data],
    at: usize,
) -> Option<&str> {
    text_opt(args, at).map_or_else(|| name(args, at), Some)
}

fn mp_const(args: &[Data]) -> Option<Data> {
    let symbol = name_or_text(args, 0)?;
    let places = digits(args, 1)?;
    let rendered = reason(0, mp::const_text(symbol, places))?;
    super::bridge::text_out(rendered)
}

fn mp_round(args: &[Data]) -> Option<Data> {
    let value = text(args, 0)?;
    let places = digits(args, 1)?;
    let rendered = reason(0, mp::round_text(value, places, mp::Round::NearestEven))?;
    super::bridge::text_out(rendered)
}

fn mp_round_mode(args: &[Data]) -> Option<Data> {
    let value = text(args, 0)?;
    let places = digits(args, 1)?;
    let mode = name_or_text(args, 2)?;
    let rounding = reason(2, mp::round_named(mode))?;
    let rendered = reason(0, mp::round_text(value, places, rounding))?;
    super::bridge::text_out(rendered)
}

fn arithmetic(
    operation: &str,
    args: &[Data],
) -> Option<Data> {
    let left = text(args, 0)?;
    let right = text(args, 1)?;
    let places = digits(args, 2)?;
    let rendered = reason(
        0,
        mp::arithmetic_text(operation, left, right, places, mp::Round::NearestEven),
    )?;
    super::bridge::text_out(rendered)
}

fn mp_add(args: &[Data]) -> Option<Data> {
    arithmetic("add", args)
}

fn mp_sub(args: &[Data]) -> Option<Data> {
    arithmetic("sub", args)
}

fn mp_mul(args: &[Data]) -> Option<Data> {
    arithmetic("mul", args)
}

fn mp_div(args: &[Data]) -> Option<Data> {
    arithmetic("div", args)
}

fn mp_compare(args: &[Data]) -> Option<Data> {
    let left = text(args, 0)?;
    let right = text(args, 1)?;
    let ordering = reason(0, mp::compare_text(left, right))?;
    num(f64::from(ordering))
}

fn mp_to_float(args: &[Data]) -> Option<Data> {
    let value = text(args, 0)?;
    let narrowed = reason(0, mp::to_float_text(value))?;
    if narrowed.is_finite() {
        return num(narrowed);
    }
    bad(
        0,
        DiagKind::Domain,
        format_args!("the value is outside the range of a double"),
    )
}

fn mp_from_float(args: &[Data]) -> Option<Data> {
    let value = super::bridge::scalar(args, 0)?;
    let places = digits(args, 1)?;
    let rendered = reason(0, mp::from_float_text(value, places))?;
    super::bridge::text_out(rendered)
}

fn mp_ball(args: &[Data]) -> Option<Data> {
    let expression = text(args, 0)?;
    let places = digits(args, 1)?;
    let (low, high) = reason(0, mp::enclosure_text(expression, places))?;
    Some(Data::Tree(vec![Data::Text(low), Data::Text(high)]))
}

fn mp_complex_eval(args: &[Data]) -> Option<Data> {
    let expression = text(args, 0)?;
    let places = digits(args, 1)?;
    let (re, im) = reason(0, mp::eval_complex_text(expression, places))?;
    Some(Data::Tree(vec![Data::Text(re), Data::Text(im)]))
}

fn mp_complex_ball(args: &[Data]) -> Option<Data> {
    let expression = text(args, 0)?;
    let places = digits(args, 1)?;
    let ((re_lo, re_hi), (im_lo, im_hi)) =
        reason(0, mp::enclosure_complex_text(expression, places))?;
    Some(Data::Tree(vec![
        Data::Tree(vec![Data::Text(re_lo), Data::Text(re_hi)]),
        Data::Tree(vec![Data::Text(im_lo), Data::Text(im_hi)]),
    ]))
}

fn mp_digits(args: &[Data]) -> Option<Data> {
    let value = text(args, 0)?;
    let places = reason(0, mp::digits_of(value))?;
    #[allow(clippy::cast_precision_loss)]
    let scalar = places as f64;
    num(scalar)
}

#[cfg(test)]
mod tests {
    use crate::api::Session;
    use crate::kernel::CallCtx;
    use crate::kernel::value::Value;

    fn call(
        name: &str,
        args: &[Value],
    ) -> Value {
        let session = Session::new();
        session
            .call_value(name, args, &CallCtx::new())
            .unwrap_or_else(|error| panic!("{name}: {error}"))
    }

    fn text(
        name: &str,
        args: &[Value],
    ) -> String {
        match call(name, args) {
            | Value::String(text) => text,
            | other => panic!("{name} returned {other:?}"),
        }
    }

    fn number(
        name: &str,
        args: &[Value],
    ) -> f64 {
        match call(name, args) {
            | Value::Float(value) => value,
            | Value::Real(number) => number.to_f64(),
            | Value::Int(value) => value.to_string().parse().unwrap_or(f64::NAN),
            | other => panic!("{name} returned {other:?}"),
        }
    }

    fn literal(text: &str) -> Value {
        Value::String(text.to_string())
    }

    #[test]
    fn named_options_accept_text_or_symbol_at_the_value_boundary() {
        for upward in [literal("up"), Value::Symbol("up".to_owned())] {
            assert_eq!(
                text("mp_round_mode", &[literal("1.25"), Value::int(2), upward]),
                "1.3"
            );
        }
        for pi in [literal("pi"), Value::Symbol("pi".to_owned())] {
            assert_eq!(text("mp_const", &[pi, Value::int(4)]), "3.142");
        }
    }

    #[test]
    fn precision_operators_reach_the_engine() {
        assert_eq!(
            text("mp_eval", &[literal("pi"), Value::int(20)]),
            "3.1415926535897932385"
        );
        assert_eq!(
            text("mp_eval", &[literal("1/3"), Value::int(20)]),
            "0.33333333333333333333"
        );
        assert_eq!(
            text("mp_const", &[literal("e"), Value::int(15)]),
            "2.71828182845905"
        );
        assert_eq!(
            text("mp_add", &[literal("0.1"), literal("0.2"), Value::int(30)]),
            "0.3"
        );
        assert_eq!(
            number(
                "mp_compare",
                &[literal("1.0000000000000000001"), literal("1")]
            ),
            1.0
        );
        assert_eq!(number("mp_to_float", &[literal("0.1")]), 0.1);
        assert_eq!(
            text("mp_from_float", &[Value::Float(0.1), Value::int(30)]),
            "0.100000000000000005551115123126"
        );
        assert_eq!(text("mp_round", &[literal("2.5"), Value::int(1)]), "2");
        assert_eq!(
            text(
                "mp_round_mode",
                &[
                    literal("1.25"),
                    Value::int(2),
                    Value::Symbol("up".to_string())
                ]
            ),
            "1.3"
        );
    }
}
