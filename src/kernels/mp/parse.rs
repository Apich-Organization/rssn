//! A small expression language over the arbitrary-precision tower.
//!
//! `mp_eval` takes one string and returns one string, so the grammar is
//! deliberately tiny and total: numbers, the constants of
//! [`super::consts`], `+ - * / ^`, parentheses, and the functions below. It
//! exists so a host can ask for `sin(1/3) + pi` at 200 digits without a term
//! per operation; it is not a general language and it does not call back into
//! the rewrite engine.
//!
//! ```text
//! expr    := term (('+' | '-') term)*
//! term    := power (('*' | '/') power)*
//! power   := unary ('^' power)?          // right associative
//! unary   := ('-' | '+') unary | primary
//! primary := number | name | name '(' expr (',' expr)* ')' | '(' expr ')'
//! ```

use super::ball::Ball;
use super::consts;
use super::decimal;
use super::funcs;
use super::real::Real;
use super::real::Round;

/// Longest expression text accepted.
const MAX_TEXT: usize = 4096;
/// Deepest parenthesised nesting accepted.
const MAX_DEPTH: usize = 64;
/// Most tokens accepted.
const MAX_TOKENS: usize = 2048;

/// One lexical token.
#[derive(Clone, Debug, PartialEq)]
enum Token {
    /// A decimal literal, as written.
    Number(String),
    /// An identifier: a constant or a function name.
    Name(String),
    /// `+`
    Plus,
    /// `-`
    Minus,
    /// `*`
    Star,
    /// `/`
    Slash,
    /// `^`
    Caret,
    /// `(`
    Open,
    /// `)`
    Close,
    /// `,`
    Comma,
}

/// Splits `text` into tokens.
fn tokenize(text: &str) -> Result<Vec<Token>, String> {
    let mut tokens = Vec::new();
    let characters: Vec<char> = text.chars().collect();
    let mut at = 0;
    while at < characters.len() {
        let character = characters[at];
        match character {
            | ' ' | '\t' | '\n' | '\r' => at += 1,
            | '+' => {
                tokens.push(Token::Plus);
                at += 1;
            },
            | '-' => {
                tokens.push(Token::Minus);
                at += 1;
            },
            | '*' => {
                tokens.push(Token::Star);
                at += 1;
            },
            | '/' => {
                tokens.push(Token::Slash);
                at += 1;
            },
            | '^' => {
                tokens.push(Token::Caret);
                at += 1;
            },
            | '(' => {
                tokens.push(Token::Open);
                at += 1;
            },
            | ')' => {
                tokens.push(Token::Close);
                at += 1;
            },
            | ',' => {
                tokens.push(Token::Comma);
                at += 1;
            },
            | '0'..='9' | '.' => {
                let start = at;
                while at < characters.len()
                    && (characters[at].is_ascii_digit() || characters[at] == '.')
                {
                    at += 1;
                }
                if at < characters.len() && (characters[at] == 'e' || characters[at] == 'E') {
                    let save = at;
                    at += 1;
                    if at < characters.len() && (characters[at] == '+' || characters[at] == '-') {
                        at += 1;
                    }
                    let digits_start = at;
                    while at < characters.len() && characters[at].is_ascii_digit() {
                        at += 1;
                    }
                    if at == digits_start {
                        at = save;
                    }
                }
                tokens.push(Token::Number(characters[start..at].iter().collect()));
            },
            | _ if character.is_alphabetic() || character == '_' => {
                let start = at;
                while at < characters.len()
                    && (characters[at].is_alphanumeric() || characters[at] == '_')
                {
                    at += 1;
                }
                tokens.push(Token::Name(characters[start..at].iter().collect()));
            },
            | _ => return Err(format!("`{character}` is not allowed in an expression")),
        }
        if tokens.len() > MAX_TOKENS {
            return Err("the expression has too many tokens".to_string());
        }
    }
    if tokens.is_empty() {
        return Err("the expression is empty".to_string());
    }
    Ok(tokens)
}

fn expression_tokens(text: &str) -> Result<Vec<Token>, String> {
    if text.len() > MAX_TEXT {
        return Err("the expression is too long".to_string());
    }
    tokenize(text)
}

/// Handles plain decimal expressions without losing exact rounding ties.
pub(super) fn round_literal(
    text: &str,
    digits: usize,
) -> Result<Option<String>, String> {
    match expression_tokens(text)?.as_slice() {
        | [Token::Number(literal)] | [Token::Plus, Token::Number(literal)] => {
            decimal::round_decimal(literal, digits, Round::NearestEven).map(Some)
        },
        | [Token::Minus, Token::Number(literal)] => {
            decimal::round_decimal(&format!("-{literal}"), digits, Round::NearestEven).map(Some)
        },
        | _ => Ok(None),
    }
}

/// A cursor over the tokens.
struct Parser<'a> {
    tokens: &'a [Token],
    at: usize,
    depth: usize,
    prec: u32,
}

impl Parser<'_> {
    /// The next token, without consuming it.
    fn peek(&self) -> Option<&Token> {
        self.tokens.get(self.at)
    }

    /// Consumes the next token when it is `wanted`.
    fn eat(
        &mut self,
        wanted: &Token,
    ) -> bool {
        if self.peek() == Some(wanted) {
            self.at += 1;
            return true;
        }
        false
    }

    /// Parses a sum or difference.
    fn expression(&mut self) -> Result<Ball, String> {
        let mut value = self.term()?;
        loop {
            if self.eat(&Token::Plus) {
                value = value.add(&self.term()?, self.prec)?;
            } else if self.eat(&Token::Minus) {
                value = value.sub(&self.term()?, self.prec)?;
            } else {
                return Ok(value);
            }
        }
    }

    /// Parses a product or quotient.
    fn term(&mut self) -> Result<Ball, String> {
        let mut value = self.power()?;
        loop {
            if self.eat(&Token::Star) {
                value = value.mul(&self.power()?, self.prec)?;
            } else if self.eat(&Token::Slash) {
                value = value.div(&self.power()?, self.prec)?;
            } else {
                return Ok(value);
            }
        }
    }

    /// Parses a power, right associatively.
    fn power(&mut self) -> Result<Ball, String> {
        let base = self.unary()?;
        if self.eat(&Token::Caret) {
            let exponent = self.power()?;
            let value = funcs::pow(&base.mid, &exponent.mid, self.prec)?;
            return widen(value, &base, &exponent, self.prec);
        }
        Ok(base)
    }

    /// Parses a unary sign.
    fn unary(&mut self) -> Result<Ball, String> {
        if self.eat(&Token::Minus) {
            return Ok(self.unary()?.neg());
        }
        if self.eat(&Token::Plus) {
            return self.unary();
        }
        self.primary()
    }

    /// Parses a literal, a name, a call or a parenthesised expression.
    fn primary(&mut self) -> Result<Ball, String> {
        let token = self
            .peek()
            .cloned()
            .ok_or_else(|| "the expression ends early".to_string())?;
        match token {
            | Token::Number(text) => {
                self.at += 1;
                decimal::enclose_decimal(&text, self.prec)
            },
            | Token::Name(name) => {
                self.at += 1;
                if self.eat(&Token::Open) {
                    let mut arguments = vec![self.expression()?];
                    while self.eat(&Token::Comma) {
                        arguments.push(self.expression()?);
                    }
                    if !self.eat(&Token::Close) {
                        return Err("missing `)`".to_string());
                    }
                    return call(&name, &arguments, self.prec);
                }
                let value = consts::named(&name, self.prec)?;
                Ok(Ball::point(value))
            },
            | Token::Open => {
                self.at += 1;
                self.depth += 1;
                if self.depth > MAX_DEPTH {
                    return Err("the expression nests too deeply".to_string());
                }
                let value = self.expression()?;
                self.depth -= 1;
                if !self.eat(&Token::Close) {
                    return Err("missing `)`".to_string());
                }
                Ok(value)
            },
            | _ => Err("expected a number, a name or `(`".to_string()),
        }
    }
}

/// Widens `value` by the radii of its operands, which `pow` cannot see.
fn widen(
    value: Ball,
    base: &Ball,
    exponent: &Ball,
    prec: u32,
) -> Result<Ball, String> {
    let inner = prec.saturating_add(16);
    let relative = base
        .rad
        .div(&base.mid.abs(), inner, super::real::Round::AwayFromZero)
        .unwrap_or_else(|_| base.rad.clone());
    let factor = exponent
        .mid
        .abs()
        .mul(&relative, inner, super::real::Round::AwayFromZero)?;
    let bound = funcs::exp(&factor, inner)?;
    let scaled = value.mul(&Ball::point(bound.mid), inner)?;
    let absolute = exponent
        .rad
        .mul(&value.mid.abs(), inner, super::real::Round::AwayFromZero)?;
    let extra = funcs::exp(&absolute, inner)?;
    Ok(Ball {
        mid: scaled
            .mid
            .with_precision(prec, super::real::Round::NearestEven)?,
        rad: scaled.rad.add(
            &extra
                .mid
                .sub(&Real::one(inner), inner, super::real::Round::AwayFromZero)?,
            prec,
            super::real::Round::AwayFromZero,
        )?,
    })
}

/// Applies the named function to `arguments`.
fn call(
    name: &str,
    arguments: &[Ball],
    prec: u32,
) -> Result<Ball, String> {
    let want = |count: usize| {
        if arguments.len() == count {
            Ok(())
        } else {
            Err(format!("`{name}` takes {count} arguments"))
        }
    };
    let first = || {
        arguments
            .first()
            .map(|ball| ball.mid.clone())
            .ok_or_else(|| format!("`{name}` needs an argument"))
    };
    match name {
        | "exp" => {
            want(1)?;
            funcs::map(&arguments[0], &funcs::exp, true, prec)
        },
        | "ln" | "log" => {
            want(1)?;
            if name == "log" {
                return funcs::map(&arguments[0], &funcs::log10, true, prec);
            }
            funcs::map(&arguments[0], &funcs::ln, true, prec)
        },
        | "log2" => {
            want(1)?;
            funcs::map(&arguments[0], &funcs::log2, true, prec)
        },
        | "log10" => {
            want(1)?;
            funcs::map(&arguments[0], &funcs::log10, true, prec)
        },
        | "sqrt" => {
            want(1)?;
            let root = first()?.sqrt(prec, super::real::Round::NearestEven)?;
            Ok(Ball {
                mid: root,
                rad: super::fixed::units(4, prec.saturating_add(8))?,
            })
        },
        | "cbrt" => {
            want(1)?;
            funcs::cbrt(&first()?, prec)
        },
        | "abs" => {
            want(1)?;
            arguments[0].abs(prec)
        },
        | "sin" => {
            want(1)?;
            funcs::map(&arguments[0], &funcs::sin, false, prec)
        },
        | "cos" => {
            want(1)?;
            funcs::map(&arguments[0], &funcs::cos, false, prec)
        },
        | "tan" => {
            want(1)?;
            funcs::map(&arguments[0], &funcs::tan, true, prec)
        },
        | "asin" => {
            want(1)?;
            funcs::map(&arguments[0], &funcs::asin, true, prec)
        },
        | "acos" => {
            want(1)?;
            funcs::map(&arguments[0], &funcs::acos, false, prec)
        },
        | "atan" => {
            want(1)?;
            funcs::map(&arguments[0], &funcs::atan, true, prec)
        },
        | "atan2" => {
            want(2)?;
            funcs::atan2(&arguments[0].mid, &arguments[1].mid, prec)
        },
        | "sinh" => {
            want(1)?;
            funcs::map(&arguments[0], &funcs::sinh, true, prec)
        },
        | "cosh" => {
            want(1)?;
            funcs::map(&arguments[0], &funcs::cosh, false, prec)
        },
        | "tanh" => {
            want(1)?;
            funcs::map(&arguments[0], &funcs::tanh, true, prec)
        },
        | "exp2" => {
            want(1)?;
            funcs::map(&arguments[0], &funcs::exp2, true, prec)
        },
        | "hypot" => {
            want(2)?;
            funcs::hypot(&arguments[0].mid, &arguments[1].mid, prec)
        },
        | "pow" => {
            want(2)?;
            funcs::pow(&arguments[0].mid, &arguments[1].mid, prec)
        },
        | _ => Err(format!("unknown function `{name}`")),
    }
}

/// Evaluates `text` at `prec` bits.
///
/// # Errors
/// Fails on malformed text, an unknown name, a domain error, or when an
/// exponent leaves the representable range.
pub fn evaluate(
    text: &str,
    prec: u32,
) -> Result<Ball, String> {
    let tokens = expression_tokens(text)?;
    let mut parser = Parser {
        tokens: &tokens,
        at: 0,
        depth: 0,
        prec: prec.saturating_add(32),
    };
    let value = parser.expression()?;
    if parser.at != tokens.len() {
        return Err("trailing text after the expression".to_string());
    }
    let mid = value.mid.with_precision(prec, Round::NearestEven)?;
    // A shorter radius can represent smaller exponents. Retain an outward
    // bound even when the displacement is too small at the target precision.
    let displacement = value
        .mid
        .sub(&mid, prec, Round::AwayFromZero)
        .or_else(|_| value.mid.sub(&mid, 1, Round::AwayFromZero))?
        .abs();
    let rad = if displacement.is_zero() {
        value.rad
    } else {
        value
            .rad
            .add(&displacement, prec, Round::TowardPositive)
            .or_else(|_| value.rad.add(&displacement, 1, Round::TowardPositive))?
    };
    Ok(Ball { mid, rad })
}

/// The names [`evaluate`] accepts as functions, for diagnostics and help.
#[must_use]
pub const fn functions() -> &'static [&'static str] {
    &[
        "exp", "ln", "log", "log2", "log10", "sqrt", "cbrt", "abs", "sin", "cos", "tan", "asin",
        "acos", "atan", "atan2", "sinh", "cosh", "tanh", "exp2", "hypot", "pow",
    ]
}

#[cfg(test)]
mod tests {
    use super::super::real::MAX_BITS;
    use super::*;
    use num_bigint::BigInt;
    use num_rational::BigRational;

    fn exact(value: &Real) -> BigRational {
        let shift = usize::try_from(value.exponent().unsigned_abs()).unwrap();
        if value.exponent() < 0 {
            BigRational::new(value.mantissa().clone(), BigInt::from(1) << shift)
        } else {
            BigRational::from_integer(value.mantissa() << shift)
        }
    }

    fn contains_exact(
        ball: &Ball,
        value: &BigRational,
    ) {
        let midpoint = exact(ball.midpoint());
        let radius = exact(ball.radius());
        assert!(
            &midpoint - &radius <= *value,
            "lower bound excludes {value}"
        );
        assert!(midpoint + radius >= *value, "upper bound excludes {value}");
    }

    #[test]
    fn decimal_literals_enclose_exact_rationals() {
        let one_tenth = BigRational::new(BigInt::from(1), BigInt::from(10));
        for precision in [8, MAX_BITS] {
            let positive = evaluate("0.1", precision).unwrap();
            assert_eq!(positive.midpoint().precision(), precision);
            if precision == 8 {
                assert_eq!(positive.midpoint().mantissa(), &BigInt::from(205));
                assert_eq!(positive.midpoint().exponent(), -11);
            }
            contains_exact(&positive, &one_tenth);
            contains_exact(&evaluate("-0.1", precision).unwrap(), &-one_tenth.clone());
        }

        // The power of five is too wide to fit even at maximum binary
        // precision, so both directed denominator conversions matter.
        let tiny = BigRational::new(BigInt::from(-1), BigInt::from(10_u32).pow(15_000));
        contains_exact(
            &decimal::enclose_decimal("-1e-15000", MAX_BITS).unwrap(),
            &tiny,
        );
    }

    #[test]
    fn final_midpoint_rounding_encloses_exact_dyadic_literals() {
        for (literal, integer) in [("257", 257), ("-257", -257)] {
            let ball = evaluate(literal, 8).unwrap();
            assert_eq!(ball.midpoint().precision(), 8);
            contains_exact(&ball, &BigRational::from_integer(BigInt::from(integer)));
        }
        assert!(evaluate("0.125", 8).unwrap().radius().is_zero());
    }
}
