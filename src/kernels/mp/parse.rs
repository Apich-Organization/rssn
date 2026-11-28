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
use super::real::MAX_BITS;
use super::real::Real;
use super::real::Round;

/// Longest expression text accepted.
pub(crate) const MAX_EXPRESSION_TEXT: usize = 4096;
/// Deepest parenthesised nesting accepted.
const MAX_DEPTH: usize = 64;
/// Most tokens accepted.
const MAX_TOKENS: usize = 2048;

/// One lexical token.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Token<'a> {
    /// A decimal literal, as written.
    Number(&'a str),
    /// An identifier: a constant or a function name.
    Name(&'a str),
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
fn tokenize(text: &str) -> Result<Vec<Token<'_>>, String> {
    if text.len() > MAX_EXPRESSION_TEXT {
        return Err("the expression is too long".to_string());
    }
    let mut tokens = Vec::new();
    let mut characters = text.char_indices().peekable();
    let digit = |&(_, c): &(usize, char)| c.is_ascii_digit();
    while let Some((start, character)) = characters.next() {
        let token = match character {
            | ' ' | '\t' | '\n' | '\r' => continue,
            | '+' => Token::Plus,
            | '-' => Token::Minus,
            | '*' => Token::Star,
            | '/' => Token::Slash,
            | '^' => Token::Caret,
            | '(' => Token::Open,
            | ')' => Token::Close,
            | ',' => Token::Comma,
            | '0'..='9' | '.' => {
                while characters
                    .next_if(|&(_, c)| c.is_ascii_digit() || c == '.')
                    .is_some()
                {}
                if characters
                    .peek()
                    .is_some_and(|&(_, c)| c == 'e' || c == 'E')
                {
                    let mut exponent = characters.clone();
                    exponent.next();
                    exponent.next_if(|&(_, c)| c == '+' || c == '-');
                    if exponent.next_if(digit).is_some() {
                        while exponent.next_if(digit).is_some() {}
                        characters = exponent;
                    }
                }
                let end = characters.peek().map_or(text.len(), |&(index, _)| index);
                Token::Number(&text[start..end])
            },
            | _ if character.is_alphabetic() || character == '_' => {
                while characters
                    .next_if(|&(_, c)| c.is_alphanumeric() || c == '_')
                    .is_some()
                {}
                let end = characters.peek().map_or(text.len(), |&(index, _)| index);
                Token::Name(&text[start..end])
            },
            | _ => return Err(format!("`{character}` is not allowed in an expression")),
        };
        tokens.push(token);
        if tokens.len() > MAX_TOKENS {
            return Err("the expression has too many tokens".to_string());
        }
    }
    if tokens.is_empty() {
        return Err("the expression is empty".to_string());
    }
    Ok(tokens)
}

/// Handles plain decimal expressions without losing exact rounding ties.
pub(super) fn round_literal(
    text: &str,
    digits: usize,
) -> Result<Option<String>, String> {
    match tokenize(text)?.as_slice() {
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
    tokens: &'a [Token<'a>],
    at: usize,
    depth: usize,
    prec: u32,
}

impl<'a> Parser<'a> {
    /// The next token, without consuming it.
    fn peek(&self) -> Option<&Token<'a>> {
        self.tokens.get(self.at)
    }

    /// Consumes the next token when it is `wanted`.
    fn eat(
        &mut self,
        wanted: &Token<'_>,
    ) -> bool {
        if self.peek() == Some(wanted) {
            self.at += 1;
            return true;
        }
        false
    }

    fn nested<T>(
        &mut self,
        parse: impl FnOnce(&mut Self) -> Result<T, String>,
    ) -> Result<T, String> {
        if self.depth >= MAX_DEPTH {
            return Err("the expression nests too deeply".into());
        }
        self.depth += 1;
        let result = parse(self);
        self.depth -= 1;
        result
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
            let exponent = self.nested(Self::power)?;
            return funcs::pow_ball(&base, &exponent, self.prec);
        }
        Ok(base)
    }

    /// Parses a unary sign.
    fn unary(&mut self) -> Result<Ball, String> {
        let mut negative = false;
        loop {
            if self.eat(&Token::Minus) {
                negative = !negative;
            } else if !self.eat(&Token::Plus) {
                break;
            }
        }
        let value = self.primary()?;
        Ok(if negative {
            value.neg()
        } else {
            value
        })
    }

    /// Parses a literal, a name, a call or a parenthesised expression.
    fn primary(&mut self) -> Result<Ball, String> {
        let token = self
            .peek()
            .copied()
            .ok_or_else(|| "the expression ends early".to_string())?;
        match token {
            | Token::Number(text) => {
                self.at += 1;
                decimal::enclose_decimal(text, self.prec)
            },
            | Token::Name(name) => {
                self.at += 1;
                if self.eat(&Token::Open) {
                    return self.nested(|parser| parser.function(name));
                }
                consts::named_ball(name, self.prec)
            },
            | Token::Open => {
                self.at += 1;
                let value = self.nested(Self::expression)?;
                if !self.eat(&Token::Close) {
                    return Err("missing `)`".to_string());
                }
                Ok(value)
            },
            | _ => Err("expected a number, a name or `(`".to_string()),
        }
    }

    fn function(
        &mut self,
        name: &str,
    ) -> Result<Ball, String> {
        let first = self.expression()?;
        let second = if self.eat(&Token::Comma) {
            Some(self.expression()?)
        } else {
            None
        };
        if self.eat(&Token::Comma) {
            return Err(format!("`{name}` accepts at most two arguments"));
        }
        if !self.eat(&Token::Close) {
            return Err("missing `)`".into());
        }
        match second {
            | Some(second) => call(name, &[first, second], self.prec),
            | None => call(name, &[first], self.prec),
        }
    }
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
    match name {
        | "exp" => {
            want(1)?;
            funcs::exp_ball(&arguments[0], prec)
        },
        | "ln" | "log" => {
            want(1)?;
            if name == "log" {
                return funcs::log10_ball(&arguments[0], prec);
            }
            funcs::ln_ball(&arguments[0], prec)
        },
        | "log2" => {
            want(1)?;
            funcs::log2_ball(&arguments[0], prec)
        },
        | "log10" => {
            want(1)?;
            funcs::log10_ball(&arguments[0], prec)
        },
        | "sqrt" => {
            want(1)?;
            arguments[0].sqrt(prec)
        },
        | "cbrt" => {
            want(1)?;
            arguments[0].map_monotone(&Real::cbrt, true, prec)
        },
        | "abs" => {
            want(1)?;
            arguments[0].abs(prec)
        },
        | "sin" => {
            want(1)?;
            funcs::sin_ball(&arguments[0], prec)
        },
        | "cos" => {
            want(1)?;
            funcs::cos_ball(&arguments[0], prec)
        },
        | "tan" => {
            want(1)?;
            funcs::tan_ball(&arguments[0], prec)
        },
        | "asin" => {
            want(1)?;
            funcs::asin_ball(&arguments[0], prec)
        },
        | "acos" => {
            want(1)?;
            funcs::acos_ball(&arguments[0], prec)
        },
        | "atan" => {
            want(1)?;
            funcs::atan_ball(&arguments[0], prec)
        },
        | "atan2" => {
            want(2)?;
            funcs::atan2_ball(&arguments[0], &arguments[1], prec)
        },
        | "sinh" => {
            want(1)?;
            funcs::sinh_ball(&arguments[0], prec)
        },
        | "cosh" => {
            want(1)?;
            funcs::cosh_ball(&arguments[0], prec)
        },
        | "tanh" => {
            want(1)?;
            funcs::tanh_ball(&arguments[0], prec)
        },
        | "exp2" => {
            want(1)?;
            funcs::exp2_ball(&arguments[0], prec)
        },
        | "hypot" => {
            want(2)?;
            funcs::hypot_ball(&arguments[0], &arguments[1], prec)
        },
        | "pow" => {
            want(2)?;
            funcs::pow_ball(&arguments[0], &arguments[1], prec)
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
    let tokens = tokenize(text)?;
    let mut parser = Parser {
        tokens: &tokens,
        at: 0,
        depth: 0,
        prec: prec.saturating_add(32).min(MAX_BITS),
    };
    let value = parser.expression()?;
    if parser.at != tokens.len() {
        return Err("trailing text after the expression".to_string());
    }
    value.rounded_midpoint(prec)
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
