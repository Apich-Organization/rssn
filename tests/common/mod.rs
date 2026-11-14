//! Shared helpers of the symbolic-audit integration tests.
#![allow(dead_code)]

use rssn::api::Config;
use rssn::api::Session;

/// What one request answered.
#[derive(Clone, Debug)]
pub struct Out {
    pub text: String,
    pub reduced: bool,
    pub value: Option<f64>,
    pub latex: String,
    pub typst: String,
}

/// A symbolic request with the standard rules.
pub fn sym(src: &str) -> Out {
    run(src, &Config::new())
}

/// A numeric request.
pub fn num(
    src: &str,
    tolerance: f64,
) -> Out {
    run(src, &Config::new().numeric(tolerance))
}

pub fn run(
    src: &str,
    config: &Config,
) -> Out {
    let s = Session::new();
    let t = s
        .parse(src)
        .unwrap_or_else(|e| panic!("parse `{src}`: {e}"));
    let a = s
        .compute(t, config)
        .unwrap_or_else(|e| panic!("compute `{src}`: {e}"));
    Out {
        text: a.term.to_string(),
        reduced: a.reduced,
        value: a.value,
        latex: a.term.to_latex(),
        typst: a.term.to_typst(),
    }
}

/// Evaluates `src` (no rules) at the bindings.
pub fn eval(
    src: &str,
    bindings: &[(&str, f64)],
) -> Option<f64> {
    let s = Session::new();
    let t = s.parse(src).ok()?;
    t.eval(bindings)
}

/// Whether the text contains a float literal (`1.5`, `1e-5`).
pub fn has_float(text: &str) -> bool {
    let b: Vec<char> = text.chars().collect();
    b.windows(3).any(|w| {
        (w[0].is_ascii_digit() && w[1] == '.' && w[2].is_ascii_digit())
            || (w[0].is_ascii_digit() && w[1] == 'e' && (w[2] == '-' || w[2].is_ascii_digit()))
    })
}
