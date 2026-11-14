//! Checked automatic differentiation of dimensionless elementary expressions.
//!
//! This is a host-facing bridge to the existing SSA tape interpreters. It
//! does not differentiate opaque host callbacks: such callbacks cannot carry
//! dual numbers through the ABI. Expressions use RSSN's arithmetic and
//! elementary syntax, with independent variables supplied in column order.
//! Opaque calls, branches, comparisons, gathers, and non-smooth intrinsics are
//! rejected rather than using the tape's finite-difference fallback rules.

use super::Intrinsic;
use super::tape::Inst;
use super::tape::Tape;
use super::tape::{
    self,
};
use crate::graph::Graph;
use crate::graph::rule::Program;
use crate::kernels::dense::Mat;
use std::collections::HashSet;

fn checked_tape(
    expressions: &[&str],
    variables: &[String],
    point: &[f64],
) -> Result<Tape, String> {
    if expressions.is_empty() || variables.is_empty() || point.len() != variables.len() {
        return Err(
            "AD requires nonempty expressions and variables, and one point entry per variable"
                .into(),
        );
    }
    if point.iter().any(|value| !value.is_finite()) {
        return Err("AD evaluation points must be finite".into());
    }
    let mut names = HashSet::new();
    for name in variables {
        let mut characters = name.chars();
        if !characters
            .next()
            .is_some_and(|c| c == '_' || c.is_alphabetic())
            || !characters.all(|c| c == '_' || c.is_alphanumeric())
            || !names.insert(name)
        {
            return Err("AD variable names must be unique identifiers".into());
        }
    }
    let mut graph = Graph::new();
    let mut program = Program::default();
    crate::rules::elementary()
        .install_into(&mut graph, &mut program)
        .map_err(|error| error.to_string())?;
    let inputs: Vec<_> = variables
        .iter()
        .map(|name| graph.interner_mut().symbol(name))
        .collect();
    let roots: Vec<_> = expressions
        .iter()
        .map(|expression| graph.parse(expression).map_err(|error| error.to_string()))
        .collect::<Result<_, _>>()?;
    let lowered = tape::lower(&graph, &roots, &inputs, &[]).map_err(|error| error.to_string())?;
    for instruction in &lowered.insts {
        if matches!(
            instruction,
            Inst::Call(..)
                | Inst::Extern1(..)
                | Inst::Extern2(..)
                | Inst::Cmp(..)
                | Inst::Select(..)
                | Inst::Gather(..)
                | Inst::Unary(Intrinsic::Abs | Intrinsic::Floor | Intrinsic::Ceil, _)
        ) {
            return Err("exact expression AD requires smooth elementary tape instructions; opaque or non-smooth operations are unsupported".into());
        }
    }
    let lowered = lowered.optimise();
    let mut values = vec![0.0; expressions.len()];
    lowered.eval(point, &[], &mut values);
    finite(&values)?;
    Ok(lowered)
}

fn finite(values: &[f64]) -> Result<(), String> {
    if values.iter().all(|value| value.is_finite()) {
        Ok(())
    } else {
        Err("expression value or derivative is not finite at the supplied point".into())
    }
}

/// Forward-mode gradient of an elementary scalar expression at a finite point.
/// Variable order determines output order; no step-size parameter is used.
///
/// # Errors
/// Invalid or duplicate variables, mismatched dimensions, unbound symbols,
/// unsupported operations, or non-finite values or derivatives.
pub fn forward_gradient(
    expression: &str,
    variables: &[String],
    point: &[f64],
) -> Result<Vec<f64>, String> {
    let tape = checked_tape(&[expression], variables, point)?;
    let gradient = tape.gradient(0, point, &[]);
    finite(&gradient)?;
    Ok(gradient)
}

/// Reverse-mode gradient of an elementary scalar expression in one reverse sweep.
///
/// # Errors
/// As [`forward_gradient`].
pub fn reverse_gradient(
    expression: &str,
    variables: &[String],
    point: &[f64],
) -> Result<Vec<f64>, String> {
    let tape = checked_tape(&[expression], variables, point)?;
    let gradient = tape.gradient_reverse(0, point, &[]);
    finite(&gradient)?;
    Ok(gradient)
}

/// Hessian of an elementary scalar expression by the existing hyper-dual tape.
/// Returns a square matrix in variable order, without finite differences.
///
/// # Errors
/// As [`forward_gradient`].
pub fn hessian(
    expression: &str,
    variables: &[String],
    point: &[f64],
) -> Result<Mat, String> {
    let tape = checked_tape(&[expression], variables, point)?;
    let data = tape.hessian(0, point, &[]);
    finite(&data)?;
    Ok(Mat {
        rows: point.len(),
        cols: point.len(),
        data,
    })
}

/// Jacobian of elementary expressions, with one row per expression and one
/// column per independent variable, using reverse accumulation per output.
///
/// # Errors
/// As [`forward_gradient`], including empty output lists.
pub fn jacobian(
    expressions: &[String],
    variables: &[String],
    point: &[f64],
) -> Result<Mat, String> {
    let sources: Vec<_> = expressions.iter().map(String::as_str).collect();
    let tape = checked_tape(&sources, variables, point)?;
    let data: Vec<_> = (0..expressions.len())
        .flat_map(|output| tape.gradient_reverse(output, point, &[]))
        .collect();
    finite(&data)?;
    Ok(Mat {
        rows: expressions.len(),
        cols: point.len(),
        data,
    })
}

/// Vector-Jacobian product `weights^T J` in one shared-tape reverse sweep.
///
/// # Errors
/// As [`forward_gradient`], plus mismatched or non-finite output weights.
pub fn vjp(
    expressions: &[String],
    variables: &[String],
    point: &[f64],
    weights: &[f64],
) -> Result<Vec<f64>, String> {
    if weights.len() != expressions.len() || weights.iter().any(|value| !value.is_finite()) {
        return Err("VJP requires one finite weight per expression".into());
    }
    let sources: Vec<_> = expressions.iter().map(String::as_str).collect();
    let tape = checked_tape(&sources, variables, point)?;
    let product = tape.vjp(point, &[], weights);
    finite(&product)?;
    Ok(product)
}
