//! Quadratic, cone and semidefinite programs.
//!
//! | operator | value |
//! |---|---|
//! | `qp_solve(q, c, a, rel, b, bounds)` | `list(x, objective, iterations, active)`: the minimiser of `0.5 x'Qx + c'x` subject to `a x (rel) b` and the bounds |
//! | `socp_solve(c, a, rel, b, cones, bounds, max_iter)` | `list(x, objective, iterations, primal_residual, dual_residual, gap)` |
//! | `sdp_solve(c, f0, fs, a, rel, b, max_iter)` | `list(x, objective, iterations, primal_residual, dual_residual, min_eigenvalue)` |
//!
//! `rel` is a list of the symbols `le`, `ge` and `eq`, one per constraint row.
//! `bounds` is a list of pairs `list(lower, upper)` with the symbol `inf` for a
//! variable without an upper bound, e.g. `list(list(0, inf), list(-1, 1))`.
//! `cones` is a list of `list(a, b, c, d)` entries, each standing for
//! `||A x + b|| <= c . x + d`; `fs` is the list of matrices `F_i` of the linear
//! matrix inequality `F_0 + sum x_i F_i >= 0`.
//!
//! A program the solver cannot solve (an infeasible or unbounded one, a matrix
//! that is not positive semidefinite where the method needs one) is a
//! [`DiagKind::Domain`] finding: the call stays a request instead of returning
//! a point that is not a minimiser.
//!
//! The sparse sequential-quadratic-programming subproblem of
//! [`crate::kernels::milp::sqp`] has no operator form: it takes a Rust oracle
//! for the objective, its gradient and its sparse Hessian, which a list of
//! numbers cannot carry. It is a library API for a host that owns the model.

use super::bridge::Data;
use super::bridge::Entry;
use super::bridge::bad;
use super::bridge::count;
use super::bridge::require;
use super::bridge::rows;
use super::bridge::vector;
use crate::kernel::diag::DiagKind;
use crate::kernels::convex::Relation;
use crate::kernels::convex::SocCone;

/// Most variables a program accepts.
const MAX_VARIABLES: usize = 2_000;
/// Most constraint rows a program accepts.
const MAX_ROWS: usize = 20_000;

/// Every operator of this family.
pub(super) const ENTRIES: &[Entry] = &[
    Entry {
        name: "qp_solve",
        params: "tttttt",
        eval: qp_solve,
        doc: "`qp_solve(q, c, a, rel, b, bounds)`: `list(x, objective, iterations, active)` minimising `0.5 x'Qx + c'x` subject to `a x (rel) b` and `bounds` (`list(list(lo, hi), ..)`, `inf` for no upper bound).",
    },
    Entry {
        name: "socp_solve",
        params: "ttttttd",
        eval: socp_solve,
        doc: "`socp_solve(c, a, rel, b, cones, bounds, max_iter)`: `list(x, objective, iterations, primal_residual, dual_residual, gap)` of a second-order cone program; each cone is `list(a, b, c, d)` standing for `||A x + b|| <= c . x + d`.",
    },
    Entry {
        name: "sdp_solve",
        params: "ttttttd",
        eval: sdp_solve,
        doc: "`sdp_solve(c, f0, fs, a, rel, b, max_iter)`: `list(x, objective, iterations, primal_residual, dual_residual, min_eigenvalue)` minimising `c . x` subject to `F_0 + sum x_i F_i >= 0` and the equality rows `a x = b`.",
    },
];

/// A vector argument, copied.
pub(super) fn numbers(
    args: &[Data],
    at: usize,
) -> Option<Vec<f64>> {
    let values = vector(args, at)?;
    require(
        !values.is_empty() && values.len() <= MAX_VARIABLES,
        Some(at),
        DiagKind::Shape,
        format_args!(
            "expected 1 to {MAX_VARIABLES} numbers, got {}",
            values.len()
        ),
    )?;
    Some(values.to_vec())
}

/// A matrix argument with `width` columns when `width` is known.
pub(super) fn matrix(
    args: &[Data],
    at: usize,
    width: Option<usize>,
) -> Option<Vec<Vec<f64>>> {
    let Some(entries) = args.get(at).and_then(rows) else {
        return bad(
            at,
            DiagKind::Shape,
            format_args!("expected a list of rows of numbers"),
        );
    };
    if entries.is_empty() || entries.len() > MAX_ROWS {
        return bad(
            at,
            DiagKind::Shape,
            format_args!("expected 1 to {MAX_ROWS} rows, got {}", entries.len()),
        );
    }
    if let Some(expected) = width
        && entries.iter().any(|row| row.len() != expected)
    {
        return bad(
            at,
            DiagKind::Shape,
            format_args!("every row must have {expected} entries"),
        );
    }
    Some(entries)
}

/// The symbols of a list argument.
pub(super) fn symbols<'a>(
    args: &'a [Data],
    at: usize,
) -> Option<Vec<&'a str>> {
    let Some(Data::Tree(parts)) = args.get(at) else {
        return bad(
            at,
            DiagKind::Shape,
            format_args!("expected a list of symbols"),
        );
    };
    let mut names = Vec::with_capacity(parts.len());
    for part in parts {
        match part {
            | Data::Name(name) => names.push(name.as_str()),
            | _ => {
                return bad(
                    at,
                    DiagKind::Shape,
                    format_args!("every entry must be a symbol"),
                );
            },
        }
    }
    Some(names)
}

/// The relation symbols of `rel`, one per constraint row.
pub(super) fn relations(
    args: &[Data],
    at: usize,
    rows_count: usize,
) -> Option<Vec<Relation>> {
    let names = symbols(args, at)?;
    require(
        names.len() == rows_count,
        Some(at),
        DiagKind::Shape,
        format_args!(
            "expected one relation per constraint row ({rows_count}), got {}",
            names.len()
        ),
    )?;
    let mut relations = Vec::with_capacity(names.len());
    for name in names {
        match name {
            | "le" => relations.push(Relation::Le),
            | "ge" => relations.push(Relation::Ge),
            | "eq" => relations.push(Relation::Eq),
            | _ => {
                return bad(
                    at,
                    DiagKind::Domain,
                    format_args!("`{name}` is not a relation (le, ge, eq)"),
                );
            },
        }
    }
    Some(relations)
}

/// The `(lower, upper)` bounds of `variables` variables.
pub(super) fn bounds(
    args: &[Data],
    at: usize,
    variables: usize,
) -> Option<(Vec<f64>, Vec<Option<f64>>)> {
    let Some(Data::Tree(entries)) = args.get(at) else {
        return bad(
            at,
            DiagKind::Shape,
            format_args!("expected a list of pairs list(lower, upper)"),
        );
    };
    require(
        entries.len() == variables,
        Some(at),
        DiagKind::Shape,
        format_args!(
            "expected one bound pair per variable ({variables}), got {}",
            entries.len()
        ),
    )?;
    let mut lower = Vec::with_capacity(variables);
    let mut upper = Vec::with_capacity(variables);
    for (index, entry) in entries.iter().enumerate() {
        let Data::Tree(pair) = entry else {
            return bad(
                at,
                DiagKind::Shape,
                format_args!("bound {index} must be a pair list(lower, upper)"),
            );
        };
        let [low, high] = pair.as_slice() else {
            return bad(
                at,
                DiagKind::Shape,
                format_args!("bound {index} must have exactly two entries"),
            );
        };
        let Data::Scalar(low) = low else {
            return bad(
                at,
                DiagKind::Shape,
                format_args!("the lower bound {index} must be a number"),
            );
        };
        lower.push(*low);
        match high {
            | Data::Scalar(value) => upper.push(Some(*value)),
            | Data::Name(name) if name == "inf" || name == "infinity" => upper.push(None),
            | _ => {
                return bad(
                    at,
                    DiagKind::Shape,
                    format_args!("the upper bound {index} must be a number or `inf`"),
                );
            },
        }
    }
    Some((lower, upper))
}

fn qp_solve(args: &[Data]) -> Option<Data> {
    let c = numbers(args, 1)?;
    let variables = c.len();
    let q = matrix(args, 0, Some(variables))?;
    let a = matrix(args, 2, Some(variables))?;
    let rel = relations(args, 3, a.len())?;
    let b = numbers(args, 4)?;
    require(
        b.len() == a.len(),
        Some(4),
        DiagKind::Shape,
        format_args!("expected {} right-hand sides, got {}", a.len(), b.len()),
    )?;
    let (lower, upper) = bounds(args, 5, variables)?;
    let problem = crate::kernels::convex::QpProblem {
        q,
        c,
        a,
        rel,
        b,
        lower,
        upper,
    };
    match crate::kernels::convex::solve_qp(&problem, 10_000) {
        | Ok(solution) => {
            #[allow(clippy::cast_precision_loss)]
            let iterations = solution.iterations as f64;
            Some(Data::Tree(vec![
                Data::vector(solution.x),
                Data::Scalar(solution.objective),
                Data::Scalar(iterations),
                Data::vector(
                    solution
                        .active
                        .into_iter()
                        .map(|index| {
                            #[allow(clippy::cast_precision_loss)]
                            let value = index as f64;
                            value
                        })
                        .collect(),
                ),
            ]))
        },
        | Err(message) => bad(0, DiagKind::Domain, format_args!("{message}")),
    }
}

fn socp_solve(args: &[Data]) -> Option<Data> {
    let c = numbers(args, 0)?;
    let variables = c.len();
    let a = matrix(args, 1, Some(variables))?;
    let rel = relations(args, 2, a.len())?;
    let b = numbers(args, 3)?;
    require(
        b.len() == a.len(),
        Some(3),
        DiagKind::Shape,
        format_args!("expected {} right-hand sides, got {}", a.len(), b.len()),
    )?;
    let Some(Data::Tree(entries)) = args.get(4) else {
        return bad(
            4,
            DiagKind::Shape,
            format_args!("expected a list of cones list(a, b, c, d)"),
        );
    };
    let mut cones = Vec::with_capacity(entries.len());
    for (index, entry) in entries.iter().enumerate() {
        let Data::Tree(parts) = entry else {
            return bad(
                4,
                DiagKind::Shape,
                format_args!("cone {index} must be list(a, b, c, d)"),
            );
        };
        let [inner, offset, weights, constant] = parts.as_slice() else {
            return bad(
                4,
                DiagKind::Shape,
                format_args!("cone {index} must have exactly four entries"),
            );
        };
        let Some(inner) = rows(inner) else {
            return bad(
                4,
                DiagKind::Shape,
                format_args!("cone {index} needs a matrix of rows inside the norm"),
            );
        };
        if inner.iter().any(|row| row.len() != variables) {
            return bad(
                4,
                DiagKind::Shape,
                format_args!("cone {index} has a row whose length is not {variables}"),
            );
        }
        let Some(offset) = numbers_of(offset) else {
            return bad(
                4,
                DiagKind::Shape,
                format_args!("cone {index} needs a vector inside the norm"),
            );
        };
        let Some(weights) = numbers_of(weights) else {
            return bad(
                4,
                DiagKind::Shape,
                format_args!("cone {index} needs a vector on the right-hand side"),
            );
        };
        require(
            weights.len() == variables,
            Some(4),
            DiagKind::Shape,
            format_args!("cone {index} right-hand side must have {variables} entries"),
        )?;
        let Data::Scalar(constant) = constant else {
            return bad(
                4,
                DiagKind::Shape,
                format_args!("cone {index} needs a number as its constant"),
            );
        };
        cones.push(SocCone {
            a: inner,
            b: offset,
            c: weights,
            d: *constant,
        });
    }
    let (lower, upper) = bounds(args, 5, variables)?;
    let max_iter = count(args, 6, 100_000)?;
    let problem = crate::kernels::convex::SocpProblem {
        c,
        a,
        rel,
        b,
        cones,
        lower,
        upper,
    };
    match crate::kernels::convex::solve_socp(&problem, max_iter.max(1)) {
        | Ok(solution) => {
            #[allow(clippy::cast_precision_loss)]
            let iterations = solution.iterations as f64;
            Some(Data::Tree(vec![
                Data::vector(solution.x),
                Data::Scalar(solution.objective),
                Data::Scalar(iterations),
                Data::Scalar(solution.primal_residual),
                Data::Scalar(solution.dual_residual),
                Data::Scalar(solution.gap),
            ]))
        },
        | Err(message) => bad(0, DiagKind::Domain, format_args!("{message}")),
    }
}

/// The numbers of one part of a nested list.
pub(super) fn numbers_of(data: &Data) -> Option<Vec<f64>> {
    match data {
        | Data::Array { shape, data } if shape.len() == 1 => Some(data.clone()),
        | _ => None,
    }
}

fn sdp_solve(args: &[Data]) -> Option<Data> {
    let c = numbers(args, 0)?;
    let variables = c.len();
    let Some(f0) = args.get(1).and_then(rows) else {
        return bad(
            1,
            DiagKind::Shape,
            format_args!("expected the constant matrix F0 as a list of rows"),
        );
    };
    let size = f0.len();
    require(
        size > 0 && f0.iter().all(|row| row.len() == size),
        Some(1),
        DiagKind::Shape,
        format_args!("F0 must be a non-empty square matrix"),
    )?;
    let Some(Data::Tree(entries)) = args.get(2) else {
        return bad(
            2,
            DiagKind::Shape,
            format_args!("expected a list of matrices F_i"),
        );
    };
    require(
        entries.len() == variables,
        Some(2),
        DiagKind::Shape,
        format_args!(
            "expected one matrix per variable ({variables}), got {}",
            entries.len()
        ),
    )?;
    let mut fs = Vec::with_capacity(variables);
    for (index, entry) in entries.iter().enumerate() {
        let Some(matrix) = rows(entry) else {
            return bad(
                2,
                DiagKind::Shape,
                format_args!("F_{index} must be a list of rows"),
            );
        };
        if matrix.len() != size || matrix.iter().any(|row| row.len() != size) {
            return bad(
                2,
                DiagKind::Shape,
                format_args!("F_{index} must be {size} x {size}"),
            );
        }
        fs.push(matrix);
    }
    let a = matrix(args, 3, Some(variables))?;
    let rel = relations(args, 4, a.len())?;
    let b = numbers(args, 5)?;
    require(
        b.len() == a.len(),
        Some(5),
        DiagKind::Shape,
        format_args!("expected {} right-hand sides, got {}", a.len(), b.len()),
    )?;
    let max_iter = count(args, 6, 100_000)?;
    let problem = crate::kernels::convex::SdpProblem { c, f0, fs, a, rel, b };
    match crate::kernels::convex::solve_sdp(&problem, max_iter.max(1)) {
        | Ok(solution) => {
            #[allow(clippy::cast_precision_loss)]
            let iterations = solution.iterations as f64;
            Some(Data::Tree(vec![
                Data::vector(solution.x),
                Data::Scalar(solution.objective),
                Data::Scalar(iterations),
                Data::Scalar(solution.primal_residual),
                Data::Scalar(solution.dual_residual),
                Data::Scalar(solution.min_eigenvalue),
            ]))
        },
        | Err(message) => bad(0, DiagKind::Domain, format_args!("{message}")),
    }
}
