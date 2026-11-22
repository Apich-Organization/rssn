//! Mixed-integer linear programs.
//!
//! | operator | value |
//! |---|---|
//! | `milp_solve(c, a, rel, b, bounds, kinds, max_nodes)` | `list(x, objective, nodes, gap, optimal)`: the best feasible integer point found by branch and bound, its objective, the number of nodes explored, the relative gap to the incumbent and whether it is proven optimal |
//!
//! `rel` is a list of the symbols `le`, `ge` and `eq`; `bounds` is a list of
//! pairs `list(lower, upper)` with `inf` for a variable without an upper bound;
//! `kinds` is a list of the symbols `continuous`, `integer` and `binary`, one
//! per variable (`binary` also fixes the bounds to `[0, 1]`).
//!
//! Branch and bound relaxes the integrality to a linear program solved by the
//! tested two-phase simplex of [`crate::kernels::optim`], branches on the most
//! fractional variable, and prunes by bound. When the node budget runs out the
//! best incumbent is returned with `optimal` zero and the gap it still has to
//! close, which is an honest answer for an anytime solver: a host that needs a
//! certificate asks again with a larger budget rather than trusting a point
//! the search did not finish proving.

use super::bridge::Data;
use super::bridge::Entry;
use super::bridge::bad;
use super::bridge::count;
use super::bridge::require;
use super::conic::bounds;
use super::conic::matrix;
use super::conic::numbers;
use super::conic::relations;
use super::conic::symbols;
use crate::kernel::diag::DiagKind;
use crate::kernels::milp::VarKind;

/// Every operator of this family.
pub(super) const ENTRIES: &[Entry] = &[Entry {
    name: "milp_solve",
    params: "ttttttd",
    eval: milp_solve,
    doc: "`milp_solve(c, a, rel, b, bounds, kinds, max_nodes)`: `list(x, objective, nodes, gap, optimal)` of a mixed-integer linear program by branch and bound; `kinds` is a list of `continuous`, `integer` or `binary`.",
}];

fn milp_solve(args: &[Data]) -> Option<Data> {
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
    let (mut lower, mut upper) = bounds(args, 4, variables)?;
    let names = symbols(args, 5)?;
    require(
        names.len() == variables,
        Some(5),
        DiagKind::Shape,
        format_args!(
            "expected one kind per variable ({variables}), got {}",
            names.len()
        ),
    )?;
    let mut kinds = Vec::with_capacity(variables);
    for (index, name) in names.into_iter().enumerate() {
        match name {
            | "continuous" => kinds.push(VarKind::Continuous),
            | "integer" => kinds.push(VarKind::Integer),
            | "binary" => {
                if let Some(slot) = lower.get_mut(index) {
                    *slot = slot.max(0.0);
                }
                if let Some(slot) = upper.get_mut(index) {
                    *slot = Some(slot.map_or(1.0, |value| value.min(1.0)));
                }
                kinds.push(VarKind::Binary);
            },
            | _ => {
                return bad(
                    5,
                    DiagKind::Domain,
                    format_args!("`{name}` is not a variable kind (continuous, integer, binary)"),
                );
            },
        }
    }
    let max_nodes = count(args, 6, 1_000_000)?;
    let problem = crate::kernels::milp::MilpProblem {
        c,
        a,
        rel,
        b,
        lower,
        upper,
        kinds,
    };
    match crate::kernels::milp::solve_milp(&problem, max_nodes.max(1)) {
        | Ok(solution) => {
            #[allow(clippy::cast_precision_loss)]
            let nodes = solution.nodes as f64;
            let optimal = if solution.optimal { 1.0 } else { 0.0 };
            Some(Data::Tree(vec![
                Data::vector(solution.x),
                Data::Scalar(solution.objective),
                Data::Scalar(nodes),
                Data::Scalar(solution.gap),
                Data::Scalar(optimal),
            ]))
        },
        | Err(message) => bad(0, DiagKind::Domain, format_args!("{message}")),
    }
}
