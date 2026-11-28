//! Dense linear programming with two-sided variable bounds.
//!
//! The solver is a thin layer over the tested two-phase tableau simplex in
//! [`crate::kernels::optim`]. That simplex only accepts non-negative
//! variables, so the problem is first shifted by its finite lower bounds,
//! `x = lower + y` with `y >= 0`, and every finite upper bound is turned into
//! one extra `y_j <= upper_j - lower_j` row. The shifted problem is handed to
//! [`crate::kernels::optim::simplex_with_max_iter`] with steepest-edge pricing
//! and the caller's pivot budget, and the optimum is mapped back to the
//! original variables.

use crate::kernels::optim::LpStatus as OptLpStatus;
use crate::kernels::optim::Pricing;
use crate::kernels::optim::simplex_with_max_iter;

/// Sense of a linear-programming row.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Relation {
    /// `a . x <= b`
    Le,
    /// `a . x >= b`
    Ge,
    /// `a . x == b`
    Eq,
}

/// A dense linear program `minimise c . x` subject to `a_i . x (rel_i) b_i`
/// and finite lower bounds with optional finite upper bounds.
#[derive(Debug, Clone, PartialEq)]
pub struct LpProblem {
    /// Objective coefficients, `minimise c . x`.
    pub c: Vec<f64>,
    /// Constraint rows; each row has one entry per variable.
    pub a: Vec<Vec<f64>>,
    /// Sense of each constraint row.
    pub rel: Vec<Relation>,
    /// Right-hand side of each constraint row.
    pub b: Vec<f64>,
    /// Finite lower bound of each variable.
    pub lower: Vec<f64>,
    /// Optional finite upper bound of each variable (`None` means `+inf`).
    pub upper: Vec<Option<f64>>,
}

/// Outcome of a linear program.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LpStatus {
    /// An optimal vertex was found.
    Optimal,
    /// No point satisfies every constraint and bound.
    Infeasible,
    /// The objective is unbounded below on the feasible set.
    Unbounded,
}

/// Result of [`solve_lp`].
#[derive(Debug, Clone, PartialEq)]
pub struct LpSolution {
    /// Outcome of the solve.
    pub status: LpStatus,
    /// Optimal point (empty unless `status` is [`LpStatus::Optimal`]).
    pub x: Vec<f64>,
    /// Optimal objective value (`NaN` when infeasible, `-inf` when unbounded).
    pub objective: f64,
    /// Number of simplex pivots performed.
    pub iterations: usize,
}

/// Solves a bounded dense linear program.
///
/// `max_iter` limits the total number of simplex pivots performed across
/// Phase I and Phase II.
///
/// # Errors
/// Returns an error string when the dimensions disagree, when any input entry
/// is not finite, when `max_iter` is zero, or when the simplex exceeds
/// `max_iter` pivots.
#[allow(
    clippy::arithmetic_side_effects,
    reason = "floating coordinate shifts are checked by the simplex input validation"
)]
pub fn solve_lp(
    p: &LpProblem,
    max_iter: usize,
) -> Result<LpSolution, String> {
    let n = validate(p)?;
    if max_iter == 0 {
        return Err("iteration budget is zero".to_string());
    }

    // A finite upper bound below its lower bound is an immediate infeasibility.
    let finite_upper = p
        .upper
        .iter()
        .zip(&p.lower)
        .enumerate()
        .filter_map(|(j, (&upper, &lower))| upper.map(|upper| (j, upper, lower)));
    if finite_upper.clone().any(|(_, upper, lower)| upper < lower) {
        return Ok(infeasible_solution(0));
    }

    let capacity = p.a.len().saturating_add(finite_upper.clone().count());
    let mut rows = Vec::with_capacity(capacity);
    let mut rhs = Vec::with_capacity(capacity);
    let mut mapped = Vec::with_capacity(capacity);

    // Shift the rows into the non-negative variables `y = x - lower`.
    for ((row, &value), &sense) in p.a.iter().zip(&p.b).zip(&p.rel) {
        let shifted = value - dot(row, &p.lower);
        rows.push(row.clone());
        rhs.push(shifted);
        mapped.push(map_relation(sense));
    }

    // One `y_j <= upper_j - lower_j` row per finite upper bound.
    for (j, upper, lower) in finite_upper {
        let mut row = vec![0.0; n];
        if let Some(pivot) = row.get_mut(j) {
            *pivot = 1.0;
        }
        rows.push(row);
        rhs.push(upper - lower);
        mapped.push(crate::kernels::optim::Relation::Le);
    }

    let result =
        simplex_with_max_iter(&p.c, &rows, &mapped, &rhs, Pricing::SteepestEdge, max_iter)?;

    match result.status {
        | OptLpStatus::Optimal { x, value: _ } => {
            let x_full = p
                .lower
                .iter()
                .zip(x)
                .map(|(lower, y)| lower + y)
                .collect::<Vec<_>>();
            let objective = dot(&p.c, &x_full);
            Ok(LpSolution {
                status: LpStatus::Optimal,
                x: x_full,
                objective,
                iterations: result.iterations,
            })
        },
        | OptLpStatus::Infeasible => Ok(infeasible_solution(result.iterations)),
        | OptLpStatus::Unbounded => {
            Ok(LpSolution {
                status: LpStatus::Unbounded,
                x: Vec::new(),
                objective: f64::NEG_INFINITY,
                iterations: result.iterations,
            })
        },
    }
}

const fn infeasible_solution(iterations: usize) -> LpSolution {
    LpSolution {
        status: LpStatus::Infeasible,
        x: Vec::new(),
        objective: f64::NAN,
        iterations,
    }
}

fn map_relation(r: Relation) -> crate::kernels::optim::Relation {
    match r {
        | Relation::Le => crate::kernels::optim::Relation::Le,
        | Relation::Ge => crate::kernels::optim::Relation::Ge,
        | Relation::Eq => crate::kernels::optim::Relation::Eq,
    }
}

fn dot(
    a: &[f64],
    b: &[f64],
) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// Validates shapes and finiteness, returning the number of variables.
fn validate(p: &LpProblem) -> Result<usize, String> {
    let n = p.c.len();
    if p.lower.len() != n || p.upper.len() != n {
        return Err("bound vectors must match the number of variables".to_string());
    }
    if p.rel.len() != p.a.len() || p.b.len() != p.a.len() {
        return Err("constraint rows, senses and right-hand sides must agree".to_string());
    }
    if p.c.iter().any(|v| !v.is_finite()) {
        return Err("objective coefficients must be finite".to_string());
    }
    if p.lower.iter().any(|v| !v.is_finite()) {
        return Err("lower bounds must be finite".to_string());
    }
    for v in p.upper.iter().flatten() {
        if !v.is_finite() {
            return Err("upper bounds must be finite when present".to_string());
        }
    }
    for row in &p.a {
        if row.len() != n {
            return Err("constraint rows must have one entry per variable".to_string());
        }
        if row.iter().any(|v| !v.is_finite()) {
            return Err("constraint coefficients must be finite".to_string());
        }
    }
    if p.b.iter().any(|v| !v.is_finite()) {
        return Err("right-hand sides must be finite".to_string());
    }
    Ok(n)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(
        a: f64,
        b: f64,
    ) -> bool {
        (a - b).abs() <= 1e-8
    }

    #[test]
    fn halfspace_lower_bound() {
        // min x0 + x1  s.t.  x0 + x1 >= 1,  x >= 0.
        let p = LpProblem {
            c: vec![1.0, 1.0],
            a: vec![vec![1.0, 1.0]],
            rel: vec![Relation::Ge],
            b: vec![1.0],
            lower: vec![0.0, 0.0],
            upper: vec![None, None],
        };
        let s = solve_lp(&p, 1000).unwrap();
        assert_eq!(s.status, LpStatus::Optimal);
        assert!(close(s.objective, 1.0));
        assert!(close(s.x[0] + s.x[1], 1.0));
    }

    #[test]
    fn box_constrained_lp() {
        // max x0 + 2 x1  s.t.  x0 + x1 <= 4,  x0 <= 2,  x1 <= 3,  x >= 0.
        // Optimum (1, 3) with objective -7 for the minimisation form.
        let p = LpProblem {
            c: vec![-1.0, -2.0],
            a: vec![vec![1.0, 1.0], vec![1.0, 0.0], vec![0.0, 1.0]],
            rel: vec![Relation::Le, Relation::Le, Relation::Le],
            b: vec![4.0, 2.0, 3.0],
            lower: vec![0.0, 0.0],
            upper: vec![Some(2.0), Some(3.0)],
        };
        let s = solve_lp(&p, 1000).unwrap();
        assert_eq!(s.status, LpStatus::Optimal);
        assert!(close(s.x[0], 1.0), "x0 = {}", s.x[0]);
        assert!(close(s.x[1], 3.0), "x1 = {}", s.x[1]);
        assert!(close(s.objective, -7.0));
    }

    #[test]
    fn equality_row() {
        // min x0  s.t.  x0 + x1 = 2,  x >= 0  ->  (0, 2).
        let p = LpProblem {
            c: vec![1.0, 0.0],
            a: vec![vec![1.0, 1.0]],
            rel: vec![Relation::Eq],
            b: vec![2.0],
            lower: vec![0.0, 0.0],
            upper: vec![None, None],
        };
        let s = solve_lp(&p, 1000).unwrap();
        assert_eq!(s.status, LpStatus::Optimal);
        assert!(close(s.x[0], 0.0));
        assert!(close(s.x[1], 2.0));
        assert!(close(s.objective, 0.0));
    }

    #[test]
    fn negative_rhs_is_shifted() {
        // min x0  s.t.  -x0 <= -3  (i.e. x0 >= 3),  x0 >= 0  ->  x0 = 3.
        let p = LpProblem {
            c: vec![1.0],
            a: vec![vec![-1.0]],
            rel: vec![Relation::Le],
            b: vec![-3.0],
            lower: vec![0.0],
            upper: vec![None],
        };
        let s = solve_lp(&p, 1000).unwrap();
        assert_eq!(s.status, LpStatus::Optimal);
        assert!(close(s.x[0], 3.0));
    }

    #[test]
    fn infeasible_bounds() {
        // x0 >= 5 and x0 <= 2 cannot both hold.
        let p = LpProblem {
            c: vec![1.0],
            a: vec![vec![1.0]],
            rel: vec![Relation::Ge],
            b: vec![5.0],
            lower: vec![0.0],
            upper: vec![Some(2.0)],
        };
        let s = solve_lp(&p, 1000).unwrap();
        assert_eq!(s.status, LpStatus::Infeasible);
    }

    #[test]
    fn unbounded_problem() {
        // min -x0  s.t.  x0 >= 1,  x0 >= 0, no upper bound.
        let p = LpProblem {
            c: vec![-1.0],
            a: vec![vec![1.0]],
            rel: vec![Relation::Ge],
            b: vec![1.0],
            lower: vec![0.0],
            upper: vec![None],
        };
        let s = solve_lp(&p, 1000).unwrap();
        assert_eq!(s.status, LpStatus::Unbounded);
    }

    #[test]
    fn dimension_mismatch_is_error() {
        let p = LpProblem {
            c: vec![1.0, 2.0],
            a: vec![vec![1.0]],
            rel: vec![Relation::Le],
            b: vec![1.0],
            lower: vec![0.0, 0.0],
            upper: vec![None, None],
        };
        assert!(solve_lp(&p, 100).is_err());
    }

    #[test]
    fn max_iter_budget_enforced() {
        // Requires 2 pivots: (0, 0) -> (0, 3) -> (1, 3).
        let p = LpProblem {
            c: vec![-1.0, -2.0],
            a: vec![vec![1.0, 1.0], vec![1.0, 0.0], vec![0.0, 1.0]],
            rel: vec![Relation::Le, Relation::Le, Relation::Le],
            b: vec![4.0, 2.0, 3.0],
            lower: vec![0.0, 0.0],
            upper: vec![Some(2.0), Some(3.0)],
        };
        let err = solve_lp(&p, 1).unwrap_err();
        assert!(err.contains("iteration limit"), "unexpected error: {err}");

        let ok = solve_lp(&p, 10).unwrap();
        assert_eq!(ok.status, LpStatus::Optimal);
        assert!(ok.iterations <= 10);
        assert!(close(ok.x[0], 1.0));
        assert!(close(ok.x[1], 3.0));
        assert!(close(ok.objective, -7.0));
    }
}
