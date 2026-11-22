//! Mixed-integer linear programming by branch and bound.
//!
//! The relaxation at each node is the dense linear program of
//! [`crate::kernels::convex::lp`]. The search is a best-bound enumeration: every
//! child is solved when it is created, so the heap always holds exact LP bounds
//! and its minimum is a valid global lower bound. Branching is on the most
//! fractional integer variable; pseudo-costs (initialised from the objective
//! coefficient magnitude) break ties and are refined as nodes are processed.

use std::cmp::Ordering;
use std::collections::BinaryHeap;

use crate::kernels::convex::lp::LpProblem;
use crate::kernels::convex::lp::LpStatus;
use crate::kernels::convex::lp::Relation;
use crate::kernels::convex::lp::solve_lp;

/// Kind of a decision variable.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VarKind {
    /// A real-valued variable.
    Continuous,
    /// A general integer variable.
    Integer,
    /// A binary variable taking values in `{0, 1}`.
    Binary,
}

/// A mixed-integer linear program `minimise c . x`.
#[derive(Debug, Clone, PartialEq)]
pub struct MilpProblem {
    /// Objective coefficients.
    pub c: Vec<f64>,
    /// Constraint rows.
    pub a: Vec<Vec<f64>>,
    /// Sense of each constraint row.
    pub rel: Vec<Relation>,
    /// Right-hand side of each constraint row.
    pub b: Vec<f64>,
    /// Finite lower bound of each variable.
    pub lower: Vec<f64>,
    /// Optional finite upper bound of each variable (`None` means `+inf`).
    pub upper: Vec<Option<f64>>,
    /// Kind of each variable.
    pub kinds: Vec<VarKind>,
}

/// Result of [`solve_milp`].
#[derive(Debug, Clone, PartialEq)]
pub struct MilpSolution {
    /// Best integer-feasible point found.
    pub x: Vec<f64>,
    /// Objective value at `x`.
    pub objective: f64,
    /// Number of branch-and-bound nodes processed.
    pub nodes: usize,
    /// Total number of simplex pivots over every relaxation.
    pub lp_iterations: usize,
    /// Relative gap to the best open bound (`0` when the search is complete).
    pub gap: f64,
    /// Whether the incumbent is proven optimal.
    pub optimal: bool,
}

/// A solved relaxation node of the search tree.
struct Node {
    lower: Vec<f64>,
    upper: Vec<Option<f64>>,
    bound: f64,
    sol: Vec<f64>,
    parent_bound: f64,
    branch: Option<(usize, bool, f64)>,
}

impl PartialEq for Node {
    fn eq(
        &self,
        other: &Self,
    ) -> bool {
        self.bound == other.bound
    }
}

impl Eq for Node {}

impl Ord for Node {
    fn cmp(
        &self,
        other: &Self,
    ) -> Ordering {
        // Reverse order so the binary heap pops the smallest bound first.
        other.bound.total_cmp(&self.bound)
    }
}

impl PartialOrd for Node {
    fn partial_cmp(
        &self,
        other: &Self,
    ) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Solves a mixed-integer linear program.
///
/// # Errors
/// Returns `"infeasible"` when the root relaxation is infeasible (or when no
/// integer-feasible point exists), `"unbounded"` when a relaxation is
/// unbounded, and a descriptive string for dimension or finiteness errors.
pub fn solve_milp(
    p: &MilpProblem,
    max_nodes: usize,
) -> Result<MilpSolution, String> {
    let n = validate(p)?;
    if max_nodes == 0 {
        return Err("node budget is zero".to_string());
    }

    let mut lp_iterations = 0usize;

    // Root relaxation.
    let root = make_lp(p, &p.lower, &p.upper);
    let root_sol = solve_lp(&root, 1_000_000)?;
    lp_iterations += root_sol.iterations;
    match root_sol.status {
        | LpStatus::Infeasible => return Err("infeasible".to_string()),
        | LpStatus::Unbounded => return Err("unbounded".to_string()),
        | LpStatus::Optimal => {},
    }

    let mut heap: BinaryHeap<Node> = BinaryHeap::new();
    heap.push(Node {
        lower: p.lower.clone(),
        upper: p.upper.clone(),
        bound: root_sol.objective,
        sol: root_sol.x.clone(),
        parent_bound: root_sol.objective,
        branch: None,
    });

    let mut pc_up = vec![0.0; n];
    let mut pc_dn = vec![0.0; n];
    for j in 0..n {
        let base = p.c[j].abs().max(1e-6);
        pc_up[j] = base;
        pc_dn[j] = base;
    }

    let mut incumbent: Option<(Vec<f64>, f64)> = None;
    let mut nodes = 0usize;
    let mut hit_limit = false;

    while let Some(node) = heap.pop() {
        if nodes >= max_nodes {
            hit_limit = true;
            heap.push(node);
            break;
        }
        nodes += 1;

        // Refine pseudo-costs from the branch that produced this node.
        if let Some((j, up, dist)) = node.branch {
            let observed = (node.bound - node.parent_bound) / dist.max(1e-3);
            let target = if up {
                &mut pc_up[j]
            } else {
                &mut pc_dn[j]
            };
            *target = 0.5 * *target + 0.5 * observed.max(0.0);
        }

        // Bound pruning.
        if let Some((_, inc)) = &incumbent {
            if node.bound >= *inc - 1e-9 {
                continue;
            }
        }

        let Some(j) = most_fractional(p, &node.sol, &pc_up, &pc_dn) else {
            // Integer-feasible: update the incumbent.
            let better = incumbent
                .as_ref()
                .is_none_or(|(_, inc)| node.bound < *inc - 1e-9);
            if better {
                incumbent = Some((node.sol.clone(), node.bound));
            }
            continue;
        };

        let v = node.sol[j];
        let floor_v = v.floor();
        let ceil_v = v.ceil();

        // Down branch: x_j <= floor(v).
        if floor_v >= node.lower[j] - 1e-12 {
            let mut upper = node.upper.clone();
            let capped = match upper[j] {
                | Some(u) => Some(u.min(floor_v)),
                | None => Some(floor_v),
            };
            upper[j] = capped;
            if let Some(u) = upper[j] {
                if u >= node.lower[j] - 1e-12 {
                    let dist = v - floor_v;
                    if let Some(child) = solve_child(
                        p,
                        &node.lower,
                        &upper,
                        node.bound,
                        (j, false, dist),
                        &mut lp_iterations,
                    )? {
                        heap.push(child);
                    }
                }
            }
        }

        // Up branch: x_j >= ceil(v).
        if node.upper[j].is_none_or(|u| ceil_v <= u + 1e-12) {
            let mut lower = node.lower.clone();
            lower[j] = lower[j].max(ceil_v);
            if node.upper[j].is_none_or(|u| lower[j] <= u + 1e-12) {
                let dist = ceil_v - v;
                if let Some(child) = solve_child(
                    p,
                    &lower,
                    &node.upper,
                    node.bound,
                    (j, true, dist),
                    &mut lp_iterations,
                )? {
                    heap.push(child);
                }
            }
        }
    }

    let (x, objective) = match incumbent {
        | Some(pair) => pair,
        | None => return Err("infeasible".to_string()),
    };

    let optimal = if hit_limit {
        false
    } else {
        heap.peek()
            .is_none_or(|node| node.bound >= objective - 1e-9)
    };
    let gap = if optimal {
        0.0
    } else {
        match heap.peek() {
            | Some(node) => ((objective - node.bound) / objective.abs().max(1.0)).max(0.0),
            | None => 0.0,
        }
    };

    Ok(MilpSolution {
        x,
        objective,
        nodes,
        lp_iterations,
        gap,
        optimal,
    })
}

/// Solves one child relaxation; returns `None` when it is infeasible or pruned.
fn solve_child(
    p: &MilpProblem,
    lower: &[f64],
    upper: &[Option<f64>],
    parent_bound: f64,
    branch: (usize, bool, f64),
    lp_iterations: &mut usize,
) -> Result<Option<Node>, String> {
    let lp = make_lp(p, lower, upper);
    let sol = solve_lp(&lp, 1_000_000)?;
    *lp_iterations += sol.iterations;
    match sol.status {
        | LpStatus::Infeasible => Ok(None),
        | LpStatus::Unbounded => Err("unbounded".to_string()),
        | LpStatus::Optimal => {
            Ok(Some(Node {
                lower: lower.to_vec(),
                upper: upper.to_vec(),
                bound: sol.objective,
                sol: sol.x,
                parent_bound,
                branch: Some(branch),
            }))
        },
    }
}

fn make_lp(
    p: &MilpProblem,
    lower: &[f64],
    upper: &[Option<f64>],
) -> LpProblem {
    LpProblem {
        c: p.c.clone(),
        a: p.a.clone(),
        rel: p.rel.clone(),
        b: p.b.clone(),
        lower: lower.to_vec(),
        upper: upper.to_vec(),
    }
}

/// Most fractional non-continuous variable, with pseudo-cost tie-breaking.
fn most_fractional(
    p: &MilpProblem,
    x: &[f64],
    pc_up: &[f64],
    pc_dn: &[f64],
) -> Option<usize> {
    let tol = 1e-6;
    let mut best: Option<usize> = None;
    let mut best_frac = tol;
    let mut best_pc = f64::NEG_INFINITY;
    for j in 0..p.kinds.len() {
        if p.kinds[j] == VarKind::Continuous {
            continue;
        }
        let v = x.get(j).copied().unwrap_or(0.0);
        let frac = v - v.floor();
        let f = frac.min(1.0 - frac);
        if f > tol {
            let pc = pc_up[j] + pc_dn[j];
            if f > best_frac + 1e-12 || ((f - best_frac).abs() <= 1e-12 && pc > best_pc) {
                best_frac = f;
                best_pc = pc;
                best = Some(j);
            }
        }
    }
    best
}

/// Validates shapes and finiteness, returning the number of variables.
fn validate(p: &MilpProblem) -> Result<usize, String> {
    let n = p.c.len();
    if p.lower.len() != n || p.upper.len() != n || p.kinds.len() != n {
        return Err("variable vectors must all have the same length".to_string());
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
    for u in &p.upper {
        if let Some(v) = u {
            if !v.is_finite() {
                return Err("upper bounds must be finite when present".to_string());
            }
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
    fn three_variable_knapsack() {
        // max 3 a + 4 b + 5 c  s.t.  2 a + 3 b + 4 c <= 5,  binary.
        // Best packing {a, b} with weight 5 and value 7.
        let p = MilpProblem {
            c: vec![-3.0, -4.0, -5.0],
            a: vec![vec![2.0, 3.0, 4.0]],
            rel: vec![Relation::Le],
            b: vec![5.0],
            lower: vec![0.0, 0.0, 0.0],
            upper: vec![Some(1.0), Some(1.0), Some(1.0)],
            kinds: vec![VarKind::Binary, VarKind::Binary, VarKind::Binary],
        };
        let s = solve_milp(&p, 1000).unwrap();
        assert!(s.optimal);
        assert!(close(s.objective, -7.0), "objective = {}", s.objective);
        assert!(close(s.x[0], 1.0));
        assert!(close(s.x[1], 1.0));
        assert!(close(s.x[2], 0.0));
    }

    #[test]
    fn binary_assignment() {
        // 2x2 assignment, cost [[1, 2], [4, 8]]; optimum pairs (1,2) and (2,1).
        let p = MilpProblem {
            c: vec![1.0, 2.0, 4.0, 8.0],
            a: vec![
                vec![1.0, 1.0, 0.0, 0.0],
                vec![0.0, 0.0, 1.0, 1.0],
                vec![1.0, 0.0, 1.0, 0.0],
                vec![0.0, 1.0, 0.0, 1.0],
            ],
            rel: vec![Relation::Eq, Relation::Eq, Relation::Eq, Relation::Eq],
            b: vec![1.0, 1.0, 1.0, 1.0],
            lower: vec![0.0, 0.0, 0.0, 0.0],
            upper: vec![Some(1.0), Some(1.0), Some(1.0), Some(1.0)],
            kinds: vec![
                VarKind::Binary,
                VarKind::Binary,
                VarKind::Binary,
                VarKind::Binary,
            ],
        };
        let s = solve_milp(&p, 1000).unwrap();
        assert!(s.optimal);
        assert!(close(s.objective, 6.0), "objective = {}", s.objective);
        // x12 = 1 and x21 = 1.
        assert!(close(s.x[1], 1.0));
        assert!(close(s.x[2], 1.0));
    }

    #[test]
    fn mixed_continuous_integer() {
        // max x + y  s.t.  x + y <= 3.5,  0 <= x <= 3 integer,  0 <= y <= 2.
        let p = MilpProblem {
            c: vec![-1.0, -1.0],
            a: vec![vec![1.0, 1.0]],
            rel: vec![Relation::Le],
            b: vec![3.5],
            lower: vec![0.0, 0.0],
            upper: vec![Some(3.0), Some(2.0)],
            kinds: vec![VarKind::Integer, VarKind::Continuous],
        };
        let s = solve_milp(&p, 1000).unwrap();
        assert!(s.optimal);
        assert!(close(s.objective, -3.5), "objective = {}", s.objective);
        assert!((s.x[0] - s.x[0].round()).abs() < 1e-6);
    }

    #[test]
    fn infeasible_root() {
        // x0 >= 5 and x0 <= 2 with an integer x0.
        let p = MilpProblem {
            c: vec![1.0],
            a: vec![vec![1.0]],
            rel: vec![Relation::Ge],
            b: vec![5.0],
            lower: vec![0.0],
            upper: vec![Some(2.0)],
            kinds: vec![VarKind::Integer],
        };
        let err = solve_milp(&p, 100).unwrap_err();
        assert_eq!(err, "infeasible");
    }
}
