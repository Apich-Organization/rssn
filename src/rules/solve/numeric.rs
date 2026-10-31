//! Numeric fallbacks for equations the symbolic kernel cannot solve, used
//! when the engine runs in numeric mode.
//!
//! * One equation in one unknown: the residual is sampled on a grid that
//!   is dense near zero and reaches `±200`; every sign change is refined
//!   by bisection and kept when the residual there is negligible (poles
//!   are rejected). Roots of even multiplicity are the zeros of the
//!   residual at sign changes of its (symbolic) derivative.
//! * A system: damped Newton iterations (Gauss–Newton when there are more
//!   equations than unknowns) from a deterministic spread of starting
//!   points, with a finite-difference Jacobian. A root is accepted when
//!   the residual of every equation is below `1e-9`; the search is not
//!   exhaustive.
//! * `nsolve(list(equations), list(unknowns), list(guesses))`: Newton from
//!   the given guess.

use super::as_expression;
use crate::backend::Compiled;
use crate::graph::Cx;
use crate::graph::Graph;
use crate::graph::Kernel;
use crate::graph::NodeId;
use crate::graph::OpFlags;
use crate::graph::OpId;
use crate::graph::Outcome;
use crate::graph::SymbolId;
use crate::graph::op::core;
use crate::rules::poly::best;

/// The equations as compiled functions of the unknowns (then the bound
/// parameters).
struct Compiled_ {
    functions: Vec<Box<dyn Compiled>>,
    parameters: Vec<f64>,
}

impl Compiled_ {
    fn new(
        graph: &mut Graph,
        equations: &[NodeId],
        unknowns: &[SymbolId],
        bindings: &[(SymbolId, f64)],
    ) -> Option<Self> {
        let mut inputs: Vec<SymbolId> = unknowns.to_vec();
        let mut parameters = Vec::new();
        for &(s, v) in bindings {
            if !inputs.contains(&s) {
                inputs.push(s);
                parameters.push(v);
            }
        }
        let mut functions = Vec::new();
        for &e in equations {
            let expr = as_expression(graph, e);
            let term = best(graph, expr)?;
            functions.push(crate::backend::compile(graph, term, &inputs).ok()?);
        }
        Some(Self {
            functions,
            parameters,
        })
    }

    fn residual(
        &self,
        point: &[f64],
    ) -> Vec<f64> {
        let mut args = point.to_vec();
        args.extend_from_slice(&self.parameters);
        self.functions.iter().map(|f| f.call(&args)).collect()
    }
}

fn solve_dense(
    mut a: Vec<Vec<f64>>,
    mut b: Vec<f64>,
) -> Option<Vec<f64>> {
    let n = b.len();
    for col in 0..n {
        let pivot = (col..n).max_by(|&i, &j| a[i][col].abs().total_cmp(&a[j][col].abs()))?;
        if a[pivot][col].abs() < 1e-14 {
            return None;
        }
        a.swap(col, pivot);
        b.swap(col, pivot);
        for row in col + 1..n {
            let factor = a[row][col] / a[col][col];
            let pivot_row = a[col].clone();
            for (v, p) in a[row].iter_mut().zip(&pivot_row).skip(col) {
                *v -= factor * p;
            }
            let v = b[col];
            b[row] -= factor * v;
        }
    }
    let mut x = vec![0.0; n];
    for row in (0..n).rev() {
        let tail: f64 = (row + 1..n).map(|k| a[row][k] * x[k]).sum();
        x[row] = (b[row] - tail) / a[row][row];
    }
    x.iter().all(|v| v.is_finite()).then_some(x)
}

fn norm(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

/// One damped Newton / Gauss–Newton run.
fn newton(
    system: &Compiled_,
    start: &[f64],
) -> Option<Vec<f64>> {
    let n = start.len();
    let mut x = start.to_vec();
    let mut r = system.residual(&x);
    if r.iter().any(|v| !v.is_finite()) {
        return None;
    }
    for _ in 0..80 {
        let m = r.len();
        // Jacobian by central differences.
        let mut jac = vec![vec![0.0; n]; m];
        for j in 0..n {
            let h = 1e-6 * (1.0 + x[j].abs());
            let (mut plus, mut minus) = (x.clone(), x.clone());
            plus[j] += h;
            minus[j] -= h;
            let (rp, rm) = (system.residual(&plus), system.residual(&minus));
            for i in 0..m {
                jac[i][j] = (rp[i] - rm[i]) / (2.0 * h);
            }
        }
        if jac.iter().flatten().any(|v| !v.is_finite()) {
            return None;
        }
        // Normal equations J^T J d = -J^T r (exact Newton when square).
        let mut a = vec![vec![0.0; n]; n];
        let mut b = vec![0.0; n];
        for i in 0..m {
            for j in 0..n {
                b[j] -= jac[i][j] * r[i];
                for k in 0..n {
                    a[j][k] += jac[i][j] * jac[i][k];
                }
            }
        }
        for (j, row) in a.iter_mut().enumerate() {
            row[j] += 1e-14;
        }
        let step = solve_dense(a, b)?;
        let current = norm(&r);
        let mut scale = 1.0;
        let mut accepted = false;
        for _ in 0..30 {
            let trial: Vec<f64> = x.iter().zip(&step).map(|(a, d)| a + scale * d).collect();
            let rt = system.residual(&trial);
            if rt.iter().all(|v| v.is_finite()) && norm(&rt) < current {
                x = trial;
                r = rt;
                accepted = true;
                break;
            }
            scale *= 0.5;
        }
        if !accepted {
            break;
        }
        if norm(&r) < 1e-13 {
            break;
        }
    }
    (r.iter().all(|v| v.abs() < 1e-9) && x.iter().all(|v| v.is_finite())).then_some(x)
}

/// Deterministic pseudo-random points in the box `[-radius, radius]^n`.
fn starts(
    n: usize,
    radius: f64,
    count: usize,
    seed: &mut u64,
) -> Vec<Vec<f64>> {
    (0..count)
        .map(|_| {
            (0..n)
                .map(|_| {
                    *seed = seed
                        .wrapping_mul(6_364_136_223_846_793_005)
                        .wrapping_add(1_442_695_040_888_963_407);
                    #[allow(clippy::cast_precision_loss)]
                    let unit = (*seed >> 11) as f64 / (1_u64 << 53) as f64;
                    radius * (2.0 * unit - 1.0)
                })
                .collect()
        })
        .collect()
}

/// All roots of the system found from a spread of starting points.
pub(super) fn system_roots(
    graph: &mut Graph,
    equations: &[NodeId],
    unknowns: &[SymbolId],
    bindings: &[(SymbolId, f64)],
) -> Option<Vec<Vec<f64>>> {
    let system = Compiled_::new(graph, equations, unknowns, bindings)?;
    let n = unknowns.len();
    let mut seed = 0x2545_f491_4f6c_dd1d_u64;
    let mut roots: Vec<Vec<f64>> = Vec::new();
    let mut all = vec![vec![0.0; n]];
    for (radius, count) in [(1.0, 40), (3.0, 80), (10.0, 80), (40.0, 40)] {
        all.extend(starts(n, radius, count, &mut seed));
    }
    for start in all {
        if let Some(root) = newton(&system, &start)
            && !roots.iter().any(|r| {
                r.iter()
                    .zip(&root)
                    .all(|(a, b)| (a - b).abs() < 1e-6 * (1.0 + a.abs()))
            })
        {
            roots.push(root);
        }
    }
    roots.sort_by(|a, b| {
        a.iter()
            .zip(b)
            .map(|(x, y)| x.total_cmp(y))
            .find(|o| o.is_ne())
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    Some(roots)
}

/// The abscissae of the scan grid: dense near zero, reaching `±200`.
fn grid() -> Vec<f64> {
    let steps = 8000_i32;
    (-steps..=steps)
        .map(|i| {
            let u = f64::from(i) / f64::from(steps);
            200.0 * u * u.abs() * u.abs()
        })
        .collect()
}

/// Sign changes of `g` on the grid, each refined by bisection.
fn sign_changes(
    g: &dyn Fn(f64) -> f64,
    points: &[f64],
) -> Vec<(f64, f64, f64)> {
    let mut out = Vec::new();
    let mut previous = g(points.first().copied().unwrap_or(0.0));
    for pair in points.windows(2) {
        let (a, b) = (pair[0], pair[1]);
        let fb = g(b);
        let fa = previous;
        previous = fb;
        if fa.is_nan() || fb.is_nan() {
            continue;
        }
        if fa == 0.0 {
            out.push((a, fa, fb));
            continue;
        }
        if fa * fb < 0.0 {
            let (mut lo, mut hi, mut flo) = (a, b, fa);
            for _ in 0..200 {
                let mid = f64::midpoint(lo, hi);
                let fm = g(mid);
                if fm.is_nan() {
                    break;
                }
                if flo * fm <= 0.0 {
                    hi = mid;
                } else {
                    lo = mid;
                    flo = fm;
                }
            }
            out.push((f64::midpoint(lo, hi), fa, fb));
        }
    }
    out
}

/// Real roots of one equation: sign changes of the residual on a grid, and
/// the roots of even multiplicity, which are the zeros of the residual at
/// sign changes of its derivative (a root of multiplicity `m` is a root of
/// multiplicity `m - 1` of `f'`).
pub(super) fn scalar_roots(
    graph: &mut Graph,
    equation: NodeId,
    unknown: SymbolId,
    bindings: &[(SymbolId, f64)],
) -> Option<Vec<f64>> {
    let system = Compiled_::new(graph, &[equation], &[unknown], bindings)?;
    let f = |t: f64| {
        system
            .residual(&[t])
            .first()
            .copied()
            .filter(|v| v.is_finite())
            .unwrap_or(f64::NAN)
    };
    let points = grid();
    let mut roots: Vec<f64> = Vec::new();
    for (r, fa, fb) in sign_changes(&f, &points) {
        let value = f(r).abs();
        let scale = fa.abs().max(fb.abs());
        // Exact zeros of the grid are reported as such; a sign change across
        // a pole leaves a large residual and is rejected.
        if fa == 0.0 || (value <= 1e-8 * (1.0 + scale.min(1e6)) && value < 1e-6) {
            roots.push(r);
        }
    }
    // Roots of even multiplicity through f'.
    let expr = as_expression(graph, equation);
    if let Some(term) = best(graph, expr) {
        let x_node = graph.symbol_node(unknown);
        if let Some(d) = crate::rules::calculus::derivative(graph, term, x_node)
            && let Some(dsystem) = Compiled_::new(graph, &[d], &[unknown], bindings)
        {
            let df = |t: f64| {
                dsystem
                    .residual(&[t])
                    .first()
                    .copied()
                    .filter(|v| v.is_finite())
                    .unwrap_or(f64::NAN)
            };
            // Round-off makes f' change sign many times near a root of high
            // multiplicity: nearby candidates form one cluster, represented
            // by the member where |f| is smallest.
            let mut candidates: Vec<f64> = Vec::new();
            for (r, _, _) in sign_changes(&df, &points) {
                let value = f(r).abs();
                let delta = 1e-3 * (1.0 + r.abs());
                let scale = f(r - delta).abs().max(f(r + delta).abs());
                if value.is_finite() && scale.is_finite() && value <= 1e-9 * scale + 1e-15 {
                    candidates.push(r);
                }
            }
            let mut clusters: Vec<Vec<f64>> = Vec::new();
            for r in candidates {
                match clusters.last_mut() {
                    | Some(cluster)
                        if cluster
                            .last()
                            .is_some_and(|&l| (r - l).abs() <= 1e-2 * (1.0 + r.abs())) =>
                    {
                        cluster.push(r);
                    },
                    | _ => clusters.push(vec![r]),
                }
            }
            for cluster in clusters {
                if let Some(best) = cluster
                    .into_iter()
                    .min_by(|a, b| f(*a).abs().total_cmp(&f(*b).abs()))
                {
                    roots.push(best);
                }
            }
        }
    }
    roots.sort_by(f64::total_cmp);
    roots.dedup_by(|a, b| (*a - *b).abs() < 1e-6 * (1.0 + a.abs()));
    // Round-off flips the sign of a residual that is flat to machine
    // precision many times: a dense chain of "roots" is one root of high
    // multiplicity, represented by the member with the smallest residual.
    let mut chains: Vec<Vec<f64>> = Vec::new();
    for r in roots {
        match chains.last_mut() {
            | Some(chain)
                if chain
                    .last()
                    .is_some_and(|&l| (r - l).abs() <= 1e-2 * (1.0 + r.abs())) =>
            {
                chain.push(r);
            },
            | _ => chains.push(vec![r]),
        }
    }
    let mut roots: Vec<f64> = Vec::new();
    for chain in chains {
        if chain.len() >= 4 {
            roots.extend(
                chain
                    .into_iter()
                    .min_by(|a, b| f(*a).abs().total_cmp(&f(*b).abs())),
            );
        } else {
            roots.extend(chain);
        }
    }
    Some(roots)
}

/// Numeric solutions in numeric mode.
pub(super) struct NumericSolve {
    pub(super) solve: OpId,
    pub(super) nsolve: OpId,
}

impl Kernel for NumericSolve {
    fn ops(&self) -> Vec<OpId> {
        vec![self.solve, self.nsolve]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        if !cx.env.numeric || cx.graph.approx(cx.graph.find(node)).is_some() {
            return Outcome::Pass;
        }
        let heavy = |g: &Graph, n: NodeId| g.ops().get(g.op(n)).flags.has(OpFlags::HEAVY);
        if cx
            .graph
            .enodes(cx.graph.find(node))
            .any(|n| !heavy(cx.graph, n))
        {
            return Outcome::Pass;
        }
        let bindings = cx.env.bindings().to_vec();
        let graph = &mut *cx.graph;
        let children = graph.children(node).to_vec();
        let is_nsolve = graph.op(node) == self.nsolve;
        let (equation, unknown, guess) = match children.as_slice() {
            | &[e, u] if !is_nsolve => (e, u, None),
            | &[e, u, g] if is_nsolve => (e, u, Some(g)),
            | _ => return Outcome::Pass,
        };
        let equations: Vec<NodeId> = if graph.op(equation) == core::LIST {
            graph.children(equation).to_vec()
        } else {
            vec![equation]
        };
        if graph.op(unknown) == core::LIST {
            let nodes = graph.children(unknown).to_vec();
            let Some(symbols) = nodes
                .iter()
                .map(|&u| graph.symbol_of(u))
                .collect::<Option<Vec<_>>>()
            else {
                return Outcome::Pass;
            };
            if let Some(guess) = guess {
                let Some(system) = Compiled_::new(graph, &equations, &symbols, &bindings) else {
                    return Outcome::Pass;
                };
                let guesses: Option<Vec<f64>> = graph
                    .children(guess)
                    .to_vec()
                    .iter()
                    .map(|&g| best(graph, g).and_then(|t| graph.eval(t, cx.env)))
                    .collect();
                let Some(start) = guesses.filter(|g| g.len() == symbols.len()) else {
                    return Outcome::Pass;
                };
                return match newton(&system, &start) {
                    | Some(root) => {
                        let items: Vec<NodeId> = root.into_iter().map(|v| graph.float(v)).collect();
                        Outcome::Equal(graph.node(core::LIST, &items))
                    },
                    | None => Outcome::Pass,
                };
            }
            let Some(roots) = system_roots(graph, &equations, &symbols, &bindings) else {
                return Outcome::Pass;
            };
            let tuples: Vec<NodeId> = roots
                .into_iter()
                .map(|r| {
                    let items: Vec<NodeId> = r.into_iter().map(|v| graph.float(v)).collect();
                    graph.node(core::LIST, &items)
                })
                .collect();
            return Outcome::Equal(graph.node(core::LIST, &tuples));
        }
        if is_nsolve {
            return Outcome::Pass;
        }
        let Some(symbol) = graph.symbol_of(unknown) else {
            return Outcome::Pass;
        };
        // Polynomials have the certified kernel.
        let expr = as_expression(graph, equation);
        if super::numeric_polynomial(graph, expr, unknown, cx.env).is_some() {
            return Outcome::Pass;
        }
        match scalar_roots(graph, equation, symbol, &bindings) {
            | Some(roots) if !roots.is_empty() => {
                let items: Vec<NodeId> = roots.into_iter().map(|v| graph.float(v)).collect();
                Outcome::Equal(graph.node(core::LIST, &items))
            },
            | _ => Outcome::Pass,
        }
    }

    fn revisit(&self) -> bool {
        true
    }
}
