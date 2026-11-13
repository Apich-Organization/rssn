//! Initial-value problems with a chosen integrator, as operators over
//! expressions.
//!
//! The right-hand sides are expressions in the named variables (and in the
//! independent variable `t`), compiled once with the current backend; every
//! other parameter must be a number. The answer is the final state as a list.
//! Anything else, a failed integration or an expression the backend cannot
//! compile leaves the call as a request. (`odeint` of the `ode` rule set is
//! the scalar problem with the symbolic solver's fallback; these take
//! systems and name the integrator.)
//!
//! | operator | value |
//! |---|---|
//! | `ode_solve(rhs, vars, t, y0, t0, t1, tol, method)` | `y(t1)` for `y' = rhs(t, y)`, `y(t0) = y0`, with relative tolerance `tol` (absolute `tol / 100`); `method` is `dopri5`, `rosenbrock23`, `radau5`, `bdf`, `sdirk3` or `exprb32` |
//! | `ode_symplectic(accel, vars, q0, p0, t0, t1, steps, method)` | `list(q, p)` at `t1` for `q' = p`, `p' = accel(q)` (unit masses) by a symplectic composition of Stormer-Verlet steps; `method` is `verlet`, `yoshida4`, `yoshida6` or `yoshida8` |
//! | `ode_imex(f_ex, f_im, vars, t, y0, t0, t1, steps, scheme)` | `y(t1)` for `y' = f_ex + f_im` with the first part explicit and the stiff second part implicit; `scheme` is `euler`, `ars222` or `ars343` |
//! | `ode_etdrk4(L, f, vars, t, u0, t0, t1, steps)` | `u(t1)` for `u' = L u + f(t, u)` by exponential time differencing, fourth order whatever the stiffness of the matrix `L` |

use crate::graph::Cx;
use crate::graph::Kernel;
use crate::graph::NodeId;
use crate::graph::OpId;
use crate::graph::Outcome;
use crate::graph::SymbolId;
use crate::kernels::dense::Mat;
use crate::kernels::ode_adaptive;
use crate::kernels::ode_adaptive::OdeOptions;
use crate::kernels::ode_ext;
use crate::kernels::ode_sdirk;
use crate::kernels::ode_stiff;

use super::bridge;
use super::functions::MAX_VARIABLES;
use super::functions::compile;
use super::functions::items;
use super::functions::list;
use super::functions::number;
use super::functions::numbers;
use super::functions::symbols;

/// Most fixed steps of an exponential or IMEX integration, and of a
/// symplectic one.
const MAX_STEPS: f64 = 200_000.0;
/// Most steps an adaptive integrator may attempt.
const MAX_ADAPTIVE_STEPS: usize = 200_000;

/// The operators of this module.
#[derive(Clone, Copy)]
pub(super) struct OdeOps {
    pub(super) solve: OpId,
    pub(super) symplectic: OpId,
    pub(super) imex: OpId,
    pub(super) etdrk4: OpId,
}

/// A right-hand side `f(t, y, out)` compiled from expressions in `t` and `vars`.
fn timed(
    cx: &mut Cx<'_>,
    terms: &[NodeId],
    t: SymbolId,
    vars: &[SymbolId],
) -> Option<impl Fn(f64, &[f64], &mut [f64]) + use<>> {
    (!vars.contains(&t) && terms.len() == vars.len()).then_some(())?;
    let mut all = vec![t];
    all.extend_from_slice(vars);
    let f = compile(cx.graph, terms, &all)?;
    Some(move |time: f64, y: &[f64], out: &mut [f64]| {
        let mut args = Vec::with_capacity(y.len() + 1);
        args.push(time);
        args.extend_from_slice(y);
        f(&args, out);
    })
}

/// The symbol a term names, as text.
fn symbol_text(
    cx: &Cx<'_>,
    node: NodeId,
) -> Option<String> {
    let symbol = cx.graph.symbol_of(node)?;
    Some(cx.graph.interner().symbol_name(symbol).to_owned())
}

/// A count in `1..=MAX_STEPS`.
fn steps(
    cx: &mut Cx<'_>,
    node: NodeId,
) -> Option<usize> {
    let x = number(cx, node)?;
    ((1.0..=MAX_STEPS).contains(&x) && x.fract() == 0.0).then_some(x as usize)
}

fn finite_state(state: Vec<f64>) -> Option<Vec<f64>> {
    state.iter().all(|x| x.is_finite()).then_some(state)
}

impl OdeOps {
    fn solve(
        cx: &mut Cx<'_>,
        args: &[NodeId],
    ) -> Option<NodeId> {
        let &[rhs, vars, t, y0, t0, t1, tol, method] = args else {
            return None;
        };
        let vars = symbols(cx.graph, vars)?;
        let t = cx.graph.symbol_of(t)?;
        let rhs = items(cx.graph, rhs);
        let y0 = numbers(cx, y0)?;
        let (t0, t1, tol) = (number(cx, t0)?, number(cx, t1)?, number(cx, tol)?);
        let method = symbol_text(cx, method)?;
        (y0.len() == vars.len() && (1e-13..=1e-2).contains(&tol)).then_some(())?;
        let f = timed(cx, &rhs, t, &vars)?;
        let options = OdeOptions {
            rtol: tol,
            atol: tol * 1e-2,
            max_steps: MAX_ADAPTIVE_STEPS,
            ..OdeOptions::default()
        };
        let last = |solution: ode_adaptive::OdeSolution| solution.y.last().cloned();
        let state = match method.as_str() {
            | "dopri5" => last(ode_adaptive::dopri5(&f, t0, t1, &y0, &options, &[]).ok()?),
            | "rosenbrock23" => last(ode_adaptive::rosenbrock23(&f, t0, t1, &y0, &options).ok()?),
            | "radau5" => last(ode_stiff::radau5(&f, t0, t1, &y0, &options).ok()?),
            | "bdf" => last(ode_stiff::bdf(&f, t0, t1, &y0, &options, 5).ok()?),
            | "exprb32" => last(ode_ext::exprb32(&f, None, t0, t1, &y0, &options).ok()?),
            | "sdirk3" => {
                let sdirk = ode_sdirk::SdirkOptions {
                    rtol: tol,
                    atol: tol * 1e-2,
                    max_steps: MAX_ADAPTIVE_STEPS,
                    ..ode_sdirk::SdirkOptions::default()
                };
                ode_sdirk::sdirk3_fd(&f, t0, t1, &y0, &sdirk)
                    .ok()?
                    .y
                    .last()
                    .cloned()
            },
            | _ => None,
        }?;
        Some(list(cx.graph, &finite_state(state)?))
    }

    fn symplectic(
        cx: &mut Cx<'_>,
        args: &[NodeId],
    ) -> Option<NodeId> {
        let &[accel, vars, q0, p0, t0, t1, count, method] = args else {
            return None;
        };
        let vars = symbols(cx.graph, vars)?;
        let rhs = items(cx.graph, accel);
        let (q0, p0) = (numbers(cx, q0)?, numbers(cx, p0)?);
        let (t0, t1) = (number(cx, t0)?, number(cx, t1)?);
        let count = steps(cx, count)?;
        let method = match symbol_text(cx, method)?.as_str() {
            | "verlet" => ode_ext::SplitMethod::Verlet,
            | "yoshida4" => ode_ext::SplitMethod::TripleJump(4),
            | "yoshida6" => ode_ext::SplitMethod::Yoshida6,
            | "yoshida8" => ode_ext::SplitMethod::Yoshida8,
            | _ => return None,
        };
        let n = vars.len();
        (rhs.len() == n && q0.len() == n && p0.len() == n).then_some(())?;
        let force = compile(cx.graph, &rhs, &vars)?;
        let acceleration = |q: &[f64]| {
            let mut out = vec![0.0; n];
            force(q, &mut out);
            out
        };
        let (q, p) =
            ode_ext::symplectic_split(acceleration, &q0, &p0, t0, t1, count, method).ok()?;
        let (q, p) = (finite_state(q)?, finite_state(p)?);
        let (q, p) = (list(cx.graph, &q), list(cx.graph, &p));
        Some(cx.graph.node(crate::graph::op::core::LIST, &[q, p]))
    }

    fn imex(
        cx: &mut Cx<'_>,
        args: &[NodeId],
    ) -> Option<NodeId> {
        let &[explicit, implicit, vars, t, y0, t0, t1, count, scheme] = args else {
            return None;
        };
        let vars = symbols(cx.graph, vars)?;
        let t = cx.graph.symbol_of(t)?;
        let (explicit, implicit) = (items(cx.graph, explicit), items(cx.graph, implicit));
        let y0 = numbers(cx, y0)?;
        let (t0, t1) = (number(cx, t0)?, number(cx, t1)?);
        let count = steps(cx, count)?;
        let scheme = match symbol_text(cx, scheme)?.as_str() {
            | "euler" => ode_ext::ImexScheme::Euler,
            | "ars222" => ode_ext::ImexScheme::Ars222,
            | "ars343" => ode_ext::ImexScheme::Ars343,
            | _ => return None,
        };
        (y0.len() == vars.len() && y0.len() * count <= 2_000_000).then_some(())?;
        let f_ex = timed(cx, &explicit, t, &vars)?;
        let f_im = timed(cx, &implicit, t, &vars)?;
        let options = ode_ext::ImexOptions {
            jac_im: None,
            steps: count,
            scheme,
        };
        let path = ode_ext::imex_rk_with(&f_ex, &f_im, t0, t1, &y0, &options).ok()?;
        Some(list(cx.graph, &finite_state(path.last()?.clone())?))
    }

    fn etdrk4(
        cx: &mut Cx<'_>,
        args: &[NodeId],
    ) -> Option<NodeId> {
        let &[linear, nonlinear, vars, t, u0, t0, t1, count] = args else {
            return None;
        };
        let vars = symbols(cx.graph, vars)?;
        let t = cx.graph.symbol_of(t)?;
        let nonlinear = items(cx.graph, nonlinear);
        let u0 = numbers(cx, u0)?;
        let (t0, t1) = (number(cx, t0)?, number(cx, t1)?);
        let count = steps(cx, count)?;
        let n = vars.len();
        let decoded = [bridge::decode(cx.graph, linear)?];
        let (rows, cols, data) = bridge::matrix(&decoded, 0)?;
        (rows == n && cols == n && u0.len() == n && n <= MAX_VARIABLES && n * count <= 200_000)
            .then_some(())?;
        let l = Mat {
            rows,
            cols,
            data: data.to_vec(),
        };
        let f = timed(cx, &nonlinear, t, &vars)?;
        let path = ode_ext::etdrk4(&l, &f, t0, t1, &u0, count).ok()?;
        Some(list(cx.graph, &finite_state(path.last()?.clone())?))
    }
}

impl Kernel for OdeOps {
    fn ops(&self) -> Vec<OpId> {
        vec![self.solve, self.symplectic, self.imex, self.etdrk4]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let op = cx.graph.op(node);
        let args = cx.graph.children(node).to_vec();
        let found = if op == self.solve {
            Self::solve(cx, &args)
        } else if op == self.symplectic {
            Self::symplectic(cx, &args)
        } else if op == self.imex {
            Self::imex(cx, &args)
        } else if op == self.etdrk4 {
            Self::etdrk4(cx, &args)
        } else {
            None
        };
        found.map_or(Outcome::Pass, Outcome::Equal)
    }
}

#[cfg(test)]
mod tests {
    use crate::rules::numeric::numeric;
    use crate::rules::testing::reduce_with;
    use crate::rules::testing::simplify;

    fn run(src: &str) -> Vec<f64> {
        simplify(&[numeric()], src)
            .replace(['(', ')', ','], " ")
            .split_whitespace()
            .filter(|word| *word != "list")
            .map(|word| word.parse().unwrap_or(f64::NAN))
            .collect()
    }

    fn close(
        got: &[f64],
        want: &[f64],
        tol: f64,
    ) {
        assert_eq!(got.len(), want.len(), "{got:?} vs {want:?}");
        for (g, w) in got.iter().zip(want) {
            assert!((g - w).abs() <= tol, "{got:?} vs {want:?}");
        }
    }

    #[test]
    fn every_integrator_solves_exponential_decay_and_rotation() {
        let decay = (-2.0_f64).exp();
        for (method, tol) in [
            ("dopri5", 1e-7),
            ("rosenbrock23", 1e-5),
            ("radau5", 1e-6),
            ("bdf", 1e-5),
            ("sdirk3", 1e-6),
            ("exprb32", 1e-6),
        ] {
            let got = run(&format!("ode_solve(-2*y, y, t, 1, 0, 1, 1e-9, {method})"));
            close(&got, &[decay], tol);
        }
        // y1' = y2, y2' = -y1 from (1, 0) is (cos t, -sin t).
        let got =
            run("ode_solve(list(y2, -y1), list(y1, y2), t, list(1, 0), 0, 1.5, 1e-10, dopri5)");
        close(&got, &[1.5_f64.cos(), -1.5_f64.sin()], 1e-7);
        // An explicit dependence on the independent variable: y' = t, y(0) = 0.
        for method in ["dopri5", "exprb32"] {
            let got = run(&format!("ode_solve(t, y, t, 0, 0, 2, 1e-9, {method})"));
            close(&got, &[2.0], 1e-6);
        }
    }

    #[test]
    fn stiff_integrators_follow_the_slow_manifold() {
        // y' = -k (y - cos t): y(1) = (k^2 cos 1 + k sin 1 - k^2 e^-k) / (k^2 + 1).
        let k = 1.0e4_f64;
        let exact =
            (k * k * 1.0_f64.cos() + k * 1.0_f64.sin() - k * k * (-k).exp()) / (k * k + 1.0);
        for method in ["radau5", "bdf", "sdirk3"] {
            let got = run(&format!(
                "ode_solve(-10000*(y - cos(t)), y, t, 0, 0, 1, 1e-8, {method})"
            ));
            close(&got, &[exact], 1e-5);
        }
    }

    #[test]
    fn symplectic_composition_conserves_the_oscillator() {
        // q'' = -q from (1, 0) over half a period is (-1, 0).
        for (method, steps, tol) in [
            ("verlet", 4000, 1e-6),
            ("yoshida4", 200, 1e-7),
            ("yoshida6", 100, 1e-8),
            ("yoshida8", 60, 1e-8),
        ] {
            let got = run(&format!(
                "ode_symplectic(-q, q, list(1), list(0), 0, 3.141592653589793, {steps}, {method})"
            ));
            close(&got, &[-1.0, 0.0], tol);
        }
        // Two degrees of freedom, different frequencies: q1'' = -q1, q2'' = -4 q2.
        let got = run(
            "ode_symplectic(list(-q1, -4*q2), list(q1, q2), list(1, 0), list(0, 1), 0, 1, 400, yoshida6)",
        );
        close(
            &got,
            &[
                1.0_f64.cos(),
                0.5 * 2.0_f64.sin(),
                -1.0_f64.sin(),
                2.0_f64.cos(),
            ],
            1e-8,
        );
    }

    #[test]
    fn imex_and_exponential_integrators_split_the_stiffness() {
        // y' = -y (explicit) - 2 y (implicit): y(1) = e^-3.
        let want = (-3.0_f64).exp();
        for (scheme, tol) in [("euler", 2e-2), ("ars222", 1e-4), ("ars343", 1e-6)] {
            let got = run(&format!("ode_imex(-y, -2*y, y, t, 1, 0, 1, 200, {scheme})"));
            close(&got, &[want], tol);
        }
        // u' = L u + f with a rotation L and f = -u / 2: e^{-t/2} times the rotation.
        let got = run(
            "ode_etdrk4(list(list(0, 1), list(-1, 0)), list(-0.5*u1, -0.5*u2), list(u1, u2), t, list(1, 0), 0, 1, 50)",
        );
        let damping = (-0.5_f64).exp();
        close(
            &got,
            &[damping * 1.0_f64.cos(), -damping * 1.0_f64.sin()],
            1e-7,
        );
    }

    #[test]
    fn invalid_problems_are_left_as_requests() {
        for src in [
            "ode_solve(-2*y, y, t, 1, 0, 1, 1e-9, nonesuch)",
            "ode_solve(-2*y, y, t, 1, 0, 1, 0.5, dopri5)",
            "ode_solve(-2*a*y, y, t, 1, 0, 1, 1e-9, dopri5)",
            "ode_solve(-2*y, y, t, list(1, 2), 0, 1, 1e-9, dopri5)",
            "ode_solve(-2*y, y, y, 1, 0, 1, 1e-9, dopri5)",
            "ode_solve(y^2, y, t, 1, 0, 2, 1e-9, dopri5)",
            "ode_symplectic(-q, q, list(1), list(0), 0, 1, 0, verlet)",
            "ode_symplectic(-q, q, list(1), list(0), 0, 1, 10, nonesuch)",
            "ode_imex(-y, -2*y, y, t, 1, 0, 1, 10, nonesuch)",
            "ode_etdrk4(list(list(0, 1)), -u, u, t, 1, 0, 1, 10)",
        ] {
            let (text, reduced) = reduce_with(&[numeric()], src, &[]);
            assert!(!reduced, "{text}");
        }
    }
}
