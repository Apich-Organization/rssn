//! Numeric kernels that have no symbolic rule set of their own.
//!
//! These are the routines of [`crate::kernels`] that a host otherwise has to
//! reimplement: they take numbers and give numbers, so they belong in the
//! manifest as operators rather than staying library functions.
//!
//! | operator | value |
//! |---|---|
//! | `fftconv(a, b)` | the linear convolution of two sequences, `c[k] = sum_i a[i] b[k - i]`, of length `len(a) + len(b) - 1` |
//! | `interpolate(xs, ys, x)` | the value at `x` of the **Lagrange polynomial** through the points `(xs[i], ys[i])` |
//!
//! Both compute only from literal numbers. A symbolic argument leaves the
//! call as a request: there is no symbolic form of a convolution or of an
//! interpolating polynomial's value, and inventing one (or returning zeros)
//! would be a wrong answer. `interpolate` is polynomial interpolation, not
//! the piecewise-linear interpolant some standard libraries call `interp1`;
//! the two agree only at the nodes, so a host must not substitute one for the
//! other.
//!
//! The rest of the rule set is table driven and lives in the submodules:
//!
//! | module | operators |
//! |---|---|
//! | `bridge` | robust statistics, windows and spectra, one-level wavelets, finite-difference PDEs, dense Krylov solves, low-rank, tensors, sampled functions |
//! | `special` | `lgamma_complex`, `gamma_complex`, `digamma_complex`, `lambert_w_complex`, `bessel_j_complex`, `bessel_i_complex` (complex values as `re, im` arguments and `list(re, im)` results) |
//! | `resample` | `quantile_hf`, `bootstrap_ci` |
//! | `transform` | `wavedec`, `waverec`, `dwt2`, `idwt2`, `frft_re`, `frft_im` |
//! | `decomp` | `qr_pivoted_*`, `hosvd`, `hooi`, `cp_als`, `tt_svd`, `tt_dense` |
//! | `solvers` | `sparse_solve`, `amg_solve`, `nested_dissection_order`, `polynomial_system_roots` |
//! | `functions` | expression operators: `oscillatory_integral`, `anderson_fixed_point`, `newton_cg_minimize` |
//! | `ode` | expression operators: `ode_solve`, `ode_symplectic`, `ode_imex`, `ode_etdrk4` |
//!
//! `fft` itself returns complex values, and a complex tensor has no term
//! encoding ([`crate::api`] refuses a `Value::Complex` with a nonzero
//! imaginary part), so only its real, imaginary and absolute parts are
//! operators; `odeint` is the scalar operator of the `ode` rule set.

mod bridge;
mod decomp;
mod functions;
mod ode;
mod resample;
mod solvers;
mod special;
mod transform;

use crate::graph::Arity;
use crate::graph::Cx;
use crate::graph::Graph;
use crate::graph::Kernel;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::OpDescriptor;
use crate::graph::OpFlags;
use crate::graph::OpId;
use crate::graph::Outcome;
use crate::graph::RuleError;
use crate::graph::RuleSet;
use crate::graph::Tier;
use crate::graph::dim::DimRule;
use crate::graph::op::core;
use crate::graph::rule::Installer;
use crate::kernels::fftconv;
use crate::kernels::interpolate::lagrange_interpolation;
use crate::rules::calculus::calculus;

/// The numeric-kernel rule set.
#[must_use]
pub fn numeric() -> RuleSet {
    RuleSet::new("numeric", install).needs(calculus())
}

fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    let request = |name: &str, arity: u8| {
        OpDescriptor::new(name, Arity::Fixed(arity))
            .flags(OpFlags::HEAVY)
            .cost(100)
    };
    let fftconv = i.op(request("fftconv", 2))?;
    let interpolate = i.op(request("interpolate", 3))?;
    // A convolution multiplies the two signals' dimensions; interpolation
    // carries the ordinate's dimension through. Without these the operators
    // would be dimensionally polymorphic by omission, which the coverage test
    // refuses — correctly, because a silently dimensionless convolution is a
    // wrong answer rather than a missing one.
    i.dim(fftconv, DimRule::Product);
    i.dim(interpolate, DimRule::Same(1));
    // The table-driven numeric kernels (robust statistics, signal processing,
    // wavelets and the dimension-free finite-difference solvers).
    let tables: [&[bridge::Entry]; 6] = [
        bridge::ENTRIES,
        special::ENTRIES,
        resample::ENTRIES,
        transform::ENTRIES,
        decomp::ENTRIES,
        solvers::ENTRIES,
    ];
    let mut bridged = Vec::new();
    for entry in tables.into_iter().flatten() {
        let op = i.op(request(entry.name, entry.params.len() as u8))?;
        i.dim(op, bridge_dimension(entry.name));
        bridged.push((op, entry.eval));
    }
    i.kernel(
        "numeric/bridge",
        Tier::Reduce,
        bridge::Bridge { ops: bridged },
    );
    let heavy = |name: &str, arity: u8| request(name, arity);
    let oscillatory = i.op(heavy("oscillatory_integral", 7))?;
    let anderson = i.op(heavy("anderson_fixed_point", 4))?;
    let newton_cg = i.op(heavy("newton_cg_minimize", 3))?;
    // The unknowns carry dimensions independent of the expression's.
    for op in [oscillatory, anderson, newton_cg] {
        i.dim(op, DimRule::Any);
    }
    i.kernel(
        "numeric/functions",
        Tier::Reduce,
        functions::Functions {
            oscillatory,
            anderson,
            newton_cg,
        },
    );
    // Initial-value problems with a named integrator. The state keeps the
    // initial state's dimension; the symplectic pair mixes positions and
    // momenta, so it is a container.
    let ode_solve = i.op(heavy("ode_solve", 8))?;
    let ode_symplectic = i.op(heavy("ode_symplectic", 8))?;
    let ode_imex = i.op(heavy("ode_imex", 9))?;
    let ode_etdrk4 = i.op(heavy("ode_etdrk4", 8))?;
    i.dim(ode_solve, DimRule::Same(3));
    i.dim(ode_symplectic, DimRule::Container);
    i.dim(ode_imex, DimRule::Same(4));
    i.dim(ode_etdrk4, DimRule::Same(4));
    i.kernel(
        "numeric/ode",
        Tier::Reduce,
        ode::OdeOps {
            solve: ode_solve,
            symplectic: ode_symplectic,
            imex: ode_imex,
            etdrk4: ode_etdrk4,
        },
    );
    i.kernel("numeric/fftconv", Tier::Reduce, FftConv { op: fftconv });
    i.kernel(
        "numeric/interpolate",
        Tier::Reduce,
        Interpolate { op: interpolate },
    );
    Ok(())
}

/// The dimension rule of a bridged operator: location and scale estimates
/// and transforms of a signal keep the signal's dimension, windows are pure
/// numbers, a correlation multiplies the two signals' dimensions, and the
/// remaining operators (probabilities, densities, the finite-difference
/// solvers whose grid is the unit box) do not track one.
fn bridge_dimension(name: &str) -> DimRule {
    match name {
        // Windows and special functions of a complex argument take bare
        // numbers.
        | "hann_window" | "hamming_window" | "blackman_window" | "lgamma_complex"
        | "gamma_complex" | "digamma_complex" | "lambert_w_complex" | "bessel_j_complex"
        | "bessel_i_complex" => DimRule::Dimensionless,
        | "cross_correlation" | "tensor_inner" | "tensor_outer" | "tensordot" => DimRule::Product,
        // A probability is a pure number whatever the sample's dimension.
        | "ecdf" => DimRule::NonMetric,
        // A density is per unit of the sample's dimension.
        | "kde" => DimRule::Inverse,
        // The unit box has no length scale, so the field keeps its dimension.
        | "heat_fd" | "wave_fd" | "advect_fd" => DimRule::Same(1),
        | "poisson_fd" => DimRule::Same(0),
        // Ranks and column indices are counts, not quantities.
        | "qr_pivoted_rank" | "qr_pivoted_permutation" | "nested_dissection_order" => {
            DimRule::NonMetric
        },
        // A decomposition is a tuple of tensors of different roles.
        | "hosvd" | "hooi" | "cp_als" | "tt_svd" => DimRule::Container,
        // Documented bypasses (`graph::dim` tests): the unknown is `b / A`,
        // the roots' dimensions follow from coefficients that are bare numbers
        // here, or the result's dimension depends on an abscissa or variables
        // whose dimensions are independent of the operands'.
        | "sparse_solve"
        | "amg_solve"
        | "polynomial_system_roots"
        | "minres_solve"
        | "fgmres_solve"
        | "bicgstab_l_solve" => DimRule::Any,
        | _ if name.starts_with("function_") => DimRule::Any,
        | _ => DimRule::Same(0),
    }
}

/// The numbers a `list` term holds, in order.
fn numbers(
    graph: &Graph,
    node: NodeId,
) -> Option<Vec<f64>> {
    if graph.op(node) != core::LIST {
        return None;
    }
    graph
        .children(node)
        .iter()
        .map(|&item| graph.number_of(item).map(Number::to_f64))
        .collect()
}

/// The single number a term holds.
fn number(
    graph: &Graph,
    node: NodeId,
) -> Option<f64> {
    graph.number_of(node).map(Number::to_f64)
}

/// `fftconv(a, b)`: linear convolution of two real sequences.
struct FftConv {
    op: OpId,
}

impl Kernel for FftConv {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let args = cx.graph.children(node).to_vec();
        let &[a, b] = args.as_slice() else {
            return Outcome::Pass;
        };
        let (Some(a), Some(b)) = (numbers(cx.graph, a), numbers(cx.graph, b)) else {
            return Outcome::Pass;
        };
        let out = fftconv::convolve(&a, &b);
        let items: Vec<NodeId> = out.iter().map(|&value| cx.graph.float(value)).collect();
        Outcome::Equal(cx.graph.node(core::LIST, &items))
    }
}

/// `interpolate(xs, ys, x)`: the Lagrange polynomial through the points,
/// evaluated at `x`.
struct Interpolate {
    op: OpId,
}

impl Kernel for Interpolate {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let args = cx.graph.children(node).to_vec();
        let &[xs, ys, x] = args.as_slice() else {
            return Outcome::Pass;
        };
        let (Some(xs), Some(ys), Some(x)) = (
            numbers(cx.graph, xs),
            numbers(cx.graph, ys),
            number(cx.graph, x),
        ) else {
            return Outcome::Pass;
        };
        if xs.len() != ys.len() || xs.is_empty() {
            return Outcome::Pass;
        }
        let points: Vec<(f64, f64)> = xs.iter().copied().zip(ys).collect();
        // Duplicate abscissae make the interpolant undefined; the call stays
        // a request rather than picking one of the two ordinates.
        let Ok(polynomial) = lagrange_interpolation(&points) else {
            return Outcome::Pass;
        };
        Outcome::Equal(cx.graph.float(polynomial.eval(x)))
    }
}

#[cfg(test)]
mod tests {
    use crate::graph::Engine;
    use crate::graph::Graph;
    use crate::graph::RuleSet;
    use crate::kernel::SignatureSource;
    use crate::kernel::infer::signatures;
    use crate::kernel::overload::Best;
    use crate::kernel::overload::best_types;
    use crate::kernel::signature::NatExpr;
    use crate::kernel::signature::TypeExpr;
    use crate::rules::testing::reduce_with;
    use crate::rules::testing::simplify;

    use super::numeric;

    fn run(src: &str) -> String {
        simplify(&[numeric()], src)
    }

    fn vector() -> TypeExpr {
        TypeExpr::Tensor {
            rank: NatExpr::Lit(1),
            shape: None,
            elem: Box::new(TypeExpr::Real),
        }
    }

    #[test]
    fn convolution_is_the_linear_one() {
        // [1, 2, 3] * [0, 1, 1/2] = [0, 1, 5/2, 4, 3/2].
        assert_eq!(
            run("fftconv(list(1, 2, 3), list(0, 1, 0.5))"),
            "list(0, 1, 2.5, 4, 1.5)"
        );
        // The length is len(a) + len(b) - 1, not either input's.
        assert_eq!(run("fftconv(list(1), list(2))"), "list(2)");
        // The empty case is empty, not zeros of some length.
        assert_eq!(run("fftconv(list(), list(1, 2))"), "list()");
    }

    #[test]
    fn lagrange_interpolation_hits_the_nodes_and_the_midpoint() {
        // The parabola through (0, 0), (1, 1), (2, 4) is x^2.
        assert_eq!(
            run("interpolate(list(0, 1, 2), list(0, 1, 4), 3/2)"),
            "2.25"
        );
        assert_eq!(run("interpolate(list(0, 1, 2), list(0, 1, 4), 1)"), "1");
    }

    /// The numbers of a printed (possibly nested) list, in order.
    fn flat(text: &str) -> Vec<f64> {
        text.replace(['(', ')', ','], " ")
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
    fn robust_statistics_resist_an_outlier() {
        assert_eq!(run("mad(list(1, 2, 3, 4, 100))"), "1");
        assert!((flat(&run("mad_sigma(list(1, 2, 3, 4, 100))"))[0] - 1.4826).abs() < 1e-3);
        assert_eq!(run("trimmed_mean(list(1, 2, 3, 4, 100), 0.2)"), "3");
        assert_eq!(run("ecdf(list(1, 2, 3, 4), 2.5)"), "0.5");
        let huber = flat(&run("huber_location(list(1, 2, 3, 4, 100), 1.345)"))[0];
        assert!((2.5..4.5).contains(&huber), "{huber}");
        let h = flat(&run("bandwidth_silverman(list(1, 2, 3, 4, 5, 6, 7, 8))"))[0];
        assert!(h > 0.0 && h < 3.0);
        let density = flat(&run("kde(list(0, 0, 0), 0, 1)"))[0];
        assert!((density - 0.398_942_28).abs() < 1e-6);
        // Empty or symbolic samples stay requests.
        for src in [
            "mad(list())",
            "mad(list(1, a, 3))",
            "trimmed_mean(list(1, 2), x)",
        ] {
            let (text, reduced) = reduce_with(&[numeric()], src, &[]);
            assert!(!reduced, "{text}");
        }
    }

    #[test]
    fn windows_and_spectra() {
        close(
            &flat(&run("hann_window(5)")),
            &[0.0, 0.5, 1.0, 0.5, 0.0],
            1e-12,
        );
        assert_eq!(flat(&run("hamming_window(4)")).len(), 4);
        close(
            &flat(&run("fft_abs(list(1, 1, 1, 1))")),
            &[4.0, 0.0, 0.0, 0.0],
            1e-12,
        );
        close(&flat(&run("fft_re(list(1, 0, 0, 0))")), &[1.0; 4], 1e-12);
        close(
            &flat(&run(
                "ifft_re(fft_re(list(1, 2, 3)), fft_im(list(1, 2, 3)))",
            )),
            &[1.0, 2.0, 3.0],
            1e-12,
        );
        close(
            &flat(&run("cross_correlation(list(1, 2, 3), list(0, 1, 0.5))")),
            &[0.5, 2.0, 3.5, 3.0, 0.0],
            1e-12,
        );
        let (text, reduced) = reduce_with(&[numeric()], "hann_window(-1)", &[]);
        assert!(!reduced, "{text}");
    }

    #[test]
    fn wavelets_by_name_round_trip() {
        close(
            &flat(&run("dwt_approx(list(1, 1, 1, 1), haar)")),
            &[std::f64::consts::SQRT_2; 2],
            1e-12,
        );
        close(
            &flat(&run("dwt_detail(list(1, 1, 1, 1), haar)")),
            &[0.0; 2],
            1e-12,
        );
        for name in ["haar", "db2", "db4", "sym4", "bior2_2"] {
            let x = "list(3, 1, 4, 1, 5, 9, 2, 6)";
            close(
                &flat(&run(&format!(
                    "idwt(dwt_approx({x}, {name}), dwt_detail({x}, {name}), {name})"
                ))),
                &[3.0, 1.0, 4.0, 1.0, 5.0, 9.0, 2.0, 6.0],
                1e-9,
            );
        }
        // An odd signal or an unknown wavelet is not transformed.
        for src in [
            "dwt_approx(list(1, 2, 3), haar)",
            "dwt_approx(list(1, 2), nonesuch)",
        ] {
            let (text, reduced) = reduce_with(&[numeric()], src, &[]);
            assert!(!reduced, "{text}");
        }
    }

    /// One operator per equation whatever the dimension: the nesting depth of
    /// the data chooses 1D, 2D or 3D.
    #[test]
    fn finite_difference_operators_are_dimension_free() {
        // -u'' = 2, u(0) = u(1) = 0: u = x (1 - x), exact on the grid.
        close(
            &flat(&run("poisson_fd(list(2, 2, 2, 2, 2), 0)")),
            &[0.0, 0.1875, 0.25, 0.1875, 0.0],
            1e-12,
        );
        // The same operator on a square and a cube: constant data and
        // boundary value give the constant.
        let row = "list(0, 0, 0, 0, 0)";
        let square = format!("list({row}, {row}, {row}, {row}, {row})");
        let out = run(&format!("poisson_fd({square}, 1)"));
        assert_eq!(flat(&out).len(), 25);
        assert!(flat(&out).iter().all(|v| (v - 1.0).abs() < 1e-12), "{out}");
        assert_eq!(out.matches("list(").count(), 6, "{out}");
        let one = "list(1, 1, 1, 1, 1)";
        let layer = format!("list({one}, {one}, {one}, {one}, {one})");
        let cube = format!("list({layer}, {layer}, {layer}, {layer})");
        let out = run(&format!("heat_fd(0.1, {cube}, 1, 0.01, 4)"));
        assert_eq!(flat(&out).len(), 100);
        assert!(flat(&out).iter().all(|v| (v - 1.0).abs() < 1e-9), "{out}");
        // Waves and advection.
        let out = run("wave_fd(1, list(0, 0, 0, 0, 0), list(0, 0, 0, 0, 0), 0, 0.1, 10)");
        close(&flat(&out), &[0.0; 5], 1e-12);
        // A Courant number of one shifts one cell per step: a full turn is the identity.
        close(
            &flat(&run(
                "advect_fd(list(1), list(0, 1, 2, 3, 4, 5, 6, 7), 1, 8)",
            )),
            &[0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0],
            1e-9,
        );
        // Mismatched velocity length, ragged data and too few nodes stay requests.
        for src in [
            "advect_fd(list(1, 1), list(0, 1, 2, 3), 1, 4)",
            "poisson_fd(list(list(1, 2, 3), list(1, 2)), 0)",
            "poisson_fd(list(1, 2), 0)",
            "heat_fd(1, list(1, 2, 3), 0, 1, 0)",
        ] {
            let (text, reduced) = reduce_with(&[numeric()], src, &[]);
            assert!(!reduced, "{text}");
        }
    }

    #[test]
    fn linear_solvers_tensors_functions_and_low_rank() {
        close(
            &flat(&run(
                "minres_solve(list(list(2, 1), list(1, 3)), list(1, 2), 1e-10)",
            )),
            &[0.2, 0.6],
            1e-8,
        );
        close(
            &flat(&run(
                "fgmres_solve(list(list(2, 1), list(1, 3)), list(1, 2), 1e-10)",
            )),
            &[0.2, 0.6],
            1e-8,
        );
        assert_eq!(run("tensor_norm(list(list(3, 0), list(0, 4)))"), "5");
        assert_eq!(
            run("tensor_inner(list(list(3, 0), list(0, 4)), list(list(3, 0), list(0, 4)))"),
            "25"
        );
        close(
            &flat(&run(
                "tensordot(list(list(1, 2), list(3, 4)), list(list(0, 1), list(1, 0)), list(1), list(0))",
            )),
            &[2.0, 1.0, 4.0, 3.0],
            1e-12,
        );
        assert_eq!(
            run("function_l1_norm(list(list(0, 1), list(1, 1), list(2, 1)))"),
            "2"
        );
        close(
            &flat(&run(
                "randomized_singular_values(list(list(3, 0), list(0, 1)), 1)",
            )),
            &[3.0],
            1e-9,
        );
        close(
            &flat(&run(
                "lobpcg_eigenvalues(list(list(1, 0, 0), list(0, 2, 0), list(0, 0, 3)), 1, 1e-8)",
            )),
            &[1.0],
            1e-6,
        );
        // A non-square or asymmetric input is not solved.
        for src in [
            "minres_solve(list(list(1, 2, 3)), list(1), 1e-8)",
            "lobpcg_eigenvalues(list(list(1, 2, 0), list(0, 1, 0), list(0, 0, 1)), 1, 1e-8)",
        ] {
            let (text, reduced) = reduce_with(&[numeric()], src, &[]);
            assert!(!reduced, "{text}");
        }
    }

    #[test]
    fn expression_operators_integrate_iterate_and_minimise() {
        // int_0^1 exp(10 i x) dx = (exp(10 i) - 1) / (10 i).
        let got = flat(&run("oscillatory_integral(1, x, x, 0, 1, 10, 24)"));
        close(
            &got,
            &[10.0_f64.sin() / 10.0, (1.0 - 10.0_f64.cos()) / 10.0],
            1e-8,
        );
        // The Dottie number.
        close(
            &flat(&run("anderson_fixed_point(cos(x), x, 1, 1e-12)")),
            &[0.739_085_133_215_160_6],
            1e-9,
        );
        close(
            &flat(&run(
                "newton_cg_minimize((x - 1)^2 + 2*(y + 2)^2, list(x, y), list(0, 0))",
            )),
            &[1.0, -2.0],
            1e-6,
        );
        // An unbound parameter leaves the call a request.
        let (text, reduced) = reduce_with(
            &[numeric()],
            "anderson_fixed_point(cos(a*x), x, 1, 1e-12)",
            &[],
        );
        assert!(!reduced, "{text}");
    }

    #[test]
    fn a_symbolic_argument_is_left_as_a_request() {
        let (text, reduced) = reduce_with(
            &[numeric()],
            "interpolate(list(0, 1), list(a, b), 1/2)",
            &[],
        );
        assert!(!reduced, "{text}");
        assert_eq!(text, "interpolate(list(0, 1), list(a, b), 1/2)");
    }

    /// The tensor overload must win over the generic `Term` fallback.
    #[test]
    fn the_tensor_signature_is_the_one_that_resolves() {
        let mut graph = Graph::new();
        let sets: Vec<RuleSet> = vec![numeric()];
        assert!(Engine::install(&mut graph, &sets).is_ok());
        for (name, args) in [
            ("fftconv", vec![vector(), vector()]),
            ("interpolate", vec![vector(), vector(), TypeExpr::Real]),
        ] {
            let op = graph
                .ops()
                .lookup(name)
                .unwrap_or_else(|| panic!("{name} is registered"));
            let (overloads, source) = signatures(graph.ops().get(op));
            assert_eq!(source, SignatureSource::Overlay, "{name}");
            let Best::Unique(chosen) = best_types(&overloads, &args) else {
                panic!("{name}: {args:?} did not resolve to one overload");
            };
            assert_eq!(overloads[chosen].params, args, "{name}");
        }
    }

    /// The path a host takes: resolve the overload from values, run it, and
    /// read the answer back as a value.
    #[test]
    fn the_session_can_delegate_the_numeric_kernels() {
        use crate::api::Session;
        use crate::kernel::CallCtx;
        use crate::kernel::value::TensorValue;
        use crate::kernel::value::Value;

        let session = Session::new();
        let a = Value::Tensor(TensorValue::vector(vec![
            Value::int(1_i32),
            Value::int(2_i32),
            Value::int(3_i32),
        ]));
        let b = Value::Tensor(TensorValue::vector(vec![
            Value::int(0_i32),
            Value::int(1_i32),
            Value::Float(0.5),
        ]));
        let out = session
            .call_value("fftconv", &[a, b], &CallCtx::new())
            .unwrap_or_else(|e| panic!("{e}"));
        let Value::Tensor(convolution) = out else {
            panic!("expected a tensor, got {out:?}");
        };
        assert_eq!(convolution.shape, vec![5]);
        let values: Vec<f64> = convolution
            .data
            .iter()
            .map(|value| {
                match value {
                    | Value::Int(n) => n.to_string().parse().unwrap_or(f64::NAN),
                    | Value::Real(number) => number.to_f64(),
                    | other => panic!("expected a number, got {other:?}"),
                }
            })
            .collect();
        assert_eq!(values, vec![0.0, 1.0, 2.5, 4.0, 1.5]);

        let xs = Value::Tensor(TensorValue::vector(vec![
            Value::int(0_i32),
            Value::int(1_i32),
            Value::int(2_i32),
        ]));
        let ys = Value::Tensor(TensorValue::vector(vec![
            Value::int(0_i32),
            Value::int(1_i32),
            Value::int(4_i32),
        ]));
        let out = session
            .call_value("interpolate", &[xs, ys, Value::Float(1.5)], &CallCtx::new())
            .unwrap_or_else(|e| panic!("{e}"));
        let Value::Real(answer) = out else {
            panic!("expected a number, got {out:?}");
        };
        assert_eq!(answer.to_f64(), 2.25);
    }
}
