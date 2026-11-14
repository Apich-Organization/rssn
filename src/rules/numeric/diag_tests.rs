//! Tests of the pure-Neumann Poisson solves (`poisson_fd_bc`) and of the
//! diagnostics channel that explains why a numeric call did not reduce.

use std::f64::consts::PI;

use super::explain;
use super::numeric;
use crate::api::Session;
use crate::graph::Budget;
use crate::graph::Engine;
use crate::graph::Env;
use crate::graph::Graph;
use crate::graph::Saturate;
use crate::kernel::call::CallCtx;
use crate::kernel::diag::Diag;
use crate::kernel::diag::DiagKind;
use crate::kernel::signature::TypeExpr;
use crate::kernel::value::TensorValue;
use crate::kernel::value::Value;
use crate::rules::testing::reduce_with;
use crate::rules::testing::simplify;

fn run(src: &str) -> String {
    simplify(&[numeric()], src)
}

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

/// A function sampled on the unit box with `dims` nodes per axis, as nested
/// list text and as the flat row-major values.
fn sample(
    dims: &[usize],
    f: &dyn Fn(&[f64]) -> f64,
) -> (String, Vec<f64>) {
    fn go(
        dims: &[usize],
        at: &mut Vec<f64>,
        f: &dyn Fn(&[f64]) -> f64,
        flat: &mut Vec<f64>,
    ) -> String {
        let depth = at.len();
        let n = dims[depth];
        let items: Vec<String> = (0..n)
            .map(|i| {
                at.push(i as f64 / (n - 1) as f64);
                let item = if depth + 1 == dims.len() {
                    let v = f(at);
                    flat.push(v);
                    format!("{v:?}")
                } else {
                    go(dims, at, f, flat)
                };
                at.pop();
                item
            })
            .collect();
        format!("list({})", items.join(", "))
    }
    let mut flat = Vec::new();
    let text = go(dims, &mut Vec::new(), f, &mut flat);
    (text, flat)
}

/// The eigenvalue of the second difference with ghost-node Neumann ends of the
/// cosine mode `k` on `n` nodes.
fn neumann_eigenvalue(
    k: f64,
    n: usize,
) -> f64 {
    let h = 1.0 / (n - 1) as f64;
    (2.0 - 2.0 * (PI * k * h).cos()) / (h * h)
}

/// The eigenvalue of the Dirichlet second difference of the sine mode on `n` nodes.
fn dirichlet_eigenvalue(n: usize) -> f64 {
    let h = 1.0 / (n - 1) as f64;
    (2.0 - 2.0 * (PI * h).cos()) / (h * h)
}

/// The trapezoid-weighted mean of values on a 1D grid.
fn weighted_mean(values: &[f64]) -> f64 {
    let n = values.len();
    let w = |i: usize| {
        if i == 0 || i + 1 == n {
            0.5
        } else {
            1.0
        }
    };
    let total: f64 = (0..n).map(w).sum();
    (0..n).map(|i| w(i) * values[i]).sum::<f64>() / total
}

#[test]
fn neumann_cosine_modes_are_inverted_exactly_in_one_dimension() {
    for (n, k) in [(9_usize, 1.0_f64), (9, 3.0), (17, 2.0)] {
        let (text, f) = sample(&[n], &|p| (PI * k * p[0]).cos());
        let out = flat(&run(&format!("poisson_fd_bc({text}, neumann, 0)")));
        let want: Vec<f64> = f.iter().map(|v| v / neumann_eigenvalue(k, n)).collect();
        close(&out, &want, 1e-9);
        // The answer is the zero-mean one.
        assert!(weighted_mean(&out).abs() < 1e-10, "{out:?}");
    }
}

#[test]
fn neumann_cosine_modes_are_inverted_exactly_in_two_and_three_dimensions() {
    let n = 7;
    let (text, f) = sample(&[n, n], &|p| (PI * p[0]).cos() * (2.0 * PI * p[1]).cos());
    let out = flat(&run(&format!("poisson_fd_bc({text}, neumann, 0)")));
    let lambda = neumann_eigenvalue(1.0, n) + neumann_eigenvalue(2.0, n);
    let want: Vec<f64> = f.iter().map(|v| v / lambda).collect();
    close(&out, &want, 1e-9);
    let n = 5;
    let (text, f) = sample(&[n, n, n], &|p| {
        (PI * p[0]).cos() * (PI * p[1]).cos() * (2.0 * PI * p[2]).cos()
    });
    let out = flat(&run(&format!("poisson_fd_bc({text}, neumann, 0)")));
    let lambda = 2.0 * neumann_eigenvalue(1.0, n) + neumann_eigenvalue(2.0, n);
    let want: Vec<f64> = f.iter().map(|v| v / lambda).collect();
    close(&out, &want, 1e-9);
}

/// A source that is not a cosine series: a constant balanced by boundary
/// fluxes. The zero-mean shift of the parabola is checked against the
/// trapezoid mean.
#[test]
fn neumann_fluxes_balance_a_constant_source() {
    let n = 9;
    let c = 2.0;
    let (text, _) = sample(&[n], &|_| c);
    // -u'' = c with du/dn = -c/2 at both ends: u = c (x - x^2) / 2 + const.
    let out = flat(&run(&format!(
        "poisson_fd_bc({text}, neumann, list({}, {}))",
        -c / 2.0,
        -c / 2.0
    )));
    let shape: Vec<f64> = (0..n)
        .map(|i| {
            let x = i as f64 / (n - 1) as f64;
            c * (x - x * x) / 2.0
        })
        .collect();
    let mean = weighted_mean(&shape);
    let want: Vec<f64> = shape.iter().map(|v| v - mean).collect();
    close(&out, &want, 1e-9);
    assert!(weighted_mean(&out).abs() < 1e-10, "{out:?}");
}

#[test]
fn mixed_conditions_solve_in_every_dimension() {
    // 1D, one letter per face: Dirichlet 2 at x = 0, outward flux 1/2 at x = 1.
    // -u'' = 1: u = -x^2/2 + 3x/2 + 2, exact on the grid.
    let n = 9;
    let (text, _) = sample(&[n], &|_| 1.0);
    let out = flat(&run(&format!("poisson_fd_bc({text}, dn, list(2, 0.5))")));
    let want: Vec<f64> = (0..n)
        .map(|i| {
            let x = i as f64 / (n - 1) as f64;
            -x * x / 2.0 + 1.5 * x + 2.0
        })
        .collect();
    close(&out, &want, 1e-9);
    // 2D: Dirichlet on the x faces, Neumann on the y faces; a source that does
    // not depend on y has the y-independent solution sin(pi x) / lambda. The
    // per-axis and per-face spellings agree.
    let (text, f) = sample(&[n, n], &|p| (PI * p[0]).sin());
    let lambda = dirichlet_eigenvalue(n);
    let want: Vec<f64> = f.iter().map(|v| v / lambda).collect();
    let by_axis = run(&format!("poisson_fd_bc({text}, dn, 0)"));
    let by_face = run(&format!("poisson_fd_bc({text}, ddnn, 0)"));
    close(&flat(&by_axis), &want, 1e-9);
    assert_eq!(by_axis, by_face);
    // 3D: Dirichlet in x only.
    let m = 5;
    let (text, f) = sample(&[m, m, m], &|p| (PI * p[0]).sin());
    let out = flat(&run(&format!("poisson_fd_bc({text}, dnn, 0)")));
    let want: Vec<f64> = f.iter().map(|v| v / dirichlet_eigenvalue(m)).collect();
    close(&out, &want, 1e-9);
    // All Dirichlet is the plain operator.
    let (text, _) = sample(&[n], &|_| 2.0);
    assert_eq!(
        run(&format!("poisson_fd_bc({text}, dirichlet, 1)")),
        run(&format!("poisson_fd({text}, 1)"))
    );
}

fn why(src: &str) -> Option<Diag> {
    let mut g = Graph::new();
    let engine = Engine::install(&mut g, &[numeric()]).unwrap_or_else(|e| panic!("{e}"));
    let root = g
        .parse(src)
        .unwrap_or_else(|e| panic!("cannot parse `{src}`: {e}"));
    engine.run(
        &mut g,
        &[root],
        &Env::symbolic(),
        &Saturate,
        &Budget::default(),
    );
    explain(&g, root)
}

/// The finding for a call that must not reduce.
fn refused(src: &str) -> Diag {
    let (text, reduced) = reduce_with(&[numeric()], src, &[]);
    assert!(!reduced, "`{src}` reduced to {text}");
    why(src).unwrap_or_else(|| panic!("`{src}` gave no reason"))
}

#[test]
fn incompatible_neumann_data_say_so() {
    let diag = refused("poisson_fd_bc(list(1, 1, 1, 1, 1), neumann, 0)");
    assert_eq!(diag.kind, DiagKind::Incompatible, "{diag}");
    assert_eq!(diag.op, "poisson_fd_bc");
    assert!(diag.message.contains("imbalance"), "{diag}");
    // A source that integrates to zero against zero flux is accepted.
    let ok = run("poisson_fd_bc(list(1, -1, 0, 1, -1), neumann, 0)");
    assert_eq!(flat(&ok).len(), 5);
    // Wrong fluxes for a balanced source are incompatible again.
    let diag = refused("poisson_fd_bc(list(1, -1, 0, 1, -1), neumann, 0.5)");
    assert_eq!(diag.kind, DiagKind::Incompatible, "{diag}");
}

#[test]
fn boundary_specifications_are_checked() {
    let diag = refused("poisson_fd_bc(list(1, 2, 3, 4), robin, 0)");
    assert_eq!((diag.kind, diag.arg), (DiagKind::Domain, Some(1)), "{diag}");
    assert!(diag.message.contains("robin"), "{diag}");
    // Three letters for a 2D field are neither per axis (2) nor per face (4).
    let diag = refused("poisson_fd_bc(list(list(1, 2, 3), list(1, 2, 3), list(1, 2, 3)), dnd, 0)");
    assert_eq!((diag.kind, diag.arg), (DiagKind::Domain, Some(1)), "{diag}");
    // Three boundary values for a 1D field.
    let diag = refused("poisson_fd_bc(list(1, 2, 3, 4), dirichlet, list(0, 0, 0))");
    assert_eq!((diag.kind, diag.arg), (DiagKind::Shape, Some(2)), "{diag}");
}

#[test]
fn numeric_operators_name_the_argument_at_fault() {
    // tensordot: the paired axes have extents 3 and 2.
    let diag = refused(
        "tensordot(list(list(1, 2, 3), list(4, 5, 6)), list(list(1, 2), list(3, 4)), list(1), list(0))",
    );
    assert_eq!(diag.kind, DiagKind::Shape, "{diag}");
    assert!(diag.message.contains("equal extents"), "{diag}");
    let diag = refused(
        "tensordot(list(list(1, 2), list(3, 4)), list(list(1, 2), list(3, 4)), list(5), list(0))",
    );
    assert_eq!((diag.kind, diag.arg), (DiagKind::Domain, Some(2)), "{diag}");
    // tensor_contract.
    let diag = refused("tensor_contract(list(list(1, 2, 3), list(4, 5, 6)), 0, 1)");
    assert_eq!((diag.kind, diag.arg), (DiagKind::Shape, Some(1)), "{diag}");
    let diag = refused("tensor_contract(list(list(1, 2), list(3, 4)), 0, 0)");
    assert_eq!((diag.kind, diag.arg), (DiagKind::Domain, Some(2)), "{diag}");
    // Decompositions.
    let diag = refused("hosvd(list(list(1, 2), list(3, 4)), list(3, 1))");
    assert_eq!((diag.kind, diag.arg), (DiagKind::Domain, Some(1)), "{diag}");
    assert!(diag.message.contains("mode 0"), "{diag}");
    let diag = refused("cp_als(list(list(1, 2), list(3, 4)), 0, 10, 1)");
    assert_eq!((diag.kind, diag.arg), (DiagKind::Domain, Some(1)), "{diag}");
    let diag = refused("tt_svd(list(list(1, 2), list(3, 4)), 1.5, 0)");
    assert_eq!((diag.kind, diag.arg), (DiagKind::Domain, Some(1)), "{diag}");
    let diag = refused("tt_dense(list(list(1, 2)))");
    assert_eq!(diag.kind, DiagKind::Shape, "{diag}");
    // Wavelets.
    let diag = refused("dwt_approx(list(1, 2, 3), haar)");
    assert_eq!(diag.kind, DiagKind::Shape, "{diag}");
    assert!(diag.message.contains("even"), "{diag}");
    let diag = refused("dwt_approx(list(1, 2), nonesuch)");
    assert_eq!((diag.kind, diag.arg), (DiagKind::Domain, Some(1)), "{diag}");
    assert!(diag.message.contains("nonesuch"), "{diag}");
    // Polynomial systems.
    let diag = refused("polynomial_system_roots(list(list(list(1, 0))))");
    assert_eq!(diag.kind, DiagKind::Domain, "{diag}");
    assert!(diag.message.contains("constant"), "{diag}");
}

#[test]
fn solvers_distinguish_bad_input_from_failure_to_converge() {
    // A tolerance out of range is the argument's fault.
    let diag = refused("sparse_solve(list(list(0, 0, 1)), list(1), 2, fgmres)");
    assert_eq!((diag.kind, diag.arg), (DiagKind::Domain, Some(2)), "{diag}");
    let diag = refused("sparse_solve(list(list(0, 0, 1)), list(1), 1e-8, nonesuch)");
    assert_eq!((diag.kind, diag.arg), (DiagKind::Domain, Some(3)), "{diag}");
    // An index outside the order of the system.
    let diag = refused("sparse_solve(list(list(0, 5, 1)), list(1), 1e-8, minres)");
    assert_eq!(diag.kind, DiagKind::Domain, "{diag}");
    // The matrix does not match b.
    let diag = refused("minres_solve(list(list(1, 0), list(0, 1)), list(1, 2, 3), 1e-8)");
    assert_eq!((diag.kind, diag.arg), (DiagKind::Shape, Some(0)), "{diag}");
    // A singular system whose right-hand side is outside the range does not converge.
    let diag = refused("minres_solve(list(list(1, 0), list(0, 0)), list(1, 1), 1e-10)");
    assert_eq!(diag.kind, DiagKind::NotConverged, "{diag}");
    // A multigrid solve needs a positive diagonal.
    let diag = refused("amg_solve(list(list(0, 0, -1)), list(1), 1e-8, ruge_stuben)");
    assert_eq!(diag.kind, DiagKind::Domain, "{diag}");
}

#[test]
fn finite_difference_failures_carry_their_reason() {
    // A symbol where the field belongs, and a field with an unevaluated entry.
    let diag = refused("heat_fd(0.1, a, 0, 0.01, 4)");
    assert_eq!((diag.kind, diag.arg), (DiagKind::Shape, Some(1)), "{diag}");
    let diag = refused("heat_fd(0.1, list(1, 2, a + 1), 0, 0.01, 4)");
    assert_eq!(
        (diag.kind, diag.arg),
        (DiagKind::Symbolic, Some(1)),
        "{diag}"
    );
    // Too few nodes.
    let diag = refused("poisson_fd(list(1, 2), 0)");
    assert_eq!((diag.kind, diag.arg), (DiagKind::Domain, Some(0)), "{diag}");
    // A step count that is not a positive integer.
    let diag = refused("heat_fd(1, list(1, 2, 3), 0, 1, 0)");
    assert_eq!((diag.kind, diag.arg), (DiagKind::Domain, Some(4)), "{diag}");
    // The Courant number of this wave is 4.
    let diag = refused("wave_fd(1, list(0, 0, 0, 0, 0), list(0, 0, 0, 0, 0), 0, 1, 1)");
    assert_eq!(diag.kind, DiagKind::Domain, "{diag}");
    assert!(diag.message.contains("Courant"), "{diag}");
    // Fields that do not match.
    let diag = refused("wave_fd_bc(1, list(1, 1, 1), list(0, 0, 0, 0), neumann, 0, 0.01, 4)");
    assert_eq!((diag.kind, diag.arg), (DiagKind::Shape, Some(2)), "{diag}");
    // A velocity needs one component per space dimension.
    let diag = refused("advect_fd(list(1, 1), list(0, 1, 2, 3), 1, 4)");
    assert_eq!((diag.kind, diag.arg), (DiagKind::Shape, Some(0)), "{diag}");
}

#[test]
fn a_call_that_reduces_has_nothing_to_explain() {
    assert!(why("tensor_norm(list(list(3, 0), list(0, 4)))").is_none());
    assert!(why("poisson_fd(list(2, 2, 2, 2, 2), 0)").is_none());
    // An operator outside the table-driven families keeps no reasons.
    assert!(why("fftconv(list(1, 2), list(1))").is_none());
}

fn tensor(
    shape: &[usize],
    data: &[f64],
) -> Value {
    Value::Tensor(TensorValue::of(
        shape.to_vec(),
        data.iter().map(|&x| Value::Float(x)).collect(),
        TypeExpr::Real,
    ))
}

#[test]
fn the_session_reports_why_an_operator_did_not_reduce() {
    let session = Session::new();
    let ctx = CallCtx::new();
    // Success is the value `call_value` gives.
    let norm = tensor(&[2, 2], &[3.0, 0.0, 0.0, 4.0]);
    let diagnosed = session.apply_diag("tensor_norm", std::slice::from_ref(&norm), &ctx);
    let plain = session.call_value("tensor_norm", std::slice::from_ref(&norm), &ctx);
    assert_eq!(diagnosed.ok(), plain.ok());
    // The reason of a refusal, with the argument that caused it.
    let x = tensor(&[2, 3], &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
    let y = tensor(&[2, 2], &[1.0, 2.0, 3.0, 4.0]);
    let (one, zero) = (tensor(&[1], &[1.0]), tensor(&[1], &[0.0]));
    let args = [x, y, one, zero];
    let found = session.apply_diag("tensordot", &args, &ctx);
    let diag = found.expect_err("mismatched axes");
    assert_eq!(diag.op, "tensordot");
    assert_eq!(diag.kind, DiagKind::Shape);
    assert!(diag.message.contains("equal extents"), "{diag}");
    // Incompatible Neumann data through the typed path.
    let ones = tensor(&[5], &[1.0; 5]);
    let args = [ones, Value::Symbol("neumann".to_owned()), Value::Float(0.0)];
    let diag = session
        .apply_diag("poisson_fd_bc", &args, &ctx)
        .expect_err("incompatible data");
    assert_eq!(
        (diag.op.as_str(), diag.kind),
        ("poisson_fd_bc", DiagKind::Incompatible)
    );
    // A numeric target loses the term in a `NotNumeric` error; the reason is
    // still found.
    let numeric_ctx = CallCtx::new().numeric(1e-9);
    let ones = tensor(&[5], &[1.0; 5]);
    let args = [ones, Value::Symbol("neumann".to_owned()), Value::Float(0.0)];
    let diag = session
        .apply_diag("poisson_fd_bc", &args, &numeric_ctx)
        .expect_err("incompatible data");
    assert_eq!(diag.kind, DiagKind::Incompatible, "{diag}");
    // The same field with a matching flux reduces.
    let balanced = tensor(&[5], &[1.0, -1.0, 0.0, 1.0, -1.0]);
    let args = [
        balanced,
        Value::Symbol("neumann".to_owned()),
        Value::Float(0.0),
    ];
    assert!(session.apply_diag("poisson_fd_bc", &args, &ctx).is_ok());
    // A wrong argument count is an arity finding, an unknown name unsupported.
    let diag = session
        .apply_diag("tensor_norm", &[], &ctx)
        .expect_err("no arguments");
    assert_eq!(diag.kind, DiagKind::Arity, "{diag}");
    let diag = session
        .apply_diag("no_such_operator", &[Value::Float(1.0)], &ctx)
        .expect_err("unknown");
    assert_eq!(diag.kind, DiagKind::Unsupported, "{diag}");
    // `explain` works on a term that is already built.
    let term = session
        .parse("tensor_contract(list(list(1, 2, 3), list(4, 5, 6)), 0, 1)")
        .expect("parses");
    let diag = session.explain(term).expect("a reason");
    assert_eq!(diag.kind, DiagKind::Shape, "{diag}");
}
