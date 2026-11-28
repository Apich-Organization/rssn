//! Numerical modal solutions through the public symbolic PDE operator.

use rssn::graph::Budget;
use rssn::graph::ClosedForm;
use rssn::graph::Engine;
use rssn::graph::Env;
use rssn::graph::Extractor;
use rssn::graph::Graph;
use rssn::graph::NodeId;
use rssn::graph::Saturate;
use std::f64::consts::PI;

fn solve(source: &str) -> (Graph, NodeId) {
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &[rssn::rules::pde()]).unwrap();
    let query = graph.parse(source).unwrap();
    engine.run(
        &mut graph,
        &[query],
        &Env::symbolic(),
        &Saturate,
        &Budget::default(),
    );
    let answer = Extractor::new(&graph, &[query], &ClosedForm)
        .build(&mut graph, query)
        .unwrap();
    let text = graph.display(answer);
    assert!(!text.contains("defint"), "{text}");
    let &[_, rhs] = graph.children(answer) else {
        panic!("expected a solved equation, got {text}");
    };
    (graph, rhs)
}

fn evaluate(
    graph: &mut Graph,
    rhs: NodeId,
    bindings: &[(&str, f64)],
) -> f64 {
    let mut env = Env::numeric(0.0);
    for &(name, value) in bindings {
        let node = graph.parse(name).unwrap();
        env.bind(graph.symbol_of(node).unwrap(), value);
    }
    graph.eval(rhs, &env).unwrap()
}

fn integrate(f: impl Fn(f64) -> f64) -> f64 {
    let n = 2048;
    let step = 1.0 / f64::from(n);
    let mut total = f(0.0) + f(1.0);
    for i in 1..n {
        total += if i % 2 == 0 { 2.0 } else { 4.0 } * f(f64::from(i) * step);
    }
    total * step / 3.0
}

fn positive_root(
    equation: impl Fn(f64) -> f64,
    index: u32,
) -> f64 {
    let mut lo = f64::from(index) * PI;
    let mut hi = (f64::from(index) + 0.5) * PI;
    let first_sign = equation(lo).is_sign_positive();
    assert_ne!(first_sign, equation(hi).is_sign_positive());
    for _ in 0..60 {
        let mid = f64::midpoint(lo, hi);
        if equation(mid).is_sign_positive() == first_sign {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    f64::midpoint(lo, hi)
}

fn spherical_j0(x: f64) -> f64 {
    if x > 0.0 { x.sin() / x } else { 1.0 }
}

fn spherical_j1(x: f64) -> f64 {
    if x > 0.0 {
        (x.sin() / x - x.cos()) / x
    } else {
        0.0
    }
}

fn close(
    got: f64,
    expected: f64,
) {
    assert!(
        (got - expected).abs() < 3e-8,
        "got {got}, expected {expected}"
    );
}

#[test]
fn general_radial_neumann_data_on_disks_and_balls() {
    use rssn::kernels::special::bessel_j0;
    use rssn::kernels::special::bessel_j1;

    for (dimension, coefficient) in [(2_i32, "1"), (3_i32, "2")] {
        let source = format!(
            "pdsolve(diff(u(r,t),t) = diff(diff(u(r,t),r),r) + {coefficient}*diff(u(r,t),r)/r, u(r,t), list(at(diff(u(r,t),r),r,1) = 0, u(r,0) = 1.25 + cos(pi*r^2)))"
        );
        let (mut graph, rhs) = solve(&source);
        let initial = |r: f64| 1.25 + (PI * r * r).cos();
        let weight = |r: f64| r.powi(dimension - 1);
        let mean = f64::from(dimension) * integrate(|r| weight(r) * initial(r));
        let (r, t) = (0.42_f64, 0.055_f64);
        let mut expected = mean;
        for index in 1..=8 {
            let (root, mode): (f64, fn(f64) -> f64) = if dimension == 2 {
                (positive_root(bessel_j1, index), bessel_j0)
            } else {
                (
                    positive_root(|x| x * x.cos() - x.sin(), index),
                    spherical_j0,
                )
            };
            let c = integrate(|q| weight(q) * initial(q) * mode(root * q))
                / integrate(|q| weight(q) * mode(root * q).powi(2));
            expected += c * (-root * root * t).exp() * mode(root * r);
        }
        close(evaluate(&mut graph, rhs, &[("r", r), ("t", t)]), expected);
        close(evaluate(&mut graph, rhs, &[("r", r), ("t", 6.0)]), mean);
    }
}

#[test]
fn ball_heat_preserves_non_axisymmetric_angular_projection() {
    let (mut graph, rhs) = solve(concat!(
        "pdsolve(diff(u(r,th,ph,t),t) = diff(diff(u(r,th,ph,t),r),r)",
        " + 2*diff(u(r,th,ph,t),r)/r + diff(diff(u(r,th,ph,t),th),th)/r^2",
        " + cos(th)*diff(u(r,th,ph,t),th)/(r^2*sin(th))",
        " + diff(diff(u(r,th,ph,t),ph),ph)/(r^2*sin(th)^2), u(r,th,ph,t),",
        " list(u(1,th,ph,t) = 0, u(r,th,ph,0) = cos(pi*r^2)*(2 + sin(th)*(cos(ph) + 2*sin(ph)))))"
    ));
    let (r, theta, phi, t) = (0.42_f64, 1.1_f64, 0.8_f64, 0.055_f64);
    let angular = theta.sin() * (phi.cos() + 2.0 * phi.sin());
    let mut expected = 0.0;
    for index in 1..=8 {
        for (root, mode, scale) in [
            (f64::from(index) * PI, spherical_j0 as fn(f64) -> f64, 2.0),
            (
                positive_root(|x| x * x.cos() - x.sin(), index),
                spherical_j1,
                angular,
            ),
        ] {
            let c = integrate(|q| q * q * (PI * q * q).cos() * mode(root * q))
                / integrate(|q| q * q * mode(root * q).powi(2));
            expected += scale * c * (-root * root * t).exp() * mode(root * r);
        }
    }
    close(
        evaluate(
            &mut graph,
            rhs,
            &[("r", r), ("th", theta), ("ph", phi), ("t", t)],
        ),
        expected,
    );
}
