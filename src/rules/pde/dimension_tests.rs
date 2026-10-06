//! Tests of the dimension-uniform solvers.

use crate::rules::testing::simplify;

fn show(src: &str) {
    let rules = crate::rules::standard();
    let (text, ok) = crate::rules::testing::reduce_with(&rules, src, &[]);
    println!("{src}\n  => [{ok}] {text}\n");
}

#[test]
fn explore() {
    let _ = simplify;
    let heat = "diff(u(x, t), t) = diff(diff(u(x, t), x), x)";
    show(&format!("pdsolve({heat}, u(x, t), list(u(0, t) = u(2*pi, t), u(x, 0) = 1 + cos(x)))"));
    show("pdsolve(diff(diff(u(x, t), t), t) + diff(diff(diff(diff(u(x, t), x), x), x), x) = 0, u(x, t), list(u(0, t) = 0, u(pi, t) = 0, at(diff(diff(u(x, t), x), x), x, 0) = 0, at(diff(diff(u(x, t), x), x), x, pi) = 0, u(x, 0) = sin(2*x), at(diff(u(x, t), t), t, 0) = 0))");
    show("pdsolve(diff(diff(u(x, y), x), x) + diff(diff(u(x, y), y), y) = 1, u(x, y), list(u(0, y) = 0, u(1, y) = 0, u(x, 0) = 0, u(x, 1) = 0))");
    show("pdsolve(diff(u(r, t), t) = diff(diff(u(r, t), r), r) + diff(u(r, t), r)/r, u(r, t), list(u(1, t) = 0, u(r, 0) = besselj(0, bessel_zero(0, 1)*r)))");
    show("pdsolve(diff(u(r, t), t) = diff(diff(u(r, t), r), r) + 2*diff(u(r, t), r)/r, u(r, t), list(u(pi, t) = 0, u(r, 0) = sin(r)/r))");
    // systems
    show("pdsolve(list(diff(u(x, t), t) + diff(v(x, t), x) = 0, diff(v(x, t), t) + diff(u(x, t), x) = 0), list(u(x, t), v(x, t)))");
    show("pdsolve(list(diff(u(x, t), t) + diff(v(x, t), x) = 0, diff(v(x, t), t) + diff(u(x, t), x) = 0), list(u(x, t), v(x, t)), list(u(x, 0) = sin(x), v(x, 0) = 0))");
    show("pdsolve(list(diff(u(x, t), t) = 2*diff(diff(u(x, t), x), x) + diff(diff(v(x, t), x), x), diff(v(x, t), t) = diff(diff(u(x, t), x), x) + 2*diff(diff(v(x, t), x), x)), list(u(x, t), v(x, t)), list(u(0, t) = 0, u(pi, t) = 0, v(0, t) = 0, v(pi, t) = 0, u(x, 0) = sin(x), v(x, 0) = 0))");
    show("pdsolve(list(diff(diff(u(x, t), t), t) = 3*diff(diff(u(x, t), x), x) + diff(diff(v(x, t), x), x), diff(diff(v(x, t), t), t) = diff(diff(u(x, t), x), x) + 3*diff(diff(v(x, t), x), x)), list(u(x, t), v(x, t)), list(u(x, 0) = sin(x), v(x, 0) = 0, at(diff(u(x, t), t), t, 0) = 0, at(diff(v(x, t), t), t, 0) = 0))");
    // waves
    show("pdsolve(diff(diff(u(x, t), t), t) = diff(diff(u(x, t), x), x) + x, u(x, t), list(u(x, 0) = 0, at(diff(u(x, t), t), t, 0) = 0))");
    show("pdsolve(diff(diff(u(x, y, t), t), t) = diff(diff(u(x, y, t), x), x) + diff(diff(u(x, y, t), y), y), u(x, y, t), list(u(x, y, 0) = x^2 + y^2, at(diff(u(x, y, t), t), t, 0) = 0))");
    show("pdsolve(diff(diff(u(x, y, t), t), t) = diff(diff(u(x, y, t), x), x) + diff(diff(u(x, y, t), y), y), u(x, y, t), list(u(x, y, 0) = 0, at(diff(u(x, y, t), t), t, 0) = x*y))");
    show("pdsolve(diff(diff(u(x, t), t), t) + 2*diff(u(x, t), t) = diff(diff(u(x, t), x), x), u(x, t), list(u(x, 0) = exp(-x^2), at(diff(u(x, t), t), t, 0) = 0))");
    show("pdsolve(diff(diff(u(x, t), t), t) - diff(diff(u(x, t), x), x) + u(x, t) = 0, u(x, t), list(u(x, 0) = 0, at(diff(u(x, t), t), t, 0) = exp(-x^2)))");
    // elliptic half-space
    show("pdsolve(diff(diff(u(x, y), x), x) + diff(diff(u(x, y), y), y) = 0, u(x, y), list(u(x, 0) = 1/(1 + x^2)))");
    show("pdsolve(diff(diff(u(x, y), x), x) + diff(diff(u(x, y), y), y) = exp(-x^2 - y^2), u(x, y), list(u(x, 0) = 0))");
    // Poisson in the disk
    show("pdsolve(diff(diff(u(r, th), r), r) + diff(u(r, th), r)/r + diff(diff(u(r, th), th), th)/r^2 = 1, u(r, th), list(u(1, th) = 0))");
}
