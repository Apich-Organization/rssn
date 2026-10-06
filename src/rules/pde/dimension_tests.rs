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
    show(&format!("pdsolve({heat}, u(x, t), list(u(0, t) = 0, u(1, t) = 0, u(x, 0) = x*(1 - x)))"));
    show("pdsolve(diff(diff(u(x, t), t), t) + 2*diff(u(x, t), t) = diff(diff(u(x, t), x), x), u(x, t), list(u(0, t) = 0, u(pi, t) = 0, u(x, 0) = sin(4*x), at(diff(u(x, t), t), t, 0) = 0))");
    show("pdsolve(diff(diff(u(x, t), t), t) = diff(diff(u(x, t), x), x) - u(x, t), u(x, t), list(u(0, t) = 0, u(pi, t) = 0, u(x, 0) = sin(2*x), at(diff(u(x, t), t), t, 0) = 0))");
    show("pdsolve(I*diff(u(x, t), t) = -diff(diff(u(x, t), x), x), u(x, t), list(u(0, t) = 0, u(pi, t) = 0, u(x, 0) = sin(x)))");
    show("pdsolve(diff(diff(u(x, y, t), t), t) = diff(diff(u(x, y, t), x), x) + diff(diff(u(x, y, t), y), y), u(x, y, t), list(u(0, y, t) = 0, u(pi, y, t) = 0, u(x, 0, t) = 0, u(x, pi, t) = 0, u(x, y, 0) = sin(x)*sin(2*y), at(diff(u(x, y, t), t), t, 0) = 0))");
    show("pdsolve(diff(diff(u(x, y), x), x) + diff(diff(u(x, y), y), y) = 1, u(x, y), list(u(0, y) = 0, u(1, y) = 0, u(x, 0) = 0, u(x, 1) = 0))");
    show("pdsolve(diff(diff(u(x, y), x), x) + diff(diff(u(x, y), y), y) = 0, u(x, y), list(u(0, y) = y, u(1, y) = 0, u(x, 0) = 0, u(x, 1) = sin(pi*x)))");
    show("pdsolve(diff(diff(u(x, y), x), x) + diff(diff(u(x, y), y), y) + 5*u(x, y) = sin(x)*sin(y), u(x, y), list(u(0, y) = 0, u(pi, y) = 0, u(x, 0) = 0, u(x, pi) = 0))");
    show("pdsolve(diff(diff(u(x, t), t), t) + diff(diff(diff(diff(u(x, t), x), x), x), x) = 0, u(x, t), list(u(0, t) = 0, u(pi, t) = 0, at(diff(diff(u(x, t), x), x), x, 0) = 0, at(diff(diff(u(x, t), x), x), x, pi) = 0, u(x, 0) = sin(2*x), at(diff(u(x, t), t), t, 0) = 0))");
    show("pdsolve(diff(u(x, t), t) = diff(diff(u(x, t), x), x), u(x, t), list(u(0, t) = sin(t), u(1, t) = 0, u(x, 0) = 0))");
    // curvilinear
    show("pdsolve(diff(u(r, t), t) = diff(diff(u(r, t), r), r) + diff(u(r, t), r)/r, u(r, t), list(u(1, t) = 0, u(r, 0) = besselj(0, bessel_zero(0, 1)*r)))");
    show("pdsolve(diff(u(r, t), t) = diff(diff(u(r, t), r), r) + diff(u(r, t), r)/r, u(r, t), list(u(1, t) = 0, u(r, 0) = 1 - r^2))");
    show("pdsolve(diff(u(r, t), t) = diff(diff(u(r, t), r), r) + 2*diff(u(r, t), r)/r, u(r, t), list(u(1, t) = 0, u(r, 0) = 1 - r^2))");
    show("pdsolve(diff(u(r, t), t) = diff(diff(u(r, t), r), r) + 2*diff(u(r, t), r)/r, u(r, t), list(u(pi, t) = 0, u(r, 0) = sin(r)/r))");
    show("pdsolve(diff(u(r, th, t), t) = diff(diff(u(r, th, t), r), r) + diff(u(r, th, t), r)/r + diff(diff(u(r, th, t), th), th)/r^2, u(r, th, t), list(u(1, th, t) = 0, u(r, th, 0) = r*cos(th)*(1 - r^2)))");
    show("pdsolve(diff(diff(u(r, z), r), r) + diff(u(r, z), r)/r + diff(diff(u(r, z), z), z) = 0, u(r, z), list(u(1, z) = 0, u(r, 0) = 0, u(r, 1) = 1 - r^2))");
    show("pdsolve(diff(diff(u(r, th), r), r) + diff(u(r, th), r)/r + diff(diff(u(r, th), th), th)/r^2 = 0, u(r, th), list(u(2, th) = cos(3*th)))");
    show("pdsolve(diff(diff(u(r, th), r), r) + diff(u(r, th), r)/r + diff(diff(u(r, th), th), th)/r^2 = 1, u(r, th), list(u(1, th) = 0))");
    // diffusion
    show("pdsolve(diff(u(x, t), t) = diff(diff(u(x, t), x), x), u(x, t), list(u(0, t) = 1, u(x, 0) = 0))");
    show("pdsolve(diff(u(x, t), t) = diff(diff(u(x, t), x), x), u(x, t), list(at(diff(u(x, t), x), x, 0) = -1, u(x, 0) = 0))");
    show("pdsolve(diff(u(x, y, t), t) = diff(diff(u(x, y, t), x), x) + diff(diff(u(x, y, t), y), y), u(x, y, t), list(u(x, y, 0) = exp(-x^2 - y^2), u(0, y, t) = 0, u(x, 0, t) = 0))");
    show("pdsolve(diff(u(x, y, t), t) = diff(diff(u(x, y, t), x), x) + diff(diff(u(x, y, t), y), y) + 1, u(x, y, t), list(u(x, y, 0) = 0))");
}
