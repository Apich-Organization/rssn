//! Round-2 coverage: rigorous (interval) machinery and exact D-finite
//! functions -- `kernels::ball`, `kernels::homotopy_cert` (the
//! complex-coefficient entry points) and `kernels::holonomic`.
//! Public API only; every input is deterministic.

#![allow(clippy::needless_range_loop)]

use num_bigint::BigInt;
use num_complex::Complex64 as Cx;
use num_rational::BigRational;
use rssn::kernels::ball::CTerm;
use rssn::kernels::ball::Cball;
use rssn::kernels::ball::eval_jacobian;
use rssn::kernels::ball::eval_jacobian_terms;
use rssn::kernels::ball::eval_polys;
use rssn::kernels::ball::eval_terms;
use rssn::kernels::ball::jacobian_at;
use rssn::kernels::ball::jacobian_at_terms;
use rssn::kernels::ball::real_terms;
use rssn::kernels::ball::values_at;
use rssn::kernels::ball::values_at_terms;
use rssn::kernels::holonomic::Atom;
use rssn::kernels::holonomic::Expr;
use rssn::kernels::holonomic::Limits;
use rssn::kernels::holonomic::Poly;
use rssn::kernels::holonomic::RatFn;
use rssn::kernels::holonomic::annihilator;
use rssn::kernels::holonomic::deciding_terms;
use rssn::kernels::holonomic::is_zero;
use rssn::kernels::holonomic::recurrence;
use rssn::kernels::homotopy_cert::certify_homotopy_node;
use rssn::kernels::homotopy_cert::certify_zero_terms;
use rssn::kernels::homotopy_cert::homotopy_node_terms;
use rssn::kernels::homotopy_cert::verify_zero_terms;
use rssn::kernels::rootfind::PolyTerm;

type Q = BigRational;

fn qi(n: i64) -> Q {
    Q::from(BigInt::from(n))
}

fn q(
    n: i64,
    d: i64,
) -> Q {
    Q::new(BigInt::from(n), BigInt::from(d))
}

// ------------------------------------------------------------------ ball

fn grid(b: &Cball) -> Vec<Cx> {
    // corners, edge midpoints and the centre of the rectangle.
    let re = [b.re.lo, 0.5 * (b.re.lo + b.re.hi), b.re.hi];
    let im = [b.im.lo, 0.5 * (b.im.lo + b.im.hi), b.im.hi];
    let mut out = Vec::new();
    for &r in &re {
        for &i in &im {
            out.push(Cx::new(r, i));
        }
    }
    out
}

#[test]
fn ball_modulus_bounds_bracket_every_sample() {
    let boxes = [
        Cball::rectangle((1.0, 2.0), (3.0, 4.0)),
        Cball::rectangle((-2.0, -1.0), (0.5, 1.5)),
        Cball::rectangle((-1.0, 1.0), (-1.0, 1.0)),
        Cball::rectangle((-3.0, 2.0), (4.0, 5.0)),
        Cball::point(3.0, -4.0),
    ];
    for b in &boxes {
        let (lo, hi) = (b.modulus_lower(), b.modulus_upper());
        assert!(lo <= hi);
        for z in grid(b) {
            assert!(lo <= z.norm() + 1e-15, "{lo} > |{z}|");
            assert!(z.norm() <= hi + 1e-15, "|{z}| > {hi}");
        }
    }
    // |3 - 4i| = 5 exactly, a point box has no slack.
    let p = Cball::point(3.0, -4.0);
    assert!((p.modulus_lower() - 5.0).abs() < 1e-15 && (p.modulus_upper() - 5.0).abs() < 1e-15);
    // a box containing the origin has lower bound 0.
    assert_eq!(boxes[2].modulus_lower(), 0.0);
}

#[test]
fn ball_square_is_outward_and_geometry_is_consistent() {
    let s = Cball::square(0.25, -0.75, 0.125);
    // exact edges are strictly inside because of the one-ulp widening
    assert!(s.contains(Cx::new(0.25 - 0.125, -0.75)));
    assert!(s.contains(Cx::new(0.25 + 0.125, -0.75 + 0.125)));
    assert!(s.half_width() >= 0.125);
    assert!(s.half_width() <= 0.125 * (1.0 + 1e-14));
    assert!((s.centre() - Cx::new(0.25, -0.75)).norm() < 1e-15);
    // zero-radius square is a (widened) neighbourhood of the point
    let t = Cball::square(1.0, 1.0, 0.0);
    assert!(t.contains(Cx::new(1.0, 1.0)) && t.half_width() < 1e-15);
    // a rectangle's half width is the larger side
    let r = Cball::rectangle((0.0, 1.0), (0.0, 4.0));
    assert!((r.half_width() - 2.0).abs() < 1e-15);
    assert!((r.centre() - Cx::new(0.5, 2.0)).norm() < 1e-15);
    // non-finite boxes are flagged
    assert!(s.is_finite());
    assert!(!Cball::rectangle((0.0, f64::INFINITY), (0.0, 1.0)).is_finite());
    assert!(!Cball::point(f64::NAN, 0.0).is_finite());
    assert!(!Cball::square(0.0, 0.0, f64::INFINITY).is_finite());
}

#[test]
fn ball_intersection_is_symmetric_and_inside_strictly_is_transitive() {
    let a = Cball::rectangle((0.0, 4.0), (0.0, 4.0));
    let b = Cball::rectangle((2.0, 6.0), (-1.0, 3.0));
    let ab = a.intersect(&b).unwrap();
    let ba = b.intersect(&a).unwrap();
    assert_eq!(ab, ba);
    assert_eq!(
        (ab.re.lo, ab.re.hi, ab.im.lo, ab.im.hi),
        (2.0, 4.0, 0.0, 3.0)
    );
    // touching boxes meet in a degenerate rectangle, separated boxes do not
    let touch = Cball::rectangle((4.0, 5.0), (0.0, 1.0));
    let meet = a.intersect(&touch).unwrap();
    assert_eq!((meet.re.lo, meet.re.hi), (4.0, 4.0));
    assert!(
        a.intersect(&Cball::rectangle((4.5, 5.0), (0.0, 1.0)))
            .is_none()
    );
    // separation in the imaginary direction alone suffices
    assert!(
        a.intersect(&Cball::rectangle((1.0, 2.0), (4.5, 6.0)))
            .is_none()
    );
    // strict containment chain
    let small = Cball::rectangle((1.0, 2.0), (1.0, 2.0));
    let mid = Cball::rectangle((0.5, 3.0), (0.5, 3.0));
    let big = Cball::rectangle((0.0, 4.0), (0.0, 4.0));
    assert!(small.inside_strictly(&mid) && mid.inside_strictly(&big));
    assert!(small.inside_strictly(&big));
    assert!(!big.inside_strictly(&small));
    // sharing one face is not strict
    let flush = Cball::rectangle((0.0, 2.0), (1.0, 2.0));
    assert!(!flush.inside_strictly(&big));
}

#[test]
fn ball_arithmetic_is_inclusion_monotone_and_division_inverts_multiplication() {
    let outer_x = Cball::rectangle((0.5, 2.0), (-1.0, 1.0));
    let inner_x = Cball::rectangle((0.75, 1.5), (-0.5, 0.25));
    let outer_y = Cball::rectangle((-2.0, 1.0), (0.5, 2.0));
    let inner_y = Cball::rectangle((-1.5, 0.0), (0.75, 1.25));
    let po = outer_x * outer_y;
    let pi = inner_x * inner_y;
    for z in grid(&pi) {
        assert!(po.contains(z), "monotonicity lost {z}");
    }
    let so = outer_x - outer_y;
    for z in grid(&(inner_x - inner_y)) {
        assert!(so.contains(z));
    }
    // (a / b) * b contains a for a denominator bounded away from zero
    let a = Cball::rectangle((1.0, 1.5), (0.5, 1.0));
    let b = Cball::rectangle((2.0, 3.0), (1.0, 1.25));
    let back = a.div(&b).unwrap() * b;
    for z in grid(&a) {
        assert!(back.contains(z), "{z} lost in (a/b)*b");
    }
    // negation reflects through the origin and is an involution
    let n = -a;
    assert_eq!((n.re.lo, n.re.hi), (-1.5, -1.0));
    assert_eq!(-n, a);
    // x + (-x) contains zero
    assert!((a + n).contains(Cx::new(0.0, 0.0)));
    // powi: z^0 = 1, z^1 contains z, z^4 contains (z^2)^2 samples
    assert!(a.powi(0).contains(Cx::new(1.0, 0.0)));
    for z in grid(&a) {
        assert!(a.powi(1).contains(z));
        assert!(a.powi(4).contains(z.powi(4)));
    }
    // i^2 = -1
    let i2 = Cball::point(0.0, 1.0).powi(2);
    assert!(i2.contains(Cx::new(-1.0, 0.0)) && i2.half_width() < 1e-15);
    // a denominator box that touches zero gives no enclosure
    assert!(a.div(&Cball::zero()).is_none());
    assert!(a.div(&Cball::rectangle((-1.0, 1.0), (0.0, 0.0))).is_none());
}

fn t(
    coef: f64,
    exps: &[u32],
) -> PolyTerm {
    PolyTerm {
        coef,
        exps: exps.to_vec(),
    }
}

fn cx3_system() -> Vec<Vec<PolyTerm>> {
    // f0 = x^2 y - z^3 + 2, f1 = x y z - 1 + 0.5 x, f2 = y^3 + x z^2
    vec![
        vec![t(1.0, &[2, 1, 0]), t(-1.0, &[0, 0, 3]), t(2.0, &[0, 0, 0])],
        vec![t(1.0, &[1, 1, 1]), t(-1.0, &[0, 0, 0]), t(0.5, &[1, 0, 0])],
        vec![t(1.0, &[0, 3, 0]), t(1.0, &[1, 0, 2])],
    ]
}

#[test]
fn ball_jacobian_matches_central_differences() {
    // Second implementation: a complex-step-free central difference of
    // `values_at`, h = 1e-6, versus the analytic `jacobian_at`.
    let system = cx3_system();
    let p = [Cx::new(0.7, 0.2), Cx::new(-0.4, 0.9), Cx::new(1.1, -0.3)];
    let j = jacobian_at(&system, &p);
    assert_eq!((j.len(), j[0].len()), (3, 3));
    let h = 1e-6;
    for col in 0..3 {
        let (mut plus, mut minus) = (p, p);
        plus[col] += h;
        minus[col] -= h;
        let (fp, fm) = (values_at(&system, &plus), values_at(&system, &minus));
        for row in 0..3 {
            let fd = (fp[row] - fm[row]) / (2.0 * h);
            assert!(
                (fd - j[row][col]).norm() < 1e-7,
                "J[{row}][{col}]: {fd} vs {}",
                j[row][col]
            );
        }
    }
}

#[test]
fn ball_point_boxes_enclose_point_evaluation() {
    let system = cx3_system();
    let p = [Cx::new(0.7, 0.2), Cx::new(-0.4, 0.9), Cx::new(1.1, -0.3)];
    let boxes: Vec<Cball> = p.iter().map(|&z| Cball::of(z)).collect();
    let (vb, vp) = (eval_polys(&system, &boxes), values_at(&system, &p));
    let (jb, jp) = (eval_jacobian(&system, &boxes), jacobian_at(&system, &p));
    for i in 0..3 {
        assert!(vb[i].contains(vp[i]), "value {i}");
        assert!(vb[i].half_width() < 1e-13, "point box should stay tight");
        for k in 0..3 {
            assert!(jb[i][k].contains(jp[i][k]), "jac {i},{k}");
            assert!(jb[i][k].half_width() < 1e-13);
        }
    }
    // widening the inputs widens (never shrinks) the enclosures
    let wide: Vec<Cball> = p.iter().map(|z| Cball::square(z.re, z.im, 1e-3)).collect();
    let vw = eval_polys(&system, &wide);
    for i in 0..3 {
        assert!(vw[i].half_width() >= vb[i].half_width());
        assert!(vw[i].contains(vp[i]));
    }
}

#[test]
fn ball_terms_edge_cases() {
    // empty system, empty equation, zero-variable evaluation
    assert!(eval_terms(&[], &[]).is_empty());
    assert!(values_at_terms(&[], &[]).is_empty());
    assert!(jacobian_at_terms(&[], &[Cx::new(1.0, 0.0)]).is_empty());
    let empty_eq: Vec<Vec<CTerm>> = vec![vec![]];
    let v = eval_terms(&empty_eq, &[Cball::one()]);
    assert_eq!(v[0], Cball::zero());
    assert_eq!(
        values_at_terms(&empty_eq, &[Cx::new(2.0, 0.0)])[0],
        Cx::new(0.0, 0.0)
    );
    let j = eval_jacobian_terms(&empty_eq, &[Cball::one(), Cball::one()]);
    assert_eq!(j[0].len(), 2);
    assert!(j[0].iter().all(|c| *c == Cball::zero()));
    // a constant has a zero derivative
    let constant = vec![vec![CTerm::real(7.0, vec![0])]];
    let jc = jacobian_at_terms(&constant, &[Cx::new(3.0, 1.0)]);
    assert_eq!(jc[0][0], Cx::new(0.0, 0.0));
    // real_terms preserves shape and coefficients
    let sys = cx3_system();
    let lifted = real_terms(&sys);
    assert_eq!(lifted.len(), 3);
    for (a, b) in sys.iter().zip(&lifted) {
        assert_eq!(a.len(), b.len());
        for (p, c) in a.iter().zip(b) {
            assert_eq!(c.coef, Cx::new(p.coef, 0.0));
            assert_eq!(c.exps, p.exps);
        }
    }
    // an infinite coefficient propagates to a non-finite enclosure rather than a panic
    let bad = vec![vec![CTerm::real(f64::INFINITY, vec![1])]];
    assert!(!eval_terms(&bad, &[Cball::one()])[0].is_finite());
}

// --------------------------------------------------- homotopy_cert (terms)

#[test]
fn rigorous_zero_of_a_complex_coefficient_quadratic() {
    // z^2 - (3 + 4i) = 0 has the roots +-(2 + i).
    let system = vec![vec![
        CTerm::real(1.0, vec![2]),
        CTerm {
            coef: Cx::new(-3.0, -4.0),
            exps: vec![0],
        },
    ]];
    let ball = certify_zero_terms(&system, &[(2.000_000_1, 0.999_999_9)]).unwrap();
    assert!(
        ball.enclosure[0].contains(Cx::new(2.0, 1.0)),
        "root must be enclosed"
    );
    assert!(
        !ball.enclosure[0].contains(Cx::new(-2.0, -1.0)),
        "other root excluded"
    );
    assert!(ball.radius < 1e-6);
    assert!(ball.contains(2.0, 1.0));
    assert!(ball.distance(2.0, 1.0) < 1e-6);
    // the negative root is certified independently and the boxes are disjoint
    let other = certify_zero_terms(&system, &[(-2.0, -1.000_000_1)]).unwrap();
    assert!(other.enclosure[0].contains(Cx::new(-2.0, -1.0)));
    assert!(rssn::kernels::homotopy_cert::disjoint(&ball, &other));
    assert!(!rssn::kernels::homotopy_cert::disjoint(&ball, &ball));
}

#[test]
fn rigorous_zero_of_a_two_variable_complex_system() {
    // f0 = z0^2 - z1, f1 = z1 - (1 + i)  => z1 = 1 + i, z0 = +-sqrt(1 + i).
    let root = Cx::new(1.0, 1.0).sqrt();
    let system = vec![
        vec![CTerm::real(1.0, vec![2, 0]), CTerm::real(-1.0, vec![0, 1])],
        vec![
            CTerm::real(1.0, vec![0, 1]),
            CTerm {
                coef: Cx::new(-1.0, -1.0),
                exps: vec![0, 0],
            },
        ],
    ];
    let ball = verify_zero_terms(
        &system,
        &[(root.re + 1e-9, root.im - 1e-9), (1.0, 1.0)],
        1e-6,
    )
    .unwrap();
    assert!(ball.enclosure[0].contains(root));
    assert!(ball.enclosure[1].contains(Cx::new(1.0, 1.0)));
    assert!(!ball.enclosure[0].contains(-root));
    assert_eq!(ball.centre.len(), 2);
    // the reported centre is a near-exact zero
    let c: Vec<Cx> = ball.centre.iter().map(|&(a, b)| Cx::new(a, b)).collect();
    let res = values_at_terms(&system, &c);
    assert!(res.iter().all(|r| r.norm() < 1e-14), "{res:?}");
}

#[test]
fn rigorous_verification_refuses_what_it_cannot_prove() {
    // double zero z^2 = 0: the Jacobian vanishes at the root, no proof possible.
    let double = vec![vec![CTerm::real(1.0, vec![2])]];
    assert!(certify_zero_terms(&double, &[(1e-9, 0.0)]).is_none());
    assert!(verify_zero_terms(&double, &[(0.0, 0.0)], 1e-3).is_none());
    // a guess far from any zero of z^2 + 1 inside a tiny box: no zero inside, no proof.
    let circle = vec![vec![CTerm::real(1.0, vec![2]), CTerm::real(1.0, vec![0])]];
    assert!(verify_zero_terms(&circle, &[(3.0, 0.0)], 1e-6).is_none());
    // ... but near i the proof succeeds
    let ok = verify_zero_terms(&circle, &[(1e-8, 1.0)], 1e-6).unwrap();
    assert!(ok.enclosure[0].contains(Cx::new(0.0, 1.0)));
    // hostile arguments return None, never panic
    assert!(verify_zero_terms(&[], &[], 1e-6).is_none());
    assert!(verify_zero_terms(&circle, &[], 1e-6).is_none());
    assert!(verify_zero_terms(&circle, &[(0.0, 1.0), (0.0, 1.0)], 1e-6).is_none());
    for bad in [0.0, -1e-6, f64::NAN, f64::INFINITY] {
        assert!(
            verify_zero_terms(&circle, &[(0.0, 1.0)], bad).is_none(),
            "relative = {bad}"
        );
    }
    assert!(verify_zero_terms(&circle, &[(f64::NAN, 1.0)], 1e-6).is_none());
    assert!(certify_zero_terms(&circle, &[(f64::INFINITY, 0.0)]).is_none());
    assert!(certify_zero_terms(&[], &[]).is_none());
}

#[test]
fn homotopy_node_system_matches_the_closed_form() {
    // n = 1, f = x^2 - 2 homogenised: X1^2 - 2 X0^2, degree 2.
    let hom = vec![vec![t(1.0, &[0, 2]), t(-2.0, &[2, 0])]];
    let gamma = Cx::new(0.3, 0.8);
    let patch = [Cx::new(1.0, 0.0), Cx::new(0.5, -0.25)];
    let x = [Cx::new(1.0, 0.1), Cx::new(1.4, -0.2)];
    for tt in [Cx::new(0.0, 0.0), Cx::new(1.0, 0.0), Cx::new(0.37, 0.21)] {
        let sys = homotopy_node_terms(&hom, &[2], gamma, &patch, tt).unwrap();
        assert_eq!(sys.len(), 2, "n equations plus the patch row");
        let v = values_at_terms(&sys, &x);
        // H0 = gamma (1 - t)(X1^2 - X0^2) + t (X1^2 - 2 X0^2)
        let h0 = gamma * (Cx::new(1.0, 0.0) - tt) * (x[1] * x[1] - x[0] * x[0])
            + tt * (x[1] * x[1] - 2.0 * x[0] * x[0]);
        assert!((v[0] - h0).norm() < 1e-13, "t = {tt}: {} vs {h0}", v[0]);
        // patch row: -1 + a . X
        let p = -1.0 + patch[0] * x[0] + patch[1] * x[1];
        assert!((v[1] - p).norm() < 1e-13);
    }
    // t = 1 kills the start system exactly: H0 reduces to f
    let at_one = homotopy_node_terms(&hom, &[2], gamma, &patch, Cx::new(1.0, 0.0)).unwrap();
    let v = values_at_terms(&at_one, &[Cx::new(1.0, 0.0), Cx::new(2.0_f64.sqrt(), 0.0)]);
    assert!(v[0].norm() < 1e-14);
}

#[test]
fn homotopy_node_malformed_inputs_are_none() {
    let hom = vec![vec![t(1.0, &[0, 2]), t(-2.0, &[2, 0])]];
    let g = Cx::new(0.3, 0.8);
    let patch = [Cx::new(1.0, 0.0), Cx::new(0.0, 0.0)];
    let one = Cx::new(1.0, 0.0);
    assert!(homotopy_node_terms(&[], &[], g, &[Cx::new(1.0, 0.0)], one).is_none());
    assert!(
        homotopy_node_terms(&hom, &[], g, &patch, one).is_none(),
        "degree count"
    );
    assert!(homotopy_node_terms(&hom, &[2, 2], g, &patch, one).is_none());
    assert!(
        homotopy_node_terms(&hom, &[2], g, &patch[..1], one).is_none(),
        "patch length"
    );
    assert!(
        homotopy_node_terms(&[vec![]], &[2], g, &patch, one).is_none(),
        "empty equation"
    );
    let wrong_arity = vec![vec![t(1.0, &[2])]];
    assert!(
        homotopy_node_terms(&wrong_arity, &[2], g, &patch, one).is_none(),
        "exponent arity"
    );
    // certify_homotopy_node propagates the failure
    assert!(
        certify_homotopy_node(
            &wrong_arity,
            &[2],
            g,
            &patch,
            &[(1.0, 0.0), (1.0, 0.0)],
            (1.0, 0.0)
        )
        .is_none()
    );
}

#[test]
fn homotopy_node_is_certified_at_the_end_of_the_path() {
    // At t = 1 the node equations are {f(X) = 0, a.X = 1}: a unique simple zero
    // (1, sqrt 2) for a = (1, 0). The Krawczyk operator must prove it.
    let hom = vec![vec![t(1.0, &[0, 2]), t(-2.0, &[2, 0])]];
    let gamma = Cx::new(0.3, 0.8);
    let patch = [Cx::new(1.0, 0.0), Cx::new(0.0, 0.0)];
    let s2 = 2.0_f64.sqrt();
    let b = certify_homotopy_node(
        &hom,
        &[2],
        gamma,
        &patch,
        &[(1.0, 0.0), (s2 + 1e-10, 0.0)],
        (1.0, 0.0),
    )
    .expect("simple zero at t = 1");
    assert!(b.enclosure[0].contains(Cx::new(1.0, 0.0)));
    assert!(b.enclosure[1].contains(Cx::new(s2, 0.0)));
    assert!(!b.enclosure[1].contains(Cx::new(-s2, 0.0)));
    // the same node at a wrong point (outside any contractible box) is not certified
    assert!(
        certify_homotopy_node(
            &hom,
            &[2],
            gamma,
            &patch,
            &[(1.0, 0.0), (5.0, 5.0)],
            (1.0, 0.0)
        )
        .is_none()
    );
}

// ------------------------------------------------------------- holonomic

fn x_poly() -> Poly {
    vec![qi(0), qi(1)]
}

fn rat(p: Poly) -> RatFn {
    RatFn::poly(p)
}

fn atom(a: Atom) -> Expr {
    Expr::Atom(a)
}

fn add(
    a: Expr,
    b: Expr,
) -> Expr {
    Expr::Add(Box::new(a), Box::new(b))
}

fn mul(
    a: Expr,
    b: Expr,
) -> Expr {
    Expr::Mul(Box::new(a), Box::new(b))
}

fn scale(
    c: i64,
    a: Expr,
) -> Expr {
    Expr::Scale(qi(c), Box::new(a))
}

fn sq(a: Expr) -> Expr {
    Expr::Pow(Box::new(a), 2)
}

fn diff(a: Expr) -> Expr {
    Expr::Diff(Box::new(a))
}

fn sub(
    a: Expr,
    b: Expr,
) -> Expr {
    add(a, scale(-1, b))
}

fn one() -> Expr {
    atom(Atom::Rational(RatFn::constant(qi(1))))
}

fn sin() -> Expr {
    atom(Atom::Sin(x_poly()))
}

fn cos() -> Expr {
    atom(Atom::Cos(x_poly()))
}

fn decide(e: &Expr) -> Option<bool> {
    is_zero(e, &Limits::default())
}

/// Applies `sum p_i(x) D^i` to a truncated series; all coefficients must vanish.
fn operator_kills(
    op: &[Poly],
    e: &Expr,
) -> bool {
    let n = 20;
    let mut d = e.series(n + op.len() + 1).unwrap();
    let mut total = vec![Q::from(BigInt::from(0)); n];
    for p in op {
        for (k, total_k) in total.iter_mut().enumerate() {
            for (j, c) in p.iter().enumerate() {
                if j <= k {
                    *total_k += c * &d[k - j];
                }
            }
        }
        // differentiate the series
        let next: Vec<Q> = (1..d.len()).map(|m| &d[m] * qi(m as i64)).collect();
        d = next;
    }
    total.iter().all(|c| c == &Q::from(BigInt::from(0)))
}

#[test]
fn holonomic_decides_classical_identities() {
    // sin^2 + cos^2 = 1
    assert_eq!(decide(&sub(add(sq(sin()), sq(cos())), one())), Some(true));
    // ... and is not identically 0 without the "- 1"
    assert_eq!(decide(&add(sq(sin()), sq(cos()))), Some(false));
    // cosh^2 - sinh^2 = 1
    let ch = atom(Atom::Cosh(x_poly()));
    let sh = atom(Atom::Sinh(x_poly()));
    assert_eq!(decide(&sub(sub(sq(ch), sq(sh)), one())), Some(true));
    // sin(2x) = 2 sin x cos x
    let sin2x = atom(Atom::Sin(vec![qi(0), qi(2)]));
    assert_eq!(decide(&sub(sin2x, scale(2, mul(sin(), cos())))), Some(true));
    // cos(2x) = cos^2 - sin^2 ; and the wrong sign is detected
    let cos2x = atom(Atom::Cos(vec![qi(0), qi(2)]));
    assert_eq!(
        decide(&sub(cos2x.clone(), sub(sq(cos()), sq(sin())))),
        Some(true)
    );
    assert_eq!(decide(&sub(cos2x, add(sq(cos()), sq(sin())))), Some(false));
    // exp(x) exp(-x) = 1
    let e = atom(Atom::Exp(rat(x_poly())));
    let em = atom(Atom::Exp(rat(vec![qi(0), qi(-1)])));
    assert_eq!(decide(&sub(mul(e, em), one())), Some(true));
    // exp(x) - exp(2x) is a non-zero function
    let e2 = atom(Atom::Exp(rat(vec![qi(0), qi(2)])));
    assert_eq!(
        decide(&sub(atom(Atom::Exp(rat(x_poly()))), e2)),
        Some(false)
    );
}

#[test]
fn holonomic_decides_algebraic_and_inverse_function_identities() {
    // ((1+x)^(1/2))^2 = 1 + x
    let onepx = RatFn::poly(vec![qi(1), qi(1)]);
    let root = atom(Atom::Power(onepx.clone(), q(1, 2)));
    assert_eq!(
        decide(&sub(sq(root), atom(Atom::Rational(onepx.clone())))),
        Some(true)
    );
    // (1+x) * d/dx log(1+x) = 1
    let log = atom(Atom::Log(onepx.clone()));
    let lhs = mul(atom(Atom::Rational(onepx)), diff(log));
    assert_eq!(decide(&sub(lhs, one())), Some(true));
    // d/dx atan x = 1/(1+x^2)
    let atan = atom(Atom::Atan(rat(x_poly())));
    let denom = RatFn::new(vec![qi(1)], vec![qi(1), qi(0), qi(1)]).unwrap();
    assert_eq!(
        decide(&sub(diff(atan.clone()), atom(Atom::Rational(denom)))),
        Some(true)
    );
    // d/dx atanh x = 1/(1-x^2)
    let atanh = atom(Atom::Atanh(rat(x_poly())));
    let denom = RatFn::new(vec![qi(1)], vec![qi(1), qi(0), qi(-1)]).unwrap();
    assert_eq!(
        decide(&sub(diff(atanh), atom(Atom::Rational(denom)))),
        Some(true)
    );
    // d/dx asin x = (1-x^2)^(-1/2)
    let asin = atom(Atom::Asin(rat(x_poly())));
    let w = atom(Atom::Power(
        RatFn::poly(vec![qi(1), qi(0), qi(-1)]),
        q(-1, 2),
    ));
    assert_eq!(decide(&sub(diff(asin), w)), Some(true));
    // but d/dx atan x is not (1-x^2)^(-1)
    let wrong = RatFn::new(vec![qi(1)], vec![qi(1), qi(0), qi(-1)]).unwrap();
    assert_eq!(
        decide(&sub(diff(atan), atom(Atom::Rational(wrong)))),
        Some(false)
    );
}

#[test]
fn holonomic_diff_and_pow_agree_with_series_arithmetic() {
    // The series of D(sin) equals the series of cos, coefficient by coefficient.
    let a = diff(sin()).series(12).unwrap();
    let b = cos().series(12).unwrap();
    assert_eq!(a, b);
    // series of the closure expressions match direct rational arithmetic:
    // (sin + cos)^2 = 1 + 2 sin cos -> coefficients of sin(2x) at odd n.
    let lhs = Expr::Pow(Box::new(add(sin(), cos())), 2)
        .series(10)
        .unwrap();
    let sin2x = atom(Atom::Sin(vec![qi(0), qi(2)])).series(10).unwrap();
    for n in 0..10 {
        let expect = if n == 0 {
            qi(1) + sin2x[0].clone()
        } else {
            sin2x[n].clone()
        };
        assert_eq!(lhs[n], expect, "coefficient {n}");
    }
    // Pow(.., 0) is the constant 1
    let z = Expr::Pow(Box::new(sin()), 0).series(4).unwrap();
    assert_eq!(z, vec![qi(1), qi(0), qi(0), qi(0)]);
    // series(0) is empty
    assert!(sin().series(0).unwrap().is_empty());
}

#[test]
fn holonomic_annihilators_kill_their_functions_and_stay_small() {
    let cases: Vec<(Expr, usize)> = vec![
        (sin(), 2),
        (atom(Atom::Exp(rat(x_poly()))), 1),
        // sin * cos = sin(2x)/2 has a minimal order of 2, but the closure construction (documented as
        // "usually but not always minimal") returns an order-3 operator; any order <= 4 is acceptable.
        (mul(sin(), cos()), 4),
        (add(sin(), atom(Atom::Exp(rat(x_poly())))), 3),
        (mul(sin(), atom(Atom::Exp(rat(x_poly())))), 2),
        (diff(sin()), 2),
        (sq(sin()), 3),
        (atom(Atom::Log(RatFn::poly(vec![qi(1), qi(1)]))), 2),
    ];
    for (idx, (e, max_order)) in cases.iter().enumerate() {
        let op = annihilator(e, &Limits::default()).expect("operator");
        assert!(
            op.len() - 1 <= *max_order,
            "case {idx}: order {} > {max_order}",
            op.len() - 1
        );
        assert!(
            !op.last().unwrap().is_empty(),
            "leading coefficient must be non-zero"
        );
        assert!(operator_kills(&op, e), "operator does not annihilate");
    }
    // the zero function gets the operator [1]
    let zero = scale(0, sin());
    assert_eq!(
        annihilator(&zero, &Limits::default()),
        Some(vec![vec![qi(1)]])
    );
    // a rational function: the operator is of order 1
    let r = atom(Atom::Rational(
        RatFn::new(vec![qi(1)], vec![qi(1), qi(-1)]).unwrap(),
    ));
    let op = annihilator(&r, &Limits::default()).unwrap();
    assert_eq!(op.len(), 2);
    assert!(operator_kills(&op, &r));
}

#[test]
fn holonomic_atom_annihilators_and_validity() {
    let valid = [
        Atom::Sin(x_poly()),
        Atom::Cos(x_poly()),
        Atom::Exp(rat(x_poly())),
        Atom::Power(RatFn::poly(vec![qi(1), qi(2)]), q(-3, 2)),
        Atom::Log(RatFn::poly(vec![qi(1), qi(0), qi(1)])),
        Atom::Rational(RatFn::new(vec![qi(1)], vec![qi(2), qi(1)]).unwrap()),
    ];
    for a in &valid {
        assert!(a.is_valid());
        let op = a.annihilator().expect("operator");
        assert!(operator_kills(&op, &atom(a.clone())), "{a:?}");
    }
    let invalid = [
        Atom::Sin(vec![qi(1), qi(1)]),
        Atom::Exp(rat(vec![qi(1), qi(1)])),
        Atom::Power(RatFn::constant(qi(2)), q(1, 2)),
        Atom::Log(rat(x_poly())),
        Atom::Atan(rat(vec![qi(1)])),
        Atom::Rational(RatFn::new(vec![qi(1)], x_poly()).unwrap()), // pole at 0
    ];
    for a in &invalid {
        assert!(!a.is_valid(), "{a:?}");
        assert!(a.annihilator().is_none() && a.series(3).is_none());
        let e = atom(a.clone());
        assert_eq!(decide(&e), None, "invalid atoms are undecidable, not wrong");
        assert!(annihilator(&e, &Limits::default()).is_none());
        // and the invalidity infects any expression containing it
        assert_eq!(decide(&add(sin(), e)), None);
    }
}

#[test]
fn holonomic_recurrence_of_exp_and_sin_is_the_textbook_one() {
    // (D - 1) e^x = 0  =>  (m + 1) a_{m+1} - a_m = 0
    let rec = recurrence(&[vec![qi(-1)], vec![qi(1)]]);
    assert_eq!(rec, vec![(0, vec![qi(-1)]), (1, vec![qi(1), qi(1)])]);
    // D^2 + 1 : (m + 1)(m + 2) a_{m+2} + a_m = 0, verified on the sin / cos series
    let op: Vec<Poly> = vec![vec![qi(1)], vec![], vec![qi(1)]];
    let rec = recurrence(&op);
    assert_eq!(rec, vec![(0, vec![qi(1)]), (2, vec![qi(2), qi(3), qi(1)])]);
    for series in [sin().series(14).unwrap(), cos().series(14).unwrap()] {
        for m in 0..12usize {
            let mm = qi(m as i64);
            let lhs = (mm.clone() + qi(1)) * (mm + qi(2)) * &series[m + 2] + &series[m];
            assert_eq!(lhs, qi(0), "m = {m}");
        }
    }
    // the zero operator has no recurrence
    assert!(recurrence(&[]).is_empty());
    assert!(recurrence(&[vec![]]).is_empty());
}

#[test]
fn holonomic_deciding_terms_counts_undetermined_coefficients() {
    let lim = Limits::default();
    // empty operator: nothing to decide
    assert_eq!(deciding_terms(&[], &lim), Some(0));
    // exp: a_0 determines everything
    assert_eq!(deciding_terms(&[vec![qi(-1)], vec![qi(1)]], &lim), Some(1));
    // sin / cos: two free initial values
    assert_eq!(
        deciding_terms(&[vec![qi(1)], vec![], vec![qi(1)]], &lim),
        Some(2)
    );
    // x f' - 3 f = 0 has the solution x^3 (recurrence (m - 3) a_m = 0), so the
    // coefficient at the integer root m = 3 is free: four coefficients decide.
    assert_eq!(
        deciding_terms(&[vec![qi(-3)], vec![qi(0), qi(1)]], &lim),
        Some(4)
    );
    // too small a budget gives None, never a wrong answer
    let tight = Limits { max_terms: 2, ..lim };
    assert_eq!(
        deciding_terms(&[vec![qi(-3)], vec![qi(0), qi(1)]], &tight),
        None
    );
}

#[test]
fn holonomic_limits_never_produce_a_wrong_answer() {
    let identity = sub(add(sq(sin()), sq(cos())), one());
    let non_identity = add(sq(sin()), sq(cos()));
    for order in [0usize, 1, 2, 3, 6] {
        for degree in [0usize, 2, 8, 50] {
            let l = Limits {
                max_order: order,
                max_degree: degree,
                max_terms: 4000,
            };
            let a = is_zero(&identity, &l);
            let b = is_zero(&non_identity, &l);
            assert_ne!(a, Some(false), "identity refuted at {order}/{degree}");
            assert_ne!(b, Some(true), "non-identity proved at {order}/{degree}");
        }
    }
    // generous limits decide both
    assert_eq!(is_zero(&identity, &Limits::default()), Some(true));
    assert_eq!(is_zero(&non_identity, &Limits::default()), Some(false));
    // a deep power exceeds a small order budget rather than looping
    let big = Expr::Pow(Box::new(add(sin(), cos())), 9);
    let small = Limits {
        max_order: 3,
        ..Limits::default()
    };
    assert!(annihilator(&big, &small).is_none());
}

#[test]
fn holonomic_ratfn_is_a_field_of_fractions() {
    assert!(
        RatFn::new(vec![qi(1)], vec![]).is_none(),
        "zero denominator"
    );
    // (x^2 - 1)/(x - 1) = x + 1 in lowest terms, monic denominator
    let r = RatFn::new(vec![qi(-1), qi(0), qi(1)], vec![qi(-1), qi(1)]).unwrap();
    assert_eq!(r, RatFn::poly(vec![qi(1), qi(1)]));
    // scaling the denominator is normalised away
    let a = RatFn::new(vec![qi(2)], vec![qi(4), qi(2)]).unwrap();
    assert_eq!(a, RatFn::new(vec![qi(1)], vec![qi(2), qi(1)]).unwrap());
    // zero is canonical
    assert!(RatFn::new(vec![], vec![qi(5), qi(3)]).unwrap().is_zero());
    assert!(RatFn::zero().is_zero());
    assert!(!RatFn::constant(qi(1)).is_zero());
    // f - f = 0 and f / f = 1
    let f = RatFn::new(vec![qi(1), qi(1)], vec![qi(2), qi(0), qi(1)]).unwrap();
    assert!(f.sub(&f).is_zero());
    assert_eq!(f.div(&f).unwrap(), RatFn::constant(qi(1)));
    assert!(f.div(&RatFn::zero()).is_none());
    assert_eq!(f.add(&f.neg()), RatFn::zero());
    // product / sum distribute
    let g = RatFn::new(vec![qi(3)], vec![qi(1), qi(1)]).unwrap();
    let h = RatFn::poly(vec![qi(1), qi(2)]);
    assert_eq!(f.add(&g).mul(&h), f.mul(&h).add(&g.mul(&h)));
    // derivative of 1/x-like: d/dx 1/(1 - x) = 1/(1 - x)^2 ; series 1,1,1,... and 1,2,3,...
    let geo = RatFn::new(vec![qi(1)], vec![qi(1), qi(-1)]).unwrap();
    assert_eq!(geo.series(5).unwrap(), vec![qi(1); 5]);
    assert_eq!(
        geo.deriv().series(5).unwrap(),
        (1..=5).map(qi).collect::<Vec<_>>()
    );
    // product rule (fg)' = f'g + fg'
    assert_eq!(f.mul(&g).deriv(), f.deriv().mul(&g).add(&f.mul(&g.deriv())));
    // a pole at the origin has no Taylor series
    assert!(
        RatFn::new(vec![qi(1)], vec![qi(0), qi(1)])
            .unwrap()
            .series(3)
            .is_none()
    );
}
