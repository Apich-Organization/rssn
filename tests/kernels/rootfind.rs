//! Tests for the root-finding kernels.

use rssn::kernels::dense::Mat;
use rssn::kernels::rootfind::{
    RootError, brent, broyden, halley, homotopy_univariate, levenberg_marquardt,
    newton_bracketed, newton_system, polynomial_roots_aberth, polynomial_roots_companion,
    polynomial_roots_jenkins_traub, polynomial_roots_jenkins_traub_complex, homotopy_system,
    PolyTerm,
};

#[test]
fn brent_roots() {
    let r = brent(|x| x * x - 2.0, 0.0, 2.0, 1e-14, 100).unwrap();
    assert!((r - 2.0_f64.sqrt()).abs() < 1e-13);
    let r = brent(|x| x.cos() - x, 0.0, 1.0, 1e-14, 100).unwrap();
    assert!((r - 0.739_085_133_215_160_6).abs() < 1e-13);
    assert_eq!(brent(|x| x * x + 1.0, -1.0, 1.0, 1e-12, 50), Err(RootError::NotBracketed));
}

#[test]
fn halley_cubic_convergence() {
    let r = halley(|x| x * x * x - 2.0, |x| 3.0 * x * x, |x| 6.0 * x, 1.0, 1e-15, 20).unwrap();
    assert!((r - 2.0_f64.cbrt()).abs() < 1e-14);
}

#[test]
fn bracketed_newton() {
    let r = newton_bracketed(|x| x * x * x - x - 1.0, |x| 3.0 * x * x - 1.0, 1.0, 2.0, 1e-14, 100)
        .unwrap();
    assert!((r - 1.324_717_957_244_746).abs() < 1e-12);
    let r = newton_bracketed(f64::atan, |x| 1.0 / (1.0 + x * x), -10.0, 15.0, 1e-13, 200).unwrap();
    assert!(r.abs() < 1e-10);
}

fn sorted_re(mut v: Vec<(f64, f64)>) -> Vec<(f64, f64)> {
    v.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.total_cmp(&b.1)));
    v
}

#[test]
fn polynomial_root_methods_agree() {
    let c = [-24.0, 38.0, -13.0, -2.0, 1.0];
    let expect = [-4.0, 1.0, 2.0, 3.0];
    for roots in [
        polynomial_roots_aberth(&c).unwrap(),
        polynomial_roots_companion(&c).unwrap(),
        homotopy_univariate(&c).unwrap(),
    ] {
        let r = sorted_re(roots);
        for (v, e) in r.iter().zip(expect) {
            assert!((v.0 - e).abs() < 1e-9 && v.1.abs() < 1e-9, "{r:?}");
        }
    }
}

#[test]
fn complex_polynomial_roots() {
    let c = [-1.0, 0.0, 0.0, 1.0];
    for roots in [polynomial_roots_aberth(&c).unwrap(), homotopy_univariate(&c).unwrap()] {
        for (re, im) in roots {
            assert!(((re * re + im * im) - 1.0).abs() < 1e-10);
            let cube_re = re * re * re - 3.0 * re * im * im;
            assert!((cube_re - 1.0).abs() < 1e-9);
        }
    }
    assert!(polynomial_roots_aberth(&[1.0]).is_err());
    let r = polynomial_roots_companion(&[1.0, -2.0, 1.0]).unwrap();
    assert!(r.iter().all(|z| (z.0 - 1.0).abs() < 1e-6));
}

fn system(x: &[f64]) -> Vec<f64> {
    vec![x[0] * x[0] + x[1] * x[1] - 4.0, x[0] * x[1] - 1.0]
}

#[test]
fn newton_broyden_system() {
    let s = newton_system(system, None::<fn(&[f64]) -> Mat>, &[2.0, 0.5], 1e-12, 50).unwrap();
    assert!(system(&s.x).iter().all(|v| v.abs() < 1e-11));
    let jac = |x: &[f64]| Mat { rows: 2, cols: 2, data: vec![2.0 * x[0], 2.0 * x[1], x[1], x[0]] };
    let s2 = newton_system(system, Some(jac), &[2.0, 0.5], 1e-12, 50).unwrap();
    assert!((s.x[0] - s2.x[0]).abs() < 1e-8);
    let b = broyden(system, &[2.0, 0.5], 1e-10, 100).unwrap();
    assert!(system(&b.x).iter().all(|v| v.abs() < 1e-9));
}

#[test]
fn newton_line_search_helps() {
    let f = |x: &[f64]| vec![x[0].atan()];
    let s = newton_system(f, None::<fn(&[f64]) -> Mat>, &[3.0], 1e-12, 100).unwrap();
    assert!(s.x[0].abs() < 1e-10);
}

#[test]
fn levenberg_marquardt_fit() {
    let xs: Vec<f64> = (0..10).map(|i| f64::from(i) * 0.5).collect();
    let ys: Vec<f64> = xs.iter().map(|x| 2.0 * (-0.5 * x).exp()).collect();
    let res = move |p: &[f64]| -> Vec<f64> {
        xs.iter().zip(&ys).map(|(x, y)| p[0] * (p[1] * x).exp() - y).collect()
    };
    let s = levenberg_marquardt(res, &[1.0, 0.0], 1e-12, 200).unwrap();
    assert!((s.x[0] - 2.0).abs() < 1e-6 && (s.x[1] + 0.5).abs() < 1e-6, "{s:?}");
}

fn poly_from_roots(roots: &[f64]) -> Vec<f64> {
    let mut c = vec![1.0];
    for &r in roots {
        let mut n = vec![0.0; c.len() + 1];
        for (i, &v) in c.iter().enumerate() {
            n[i + 1] += v;
            n[i] -= r * v;
        }
        c = n;
    }
    c
}

fn sorted_by_re(mut r: Vec<(f64, f64)>) -> Vec<(f64, f64)> {
    r.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap().then(a.1.partial_cmp(&b.1).unwrap()));
    r
}

/// Backward error of a root: |p(z)| / sum |a_i| |z|^i.
fn backward_err(c: &[f64], z: (f64, f64)) -> f64 {
    let (v, _) = rssn::kernels::rootfind::poly_eval_complex(c, z);
    let m = z.0.hypot(z.1);
    let s: f64 = c.iter().enumerate().map(|(i, a)| a.abs() * m.powi(i as i32)).sum();
    v.0.hypot(v.1) / s
}

#[test]
fn jenkins_traub_wilkinson_12() {
    let roots: Vec<f64> = (1..=12).map(f64::from).collect();
    let c = poly_from_roots(&roots);
    let jt = sorted_by_re(polynomial_roots_jenkins_traub(&c).unwrap());
    let ab = sorted_by_re(polynomial_roots_aberth(&c).unwrap());
    for (k, (j, a)) in jt.iter().zip(&ab).enumerate() {
        let want = (k + 1) as f64;
        assert!((j.0 - want).abs() < 1e-6 && j.1.abs() < 1e-6, "jt {j:?} vs {want}");
        assert!((j.0 - a.0).abs() < 1e-5, "jt {j:?} aberth {a:?}");
    }
}

#[test]
fn jenkins_traub_wilkinson_20() {
    let roots: Vec<f64> = (1..=20).map(f64::from).collect();
    let c = poly_from_roots(&roots);
    let jt = polynomial_roots_jenkins_traub(&c).unwrap();
    assert_eq!(jt.len(), 20);
    // backward stable: residuals tiny relative to the coefficient size
    for z in &jt {
        assert!(backward_err(&c, *z) < 1e-12, "jt {z:?} {}", backward_err(&c, *z));
    }
    // the Aberth iteration is much less reliable on this polynomial
    // (it can return non-finite or huge values); JT must beat or equal it
    let jt_worst = jt.iter().map(|z| backward_err(&c, *z)).fold(0.0, f64::max);
    if let Ok(ab) = polynomial_roots_aberth(&c) {
        let ab_worst = ab.iter().map(|z| backward_err(&c, *z)).fold(0.0, f64::max);
        assert!(!(ab_worst < jt_worst) || ab_worst < 1e-12);
    }
    // QR eigenvalues of the companion matrix agree up to the conditioning
    let cq = polynomial_roots_companion(&c).unwrap();
    for z in &jt {
        let d = cq.iter().map(|w| (z.0 - w.0).hypot(z.1 - w.1)).fold(f64::INFINITY, f64::min);
        assert!(d < 0.1, "jt root {z:?} has no companion partner (d = {d})");
    }
    // and approximate the true integer roots (perturbed pairs allowed)
    for z in &jt {
        let d = roots.iter().map(|r| (z.0 - r).hypot(z.1)).fold(f64::INFINITY, f64::min);
        assert!(d < 0.5, "{z:?}");
    }
}

#[test]
fn jenkins_traub_clustered_and_multiple_roots() {
    // five roots clustered within 4e-3, plus a well-separated set
    let mut roots: Vec<f64> = (0..5).map(|k| 1.0 + 1e-3 * f64::from(k)).collect();
    roots.extend([-2.0, 3.0, 5.5]);
    let c = poly_from_roots(&roots);
    let jt = sorted_by_re(polynomial_roots_jenkins_traub(&c).unwrap());
    for z in &jt {
        assert!(backward_err(&c, *z) < 1e-12, "{z:?}");
    }
    // the isolated roots are recovered accurately
    for want in [-2.0, 3.0, 5.5] {
        assert!(jt.iter().any(|z| (z.0 - want).abs() < 1e-8 && z.1.abs() < 1e-8), "{want}");
    }
    // the cluster is recovered to the conditioning of the problem
    let cluster: Vec<&(f64, f64)> = jt.iter().filter(|z| z.0 > 0.99 && z.0 < 1.01).collect();
    assert_eq!(cluster.len(), 5);
    let cx = cluster.iter().map(|z| z.0).sum::<f64>() / 5.0;
    assert!((cx - 1.002).abs() < 1e-8, "centroid {cx}");
    // multiple root (z-1)^4 (z+2): accuracy ~ eps^(1/4)
    let c = poly_from_roots(&[1.0, 1.0, 1.0, 1.0, -2.0]);
    let jt = polynomial_roots_jenkins_traub(&c).unwrap();
    let near1 = jt.iter().filter(|z| (z.0 - 1.0).hypot(z.1) < 1e-2).count();
    assert_eq!(near1, 4);
    assert!(jt.iter().any(|z| (z.0 + 2.0).abs() < 1e-8));
}

#[test]
fn jenkins_traub_complex_zero_roots_and_unity() {
    // (z - (1+2i)) (z - (3-i)) = z^2 - (4+i) z + (5+5i)
    let r = polynomial_roots_jenkins_traub_complex(&[(5.0, 5.0), (-4.0, -1.0), (1.0, 0.0)]).unwrap();
    for want in [(1.0, 2.0), (3.0, -1.0)] {
        assert!(r.iter().any(|z| (z.0 - want.0).hypot(z.1 - want.1) < 1e-12), "{want:?}");
    }
    // z^3 (z - 2)
    let r = polynomial_roots_jenkins_traub(&[0.0, 0.0, 0.0, -2.0, 1.0]).unwrap();
    assert_eq!(r.iter().filter(|z| z.0.hypot(z.1) < 1e-12).count(), 3);
    assert!(r.iter().any(|z| (z.0 - 2.0).abs() < 1e-12));
    // z^30 - 1
    let mut c = vec![0.0; 31];
    c[0] = -1.0;
    c[30] = 1.0;
    let r = polynomial_roots_jenkins_traub(&c).unwrap();
    assert_eq!(r.len(), 30);
    for z in &r {
        assert!((z.0.hypot(z.1) - 1.0).abs() < 1e-12);
    }
    for k in 0..30 {
        let a = 2.0 * std::f64::consts::PI * f64::from(k) / 30.0;
        assert!(r.iter().any(|z| (z.0 - a.cos()).hypot(z.1 - a.sin()) < 1e-10));
    }
    assert!(polynomial_roots_jenkins_traub(&[3.0]).is_err());
}

fn term(coef: f64, exps: &[u32]) -> PolyTerm {
    PolyTerm { coef, exps: exps.to_vec() }
}

fn finite(sols: &[rssn::kernels::rootfind::PathResult]) -> Vec<Vec<(f64, f64)>> {
    sols.iter().filter(|p| p.finite).map(|p| p.x.clone()).collect()
}

#[test]
fn homotopy_system_2x2() {
    // x^2 + y^2 = 5, x y = 2  ->  (1,2) (2,1) (-1,-2) (-2,-1)
    let sys = vec![
        vec![term(1.0, &[2, 0]), term(1.0, &[0, 2]), term(-5.0, &[0, 0])],
        vec![term(1.0, &[1, 1]), term(-2.0, &[0, 0])],
    ];
    let sols = homotopy_system(&sys).unwrap();
    assert_eq!(sols.len(), 4);
    let fin = finite(&sols);
    assert_eq!(fin.len(), 4);
    for want in [(1.0, 2.0), (2.0, 1.0), (-1.0, -2.0), (-2.0, -1.0)] {
        assert!(
            fin.iter().any(|s| (s[0].0 - want.0).hypot(s[0].1) < 1e-9 && (s[1].0 - want.1).hypot(s[1].1) < 1e-9),
            "{want:?}"
        );
    }
    // circle and line: x^2 + y^2 = 1, x = y
    let sys = vec![
        vec![term(1.0, &[2, 0]), term(1.0, &[0, 2]), term(-1.0, &[0, 0])],
        vec![term(1.0, &[1, 0]), term(-1.0, &[0, 1])],
    ];
    let fin = finite(&homotopy_system(&sys).unwrap());
    assert_eq!(fin.len(), 2);
    let h = 0.5_f64.sqrt();
    for s in &fin {
        assert!((s[0].0.abs() - h).abs() < 1e-10 && (s[0].0 - s[1].0).abs() < 1e-10);
    }
}

#[test]
fn homotopy_system_solutions_at_infinity() {
    // x^2 = 1, x y = 1: Bezout number 4 but only (1,1) and (-1,-1) are finite
    let sys = vec![
        vec![term(1.0, &[2, 0]), term(-1.0, &[0, 0])],
        vec![term(1.0, &[1, 1]), term(-1.0, &[0, 0])],
    ];
    let sols = homotopy_system(&sys).unwrap();
    assert_eq!(sols.len(), 4);
    let fin = finite(&sols);
    assert_eq!(fin.len(), 2, "{sols:?}");
    for s in &fin {
        assert!((s[0].0 * s[1].0 - 1.0).abs() < 1e-9 && s[0].0.abs() - 1.0 < 1e-9);
    }
}

#[test]
fn homotopy_system_3x3() {
    // elementary symmetric: roots are permutations of (1, 2, 3)
    let sys = vec![
        vec![term(1.0, &[1, 0, 0]), term(1.0, &[0, 1, 0]), term(1.0, &[0, 0, 1]), term(-6.0, &[0, 0, 0])],
        vec![term(1.0, &[1, 1, 0]), term(1.0, &[0, 1, 1]), term(1.0, &[1, 0, 1]), term(-11.0, &[0, 0, 0])],
        vec![term(1.0, &[1, 1, 1]), term(-6.0, &[0, 0, 0])],
    ];
    let sols = homotopy_system(&sys).unwrap();
    assert_eq!(sols.len(), 6);
    let fin = finite(&sols);
    assert_eq!(fin.len(), 6);
    let perms = [[1.0, 2.0, 3.0], [1.0, 3.0, 2.0], [2.0, 1.0, 3.0], [2.0, 3.0, 1.0], [3.0, 1.0, 2.0], [3.0, 2.0, 1.0]];
    for p in perms {
        assert!(
            fin.iter().any(|s| (0..3).all(|i| (s[i].0 - p[i]).hypot(s[i].1) < 1e-8)),
            "{p:?}"
        );
    }
    // complex solutions: e1 = e2 = 0, e3 = 1 -> permutations of the cube roots of unity
    let sys = vec![
        vec![term(1.0, &[1, 0, 0]), term(1.0, &[0, 1, 0]), term(1.0, &[0, 0, 1])],
        vec![term(1.0, &[1, 1, 0]), term(1.0, &[0, 1, 1]), term(1.0, &[1, 0, 1])],
        vec![term(1.0, &[1, 1, 1]), term(-1.0, &[0, 0, 0])],
    ];
    let fin = finite(&homotopy_system(&sys).unwrap());
    assert_eq!(fin.len(), 6);
    for s in &fin {
        for v in s {
            let m = v.0.hypot(v.1);
            assert!((m - 1.0).abs() < 1e-8);
            // cube of the root is 1
            let (re, im) = (v.0 * v.0 - v.1 * v.1, 2.0 * v.0 * v.1);
            let (c3r, c3i) = (re * v.0 - im * v.1, re * v.1 + im * v.0);
            assert!((c3r - 1.0).abs() < 1e-8 && c3i.abs() < 1e-8);
        }
    }
    assert!(homotopy_system(&[]).is_err());
}
