//! Tests for the root-finding kernels.

use rssn::kernels::dense::Mat;
use rssn::kernels::rootfind::{
    RootError, brent, broyden, halley, homotopy_univariate, levenberg_marquardt,
    newton_bracketed, newton_system, polynomial_roots_aberth, polynomial_roots_companion,
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
