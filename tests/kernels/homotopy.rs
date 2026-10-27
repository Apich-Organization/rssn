//! Tests for the projective homotopy solver with endgame.

use rssn::kernels::homotopy::{HomotopyOptions, PathEnd, Precision, homotopy_solve};
use rssn::kernels::rootfind::PolyTerm;

fn t(coef: f64, exps: &[u32]) -> PolyTerm {
    PolyTerm { coef, exps: exps.to_vec() }
}

fn finite_regular(v: &[PathEnd]) -> Vec<&PathEnd> {
    v.iter().filter(|p| !p.at_infinity && !p.singular).collect()
}

#[test]
fn solutions_at_infinity_in_projective_coordinates() {
    // x^2 - 1 = 0, x y - 1 = 0: solutions (1,1), (-1,-1) and a double root at infinity [0:1:0]
    let sys = vec![vec![t(1.0, &[2, 0]), t(-1.0, &[0, 0])], vec![t(1.0, &[1, 1]), t(-1.0, &[0, 0])]];
    let r = homotopy_solve(&sys, &HomotopyOptions::default()).unwrap();
    assert_eq!(r.len(), 4);
    let fin = finite_regular(&r);
    assert_eq!(fin.len(), 2);
    for p in &fin {
        let (x, y) = (p.x[0], p.x[1]);
        assert!((x.0.abs() - 1.0).abs() < 1e-10 && x.1.abs() < 1e-10);
        assert!((x.0 - y.0).abs() < 1e-10 && p.residual < 1e-12);
        assert!(!p.suspect);
    }
    let inf: Vec<&PathEnd> = r.iter().filter(|p| p.at_infinity).collect();
    assert_eq!(inf.len(), 2);
    for p in inf {
        assert!(p.singular, "winding {}", p.winding);
        // the point at infinity is [0 : 0 : 1] up to the patch scaling: x_1 = 0
        assert!(p.homogeneous[1].0.hypot(p.homogeneous[1].1) < 1e-6 * p.homogeneous[2].0.hypot(p.homogeneous[2].1));
    }
}

#[test]
fn singular_finite_endpoints_via_endgame() {
    // (x - 1)^3 = 0, y - x^2 = 0 : triple root (1, 1)
    let sys = vec![
        vec![t(1.0, &[3, 0]), t(-3.0, &[2, 0]), t(3.0, &[1, 0]), t(-1.0, &[0, 0])],
        vec![t(1.0, &[0, 1]), t(-1.0, &[2, 0])],
    ];
    let r = homotopy_solve(&sys, &HomotopyOptions::default()).unwrap();
    assert_eq!(r.len(), 6);
    let mut tri = 0;
    let mut at_inf = 0;
    for p in &r {
        if p.at_infinity {
            at_inf += 1;
            continue;
        }
        assert!(p.singular && (1..=2).contains(&p.winding), "winding {}", p.winding);
        assert!((p.x[0].0 - 1.0).abs() < 1e-7 && p.x[0].1.abs() < 1e-7, "{:?}", p.x[0]);
        assert!((p.x[1].0 - 1.0).abs() < 1e-6 && p.x[1].1.abs() < 1e-6);
        tri += 1;
    }
    assert_eq!(tri, 3);
    assert_eq!(at_inf, 3);
    // double root: (x - 2)^2 (x + 1) = 0 in one variable
    let sys = vec![vec![t(1.0, &[3]), t(-3.0, &[2]), t(0.0, &[1]), t(4.0, &[0])]];
    let r = homotopy_solve(&sys, &HomotopyOptions::default()).unwrap();
    let mut roots: Vec<(f64, usize)> = r.iter().map(|p| (p.x[0].0, p.winding)).collect();
    roots.sort_by(|a, b| a.0.total_cmp(&b.0));
    assert!((roots[0].0 + 1.0).abs() < 1e-10 && roots[0].1 == 1);
    for k in 1..3 {
        assert!((roots[k].0 - 2.0).abs() < 1e-7 && roots[k].1 == 2, "{:?}", roots[k]);
    }
}

#[test]
fn regular_systems_and_path_jumping() {
    // cyclic 3-roots: x+y+z = 0, xy+yz+zx = 0, xyz - 1 = 0 -> 6 regular solutions
    let sys = vec![
        vec![t(1.0, &[1, 0, 0]), t(1.0, &[0, 1, 0]), t(1.0, &[0, 0, 1])],
        vec![t(1.0, &[1, 1, 0]), t(1.0, &[0, 1, 1]), t(1.0, &[1, 0, 1])],
        vec![t(1.0, &[1, 1, 1]), t(-1.0, &[0, 0, 0])],
    ];
    let r = homotopy_solve(&sys, &HomotopyOptions::default()).unwrap();
    assert_eq!(r.len(), 6);
    assert!(r.iter().all(|p| !p.at_infinity && !p.singular && !p.suspect && p.residual < 1e-12));
    for i in 0..6 {
        for j in i + 1..6 {
            let d = r[i].x.iter().zip(&r[j].x).map(|(a, b)| (a.0 - b.0).hypot(a.1 - b.1)).fold(0.0, f64::max);
            assert!(d > 1e-3, "paths {i} and {j} ended at the same point");
        }
    }
    // moderately clustered roots (spacing 0.25) of a degree-8 polynomial: all paths must stay distinct
    let roots: Vec<f64> = (0..8).map(|k| 1.0 + 0.25 * f64::from(k)).collect();
    let mut poly = vec![1.0];
    for &r0 in &roots {
        let mut next = vec![0.0; poly.len() + 1];
        for (i, c) in poly.iter().enumerate() {
            next[i] += c;
            next[i + 1] -= c * r0;
        }
        poly = next;
    }
    let sys = vec![poly.iter().enumerate().map(|(i, &c)| t(c, &[(poly.len() - 1 - i) as u32])).collect::<Vec<_>>()];
    let r = homotopy_solve(&sys, &HomotopyOptions::default()).unwrap();
    let mut got: Vec<f64> = r.iter().map(|p| p.x[0].0).collect();
    got.sort_by(f64::total_cmp);
    for (g, w) in got.iter().zip(&roots) {
        assert!((g - w).abs() < 1e-5, "{g} vs {w}");
    }
    assert!(r.iter().all(|p| !p.suspect));
    assert!(homotopy_solve(&[], &HomotopyOptions::default()).is_err());
}

/// `x (x - d)(x - 3d)` with `d = 2^-27`: three roots 7.5e-9 apart (centroid `4d/3`, not a root), coefficients exact in `f64`.
fn clustered_cubic() -> (Vec<Vec<PolyTerm>>, f64) {
    let d = 2.0_f64.powi(-27);
    let sys = vec![vec![t(1.0, &[3]), t(-4.0 * d, &[2]), t(3.0 * d * d, &[1]), t(0.0, &[0])]];
    (sys, d)
}

#[test]
fn double_double_resolves_tight_clusters() {
    let (sys, d) = clustered_cubic();
    for precision in [Precision::DoubleDouble, Precision::Adaptive] {
        let opts = HomotopyOptions { precision, ..HomotopyOptions::default() };
        let r = homotopy_solve(&sys, &opts).unwrap();
        assert_eq!(r.len(), 3);
        let mut roots: Vec<(f64, f64)> = r.iter().map(|p| (p.x[0].0 + p.x_lo[0].0, p.x[0].1)).collect();
        roots.sort_by(|a, b| a.0.total_cmp(&b.0));
        for (k, rt) in roots.iter().enumerate() {
            let want = [0.0, 1.0, 3.0][k] * d;
            assert!((rt.0 - want).abs() < 1e-5 * d, "{precision:?}: root {k}: {} vs {want}", rt.0);
            assert!(rt.1.abs() < 1e-5 * d);
        }
        assert!(
            r.iter().all(|p| p.winding == 1 && !p.singular && p.converged && !p.suspect),
            "{precision:?}: {r:?}"
        );
        assert!(r.iter().all(|p| p.extended), "{precision:?}");
    }
    // ordinary systems report extra digits in the extended mode
    let sys = vec![vec![t(1.0, &[2]), t(-2.0, &[0])]];
    let opts = HomotopyOptions { precision: Precision::DoubleDouble, ..HomotopyOptions::default() };
    let r = homotopy_solve(&sys, &opts).unwrap();
    let mut v: Vec<(f64, f64)> = r.iter().map(|p| (p.x[0].0, p.x_lo[0].0)).collect();
    v.sort_by(|a, b| a.0.total_cmp(&b.0));
    // sqrt(2) = 1.41421356237309504880168872420969807856967...
    let hi = 2.0_f64.sqrt();
    let lo = -9.667_293_313_452_913e-17;
    assert!((v[1].0 - hi).abs() < 1e-15);
    assert!((v[1].1 - lo).abs() < 1e-30, "{} vs {lo}", v[1].1);
    assert!((v[0].0 + hi).abs() < 1e-15);
}

#[test]
fn adaptive_leaves_easy_paths_in_double() {
    let sys = vec![vec![t(1.0, &[2, 0]), t(-1.0, &[0, 0])], vec![t(1.0, &[1, 1]), t(-1.0, &[0, 0])]];
    let opts = HomotopyOptions { precision: Precision::Adaptive, ..HomotopyOptions::default() };
    let r = homotopy_solve(&sys, &opts).unwrap();
    assert_eq!(r.len(), 4);
    // the two regular solutions stay in f64
    assert!(r.iter().filter(|p| !p.at_infinity && !p.singular).all(|p| !p.extended));
}
