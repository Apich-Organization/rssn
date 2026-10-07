//! Discretised function-space norms and orthogonalisation (ported from
//! `numerical_functional_analysis_test.rs`).

use rssn::kernels::functional_analysis::*;

#[test]

fn test_l1_norm() {
    let points = vec![(0.0, 1.0), (1.0, 1.0), (2.0, 1.0)];

    let res = l1_norm(&points);

    assert!((res - 2.0).abs() < 1e-9);
}

#[test]

fn test_l2_norm() {
    let points = vec![(0.0, 1.0), (1.0, 1.0), (2.0, 1.0)];

    let res = l2_norm(&points);

    assert!((res - 2.0f64.sqrt()).abs() < 1e-9);
}

#[test]

fn test_inner_product() {
    let f = vec![(0.0, 1.0), (1.0, 1.0)];

    let g = vec![(0.0, 2.0), (1.0, 2.0)];

    let res = inner_product(&f, &g).unwrap();

    assert!((res - 2.0).abs() < 1e-9);
}

#[test]

fn test_project() {
    let f = vec![(0.0, 1.0), (1.0, 1.0)];

    let g = vec![(0.0, 2.0), (1.0, 2.0)];

    let res = project(&f, &g).unwrap();

    // proj_g(f) = (<f, g> / <g, g>) * g
    // <f, g> = 1*2 * 1 = 2
    // <g, g> = 2*2 * 1 = 4
    // coeff = 2/4 = 0.5
    // proj = 0.5 * [(0, 2), (1, 2)] = [(0, 1), (1, 1)]
    assert_eq!(res.len(), 2);

    assert!((res[0].1 - 1.0).abs() < 1e-9);

    assert!((res[1].1 - 1.0).abs() < 1e-9);
}

#[test]

fn test_gram_schmidt() {
    let u1 = vec![(0.0, 1.0), (1.0, 1.0)];

    let u2 = vec![(0.0, 0.0), (1.0, 1.0)];

    let basis = vec![u1, u2];

    let orth = gram_schmidt(&basis).unwrap();

    assert_eq!(orth.len(), 2);

    // <v1, v2> should be 0
    let ip = inner_product(&orth[0], &orth[1]).unwrap();

    assert!(ip.abs() < 1e-15);
}

#[test]

fn test_gram_schmidt_orthonormal() {
    let u1 = vec![(0.0, 1.0), (1.0, 1.0)];

    let u2 = vec![(0.0, 0.0), (1.0, 1.0)];

    let basis = vec![u1, u2];

    let orthonorm = gram_schmidt_orthonormal(&basis).unwrap();

    assert_eq!(orthonorm.len(), 2);

    // Norms should be 1
    assert!((l2_norm(&orthonorm[0]) - 1.0).abs() < 1e-9);

    assert!((l2_norm(&orthonorm[1]) - 1.0).abs() < 1e-9);

    // Inner product should be 0
    let ip = inner_product(&orthonorm[0], &orthonorm[1]).unwrap();

    assert!(ip.abs() < 1e-15);
}

#[cfg(test)]
mod proptests {

    use proptest::prelude::*;
    use proptest::test_runner::RngSeed;

    use super::*;

    // Strategy to generate a vector of Y values.
    // We will map these to X values 0.0, 1.0, 2.0, ...
    fn fun_strategy() -> impl Strategy<Value = Vec<f64>> {
        proptest::collection::vec(-100.0..100.0f64, 2..20)
    }

    fn make_points(ys: &[f64]) -> Vec<(f64, f64)> {
        ys.iter().enumerate().map(|(i, &y)| (i as f64, y)).collect()
    }

    proptest! {
        #![proptest_config(ProptestConfig { rng_seed: RngSeed::Fixed(0x5EED), failure_persistence: None, ..ProptestConfig::default() })]

        #[test]
        fn prop_l2_norm_positive(ys in fun_strategy()) {
            let points = make_points(&ys);
            let norm = l2_norm(&points);
            prop_assert!(norm >= 0.0);
        }

        #[test]
        fn prop_inner_product_symmetry(ys1 in fun_strategy(), ys2 in fun_strategy()) {
            // Ensure same length
            let min_len = std::cmp::min(ys1.len(), ys2.len());
            let p1 = make_points(&ys1[..min_len]);
            let p2 = make_points(&ys2[..min_len]);

            let ip1 = inner_product(&p1, &p2).unwrap();
            let ip2 = inner_product(&p2, &p1).unwrap();

            // Allow for some floating point error
            prop_assert!((ip1 - ip2).abs() < 1e-9);
        }

        #[test]
        fn prop_inner_product_linearity(ys1 in fun_strategy(), ys2 in fun_strategy(), c in -10.0..10.0f64) {
             // Ensure same length
            let min_len = std::cmp::min(ys1.len(), ys2.len());
            let p1 = make_points(&ys1[..min_len]);
            let p2 = make_points(&ys2[..min_len]);

            // c * p1
            let cp1: Vec<(f64, f64)> = p1.iter().map(|&(x, y)| (x, c * y)).collect();

            let ip_scaled = inner_product(&cp1, &p2).unwrap();
            let ip_orig = inner_product(&p1, &p2).unwrap();

            // Check <c*f, g> = c * <f, g>
            // Note: inner_product can be large, so relative error might be better, or loose absolute.
            // Using a somewhat loose tolerance or checking relative difference if values are large.
            let diff = (ip_scaled - c * ip_orig).abs();
            let threshold = 1e-7 * ip_scaled.abs().max(1.0);
            prop_assert!(diff <= threshold, "diff: {}, threshold: {}", diff, threshold);
        }

        #[test]
        fn prop_projection_orthogonality(ys1 in fun_strategy(), ys2 in fun_strategy()) {
             // Ensure same length
            let min_len = std::cmp::min(ys1.len(), ys2.len());
            let f = make_points(&ys1[..min_len]);
            let g = make_points(&ys2[..min_len]);

            // Avoid division by zero if g is essentially zero
            let g_norm = l2_norm(&g);
            if g_norm > 1e-6 {
                let proj = project(&f, &g).unwrap();

                // residual = f - proj
                let residual: Vec<(f64, f64)> = f.iter().zip(proj.iter())
                    .map(|(&(x1, y1), &(_, y2))| (x1, y1 - y2))
                    .collect();

                let ip = inner_product(&residual, &g).unwrap();
                prop_assert!(ip.abs() < 1e-7 * l2_norm(&f).max(1.0) * g_norm.max(1.0));
            }
        }
    }
}

// ============================================================================
// Added: analytic reference values
// ============================================================================

fn sample(f: impl Fn(f64) -> f64, a: f64, b: f64, n: usize) -> Vec<(f64, f64)> {
    (0..=n).map(|i| {
        let x = a + (b - a) * i as f64 / n as f64;
        (x, f(x))
    }).collect()
}

#[test]
fn norms_of_sine_on_half_period() {
    let s = sample(f64::sin, 0.0, std::f64::consts::PI, 2000);
    assert!((l1_norm(&s) - 2.0).abs() < 1e-5);
    assert!((l2_norm(&s) - (std::f64::consts::PI / 2.0).sqrt()).abs() < 1e-5);
    assert!((infinity_norm(&s) - 1.0).abs() < 1e-6);
}

#[test]
fn l1_norm_of_absolute_value_uses_absolute_values() {
    // int_{-1}^{1} |x| dx = 1 (the grid contains the kink at 0).
    let s = sample(|x| x, -1.0, 1.0, 100);
    assert!((l1_norm(&s) - 1.0).abs() < 1e-12);
    assert_eq!(infinity_norm(&s), 1.0);
    assert_eq!(infinity_norm(&[]), 0.0);
}

#[test]
fn sines_of_different_frequency_are_orthogonal() {
    let n = 4000;
    let s1 = sample(|x| (x).sin(), 0.0, 2.0 * std::f64::consts::PI, n);
    let s2 = sample(|x| (2.0 * x).sin(), 0.0, 2.0 * std::f64::consts::PI, n);
    assert!(inner_product(&s1, &s2).unwrap_or(f64::NAN).abs() < 1e-8);
    // <sin, sin> = pi
    assert!((inner_product(&s1, &s1).unwrap_or(f64::NAN) - std::f64::consts::PI).abs() < 1e-5);
}

#[test]
fn inner_product_rejects_length_mismatch_and_grid_mismatch() {
    let a = [(0.0, 1.0), (1.0, 1.0)];
    let b = [(0.0, 1.0), (1.0, 1.0), (2.0, 1.0)];
    assert!(inner_product(&a, &b).is_err());
    // Different abscissae contribute nothing.
    let c = [(5.0, 1.0), (6.0, 1.0)];
    assert_eq!(inner_product(&a, &c), Ok(0.0));
}

#[test]
fn projection_edge_cases() {
    let f = [(0.0, 1.0), (1.0, 3.0)];
    let zero = [(0.0, 0.0), (1.0, 0.0)];
    assert_eq!(project(&f, &zero), Ok(vec![(0.0, 0.0), (1.0, 0.0)]));
    // Projecting onto itself returns itself.
    let p = project(&f, &f).unwrap_or_default();
    for ((_, a), (_, b)) in p.iter().zip(&f) {
        assert!((a - b).abs() < 1e-12);
    }
}

#[test]
fn normalize_gives_unit_norm_and_leaves_zero_alone() {
    let f = sample(|x| 3.0 * x, 0.0, 1.0, 100);
    assert!((l2_norm(&normalize(&f)) - 1.0).abs() < 1e-12);
    let zero = [(0.0, 0.0), (1.0, 0.0)];
    assert_eq!(normalize(&zero), zero.to_vec());
}

#[test]
fn gram_schmidt_three_polynomials_is_orthonormal() {
    let n = 400;
    let basis = vec![sample(|_| 1.0, -1.0, 1.0, n), sample(|x| x, -1.0, 1.0, n), sample(|x| x * x, -1.0, 1.0, n)];
    let q = gram_schmidt_orthonormal(&basis).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(q.len(), 3);
    for i in 0..3 {
        for j in 0..3 {
            let ip = inner_product(&q[i], &q[j]).unwrap_or(f64::NAN);
            let want = if i == j { 1.0 } else { 0.0 };
            assert!((ip - want).abs() < 1e-9, "<q{i}, q{j}> = {ip}");
        }
    }
    // Legendre P2 shape: the third vector is even and proportional to x^2 - 1/3.
    let mid = n / 2;
    assert!(q[2][mid].1 * q[2][0].1 < 0.0, "P2 changes sign between centre and ends");
}
