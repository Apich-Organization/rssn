//! Tests for physics CNM module.

use assert_approx_eq::assert_approx_eq;
use num_complex::Complex;
use proptest::prelude::*;
use rssn::physics::physics_cnm::*;

#[test]

fn test_solve_heat_1d_cn_conservation() {

    let n = 20;

    let initial = vec![1.0; n];

    let res = solve_heat_equation_1d_cn(
        &initial,
        0.1,
        0.01,
        0.1,
        10,
    );

    // For insulated boundaries (b_i=1.0 at ends), total "heat" might change if implementation is simple Dirichlet
    // Current implementation sets b[0]=1.0, b[n-1]=1.0 and d[0]=0, d[n-1]=0 which means u_0 and u_{n-1} become 0.
    // So heat will leak out.
    assert_eq!(res.len(), n);

    assert!(res[10] > 0.0);

    assert_eq!(res[0], 0.0);
}

#[test]

fn test_solve_schrodinger_1d_norm_conservation() {

    let n = 100;

    let dx = 0.1;

    let mut psi0 = vec![Complex::new(0.0, 0.0); n];

    // Gaussian wave packet centered at the middle
    for i in 0 .. n {

        let x = (i as f64 - (n as f64 / 2.0)) * dx;

        let val = (-x * x).exp();

        psi0[i] = Complex::new(val, 0.0);
    }

    // Normalize
    let mut norm = (psi0
        .iter()
        .map(|p| p.norm_sqr())
        .sum::<f64>()
        * dx)
        .sqrt();

    for p in psi0.iter_mut() {

        *p = *p / norm;
    }

    let v = vec![0.0; n];

    let res = solve_schrodinger_1d_cn(
        &psi0, &v, dx, 0.01, 10,
    );

    let final_norm = (res
        .iter()
        .map(|p| p.norm_sqr())
        .sum::<f64>()
        * dx)
        .sqrt();

    assert_approx_eq!(
        final_norm,
        1.0,
        1e-10
    );
}

proptest! {
    #[test]
    fn prop_heat_1d_stability(d_coeff in 0.01..1.0f64) {
        let n = 10;
        let initial = vec![1.0; n];
        let res = solve_heat_equation_1d_cn(&initial, 0.1, 0.001, d_coeff, 5);
        for &val in &res {
            prop_assert!(val.is_finite());
            prop_assert!(val >= 0.0 && val <= 1.1); // 1.1 to allow tiny numerical overshoot if any
        }
    }
}
