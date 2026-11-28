//! SDIRK analytic references and checked work during Newton factor reuse.

use std::cell::Cell;
use std::sync::Arc;
use std::sync::atomic::AtomicBool;
use std::sync::atomic::Ordering;

use rssn::cancellation;
use rssn::kernels::dense::Mat;
use rssn::kernels::ode_adaptive::OdeError;
use rssn::kernels::ode_sdirk::SdirkOptions;
use rssn::kernels::ode_sdirk::sdirk3;
use rssn::kernels::ode_sdirk::sdirk3_fixed;

#[test]
fn constant_dense_jacobian_keeps_the_analytic_solution_and_evaluation_counts() {
    const N: usize = 8;
    const SHIFT: f64 = 1.25;
    const DECAY: f64 = 5.0;
    let y0 = [0.7, 1.2, 0.5, 0.8, -0.9, 0.8, 1.1, 0.6];
    let direction = [1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0];
    let projection =
        |y: &[f64]| y.iter().zip(direction).map(|(a, b)| a * b).sum::<f64>() / N as f64;
    let matrix = Mat::from_rows(
        &direction
            .iter()
            .enumerate()
            .map(|(i, &row)| {
                direction
                    .iter()
                    .enumerate()
                    .map(|(j, &col)| SHIFT * (row * col) / N as f64 - f64::from(i == j) * DECAY)
                    .collect()
            })
            .collect::<Vec<Vec<f64>>>(),
    )
    .unwrap();
    let rhs_count = Cell::new(0);
    let jacobian_count = Cell::new(0);
    let f = |_: f64, y: &[f64], out: &mut [f64]| {
        rhs_count.set(rhs_count.get() + 1);
        let shared = SHIFT * projection(y);
        for ((value, &state), sign) in out.iter_mut().zip(y).zip(direction) {
            *value = -DECAY * state + sign * shared;
        }
    };
    let jac = |_: f64, _: &[f64]| {
        jacobian_count.set(jacobian_count.get() + 1);
        matrix.clone()
    };
    let options = SdirkOptions {
        h0: Some(0.05),
        ..SdirkOptions::default()
    };
    let time = 0.2;
    let result = sdirk3(f, jac, 0.0, time, &y0, &options).unwrap();
    assert!(result.rejected > 0);
    assert_eq!(result.nfev, rhs_count.get());
    assert_eq!(result.njev, jacobian_count.get());
    assert_eq!(result.newton_iterations, result.njev);
    let parallel = projection(&y0);
    for ((got, initial), sign) in result.y.last().unwrap().iter().zip(y0).zip(direction) {
        let ref_value = (-DECAY * time).exp() * (initial - parallel * sign)
            + ((SHIFT - DECAY) * time).exp() * parallel * sign;
        assert!((got - ref_value).abs() < 8e-7, "{got} != {ref_value}");
    }
}

#[test]
fn moving_jacobians_remain_correct_in_both_time_directions() {
    let f = |t: f64, y: &[f64], out: &mut [f64]| out[0] = -(1.0 + t) * y[0];
    let jac = |t: f64, _: &[f64]| Mat::from_rows(&[vec![-1.0 - t]]).unwrap();
    let exact = (-1.5_f64).exp();
    let forward = sdirk3_fixed(f, jac, 0.0, 1.0, &[1.0], 256).unwrap();
    let backward = sdirk3_fixed(f, jac, 1.0, 0.0, &[exact], 256).unwrap();
    assert!((forward[0] - exact).abs() < 8e-8);
    assert!((backward[0] - 1.0).abs() < 8e-8);
}

#[test]
fn malformed_jacobian_storage_is_rejected_before_factorization() {
    let f = |_: f64, _: &[f64], out: &mut [f64]| out[0] = 1.0;
    for data in [vec![], vec![0.0, 0.0]] {
        let jac = |_: f64, _: &[f64]| {
            Mat {
                rows: 1,
                cols: 1,
                data: data.clone(),
            }
        };
        assert_eq!(
            sdirk3_fixed(f, jac, 0.0, 0.01, &[0.0], 1).unwrap_err(),
            OdeError::Invalid
        );
    }
}

#[test]
fn new_and_reused_factors_observe_cancellation() {
    let f = |_: f64, y: &[f64], out: &mut [f64]| out[0] = -y[0];
    for cancel_on in [1, 2] {
        let token = Arc::new(AtomicBool::new(false));
        let _scope = cancellation::install(Arc::clone(&token));
        let calls = Cell::new(0);
        let jac = |_: f64, _: &[f64]| {
            calls.set(calls.get() + 1);
            if calls.get() == cancel_on {
                token.store(true, Ordering::Relaxed);
            }
            Mat::from_rows(&[vec![-1.0]]).unwrap()
        };
        assert_eq!(
            sdirk3_fixed(f, jac, 0.0, 0.01, &[1.0], 1).unwrap_err(),
            OdeError::Cancelled
        );
        assert_eq!(calls.get(), cancel_on);
    }
}
