//! Scale and fixed-noise regressions for the commutative Platen method.

use rssn::kernels::sde::srk_platen;

fn zero_drift(
    _: f64,
    _: &[f64],
    out: &mut [f64],
) {
    out.fill(0.0);
}

#[test]
fn finite_platen_correction_survives_a_nonrepresentable_product() {
    let linear = |_: f64, state: &[f64], out: &mut [f64]| out[0] = state[0];
    let path = srk_platen(
        zero_drift,
        linear,
        1,
        0.0,
        2.0_f64.powi(1000),
        &[1.0],
        &[vec![0.0]],
    )
    .unwrap();
    assert_eq!(path[1][0], -2.0_f64.powi(999));

    let tiny_affine = |_: f64, state: &[f64], out: &mut [f64]| {
        out[0] = 2.0_f64.powi(-500) + 2.0_f64.powi(500) * state[0];
    };
    let path = srk_platen(
        zero_drift,
        tiny_affine,
        1,
        0.0,
        2.0_f64.powi(-1000),
        &[0.0],
        &[vec![0.0]],
    )
    .unwrap();
    assert_eq!(path[1][0], -2.0_f64.powi(-1001));
}

#[test]
fn diagonal_brownian_weight_keeps_a_large_square_or_tiny_cancellation() {
    let affine = |_: f64, state: &[f64], out: &mut [f64]| {
        out[0] = 2.0_f64.powi(-512) * (1.0 + state[0]);
    };
    let path = srk_platen(
        zero_drift,
        affine,
        1,
        0.0,
        2.0_f64.powi(1022),
        &[0.0],
        &[vec![2.0_f64.powi(513)]],
    )
    .unwrap();
    // The linear contribution is 2 and the weighted correction is 15/8.
    assert_eq!(path[1][0], 3.875);

    let affine = |_: f64, state: &[f64], out: &mut [f64]| {
        out[0] = 1.0 + 2.0_f64.powi(537) * state[0];
    };
    let root_h = 2.0_f64.powi(-537);
    let dw = root_h.next_up();
    let path = srk_platen(
        zero_drift,
        affine,
        1,
        0.0,
        f64::from_bits(1),
        &[0.0],
        &[vec![dw]],
    )
    .unwrap();
    // The exact diagonal correction is one noise ULP plus 2^-642.
    assert_eq!(path[1][0], dw.next_up());
}

#[test]
fn cross_brownian_weight_preserves_commutative_noise_below_the_float_range() {
    // The two columns commute. Their mixed directional derivatives are
    // (0, 0, b) in both orders, and their self derivatives vanish.
    let diffusion = |_: f64, state: &[f64], out: &mut [f64]| {
        let b = 2.0_f64.powi(1000);
        out.copy_from_slice(&[1.0, 0.0, 0.0, 1.0, b * state[1], b * state[0]]);
    };
    let dw = 2.0_f64.powi(-550);
    let path = srk_platen(
        zero_drift,
        diffusion,
        2,
        0.0,
        2.0_f64.powi(-1000),
        &[0.0; 3],
        &[vec![dw, dw]],
    )
    .unwrap();
    assert_eq!(path[1], [dw, dw, 2.0_f64.powi(-100)]);
}

#[test]
fn weighted_difference_between_opposite_finite_coefficients() {
    let diffusion = |_: f64, state: &[f64], out: &mut [f64]| {
        out[0] = -f64::MAX * (2.0 * state[0].exp() - 1.0);
    };
    let path = srk_platen(
        zero_drift,
        diffusion,
        1,
        0.0,
        2.0_f64.powi(-1000),
        &[0.0],
        &[vec![0.0]],
    )
    .unwrap();
    assert_eq!(path[1][0], -f64::MAX * 2.0_f64.powi(-500));
}

#[test]
fn empty_increment_sequence_needs_no_diffusion_workspace() {
    let unused = |_: f64, _: &[f64], _: &mut [f64]| panic!("empty path called a coefficient");
    let path = srk_platen(unused, unused, usize::MAX, 0.0, 1.0, &[1.0, 2.0], &[]).unwrap();
    assert_eq!(path, [vec![1.0, 2.0]]);
}
