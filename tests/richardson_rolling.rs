//! Rolling tables and large finite inputs in Richardson extrapolation.

use rssn::kernels::convergence::richardson_extrapolation;
use rssn::kernels::convergence::richardson_extrapolation_with;

#[test]
fn handles_empty_singleton_and_two_error_orders() {
    assert!(richardson_extrapolation(&[]).is_empty());
    assert_eq!(
        richardson_extrapolation(&[-0.0])[0].to_bits(),
        (-0.0_f64).to_bits()
    );
    let values = [11.0, 4.25, 3.265625, 3.0634765625];
    let diagonal = richardson_extrapolation(&values);
    assert_eq!(diagonal.len(), values.len());
    assert_eq!(diagonal[0], values[0]);
    assert!((diagonal[2] - 3.0).abs() < 1e-14);
    assert!((diagonal[3] - 3.0).abs() < 1e-14);
}

#[test]
fn finite_quotients_survive_overflowing_fused_numerators() {
    let same = 1e308;
    let constant = richardson_extrapolation(&[same, same]);
    assert_eq!(constant, [same, same]);
    let opposite = richardson_extrapolation(&[-same, same]);
    assert!((opposite[1] / same - 5.0 / 3.0).abs() < 4e-16);
    let unrepresentable = richardson_extrapolation(&[-f64::MAX, f64::MAX]);
    assert_eq!(unrepresentable[1], f64::INFINITY);
}

#[test]
fn overflowing_growth_power_can_have_a_representable_correction() {
    let growth = 2.0_f64.powi(520);
    let large = 2.0_f64.powi(1000);
    let small = 2.0_f64.powi(480);
    // The second-order cells are A and zero. Subtracting A over 2^1040
    // from zero must retain the small nonzero result.
    let diagonal = richardson_extrapolation_with(&[large, large, small], growth);
    assert_eq!(diagonal, [large, large, -2.0_f64.powi(-40)]);

    let input = vec![1e308; 1181];
    assert_eq!(richardson_extrapolation(&input), input);
}

#[test]
fn long_diagonal_retains_every_input_when_corrections_round_to_zero() {
    let values: Vec<_> = (0..32_481)
        .map(|i| 2.0 + (i as f64 / 32.0).sin() / 1024.0)
        .collect();
    let diagonal = richardson_extrapolation_with(&values, 2.0_f64.powi(80));
    assert_eq!(diagonal, values);
}
