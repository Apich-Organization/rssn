//! Scaled reflection and preservation of untouched components.

use rssn::kernels::vector::reflect;

#[test]
fn finite_coordinates_survive_overflow_of_the_reflection_factor() {
    let tiny = f64::from_bits(1);
    let got = reflect(&[f64::MAX, tiny], &[1.0, 0.0]).unwrap();
    assert_eq!(got[0].to_bits(), (-f64::MAX).to_bits());
    assert_eq!(got[1].to_bits(), 1);
}

#[test]
fn reflection_recovers_when_the_dot_itself_overflows() {
    let tiny = f64::from_bits(1);
    let v = [f64::MAX, f64::MAX, f64::MAX, f64::MAX, tiny];
    let result = reflect(&v, &[0.5, 0.5, 0.5, 0.5, 0.0]).unwrap();
    for &value in &result[..4] {
        assert_eq!(value.to_bits(), (-f64::MAX).to_bits());
    }
    assert_eq!(result[4].to_bits(), 1);
}

#[test]
fn overflowing_partial_sums_can_cancel_to_zero() {
    let max = f64::MAX;
    let v = [max, max, max, max, -max, -max, -max, -max];
    let n = [0.5 * std::f64::consts::FRAC_1_SQRT_2; 8];
    let got = reflect(&v, &n).unwrap();
    assert_eq!(v.as_slice(), got.as_slice());
}
