//! Exact rounding of arbitrary-precision real division.

use rssn::kernels::mp::Real;
use rssn::kernels::mp::Round;

fn assert_divided_by_one(
    value: i64,
    input_bits: u32,
    output_bits: u32,
    mode: Round,
    expected: i64,
) {
    let numerator = Real::from_i64(value, input_bits).unwrap();
    let result = numerator.div(&Real::one(1), output_bits, mode).unwrap();
    let expected = Real::from_i64(expected, output_bits).unwrap();
    assert_eq!(
        (result.mantissa(), result.exponent()),
        (expected.mantissa(), expected.exponent()),
        "{value} under {mode:?}"
    );
}

#[test]
fn directed_rounding_keeps_bits_discarded_before_division() {
    assert_divided_by_one(1025, 11, 1, Round::TowardPositive, 2048);
    assert_divided_by_one(-1025, 11, 1, Round::TowardNegative, -2048);
    assert_divided_by_one(1025, 11, 1, Round::AwayFromZero, 2048);
    assert_divided_by_one(1025, 11, 1, Round::TowardZero, 1024);
    assert_divided_by_one(1024, 11, 1, Round::TowardPositive, 1024);
}

#[test]
fn nearest_even_distinguishes_discarded_bits_from_an_exact_tie() {
    assert_divided_by_one(5121, 13, 2, Round::NearestEven, 6144);
    assert_divided_by_one(5120, 13, 2, Round::NearestEven, 4096);
    assert_divided_by_one(5119, 13, 2, Round::NearestEven, 4096);
    assert_divided_by_one(-5121, 13, 2, Round::NearestEven, -6144);
}
