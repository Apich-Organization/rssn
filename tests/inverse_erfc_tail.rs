//! Multiprecision reference values around the inverse complementary error function's tail.

use rssn::kernels::special::inverse_erfc;

#[test]
fn inverse_erfc_resolves_the_smallest_positive_probability() {
    let smallest = f64::from_bits(1);
    let reference = f64::from_bits(0x403b_369a_6244_e684);
    let result = inverse_erfc(smallest);
    assert!(result.to_bits().abs_diff(reference.to_bits()) <= 2);
    assert!(inverse_erfc(f64::from_bits(2)) < result);
}

#[test]
fn inverse_erfc_matches_both_sides_of_tail_transitions() {
    // These probabilities are erfc(x) rounded from 100-digit arithmetic.
    for (p, quantile) in [
        (1.318_675_289_407_676_7e-29, 7.99),
        (1.122_429_717_298_292_6e-29, 8.0),
        (9.551_998_261_605_015e-30, 8.01),
        (3.208_299_158_374_876e-113, 15.99),
        (2.328_485_751_571_530_8e-113, 16.0),
        (1.689_606_691_894_410_7e-113, 16.01),
        (1.368_876_944_999_886_5e-273, 24.99),
        (8.300_172_571_196_523e-274, 25.0),
        (5.031_796_228_007_791e-274, 25.01),
    ] {
        let quantile: f64 = quantile;
        let result = inverse_erfc(p);
        assert!(
            result.to_bits().abs_diff(quantile.to_bits()) <= 2,
            "at p={p:e}, expected {quantile}, got {result}"
        );
    }
}
