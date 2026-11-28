//! High-range checks of the flexible-cell matrix propagators.

use rssn::sim::md_engine::mat_exp_phi1;

fn close_ulp(
    actual: f64,
    expected_bits: u64,
    max_ulp: u64,
) {
    let expected = f64::from_bits(expected_bits);
    if expected.is_infinite() {
        assert_eq!(actual.to_bits(), expected_bits);
    } else {
        assert!(
            actual.is_finite() && actual.to_bits().abs_diff(expected_bits) <= max_ulp,
            "{actual:.17e} vs {expected:.17e}"
        );
    }
}

#[test]
fn scalar_propagators_match_high_precision_references_across_the_float_range() {
    // Decimal exp and (exp(x) - 1) / x at 420 digits, rounded to f64.
    let references = [
        (-f64::MAX, 0, 0x0004_0000_0000_0000),
        (-1e308, 0, 0x0007_30d6_7819_e8d2),
        (-1000.0, 0, 0x3f50_624d_d2f1_a9fc),
        (-0.7, 0x3fdf_c80d_b9dd_5542, 0x3fe7_0363_e8f4_30d1),
        (-1e-16, 0x3fef_ffff_ffff_ffff, 0x3ff0_0000_0000_0000),
        (
            -f64::from_bits(1),
            0x3ff0_0000_0000_0000,
            0x3ff0_0000_0000_0000,
        ),
        (0.0, 0x3ff0_0000_0000_0000, 0x3ff0_0000_0000_0000),
        (
            f64::from_bits(1),
            0x3ff0_0000_0000_0000,
            0x3ff0_0000_0000_0000,
        ),
        (1e-16, 0x3ff0_0000_0000_0000, 0x3ff0_0000_0000_0000),
        (0.3, 0x3ff5_9905_8c8c_1a96, 0x3ff2_a8bd_29d3_0348),
        (25.0, 0x4230_c3d3_9209_62c9, 0x41e5_7589_b0ae_92ec),
        (709.0, 0x7fdd_422d_2be5_dc9b, 0x7f45_20fc_26ec_0a8c),
        (710.0, 0x7ff0_0000_0000_0000, 0x7f5c_ad37_ec32_8617),
        (714.0, 0x7ff0_0000_0000_0000, 0x7fb8_53ad_e482_2d24),
        (716.35, 0x7ff0_0000_0000_0000, 0x7fef_c7d1_6ebc_557e),
        (716.36, 0x7ff0_0000_0000_0000, 0x7ff0_0000_0000_0000),
        (1000.0, 0x7ff0_0000_0000_0000, 0x7ff0_0000_0000_0000),
    ];
    for (x, expected_e, expected_phi) in references {
        let (e, phi) = mat_exp_phi1(&[[x]]);
        close_ulp(e[0][0], expected_e, 2);
        close_ulp(phi[0][0], expected_phi, 3);
    }
}

#[test]
fn an_overflowing_diagonal_exponential_leaves_other_components_independent() {
    let a = [[716.35, 0.0, -0.0], [0.0, -f64::MAX, 0.0], [0.0, 0.0, 0.0]];
    let (e, phi) = mat_exp_phi1(&a);
    assert_eq!(e[0][0].to_bits(), f64::INFINITY.to_bits());
    assert_eq!(e[1][1].to_bits(), 0);
    assert_eq!(e[2][2].to_bits(), 1.0_f64.to_bits());
    close_ulp(phi[0][0], 0x7fef_c7d1_6ebc_557e, 3);
    close_ulp(phi[1][1], 0x0004_0000_0000_0000, 0);
    assert_eq!(phi[2][2].to_bits(), 1.0_f64.to_bits());
    for i in 0..3 {
        for j in 0..3 {
            if i != j {
                assert_eq!(e[i][j].to_bits(), 0);
                assert_eq!(phi[i][j].to_bits(), 0);
            }
        }
    }
}

#[test]
fn nilpotent_simple_shears_use_the_exact_truncated_series() {
    for matrix in [
        [[0.0, f64::MAX, -f64::MAX], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
        [[0.0, 0.0, 0.0], [f64::MAX, 0.0, 0.0], [-f64::MAX, 0.0, 0.0]],
        [
            [0.0, 0.03125, -f64::from_bits(1)],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
        ],
    ] {
        let (e, phi) = mat_exp_phi1(&matrix);
        for (i, row) in matrix.iter().enumerate() {
            for (j, &a) in row.iter().enumerate() {
                // A² = 0 structurally, including these multiple-shear blocks.
                let base = if i == j { 1.0_f64 } else { 0.0 };
                assert_eq!(e[i][j].to_bits(), (base + a).to_bits());
                assert_eq!(phi[i][j].to_bits(), (base + a * 0.5).to_bits());
            }
        }
    }
}
