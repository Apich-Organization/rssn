//! Scale and exterior-domain checks for Kerr-Newman characteristic radii.

use rssn::sim::relativity::KerrNewman;

fn check_ratio(
    radius: Option<f64>,
    mass: f64,
    expected: Option<f64>,
) {
    match (radius, expected) {
        | (Some(radius), Some(expected)) => {
            let ratio = radius / mass;
            assert!(
                (ratio - expected).abs() < 16.0 * f64::EPSILON,
                "r/M={ratio:.17e}, expected {expected:.17e}"
            );
        },
        | (None, None) => {},
        | _ => panic!("radius {radius:?}, mass {mass}, expected ratio {expected:?}"),
    }
}

#[test]
fn ordinary_photon_radii_match_independent_high_precision_metric_roots() {
    // Reference values from the signed null-circle equation with 200-digit arithmetic,
    // checked against g_tt + 2 Omega g_tphi + Omega^2 g_phiphi at each root.
    for (m, a, q, horizon, positive, negative) in [
        (
            3.0,
            1.25,
            -1.75,
            1.697_216_688_778_396_4,
            2.165_935_603_919_166_4,
            3.227_865_673_503_483_6,
        ),
        (
            4.0,
            2.0,
            3.0,
            1.433_012_701_892_219_4,
            1.720_804_296_176_275_5,
            3.153_698_926_392_698,
        ),
        (
            7.8125,
            -2.01,
            5.67,
            1.638_027_761_151_503_5,
            2.911_060_457_601_111,
            2.225_216_472_270_077_4,
        ),
    ] {
        for scale in [
            2.0_f64.powi(-900),
            2.0_f64.powi(-500),
            1.0,
            2.0_f64.powi(600),
            2.0_f64.powi(900),
        ] {
            let metric = KerrNewman {
                m: m * scale,
                a: a * scale,
                q: q * scale,
            };
            check_ratio(metric.horizon(), metric.m, Some(horizon));
            check_ratio(metric.photon_orbit(true), metric.m, Some(positive));
            check_ratio(metric.photon_orbit(false), metric.m, Some(negative));
            let charge_flipped = KerrNewman {
                q: -metric.q,
                ..metric
            };
            assert_eq!(metric.horizon(), charge_flipped.horizon());
            assert_eq!(metric.photon_orbit(true), charge_flipped.photon_orbit(true));
            let spin_flipped = KerrNewman {
                a: -metric.a,
                ..metric
            };
            assert_eq!(metric.photon_orbit(true), spin_flipped.photon_orbit(false));
        }
    }
}

#[test]
fn near_extremality_distinguishes_neighboring_binary_inputs() {
    for scale in [2.0_f64.powi(-900), 1.0, 2.0_f64.powi(900)] {
        for (a, q, horizon, positive, negative) in [
            (
                0.6,
                0.8_f64.next_down(),
                1.000_000_011_542_39,
                1.000_000_020_880_968_5,
                3.2,
            ),
            (
                0.6_f64.next_down(),
                0.8,
                1.000_000_009_424_321_8,
                1.000_000_017_049_239_4,
                3.199_999_999_999_999_7,
            ),
            (
                1.0_f64.next_down(),
                0.0,
                1.000_000_014_901_161_2,
                1.000_000_017_206_379,
                4.0,
            ),
            (
                0.125,
                -0.992_156_741_649_221_4,
                1.000_000_011_260_114_8,
                1.750_000_000_000_000_4,
                2.25,
            ),
        ] {
            let metric = KerrNewman {
                m: scale,
                a: a * scale,
                q: q * scale,
            };
            check_ratio(metric.horizon(), metric.m, Some(horizon));
            check_ratio(metric.photon_orbit(true), metric.m, Some(positive));
            check_ratio(metric.photon_orbit(false), metric.m, Some(negative));
            assert!(metric.photon_orbit(true).unwrap() > metric.horizon().unwrap());
        }
        let super_extremal = KerrNewman {
            m: scale,
            a: 0.8 * scale,
            q: 0.6 * scale,
        };
        assert_eq!(super_extremal.horizon(), None);
        assert_eq!(super_extremal.photon_orbit(true), None);
        assert_eq!(super_extremal.photon_orbit(false), None);
    }
}

#[test]
fn extremal_exterior_branch_is_distinct_from_the_horizon() {
    let cases: [(f64, f64, f64, Option<f64>, f64); 4] = [
        (1.0, 1.0, 0.0, None, 4.0),
        (5.0, 3.0, 4.0, None, 16.0),
        (13.0, 5.0, 12.0, Some(16.0), 36.0),
        (1.0, 0.0, 1.0, Some(2.0), 2.0),
    ];
    for (m, a, q, positive, negative) in cases {
        let metric = KerrNewman { m, a, q };
        assert_eq!(metric.horizon(), Some(m));
        check_ratio(metric.photon_orbit(true), m, positive.map(|r| r / m));
        check_ratio(metric.photon_orbit(false), m, Some(negative / m));
    }

    let quantum = f64::from_bits(1);
    let metric = KerrNewman {
        m: 13.0 * quantum,
        a: 5.0 * quantum,
        q: 12.0 * quantum,
    };
    assert_eq!(metric.horizon(), Some(13.0 * quantum));
    assert_eq!(metric.photon_orbit(true), Some(16.0 * quantum));
    assert_eq!(metric.photon_orbit(false), Some(36.0 * quantum));
}

#[test]
fn radii_are_finite_or_absent_for_extreme_or_invalid_metrics() {
    let tiny = f64::from_bits(1);
    let small = KerrNewman::schwarzschild(tiny);
    assert_eq!(small.horizon(), Some(2.0 * tiny));
    assert_eq!(small.photon_orbit(true), Some(3.0 * tiny));
    let big = KerrNewman::schwarzschild(2.0_f64.powi(1021));
    assert_eq!(big.horizon(), Some(2.0 * big.m));
    assert_eq!(big.photon_orbit(true), Some(3.0 * big.m));

    let too_big = KerrNewman::schwarzschild(2.0_f64.powi(1023));
    assert_eq!(too_big.horizon(), None);
    assert_eq!(too_big.photon_orbit(true), None);
    let super_extremal = KerrNewman {
        m: 1.0,
        a: 1.0,
        q: tiny,
    };
    assert_eq!(super_extremal.horizon(), None);
    assert_eq!(super_extremal.photon_orbit(false), None);
    for invalid in [
        KerrNewman::schwarzschild(-1.0),
        KerrNewman::schwarzschild(f64::INFINITY),
        KerrNewman {
            m: 1.0,
            a: f64::NAN,
            q: 0.0,
        },
    ] {
        assert_eq!(invalid.horizon(), None);
        assert_eq!(invalid.photon_orbit(true), None);
    }
    let flat = KerrNewman::schwarzschild(0.0);
    assert_eq!(flat.horizon(), Some(0.0));
    assert_eq!(flat.photon_orbit(true), None);
}
