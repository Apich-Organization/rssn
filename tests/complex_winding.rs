//! Scale, geometry, and sampling regressions for complex winding counts.

use std::cell::Cell;

use num_complex::Complex64;
use rssn::kernels::complex::count_zeros_poles;
use rssn::kernels::complex::winding_number;

#[test]
fn argument_principle_is_independent_of_the_callback_scale() {
    for scale in [1.0, 1e200, 1e-200, f64::from_bits(16)] {
        let count = |f: fn(Complex64) -> Complex64| {
            count_zeros_poles(|z| f(z) * scale, Complex64::ZERO, 1.0, 64)
        };
        assert_eq!(count(|z| z), Some(1));
        assert_eq!(count(|z| Complex64::ONE / z), Some(-1));
    }
}

#[test]
fn samples_once_per_node_and_rejects_the_first_invalid_value() {
    let visits = Cell::new(0);
    let wind = count_zeros_poles(
        |z| {
            visits.set(visits.get() + 1);
            z
        },
        Complex64::ZERO,
        1.0,
        1,
    );
    assert_eq!(wind, Some(1));
    assert_eq!(visits.get(), 64);

    for invalid in [Complex64::ZERO, Complex64::new(f64::INFINITY, 0.0)] {
        visits.set(0);
        let wind = count_zeros_poles(
            |_| {
                visits.set(visits.get() + 1);
                invalid
            },
            Complex64::ZERO,
            1.0,
            64,
        );
        assert_eq!(wind, None);
        assert_eq!(visits.get(), 1);
    }
}

#[test]
fn sampled_count_allows_exact_quarter_turns_and_rejects_larger_steps() {
    let cardinals = [
        Complex64::new(1.0, 0.0),
        Complex64::new(0.0, 1.0),
        Complex64::new(-1.0, 0.0),
        Complex64::new(0.0, -1.0),
    ];
    let visits = Cell::new(0);
    let wind = count_zeros_poles(
        |_| {
            let k = visits.get();
            visits.set(k + 1);
            cardinals[k / 16]
        },
        Complex64::ZERO,
        1.0,
        64,
    );
    assert_eq!(wind, Some(1));
    assert_eq!(
        count_zeros_poles(|z| z.powi(40), Complex64::ZERO, 1.0, 64),
        None
    );

    visits.set(0);
    let e = f64::EPSILON;
    let first = Complex64::new(1.0 + e, 1.0);
    let second = Complex64::new(1.0 - e, -1.0);
    // The dot product is -e^2, although the rounded products sum to zero.
    let wind = count_zeros_poles(
        |_| {
            let k = visits.get();
            visits.set(k + 1);
            if k == 0 { first } else { second }
        },
        Complex64::ZERO,
        1.0,
        64,
    );
    assert_eq!(wind, None);
    assert_eq!(visits.get(), 2);
}

#[test]
fn polygon_handles_orientation_and_degenerate_edges() {
    let triangle = [
        Complex64::new(1.0, 0.0),
        Complex64::new(-0.5, 1.0),
        Complex64::new(-0.5, -1.0),
    ];
    assert_eq!(winding_number(&triangle), Some(1));
    assert_eq!(
        winding_number(&[triangle[2], triangle[1], triangle[0]]),
        Some(-1)
    );
    assert_eq!(winding_number(&triangle.map(|p| p + 2.0)), Some(0));
    assert_eq!(winding_number(&[triangle[0], triangle[0]]), Some(0));
    assert_eq!(winding_number(&[triangle[0]]), Some(0));
    assert_eq!(winding_number(&[]), None);
    assert_eq!(winding_number(&[triangle[0], Complex64::ZERO]), None);
    assert_eq!(winding_number(&[triangle[0], -triangle[0]]), None);
    assert_eq!(winding_number(&[Complex64::new(f64::NAN, 0.0)]), None);
}

#[test]
fn polygon_uses_exact_crossings_when_float_products_tie_or_leave_range() {
    for (x, y) in [
        (2.0_f64.powi(900), 2.0_f64.powi(-900)),
        (f64::from_bits(1), f64::from_bits(1)),
        (f64::MAX, f64::MAX),
    ] {
        let rectangle = [
            Complex64::new(x, -y),
            Complex64::new(x, y),
            Complex64::new(-x, y),
            Complex64::new(-x, -y),
        ];
        assert_eq!(winding_number(&rectangle), Some(1));
        assert_eq!(winding_number(&[rectangle[1], rectangle[3]]), None);
    }
    let e = f64::EPSILON;
    // The first edge's determinant is e^2 even though both rounded
    // product magnitudes are 1. It crosses the real axis to the left.
    let near_antipodal = [
        Complex64::new(1.0 + e, 1.0),
        Complex64::new(-1.0, -(1.0 - e)),
        Complex64::new(0.0, -1.0),
    ];
    assert_eq!(winding_number(&near_antipodal), Some(1));
}
