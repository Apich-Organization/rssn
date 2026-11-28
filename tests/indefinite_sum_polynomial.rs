use rssn::kernels::indefinite_sum::PolynomialAntidifference;
use rssn::kernels::indefinite_sum::series_antidiff;

#[test]
fn normalized_quadratic_retains_near_anchor_values() {
    let coefficients = [0.0, 1.0];
    let plan = PolynomialAntidifference::new(&coefficients, 0.0);
    for x in [2.0_f64.powi(-60), -2.0_f64.powi(-60), f64::MIN_POSITIVE] {
        // The unique polynomial with F(x+1) - F(x) = x and F(0) = 0.
        let expected = 0.5 * x * (x - 1.0);
        assert_eq!(
            series_antidiff(&coefficients, 0.0, x).to_bits(),
            expected.to_bits()
        );
        assert_eq!(plan.evaluate(x).to_bits(), expected.to_bits());
        assert_eq!(plan.difference(0.0, x).to_bits(), expected.to_bits());
    }
    assert_eq!(plan.evaluate(0.0).to_bits(), 0.0_f64.to_bits());
}

#[test]
fn factored_polynomial_retains_small_values_at_distinct_scales() {
    let x = 2.0_f64.powi(-60);
    let cubic = PolynomialAntidifference::new(&[0.0, 0.0, 0.0, 1.0], 0.0);
    // The sum of cubes extends the square of the triangular polynomial.
    let expected = (0.5 * x * (x - 1.0)).powi(2);
    assert_eq!(cubic.evaluate(x).to_bits(), expected.to_bits());
    assert_eq!(cubic.difference(0.0, x).to_bits(), expected.to_bits());

    let mixed_scale = PolynomialAntidifference::new(&[f64::MAX, -f64::MAX], 0.0);
    let expected = f64::MAX * x + (-0.5 * f64::MAX * x) * (x - 1.0);
    assert!(expected.is_finite());
    assert_eq!(mixed_scale.evaluate(x).to_bits(), expected.to_bits());
    assert_eq!(mixed_scale.difference(0.0, x).to_bits(), expected.to_bits());
}

#[test]
fn nearby_endpoints_retain_their_physical_separation() {
    let plan = PolynomialAntidifference::new(&[0.0, 1.0], 0.0);
    let from = 1.5_f64;
    let to = from.next_up();
    let expected = 0.5 * (to - from) * (to + from - 1.0);
    assert_eq!(plan.difference(from, to).to_bits(), expected.to_bits());
    assert_eq!(plan.difference(to, from).to_bits(), (-expected).to_bits());

    // Subtracting a remote center rounds both transformed endpoints to
    // the same float. For a constant summand F(to) - F(from) is still to-from.
    let constant = PolynomialAntidifference::new(&[1.0], 2.0_f64.powi(80));
    assert_eq!(
        constant.difference(from, to).to_bits(),
        (to - from).to_bits()
    );
    assert_eq!(constant.difference(from, from).to_bits(), 0.0_f64.to_bits());
}

#[test]
fn normalized_unit_step_preserves_the_constant_input_amid_large_terms() {
    let coefficients = [1.25, 1e30, -8e28, 1e27];
    let center = 7.0;
    let from = center;
    let to = center + 1.0;
    let plan = PolynomialAntidifference::new(&coefficients, center);

    assert_eq!(
        series_antidiff(&coefficients, center, to).to_bits(),
        coefficients[0].to_bits()
    );
    assert_eq!(plan.evaluate(to).to_bits(), coefficients[0].to_bits());
    assert_eq!(
        plan.difference(from, to).to_bits(),
        coefficients[0].to_bits()
    );
    assert_eq!(
        plan.difference(to, from).to_bits(),
        (-coefficients[0]).to_bits()
    );
}

#[test]
fn prepared_polynomial_agrees_with_direct_lattice_sums_and_recurrence() {
    let coefficients = [0.5, -0.25, 0.0, 0.75, -0.125, 1.0 / 16.0];
    let center = -3.25;
    let plan = PolynomialAntidifference::new(&coefficients, center);
    let summand = |u: f64| 0.5 - 0.25 * u + 0.75 * u.powi(3) - 0.125 * u.powi(4) + u.powi(5) / 16.0;
    for steps in -12..=12 {
        let expected = if steps < 0 {
            -(steps..0).map(|k| summand(f64::from(k))).sum::<f64>()
        } else {
            (0..steps).map(|k| summand(f64::from(k))).sum::<f64>()
        };
        let x = center + f64::from(steps);
        let scale = 1.0_f64.max(expected.abs());
        assert!((plan.evaluate(x) - expected).abs() <= 16.0 * f64::EPSILON * scale);
        assert!(
            (series_antidiff(&coefficients, center, x) - expected).abs()
                <= 16.0 * f64::EPSILON * scale
        );
    }
    for u in [-4.25_f64, -1.25, 0.125, 0.75, 3.125, 9.5] {
        let x = center + u;
        let expected = summand(u);
        let scale = 1.0_f64.max(expected.abs());
        assert!((plan.difference(x, x + 1.0) - expected).abs() <= 16.0 * f64::EPSILON * scale);
    }

    let empty = PolynomialAntidifference::new(&[], center);
    assert_eq!(empty.evaluate(1.0).to_bits(), 0.0_f64.to_bits());
    assert_eq!(empty.difference(0.0, 1.0).to_bits(), 0.0_f64.to_bits());
}
