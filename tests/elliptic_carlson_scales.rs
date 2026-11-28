//! Scale, principal-value and continuation tests for Carlson integrals.

use std::f64::consts::FRAC_PI_2;
use std::f64::consts::LN_2;
use std::f64::consts::PI;

use rssn::kernels::special_ext::elliptic_f;
use rssn::sim::elliptic::EllipticPiPlan;
use rssn::sim::elliptic::carlson_rc;
use rssn::sim::elliptic::carlson_rj;
use rssn::sim::elliptic::elliptic_pi;
use rssn::sim::elliptic::elliptic_pi_complete;

fn relative(
    actual: f64,
    expected: f64,
    tolerance: f64,
) {
    assert!(
        (actual / expected - 1.0).abs() <= tolerance,
        "{actual:e} versus {expected:e}"
    );
}

fn power(exponent: i32) -> f64 {
    f64::from_bits(((exponent + 1023) as u64) << 52)
}

#[test]
fn rc_stays_finite_across_large_and_subnormal_inputs() {
    let tiny = f64::from_bits(1);
    let huge = 1.0e308_f64;
    relative(carlson_rc(0.0, tiny), FRAC_PI_2 / tiny.sqrt(), 2e-15);
    relative(carlson_rc(huge, huge), 1.0 / huge.sqrt(), 2e-15);
    relative(
        carlson_rc(huge, tiny),
        (LN_2 + 0.5 * (huge.ln() - tiny.ln())) / huge.sqrt(),
        4e-15,
    );
    relative(
        carlson_rc(huge, -huge),
        1.0_f64.asinh() / (2.0_f64.sqrt() * huge.sqrt()),
        4e-15,
    );
    assert_eq!(carlson_rc(0.0, -huge), 0.0);
    assert!(carlson_rc(-1.0, 2.0).is_nan());
    assert!(carlson_rc(1.0, 0.0).is_nan());
}

#[test]
fn rj_homogeneity_survives_raw_cubic_overflow_and_underflow() {
    let base = carlson_rj(0.0, 0.75, 1.25, 2.0);
    for exponent in [-600, -400, 400, 600] {
        let t = power(exponent);
        relative(
            carlson_rj(0.0, 0.75 * t, 1.25 * t, 2.0 * t),
            base * power(-3 * exponent / 2),
            5e-14,
        );
    }
    assert!(carlson_rj(0.0, 0.0, 1.0, 1.0).is_nan());
    assert!(carlson_rj(1.0, 1.0, 1.0, 0.0).is_nan());
    assert_eq!(carlson_rj(f64::MAX, f64::MAX, f64::MAX, f64::MAX), 0.0);
    let tiny = f64::from_bits(1);
    assert_eq!(carlson_rj(tiny, tiny, tiny, tiny), f64::INFINITY);
}

#[test]
fn mixed_rj_orders_keep_small_arguments_and_long_duplication_weights() {
    let x = power(-1000);
    let p = power(1000);
    // Integrating three equal arguments gives 3(1/sqrt(x)-RC(x,p))/(p-x).
    // The omitted RC term is negligible at this exponent gap.
    relative(carlson_rj(x, x, x, p), 3.0 * power(-500), 1e-12);

    let x = power(700);
    let p = power(-700);
    let subnormal = f64::from_bits(1 << 24);
    relative(
        carlson_rj(x, x, x, p),
        3.0 * (701.0 * LN_2 - 1.0) * subnormal,
        2e-10,
    );

    // With one zero and two unit arguments, the two simple poles after
    // t = s^2 give 3*pi/2 * (1/sqrt(p)-1)/(1-p).
    for p in [power(-1000), power(1000)] {
        let exact = 3.0 * FRAC_PI_2 * (1.0 / p.sqrt() - 1.0) / (1.0 - p);
        relative(carlson_rj(0.0, 1.0, 1.0, p), exact, 5e-13);
    }

    let [x, y, z, p] = [power(-400), power(400), power(500), power(-300)];
    let expected = 2.313_355_794_954_748_4e-90;
    for [a, b, c] in [[x, y, z], [z, x, y], [y, z, x]] {
        relative(carlson_rj(a, b, c, p), expected, 5e-13);
    }
}

#[test]
fn prepared_pi_reuses_the_period_over_positive_and_negative_amplitudes() {
    let (n, m) = (-0.75, 0.8);
    let plan = EllipticPiPlan::new(n, m).unwrap();
    assert_eq!(
        plan.complete().to_bits(),
        elliptic_pi_complete(n, m).to_bits()
    );
    let local = 0.44;
    for j in -100..=100 {
        let angle = local + f64::from(j) * PI;
        let expected = plan.evaluate(local) + 2.0 * f64::from(j) * plan.complete();
        assert!((plan.evaluate(angle) - expected).abs() < 2e-13);
        assert_eq!(
            plan.evaluate(angle).to_bits(),
            elliptic_pi(n, angle, m).to_bits()
        );
        assert!((plan.evaluate(angle) + plan.evaluate(-angle)).abs() < 2e-13);
    }
    assert!(EllipticPiPlan::new(1.0, m).is_err());
    assert!(EllipticPiPlan::new(n, 1.0).is_err());
    assert!(EllipticPiPlan::new(f64::NAN, m).is_err());
}

#[test]
fn quarter_period_is_exact_even_near_the_pole() {
    let n = 1.0 - power(-44);
    let m = 1.0 - power(-51);
    let plan = EllipticPiPlan::new(n, m).unwrap();
    assert_eq!(
        plan.evaluate(FRAC_PI_2).to_bits(),
        plan.complete().to_bits()
    );
    assert_eq!(
        plan.evaluate(-FRAC_PI_2).to_bits(),
        (-plan.complete()).to_bits()
    );
    assert_eq!(
        elliptic_pi(n, FRAC_PI_2, m).to_bits(),
        plan.complete().to_bits()
    );
}

#[test]
fn zero_parameter_and_characteristic_use_their_elementary_limits() {
    let n = -1.0e280_f64;
    let k = (1.0 - n).sqrt();
    let plan = EllipticPiPlan::new(n, 0.0).unwrap();
    relative(plan.complete(), FRAC_PI_2 / k, 1e-15);
    relative(plan.evaluate(0.7), FRAC_PI_2 / k, 2e-15);
    relative(plan.evaluate(-PI), -PI / k, 2e-15);
    let phi = 1e-142;
    relative(plan.evaluate(phi), (k * phi).atan() / k, 2e-15);

    for phi in [-7.25_f64, -0.4, 0.1, 4.5] {
        relative(elliptic_pi(0.0, phi, 0.7), elliptic_f(phi, 0.7), 2e-15);
    }
    let phi = 0.4_f64;
    relative(elliptic_pi(1.0, phi, 0.0), phi.tan(), 5e-15);
}

#[test]
fn large_negative_characteristic_avoids_first_kind_cancellation() {
    let n = -1.0e200_f64;
    let k = (1.0 - n).sqrt();
    for m in [-0.45, 0.7] {
        let plan = EllipticPiPlan::new(n, m).unwrap();
        relative(plan.complete(), FRAC_PI_2 / k, 2e-15);
        relative(plan.evaluate(0.7), FRAC_PI_2 / k, 2e-15);
        relative(plan.evaluate(1e-102), (k * 1e-102).atan() / k, 2e-15);
        relative(plan.evaluate(-PI), -PI / k, 2e-15);
    }
}
