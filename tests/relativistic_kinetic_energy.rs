//! Numerical boundary checks for relativistic kinetic energy in SI.

use rssn::sim::classical::SPEED_OF_LIGHT;
use rssn::sim::classical::lorentz_factor;
use rssn::sim::classical::relativistic_kinetic_energy;

fn check(
    value: f64,
    reference: f64,
) {
    let magnitude = value.abs().max(reference.abs());
    let tolerance = (8.0 * f64::EPSILON * magnitude).max(4.0 * f64::from_bits(1));
    assert!(
        value.is_finite() && (value - reference).abs() <= tolerance,
        "{value:.17e} vs {reference:.17e}"
    );
}

#[test]
fn low_speed_remains_resolved_when_lorentz_factor_rounds_to_one() {
    assert_eq!(lorentz_factor(1.0), 1.0);
    assert_eq!(relativistic_kinetic_energy(1.0, 1.0), 0.5);
    assert_eq!(relativistic_kinetic_energy(1.0, -1.0), 0.5);
    check(
        relativistic_kinetic_energy(2.5, 60_000.0),
        4_500_000_135.186_986,
    );
    check(
        relativistic_kinetic_energy(2.5, 0.6 * SPEED_OF_LIGHT),
        5.617_219_867_105_109e16,
    );
}

#[test]
fn near_light_speed_uses_the_actual_si_velocity_gap() {
    // References evaluate (1 - v²/c²)^-1/2 and (gamma - 1)mc² at 840 decimal
    // digits from the exact binary64 velocity and the exact integer value of c.
    for (velocity, gamma, energy) in [
        (
            SPEED_OF_LIGHT.next_down(),
            50_148_194.499_089_04,
            4.507_094_861_160_229e24,
        ),
        (
            SPEED_OF_LIGHT.next_down().next_down(),
            35_460_128.394_567_78,
            3.186_997_313_453_509e24,
        ),
    ] {
        check(lorentz_factor(velocity), gamma);
        check(lorentz_factor(-velocity), gamma);
        check(relativistic_kinetic_energy(1.0, velocity), energy);
        check(relativistic_kinetic_energy(1.0, -velocity), energy);
    }
}

#[test]
fn tiny_speed_and_large_mass_can_have_a_representable_energy() {
    let mass = 2.0_f64.powi(1000);
    let velocity = 2.0_f64.powi(-700);
    assert_eq!(velocity * velocity, 0.0);
    assert_eq!(
        relativistic_kinetic_energy(mass, velocity),
        2.0_f64.powi(-401)
    );
    assert_eq!(
        relativistic_kinetic_energy(-mass, velocity),
        -2.0_f64.powi(-401)
    );
    assert!((1.25_f64 * 1.25 * f64::MAX).is_infinite());
    check(
        relativistic_kinetic_energy(f64::MAX, 1.25),
        1.404_447_761_611_184_1e308,
    );

    let quantum = f64::from_bits(1);
    check(
        relativistic_kinetic_energy(quantum, 2.0_f64.powi(26)),
        1.156_174_385_898_489e-308,
    );
    assert_eq!(relativistic_kinetic_energy(1.0, quantum).to_bits(), 0);
    assert_eq!(
        relativistic_kinetic_energy(-1.0, quantum).to_bits(),
        (-0.0_f64).to_bits()
    );
}

#[test]
fn exceptional_values_retain_the_ieee_return_values() {
    for speed in [SPEED_OF_LIGHT, SPEED_OF_LIGHT.next_up(), f64::INFINITY] {
        for velocity in [speed, -speed] {
            assert_eq!(lorentz_factor(velocity), f64::INFINITY);
            assert_eq!(relativistic_kinetic_energy(1.0, velocity), f64::INFINITY);
            assert_eq!(
                relativistic_kinetic_energy(-1.0, velocity),
                f64::NEG_INFINITY
            );
            assert!(relativistic_kinetic_energy(0.0, velocity).is_nan());
        }
    }
    assert!(lorentz_factor(f64::NAN).is_nan());
    assert!(relativistic_kinetic_energy(1.0, f64::NAN).is_nan());
    assert!(relativistic_kinetic_energy(f64::NAN, 1.0).is_nan());
    assert!(relativistic_kinetic_energy(f64::INFINITY, 0.0).is_nan());
    assert_eq!(
        relativistic_kinetic_energy(-1.0, 0.0).to_bits(),
        (-0.0_f64).to_bits()
    );
}
