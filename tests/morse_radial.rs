//! Scaled and near-equilibrium Morse energy and forces.

use rssn::sim::physics_md::Particle;
use rssn::sim::physics_md::morse_interaction;
use rssn::sim::physics_md::morse_radial;

fn particle(position: Vec<f64>) -> Particle {
    let dim = position.len();
    Particle::new(0, 1.0, position, vec![0.0; dim])
}

fn close(
    actual: f64,
    expected: f64,
) {
    if expected == 0.0 {
        assert_eq!(actual.abs().to_bits(), 0);
    } else {
        assert!(
            (actual / expected - 1.0).abs() <= 4e-13,
            "observed {actual:.18e}, expected {expected:.18e}"
        );
    }
}

#[test]
fn scaled_radial_values_match_a_decimal_reference() {
    // Independent 460-digit Decimal evaluation of the exact binary64 inputs.
    let cases = [
        (
            [1.0, 0.75, 1e-7, 1.0_f64.next_down()],
            [9.244_463_733_058_732e-47, -1.665_334_536_937_734_7e-30],
        ),
        (
            [1.0_f64.next_up(), 1e308, 2.0, 1.0],
            [1.972_152_263_052_528_8e277, -1.776_356_839_400_249_4e293],
        ),
        (
            [1.0, 1e-320, 2.0, 361.0],
            [2.421_302_808_287_985_6e305, 9.685_211_233_151_942e305],
        ),
        (
            [800.0, 1e307, 1.0, 0.0],
            [1e307, -7.335_749_168_355_375e-41],
        ),
        (
            [1.0, f64::from_bits(1), 0.75, 601.0],
            [3.620_915_331_994_829e67, 5.431_372_997_992_244e67],
        ),
        (
            [1e-308, 1e308, 1e-20, 2e-308],
            [0.0, 1.999_999_999_999_999_6e-40],
        ),
        (
            [1e300, 1e100, 1e-300, 0.0],
            [3.995_764_008_937_281e99, -4.650_883_158_696_593e-201],
        ),
    ];
    for ([r, de, a, re], [energy, force]) in cases {
        let value = morse_radial(r, de, a, re).unwrap();
        close(value.0, energy);
        close(value.1, force);
    }
}

#[test]
fn signed_radial_force_is_the_negative_energy_derivative() {
    let h = 1e-5;
    for r in [0.75, 0.999, 1.001, 1.375] {
        let eval = |r| morse_radial(r, 2.0, 1.25, 1.0).unwrap();
        let (_, force) = eval(r);
        let derivative = (eval(r + h).0 - eval(r - h).0) / (2.0 * h);
        assert!((force + derivative).abs() < 1e-8 * (1.0 + force.abs()));
        assert!(force * (r - 1.0) < 0.0);
    }
    let at_rest = morse_radial(1.0, 1e308, 2.0, 1.0).unwrap();
    close(at_rest.0, 0.0);
    close(at_rest.1, 0.0);
}

#[test]
fn radial_origin_and_particle_collision_have_distinct_contracts() {
    let radial = morse_radial(0.0, 0.5, 2.0, 0.25).unwrap();
    close(radial.0, 0.210_419_643_529_394_46);
    close(radial.1, 2.139_121_115_517_834_2);

    let origin = particle(vec![0.0; 3]);
    for r in [0.0, 1e-9_f64.next_down()] {
        let near = particle(vec![r, 0.0, 0.0]);
        let (energy, force) = morse_interaction(&origin, &near, 0.5, 2.0, 0.25).unwrap();
        assert!(energy.is_infinite() && energy.is_sign_positive());
        assert_eq!(force, vec![0.0; 3]);
    }
    let threshold = particle(vec![1e-9, 0.0, 0.0]);
    assert!(
        morse_interaction(&origin, &threshold, 0.5, 2.0, 0.25)
            .unwrap()
            .0
            .is_finite()
    );
    assert!(morse_interaction(&origin, &particle(vec![1.0]), 1.0, 1.0, 1.0).is_err());
}

#[test]
fn extreme_vector_scaling_keeps_representable_force_components() {
    let distant = particle(vec![0.0, 1e300, 0.0]);
    let origin = particle(vec![0.0; 3]);
    let (energy, force) = morse_interaction(&distant, &origin, 1e100, 1e-300, 0.0).unwrap();
    let (_, opposite) = morse_interaction(&origin, &distant, 1e100, 1e-300, 0.0).unwrap();
    close(energy, 3.995_764_008_937_281e99);
    close(force[1], -4.650_883_158_696_593e-201);
    close(opposite[1], -force[1]);
    close(force[0], 0.0);
    close(force[2], 0.0);
}

#[test]
fn radial_rejects_invalid_arguments_or_truly_overflowing_results() {
    for (r, de, a, re) in [
        (-1.0, 1.0, 1.0, 1.0),
        (1.0, -1.0, 1.0, 1.0),
        (1.0, 1.0, 0.0, 1.0),
        (1.0, 1.0, 1.0, f64::NAN),
        (1.0, f64::INFINITY, 1.0, 1.0),
        (1.0, 1.0, 1e100, 2.0),
    ] {
        assert!(morse_radial(r, de, a, re).is_err());
    }
    let zero_strength = morse_radial(1.0, 0.0, 1e100, 2.0).unwrap();
    close(zero_strength.0, 0.0);
    close(zero_strength.1, 0.0);
    let long_tail = morse_radial(f64::MAX, 3.0, 2.0, 0.0).unwrap();
    close(long_tail.0, 3.0);
    close(long_tail.1, 0.0);
}
