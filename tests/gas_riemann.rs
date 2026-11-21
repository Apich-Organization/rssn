//! Scale, wave-boundary and prepared sampling checks for the ideal-gas Riemann solver.

use rssn::sim::gas_dynamics::ExactRiemann;
use rssn::sim::gas_dynamics::Limiter;
use rssn::sim::gas_dynamics::Primitive;
use rssn::sim::gas_dynamics::exact_riemann;
use rssn::sim::gas_dynamics::sod_shock_tube;

const LEFT: Primitive = Primitive {
    rho: 1.0,
    u: 0.0,
    p: 1.0,
};
const RIGHT: Primitive = Primitive {
    rho: 0.125,
    u: 0.0,
    p: 0.1,
};
const GAMMA: f64 = 1.4;

fn near(
    actual: f64,
    expected: f64,
    relative_tolerance: f64,
) {
    assert!(
        (actual - expected).abs() <= relative_tolerance * expected.abs(),
        "{actual:.17e} vs {expected:.17e}"
    );
}

fn bits(w: Primitive) -> [u64; 3] {
    [w.rho.to_bits(), w.u.to_bits(), w.p.to_bits()]
}

#[test]
fn prepared_solver_matches_the_one_shot_call() {
    let prepared = ExactRiemann::new(LEFT, RIGHT, GAMMA);
    for k in -200..=200 {
        let xi = f64::from(k) / 70.0;
        assert_eq!(
            bits(prepared.sample(xi)),
            bits(exact_riemann(LEFT, RIGHT, GAMMA, xi))
        );
    }
}

#[test]
fn prepared_batch_matches_scalar_samples_in_order_under_an_active_pool() {
    let solution = ExactRiemann::new(LEFT, RIGHT, GAMMA);
    let queries = (0..10_000)
        .map(|i| f64::from((71 * i) % 1009) / 119.0 - 4.0)
        .collect::<Vec<_>>();
    let expected = queries
        .iter()
        .map(|&x| bits(solution.sample(x)))
        .collect::<Vec<_>>();
    for workers in [1, 2] {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(workers)
            .build()
            .unwrap();
        let values = pool.install(|| solution.sample_batch(&queries));
        assert_eq!(values.into_iter().map(bits).collect::<Vec<_>>(), expected);
    }
    assert!(solution.sample_batch(&[]).is_empty());
}

#[test]
fn sod_preserves_the_wave_solution_under_common_thermodynamic_scaling() {
    // Independent bisection of the ideal-gas Riemann pressure equation.
    let reference = [
        (-2.0, 1.0, 0.0, 1.0),
        (
            -0.5,
            0.602_937_696_498_180_7,
            0.569_346_630_516_602_9,
            0.492_471_851_553_222_5,
        ),
        (
            0.5,
            0.426_319_428_178_495_2,
            0.927_452_620_048_949_8,
            0.303_130_178_050_646_8,
        ),
        (
            1.2,
            0.265_573_711_705_307_1,
            0.927_452_620_048_949_8,
            0.303_130_178_050_646_8,
        ),
        (2.0, 0.125, 0.0, 0.1),
    ];
    for scale in [1e-300, 1e-200, 1e-10, 1.0, 1e150, 1e200, 1e300] {
        let adjust = |w: Primitive| {
            Primitive {
                rho: w.rho * scale,
                p: w.p * scale,
                ..w
            }
        };
        let prepared = ExactRiemann::new(adjust(LEFT), adjust(RIGHT), GAMMA);
        for (xi, rho, u, p) in reference {
            let w = prepared.sample(xi);
            near(w.rho / scale, rho, 3e-14);
            near(w.u, u, 3e-14);
            near(w.p / scale, p, 3e-14);
        }
    }
}

#[test]
fn exact_contact_chooses_the_left_trace() {
    let l = Primitive { u: 0.25, ..LEFT };
    let r = Primitive { rho: 0.25, ..l };
    let solved = ExactRiemann::new(l, r, GAMMA);
    assert_eq!(bits(solved.sample(l.u)), bits(l));
    assert_eq!(bits(solved.sample(l.u.next_up())), bits(r));
}

#[test]
fn symmetric_shock_uses_undisturbed_traces_at_the_fronts() {
    let l = Primitive { u: 1.0, ..LEFT };
    let r = Primitive { u: -1.0, ..LEFT };
    let solved = ExactRiemann::new(l, r, 3.0);
    // Rankine-Hugoniot gives speeds −2 and +2, star density 1.5 and p = 4.
    assert_eq!(bits(solved.sample(-2.0)), bits(l));
    assert_eq!(bits(solved.sample(2.0)), bits(r));
    for xi in [-1.75, 0.0, 1.75] {
        let w = solved.sample(xi);
        near(w.rho, 1.5, 4e-15);
        near(w.p, 4.0, 4e-15);
        near(w.u, 0.0, 4e-15);
    }
}

#[test]
fn vacuum_gap_has_zero_pressure_and_density_and_finite_fans() {
    let l = Primitive { u: -10.0, ..LEFT };
    let r = Primitive { u: 10.0, ..LEFT };
    let c = GAMMA.sqrt();
    let tail_distance = 2.0 * c / (GAMMA - 1.0);
    let l_tail = l.u + tail_distance;
    let r_tail = r.u - tail_distance;
    let solved = ExactRiemann::new(l, r, GAMMA);
    assert_eq!(bits(solved.sample(l.u - c)), bits(l));
    assert_eq!(bits(solved.sample(r.u + c)), bits(r));
    for xi in [l_tail, -2.0, 0.0, 2.0, r_tail] {
        assert_eq!(bits(solved.sample(xi)), [0, xi.to_bits(), 0]);
    }
    let near_tail = solved.sample(l_tail.next_down());
    assert!(near_tail.rho.is_finite() && near_tail.p.is_finite());
    let xi = -6.0;
    let fan_velocity = (2.0 * (c + xi) + (GAMMA - 1.0) * l.u) / (GAMMA + 1.0);
    let relative_sound = (l_tail - xi) / ((GAMMA + 1.0) / (GAMMA - 1.0) * c);
    near(solved.sample(xi).u, fan_velocity, 4e-15);
    near(
        solved.sample(xi).p,
        relative_sound.powf(2.0 * GAMMA / (GAMMA - 1.0)),
        4e-14,
    );
    let mirrored = solved.sample(-xi);
    near(mirrored.p, solved.sample(xi).p, 4e-15);
    near(mirrored.u, -fan_velocity, 4e-15);

    let l = Primitive {
        u: -tail_distance,
        ..LEFT
    };
    let r = Primitive {
        u: tail_distance,
        ..LEFT
    };
    let touching_fans = ExactRiemann::new(l, r, GAMMA);
    assert_eq!(bits(touching_fans.sample(0.0)), [0; 3]);
}

#[test]
fn approaching_vacuum_from_positive_pressure_does_not_hit_an_absolute_floor() {
    let v = 5.91;
    let l = Primitive { u: -v, ..LEFT };
    let r = Primitive { u: v, ..LEFT };
    let c = GAMMA.sqrt();
    let ratio = (c - 0.5 * (GAMMA - 1.0) * v) / c;
    let w = ExactRiemann::new(l, r, GAMMA).sample(0.0);
    assert!(w.rho > 0.0 && w.p > 0.0 && w.p < 1e-12);
    near(w.rho, ratio.powf(2.0 / (GAMMA - 1.0)), 5e-14);
    near(w.p, ratio.powf(2.0 * GAMMA / (GAMMA - 1.0)), 5e-14);
    near(w.u, 0.0, 4e-15);
}

#[test]
fn initial_shock_tube_densities_agree_at_the_central_discontinuity() {
    let (x, numeric, exact) = sod_shock_tube(9, 0.0, Limiter::VanLeer);
    assert_eq!(x[4].to_bits(), 0.5f64.to_bits());
    assert_eq!(exact[4].to_bits(), RIGHT.rho.to_bits());
    assert_eq!(
        numeric.into_iter().map(f64::to_bits).collect::<Vec<_>>(),
        exact.into_iter().map(f64::to_bits).collect::<Vec<_>>()
    );
}
