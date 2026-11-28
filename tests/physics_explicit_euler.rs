//! Endpoint, callback and arithmetic regressions for the fixed-step Euler family.

use rayon::ThreadPoolBuilder;
use rssn::sim::physics_em::ExplicitEulerMethod;
use rssn::sim::physics_em::solve_explicit_euler;
use rssn::sim::physics_em::solve_forward_euler;
use rssn::sim::physics_em::solve_heun_euler;
use rssn::sim::physics_em::solve_midpoint_euler;
use rssn::sim::physics_rkm::OdeSystem;
use std::sync::Mutex;

const METHODS: [ExplicitEulerMethod; 3] = [
    ExplicitEulerMethod::Forward,
    ExplicitEulerMethod::Midpoint,
    ExplicitEulerMethod::Heun,
];

struct Recorder {
    seen: Mutex<Vec<(f64, f64)>>,
    rhs: fn(f64, f64) -> f64,
}

impl Recorder {
    fn new(rhs: fn(f64, f64) -> f64) -> Self {
        Self {
            seen: Mutex::new(Vec::new()),
            rhs,
        }
    }
}

impl OdeSystem for Recorder {
    fn dim(&self) -> usize {
        1
    }

    fn eval(
        &self,
        t: f64,
        y: &[f64],
        out: &mut [f64],
    ) {
        self.seen.lock().unwrap().push((t, y[0]));
        out[0] = (self.rhs)(t, y[0]);
    }
}

#[test]
fn nominal_grid_clips_the_last_step_and_preserves_stage_order() {
    let dt: f64 = 0.3;
    let times = [0.0, dt, 2.0 * dt, dt.mul_add(3.0, 0.0), 1.0];
    for method in METHODS {
        let system = Recorder::new(|t, _| t);
        let history = solve_explicit_euler(&system, &[0.0], (0.0, 1.0), dt, method).unwrap();
        let reported: Vec<_> = history.iter().map(|row| row.0.to_bits()).collect();
        assert_eq!(reported, times.map(f64::to_bits));
        let expected = if method == ExplicitEulerMethod::Forward {
            0.36
        } else {
            0.5
        };
        assert!((history[4].1[0] - expected).abs() <= 3e-16);

        let evaluated: Vec<_> = system
            .seen
            .lock()
            .unwrap()
            .iter()
            .map(|row| row.0.to_bits())
            .collect();
        let mut expected_calls = Vec::new();
        for pair in times.windows(2) {
            expected_calls.push(pair[0].to_bits());
            match method {
                | ExplicitEulerMethod::Forward => (),
                | ExplicitEulerMethod::Midpoint => {
                    expected_calls.push(pair[0].midpoint(pair[1]).to_bits())
                },
                | ExplicitEulerMethod::Heun => expected_calls.push(pair[1].to_bits()),
            }
        }
        assert_eq!(evaluated, expected_calls);
    }
}

#[test]
fn aligned_decimal_steps_finish_without_a_duplicate_or_extra_callback() {
    let system = Recorder::new(|_, _| 1.0);
    let history =
        solve_explicit_euler(&system, &[1.0], (0.0, 1.0), 0.1, ExplicitEulerMethod::Heun).unwrap();
    assert_eq!(history.len(), 11);
    assert_eq!(system.seen.lock().unwrap().len(), 20);
    assert_eq!(history[10].0.to_bits(), 1.0f64.to_bits());
    assert!((history[10].1[0] - 2.0).abs() <= 2e-15);
    assert!(history.windows(2).all(|w| w[1].0 > w[0].0));
}

#[test]
fn legacy_history_keeps_complete_steps_and_accumulation_order() {
    let dt = 0.3;
    for solver in [
        solve_forward_euler::<Recorder>,
        solve_midpoint_euler::<Recorder>,
        solve_heun_euler::<Recorder>,
    ] {
        let system = Recorder::new(|_, _| 1.0);
        let history = solver(&system, &[1.0], (0.0, 1.0), dt);
        assert_eq!(history.len(), 5);
        let mut expected_t = 0.0f64;
        let mut expected_y = 1.0f64;
        for (t, y) in history {
            assert_eq!(t.to_bits(), expected_t.to_bits());
            assert_eq!(y[0].to_bits(), expected_y.to_bits());
            expected_t += dt;
            expected_y += dt;
        }
    }
}

#[test]
fn checked_backward_and_wide_finite_spans_are_supported() {
    for method in METHODS {
        let system = Recorder::new(|_, _| 1.0);
        let backward = solve_explicit_euler(&system, &[2.0], (1.0, 0.0), -0.3, method).unwrap();
        let (t, y) = backward.last().unwrap();
        assert_eq!(t.to_bits(), 0);
        assert!((y[0] - 1.0).abs() <= 1e-15);
        assert!(backward.windows(2).all(|w| w[1].0 < w[0].0));
        let system = Recorder::new(|_, _| 1e-308);
        let wide = solve_explicit_euler(&system, &[0.0], (-1e308, 1e308), 1e308, method).unwrap();
        assert_eq!(wide.len(), 3);
        assert_eq!(wide[2].0.to_bits(), 1e308f64.to_bits());
        assert!((wide[2].1[0] - 2.0).abs() <= 4e-16);
    }
}

#[test]
fn invalid_geometry_fails_before_any_rhs_call_and_zero_span_is_empty_work() {
    let system = Recorder::new(|_, _| 1.0);
    let method = ExplicitEulerMethod::Forward;
    let large = 9_007_199_254_740_992.0_f64;
    for (span, dt) in [
        ((0.0, 1.0), 0.0),
        ((0.0, 1.0), -0.5),
        ((0.0, 1.0), f64::INFINITY),
        ((f64::NAN, 1.0), 1.0),
        ((large, large.next_up()), 0.5),
        ((0.0, 1.0), f64::from_bits(1)),
    ] {
        assert!(solve_explicit_euler(&system, &[0.0], span, dt, method).is_err());
    }
    assert!(solve_explicit_euler(&system, &[0.0, 1.0], (0.0, 1.0), 0.25, method).is_err());
    assert!(solve_explicit_euler(&system, &[f64::NAN], (0.0, 1.0), 0.25, method).is_err());
    assert_eq!(
        solve_forward_euler(&system, &[0.0], (0.0, 1.0), 0.0).len(),
        1
    );
    assert_eq!(
        solve_forward_euler(&system, &[0.0], (large, large.next_up()), 0.5).len(),
        1
    );
    assert_eq!(
        solve_explicit_euler(&system, &[0.0], (0.5, 0.5), 0.25, method)
            .unwrap()
            .len(),
        1
    );
    assert!(system.seen.lock().unwrap().is_empty());
}

#[test]
fn checked_solver_stops_before_propagating_nonfinite_rhs() {
    for method in METHODS {
        let system = Recorder::new(|_, _| f64::NAN);
        assert!(solve_explicit_euler(&system, &[0.0], (0.0, 1.0), 0.3, method).is_err());
        assert_eq!(system.seen.lock().unwrap().len(), 1);
    }

    for method in [ExplicitEulerMethod::Midpoint, ExplicitEulerMethod::Heun] {
        let system = Recorder::new(|t, _| {
            if t == 0.0 {
                1.0
            } else {
                f64::INFINITY
            }
        });
        assert!(solve_explicit_euler(&system, &[0.0], (0.0, 1.0), 0.3, method).is_err());
        assert_eq!(system.seen.lock().unwrap().len(), 2);
    }
}

#[test]
fn checked_stage_arithmetic_retains_representable_extreme_updates() {
    let large = Recorder::new(|_, _| 1e308);
    let heun =
        solve_explicit_euler(&large, &[0.0], (0.0, 0.25), 0.25, ExplicitEulerMethod::Heun).unwrap();
    assert_eq!(heun[1].1[0].to_bits(), 2.5e307f64.to_bits());
    let forward = solve_explicit_euler(
        &large,
        &[-1e308],
        (0.0, 2.0),
        2.0,
        ExplicitEulerMethod::Forward,
    )
    .unwrap();
    assert_eq!(forward[1].1[0].to_bits(), 1e308f64.to_bits());

    for h_bits in [1, 3] {
        let tiny = f64::from_bits(h_bits);
        for method in METHODS {
            let system = Recorder::new(|_, _| 1e308);
            let trace = solve_explicit_euler(&system, &[0.0], (0.0, tiny), tiny, method).unwrap();
            let exact = tiny * 1e308;
            assert_eq!(trace[1].1[0].to_bits(), exact.to_bits());
            assert_eq!(trace[1].0.to_bits(), h_bits);
            if method == ExplicitEulerMethod::Midpoint {
                assert_eq!(
                    system.seen.lock().unwrap()[1].1.to_bits(),
                    (0.5 * exact).to_bits()
                );
            }
        }
    }
}

struct Coupled(usize);

impl OdeSystem for Coupled {
    fn dim(&self) -> usize {
        self.0
    }

    fn eval(
        &self,
        t: f64,
        y: &[f64],
        out: &mut [f64],
    ) {
        for (dest, &partner) in out.iter_mut().zip(y.iter().rev()) {
            *dest = 0.3 * t - 0.15 * partner;
        }
    }
}

#[test]
fn active_worker_count_does_not_change_numerical_bits() {
    let n = 65_536;
    let start: Vec<_> = (0..n).map(|i| f64::from(i % 17) / 8.0).collect();
    for method in METHODS {
        let run = |workers| {
            ThreadPoolBuilder::new()
                .num_threads(workers)
                .build()
                .unwrap()
                .install(|| {
                    solve_explicit_euler(&Coupled(start.len()), &start, (0.0, 1.0), 0.3, method)
                        .unwrap()
                })
        };
        let one = run(1);
        let two = run(2);
        for ((ta, a), (tb, b)) in one.iter().zip(&two) {
            assert_eq!(ta.to_bits(), tb.to_bits());
            assert!(a.iter().zip(b).all(|(&x, &y)| x.to_bits() == y.to_bits()));
        }
    }
}
