use std::cell::RefCell;

use rssn::kernels::ode_adaptive::Event;
use rssn::kernels::ode_adaptive::OdeError;
use rssn::kernels::ode_adaptive::OdeOptions;
use rssn::kernels::ode_adaptive::OdeSolution;
use rssn::kernels::ode_adaptive::dopri5;
use rssn::kernels::ode_adaptive::rosenbrock23_events;

type Rhs<'a> = &'a dyn Fn(f64, &[f64], &mut [f64]);
type Solve =
    fn(Rhs<'_>, f64, f64, &[f64], &OdeOptions, &[Event<'_>]) -> Result<OdeSolution, OdeError>;

fn solvers() -> [(&'static str, Solve); 2] {
    [
        ("DOPRI5", |f, a, b, y, opts, events| {
            dopri5(f, a, b, y, opts, events)
        }),
        ("Rosenbrock23", |f, a, b, y, opts, events| {
            rosenbrock23_events(f, a, b, y, opts, events)
        }),
    ]
}

#[test]
fn curvature_heuristic_preserves_progress_at_large_times_in_both_directions() {
    let time = 2.0_f64.powi(50);
    for (name, solve) in solvers() {
        for (start, end) in [
            (time, time + 16.0),
            (time + 16.0, time),
            (-time, -time + 16.0),
            (-time + 16.0, -time),
        ] {
            let times = RefCell::new(Vec::new());
            let rhs = |t: f64, _: &[f64], dy: &mut [f64]| {
                assert!((start.min(end)..=start.max(end)).contains(&t));
                times.borrow_mut().push(t);
                dy.fill(0.0);
            };
            let sol = solve(&rhs, start, end, &[3.0], &OdeOptions::default(), &[])
                .unwrap_or_else(|error| panic!("{name}: {error:?} at {start},{end}"));
            assert_eq!(sol.t.last().unwrap().to_bits(), end.to_bits());
            assert_eq!(sol.nfev, times.borrow().len());
            assert!(sol.y.iter().all(|y| y == &[3.0]));
            let advance = if end > start {
                start.next_up()
            } else {
                start.next_down()
            };
            assert!((sol.t[1] - start).abs() >= (advance - start).abs());
        }
    }
}

#[test]
fn large_time_startup_still_reaches_a_terminal_event_on_the_local_lattice() {
    let origin = 2.0_f64.powi(50);
    let target = origin + 7.0;
    let event = Event {
        g: &|t, _| t - target,
        terminal: true,
        direction: 0,
    };
    for (start, end) in [(origin, origin + 16.0), (origin + 16.0, origin)] {
        let sol = dopri5(
            |_, _, out| out.fill(0.0),
            start,
            end,
            &[3.0],
            &OdeOptions::default(),
            std::slice::from_ref(&event),
        )
        .unwrap();
        assert_eq!(sol.events.len(), 1);
        let location = sol.events[0].t;
        assert!((location - target).abs() <= target.next_up() - target);
        assert_eq!(sol.t.last().unwrap().to_bits(), location.to_bits());
        assert_eq!(sol.events[0].y, [3.0]);
    }
}

#[test]
fn automatic_probe_stays_within_a_maximum_between_time_lattice_points() {
    let start = 2.0_f64.powi(50);
    let times = RefCell::new(Vec::new());
    let rhs = |t: f64, _: &[f64], dy: &mut [f64]| {
        times.borrow_mut().push(t);
        dy.fill(0.0);
    };
    let tighter = OdeOptions {
        hmax: 0.24,
        ..Default::default()
    };
    assert!(matches!(
        dopri5(rhs, start, start + 1.0, &[1.0], &tighter, &[]),
        Err(OdeError::StepLimit)
    ));
    assert_eq!(*times.borrow(), [start]);
    times.borrow_mut().clear();

    let rounded = OdeOptions {
        hmax: 0.4,
        ..Default::default()
    };
    let sol = dopri5(rhs, start, start + 1.0, &[1.0], &rounded, &[]).unwrap();
    assert_eq!(times.borrow()[1].to_bits(), start.next_up().to_bits());
    assert_eq!(sol.nfev, times.borrow().len());
    assert_eq!(sol.t.last().unwrap().to_bits(), (start + 1.0).to_bits());
}

#[test]
fn startup_uses_the_actual_inward_neighbor_at_the_finite_time_limits() {
    for (start, end) in [
        (f64::MAX, f64::MAX.next_down()),
        (-f64::MAX, (-f64::MAX).next_up()),
        (0.0, f64::from_bits(1)),
    ] {
        let times = RefCell::new(Vec::new());
        let solution = dopri5(
            |t, _, dy| {
                times.borrow_mut().push(t);
                dy.fill(0.0);
            },
            start,
            end,
            &[3.0],
            &OdeOptions::default(),
            &[],
        )
        .unwrap();
        assert_eq!(solution.t, [start, end]);
        assert_eq!(solution.y, [[3.0], [3.0]]);
        assert_eq!(times.borrow()[1].to_bits(), end.to_bits());
        assert_eq!(solution.nfev, times.borrow().len());
        assert!(
            times
                .borrow()
                .iter()
                .all(|t| (start.min(end)..=start.max(end)).contains(t))
        );
    }
}
