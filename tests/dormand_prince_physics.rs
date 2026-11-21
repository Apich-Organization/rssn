//! Scope checks for the Dormand-Prince physics driver and its callback contract.

use rssn::sim::physics_rkm::DormandPrince54;
use rssn::sim::physics_rkm::OdeSystem;
use std::sync::Mutex;

struct Traced<F> {
    times: Mutex<Vec<f64>>,
    slope: F,
}

impl<F> Traced<F> {
    fn new(slope: F) -> Self {
        Self {
            times: Mutex::new(Vec::new()),
            slope,
        }
    }

    fn times(&self) -> Vec<f64> {
        self.times.lock().unwrap().clone()
    }
}

impl<F: Fn(f64, usize) -> f64 + Send + Sync> OdeSystem for Traced<F> {
    fn dim(&self) -> usize {
        1
    }

    fn eval(
        &self,
        t: f64,
        y: &[f64],
        out: &mut [f64],
    ) {
        assert!(t.is_finite() && y.iter().all(|v| v.is_finite()));
        let mut times = self.times.lock().unwrap();
        out[0] = (self.slope)(t, times.len());
        times.push(t);
    }
}

#[test]
fn short_and_subnormal_spans_reach_the_endpoint() {
    let dp = DormandPrince54::new();
    let fast = Traced::new(|_, _| 1e300);
    let rows = dp.solve(&fast, &[0.0], (0.0, 1e-300), 1e-301, (1e-10, 1e-12));
    let (t, y) = rows.last().unwrap();
    assert_eq!(t.to_bits(), 1e-300_f64.to_bits());
    assert!((y[0] - 1.0).abs() < 1e-14);
    assert!(fast.times().len() >= 7);

    let tiny = f64::from_bits(1);
    let rate = f64::from_bits(0x7fe0_0000_0000_0000);
    let fastest = Traced::new(move |_, _| rate);
    let rows = dp.solve(&fastest, &[0.0], (0.0, tiny), tiny, (1e-8, 1e-24));
    let (t, y) = rows.last().unwrap();
    assert_eq!(t.to_bits(), tiny.to_bits());
    assert!((y[0] - tiny * rate).abs() < 1e-28);
    assert_eq!(fastest.times().len(), 7);
}

#[test]
fn clock_checks_do_no_rhs_work_and_clipped_endpoint_stages_share_the_exact_time() {
    let dp = DormandPrince54::new();
    let traced = Traced::new(|_, _| 1.0);
    let t = 2.0_f64.powi(54);
    let rows = dp.solve(&traced, &[0.0], (t, t.next_up()), 1.0, (1e-8, 1e-9));
    assert_eq!(rows.len(), 1);
    assert!(traced.times().is_empty());

    for (span, step, tol, init) in [
        ((1.0, 0.0), 0.1, (1e-8, 1e-9), vec![0.0]),
        ((0.0, f64::INFINITY), 0.1, (1e-8, 1e-9), vec![0.0]),
        ((0.0, 1.0), 0.0, (1e-8, 1e-9), vec![0.0]),
        ((0.0, 1.0), 0.1, (1e-8, f64::NAN), vec![0.0]),
        ((0.0, 1.0), 0.1, (0.0, 0.0), vec![0.0]),
        ((0.0, 1.0), 0.1, (1e-8, 1e-9), vec![0.0, 0.0]),
        ((0.0, 1.0), 0.1, (1e-8, 1e-9), vec![f64::NAN]),
    ] {
        assert_eq!(dp.solve(&traced, &init, span, step, tol).len(), 1);
    }
    assert!(traced.times().is_empty());

    let end = 2.0_f64.powi(-54);
    let rows = dp.solve(&traced, &[0.0], (-1.0, end), 2.0, (1e-8, 1e-9));
    assert_eq!(rows.len(), 2);
    assert_eq!(rows[1].0.to_bits(), end.to_bits());
    assert!((rows[1].1[0] - 1.0).abs() < 1e-14);
    let times = traced.times();
    assert_eq!(times.len(), 7);
    assert_eq!(times[5].to_bits(), end.to_bits());
    assert_eq!(times[6].to_bits(), end.to_bits());
    assert!(times.into_iter().all(|t| (-1.0..=end).contains(&t)));
}

#[test]
fn an_overflowing_error_square_rejects_and_recovers() {
    // The isolated nonzero sample has zero integral. The first large step
    // samples it at c=0.3. The smaller trials do not see the point.
    let spike = Traced::new(|t: f64, _| {
        if t.to_bits() == 0.3_f64.to_bits() {
            1e100
        } else {
            0.0
        }
    });
    let rows = DormandPrince54::new().solve(&spike, &[2.0], (0.0, 1.0), 1.0, (0.0, 1e-100));
    let times = spike.times();
    assert!(times.len() >= 14);
    assert_eq!(times[2].to_bits(), 0.3_f64.to_bits());
    assert_eq!(times[7].to_bits(), 0.0_f64.to_bits());
    assert_eq!(rows[1].0.to_bits(), 0.1_f64.to_bits());
    assert_eq!(rows.last().unwrap().0.to_bits(), 1.0_f64.to_bits());
    assert!(
        rows.iter()
            .all(|(_, y)| y[0].to_bits() == 2.0_f64.to_bits())
    );
}

#[test]
fn bad_oracle_outputs_preserve_the_last_accepted_state_without_retry() {
    let first = Traced::new(|_, _| f64::NAN);
    let result = DormandPrince54::new().solve(&first, &[1.0], (0.0, 1.0), 0.125, (1e-8, 1e-9));
    assert_eq!(result.len(), 1);
    assert_eq!(first.times().len(), 1);

    let next = Traced::new(|_, i| if i == 7 { f64::INFINITY } else { 1.0 });
    let rows = DormandPrince54::new().solve(&next, &[0.0], (0.0, 1.0), 0.125, (1e-8, 1e-9));
    assert_eq!(rows.len(), 2);
    assert_eq!(rows[1].0.to_bits(), 0.125_f64.to_bits());
    assert!((rows[1].1[0] - 0.125).abs() < 1e-16);
    assert_eq!(next.times().len(), 8);
}

#[test]
fn relative_only_zero_state_has_zero_defect() {
    let zero = Traced::new(|_, _| 0.0);
    let rows = DormandPrince54::new().solve(&zero, &[0.0], (0.0, 1.0), 0.0625, (1e-8, 0.0));
    assert_eq!(rows.last().unwrap().0.to_bits(), 1.0_f64.to_bits());
    assert!(
        rows.iter()
            .all(|(_, y)| y[0].to_bits() == 0.0_f64.to_bits())
    );
}

struct LinearSystem(usize);

impl OdeSystem for LinearSystem {
    fn dim(&self) -> usize {
        self.0
    }

    fn eval(
        &self,
        _t: f64,
        y: &[f64],
        out: &mut [f64],
    ) {
        for (v, &y) in out.iter_mut().zip(y) {
            *v = -0.5 * y;
        }
    }
}

#[test]
fn fused_stage_rows_give_the_same_bits_on_active_worker_pools() {
    for dim in [4096, 16_385] {
        let y = vec![2.0; dim];
        let solve = |workers| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(workers)
                .build()
                .unwrap()
                .install(|| {
                    DormandPrince54::new().solve(
                        &LinearSystem(y.len()),
                        &y,
                        (0.0, 2.0),
                        0.8,
                        (1e-9, 1e-10),
                    )
                })
        };
        let a = solve(1);
        let b = solve(2);
        assert_eq!(a.len(), b.len());
        for ((t1, y1), (t2, y2)) in a.iter().zip(&b) {
            assert_eq!(t1.to_bits(), t2.to_bits());
            assert!(
                y1.iter()
                    .zip(y2)
                    .all(|(v1, v2)| v1.to_bits() == v2.to_bits())
            );
        }
        let (t, y) = a.last().unwrap();
        assert_eq!(t.to_bits(), 2.0_f64.to_bits());
        assert!((y[0] - 2.0 * (-1.0_f64).exp()).abs() < 1e-9);
    }
}
