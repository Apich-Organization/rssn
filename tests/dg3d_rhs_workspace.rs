//! Three-dimensional DG spatial workspace and face callback regressions.

use std::sync::Arc;
use std::sync::Mutex;
use std::sync::atomic::AtomicUsize;
use std::sync::atomic::Ordering;

use rssn::sim::dg::DgLaw;
use rssn::sim::dg3d::Dg3Bc;
use rssn::sim::dg3d::Dg3Mesh;
use rssn::sim::dg3d::Dg3Options;
use rssn::sim::dg3d::LinearAdvection3d;
use rssn::sim::dg3d::solve_dg3d;

fn mesh(
    ne: [usize; 3],
    p: usize,
) -> Dg3Mesh {
    Dg3Mesh {
        ne,
        lo: [-1.0, 0.0, -2.0],
        hi: [1.0, 2.0, 2.0],
        p,
    }
}

fn options(entropy_stable: bool) -> Dg3Options {
    Dg3Options {
        entropy_stable,
        ..Dg3Options::default()
    }
}

#[test]
fn affine_advection_preserves_fixed_face_stage_arguments_and_order() {
    let mesh = mesh([2, 1, 1], 2);
    let a = [0.75, -0.5, 0.25];
    let exact = |t: f64, x: [f64; 3]| 2.0 + 0.125 * x[0] - 0.25 * x[1] + 0.375 * x[2] - 0.3125 * t;
    let dt = 0.015625;
    for entropy_stable in [false, true] {
        let visited = Arc::new(Mutex::new(Vec::new()));
        let bcs: [Dg3Bc; 6] = std::array::from_fn(|face| {
            let visited = Arc::clone(&visited);
            Dg3Bc::Fixed(Arc::new(move |t, x| {
                visited.lock().unwrap().push((face, t, x));
                vec![exact(t, x)]
            }))
        });
        let s = solve_dg3d(
            &mesh,
            &LinearAdvection3d { a },
            &bcs,
            &|x| vec![exact(0.0, x)],
            dt,
            0.2,
            &options(entropy_stable),
        )
        .unwrap();
        assert_eq!(s.steps, 1);
        assert_eq!(s.t, dt);
        for e in 0..2 {
            for nd in 0..mesh.nodes_per_element() {
                assert!((s.value(e, nd, 0) - exact(dt, s.coord(e, nd))).abs() < 2e-14);
            }
        }

        let mut expected = Vec::new();
        for t in [0.0, dt, 0.5 * dt] {
            for e in 0..2 {
                for axis in 0..3 {
                    let transverse: Vec<_> = (0..3).filter(|&i| i != axis).collect();
                    for slow in [0.0, 0.5, 1.0] {
                        for fast in [0.0, 0.5, 1.0] {
                            for side in 0..2 {
                                if axis == 0 && side != e {
                                    continue;
                                }
                                let mut fraction = [0.0; 3];
                                fraction[axis] = side as f64;
                                fraction[transverse[0]] = fast;
                                fraction[transverse[1]] = slow;
                                let x = [
                                    -1.0 + e as f64 + fraction[0],
                                    2.0 * fraction[1],
                                    -2.0 + 4.0 * fraction[2],
                                ];
                                expected.push((2 * axis + side, t, x));
                            }
                        }
                    }
                }
            }
        }
        assert_eq!(*visited.lock().unwrap(), expected);
    }
}

#[test]
fn malformed_fixed_face_is_an_error_in_initial_and_later_rk_stages() {
    let mesh = mesh([1; 3], 1);
    let law = LinearAdvection3d { a: [1.0, 0.0, 0.0] };
    for (malformed_at_start, entropy_stable) in [(true, false), (false, true)] {
        let mut bcs = std::array::from_fn(|_| Dg3Bc::Outflow);
        bcs[0] = Dg3Bc::Fixed(Arc::new(move |t, _| {
            if malformed_at_start || t > 0.0 {
                vec![1.0, 2.0]
            } else {
                vec![1.0]
            }
        }));
        let error = solve_dg3d(
            &mesh,
            &law,
            &bcs,
            &|_| vec![1.0],
            0.01,
            0.1,
            &options(entropy_stable),
        )
        .unwrap_err();
        assert!(error.contains("fixed DG boundary state"));
    }
}

#[test]
fn a_malformed_face_stops_later_serial_callbacks_and_reports_parallel_errors() {
    let mesh = mesh([4, 2, 2], 1);
    let law = LinearAdvection3d { a: [1.0, 0.0, 0.0] };
    for threads in [1, 2] {
        let later = Arc::new(AtomicUsize::new(0));
        let mut bcs = std::array::from_fn(|_| Dg3Bc::Outflow);
        bcs[0] = Dg3Bc::Fixed(Arc::new(|_, _| Vec::new()));
        let later_calls = Arc::clone(&later);
        bcs[1] = Dg3Bc::Fixed(Arc::new(move |_, _| {
            later_calls.fetch_add(1, Ordering::Relaxed);
            vec![1.0]
        }));
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        let error = pool.install(|| {
            solve_dg3d(
                &mesh,
                &law,
                &bcs,
                &|_| vec![1.0],
                0.01,
                0.1,
                &options(false),
            )
            .unwrap_err()
        });
        assert!(error.contains("fixed DG boundary state"));
        if threads == 1 {
            assert_eq!(later.load(Ordering::Relaxed), 0);
        }
    }
}

#[test]
fn correctly_sized_nonfinite_fixed_values_reach_the_flux_unchanged() {
    let mesh = mesh([1; 3], 1);
    let law = LinearAdvection3d { a: [1.0, 0.0, 0.0] };
    for entropy_stable in [false, true] {
        let mut bcs = std::array::from_fn(|_| Dg3Bc::Outflow);
        bcs[0] = Dg3Bc::Fixed(Arc::new(|_, _| vec![f64::NAN]));
        let s = solve_dg3d(
            &mesh,
            &law,
            &bcs,
            &|_| vec![1.0],
            0.01,
            0.1,
            &options(entropy_stable),
        )
        .unwrap();
        assert_eq!(s.steps, 1);
        assert!(s.q.iter().any(|v| v.is_nan()));
    }
}

struct WallLaw {
    walls: AtomicUsize,
}

impl DgLaw for WallLaw {
    fn nvar(&self) -> usize {
        2
    }

    fn flux(
        &self,
        q: &[f64],
        dir: usize,
        out: &mut [f64],
    ) {
        if dir == 0 {
            out.copy_from_slice(q);
        } else {
            out.fill(0.0);
        }
    }

    fn max_speed(
        &self,
        _: &[f64],
        dir: usize,
    ) -> f64 {
        f64::from(u8::from(dir == 0))
    }

    fn numerical_flux(
        &self,
        left: &[f64],
        _: &[f64],
        dir: usize,
        out: &mut [f64],
    ) {
        assert_eq!(out, [0.0, 0.0]);
        if dir == 0 {
            out[0] = left[0];
        }
    }

    fn wall_state(
        &self,
        inside: &[f64],
        _: usize,
        outside: &mut [f64],
    ) {
        assert_eq!(outside, [0.0, 0.0]);
        outside[0] = inside[0];
        self.walls.fetch_add(1, Ordering::Relaxed);
    }
}

#[test]
fn reusable_surface_buffers_are_fresh_and_worker_private() {
    let mesh = mesh([4, 2, 2], 2);
    let bcs = std::array::from_fn(|_| Dg3Bc::Wall);
    let run = |threads| {
        let law = WallLaw {
            walls: AtomicUsize::new(0),
        };
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        let s = pool
            .install(|| {
                solve_dg3d(
                    &mesh,
                    &law,
                    &bcs,
                    &|x| vec![1.0 + 0.25 * x[0], 0.0],
                    0.025,
                    0.05,
                    &options(false),
                )
            })
            .unwrap();
        (
            s.q.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
            s.steps,
            law.walls.into_inner(),
        )
    };
    let one = run(1);
    let two = run(2);
    assert_eq!(one, two);
    let faces = 2 * (mesh.ne[0] * mesh.ne[1] + mesh.ne[1] * mesh.ne[2] + mesh.ne[0] * mesh.ne[2]);
    assert_eq!(one.2, 3 * one.1 * faces * (mesh.p + 1).pow(2));
}
