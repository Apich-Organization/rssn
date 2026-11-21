//! Boundary and parallel execution checks for three-dimensional finite-volume advection.

use rssn::sim::physics_fvm::Mesh3D;
use rssn::sim::physics_fvm::solve_advection_3d;

#[test]
fn advection_3d_unmarked_faces_translate_with_zero_inflow() {
    let shape = [5, 3, 4];
    let [w, h, d] = shape;
    let mut mesh = Mesh3D::new(w, h, d, (5.0, 3.0, 4.0), |x, y, z| x + 10.0 * y + 100.0 * z);
    let before: Vec<f64> = mesh.cells.iter().map(|cell| cell.value).collect();
    let strides = [1, w, w * h];
    for (axis, stride) in strides.into_iter().enumerate() {
        for speed in [-1.0, 1.0] {
            let mut velocity = [0.0; 3];
            velocity[axis] = speed;
            let [vx, vy, vz] = velocity;
            let result =
                solve_advection_3d(&mut mesh, (vx, vy, vz), 1.0, 2, |_, _, _, _, _, _| false);
            for (idx, actual) in result.into_iter().enumerate() {
                let coordinate = idx / stride % shape[axis];
                let expected = if speed > 0.0 && coordinate >= 2 {
                    before[idx - 2 * stride]
                } else if speed < 0.0 && coordinate + 2 < shape[axis] {
                    before[idx + 2 * stride]
                } else {
                    0.0
                };
                assert_eq!(actual, expected, "axis {axis}, speed {speed}, cell {idx}");
            }
        }
    }
    assert_eq!(
        mesh.cells.iter().map(|cell| cell.value).collect::<Vec<_>>(),
        before
    );
}

#[test]
fn advection_3d_masked_cells_take_effect_at_the_end_of_each_step() {
    use std::sync::atomic::AtomicUsize;
    use std::sync::atomic::Ordering;

    let (w, h, d) = (67, 41, 8);
    let mut mesh = Mesh3D::new(w, h, d, (67.0, 41.0, 8.0), |_, _, _| 3.0);
    let callbacks: Vec<AtomicUsize> = (0..mesh.cells.len()).map(|_| AtomicUsize::new(0)).collect();
    let result = solve_advection_3d(
        &mut mesh,
        (1.0, 0.0, 0.0),
        1.0,
        2,
        |x, y, z, width, height, depth| {
            assert_eq!((width, height, depth), (w, h, d));
            callbacks[x + w * (y + h * z)].fetch_add(1, Ordering::Relaxed);
            x == 6
        },
    );
    assert!(
        callbacks
            .iter()
            .all(|calls| calls.load(Ordering::Relaxed) == 2)
    );
    for row in result.chunks_exact(w) {
        for (x, &value) in row.iter().enumerate() {
            let expected = if x < 2 || x == 6 || x == 7 {
                0.0
            } else {
                3.0
            };
            assert_eq!(value, expected);
        }
    }
}

#[test]
fn advection_3d_handles_empty_and_single_cell_meshes() {
    let mut empty = Mesh3D::new(3, 0, 3, (3.0, 1.0, 3.0), |_, _, _| unreachable!());
    assert!(
        solve_advection_3d(
            &mut empty,
            (1.0, 1.0, 1.0),
            0.25,
            2,
            |_, _, _, _, _, _| unreachable!()
        )
        .is_empty()
    );
    let mut mesh = Mesh3D::new(1, 1, 1, (1.0, 1.0, 1.0), |_, _, _| 8.0);
    let velocity = (0.25, -0.25, 0.25);
    assert_eq!(
        solve_advection_3d(
            &mut mesh,
            velocity,
            1.0,
            0,
            |_, _, _, _, _, _| unreachable!()
        ),
        [8.0]
    );
    assert_eq!(
        solve_advection_3d(&mut mesh, velocity, 1.0, 2, |_, _, _, _, _, _| false),
        [0.5]
    );
}

#[test]
fn advection_3d_is_bit_identical_across_worker_counts() {
    let solve = || {
        let mut mesh = Mesh3D::new(37, 17, 23, (1.0, 1.5, 2.0), |x, y, z| (x - y + z).sin());
        solve_advection_3d(
            &mut mesh,
            (-0.25, 0.125, -0.5),
            0.0125,
            7,
            |x, y, z, _, _, _| (x + y + z) % 29 == 0,
        )
        .into_iter()
        .map(f64::to_bits)
        .collect::<Vec<_>>()
    };
    let one = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .unwrap();
    let four = rayon::ThreadPoolBuilder::new()
        .num_threads(4)
        .build()
        .unwrap();
    assert_eq!(one.install(solve), four.install(solve));
}
