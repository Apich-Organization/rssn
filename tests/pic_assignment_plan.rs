//! Independent particle-mesh transfer and conservation checks.

use rssn::sim::pic::Grid;
use rssn::sim::pic::Shape;

const SHAPES: [Shape; 2] = [Shape::Linear, Shape::Cubic];

fn independent_stencil(
    shape: Shape,
    grid: &Grid,
    x: f64,
) -> Vec<(usize, f64)> {
    let nodes: &[i32] = match shape {
        | Shape::Linear => &[0, 1],
        | Shape::Cubic => &[-1, 0, 1, 2],
    };
    let u = (x.rem_euclid(grid.length()) / grid.length()) * grid.cells() as f64;
    let cell = u.floor() as i32;
    let fraction = u - f64::from(cell);
    nodes
        .iter()
        .map(|&a| {
            let w = nodes
                .iter()
                .filter(|&&b| b != a)
                .map(|&b| (fraction - f64::from(b)) / f64::from(a - b))
                .product();
            ((cell + a).rem_euclid(grid.cells() as i32) as usize, w)
        })
        .collect()
}

fn assert_bits(
    x: &[f64],
    y: &[f64],
) {
    assert_eq!(x.len(), y.len());
    for (i, (&x, &y)) in x.iter().zip(y).enumerate() {
        assert_eq!(x.to_bits(), y.to_bits(), "index {i}");
    }
}

#[test]
fn direct_and_prepared_transfers_match_lagrange_references_at_periodic_seams() {
    let grid = Grid::new(16, 2.0).unwrap();
    let winding = (1_u64 << 40) as f64;
    let pos = [
        -0.03125,
        -2.0,
        -0.0,
        0.09375,
        0.625,
        1.96875,
        winding + 0.0625,
        -winding - 0.0625,
    ];
    let q = [0.4, -0.2, 0.6, -0.25, 1.0, -0.8, 0.75, -1.25];
    let field: Vec<f64> = (0..16).map(|i| (f64::from(i) * 0.37).sin()).collect();
    for shape in SHAPES {
        let plan = grid.prepare_assignment(&pos, shape).unwrap();
        let got = grid.deposit(&pos, &q, shape);
        let gathered = grid.gather(&field, &pos, shape).unwrap();
        assert_bits(&plan.deposit(&q).unwrap(), &got);
        assert_bits(&plan.gather(&field).unwrap(), &gathered);
        let mut want = [0.0; 16];
        for (i, &x) in pos.iter().enumerate() {
            let mut expected_gather = 0.0;
            for (cell, weight) in independent_stencil(shape, &grid, x) {
                want[cell] += q[i] * weight;
                expected_gather += field[cell] * weight;
            }
            assert!((gathered[i] - expected_gather).abs() < 1e-15);
        }
        for (&actual, expected) in got.iter().zip(want) {
            assert!((actual * grid.spacing() - expected).abs() < 1e-15);
        }
        assert!((got.iter().sum::<f64>() * grid.spacing() - q.iter().sum::<f64>()).abs() < 2e-15);
    }
}

#[test]
fn deposit_and_gather_are_adjoint_for_reused_particle_positions() {
    let grid = Grid::new(64, 3.7).unwrap();
    let pos: Vec<f64> = (0..101)
        .map(|i| f64::from(i * 73 % 211) / 23.0 - 4.0)
        .collect();
    for shape in SHAPES {
        let plan = grid.prepare_assignment(&pos, shape).unwrap();
        for sequence in [7, 11, 15] {
            let charges: Vec<f64> = (0..101)
                .map(|i| f64::from(i * sequence % 37) / 16.0 - 1.0)
                .collect();
            let field: Vec<f64> = (0..64)
                .map(|i| (f64::from(i * sequence) / 11.0).sin())
                .collect();
            let rho = plan.deposit(&charges).unwrap();
            let gather = plan.gather(&field).unwrap();
            let mesh: f64 = rho
                .iter()
                .zip(&field)
                .map(|(r, f)| r * f * grid.spacing())
                .sum();
            let particles: f64 = charges.iter().zip(&gather).map(|(q, f)| q * f).sum();
            assert!((mesh - particles).abs() < 2e-14);
        }
    }
}

#[test]
fn subnormal_mesh_assigns_the_last_particle_before_the_periodic_edge() {
    let tiny = f64::from_bits(1);
    assert!(Grid::new(8, tiny).is_err());
    let grid = Grid::new(8, f64::from_bits(9)).unwrap();
    assert_eq!(grid.spacing(), tiny);
    let pos = [f64::from_bits(8)];
    let q = [tiny];
    let rho = grid.deposit(&pos, &q, Shape::Linear);
    assert!((rho[7] - 64.0 / 81.0).abs() < 3e-15);
    assert!((rho[0] - 8.0 / 81.0).abs() < 3e-15);
    assert_eq!(&rho[1..7], &[0.0; 6]);

    let field = [2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0];
    for shape in SHAPES {
        let plan = grid.prepare_assignment(&pos, shape).unwrap();
        let expected: f64 = independent_stencil(shape, &grid, pos[0])
            .into_iter()
            .map(|(i, w)| w * field[i])
            .sum();
        assert!((plan.gather(&field).unwrap()[0] - expected).abs() < 2e-15);
        let rho = plan.deposit(&q).unwrap();
        assert!((rho.iter().sum::<f64>() * (9.0 / 8.0) - 1.0).abs() < 4e-15);
    }
}

#[test]
fn repeated_and_parallel_gathers_keep_particle_order_and_bits() {
    let grid = Grid::new(32, 1.0).unwrap();
    let pos: Vec<f64> = (0..32_771)
        .map(|i| f64::from(i * 97 % 109_817) / 1024.0 - 10.0)
        .collect();
    let field: Vec<f64> = (0..32).map(|i| (f64::from(i) * 0.9).sin()).collect();
    for shape in SHAPES {
        let run = |workers| {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(workers)
                .build()
                .unwrap();
            pool.install(|| {
                let plan = grid.prepare_assignment(&pos, shape).unwrap();
                (
                    grid.gather(&field, &pos, shape).unwrap(),
                    plan.gather(&field).unwrap(),
                )
            })
        };
        let (a, b) = run(1);
        let (c, d) = run(4);
        assert_bits(&a, &b);
        assert_bits(&a, &c);
        assert_bits(&a, &d);
    }
}

#[test]
fn prepared_shapes_validate_counts_without_changing_legacy_zip_or_signed_zero() {
    let grid = Grid::new(8, 1.0).unwrap();
    let shape = Shape::Linear;
    let plan = grid.prepare_assignment(&[0.05], shape).unwrap();
    assert!(plan.deposit(&[]).is_err());
    assert!(plan.gather(&[1.0; 7]).is_err());
    assert!(grid.prepare_assignment(&[f64::NAN], shape).is_err());
    assert_bits(
        &grid.deposit(&[0.05, 0.7], &[2.0], shape),
        &plan.deposit(&[2.0]).unwrap(),
    );
    assert_bits(
        &grid.deposit(&[0.05], &[2.0, 3.0], shape),
        &plan.deposit(&[2.0]).unwrap(),
    );
    let field = [-0.0; 8];
    assert_eq!(
        plan.gather(&field).unwrap()[0].to_bits(),
        (-0.0_f64).to_bits()
    );
    assert_bits(
        &grid.gather(&field, &[0.05], shape).unwrap(),
        &plan.gather(&field).unwrap(),
    );
}
