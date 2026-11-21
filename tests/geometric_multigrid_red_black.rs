//! Geometric multigrid cycles against independently assembled grid operators.

use std::f64::consts::PI;

use rssn::kernels::multigrid::Cycle;
use rssn::kernels::multigrid::GeometricMg;
use rssn::kernels::multigrid::MgOptions;
use rssn::kernels::multigrid::Smoother;

fn coordinates(
    dim: usize,
    n: usize,
) -> Vec<Vec<usize>> {
    (0..n.pow(dim as u32))
        .map(|mut index| {
            (0..dim)
                .map(|_| {
                    let coordinate = index % n;
                    index /= n;
                    coordinate
                })
                .collect()
        })
        .collect()
}

fn stencil(
    points: &[Vec<usize>],
    n: usize,
    shift: f64,
) -> Vec<Vec<f64>> {
    let inv_h2 = ((n + 1) as f64).powi(2);
    points
        .iter()
        .map(|row| {
            points
                .iter()
                .map(|column| {
                    match row
                        .iter()
                        .zip(column)
                        .map(|(a, b)| a.abs_diff(*b))
                        .sum::<usize>()
                    {
                        | 0 => 2.0 * row.len() as f64 * inv_h2 + shift,
                        | 1 => -inv_h2,
                        | _ => 0.0,
                    }
                })
                .collect()
        })
        .collect()
}

fn sweep(
    matrix: &[Vec<f64>],
    points: &[Vec<usize>],
    u: &mut [f64],
    rhs: &[f64],
    count: usize,
) {
    for _ in 0..count {
        for color in 0..2 {
            for (i, point) in points.iter().enumerate() {
                if point.iter().sum::<usize>() % 2 == color {
                    let row = &matrix[i];
                    let off_diagonal: f64 = row
                        .iter()
                        .zip(u.iter())
                        .enumerate()
                        .filter(|(j, _)| *j != i)
                        .map(|(_, (a, u))| a * u)
                        .sum();
                    u[i] = (rhs[i] - off_diagonal) / row[i];
                }
            }
        }
    }
}

fn dot(
    left: &[f64],
    right: &[f64],
) -> f64 {
    left.iter().zip(right).map(|(a, b)| a * b).sum()
}

fn reference_cycle(
    dim: usize,
    n: usize,
    shift: f64,
    u: &mut [f64],
    rhs: &[f64],
    options: &MgOptions,
) {
    let points = coordinates(dim, n);
    let matrix = stencil(&points, n, shift);
    if n == 1 {
        u[0] = rhs[0] / matrix[0][0];
        return;
    }
    sweep(&matrix, &points, u, rhs, options.pre);
    let residual: Vec<_> = rhs
        .iter()
        .zip(&matrix)
        .map(|(&rhs, row)| rhs - dot(row, u))
        .collect();
    let coarse_n = n / 2;
    let coarse = coordinates(dim, coarse_n);
    let restriction: Vec<Vec<f64>> = coarse
        .iter()
        .map(|center| {
            points
                .iter()
                .map(|point| {
                    center
                        .iter()
                        .zip(point)
                        .map(|(&c, &x)| {
                            match x.abs_diff(2 * c + 1) {
                                | 0 => 0.5,
                                | 1 => 0.25,
                                | _ => 0.0,
                            }
                        })
                        .product()
                })
                .collect()
        })
        .collect();
    let coarse_rhs: Vec<_> = restriction.iter().map(|row| dot(row, &residual)).collect();
    let mut correction = vec![0.0; coarse.len()];
    let visits = if options.cycle == Cycle::W {
        2
    } else {
        1
    };
    for _ in 0..visits {
        reference_cycle(dim, coarse_n, shift, &mut correction, &coarse_rhs, options);
    }
    // Interpolation is the transpose of full weighting times two per axis.
    let factor = 2.0_f64.powi(dim as i32);
    for (i, value) in u.iter_mut().enumerate() {
        let interpolated: f64 = restriction
            .iter()
            .zip(&correction)
            .map(|(row, &c)| factor * row[i] * c)
            .sum();
        *value += interpolated;
    }
    sweep(&matrix, &points, u, rhs, options.post);
}

#[test]
fn cycles_agree_with_dense_stencils_for_all_dimensions_and_colors() {
    for dim in 1..=3 {
        for n in [1, 3, 7] {
            let shift = 2.5;
            let grid = GeometricMg::new(dim, n, shift).unwrap();
            let initial: Vec<_> = (0..grid.size())
                .map(|i| ((29 * i % 23) as f64 - 11.0) / 32.0)
                .collect();
            let rhs: Vec<_> = (0..grid.size())
                .map(|i| (8 * i % 13) as f64 - 6.0)
                .collect();
            for cycle in [Cycle::V, Cycle::W] {
                for (pre, post) in [(0, 0), (0, 1), (1, 0), (2, 3)] {
                    let options = MgOptions {
                        cycle,
                        pre,
                        post,
                        smoother: Smoother::RedBlackGaussSeidel,
                        ..MgOptions::default()
                    };
                    let mut actual = initial.clone();
                    let mut expected = initial.clone();
                    grid.cycle_once(&mut actual, &rhs, &options);
                    reference_cycle(dim, n, shift, &mut expected, &rhs, &options);
                    for (index, (&a, &b)) in actual.iter().zip(&expected).enumerate() {
                        assert!(
                            (a - b).abs() < 5e-14,
                            "dim={dim}, n={n}, {cycle:?}, pre={pre}, post={post}, at={index}: {a} vs {b}"
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn solves_shifted_discrete_eigenmodes() {
    for (dim, n) in [(1, 127), (2, 31), (3, 15)] {
        let shift = 11.0;
        let grid = GeometricMg::new(dim, n, shift).unwrap();
        let h = 1.0 / (n + 1) as f64;
        let eigenvalue = shift
            + 4.0 / h.powi(2)
                * (1..=dim)
                    .map(|mode| (mode as f64 * PI * h * 0.5).sin().powi(2))
                    .sum::<f64>();
        let exact: Vec<f64> = coordinates(dim, n)
            .iter()
            .map(|point| {
                point
                    .iter()
                    .enumerate()
                    .map(|(axis, x)| (PI * h * ((axis + 1) * (x + 1)) as f64).sin())
                    .product()
            })
            .collect();
        let rhs: Vec<_> = exact.iter().map(|&u| eigenvalue * u).collect();
        let result = grid.solve(&rhs, None, &MgOptions::default());
        assert!(result.converged, "{dim}D residuals {:?}", result.history);
        for (actual, expected) in result.x.iter().zip(exact) {
            assert!((actual - expected).abs() < 2e-10);
        }
    }
}
