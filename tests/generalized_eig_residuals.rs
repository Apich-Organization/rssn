//! Generalized shift-invert residuals against analytical and exact pencils.

use std::cell::Cell;

use num_rational::BigRational;
use num_traits::ToPrimitive;
use rssn::kernels::eigs::EigsOptions;
use rssn::kernels::eigs::EigsResult;
use rssn::kernels::eigs_ext::GenEigsOptions;
use rssn::kernels::eigs_ext::eigs_gen_shift_invert;
use rssn::kernels::eigs_ext::eigs_gen_sym_shift_invert;

#[test]
fn real_generalized_eigenvectors_need_one_residual_matvec_each() {
    let diagonal = [2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0];
    let weight = [2.0, 1.0, 4.0, 2.0, 1.0, 4.0, 2.0, 1.0, 4.0];
    let options = GenEigsOptions {
        sigma: 3.0625,
        block: 2,
        eigs: EigsOptions::default(),
    };
    let a_calls = Cell::new(0);
    let result = eigs_gen_sym_shift_invert(
        |x, y| {
            a_calls.set(a_calls.get() + 1);
            assert!(x.iter().any(|&x| x != 0.0));
            for (i, value) in y.iter_mut().enumerate() {
                *value = (weight[i] * diagonal[i]) * x[i];
            }
        },
        |x, y| {
            for (i, value) in y.iter_mut().enumerate() {
                *value = weight[i] * x[i];
            }
        },
        |x, y| {
            for (i, value) in y.iter_mut().enumerate() {
                *value = x[i] / (weight[i] * (diagonal[i] - options.sigma));
            }
        },
        diagonal.len(),
        3,
        &options,
    )
    .unwrap();
    assert!(result.converged);
    assert_eq!(a_calls.get(), result.values.len());
    for (i, &expected) in [3.0, 3.25, 2.75].iter().enumerate() {
        let (eigenvalue, imaginary) = result.values[i];
        assert_eq!(imaginary, 0.0);
        assert!((eigenvalue - expected).abs() < 1e-12);
        let expected_residual = result.vectors_re[i]
            .iter()
            .enumerate()
            .map(|(j, &v)| weight[j] * v * (diagonal[j] - eigenvalue))
            .fold(0.0_f64, f64::hypot);
        assert!((result.residuals[i] - expected_residual).abs() < 1e-15);
        for j in 0..result.values.len() {
            let inner: f64 = result.vectors_re[i]
                .iter()
                .zip(&result.vectors_re[j])
                .zip(weight)
                .map(|((&x, &y), weight)| weight * x * y)
                .sum();
            assert!((inner - f64::from(u8::from(i == j))).abs() < 2e-14);
        }
    }
}

const BLOCKS: [(f64, f64); 7] = [
    (2.0, 1.0),
    (3.0, 2.0),
    (4.0, 3.0),
    (5.0, 1.0),
    (6.0, 3.0),
    (7.0, 4.0),
    (8.0, 2.0),
];

fn exact(value: f64) -> BigRational {
    BigRational::from_float(value).unwrap()
}

fn exact_complex_residual(
    solution: &EigsResult,
    index: usize,
) -> f64 {
    let (real, imag) = solution.values[index];
    let (real, imag) = (exact(real), exact(imag));
    let mut sum = exact(0.0);
    for (i, (a, b)) in BLOCKS.into_iter().enumerate() {
        let (a, b) = (exact(a), exact(b));
        let xr = exact(solution.vectors_re[index][2 * i]);
        let yr = exact(solution.vectors_re[index][2 * i + 1]);
        let xi = exact(solution.vectors_im[index][2 * i]);
        let yi = exact(solution.vectors_im[index][2 * i + 1]);
        let residuals = [
            (&a - &real) * &xr - &b * &yr + &imag * &xi,
            (&a - &real) * &xi - &b * &yi - &imag * &xr,
            (&a - &real) * &yr + &b * &xr + &imag * &yi,
            (&a - &real) * &yi + &b * &xi - &imag * &yr,
        ];
        // B has diagonal (1, 2) in every block.
        for (r, weight) in residuals.into_iter().zip([1.0, 1.0, 2.0, 2.0]) {
            let r = r * exact(weight);
            sum += &r * &r;
        }
    }
    sum.to_f64().unwrap().sqrt()
}

fn complex_pencil(
    scale: f64,
    ncv: usize,
) -> EigsResult {
    let a_calls = Cell::new(0);
    let options = GenEigsOptions {
        sigma: 2.0,
        block: 1,
        eigs: EigsOptions {
            ncv,
            max_restarts: 0,
            ..EigsOptions::default()
        },
    };
    let result = eigs_gen_shift_invert(
        |x, y| {
            a_calls.set(a_calls.get() + 1);
            for (i, (a, b)) in BLOCKS.into_iter().enumerate() {
                let (xr, yr) = (x[2 * i], x[2 * i + 1]);
                y[2 * i] = scale * (a * xr - b * yr);
                y[2 * i + 1] = (2.0 * scale) * (b * xr + a * yr);
            }
        },
        |x, y| {
            for (x, y) in x.chunks_exact(2).zip(y.chunks_exact_mut(2)) {
                y[0] = scale * x[0];
                y[1] = (2.0 * scale) * x[1];
            }
        },
        |x, y| {
            for (i, (a, b)) in BLOCKS.into_iter().enumerate() {
                let diagonal = a - options.sigma;
                let determinant = diagonal * diagonal + b * b;
                let (xr, yr) = (x[2 * i] / scale, x[2 * i + 1] / (2.0 * scale));
                y[2 * i] = (diagonal * xr + b * yr) / determinant;
                y[2 * i + 1] = (-b * xr + diagonal * yr) / determinant;
            }
        },
        2 * BLOCKS.len(),
        2,
        &options,
    )
    .unwrap();
    assert!(result.values.iter().all(|&(_, imag)| imag.abs() > 0.5));
    assert_eq!(a_calls.get(), 2 * result.values.len());
    result
}

#[test]
fn complex_generalized_residuals_remain_finite_at_both_exponent_extremes() {
    let full = complex_pencil(1.0, 2 * BLOCKS.len());
    assert!(full.converged);
    for ((real, imag), residual) in full.values.into_iter().zip(full.residuals) {
        assert!((real - 2.0).abs() < 1e-13);
        assert!((imag.abs() - 1.0).abs() < 1e-13);
        assert!(residual < 2e-13);
    }

    // Four Krylov vectors leave a nonzero residual from the other blocks,
    // allowing the same approximate eigenpairs to exercise the scaled norm.
    let ordinary = complex_pencil(1.0, 4);
    assert!(!ordinary.converged);
    for exp in [-600, 600] {
        let scale = 2.0_f64.powi(exp);
        let scaled = complex_pencil(scale, 4);
        assert_eq!(scaled.values, ordinary.values);
        assert_eq!(scaled.vectors_re, ordinary.vectors_re);
        assert_eq!(scaled.vectors_im, ordinary.vectors_im);
        for (i, &residual) in scaled.residuals.iter().enumerate() {
            let expected = exact_complex_residual(&scaled, i);
            assert!(expected > 1e-4);
            assert!(residual.is_finite() && residual > 0.0);
            assert!((residual / scale - expected).abs() < 4e-15);
            assert!((residual / scale - ordinary.residuals[i]).abs() < 4e-15);
        }
    }
}
