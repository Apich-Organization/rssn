// Dimensionless polynomial coordinates and coefficients in ascending order.
use crate::kernels::complex_roots as polynomial_roots;
use crate::kernels::reuse::Scratch as PolynomialScratch;
use crate::kernels::reuse::global_pool as polynomial_reuse_pool;
use num_complex::Complex64 as PolynomialComplex;

struct ComplexPolynomialRootSet {
    roots: Vec<PolynomialComplex>,
    backward_errors: Vec<f64>,
    coefficient_error: f64,
    degree: usize,
}
record!(ComplexPolynomialRootSet {
    roots: Vec<PolynomialComplex>,
    backward_errors: Vec<f64>,
    coefficient_error: f64,
    degree: usize,
});

struct ComplexPolynomialRoot {
    root: PolynomialComplex,
    backward_error: f64,
}
record!(ComplexPolynomialRoot {
    root: PolynomialComplex,
    backward_error: f64,
});

fn polynomial_problem(
    coefficients: &[PolynomialComplex],
    iterations: usize,
    tolerance: f64,
    simultaneous: bool,
) -> Result<PolynomialScratch<'static, PolynomialComplex>, OpError> {
    const MAX_COEFFICIENTS: usize = 1_025;
    const MAX_WORK: usize = 4_194_304;
    crate::cancellation::checkpoint()?;
    if coefficients.len() > MAX_COEFFICIENTS {
        return Err(OpError::BudgetExceeded);
    }
    if iterations == 0 || !tolerance.is_finite() || !(0.0 < tolerance && tolerance < 1.0) {
        return Err(OpError::Domain(
            "expected a positive iteration count and a relative tolerance in (0, 1)".into(),
        ));
    }
    let normalized =
        polynomial_roots::normalized_coefficients(coefficients).map_err(OpError::Domain)?;
    let degree = normalized.len().saturating_sub(1).max(1);
    let cost = if simultaneous {
        degree.saturating_mul(degree)
    } else {
        degree
    };
    if cost.saturating_mul(iterations) > MAX_WORK {
        return Err(OpError::BudgetExceeded);
    }
    Ok(normalized)
}

// Multiply the linear factors in place, starting with the normalized leading coefficient.
// Checking the polynomial reconstructed from all roots detects duplicate simple roots that
// would each pass a pointwise evaluation at the same root.
fn coefficient_error(
    coefficients: &[PolynomialComplex],
    roots: &[PolynomialComplex],
) -> f64 {
    let mut reconstructed = polynomial_reuse_pool().scratch_complex64(coefficients.len());
    reconstructed.extend(coefficients.last().copied());
    for &root in roots {
        let old_degree = reconstructed.len() - 1;
        let leading = reconstructed[old_degree];
        reconstructed.push(leading);
        for degree in (1..=old_degree).rev() {
            reconstructed[degree] = reconstructed[degree - 1] - root * reconstructed[degree];
        }
        reconstructed[0] *= -root;
    }
    reconstructed
        .iter()
        .zip(coefficients)
        .try_fold(0.0_f64, |maximum, (&found, &wanted)| {
            let error = (found - wanted).norm();
            error.is_finite().then(|| maximum.max(error))
        })
        .unwrap_or(f64::NAN)
}

type ComplexPolynomialSolver =
    fn(&[PolynomialComplex], usize, f64) -> Result<Vec<PolynomialComplex>, String>;

fn checked_polynomial_roots(
    coefficients: &[PolynomialComplex],
    iterations: usize,
    tolerance: f64,
    solve: ComplexPolynomialSolver,
) -> Result<ComplexPolynomialRootSet, OpError> {
    let normalized = polynomial_problem(coefficients, iterations, tolerance, true)?;
    let raw = solve(&normalized, iterations, tolerance);
    crate::cancellation::checkpoint()?;
    let roots = raw.map_err(OpError::NotReduced)?;
    let degree = normalized.len() - 1;
    if roots.len() != degree {
        return Err(OpError::NotReduced(
            "the iterate count differs from polynomial degree".into(),
        ));
    }
    let backward_errors: Vec<_> = roots
        .iter()
        .map(|&root| polynomial_roots::backward_error(&normalized, root))
        .collect();
    let full_error = coefficient_error(&normalized, &roots);
    let reconstruction_limit = (16.0 * tolerance).max(32.0 * f64::EPSILON * degree as f64);
    if backward_errors
        .iter()
        .any(|&error| !error.is_finite() || error > tolerance)
        || !full_error.is_finite()
        || full_error > reconstruction_limit
    {
        return Err(OpError::NotReduced(
            "the candidate roots failed the backward error or coefficient reconstruction test"
                .into(),
        ));
    }
    Ok(ComplexPolynomialRootSet {
        roots,
        degree,
        backward_errors,
        coefficient_error: full_error,
    })
}

fn checked_laguerre_root(
    coefficients: &[PolynomialComplex],
    seed: PolynomialComplex,
    iterations: usize,
    tolerance: f64,
) -> Result<Option<ComplexPolynomialRoot>, OpError> {
    let normalized = polynomial_problem(coefficients, iterations, tolerance, false)?;
    if !seed.is_finite() {
        return Err(OpError::Domain(
            "Laguerre requires a finite complex seed".into(),
        ));
    }
    if normalized.len() == 1 {
        return Ok(None);
    }
    let raw = polynomial_roots::laguerre_root(&normalized, seed, iterations, tolerance);
    crate::cancellation::checkpoint()?;
    let root = raw.map_err(OpError::NotReduced)?;
    let backward_error = polynomial_roots::backward_error(&normalized, root);
    Ok(Some(ComplexPolynomialRoot { root, backward_error }))
}

binding!(
    polynomial_durand_kerner,
    "num_complex_polynomial_durand_kerner",
    "numerical_vector",
    "Finds approximate roots with Durand–Kerner in dimensionless coordinates. Coefficients are ascending. Reports the degree, root candidates including multiplicities, their coefficientwise relative backward errors and the largest coefficient error of the reconstructed polynomial after normalization. A nonzero constant has an empty root set. The zero polynomial, exceeded work budget and nonconvergence return errors.",
    "crate::kernels::complex_roots::durand_kerner",
    DimensionPolicy::Dimensionless,
    (coefficients: Vec<num_complex::Complex64>, max_steps: usize, tolerance: f64) -> Result<super::ComplexPolynomialRootSet, OpError>,
    |coefficients: Vec<num_complex::Complex64>, max_steps: usize, tolerance: f64| super::checked_polynomial_roots(&coefficients, max_steps, tolerance, super::polynomial_roots::durand_kerner)
);
binding!(
    polynomial_aberth_ehrlich,
    "num_complex_polynomial_aberth_ehrlich",
    "numerical_vector",
    "Finds approximate roots with Aberth–Ehrlich in dimensionless coordinates. Coefficients are ascending. Reports the degree, root candidates including multiplicities, their coefficientwise relative backward errors and the largest coefficient error of the reconstructed polynomial after normalization. A nonzero constant has an empty root set. The zero polynomial, exceeded work budget and nonconvergence return errors.",
    "crate::kernels::complex_roots::aberth_ehrlich",
    DimensionPolicy::Dimensionless,
    (coefficients: Vec<num_complex::Complex64>, max_steps: usize, tolerance: f64) -> Result<super::ComplexPolynomialRootSet, OpError>,
    |coefficients: Vec<num_complex::Complex64>, max_steps: usize, tolerance: f64| super::checked_polynomial_roots(&coefficients, max_steps, tolerance, super::polynomial_roots::aberth_ehrlich)
);
binding!(
    polynomial_laguerre,
    "num_complex_polynomial_laguerre",
    "numerical_vector",
    "Finds an approximate root of an ascending dimensionless polynomial from the supplied complex seed. The result includes its coefficientwise relative backward error. A nonzero constant returns None. The zero polynomial, exceeded work budget and nonconvergence return errors.",
    "crate::kernels::complex_roots::laguerre_root",
    DimensionPolicy::Dimensionless,
    (coefficients: Vec<num_complex::Complex64>, seed: num_complex::Complex64, max_steps: usize, tolerance: f64) -> Result<Option<super::ComplexPolynomialRoot>, OpError>,
    |coefficients: Vec<num_complex::Complex64>, seed: num_complex::Complex64, max_steps: usize, tolerance: f64| super::checked_laguerre_root(&coefficients, seed, max_steps, tolerance)
);

struct VerifiedPolynomialRootSet {
    disks: Vec<[f64; 3]>,
    complete: bool,
    degree: usize,
}
record!(VerifiedPolynomialRootSet {
    disks: Vec<[f64; 3]>,
    complete: bool,
    degree: usize
});

#[allow(
    clippy::single_call_fn,
    reason = "the bounded certificate workflow is separate from the operator registration"
)]
fn verified_polynomial_root_set(
    coefficients: &[PolynomialComplex],
    tolerance: f64,
) -> Result<VerifiedPolynomialRootSet, OpError> {
    if coefficients.len() > 129 {
        return Err(OpError::BudgetExceeded);
    }
    if !tolerance.is_finite()
        || tolerance <= 0.0
        || coefficients
            .iter()
            .any(|coefficient| !coefficient.is_finite())
    {
        return Err(OpError::Domain(
            "verified isolation needs finite coefficients and a positive tolerance".into(),
        ));
    }
    let degree = coefficients
        .iter()
        .rposition(|&coefficient| coefficient != PolynomialComplex::new(0.0, 0.0))
        .ok_or_else(|| OpError::Domain("the zero polynomial has no isolated roots".into()))?;
    crate::cancellation::checkpoint()?;
    let roots = polynomial_roots::isolate_roots_verified(coefficients, tolerance);
    crate::cancellation::checkpoint()?;
    let complete = roots.len() == degree;
    let disks = roots
        .into_iter()
        .map(|disk| {
            let center = disk.center();
            [center.re, center.im, disk.radius()]
        })
        .collect();
    Ok(VerifiedPolynomialRootSet {
        disks,
        complete,
        degree,
    })
}

binding!(
    polynomial_isolate_verified,
    "num_complex_polynomial_isolate_verified",
    "numerical_vector",
    "Certifies disjoint simple-root disks for an ascending dimensionless polynomial of degree at most 128. Each disk is [real, imaginary, radius]. complete is true only if the certificates account for every root; unresolved or multiple roots may produce a partial set.",
    "crate::kernels::complex_roots::isolate_roots_verified",
    DimensionPolicy::Dimensionless,
    (coefficients: Vec<num_complex::Complex64>, tolerance: f64) -> Result<super::VerifiedPolynomialRootSet, OpError>,
    |coefficients: Vec<num_complex::Complex64>, tolerance: f64| super::verified_polynomial_root_set(&coefficients, tolerance)
);
