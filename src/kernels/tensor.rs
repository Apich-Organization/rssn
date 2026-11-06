//! # Numerical Tensor Operations
//!
//! This module provides numerical tensor operations, primarily using `ndarray`
//! for efficient multi-dimensional array manipulation. It includes functions
//! for tensor contraction (tensordot), outer product, and Einstein summation (`einsum`).
//!
//! [`cp_als`] adds a CP (CANDECOMP/PARAFAC) decomposition by alternating
//! least squares: `T ~= sum_r lambda_r a_r o b_r o c_r`, with unit-norm factor
//! columns and the scale carried by `lambda`. [`hosvd`] adds the higher-order
//! SVD (Tucker decomposition) `T = C x_1 U_1 ... x_N U_N` with orthonormal
//! factors and an all-orthogonal core.
//!
//! # Limitations
//!
//! * Tensors are dense `ArrayD<f64>` of any rank. There is no sparse tensor
//!   type and no complex element type; a complex tensor has to be split into
//!   real and imaginary parts.
//! * CP is ill-posed:
//!   the best rank-`R` approximation need not exist (a tensor can have border
//!   rank below its rank), so a small residual after a sweep is not a
//!   certificate of optimality. The solver reports `converged`, `stagnated`,
//!   `diverged` and `rank_deficient` rather than hiding that. [`hosvd`] is
//!   exact before truncation and gives orthonormal factors, but its truncated
//!   form is not the best multilinear-rank approximation. [`hooi`] refines it
//!   by alternating orthogonal subspaces, but only local stationarity is
//!   reported, not global optimality; tensor trains are available in [`super::tensor_train`]; there is no CP
//!   with missing entries and no non-negativity constraint.
//! * [`hosvd`] uses a thin bidiagonal SVD for mode unfoldings, including
//!   rank-deficient ones, without forming a Gram matrix.
//! * The CP solver is plain ALS: it is not globally convergent and can be slow
//!   on ill-conditioned problems. Each mode uses one direct SVD of its
//!   Khatri-Rao design, followed by all right-hand sides together; no Gram
//!   matrix is formed. Mode unfoldings are cached, so dense storage dominates
//!   for large tensors.
//! * Factor columns are returned unit-norm with the weights in
//!   [`CpDecomposition::weights`]; a single term is only identified up to the
//!   usual permutation and scaling ambiguity, so compare reconstructions, not
//!   factors.

use ndarray::ArrayD;
use ndarray::IxDyn;

use crate::kernels::dense::Mat;
use crate::kernels::dense::svd;

/// Performs tensor contraction between two N-dimensional arrays (tensordot).
///
/// # Arguments
/// * `a`, `b` - The two tensors (`ndarray::ArrayD<f64>`) to contract.
/// * `axes_a`, `axes_b` - The axes to contract for tensor `a` and `b` respectively.
///
/// # Returns
/// The resulting contracted tensor as an `ndarray::ArrayD<f64>`.
///
/// # Errors
/// Returns an error if the number of axes to contract mismatch, or if the dimensions along these axes do not match.
pub fn tensordot(
    a: &ArrayD<f64>,
    b: &ArrayD<f64>,
    axes_a: &[usize],
    axes_b: &[usize],
) -> Result<ArrayD<f64>, String> {
    if axes_a.len() != axes_b.len() {
        return Err("Contracted axes \
                    must have the \
                    same length."
            .to_string());
    }

    for (&ax_a, &ax_b) in axes_a.iter().zip(axes_b.iter()) {
        if a.shape()[ax_a] != b.shape()[ax_b] {
            return Err(format!(
                "Dimension mismatch \
                 on contracted axes: \
                 {} != {}",
                a.shape()[ax_a],
                b.shape()[ax_b]
            ));
        }
    }

    let free_axes_a: Vec<_> = (0..a.ndim()).filter(|i| !axes_a.contains(i)).collect();

    let free_axes_b: Vec<_> = (0..b.ndim()).filter(|i| !axes_b.contains(i)).collect();

    let perm_a: Vec<_> = free_axes_a.iter().chain(axes_a.iter()).copied().collect();

    let perm_b: Vec<_> = axes_b.iter().chain(free_axes_b.iter()).copied().collect();

    let a_perm = a.clone().permuted_axes(perm_a);

    let b_perm = b.clone().permuted_axes(perm_b);

    let free_dim_a = free_axes_a.iter().map(|&i| a.shape()[i]).product::<usize>();

    let free_dim_b = free_axes_b.iter().map(|&i| b.shape()[i]).product::<usize>();

    let contracted_dim = axes_a.iter().map(|&i| a.shape()[i]).product::<usize>();

    let a_mat = a_perm
        .to_shape((free_dim_a, contracted_dim))
        .map_err(|e| e.to_string())?
        .to_owned();

    let b_mat = b_perm
        .to_shape((contracted_dim, free_dim_b))
        .map_err(|e| e.to_string())?
        .to_owned();

    let result_mat = a_mat.dot(&b_mat);

    let mut final_shape_dims = Vec::new();

    final_shape_dims.extend(free_axes_a.iter().map(|&i| a.shape()[i]));

    final_shape_dims.extend(free_axes_b.iter().map(|&i| b.shape()[i]));

    Ok(result_mat
        .to_shape(IxDyn(&final_shape_dims))
        .map_err(|e| e.to_string())?
        .to_owned())
}

/// Computes the outer product of two tensors.
///
/// The outer product of two tensors `A` (rank `r`) and `B` (rank `s`)
/// results in a new tensor `C` of rank `r + s`. Each component of `C`
/// is the product of a component from `A` and a component from `B`.
///
/// # Arguments
/// * `a` - The first tensor (`ndarray::ArrayD<f64>`).
/// * `b` - The second tensor (`ndarray::ArrayD<f64>`).
///
/// # Returns
/// The resulting outer product tensor as an `ndarray::ArrayD<f64>`.
///
/// # Errors
/// Returns an error if input tensors are not contiguous.
pub fn outer_product(
    a: &ArrayD<f64>,
    b: &ArrayD<f64>,
) -> Result<ArrayD<f64>, String> {
    let mut new_shape = a.shape().to_vec();

    new_shape.extend_from_slice(b.shape());

    let a_flat = a.as_slice().ok_or_else(|| {
        "Input tensor 'a' is not \
             contiguous"
            .to_string()
    })?;

    let b_flat = b.as_slice().ok_or_else(|| {
        "Input tensor 'b' is not \
             contiguous"
            .to_string()
    })?;

    let mut result_data = Vec::with_capacity(a.len() * b.len());

    for val_a in a_flat {
        for val_b in b_flat {
            result_data.push(val_a * val_b);
        }
    }

    ArrayD::from_shape_vec(IxDyn(&new_shape), result_data).map_err(|e| e.to_string())
}

/// Performs tensor-vector multiplication.
///
/// # Errors
/// Returns an error if the tensor has zero dimensions or if the last dimension mismatches the vector size.
pub fn tensor_vec_mul(
    tensor: &ArrayD<f64>,
    vector: &[f64],
) -> Result<ArrayD<f64>, String> {
    if tensor.ndim() < 1 {
        return Err("Tensor must \
                    have at least \
                    one dimension."
            .to_string());
    }

    let last_dim = tensor.shape()[tensor.ndim() - 1];

    if last_dim != vector.len() {
        return Err(format!(
            "Dimension mismatch: last \
             tensor dim {} != vector \
             length {}",
            last_dim,
            vector.len()
        ));
    }

    let vec_arr = ndarray::Array1::from_vec(vector.to_vec());

    let res = tensordot(tensor, &vec_arr.into_dyn(), &[tensor.ndim() - 1], &[0])?;

    Ok(res)
}

/// Computes the inner product of two tensors of the same shape.
///
/// # Errors
/// Returns an error if the shapes mismatch or if tensors are not contiguous.
pub fn inner_product(
    a: &ArrayD<f64>,
    b: &ArrayD<f64>,
) -> Result<f64, String> {
    if a.shape() != b.shape() {
        return Err("Tensors must \
                    have the same \
                    shape for inner \
                    product."
            .to_string());
    }

    let a_flat = a.as_slice().ok_or("Tensor 'a' is not contiguous")?;

    let b_flat = b.as_slice().ok_or("Tensor 'b' is not contiguous")?;

    Ok(a_flat.iter().zip(b_flat.iter()).map(|(x, y)| x * y).sum())
}

/// Contracts a single tensor along two specified axes (a generalised
/// trace): the result has rank `ndim - 2` and entries
/// `sum_t a[..., t, ..., t, ...]`, the remaining axes in their original
/// order.
///
/// # Errors
/// Returns an error if the axes are equal, out of range, or of different
/// lengths.
pub fn contract(
    a: &ArrayD<f64>,
    axis1: usize,
    axis2: usize,
) -> Result<ArrayD<f64>, String> {
    if axis1 == axis2 {
        return Err("Axes must be different for contraction.".to_string());
    }
    let (Some(&n), Some(&m)) = (a.shape().get(axis1), a.shape().get(axis2)) else {
        return Err("Contraction axis out of range.".to_string());
    };
    if n != m {
        return Err("Dimensions along contraction axes must be equal.".to_string());
    }
    // Remove the higher axis first so the lower index stays valid.
    let (low, high) = if axis1 < axis2 {
        (axis1, axis2)
    } else {
        (axis2, axis1)
    };
    let shape: Vec<usize> = a
        .shape()
        .iter()
        .enumerate()
        .filter(|&(i, _)| i != low && i != high)
        .map(|(_, &d)| d)
        .collect();
    let mut out = ArrayD::<f64>::zeros(IxDyn(&shape));
    for t in 0..n {
        let diagonal = a.index_axis(ndarray::Axis(high), t);
        let diagonal = diagonal.index_axis(ndarray::Axis(low), t);
        out += &diagonal;
    }
    Ok(out)
}

/// Computes the Frobenius norm of a tensor.
#[must_use]
pub fn norm(a: &ArrayD<f64>) -> f64 {
    a.iter().fold(0.0_f64, |length, &x| length.hypot(x))
}

use serde::Deserialize;
use serde::Serialize;

/// A serializable representation of an N-dimensional tensor.
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct TensorData {
    /// Dimensions of the tensor.
    pub shape: Vec<usize>,
    /// Flat vector of tensor data.
    pub data: Vec<f64>,
}

impl From<&ArrayD<f64>> for TensorData {
    fn from(arr: &ArrayD<f64>) -> Self {
        Self {
            shape: arr.shape().to_vec(),
            data: arr.clone().into_raw_vec_and_offset().0,
        }
    }
}

impl TensorData {
    /// Converts back to an `ndarray::ArrayD`.
    ///
    /// # Errors
    /// Returns an error if the shape and data are inconsistent.
    pub fn to_arrayd(&self) -> Result<ArrayD<f64>, String> {
        ArrayD::from_shape_vec(IxDyn(&self.shape), self.data.clone()).map_err(|e| e.to_string())
    }
}

/// Configuration of the CP-ALS solver [`cp_als`].
#[derive(Debug, Clone)]
pub struct CpOptions {
    /// Number of rank-one components `R`; must be at least 1.
    pub rank: usize,
    /// Maximum number of alternating-least-squares sweeps.
    pub max_sweeps: usize,
    /// Tolerance, relative to `||T||_F`, for the fitted residual and for the
    /// per-sweep improvement.
    pub tol: f64,
    /// Seed of the uniform initialisation; the same seed reproduces the run.
    pub seed: u64,
    /// A factor weight above `divergence_limit * (1 + ||T||_F)` is reported as
    /// divergence.
    pub divergence_limit: f64,
}

impl CpOptions {
    /// Options for rank `rank` with `500` sweeps, `tol = 1e-12`, seed `0` and
    /// a divergence limit of `1e6`.
    #[must_use]
    pub const fn new(rank: usize) -> Self {
        Self {
            rank,
            max_sweeps: 500,
            tol: 1e-12,
            seed: 0,
            divergence_limit: 1e6,
        }
    }
}

/// Outcome of [`cp_als`].
#[allow(clippy::struct_excessive_bools)]
#[derive(Debug, Clone)]
pub struct CpDecomposition {
    /// Factor matrices `A^(k)`, one per mode, each `n_k x rank` with
    /// unit-norm columns.
    pub factors: Vec<Mat>,
    /// Column norms removed from the factors, one per component.
    pub weights: Vec<f64>,
    /// `sum_r weights[r] * a_r o b_r o c_r`, same shape as the input.
    pub reconstruction: ArrayD<f64>,
    /// Frobenius norm of `T - reconstruction`.
    pub residual: f64,
    /// `residual / ||T||_F`.
    pub relative_residual: f64,
    /// Number of completed sweeps.
    pub sweeps: usize,
    /// The relative residual reached `tol`.
    pub converged: bool,
    /// The per-sweep improvement fell to `tol` while the relative residual
    /// stayed above it: ALS stopped making progress without fitting the data.
    pub stagnated: bool,
    /// A factor weight left `divergence_limit * (1 + ||T||_F)` (or a
    /// non-finite value appeared), the signature of a degenerate CP problem.
    pub diverged: bool,
    /// A Khatri-Rao design matrix was numerically rank deficient and its
    /// update used a minimum-norm SVD solve; components may not be identifiable.
    pub rank_deficient: bool,
}

/// CP (CANDECOMP/PARAFAC) decomposition by alternating least squares.
///
/// Factorises a dense tensor as `T ~= sum_r lambda_r a_r o b_r o c_r` for
/// `options.rank` components, by solving one matricised least-squares problem
/// per mode per sweep (the standard ALS update on the mode-`k` unfolding
/// against the Khatri-Rao product of the other factors). Factors are
/// initialised uniformly at random from `options.seed` and normalised after
/// every update, with the removed norms kept in
/// [`CpDecomposition::weights`].
///
/// A rank-`R` tensor whose factors span a full-rank problem is recovered to
/// rounding error, but the factors themselves are only determined up to
/// permutation and scaling: compare reconstructions, not factors.
///
/// CP is ill-posed, so the result says what happened rather than assuming
/// success. `converged` means the residual reached `tol`; `stagnated` means
/// sweeps stopped improving while the residual stayed above `tol`; `diverged`
/// means a weight blew up, which is how a degenerate problem (border rank
/// below rank) shows itself; `rank_deficient` means the least-squares system
/// was numerically rank deficient in its direct SVD solve.
///
/// # Errors
/// Returns an error for a tensor of rank below 2, an empty dimension,
/// non-finite entries/options, a zero rank, zero sweeps, or rescaling overflow.
/// Updates solve the Khatri-Rao least-squares system directly by SVD, once
/// per mode, avoiding normal equations and squared condition numbers. Input
/// normalization protects uniformly very small or very large tensors.
pub fn cp_als(
    tensor: &ArrayD<f64>,
    options: &CpOptions,
) -> Result<CpDecomposition, String> {
    let ndim = tensor.ndim();
    if ndim < 2 {
        return Err("CP decomposition needs a tensor of rank at least 2.".to_string());
    }
    if options.rank == 0 {
        return Err("CP rank must be at least 1.".to_string());
    }
    if options.max_sweeps == 0 {
        return Err("max_sweeps must be at least 1.".to_string());
    }
    if tensor.shape().contains(&0) {
        return Err("CP decomposition needs non-empty dimensions.".to_string());
    }
    if tensor.iter().any(|v| !v.is_finite()) {
        return Err("Tensor entries must be finite.".to_string());
    }

    if !options.tol.is_finite()
        || options.tol < 0.0
        || !options.divergence_limit.is_finite()
        || options.divergence_limit <= 0.0
    {
        return Err(
            "CP tolerance must be finite/nonnegative and divergence limit finite/positive.".into(),
        );
    }
    let scale = tensor.iter().fold(0.0_f64, |s, x| s.max(x.abs()));
    let scale = if scale == 0.0 { 1.0 } else { scale };
    let normalized = tensor.mapv(|x| x / scale);
    let rank = options.rank;
    let shape = tensor.shape().to_vec();
    let unfoldings = mode_unfoldings(&normalized)?;
    let t_norm = {
        let n = norm(&normalized);
        if n > 0.0 { n } else { 1.0 }
    };
    let limit = options.divergence_limit * (scale.recip() + t_norm);

    let mut factors = init_factors(&shape, rank, options.seed);
    let mut weights = vec![1.0; rank];

    let mut sweeps = 0_usize;
    let mut residual = f64::INFINITY;
    let mut previous = f64::INFINITY;
    let mut last_improvement = f64::INFINITY;
    let mut converged = false;
    let mut diverged = false;
    let mut rank_deficient = false;

    for sweep in 1..=options.max_sweeps {
        for (mode, unfold) in unfoldings.iter().enumerate() {
            let kr = khatri_rao(&factors, mode, &shape);
            let decomposition = svd(&kr).map_err(|e| e.to_string())?;
            let cutoff = kr.rows.max(kr.cols) as f64
                * f64::EPSILON
                * decomposition.s.first().copied().unwrap_or(0.0);
            let active = decomposition.s.iter().filter(|&&s| s > cutoff).count();
            rank_deficient |= active < rank;
            let projected = unfold.matmul(&decomposition.u).map_err(|e| e.to_string())?;
            let rows = unfold.rows;
            let mut updated = Mat::zeros(rows, rank);
            for row in 0..rows {
                for col in 0..rank {
                    let mut value = 0.0;
                    for j in 0..active {
                        value +=
                            projected.at(row, j) / decomposition.s[j] * decomposition.v.at(col, j);
                    }
                    updated.set(row, col, value);
                }
            }
            for col in 0..rank {
                let length = column_norm(&updated, col);
                if length > 0.0 {
                    for row in 0..rows {
                        updated.set(row, col, updated.at(row, col) / length);
                    }
                }
                if let Some(w) = weights.get_mut(col) {
                    *w = length;
                }
            }
            if let Some(slot) = factors.get_mut(mode) {
                *slot = updated;
            }
        }
        sweeps = sweep;
        let flat = reconstruct(&factors, &weights, &shape);
        residual = normalized
            .iter()
            .zip(flat.iter())
            .fold(0.0_f64, |length, (t, r)| length.hypot(t - r));
        if !residual.is_finite() {
            diverged = true;
            break;
        }
        if previous.is_finite() && previous > 0.0 {
            last_improvement = (previous - residual) / previous;
        }
        if residual <= options.tol * t_norm {
            converged = true;
            break;
        }
        if weights.iter().any(|w| !w.is_finite() || w.abs() > limit) {
            diverged = true;
            break;
        }
        previous = residual;
    }

    let relative_residual = residual / t_norm;
    let flat = reconstruct(&factors, &weights, &shape)
        .into_iter()
        .map(|x| x * scale)
        .collect::<Vec<_>>();
    for weight in &mut weights {
        *weight *= scale;
    }
    let residual = residual * scale;
    if !flat.iter().chain(weights.iter()).all(|x| x.is_finite()) || !residual.is_finite() {
        return Err("CP result overflowed while rescaling.".into());
    }
    let reconstruction = ArrayD::from_shape_vec(IxDyn(&shape), flat).map_err(|e| e.to_string())?;
    Ok(CpDecomposition {
        factors,
        weights,
        reconstruction,
        residual,
        relative_residual,
        sweeps,
        converged,
        stagnated: !converged && !diverged && last_improvement <= options.tol,
        diverged,
        rank_deficient,
    })
}

/// A Tucker decomposition computed by the higher-order SVD (HOSVD).
///
/// `T = C x_1 U_1 x_2 ... x_N U_N`, where `factors[k]` is the `n_k x r_k`
/// matrix of the leading left singular vectors of the mode-`k` unfolding and
/// `core` is the `r_1 x ... x r_N` core tensor. The factors have orthonormal
/// columns and the core is all-orthogonal (its mode-`k` slices are mutually
/// orthogonal); the truncated HOSVD is *not* the best multilinear-rank
/// approximation — [`hooi`] can refine it locally — but it is the canonical starting
/// point and its truncation error is bounded by the discarded singular values.
#[derive(Debug, Clone, PartialEq)]
pub struct Hosvd {
    /// Core tensor, `r_1 x ... x r_N`.
    pub core: ArrayD<f64>,
    /// One `n_k x r_k` factor per mode, with orthonormal columns.
    pub factors: Vec<Mat>,
    /// The retained singular values of each mode unfolding, descending.
    pub singular_values: Vec<Vec<f64>>,
    /// Retained rank per mode.
    pub ranks: Vec<usize>,
}

/// Tucker approximation refined by higher-order orthogonal iteration (HOOI).
///
/// Factors have orthonormal columns. Unlike HOSVD, the core need not be
/// all-orthogonal. Convergence denotes stabilization of the local objective,
/// not global optimality or a small approximation residual.
#[derive(Debug, Clone, PartialEq)]
pub struct Tucker {
    /// Reduced core tensor.
    pub core: ArrayD<f64>,
    /// Orthonormal factor for each mode.
    pub factors: Vec<Mat>,
    /// Actual retained ranks (bounded by the requested ranks).
    pub ranks: Vec<usize>,
    /// Completed alternating sweeps.
    pub sweeps: usize,
    /// Directly measured relative Frobenius reconstruction residual.
    pub relative_error: f64,
    /// Initial HOSVD residual followed by each completed sweep's residual.
    pub error_history: Vec<f64>,
    /// Whether the objective stabilized within the requested tolerance.
    pub converged: bool,
}

impl Tucker {
    /// Reconstruct the dense approximation.
    ///
    /// # Errors
    /// Returns an error if the public core/factor dimensions are inconsistent.
    pub fn reconstruct(&self) -> Result<ArrayD<f64>, String> {
        let mut out = self.core.clone();
        for (mode, factor) in self.factors.iter().enumerate() {
            out = mode_product(&out, mode, factor)?;
        }
        Ok(out)
    }
}

/// Refine a truncated HOSVD using higher-order orthogonal iteration.
///
/// Each coordinate step computes the leading singular subspace after
/// projecting all other modes. Ranks are upper bounds and may be clamped to
/// the available subspace. Initialization is deterministic HOSVD; HOOI can
/// converge to a local optimum or stationary point. Inspect `relative_error`
/// independently from `converged`. Exhausting the sweep budget returns the
/// best observed approximation with `converged = false`. Residuals near
/// floating-point roundoff also terminate the iteration.
///
/// Input normalization protects intermediate norm calculations at uniformly
/// small/large scales. No Gram matrices or condition-number squaring are used.
/// Zero tensors use a rank-one zero core.
///
/// # Errors
/// Non-finite data, empty modes, invalid
/// ranks/tolerance/budget, SVD failure, or arithmetic overflow return an error.
pub fn hooi(
    tensor: &ArrayD<f64>,
    ranks: &[usize],
    max_sweeps: usize,
    tolerance: f64,
) -> Result<Tucker, String> {
    let shape = tensor.shape();
    if ranks.len() != shape.len()
        || ranks.contains(&0)
        || shape.contains(&0)
        || max_sweeps == 0
        || !tolerance.is_finite()
        || tolerance < 0.0
    {
        return Err("invalid Tucker ranks, tolerance or sweep budget".into());
    }
    if !tensor.iter().all(|x| x.is_finite()) {
        return Err("Tucker input contains non-finite values".into());
    }
    let scale = tensor.iter().fold(0.0_f64, |s, x| s.max(x.abs()));
    if scale == 0.0 || shape.is_empty() {
        let factors = shape
            .iter()
            .map(|&n| {
                let mut u = Mat::zeros(n, 1);
                u.set(0, 0, 1.0);
                u
            })
            .collect();
        return Ok(Tucker {
            core: if shape.is_empty() {
                tensor.clone()
            } else {
                ArrayD::zeros(IxDyn(&vec![1; shape.len()]))
            },
            factors,
            ranks: vec![1; shape.len()],
            sweeps: 0,
            relative_error: 0.0,
            error_history: vec![0.0],
            converged: true,
        });
    }
    let input = tensor.mapv(|x| x / scale);
    let input_norm = input.iter().fold(0.0_f64, |s, &x| s.hypot(x));
    let initial = hosvd(&input, Some(ranks))?;
    let mut factors = initial.factors;
    let mut actual_ranks = initial.ranks;
    let mut core = initial.core;
    let measure = |core: &ArrayD<f64>, factors: &[Mat]| -> Result<f64, String> {
        let mut reconstructed = core.clone();
        for (mode, u) in factors.iter().enumerate() {
            reconstructed = mode_product(&reconstructed, mode, u)?;
        }
        let residual = input
            .iter()
            .zip(reconstructed.iter())
            .fold(0.0_f64, |s, (a, b)| s.hypot(a - b))
            / input_norm;
        if !residual.is_finite() {
            return Err("Tucker reconstruction overflow".into());
        }
        Ok(residual)
    };
    let mut previous = measure(&core, &factors)?;
    let mut best_error = previous;
    let mut best_core = core.clone();
    let mut best_factors = factors.clone();
    let mut best_ranks = actual_ranks.clone();
    let mut history = vec![previous];
    let noise_floor = 32.0 * f64::EPSILON * shape.len().max(1) as f64;
    let mut converged = previous <= noise_floor;
    let mut sweeps = 0;
    while sweeps < max_sweeps && !converged {
        for mode in 0..shape.len() {
            let mut projected = input.clone();
            for (other, u) in factors.iter().enumerate() {
                if other != mode {
                    projected = mode_product(&projected, other, &u.transpose())?;
                }
            }
            let mut permutation = vec![mode];
            permutation.extend((0..shape.len()).filter(|&k| k != mode));
            let unfolded = projected
                .permuted_axes(permutation)
                .as_standard_layout()
                .into_owned();
            let rows = shape[mode];
            let matrix = Mat {
                rows,
                cols: unfolded.len() / rows,
                data: unfolded.into_raw_vec_and_offset().0,
            };
            let d = svd(&matrix).map_err(|e| format!("Tucker mode {mode}: {e}"))?;
            let rank = actual_ranks[mode].min(d.s.len());
            if rank == 0 {
                return Err("Tucker projection has an empty singular subspace".into());
            }
            let mut u = Mat::zeros(rows, rank);
            for i in 0..rows {
                for j in 0..rank {
                    u.set(i, j, d.u.at(i, j));
                }
            }
            factors[mode] = u;
            actual_ranks[mode] = rank;
        }
        core.clone_from(&input);
        for (mode, u) in factors.iter().enumerate() {
            core = mode_product(&core, mode, &u.transpose())?;
        }
        let error = measure(&core, &factors)?;
        history.push(error);
        sweeps += 1;
        if error < best_error {
            best_error = error;
            best_core.clone_from(&core);
            best_factors.clone_from(&factors);
            best_ranks.clone_from(&actual_ranks);
        }
        converged = error <= noise_floor
            || (previous - error).abs() <= tolerance * previous.max(f64::EPSILON);
        previous = error;
    }
    best_core.mapv_inplace(|x| x * scale);
    if !best_core.iter().all(|x| x.is_finite()) {
        return Err("Tucker core rescaling overflow".into());
    }
    Ok(Tucker {
        core: best_core,
        factors: best_factors,
        ranks: best_ranks,
        sweeps,
        relative_error: best_error,
        error_history: history,
        converged,
    })
}

impl Hosvd {
    /// Reconstructs the tensor: `C x_1 U_1 ... x_N U_N`.
    ///
    /// # Errors
    /// Fails only on an inconsistent shape, which a value produced by
    /// [`hosvd`] cannot have.
    pub fn reconstruct(&self) -> Result<ArrayD<f64>, String> {
        let mut out = self.core.clone();
        for (mode, factor) in self.factors.iter().enumerate() {
            out = mode_product(&out, mode, factor)?;
        }
        Ok(out)
    }
}

/// The mode-`k` product `T x_k M`: `M` is `r x n_k` and the result has `r` in
/// place of mode `k`.
fn mode_product(
    tensor: &ArrayD<f64>,
    mode: usize,
    matrix: &Mat,
) -> Result<ArrayD<f64>, String> {
    let ndim = tensor.ndim();
    if mode >= ndim {
        return Err(String::from("mode out of range"));
    }
    let shape = tensor.shape().to_vec();
    let rows = *shape.get(mode).ok_or("mode out of range")?;
    if matrix.cols != rows {
        return Err(String::from("mode-product shape mismatch"));
    }
    if rows == 0 {
        return Err(String::from("cannot contract an empty mode"));
    }
    let mut permutation = Vec::with_capacity(ndim);
    permutation.push(mode);
    permutation.extend((0..ndim).filter(|&axis| axis != mode));
    let standard = tensor
        .clone()
        .permuted_axes(permutation.clone())
        .as_standard_layout()
        .into_owned();
    let data = standard.into_raw_vec_and_offset().0;
    let columns = data.len() / rows;
    let mut product = vec![0.0; matrix.rows * columns];
    for row in 0..matrix.rows {
        for slice in 0..rows {
            let coefficient = matrix.at(row, slice);
            for column in 0..columns {
                product[row * columns + column] += coefficient * data[slice * columns + column];
            }
        }
    }
    let mut contracted = Vec::with_capacity(ndim);
    contracted.push(matrix.rows);
    for (axis, &extent) in shape.iter().enumerate() {
        if axis != mode {
            contracted.push(extent);
        }
    }
    let result = ArrayD::from_shape_vec(contracted, product).map_err(|e| e.to_string())?;
    // `result` has the contracted mode first; undo the permutation.
    let mut inverse = vec![0_usize; ndim];
    for (position, &axis) in permutation.iter().enumerate() {
        inverse[axis] = position;
    }
    Ok(result.permuted_axes(inverse).into_owned())
}

/// Appends a column to `u`.
#[allow(clippy::single_call_fn)]
fn push_column(
    u: &mut Mat,
    column: &[f64],
) {
    let next = u.cols;
    let mut grown = Mat::zeros(u.rows, next + 1);
    for row in 0..u.rows {
        for old in 0..next {
            grown.set(row, old, u.at(row, old));
        }
        grown.set(row, next, column.get(row).copied().unwrap_or(0.0));
    }
    *u = grown;
}

/// The non-negligible leading left singular data of a dense matrix.
///
/// Numerical left singular subspace without forming a Gram matrix.
fn reduced_svd(a: &Mat) -> Result<(Mat, Vec<f64>), String> {
    let decomposition = svd(a).map_err(|error| error.to_string())?;
    let sigma_max = decomposition.s.first().copied().unwrap_or(0.0);
    let tolerance = a.rows.max(a.cols) as f64 * f64::EPSILON * sigma_max;
    let mut u = Mat::zeros(a.rows, 0);
    let mut singular = Vec::new();
    for (column, &sigma) in decomposition.s.iter().enumerate() {
        if sigma <= tolerance {
            continue;
        }
        let entries: Vec<f64> = (0..a.rows)
            .map(|row| decomposition.u.at(row, column))
            .collect();
        push_column(&mut u, &entries);
        singular.push(sigma);
    }
    if singular.is_empty() {
        return Err("the matrix has no nonzero singular value".into());
    }
    Ok((u, singular))
}

/// The higher-order SVD (Tucker decomposition) of a dense tensor.
///
/// For each mode the tensor is unfolded to `n_k x prod(rest)`, the leading left
/// singular vectors of the unfolding become the factor, and the core is the
/// tensor contracted with every factor. `ranks` may be `None` for the numerical
/// multilinear rank of the unfoldings (at most `min(n_k, prod(rest))`, and less
/// when an unfolding is rank deficient in floating point); a requested rank
/// above that is clamped to it. The decomposition is exact before truncation
/// and gives `||T - C x U||_F <= sqrt(sum_k sum_{j > r_k} sigma_{k,j}^2)`.
///
/// # Errors
/// Fails when `ranks` does not have one entry per mode, when a requested rank
/// is zero, when a mode is empty, or when a mode unfolding has no nonzero
/// singular value.
pub fn hosvd(
    tensor: &ArrayD<f64>,
    ranks: Option<&[usize]>,
) -> Result<Hosvd, String> {
    let shape = tensor.shape().to_vec();
    if shape.is_empty() {
        return Ok(Hosvd {
            core: tensor.clone(),
            factors: Vec::new(),
            singular_values: Vec::new(),
            ranks: Vec::new(),
        });
    }
    if let Some(requested) = ranks {
        if requested.len() != shape.len() {
            return Err(String::from("ranks must have one entry per mode"));
        }
        if requested.contains(&0) {
            return Err(String::from("a requested rank is zero"));
        }
    }
    let unfoldings = mode_unfoldings(tensor)?;
    let mut factors = Vec::with_capacity(shape.len());
    let mut singular_values = Vec::with_capacity(shape.len());
    let mut kept = Vec::with_capacity(shape.len());
    for (mode, unfolding) in unfoldings.iter().enumerate() {
        let (left, singular) =
            reduced_svd(unfolding).map_err(|error| format!("mode {mode} unfolding: {error}"))?;
        let available = singular.len();
        let rank = ranks.map_or(available, |requested| requested[mode].min(available));
        if rank == 0 {
            return Err(String::from("a mode unfolding has no singular value"));
        }
        let mut factor = Mat::zeros(unfolding.rows, rank);
        for row in 0..unfolding.rows {
            for column in 0..rank {
                factor.set(row, column, left.at(row, column));
            }
        }
        singular_values.push(singular.iter().take(rank).copied().collect());
        kept.push(rank);
        factors.push(factor);
    }
    let mut core = tensor.clone();
    for (mode, factor) in factors.iter().enumerate() {
        core = mode_product(&core, mode, &factor.transpose())?;
    }
    Ok(Hosvd {
        core,
        factors,
        singular_values,
        ranks: kept,
    })
}

/// Mode-`k` unfoldings, `n_k x prod(rest)`, with the mode first and the
/// remaining modes in their original order.
#[allow(clippy::single_call_fn)]
fn mode_unfoldings(tensor: &ArrayD<f64>) -> Result<Vec<Mat>, String> {
    let ndim = tensor.ndim();
    let total = tensor.len();
    let mut out = Vec::with_capacity(ndim);
    for mode in 0..ndim {
        let mut perm = Vec::with_capacity(ndim);
        perm.push(mode);
        perm.extend((0..ndim).filter(|&i| i != mode));
        let permuted = tensor.clone().permuted_axes(perm);
        let standard = permuted.as_standard_layout().into_owned();
        let rows = tensor
            .shape()
            .get(mode)
            .copied()
            .ok_or("Mode out of range.")?;
        let cols = total.checked_div(rows).ok_or("Empty tensor dimension.")?;
        out.push(Mat {
            rows,
            cols,
            data: standard.into_raw_vec_and_offset().0,
        });
    }
    Ok(out)
}

/// Uniform `[-1, 1)` initial factors, one `n_k x rank` matrix per mode.
#[allow(clippy::single_call_fn)]
fn init_factors(
    shape: &[usize],
    rank: usize,
    seed: u64,
) -> Vec<Mat> {
    use rand_v09::Rng;
    use rand_v09::SeedableRng;

    let mut rng = rand_v09::rngs::StdRng::seed_from_u64(seed);
    shape
        .iter()
        .map(|&n| {
            let mut m = Mat::zeros(n, rank);
            for i in 0..n {
                for j in 0..rank {
                    m.set(i, j, rng.random_range(-1.0..1.0));
                }
            }
            m
        })
        .collect()
}

/// Khatri-Rao product of every factor except `skip`, with the rows ordered
/// like the corresponding mode unfolding.
#[allow(clippy::single_call_fn)]
fn khatri_rao(
    factors: &[Mat],
    skip: usize,
    shape: &[usize],
) -> Mat {
    let dims: Vec<usize> = shape
        .iter()
        .enumerate()
        .filter(|&(i, _)| i != skip)
        .map(|(_, &d)| d)
        .collect();
    let rows = dims.iter().product::<usize>();
    let rank = factors.first().map_or(0, |f| f.cols);
    let mut out = Mat::zeros(rows, rank);
    let mut index = vec![0_usize; dims.len()];
    for row in 0..rows {
        let mut rest = row;
        for d in (0..dims.len()).rev() {
            let dim = dims.get(d).copied().unwrap_or(1);
            if let Some(slot) = index.get_mut(d) {
                *slot = rest % dim;
            }
            rest /= dim;
        }
        for r in 0..rank {
            let mut product = 1.0;
            let mut position = 0_usize;
            for (mode, f) in factors.iter().enumerate() {
                if mode == skip {
                    continue;
                }
                let i = index.get(position).copied().unwrap_or(0);
                product *= f.at(i, r);
                position += 1;
            }
            out.set(row, r, product);
        }
    }
    out
}

/// Euclidean norm of column `col`.
#[allow(clippy::single_call_fn)]
fn column_norm(
    a: &Mat,
    col: usize,
) -> f64 {
    (0..a.rows).fold(0.0_f64, |length, i| length.hypot(a.at(i, col)))
}

/// Reconstruction in logical (row-major) order.
#[allow(clippy::single_call_fn)]
fn reconstruct(
    factors: &[Mat],
    weights: &[f64],
    shape: &[usize],
) -> Vec<f64> {
    let total = shape.iter().product::<usize>();
    let rank = weights.len();
    let mut data = vec![0.0; total];
    let mut index = vec![0_usize; shape.len()];
    for (flat, slot) in data.iter_mut().enumerate() {
        let mut rest = flat;
        for d in (0..shape.len()).rev() {
            let dim = shape.get(d).copied().unwrap_or(1);
            if let Some(entry) = index.get_mut(d) {
                *entry = rest % dim;
            }
            rest /= dim;
        }
        let mut sum = 0.0;
        for r in 0..rank {
            let mut product = weights.get(r).copied().unwrap_or(0.0);
            for (mode, f) in factors.iter().enumerate() {
                let i = index.get(mode).copied().unwrap_or(0);
                product *= f.at(i, r);
            }
            sum += product;
        }
        *slot = sum;
    }
    data
}

#[cfg(test)]
mod tests {

    use ndarray::array;

    use super::*;

    #[test]
    fn test_tensordot() {
        let a = array![[1.0, 2.0], [3.0, 4.0]].into_dyn();

        let b = array![[5.0, 6.0], [7.0, 8.0]].into_dyn();

        let res = tensordot(&a, &b, &[1], &[0]).unwrap();

        // Standard matrix multiplication
        assert_eq!(res.shape(), &[2, 2]);

        assert!((res[[0, 0]] - (1.0 * 5.0 + 2.0 * 7.0)).abs() < 1e-9);
    }

    #[test]
    fn test_outer_product() {
        let a = array![1.0, 2.0].into_dyn();

        let b = array![3.0, 4.0].into_dyn();

        let res = outer_product(&a, &b).unwrap();

        assert_eq!(res.shape(), &[2, 2]);

        assert!((res[[0, 0]] - 3.0).abs() < 1e-9);

        assert!((res[[1, 1]] - 8.0).abs() < 1e-9);
    }

    /// One deterministic `n x rank` factor.
    fn factor(
        n: usize,
        rank: usize,
        seed: u64,
    ) -> Mat {
        init_factors(&[n], rank, seed).into_iter().next().unwrap()
    }

    /// Dense tensor `sum_r weights[r] * a_r o b_r o ...`.
    fn build(
        factors: &[Mat],
        weights: &[f64],
        shape: &[usize],
    ) -> ArrayD<f64> {
        ArrayD::from_shape_vec(IxDyn(shape), reconstruct(factors, weights, shape)).unwrap()
    }

    /// Independent recomputation of `||T - reconstruction|| / ||T||`.
    fn measured_relative_residual(
        t: &ArrayD<f64>,
        d: &CpDecomposition,
    ) -> f64 {
        let num = t
            .iter()
            .zip(d.reconstruction.iter())
            .map(|(a, b)| (a - b) * (a - b))
            .sum::<f64>()
            .sqrt();
        num / norm(t)
    }

    /// A rank-`R` tensor built from known factors is recovered to rounding
    /// error. The factors themselves are only identified up to a column
    /// permutation and scaling, so this compares the reconstruction.
    #[test]
    fn test_cp_exact_recovery() {
        let shape = [3_usize, 4, 5];
        let rank = 3;
        let factors: Vec<Mat> = shape
            .iter()
            .enumerate()
            .map(|(mode, &n)| factor(n, rank, 11 + mode as u64))
            .collect();
        let weights = [1.0, -0.7, 0.4];
        let t = build(&factors, &weights, &shape);

        let mut options = CpOptions::new(rank);
        options.max_sweeps = 3000;
        options.tol = 1e-13;
        options.seed = 5;
        let d = cp_als(&t, &options).unwrap();

        assert!(d.converged, "not converged: {d:?}");
        assert!(!d.diverged);
        assert!(
            d.relative_residual < 1e-10,
            "relative residual {:e}",
            d.relative_residual
        );
        let max_entry_error = t
            .iter()
            .zip(d.reconstruction.iter())
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);
        assert!(
            max_entry_error < 1e-10,
            "max entry error {max_entry_error:e}"
        );
        for f in &d.factors {
            for col in 0..rank {
                assert!((column_norm(f, col) - 1.0).abs() < 1e-12);
            }
        }
    }

    /// More components never fit worse: the residual falls monotonically with
    /// `R` on a fixed tensor.
    #[test]
    fn test_cp_error_decreases_with_rank() {
        let shape = [5_usize, 5, 5];
        let rank = 4;
        let factors: Vec<Mat> = shape
            .iter()
            .enumerate()
            .map(|(mode, &n)| factor(n, rank, 31 + mode as u64))
            .collect();
        let weights = [1.0, 0.3, 0.1, 0.03];
        let t = build(&factors, &weights, &shape);

        let mut errors = Vec::new();
        for r in 1..=rank {
            let mut options = CpOptions::new(r);
            options.max_sweeps = 2000;
            options.tol = 1e-13;
            options.seed = 21;
            let d = cp_als(&t, &options).unwrap();
            errors.push(d.relative_residual);
            println!("R = {r}: relative residual = {:.6e}", d.relative_residual);
        }
        for (i, pair) in errors.windows(2).enumerate() {
            assert!(
                pair[1] < pair[0],
                "R = {} error {:e} is not below R = {} error {:e}",
                i + 2,
                pair[1],
                i + 1,
                pair[0]
            );
        }
        assert!(errors[3] < 1e-10, "rank-4 fit {:e}", errors[3]);
    }

    /// `T = e2 o e1 o e1 + e1 o e3 o e1 + e1 o e1 o e4` has rank 3 and border
    /// rank 2: `(1/eps)[(e1+eps e2) o (e1+eps e3) o (e1+eps e4) - e1 o e1 o e1]`
    /// tends to `T` with factors of norm `~1/eps`. The best rank-2
    /// approximation therefore does not exist, and the solver must say so
    /// instead of reporting a fit it does not have.
    #[test]
    fn test_cp_degenerate_reported() {
        // C order in a 2 x 3 x 4 tensor: (1,0,0), (0,2,0), (0,0,3).
        let mut data = vec![0.0_f64; 2 * 3 * 4];
        data[12] = 1.0;
        data[8] = 1.0;
        data[3] = 1.0;
        let t = ArrayD::from_shape_vec(IxDyn(&[2, 3, 4]), data).unwrap();

        let mut options = CpOptions::new(2);
        options.max_sweeps = 400;
        options.seed = 3;
        // The border-rank degeneracy grows the weights only like sqrt(sweeps),
        // so the default 1e6 limit never fires; a tight limit makes the solver
        // report the blow-up it is actually exhibiting.
        options.divergence_limit = 1.0;
        let d = cp_als(&t, &options).unwrap();

        assert!(
            !d.converged,
            "degenerate problem claimed to converge: {d:?}"
        );
        assert!(
            d.stagnated || d.diverged,
            "no diagnosis of the degeneracy: {d:?}"
        );
        let max_weight = d.weights.iter().fold(0.0_f64, |m, w| m.max(w.abs()));
        assert!(max_weight > 1.0, "weights did not grow: {max_weight:e}");
        let honest = measured_relative_residual(&t, &d);
        assert!(
            (honest - d.relative_residual).abs() < 1e-12,
            "reported {:e} vs measured {:e}",
            d.relative_residual,
            honest
        );
    }

    /// ALS on a well-posed problem is insensitive to the initialisation: two
    /// seeds reach the same reconstruction error.
    #[test]
    fn test_cp_seeds_agree() {
        let shape = [4_usize, 4, 4];
        let rank = 3;
        let factors: Vec<Mat> = shape
            .iter()
            .enumerate()
            .map(|(mode, &n)| factor(n, rank, 41 + mode as u64))
            .collect();
        let weights = [1.0, 0.5, 0.25];
        let t = build(&factors, &weights, &shape);

        let mut a = CpOptions::new(rank);
        // ALS can spend thousands of sweeps in a swamp before it converges;
        // seed 1 needs ~7.4e3 here, so the budget must not cut it short.
        a.max_sweeps = 20_000;
        a.tol = 1e-13;
        a.seed = 1;
        let mut b = a.clone();
        b.seed = 2;
        let da = cp_als(&t, &a).unwrap();
        let db = cp_als(&t, &b).unwrap();

        assert!(da.converged && db.converged);
        assert!((da.relative_residual - db.relative_residual).abs() < 1e-10);
        assert!(measured_relative_residual(&t, &da) < 1e-10);
        assert!(measured_relative_residual(&t, &db) < 1e-10);
    }

    /// Cost of one sweep on a `20 x 20 x 20` rank-3 tensor.
    #[test]
    fn test_cp_cost_20x20x20_rank3() {
        let shape = [20_usize, 20, 20];
        let rank = 3;
        let factors: Vec<Mat> = shape
            .iter()
            .enumerate()
            .map(|(mode, &n)| factor(n, rank, 51 + mode as u64))
            .collect();
        let weights = [1.0, 0.6, 0.3];
        let t = build(&factors, &weights, &shape);

        let mut options = CpOptions::new(rank);
        // This tensor sits in an ALS swamp for ~100 sweeps; 500 is enough for
        // the fit to reach the rounding floor while still measuring per-sweep
        // cost.
        options.max_sweeps = 500;
        options.tol = 0.0;
        options.seed = 0;
        let start = std::time::Instant::now();
        let d = cp_als(&t, &options).unwrap();
        let elapsed = start.elapsed();
        assert_eq!(d.sweeps, 500);
        println!(
            "20x20x20 rank-3: 500 sweeps in {:?} ({:?} per sweep), relative residual {:.3e}",
            elapsed,
            elapsed / 500,
            d.relative_residual
        );
        assert!(d.relative_residual < 1e-3, "{:e}", d.relative_residual);
    }

    /// Invalid configurations are rejected instead of silently doing nothing.
    #[test]
    fn test_cp_rejects_invalid_input() {
        let t = array![[1.0, 2.0], [3.0, 4.0]].into_dyn();
        assert!(cp_als(&t, &CpOptions::new(0)).is_err());
        let mut options = CpOptions::new(1);
        options.max_sweeps = 0;
        assert!(cp_als(&t, &options).is_err());
        let vector = array![1.0, 2.0].into_dyn();
        assert!(cp_als(&vector, &CpOptions::new(1)).is_err());
    }

    /// ALS is a descent method: the fitted residual must fall at every sweep.
    /// A diagonal-only normal matrix (the Hadamard accumulator seeded with the
    /// identity) leaves the off-diagonal component correlations out and makes
    /// the residual grow instead, so this guards the normal equations.
    #[test]
    fn test_cp_residual_decreases_each_sweep() {
        let shape = [10_usize, 10, 10];
        let rank = 3;
        let factors: Vec<Mat> = shape
            .iter()
            .enumerate()
            .map(|(mode, &n)| factor(n, rank, 61 + mode as u64))
            .collect();
        let weights = [1.0, -0.7, 0.4];
        let t = build(&factors, &weights, &shape);

        let mut previous = f64::INFINITY;
        for sweeps in 1..=10 {
            let mut options = CpOptions::new(rank);
            options.max_sweeps = sweeps;
            options.tol = 0.0;
            options.seed = 7;
            let d = cp_als(&t, &options).unwrap();
            println!(
                "sweep {sweeps:2}: relative residual {:.6e}",
                d.relative_residual
            );
            assert!(
                d.relative_residual < previous,
                "sweep {sweeps}: {:e} is not below {:e}",
                d.relative_residual,
                previous
            );
            previous = d.relative_residual;
        }
        assert!(previous < 1e-3, "ten sweeps left {previous:e}");
    }
}
