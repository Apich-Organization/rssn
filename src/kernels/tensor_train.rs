//! Tensor-train compression by successive thin SVDs (TT-SVD).
//!
//! Storage is `sum_k r[k] n[k] r[k+1]`. Reconstruction is optional: entries,
//! inner products and norms can be evaluated directly on the compressed cores.
//! The discarded singular-value norms give an a posteriori Frobenius error
//! estimate; floating-point SVD roundoff is additional to this estimate.

use super::dense::Mat;
use super::dense::svd;
use super::krylov::dot_exact;
use super::reuse::f64_scratch;

/// Invalid input, arithmetic overflow, or failure of an SVD.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TrainError {
    /// Shape, rank limit, tolerance, or core dimensions are invalid.
    InvalidInput,
    /// Input contains a non-finite value.
    NonFinite,
    /// The SVD failed to converge.
    NoConvergence,
    /// A size or floating-point calculation overflowed.
    Overflow,
}

impl std::fmt::Display for TrainError {
    fn fmt(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        f.write_str(match self {
            | Self::InvalidInput => "invalid tensor-train shape, rank or tolerance",
            | Self::NonFinite => "tensor-train input contains non-finite values",
            | Self::NoConvergence => "tensor-train singular-value decomposition did not converge",
            | Self::Overflow => "tensor-train size or arithmetic overflow",
        })
    }
}

impl std::error::Error for TrainError {}

/// A row-major core indexed as `(left * mode + index) * right + right_index`.
#[derive(Clone, Debug)]
pub struct TrainCore {
    left: usize,
    mode: usize,
    right: usize,
    data: Vec<f64>,
}

impl TrainCore {
    /// Construct a finite, nonempty row-major core, validating all extents.
    ///
    /// # Errors
    /// Returns an error for zero/inconsistent extents, non-finite data, or size overflow.
    pub fn new(
        shape: [usize; 3],
        data: Vec<f64>,
    ) -> Result<Self, TrainError> {
        if shape.contains(&0) || size(&shape)? != data.len() {
            return Err(TrainError::InvalidInput);
        }
        if !data.iter().all(|v| v.is_finite()) {
            return Err(TrainError::NonFinite);
        }
        Ok(Self {
            left: shape[0],
            mode: shape[1],
            right: shape[2],
            data,
        })
    }

    /// The two bond dimensions and physical dimension.
    #[must_use]
    pub const fn shape(&self) -> [usize; 3] {
        [self.left, self.mode, self.right]
    }

    /// Row-major core elements.
    #[must_use]
    pub fn data(&self) -> &[f64] {
        &self.data
    }
}

/// A compressed tensor with boundary bond dimensions one.
#[derive(Clone, Debug)]
pub struct TensorTrain {
    shape: Vec<usize>,
    cores: Vec<TrainCore>,
    error: f64,
    tolerance_met: bool,
}

impl TensorTrain {
    /// Import compressed cores without materializing their dense tensor.
    ///
    /// Boundary bonds must be one and neighboring bond dimensions must agree.
    /// A left orthogonalization restores the invariant used by norm and
    /// rounding. The truncation error is zero relative to these supplied
    /// cores; it is not a certificate about their original source data.
    /// Logical tensor size may exceed usize, in which case dense conversion
    /// fails while entry evaluation and compressed contractions remain usable.
    ///
    /// # Errors
    /// Returns an error for empty/mismatched bonds, SVD failure, or arithmetic overflow.
    pub fn from_cores(cores: Vec<TrainCore>) -> Result<Self, TrainError> {
        if cores.first().is_none_or(|c| c.left != 1)
            || cores.last().is_none_or(|c| c.right != 1)
            || cores.windows(2).any(|p| p[0].right != p[1].left)
        {
            return Err(TrainError::InvalidInput);
        }
        let shape = cores.iter().map(|c| c.mode).collect();
        let mut train = Self {
            shape,
            cores,
            error: 0.0,
            tolerance_met: true,
        };
        train.orthogonalize()?;
        Ok(train)
    }

    /// Add two compressed tensors by direct sums of their internal bonds.
    /// Use [`Self::round`] afterwards to reduce the resulting bond ranks.
    ///
    /// # Errors
    /// Returns an error for different tensor shapes, SVD failure, or arithmetic overflow.
    pub fn add(
        &self,
        other: &Self,
    ) -> Result<Self, TrainError> {
        if self.shape != other.shape {
            return Err(TrainError::InvalidInput);
        }
        let order = self.cores.len();
        let mut cores = Vec::with_capacity(order);
        for (k, (a, b)) in self.cores.iter().zip(&other.cores).enumerate() {
            let left = if k == 0 {
                1
            } else {
                a.left.checked_add(b.left).ok_or(TrainError::Overflow)?
            };
            let right = if k + 1 == order {
                1
            } else {
                a.right.checked_add(b.right).ok_or(TrainError::Overflow)?
            };
            let mut data = allocate(&[left, a.mode, right])?;
            for (source, lo, ro) in [
                (a, 0, 0),
                (
                    b,
                    if k == 0 { 0 } else { a.left },
                    if k + 1 == order { 0 } else { a.right },
                ),
            ] {
                for l in 0..source.left {
                    for n in 0..source.mode {
                        for r in 0..source.right {
                            data[((l + lo) * source.mode + n) * right + r + ro] +=
                                source.data[(l * source.mode + n) * source.right + r];
                        }
                    }
                }
            }
            if !data.iter().all(|v| v.is_finite()) {
                return Err(TrainError::Overflow);
            }
            cores.push(TrainCore {
                left,
                mode: a.mode,
                right,
                data,
            });
        }
        let mut result = Self {
            shape: self.shape.clone(),
            cores,
            error: self.error + other.error,
            tolerance_met: self.tolerance_met && other.tolerance_met,
        };
        if !result.error.is_finite() {
            return Err(TrainError::Overflow);
        }
        // restore the left-orthonormal invariant required by norm evaluation
        result.orthogonalize()?;
        Ok(result)
    }

    /// Entrywise product using Kronecker products of the bond dimensions.
    ///
    /// # Errors
    /// Returns an error for different tensor shapes, SVD failure, or arithmetic overflow.
    pub fn hadamard(
        &self,
        other: &Self,
    ) -> Result<Self, TrainError> {
        if self.shape != other.shape {
            return Err(TrainError::InvalidInput);
        }
        let mut cores = Vec::with_capacity(self.cores.len());
        for (a, b) in self.cores.iter().zip(&other.cores) {
            let left = a.left.checked_mul(b.left).ok_or(TrainError::Overflow)?;
            let right = a.right.checked_mul(b.right).ok_or(TrainError::Overflow)?;
            let mut data = allocate(&[left, a.mode, right])?;
            for i in 0..a.left {
                for j in 0..b.left {
                    for n in 0..a.mode {
                        for p in 0..a.right {
                            for q in 0..b.right {
                                data[(((i * b.left + j) * a.mode + n) * right) + p * b.right + q] =
                                    a.data[(i * a.mode + n) * a.right + p]
                                        * b.data[(j * b.mode + n) * b.right + q];
                            }
                        }
                    }
                }
            }
            if !data.iter().all(|v| v.is_finite()) {
                return Err(TrainError::Overflow);
            }
            cores.push(TrainCore {
                left,
                mode: a.mode,
                right,
                data,
            });
        }
        let error = self.error * other.frobenius_norm()
            + other.error * self.frobenius_norm()
            + self.error * other.error;
        if !error.is_finite() {
            return Err(TrainError::Overflow);
        }
        let mut result = Self {
            shape: self.shape.clone(),
            cores,
            error,
            tolerance_met: self.tolerance_met && other.tolerance_met,
        };
        result.orthogonalize()?;
        Ok(result)
    }

    /// Reduce bond dimensions without constructing the dense tensor.
    /// The tolerance controls additional Frobenius error relative to this
    /// represented tensor. The returned error also includes the input's
    /// previous truncation estimate, using the triangle inequality.
    ///
    /// # Errors
    /// Returns an error for invalid tolerance/rank, SVD failure, or arithmetic overflow.
    pub fn round(
        &self,
        relative_tolerance: f64,
        max_rank: Option<usize>,
    ) -> Result<Self, TrainError> {
        if !relative_tolerance.is_finite() || relative_tolerance < 0.0 || max_rank == Some(0) {
            return Err(TrainError::InvalidInput);
        }
        let mut result = self.clone();
        result.orthogonalize()?;
        let budget = relative_tolerance * result.frobenius_norm();
        let local = budget / ((result.cores.len() - 1).max(1) as f64).sqrt();
        let mut error = 0.0_f64;
        for k in (1..result.cores.len()).rev() {
            let core = &result.cores[k];
            let cols = core
                .mode
                .checked_mul(core.right)
                .ok_or(TrainError::Overflow)?;
            let d = svd(&Mat {
                rows: core.left,
                cols,
                data: core.data.clone(),
            })
            .map_err(|_| TrainError::NoConvergence)?;
            let mut tails = vec![0.0_f64; d.s.len() + 1];
            for i in (0..d.s.len()).rev() {
                tails[i] = tails[i + 1].hypot(d.s[i]);
            }
            let rank = (1..=d.s.len())
                .find(|&r| tails[r] <= local)
                .unwrap_or(d.s.len())
                .min(max_rank.unwrap_or(usize::MAX));
            error = error.hypot(tails[rank]);
            let mut data = allocate(&[rank, cols])?;
            for i in 0..rank {
                for j in 0..cols {
                    data[i * cols + j] = d.v.at(j, i);
                }
            }
            let mode = core.mode;
            let right = core.right;
            let previous = &result.cores[k - 1];
            let mut preceding = allocate(&[previous.left, previous.mode, rank])?;
            for i in 0..previous.left * previous.mode {
                for r in 0..rank {
                    for j in 0..previous.right {
                        preceding[i * rank + r] +=
                            previous.data[i * previous.right + j] * d.u.at(j, r) * d.s[r];
                    }
                }
            }
            if !preceding.iter().all(|v| v.is_finite()) {
                return Err(TrainError::Overflow);
            }
            result.cores[k] = TrainCore {
                left: rank,
                mode,
                right,
                data,
            };
            result.cores[k - 1].right = rank;
            result.cores[k - 1].data = preceding;
        }
        result.error += error;
        if !result.error.is_finite() {
            return Err(TrainError::Overflow);
        }
        result.tolerance_met = error <= budget;
        result.orthogonalize()?;
        Ok(result)
    }

    fn orthogonalize(&mut self) -> Result<(), TrainError> {
        for k in 0..self.cores.len() - 1 {
            let a = &self.cores[k];
            let rows = a.left.checked_mul(a.mode).ok_or(TrainError::Overflow)?;
            let d = svd(&Mat {
                rows,
                cols: a.right,
                data: a.data.clone(),
            })
            .map_err(|_| TrainError::NoConvergence)?;
            let rank = d.s.len();
            let mut data = allocate(&[rows, rank])?;
            for i in 0..rows {
                for j in 0..rank {
                    data[i * rank + j] = d.u.at(i, j);
                }
            }
            let b = &self.cores[k + 1];
            let width = b.mode.checked_mul(b.right).ok_or(TrainError::Overflow)?;
            let mut next = allocate(&[rank, width])?;
            for i in 0..rank {
                for j in 0..width {
                    for r in 0..a.right {
                        next[i * width + j] += d.s[i] * d.v.at(r, i) * b.data[r * width + j];
                    }
                }
            }
            if !next.iter().all(|v| v.is_finite()) {
                return Err(TrainError::Overflow);
            }
            self.cores[k].right = rank;
            self.cores[k].data = data;
            self.cores[k + 1].left = rank;
            self.cores[k + 1].data = next;
        }
        Ok(())
    }

    /// Original tensor shape.
    #[must_use]
    pub fn shape(&self) -> &[usize] {
        &self.shape
    }

    /// Compressed cores.
    #[must_use]
    pub fn cores(&self) -> &[TrainCore] {
        &self.cores
    }

    /// Frobenius norm estimate of discarded singular values.
    #[must_use]
    pub const fn truncation_error(&self) -> f64 {
        self.error
    }

    /// Whether the requested relative tolerance was met despite a rank cap.
    #[must_use]
    pub const fn tolerance_met(&self) -> bool {
        self.tolerance_met
    }

    /// Number of stored scalar elements.
    #[must_use]
    pub fn storage_len(&self) -> usize {
        self.cores.iter().map(|c| c.data.len()).sum()
    }

    /// Frobenius norm, evaluated with scaled sum-of-squares arithmetic.
    /// TT-SVD makes every preceding core left-orthonormal, so the final
    /// core alone has the same norm as the reconstructed tensor.
    #[must_use]
    pub fn frobenius_norm(&self) -> f64 {
        self.cores.last().map_or(0.0, |core| {
            core.data.iter().fold(0.0_f64, |s, &v| s.hypot(v))
        })
    }

    /// Evaluate one element without constructing the dense tensor.
    ///
    /// # Errors
    /// Returns an error for invalid indices or an unrepresentable contraction.
    pub fn entry(
        &self,
        indices: &[usize],
    ) -> Result<f64, TrainError> {
        if indices.len() != self.shape.len() || indices.iter().zip(&self.shape).any(|(i, n)| i >= n)
        {
            return Err(TrainError::InvalidInput);
        }
        let mut state = vec![1.0];
        for (core, &index) in self.cores.iter().zip(indices) {
            let mut next = vec![0.0; core.right];
            for (left, &value) in state.iter().enumerate() {
                for (right, out) in next.iter_mut().enumerate() {
                    *out += value * core.data[(left * core.mode + index) * core.right + right];
                }
            }
            state = next;
        }
        let value = state[0];
        if !value.is_finite() {
            return Err(TrainError::Overflow);
        }
        Ok(value)
    }

    /// Reconstruct the dense, row-major tensor. Allocates prod(shape) elements.
    ///
    /// # Errors
    /// Returns an error for size/allocation overflow or unrepresentable arithmetic.
    pub fn to_dense(&self) -> Result<Vec<f64>, TrainError> {
        self.to_dense_with_limit(usize::MAX)
    }

    /// Dense conversion with an explicit limit on the output element count.
    /// Use this when loading cores supplied by another process or when the
    /// logical tensor may be exponentially larger than its compressed storage.
    ///
    /// # Errors
    /// Returns an error when size exceeds the limit, allocation fails, or arithmetic overflows.
    pub fn to_dense_with_limit(
        &self,
        max_elements: usize,
    ) -> Result<Vec<f64>, TrainError> {
        if size(&self.shape)? > max_elements {
            return Err(TrainError::Overflow);
        }
        let mut result = vec![1.0];
        let mut prefixes = 1usize;
        for core in &self.cores {
            let next_prefixes = prefixes
                .checked_mul(core.mode)
                .ok_or(TrainError::Overflow)?;
            let len = next_prefixes
                .checked_mul(core.right)
                .ok_or(TrainError::Overflow)?;
            let mut next = Vec::new();
            next.try_reserve_exact(len)
                .map_err(|_| TrainError::Overflow)?;
            next.resize(len, 0.0);
            for p in 0..prefixes {
                for index in 0..core.mode {
                    for left in 0..core.left {
                        for right in 0..core.right {
                            next[(p * core.mode + index) * core.right + right] += result
                                [p * core.left + left]
                                * core.data[(left * core.mode + index) * core.right + right];
                        }
                    }
                }
            }
            if !next.iter().all(|v| v.is_finite()) {
                return Err(TrainError::Overflow);
            }
            result = next;
            prefixes = next_prefixes;
        }
        Ok(result)
    }

    /// Inner product by contraction of bond matrices, without decompression.
    ///
    /// Each physical slice advances the bond matrix `E` by a term `A_k^T E B_k`.
    /// The contraction order is selected from the bond ranks, and the same two bond
    /// buffers are reused across the cores.
    ///
    /// # Errors
    /// Returns an error for different shapes or an unrepresentable contraction.
    pub fn inner(
        &self,
        other: &Self,
    ) -> Result<f64, TrainError> {
        const STACK_ENTRIES: usize = 64;
        let layout = self.inner_layout(other)?;
        if layout.len == 0 {
            self.contract_inner(other, layout, &mut [])
        } else if layout.len <= STACK_ENTRIES {
            self.contract_inner(other, layout, &mut [0.0; STACK_ENTRIES])
        } else {
            let mut work = f64_scratch(layout.len);
            work.resize(layout.len, 0.0);
            self.contract_inner(other, layout, &mut work)
        }
    }

    /// Number of `f64` entries required by [`Self::inner_with_workspace`].
    ///
    /// The workspace is sized by the bond ranks. Two rank-one tensor trains need no workspace.
    ///
    /// # Errors
    /// Returns an error if the shapes differ or the workspace extent is unrepresentable.
    pub fn inner_workspace_len(
        &self,
        other: &Self,
    ) -> Result<usize, TrainError> {
        Ok(self.inner_layout(other)?.len)
    }

    /// Contracts two tensor trains using a caller-owned workspace.
    ///
    /// Provide at least [`Self::inner_workspace_len`] entries. The workspace can be reused
    /// with different trains, and any entries beyond the required length are left untouched.
    ///
    /// # Errors
    /// Returns an error for different shapes, insufficient workspace, or arithmetic overflow.
    pub fn inner_with_workspace(
        &self,
        other: &Self,
        work: &mut [f64],
    ) -> Result<f64, TrainError> {
        let layout = self.inner_layout(other)?;
        if work.len() < layout.len {
            return Err(TrainError::InvalidInput);
        }
        self.contract_inner(other, layout, work)
    }

    fn inner_layout(
        &self,
        other: &Self,
    ) -> Result<InnerLayout, TrainError> {
        const MAX_WORK_ENTRIES: usize = isize::MAX as usize / std::mem::size_of::<f64>();
        if self.shape != other.shape {
            return Err(TrainError::InvalidInput);
        }
        if self
            .cores
            .iter()
            .chain(&other.cores)
            .all(|core| core.right == 1)
        {
            return Ok(InnerLayout { bond: 0, len: 0 });
        }
        let mut bond = 1;
        let mut intermediate = 0;
        for (a, b) in self.cores.iter().zip(&other.cores) {
            bond = bond.max(size(&[a.right, b.right])?);
            let (_, required) = inner_step(a, b);
            intermediate = intermediate.max(required);
        }
        let len = size(&[2, bond])?
            .checked_add(intermediate)
            .ok_or(TrainError::Overflow)?;
        if len > MAX_WORK_ENTRIES {
            return Err(TrainError::Overflow);
        }
        Ok(InnerLayout { bond, len })
    }

    #[expect(
        clippy::arithmetic_side_effects,
        clippy::indexing_slicing,
        reason = "the layout checks bond sizes and the workspace before contraction"
    )]
    fn contract_inner(
        &self,
        other: &Self,
        layout: InnerLayout,
        work: &mut [f64],
    ) -> Result<f64, TrainError> {
        if let ([a], [b]) = (self.cores.as_slice(), other.cores.as_slice()) {
            let pairs = a.data.iter().copied().zip(b.data.iter().copied());
            let fast: f64 = pairs.clone().map(|(x, y)| x * y).sum();
            let value = if fast.is_normal() {
                fast
            } else {
                dot_exact(pairs)
            };
            return value
                .is_finite()
                .then_some(value)
                .ok_or(TrainError::Overflow);
        }
        if layout.len == 0 {
            let mut state = 1.0;
            for (a, b) in self.cores.iter().zip(&other.cores) {
                state = a
                    .data
                    .iter()
                    .zip(&b.data)
                    .map(|(x, y)| state * x * y)
                    .sum::<f64>();
                if !state.is_finite() {
                    return Err(TrainError::Overflow);
                }
            }
            return Ok(state);
        }

        let (mut state, rest) = work.split_at_mut(layout.bond);
        let (mut next, aux) = rest.split_at_mut(layout.bond);
        state[0] = 1.0;
        for (a, b) in self.cores.iter().zip(&other.cores) {
            let input = &state[..a.left * b.left];
            let output = &mut next[..a.right * b.right];
            output.fill(0.0);
            let (step, required) = inner_step(a, b);
            match step {
                | InnerStep::Direct => inner_direct(a, b, input, output),
                | InnerStep::Left => inner_left(a, b, input, output, &mut aux[..required]),
                | InnerStep::Right => inner_right(a, b, input, output, &mut aux[..required]),
            }
            if output.iter().any(|v| !v.is_finite()) {
                if step == InnerStep::Direct {
                    return Err(TrainError::Overflow);
                }
                // The factored intermediate can overflow before multiplication by the other core.
                output.fill(0.0);
                inner_direct(a, b, input, output);
                if output.iter().any(|v| !v.is_finite()) {
                    return Err(TrainError::Overflow);
                }
            }
            std::mem::swap(&mut state, &mut next);
        }
        Ok(state[0])
    }
}

#[derive(Clone, Copy)]
struct InnerLayout {
    bond: usize,
    len: usize,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum InnerStep {
    Direct,
    Left,
    Right,
}

fn inner_step(
    a: &TrainCore,
    b: &TrainCore,
) -> (InnerStep, usize) {
    let mut best = (InnerStep::Direct, 0);
    let mut cost = a
        .left
        .saturating_mul(b.left)
        .saturating_mul(a.right)
        .saturating_mul(b.right);
    for (step, extent, terms) in [
        (
            InnerStep::Left,
            a.right.checked_mul(b.left),
            a.left.saturating_add(b.right),
        ),
        (
            InnerStep::Right,
            a.left.checked_mul(b.right),
            b.left.saturating_add(a.right),
        ),
    ] {
        if let Some(extent) = extent {
            let candidate = extent.saturating_mul(terms);
            if (candidate, extent) < (cost, best.1) {
                best = (step, extent);
                cost = candidate;
            }
        }
    }
    best
}

#[expect(
    clippy::arithmetic_side_effects,
    clippy::indexing_slicing,
    reason = "core construction validates the storage layout and the callers bound the indices"
)]
fn inner_row(
    core: &TrainCore,
    left: usize,
    mode: usize,
) -> &[f64] {
    let start = (left * core.mode + mode) * core.right;
    &core.data[start..start + core.right]
}

#[inline]
fn inner_axpy(
    out: &mut [f64],
    coefficient: f64,
    row: &[f64],
) {
    for (value, rhs) in out.iter_mut().zip(row) {
        *value += coefficient * rhs;
    }
}

fn inner_direct(
    a: &TrainCore,
    b: &TrainCore,
    state: &[f64],
    next: &mut [f64],
) {
    for (i, s) in state.chunks_exact(b.left).enumerate() {
        for (j, &weight) in s.iter().enumerate() {
            for k in 0..a.mode {
                let right = inner_row(b, j, k);
                for (&left, out) in inner_row(a, i, k)
                    .iter()
                    .zip(next.chunks_exact_mut(b.right))
                {
                    inner_axpy(out, weight * left, right);
                }
            }
        }
    }
}

// Accumulate (A_k^T E) B_k over the physical index k.
#[expect(
    clippy::indexing_slicing,
    clippy::single_call_fn,
    reason = "the left contraction kernel reads column j within its b.left columns"
)]
fn inner_left(
    a: &TrainCore,
    b: &TrainCore,
    state: &[f64],
    next: &mut [f64],
    aux: &mut [f64],
) {
    for k in 0..a.mode {
        aux.fill(0.0);
        for (i, s) in state.chunks_exact(b.left).enumerate() {
            for (&left, t) in inner_row(a, i, k).iter().zip(aux.chunks_exact_mut(b.left)) {
                inner_axpy(t, left, s);
            }
        }
        for j in 0..b.left {
            let right = inner_row(b, j, k);
            for (t, out) in aux.chunks_exact(b.left).zip(next.chunks_exact_mut(b.right)) {
                inner_axpy(out, t[j], right);
            }
        }
    }
}

// Accumulate A_k^T (E B_k) with the opposite order of contraction.
#[expect(
    clippy::single_call_fn,
    reason = "keeps the right-first kernel separate from dispatch"
)]
fn inner_right(
    a: &TrainCore,
    b: &TrainCore,
    state: &[f64],
    next: &mut [f64],
    aux: &mut [f64],
) {
    for k in 0..a.mode {
        aux.fill(0.0);
        for (s, t) in state
            .chunks_exact(b.left)
            .zip(aux.chunks_exact_mut(b.right))
        {
            for (j, &weight) in s.iter().enumerate() {
                inner_axpy(t, weight, inner_row(b, j, k));
            }
        }
        for (i, t) in aux.chunks_exact(b.right).enumerate() {
            for (&left, out) in inner_row(a, i, k)
                .iter()
                .zip(next.chunks_exact_mut(b.right))
            {
                inner_axpy(out, left, t);
            }
        }
    }
}

fn size(shape: &[usize]) -> Result<usize, TrainError> {
    shape
        .iter()
        .try_fold(1usize, |n, &m| n.checked_mul(m).ok_or(TrainError::Overflow))
}

fn allocate(shape: &[usize]) -> Result<Vec<f64>, TrainError> {
    let len = size(shape)?;
    let mut data = Vec::new();
    data.try_reserve_exact(len)
        .map_err(|_| TrainError::Overflow)?;
    data.resize(len, 0.0);
    Ok(data)
}

/// Compress a finite, nonempty dense tensor with TT-SVD.
///
/// Each split is assigned relative tolerance / sqrt(order - 1). A rank cap
/// may prevent satisfying that tolerance; inspect `tolerance_met` and the
/// returned truncation estimate. Zero tensors use bond rank one. A one-mode
/// tensor is stored exactly. Input is normalized before decomposition to
/// avoid squaring the input's scale in intermediate norm computations.
///
/// # Errors
/// Returns an error for invalid shape/options, non-finite input, SVD failure,
/// allocation/size overflow or unrepresentable rescaling.
pub fn tt_svd(
    data: &[f64],
    shape: &[usize],
    relative_tolerance: f64,
    max_rank: Option<usize>,
) -> Result<TensorTrain, TrainError> {
    if shape.is_empty()
        || shape.contains(&0)
        || size(shape)? != data.len()
        || !relative_tolerance.is_finite()
        || relative_tolerance < 0.0
        || max_rank == Some(0)
    {
        return Err(TrainError::InvalidInput);
    }
    if !data.iter().all(|v| v.is_finite()) {
        return Err(TrainError::NonFinite);
    }
    let scale = data.iter().fold(0.0_f64, |s, v| s.max(v.abs()));
    let scale = if scale == 0.0 { 1.0 } else { scale };
    let mut remainder: Vec<f64> = data.iter().map(|x| x / scale).collect();
    let norm = remainder.iter().fold(0.0_f64, |s, &v| s.hypot(v));
    let budget = relative_tolerance * norm;
    let local_budget = budget / ((shape.len() - 1).max(1) as f64).sqrt();
    let mut error = 0.0_f64;
    let mut left = 1usize;
    let mut cores = Vec::with_capacity(shape.len());
    for mode in 0..shape.len() - 1 {
        let rows = left.checked_mul(shape[mode]).ok_or(TrainError::Overflow)?;
        let cols = remainder.len() / rows;
        let decomposition = svd(&Mat {
            rows,
            cols,
            data: remainder,
        })
        .map_err(|_| TrainError::NoConvergence)?;
        let mut tails = vec![0.0_f64; decomposition.s.len() + 1];
        for i in (0..decomposition.s.len()).rev() {
            tails[i] = tails[i + 1].hypot(decomposition.s[i]);
        }
        let desired = (1..=decomposition.s.len())
            .find(|&r| tails[r] <= local_budget)
            .unwrap_or(decomposition.s.len());
        let right = desired.min(max_rank.unwrap_or(usize::MAX));
        error = error.hypot(tails[right]);
        let mut core = Vec::with_capacity(rows * right);
        for i in 0..rows {
            for j in 0..right {
                core.push(decomposition.u.at(i, j));
            }
        }
        cores.push(TrainCore {
            left,
            mode: shape[mode],
            right,
            data: core,
        });
        remainder = vec![0.0; right * cols];
        for i in 0..right {
            for j in 0..cols {
                remainder[i * cols + j] = decomposition.s[i] * decomposition.v.at(j, i);
            }
        }
        left = right;
    }
    for value in &mut remainder {
        *value *= scale;
    }
    if !remainder.iter().all(|x| x.is_finite()) || !(error * scale).is_finite() {
        return Err(TrainError::Overflow);
    }
    cores.push(TrainCore {
        left,
        mode: shape.last().copied().ok_or(TrainError::InvalidInput)?,
        right: 1,
        data: remainder,
    });
    Ok(TensorTrain {
        shape: shape.to_vec(),
        cores,
        error: error * scale,
        tolerance_met: error <= budget,
    })
}
