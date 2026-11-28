//! One- and N-dimensional Fourier transforms and frequency indexing.
//!
//! The discrete Fourier transform (§8.8).
//!
//! `fft` does not normalise and `ifft` divides by `N`, so `ifft(fft(x)) == x`
//! (§8.8). The other conventions are reachable through [`fft_norm`], which is
//! what `fft_norm(x, Norm.Unitary)` needs.
//!
//! `RustFFT` chooses a plan for the transform length and the available CPU
//! instructions. Calls of matching lengths can reuse plans and scratch.

use num_complex::Complex64 as Complex;
/// Complex representation used by the FFT kernels.
pub use num_complex::Complex64;
use rustfft::Fft;
use rustfft::FftDirection;

use super::par;
use super::reuse::fft_plan;
use super::reuse::fft_scratch_for;

/// The normalisation convention of [`fft_norm`] (§8.8).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Norm {
    /// No scaling on the forward transform, `1/N` on the inverse (the default).
    Backward,
    /// `1/N` on the forward transform, no scaling on the inverse.
    Forward,
    /// `1/sqrt(N)` on both.
    Unitary,
}

/// The unnormalised forward transform.
#[must_use]
pub fn fft(input: &[Complex]) -> Vec<Complex> {
    fft_owned(input.to_vec())
}

pub(crate) fn fft_owned(input: Vec<Complex>) -> Vec<Complex> {
    transform(input, false, Norm::Backward)
}

/// The inverse transform, scaled by `1/N`.
#[must_use]
pub fn ifft(input: &[Complex]) -> Vec<Complex> {
    ifft_owned(input.to_vec())
}

pub(crate) fn ifft_owned(input: Vec<Complex>) -> Vec<Complex> {
    transform(input, true, Norm::Backward)
}

fn transform(
    mut input: Vec<Complex>,
    inverse: bool,
    norm: Norm,
) -> Vec<Complex> {
    transform_in_place(&mut input, inverse, norm);
    input
}

/// Applies a one-dimensional FFT to an existing complex buffer.
///
/// Uses the direction and normalization conventions of the allocating FFT
/// functions. Empty and singleton inputs are left unchanged.
pub fn transform_in_place(
    data: &mut [Complex],
    inverse: bool,
    norm: Norm,
) {
    if data.len() <= 1 {
        return;
    }
    let plan = fft_plan(data.len(), direction(inverse));
    let scratch_len = plan.get_inplace_scratch_len();
    let mut scratch = fft_scratch_for(scratch_len);
    scratch.resize(scratch_len, Complex64::ZERO);
    plan.process_with_scratch(data, &mut scratch);
    scale_transform(data, norm, inverse, data.len());
}

fn scale_transform(
    data: &mut [Complex],
    norm: Norm,
    inverse: bool,
    count: usize,
) {
    let scale = scale_for(norm, inverse, count);
    if scale < 1.0 {
        super::par::for_chunks_mut(data, 4_096, 4, |_, chunk| {
            super::simd::complex_scale_in_place(scale, chunk);
        });
    }
}

const fn direction(inverse: bool) -> FftDirection {
    if inverse {
        FftDirection::Inverse
    } else {
        FftDirection::Forward
    }
}

/// The scale a transform over `count` points carries under `norm`.
fn scale_for(
    norm: Norm,
    inverse: bool,
    count: usize,
) -> f64 {
    let count = count.max(1) as f64;
    match (norm, inverse) {
        | (Norm::Unitary, _) => 1.0 / count.sqrt(),
        | (Norm::Backward, false) | (Norm::Forward, true) => 1.0,
        | (Norm::Backward, true) | (Norm::Forward, false) => 1.0 / count,
    }
}

/// The forward transform under a chosen normalisation.
///
/// `Backward` is the plain forward transform, `Forward` scales it by `1/N`
/// (so the inverse needs none) and `Unitary` by `1/sqrt(N)`.
#[must_use]
pub fn fft_norm(
    input: &[Complex],
    norm: Norm,
) -> Vec<Complex> {
    transform(input.to_vec(), false, norm)
}

/// The inverse transform under a chosen normalisation: `Backward` scales by
/// `1/N`, `Forward` not at all, `Unitary` by `1/sqrt(N)`.
#[must_use]
pub fn ifft_norm(
    input: &[Complex],
    norm: Norm,
) -> Vec<Complex> {
    transform(input.to_vec(), true, norm)
}

fn shape_len(shape: &[usize]) -> Option<usize> {
    if shape.contains(&0) {
        return Some(0);
    }
    shape.iter().try_fold(1_usize, |n, &d| n.checked_mul(d))
}

fn check_shape(
    len: usize,
    shape: &[usize],
    axes: &[usize],
) -> Result<(), String> {
    if shape_len(shape) != Some(len) {
        return Err(format!("the shape {shape:?} does not match {len} elements"));
    }
    check_axes(shape.len(), axes)
}

fn check_axes(
    rank: usize,
    axes: &[usize],
) -> Result<(), String> {
    for (position, axis) in axes.iter().enumerate() {
        if *axis >= rank {
            return Err(format!(
                "axis {axis} is out of range for a tensor of rank {rank}"
            ));
        }
        if axes[..position].contains(axis) {
            return Err(format!("axis {axis} is listed twice"));
        }
    }
    Ok(())
}

/// The transform over `axes` of a row-major tensor of the given `shape`,
/// computed one axis at a time (row-column decomposition) with the 1-D
/// transform above. An empty `axes` is the identity.
///
/// # Errors
/// An axis out of range or repeated, or a `shape` that disagrees with the
/// number of elements.
pub fn fft_axes(
    input: &[Complex],
    shape: &[usize],
    axes: &[usize],
    inverse: bool,
    norm: Norm,
) -> Result<Vec<Complex>, String> {
    check_shape(input.len(), shape, axes)?;
    Ok(transform_axes(input.to_vec(), shape, axes, inverse, norm))
}

/// Applies the selected-direction FFT to an existing row-major complex buffer.
///
/// # Errors
/// As for an allocating axis transform. Validation completes before the
/// buffer is modified.
pub fn fft_axes_in_place(
    data: &mut [Complex],
    shape: &[usize],
    axes: &[usize],
    inverse: bool,
    norm: Norm,
) -> Result<(), String> {
    check_shape(data.len(), shape, axes)?;
    transform_axes_in_place(data, shape, axes, inverse, norm);
    Ok(())
}

fn transform_axes(
    mut data: Vec<Complex>,
    shape: &[usize],
    axes: &[usize],
    inverse: bool,
    norm: Norm,
) -> Vec<Complex> {
    transform_axes_in_place(&mut data, shape, axes, inverse, norm);
    data
}

/// Applies an axis transform to a buffer whose shape and axes are validated.
pub(crate) fn transform_axes_in_place(
    data: &mut [Complex],
    shape: &[usize],
    axes: &[usize],
    inverse: bool,
    norm: Norm,
) {
    if data.is_empty() || axes.is_empty() {
        return;
    }
    for axis in axes {
        let extent = shape[*axis];
        if extent <= 1 {
            continue;
        }
        let inner: usize = shape[axis + 1..].iter().product();
        transform_axis(data, extent, inner, inverse);
    }
    let count: usize = axes.iter().map(|axis| shape[*axis]).product();
    scale_transform(data, norm, inverse, count);
}

/// Independent lanes are batched into disjoint blocks of the tensor.
fn transform_axis(
    data: &mut [Complex64],
    extent: usize,
    inner: usize,
    inverse: bool,
) {
    // The shape has already been checked. Each chunk is a whole number of
    // outer blocks and therefore contains only complete transforms.
    if extent <= 1 {
        return;
    }
    let block_len = extent * inner;
    let plan = fft_plan(extent, direction(inverse));
    // A single wide strided block has many independent lanes, but the outer
    // block scheduler below cannot split it. Pack its lanes contiguously to
    // give workers independent FFT batches and avoid repeated strided reads.
    if inner >= 16
        && block_len >= 262_144
        && data.len() / block_len < 4
        && crate::execution::should_split(inner, 16)
    {
        let mut transposed = fft_scratch_for(block_len);
        transposed.resize(block_len, Complex64::ZERO);
        for block in data.chunks_exact_mut(block_len) {
            transform_wide_block(block, &mut transposed, extent, inner, plan.as_ref());
        }
        return;
    }
    let work = |_: usize, data: &mut [Complex64]| {
        let scratch_len = plan.get_inplace_scratch_len();
        let mut scratch = fft_scratch_for(scratch_len);
        scratch.resize(scratch_len, Complex64::ZERO);
        if inner == 1 {
            plan.process_with_scratch(data, &mut scratch);
            return;
        }
        let mut lane = fft_scratch_for(extent);
        lane.resize(extent, Complex64::ZERO);
        for block in data.chunks_exact_mut(block_len) {
            for offset in 0..inner {
                for (slot, value) in lane
                    .iter_mut()
                    .zip(block.iter().skip(offset).step_by(inner))
                {
                    *slot = *value;
                }
                plan.process_with_scratch(&mut lane, &mut scratch);
                for (slot, value) in block
                    .iter_mut()
                    .skip(offset)
                    .step_by(inner)
                    .zip(lane.iter())
                {
                    *slot = *value;
                }
            }
        }
    };
    let chunk_len = block_len * (32_768 / block_len).max(1);
    if data.len() >= 262_144 {
        par::for_chunks_mut(data, chunk_len, 2, work);
    } else {
        work(0, data);
    }
}

fn transform_wide_block(
    block: &mut [Complex64],
    transposed: &mut [Complex64],
    extent: usize,
    inner: usize,
    plan: &dyn Fft<f64>,
) {
    let lanes_per_batch = (32_768 / extent).clamp(8, inner.div_ceil(2));
    let lane_chunk = lanes_per_batch * extent;
    par::for_chunks_mut(transposed, lane_chunk, 2, |chunk, lanes| {
        let first = chunk * lanes_per_batch;
        // Work on a tile of adjacent columns so each source row is read
        // contiguously and a tile remains resident across the gather pass.
        for (row, source) in block.chunks_exact(inner).enumerate() {
            for (value, lane) in source[first..].iter().zip(lanes.chunks_exact_mut(extent)) {
                lane[row] = *value;
            }
        }
        let scratch_len = plan.get_inplace_scratch_len();
        let mut scratch = fft_scratch_for(scratch_len);
        scratch.resize(scratch_len, Complex64::ZERO);
        plan.process_with_scratch(lanes, &mut scratch);
    });

    let rows_per_batch = (32_768 / inner).max(1);
    par::for_chunks_mut(block, rows_per_batch * inner, 2, |chunk, rows| {
        let first = chunk * rows_per_batch;
        for (column, lane) in transposed.chunks_exact(extent).enumerate() {
            for (row, value) in rows.chunks_exact_mut(inner).zip(&lane[first..]) {
                row[column] = *value;
            }
        }
    });
}

/// The half-spectrum transform of a real tensor over `axes`: the last listed
/// axis keeps `n/2 + 1` bins. Returns the data and its shape.
///
/// # Errors
/// As [`fft_axes`], and an empty `axes`.
pub fn rfft_axes(
    input: &[f64],
    shape: &[usize],
    axes: &[usize],
) -> Result<(Vec<Complex>, Vec<usize>), String> {
    let Some(last) = axes.last().copied() else {
        return Err("a real transform needs at least one axis".to_string());
    };
    check_shape(input.len(), shape, axes)?;
    if shape[last] == 0 {
        return Err("real transform axis must be nonempty".to_string());
    }
    let keep = shape[last] / 2 + 1;
    let mut out_shape = shape.to_vec();
    out_shape[last] = keep;
    if input.is_empty() {
        return Ok((Vec::new(), out_shape));
    }
    if input.len() == shape[last] {
        // All other extents are one. Reuse the packed real transform and
        // its cached twiddles for a single lane.
        return Ok((rfft(input), out_shape));
    }
    let mut complex: Vec<_> = input
        .iter()
        .map(|&value| Complex::new(value, 0.0))
        .collect();
    // Hermitian symmetry holds along this axis while all other axes still
    // contain real-valued inputs. Transform the remaining axes on half as
    // many lanes by truncating before processing them.
    let inner: usize = shape[last + 1..].iter().product();
    transform_axis(&mut complex, shape[last], inner, false);
    let mut out = Vec::with_capacity(out_shape.iter().product());
    for block in complex.chunks_exact(shape[last] * inner) {
        out.extend_from_slice(&block[..keep * inner]);
    }
    let out = transform_axes(
        out,
        &out_shape,
        &axes[..axes.len() - 1],
        false,
        Norm::Backward,
    );
    Ok((out, out_shape))
}

/// Rebuilds a real tensor from the half spectrum of [`rfft_axes`]. `shape` is
/// the half-spectrum shape; `n_last` the length of the last listed axis in the
/// result (`2 * (m - 1)` when the original was even).
///
/// # Errors
/// As [`fft_axes`], and an `n_last` the half spectrum cannot belong to.
pub fn irfft_axes(
    input: &[Complex],
    shape: &[usize],
    axes: &[usize],
    n_last: usize,
) -> Result<(Vec<f64>, Vec<usize>), String> {
    let Some(last) = axes.last().copied() else {
        return Err("a real transform needs at least one axis".to_string());
    };
    check_shape(input.len(), shape, axes)?;
    if n_last == 0 || n_last / 2 + 1 != shape[last] {
        return Err(format!(
            "a half spectrum of {} bins cannot come from a length-{n_last} axis",
            shape[last]
        ));
    }
    let mut full_shape = shape.to_vec();
    full_shape[last] = n_last;
    let total = shape_len(&full_shape).ok_or_else(|| "tensor shape overflows usize".to_string())?;
    if total == 0 {
        return Ok((Vec::new(), full_shape));
    }
    if total == n_last {
        let plan = super::transforms::IrfftPlan::new(n_last)
            .map_err(|_| "tensor is too large".to_string())?;
        let mut samples = Vec::new();
        samples
            .try_reserve_exact(n_last)
            .map_err(|_| "tensor is too large".to_string())?;
        samples.resize(n_last, 0.0);
        plan.apply_into(input, &mut samples)
            .map_err(|_| "tensor shape disagrees with inverse FFT plan".to_string())?;
        return Ok((samples, full_shape));
    }
    let mut full = Vec::new();
    full.try_reserve_exact(total)
        .map_err(|_| "tensor is too large".to_string())?;
    // Undo the other axes before expanding the half spectrum. The inverse
    // across those axes reduces the mirror relation to one dimension.
    let half = transform_axes(
        input.to_vec(),
        shape,
        &axes[..axes.len() - 1],
        true,
        Norm::Backward,
    );
    let inner = shape[last + 1..].iter().product::<usize>();
    let keep = shape[last];
    for block in half.chunks_exact(keep * inner) {
        full.extend_from_slice(block);
        for bin in keep..n_last {
            let mirrored = (n_last - bin) * inner;
            full.extend(
                block[mirrored..mirrored + inner]
                    .iter()
                    .map(Complex64::conj),
            );
        }
    }
    let real = transform_axes(full, &full_shape, &[last], true, Norm::Backward);
    Ok((real.into_iter().map(|value| value.re).collect(), full_shape))
}

/// Where each element goes under `fftshift` (`inverse = false`) or
/// `ifftshift` of a row-major tensor over `axes`: result element `i` is input
/// element `permutation[i]`.
///
/// # Errors
/// An axis out of range or repeated.
pub fn shift_permutation(
    shape: &[usize],
    axes: &[usize],
    inverse: bool,
) -> Result<Vec<usize>, String> {
    let total = shape_len(shape).ok_or_else(|| "tensor shape overflows usize".to_string())?;
    check_axes(shape.len(), axes)?;
    let mut out = Vec::new();
    out.try_reserve_exact(total)
        .map_err(|_| "tensor is too large".to_string())?;
    out.extend(0..total);
    roll_axes(&mut out, shape, axes, inverse);
    Ok(out)
}

/// Shifts zero frequency along the selected axes of a row-major buffer.
///
/// Forward shifts roll each axis to the right by half its length, rounded
/// down. Inverse shifts roll left by the same amount.
///
/// # Errors
/// Returns an error for a shape that disagrees with the buffer, or for
/// repeated or out-of-range axes. The buffer is left unchanged on error.
pub fn shift_axes_in_place<T: Send>(
    data: &mut [T],
    shape: &[usize],
    axes: &[usize],
    inverse: bool,
) -> Result<(), String> {
    check_shape(data.len(), shape, axes)?;
    roll_axes(data, shape, axes, inverse);
    Ok(())
}

fn roll_axes<T: Send>(
    data: &mut [T],
    shape: &[usize],
    axes: &[usize],
    inverse: bool,
) {
    if data.is_empty() || axes.is_empty() {
        return;
    }
    let mut inner = 1;
    for (axis, &extent) in shape.iter().enumerate().rev() {
        let block_len = inner * extent;
        if extent > 1 && axes.contains(&axis) {
            let offset = (extent / 2) * inner;
            let roll = |_: usize, blocks: &mut [T]| {
                for block in blocks.chunks_exact_mut(block_len) {
                    if inverse {
                        block.rotate_left(offset);
                    } else {
                        block.rotate_right(offset);
                    }
                }
            };
            // Group neighboring blocks to avoid scheduling tiny rotations.
            if data.len() >= 262_144 {
                let chunk_len = block_len * (32_768 / block_len).max(1);
                par::for_chunks_mut(data, chunk_len, 4, roll);
            } else {
                roll(0, data);
            }
        }
        inner = block_len;
    }
}

/// The non-negative sample frequencies of a length-`n` real transform,
/// `n/2 + 1` of them.
#[must_use]
pub fn rfftfreq(
    n: usize,
    d: f64,
) -> Vec<f64> {
    let mut out = Vec::new();
    if n == 0 || d == 0.0 || out.try_reserve_exact(n / 2 + 1).is_err() {
        return out;
    }
    out.extend((0..=n / 2).map(|index| (index as f64 / n as f64) / d));
    out
}

/// The transform of a real signal: the non-redundant half of the spectrum,
/// `n/2 + 1` bins.
#[must_use]
pub fn rfft(input: &[f64]) -> Vec<Complex> {
    super::transforms::rfft(input)
}

/// Rebuilds a real signal of length `n` from the half spectrum.
#[must_use]
pub fn irfft(
    input: &[Complex],
    n: usize,
) -> Vec<f64> {
    if n == 0 || input.len() == n / 2 + 1 {
        return super::transforms::irfft(input, n);
    }
    // Retain this entry point's explicit zero-fill and overwrite behavior for
    // incomplete or overlong spectra.
    let mut spectrum = Vec::new();
    if spectrum.try_reserve_exact(n).is_err() {
        return Vec::new();
    }
    spectrum.resize(n, Complex64::default());
    for (index, value) in input.iter().enumerate() {
        if index < n {
            spectrum[index] = *value;
            if index != 0 && index != n - index && n - index < n {
                spectrum[n - index] = value.conj();
            }
        }
    }
    ifft_owned(spectrum)
        .into_iter()
        .map(|value| value.re)
        .collect()
}

/// The sample frequencies of a length-`n` transform with spacing `d`.
#[must_use]
pub fn fftfreq(
    n: usize,
    d: f64,
) -> Vec<f64> {
    let mut out = Vec::new();
    if n == 0 || d == 0.0 || out.try_reserve_exact(n).is_err() {
        return out;
    }
    let half = n.div_ceil(2);
    for index in 0..n {
        if index < half {
            out.push((index as f64 / n as f64) / d);
        } else {
            out.push(-((n - index) as f64 / n as f64) / d);
        }
    }
    out
}

/// Forward two-dimensional FFT of a row-major matrix.
///
/// # Errors
/// Returns an error if the shape overflows or does not match the input length.
pub fn fft2(
    input: &[Complex],
    rows: usize,
    cols: usize,
) -> Result<Vec<Complex>, String> {
    fft_axes(input, &[rows, cols], &[0, 1], false, Norm::Backward)
}
/// Inverse two-dimensional FFT of a row-major matrix.
///
/// # Errors
/// Returns an error if the shape overflows or does not match the input length.
pub fn ifft2(
    input: &[Complex],
    rows: usize,
    cols: usize,
) -> Result<Vec<Complex>, String> {
    fft_axes(input, &[rows, cols], &[0, 1], true, Norm::Backward)
}
/// Forward FFT over every tensor axis.
///
/// # Errors
/// Returns an error if the shape overflows or does not match the input length.
pub fn fftn(
    input: &[Complex],
    shape: &[usize],
) -> Result<Vec<Complex>, String> {
    fft_axes(
        input,
        shape,
        &(0..shape.len()).collect::<Vec<_>>(),
        false,
        Norm::Backward,
    )
}
/// Center the zero frequency in a one-dimensional spectrum.
#[must_use]
pub fn fftshift(input: &[Complex]) -> Vec<Complex> {
    let mut output = input.to_vec();
    output.rotate_right(input.len() / 2);
    output
}
/// Undo a centered one-dimensional spectrum.
#[must_use]
pub fn ifftshift(input: &[Complex]) -> Vec<Complex> {
    let mut output = input.to_vec();
    output.rotate_left(input.len() / 2);
    output
}

/// Forward two-dimensional real-to-complex FFT of a row-major real matrix.
///
/// Output has shape `(rows, cols / 2 + 1)`.
///
/// # Errors
/// Returns an error if the shape overflows or does not match the input length.
pub fn rfft2(
    input: &[f64],
    rows: usize,
    cols: usize,
) -> Result<(Vec<Complex>, Vec<usize>), String> {
    rfft_axes(input, &[rows, cols], &[0, 1])
}

/// Inverse two-dimensional complex-to-real FFT returning a `(rows, cols)` real matrix.
///
/// # Errors
/// Returns an error if the input shape is inconsistent with `(rows, cols / 2 + 1)`.
pub fn irfft2(
    input: &[Complex],
    rows: usize,
    cols: usize,
) -> Result<Vec<f64>, String> {
    irfft_axes(input, &[rows, cols / 2 + 1], &[0, 1], cols).map(|(data, _)| data)
}

/// Forward N-dimensional real-to-complex FFT over every tensor axis.
///
/// # Errors
/// Returns an error if `shape` is empty or does not match `input.len()`.
pub fn rfftn(
    input: &[f64],
    shape: &[usize],
) -> Result<(Vec<Complex>, Vec<usize>), String> {
    let axes: Vec<usize> = (0..shape.len()).collect();
    rfft_axes(input, shape, &axes)
}

/// Inverse N-dimensional complex-to-real FFT over every tensor axis.
///
/// `out_last` is the length of the last axis in the real output tensor.
///
/// # Errors
/// Returns an error if `shape` is empty or does not match `input.len()`.
pub fn irfftn(
    input: &[Complex],
    shape: &[usize],
    out_last: usize,
) -> Result<Vec<f64>, String> {
    let axes: Vec<usize> = (0..shape.len()).collect();
    irfft_axes(input, shape, &axes, out_last).map(|(data, _)| data)
}

fn shift_nd_inner(
    input: &[Complex],
    shape: &[usize],
    inverse: bool,
) -> Result<Vec<Complex>, String> {
    check_shape(input.len(), shape, &[])?;
    let mut out = input.to_vec();
    roll_axes(
        &mut out,
        shape,
        &(0..shape.len()).collect::<Vec<_>>(),
        inverse,
    );
    Ok(out)
}

/// Centers the zero-frequency component along every axis of an N-D spectrum.
///
/// # Errors
/// Returns an error if `shape` overflows or does not match `input.len()`.
pub fn fftshift_nd(
    input: &[Complex],
    shape: &[usize],
) -> Result<Vec<Complex>, String> {
    shift_nd_inner(input, shape, false)
}

/// Undoes [`fftshift_nd`] along every axis of an N-D spectrum.
///
/// # Errors
/// Returns an error if `shape` overflows or does not match `input.len()`.
pub fn ifftshift_nd(
    input: &[Complex],
    shape: &[usize],
) -> Result<Vec<Complex>, String> {
    shift_nd_inner(input, shape, true)
}

#[cfg(test)]
mod nd_tests {
    use super::*;

    #[test]
    fn shifts_agree_on_odd_even_and_partial_axes() {
        let input: Vec<_> = (0..12).map(|i| Complex::new(f64::from(i), -1.0)).collect();
        let expected = [10, 11, 6, 7, 8, 9, 4, 5, 0, 1, 2, 3];
        let shape = [2, 3, 2];
        let axes = [1, 0];
        assert_eq!(shift_permutation(&shape, &axes, false).unwrap(), expected);
        let mut shifted = input.clone();
        shift_axes_in_place(&mut shifted, &shape, &axes, false).unwrap();
        assert_eq!(shifted, expected.map(|i| input[i]));
        shift_axes_in_place(&mut shifted, &shape, &axes, true).unwrap();
        assert_eq!(shifted, input);
        assert_eq!(
            fftshift_nd(&input, &shape).unwrap(),
            [11, 10, 7, 6, 9, 8, 5, 4, 1, 0, 3, 2].map(|i| input[i])
        );
    }

    #[test]
    fn shift_handles_empty_and_invalid_shapes_before_mutating() {
        let mut values = [1, 2, 3, 4];
        assert!(shift_axes_in_place(&mut values, &[2, 2], &[0, 0], false).is_err());
        assert!(shift_axes_in_place(&mut values, &[2, 2], &[2], false).is_err());
        assert!(shift_axes_in_place(&mut values, &[3, 2], &[1], false).is_err());
        assert_eq!(values, [1, 2, 3, 4]);
        assert!(shift_permutation(&[usize::MAX, 2], &[1], false).is_err());
        assert_eq!(
            shift_permutation(&[usize::MAX, 2, 0], &[0], false).unwrap(),
            Vec::<usize>::new()
        );
        assert_eq!(
            fftshift_nd(&[], &[usize::MAX, 2, 0]).unwrap(),
            Vec::<Complex>::new()
        );
    }

    #[test]
    fn large_parallel_shift_matches_row_and_plane_offsets() {
        let shape = [4, 257, 256];
        let [planes, rows, cols] = shape;
        let mut actual = (0..planes * rows * cols).collect::<Vec<_>>();
        rayon::ThreadPoolBuilder::new()
            .num_threads(2)
            .build()
            .unwrap()
            .install(|| shift_axes_in_place(&mut actual, &shape, &[0, 2, 1], false).unwrap());
        for (destination, source) in actual.into_iter().enumerate() {
            let col = destination % cols;
            let row = destination / cols % rows;
            let plane = destination / (cols * rows);
            let expected = (((plane + planes.div_ceil(2)) % planes * rows
                + (row + rows.div_ceil(2)) % rows)
                * cols)
                + (col + cols.div_ceil(2)) % cols;
            assert_eq!(source, expected);
        }
    }

    #[test]
    fn in_place_axes_checks_inputs_and_honors_normalization() {
        let shape = [2, 3];
        let mut samples = [Complex::ZERO; 6];
        samples[0] = Complex::new(2.0, 3.0);
        let original = samples;
        assert!(fft_axes_in_place(&mut samples, &shape, &[0, 2], false, Norm::Forward).is_err());
        assert_eq!(samples, original);

        fft_axes_in_place(&mut samples, &shape, &[1, 0], false, Norm::Forward).unwrap();
        for bin in samples {
            assert!((bin - original[0] / 6.0).norm() < 1e-12);
        }
        fft_axes_in_place(&mut samples, &shape, &[0, 1], true, Norm::Forward).unwrap();
        for (sample, expected) in samples.into_iter().zip(original) {
            assert!((sample - expected).norm() < 1e-12);
        }
    }

    #[test]
    fn single_real_lane_matches_complex_reference_for_odd_and_even_lengths() {
        for len in [1, 2, 7, 16] {
            let shape = [1, len, 1];
            let signal = (0..len).map(|i| (i as f64 * 0.3).sin()).collect::<Vec<_>>();
            let embedded = signal
                .iter()
                .map(|&x| Complex::new(x, 0.0))
                .collect::<Vec<_>>();
            let reference = fft(&embedded);
            let (bins, half_shape) = rfft_axes(&signal, &shape, &[2, 0, 1]).unwrap();
            assert_eq!(half_shape, [1, len / 2 + 1, 1]);
            assert_eq!(bins.len(), len / 2 + 1);
            for (bin, expected) in bins.iter().zip(reference) {
                assert!((*bin - expected).norm() < 1e-12);
            }
            let (recovered, recovered_shape) =
                irfft_axes(&bins, &half_shape, &[2, 0, 1], len).unwrap();
            assert_eq!(recovered_shape, shape);
            for (&sample, expected) in recovered.iter().zip(signal) {
                assert!((sample - expected).abs() < 1e-12);
            }
        }

        let (samples, shape) =
            irfft_axes(&[Complex::new(f64::MAX, 0.0); 3], &[3], &[0], 4).unwrap();
        assert_eq!(shape, [4]);
        assert_eq!(samples, [f64::MAX, 0.0, 0.0, 0.0]);
    }

    #[test]
    fn rfft2_and_fftshift_nd_roundtrip() {
        let rows = 5;
        let cols = 6;
        let signal: Vec<f64> = (0..rows * cols)
            .map(|i| ((i as f64) * 0.31).sin() + 0.25 * ((i as f64) * 0.73).cos())
            .collect();
        let (spec, shape) = rfft2(&signal, rows, cols).unwrap();
        assert_eq!(shape, vec![rows, cols / 2 + 1]);
        let rec = irfft2(&spec, rows, cols).unwrap();
        for (a, b) in signal.iter().zip(&rec) {
            assert!((a - b).abs() < 1e-12);
        }

        let shifted = fftshift_nd(&spec, &shape).unwrap();
        let unshifted = ifftshift_nd(&shifted, &shape).unwrap();
        assert_eq!(spec, unshifted);
    }
}
