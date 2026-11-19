//! The discrete Fourier transform (§8.8).
//!
//! `fft` does not normalise and `ifft` divides by `N`, so `ifft(fft(x)) == x`
//! (§8.8). The other conventions are reachable through [`fft_norm`], which is
//! what `fft_norm(x, Norm.Unitary)` needs.
//!
//! The transform itself is `rustfft`'s: it handles every length (including
//! primes) through Bluestein's algorithm, and it was already in the dependency
//! graph.

use num_complex::Complex64 as Complex;
/// Complex representation used by the FFT kernels.
pub use num_complex::Complex64;
use rustfft::FftDirection;

use super::par;
use super::reuse::fft_plan;
use super::reuse::fft_scratch;

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
    if input.is_empty() {
        return input;
    }
    let plan = fft_plan(input.len(), direction(inverse));
    let mut scratch = fft_scratch();
    scratch.resize(plan.get_inplace_scratch_len(), Complex64::default());
    plan.process_with_scratch(&mut input, &mut scratch);
    let scale = scale_for(norm, inverse, input.len());
    if scale < 1.0 {
        for value in &mut input {
            *value *= scale;
        }
    }
    input
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
    for (position, axis) in axes.iter().enumerate() {
        if *axis >= shape.len() {
            return Err(format!(
                "axis {axis} is out of range for a tensor of rank {}",
                shape.len()
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

fn transform_axes(
    mut data: Vec<Complex>,
    shape: &[usize],
    axes: &[usize],
    inverse: bool,
    norm: Norm,
) -> Vec<Complex> {
    if data.is_empty() || axes.is_empty() {
        return data;
    }
    for axis in axes {
        let extent = shape[*axis];
        if extent <= 1 {
            continue;
        }
        let inner: usize = shape[axis + 1..].iter().product();
        transform_axis(&mut data, extent, inner, inverse);
    }
    let count: usize = axes.iter().map(|axis| shape[*axis]).product();
    let scale = scale_for(norm, inverse, count);
    if scale < 1.0 {
        for value in &mut data {
            *value *= scale;
        }
    }
    data
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
    let block_len = extent * inner;
    let plan = fft_plan(extent, direction(inverse));
    let work = |_: usize, data: &mut [Complex64]| {
        let mut scratch = fft_scratch();
        scratch.resize(plan.get_inplace_scratch_len(), Complex64::default());
        if inner == 1 {
            plan.process_with_scratch(data, &mut scratch);
            return;
        }
        let mut lane = vec![Complex64::default(); extent];
        for block in data.chunks_exact_mut(block_len) {
            for offset in 0..inner {
                for (slot, value) in lane
                    .iter_mut()
                    .zip(block.iter().skip(offset).step_by(inner))
                {
                    *slot = *value;
                }
                plan.process_with_scratch(&mut lane, &mut scratch);
                for (slot, value) in block.iter_mut().skip(offset).step_by(inner).zip(&lane) {
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
    check_shape(total, shape, axes)?;
    let rank = shape.len();
    let mut out = Vec::new();
    out.try_reserve_exact(total)
        .map_err(|_| "tensor is too large".to_string())?;
    let mut coordinate = vec![0usize; rank];
    for slot in 0..total {
        let mut rest = slot;
        for axis in (0..rank).rev() {
            coordinate[axis] = rest % shape[axis];
            rest /= shape[axis];
        }
        let mut source = 0usize;
        for axis in 0..rank {
            let extent = shape[axis];
            let mut at = coordinate[axis];
            if axes.contains(&axis) {
                // `fftshift` rolls by n/2 (zero frequency to the middle);
                // `ifftshift` is its inverse.
                let by = if inverse {
                    extent - extent / 2
                } else {
                    extent / 2
                };
                at = if at >= by {
                    at - by
                } else {
                    extent - (by - at)
                };
            }
            source = source * extent + at;
        }
        out.push(source);
    }
    Ok(out)
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
    let spectrum = fft_owned(
        input
            .iter()
            .map(|value| Complex::new(*value, 0.0))
            .collect(),
    );
    spectrum.into_iter().take(input.len() / 2 + 1).collect()
}

/// Rebuilds a real signal of length `n` from the half spectrum.
#[must_use]
pub fn irfft(
    input: &[Complex],
    n: usize,
) -> Vec<f64> {
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
