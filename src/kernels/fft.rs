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
use num_complex::Complex64;
use rustfft::FftPlanner;

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

fn to_internal(input: &[Complex]) -> Vec<Complex64> {
    input.to_vec()
}

const fn from_internal(input: Vec<Complex64>) -> Vec<Complex> {
    input
}

/// The unnormalised forward transform.
#[must_use]
pub fn fft(input: &[Complex]) -> Vec<Complex> {
    if input.is_empty() {
        return Vec::new();
    }
    let mut buffer = to_internal(input);
    let mut planner = FftPlanner::new();
    let transform = planner.plan_fft_forward(buffer.len());
    transform.process(&mut buffer);
    from_internal(buffer)
}

/// The inverse transform, scaled by `1/N`.
#[must_use]
pub fn ifft(input: &[Complex]) -> Vec<Complex> {
    if input.is_empty() {
        return Vec::new();
    }
    let mut buffer = to_internal(input);
    let len = buffer.len();
    let mut planner = FftPlanner::new();
    let transform = planner.plan_fft_inverse(buffer.len());
    transform.process(&mut buffer);
    let scale = 1.0 / len.max(1) as f64;
    from_internal(buffer.into_iter().map(|value| value * scale).collect())
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
    let scale = scale_for(norm, false, input.len());
    fft(input)
        .into_iter()
        .map(|value| Complex::new(value.re * scale, value.im * scale))
        .collect()
}

/// The inverse transform under a chosen normalisation: `Backward` scales by
/// `1/N`, `Forward` not at all, `Unitary` by `1/sqrt(N)`.
#[must_use]
pub fn ifft_norm(
    input: &[Complex],
    norm: Norm,
) -> Vec<Complex> {
    if input.is_empty() {
        return Vec::new();
    }
    let scale = scale_for(norm, true, input.len());
    let mut buffer = to_internal(input);
    let mut planner = FftPlanner::new();
    planner.plan_fft_inverse(buffer.len()).process(&mut buffer);
    from_internal(buffer.into_iter().map(|value| value * scale).collect())
}

fn check_shape(
    len: usize,
    shape: &[usize],
    axes: &[usize],
) -> Result<(), String> {
    if shape.iter().try_fold(1usize, |n, &d| n.checked_mul(d)) != Some(len) {
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
    if input.is_empty() {
        return Ok(Vec::new());
    }
    let mut data = to_internal(input);
    let mut planner = FftPlanner::new();
    for axis in axes {
        let extent = shape[*axis];
        if extent == 0 {
            continue;
        }
        let inner: usize = shape[axis + 1..].iter().product();
        let outer: usize = shape[..*axis].iter().product();
        let transform = if inverse {
            planner.plan_fft_inverse(extent)
        } else {
            planner.plan_fft_forward(extent)
        };
        let mut lane = vec![Complex64::new(0.0, 0.0); extent];
        for block in 0..outer {
            for offset in 0..inner {
                let base = block * extent * inner + offset;
                for (index, slot) in lane.iter_mut().enumerate() {
                    *slot = data[base + index * inner];
                }
                transform.process(&mut lane);
                for (index, slot) in lane.iter().enumerate() {
                    data[base + index * inner] = *slot;
                }
            }
        }
    }
    let count: usize = axes.iter().map(|axis| shape[*axis]).product();
    let scale = scale_for(norm, inverse, count);
    Ok(from_internal(
        data.into_iter().map(|value| value * scale).collect(),
    ))
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
    let complex: Vec<Complex> = input
        .iter()
        .map(|value| Complex::new(*value, 0.0))
        .collect();
    let full = fft_axes(&complex, shape, axes, false, Norm::Backward)?;
    if shape[last] == 0 {
        return Err("real transform axis must be nonempty".to_string());
    }
    let keep = shape[last] / 2 + 1;
    let inner: usize = shape[last + 1..].iter().product();
    let mut out_shape = shape.to_vec();
    out_shape[last] = keep;
    let mut out = Vec::with_capacity(out_shape.iter().product());
    for (index, value) in full.into_iter().enumerate() {
        if (index / inner.max(1)) % shape[last].max(1) < keep {
            out.push(value);
        }
    }
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
    let total = full_shape
        .iter()
        .try_fold(1usize, |n, &d| n.checked_mul(d))
        .ok_or_else(|| "tensor shape overflows usize".to_string())?;
    let mut full = vec![Complex::new(0.0, 0.0); total];
    let rank = shape.len();
    let mut coordinate = vec![0usize; rank];
    for (slot, target) in full.iter_mut().enumerate() {
        let mut rest = slot;
        for axis in (0..rank).rev() {
            coordinate[axis] = rest % full_shape[axis];
            rest /= full_shape[axis];
        }
        // Bins past the half spectrum are the conjugates of their mirror
        // image, mirrored on every transformed axis.
        let mirrored = coordinate[last] >= shape[last];
        let mut source = 0usize;
        for axis in 0..rank {
            let mut at = coordinate[axis];
            if mirrored && axes.contains(&axis) {
                at = (full_shape[axis] - at) % full_shape[axis];
            }
            source = source * shape[axis] + at;
        }
        let value = input[source];
        *target = if mirrored {
            value.conj()
        } else {
            value
        };
    }
    let real = fft_axes(&full, &full_shape, axes, true, Norm::Backward)?;
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
    let total = shape
        .iter()
        .try_fold(1usize, |n, &d| n.checked_mul(d))
        .ok_or_else(|| "tensor shape overflows usize".to_string())?;
    check_shape(total, shape, axes)?;
    let rank = shape.len();
    let mut out = Vec::with_capacity(total);
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
                at = (at + extent - by % extent) % extent;
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
    if n == 0 || d == 0.0 {
        return Vec::new();
    }
    (0..=n / 2)
        .map(|index| index as f64 / (n as f64 * d))
        .collect()
}

/// The transform of a real signal: the non-redundant half of the spectrum,
/// `n/2 + 1` bins.
#[must_use]
pub fn rfft(input: &[f64]) -> Vec<Complex> {
    let spectrum = fft(&input
        .iter()
        .map(|value| Complex::new(*value, 0.0))
        .collect::<Vec<_>>());
    spectrum.into_iter().take(input.len() / 2 + 1).collect()
}

/// Rebuilds a real signal of length `n` from the half spectrum.
#[must_use]
pub fn irfft(
    input: &[Complex],
    n: usize,
) -> Vec<f64> {
    let mut spectrum = vec![Complex::new(0.0, 0.0); n];
    for (index, value) in input.iter().enumerate() {
        if index < n {
            spectrum[index] = *value;
            if index != 0 && index != n - index && n - index < n {
                spectrum[n - index] = value.conj();
            }
        }
    }
    ifft(&spectrum).into_iter().map(|value| value.re).collect()
}

/// The sample frequencies of a length-`n` transform with spacing `d`.
#[must_use]
pub fn fftfreq(
    n: usize,
    d: f64,
) -> Vec<f64> {
    let mut out = Vec::with_capacity(n);
    if n == 0 || d == 0.0 {
        return out;
    }
    let half = n.div_ceil(2);
    for index in 0..n {
        if index < half {
            out.push(index as f64 / (n as f64 * d));
        } else {
            out.push((index as f64 - n as f64) / (n as f64 * d));
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
    shift_permutation(&[input.len()], &[0], false)
        .unwrap_or_default()
        .into_iter()
        .map(|i| input[i])
        .collect()
}
/// Undo a centered one-dimensional spectrum.
#[must_use]
pub fn ifftshift(input: &[Complex]) -> Vec<Complex> {
    shift_permutation(&[input.len()], &[0], true)
        .unwrap_or_default()
        .into_iter()
        .map(|i| input[i])
        .collect()
}
