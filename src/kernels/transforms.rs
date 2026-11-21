//! # Numerical integral transforms
//!
//! The discrete Fourier transform uses the convention
//! `X_k = sum_j x_j exp(-2 pi i j k / n)`; its inverse divides by `n`.
//! [`fft`] and [`ifft`] act in place. Power-of-two lengths use the scalar
//! radix-2 implementation for compatibility with earlier RSSN results;
//! other lengths share the runtime-selected plans in [`super::fft`].
//! [`fft_many_par`] and [`ifft_many_par`] schedule independent rows and
//! preserve the results and order of the individual transforms.
//!
//! [`rfft`] and [`irfft`] return or accept only the `n/2 + 1` nonredundant
//! bins. Even lengths pack two real samples into each complex sample and
//! use a half-length transform. Odd lengths use one full complex transform.

use std::alloc::Layout;
use std::f64::consts::PI;
use std::sync::Arc;

use num_complex::Complex;
use num_complex::Complex64;

use super::fft::Norm;
use super::fft::transform_in_place;
use super::reuse::PlanCache;
use super::reuse::fft_scratch_for;

/// The optimized, in-place, iterative Cooley-Tukey Fast Fourier Transform (FFT) algorithm.
///
/// This function operates directly on the slice, avoiding memory allocations and copies.
/// It assumes the input slice length is a power of two; other lengths are
/// handled by the public wrappers using shared FFT plans.
///
/// # Arguments
/// * `data` - A mutable slice of `Complex<f64>` representing the input sequence.
fn fft_cooley_tukey_in_place(data: &mut [Complex<f64>]) {
    let n = data.len();

    if n <= 1 {
        return;
    }

    let mut j = 0;

    for i in 0..n {
        if i < j {
            data.swap(i, j);
        }

        let mut bit = n >> 1;

        while j & bit != 0 {
            j ^= bit;

            bit >>= 1;
        }

        j ^= bit;
    }

    let mut len = 2;

    while len <= n {
        let half_len = len / 2;

        let w_m = Complex::from_polar(1.0, -2.0 * PI / (len as f64));

        for i in (0..n).step_by(len) {
            let mut w = Complex::new(1.0, 0.0);

            for k in 0..half_len {
                let even_idx = i + k;

                let odd_idx = i + k + half_len;

                let t = w * data[odd_idx];

                let u = data[even_idx];

                data[even_idx] = u + t;

                data[odd_idx] = u - t;

                w *= w_m;
            }
        }

        len *= 2;
    }
}

/// The inverse radix-2 transform: `conj(fft(conj(x))) / n`, in place.
///
/// The input length must be a power of two.
fn inverse_pow2_in_place(data: &mut [Complex<f64>]) {
    for val in data.iter_mut() {
        *val = val.conj();
    }

    fft_cooley_tukey_in_place(data);

    let n_f64 = data.len() as f64;

    for val in data.iter_mut() {
        *val = val.conj() / n_f64;
    }
}

/// Computes the Fast Fourier Transform (FFT) of a sequence of complex numbers.
///
/// Works at any length and leaves the length unchanged: a power-of-two `n` uses
/// the scalar radix-2 kernel. Other lengths use a shared, runtime-selected
/// FFT plan. Both paths compute the unpadded `n`-point DFT.
///
/// # Arguments
/// * `data` - A mutable `Vec<Complex<f64>>` representing the input sequence.
///
/// # Example
/// ```rust
/// use num_complex::Complex;
/// use rssn::kernels::transforms::fft;
///
/// let mut data = vec![
///     Complex::new(1.0, 0.0),
///     Complex::new(1.0, 0.0),
///     Complex::new(0.0, 0.0),
///     Complex::new(0.0, 0.0),
/// ];
///
/// fft(&mut data);
///
/// assert!((data[0].re - 2.0).abs() < 1e-9);
/// ```
///
/// A length that is not a power of two keeps its length:
/// ```rust
/// use num_complex::Complex;
/// use rssn::kernels::transforms::fft;
///
/// let mut data = vec![Complex::new(1.0, 0.0); 3];
///
/// fft(&mut data);
///
/// assert_eq!(data.len(), 3);
/// // [3, 0, 0] for a constant signal.
/// assert!((data[0].re - 3.0).abs() < 1e-12);
/// ```
pub fn fft(data: &mut [Complex<f64>]) {
    if data.len().is_power_of_two() {
        fft_cooley_tukey_in_place(data);
    } else {
        transform_in_place(data, false, Norm::Backward);
    }
}

/// Computes the Inverse Fast Fourier Transform (IFFT).
///
/// Works at any length and leaves the length unchanged, so `ifft(fft(x))`
/// returns `x` at every length. The `1/n` normalisation is included.
///
/// # Arguments
/// * `data` - A mutable `Vec<Complex<f64>>` representing the input frequency-domain sequence.
///
/// # Example
/// ```rust
/// use num_complex::Complex;
/// use rssn::kernels::transforms::fft;
/// use rssn::kernels::transforms::ifft;
///
/// let mut data = vec![Complex::new(1.0, 2.0), Complex::new(3.0, 4.0)];
///
/// let original = data.clone();
///
/// fft(&mut data);
///
/// ifft(&mut data);
///
/// assert!((data[0].re - original[0].re).abs() < 1e-9);
///
/// assert!((data[0].im - original[0].im).abs() < 1e-9);
/// ```
pub fn ifft(data: &mut [Complex<f64>]) {
    if data.len() <= 1 {
        return;
    }
    if data.len().is_power_of_two() {
        inverse_pow2_in_place(data);
    } else {
        transform_in_place(data, true, Norm::Backward);
    }
}

/// Computes the real-input FFT: the `n/2 + 1` non-redundant bins of the DFT of
/// a real sequence.
///
/// For real `x`, the full spectrum obeys `X_{n-k} = conj(X_k)`, so only the
/// bins `k = 0 ..= n/2` are returned (length `n/2 + 1`, unnormalised like
/// [`fft`]). At even `n = 2m` the sequence is packed as the `m` complex values
/// `x_{2j} + i x_{2j+1}`, so one `m`-point complex transform replaces the
/// `n`-point one. Odd `n` has no such packing and runs a full `n`-point complex
/// transform, keeping the first `n/2 + 1` bins.
///
/// # Arguments
/// * `x` - The real input sequence.
///
/// # Example
/// ```rust
/// use rssn::kernels::transforms::rfft;
///
/// // An impulse has a flat half spectrum.
/// let spectrum = rfft(&[1.0, 0.0, 0.0, 0.0]);
///
/// assert_eq!(spectrum.len(), 3);
/// assert!((spectrum[0].re - 1.0).abs() < 1e-12);
/// assert!((spectrum[2].re - 1.0).abs() < 1e-12);
/// ```
///
/// A constant signal concentrates all the energy in the DC bin:
/// ```rust
/// use rssn::kernels::transforms::rfft;
///
/// let spectrum = rfft(&[2.0; 5]);
///
/// assert!((spectrum[0].re - 10.0).abs() < 1e-12);
/// assert!(spectrum[1..].iter().all(|bin| bin.norm() < 1e-12));
/// ```
#[must_use]
pub fn rfft(x: &[f64]) -> Vec<Complex64> {
    let n = x.len();

    if n == 0 {
        return Vec::new();
    }

    if n == 1 {
        return vec![Complex64::new(x[0], 0.0)];
    }

    if n.is_multiple_of(2) {
        let m = n / 2;
        let mut packed = fft_scratch_for(m);
        packed.extend(
            x.as_chunks::<2>()
                .0
                .iter()
                .map(|pair| Complex64::new(pair[0], pair[1])),
        );
        transform_in_place(&mut packed, false, Norm::Backward);

        let twiddles = real_twiddles(m);
        let mut spectrum = Vec::with_capacity(m + 1);
        let dc = packed[0];
        spectrum.push(Complex64::new(dc.re + dc.im, 0.0));
        for (k, twiddle) in (1..m).zip(twiddles.iter()) {
            let a = packed[k] * 0.5;
            let b = packed[m - k].conj() * 0.5;
            let even = a + b;
            let odd = Complex64::new(0.0, -1.0) * (a - b);
            spectrum.push(even + twiddle * odd);
        }
        spectrum.push(Complex64::new(dc.re - dc.im, 0.0));
        spectrum
    } else {
        let mut full: Vec<_> = x.iter().map(|&value| Complex64::new(value, 0.0)).collect();
        transform_in_place(&mut full, false, Norm::Backward);
        full.truncate(n / 2 + 1);
        full
    }
}

/// Inverts [`rfft`]: reconstructs the real signal of length `n` from its
/// `n/2 + 1` non-redundant bins, with the `1/n` normalisation of [`ifft`].
///
/// The imaginary parts of the DC bin (`k = 0`) and, for even `n`, of the
/// Nyquist bin (`k = n/2`) are ignored, because a real signal has purely real
/// values there. Even `n = 2m` is unpacked into one `m`-point inverse complex
/// transform; odd `n` reconstructs the full conjugate-symmetric spectrum and
/// takes one `n`-point inverse transform.
///
/// Returns an empty vector when `n` is zero or when `spectrum.len()` is not
/// `n/2 + 1`. Extreme finite magnitudes are rescaled in the FFT workspace
/// to protect the intermediate sums and the final normalisation.
///
/// For repeated transforms of the same length, see [`IrfftPlan`].
///
/// # Arguments
/// * `spectrum` - The half spectrum produced by [`rfft`] (or any conjugate-
///   symmetric data).
/// * `n` - The length of the reconstructed real signal.
///
/// # Example
/// ```rust
/// use rssn::kernels::transforms::irfft;
/// use rssn::kernels::transforms::rfft;
///
/// let x = [0.5, -1.25, 2.0, 0.0, -3.5, 1.0];
/// let back = irfft(&rfft(&x), x.len());
///
/// for (got, want) in back.iter().zip(x) {
///     assert!((got - want).abs() < 1e-12);
/// }
/// ```
#[must_use]
pub fn irfft(
    spectrum: &[Complex64],
    n: usize,
) -> Vec<f64> {
    if n == 0 || spectrum.len() != n / 2 + 1 {
        return Vec::new();
    }
    let Ok(plan) = IrfftPlan::new(n) else {
        return Vec::new();
    };
    let mut output = vec![0.0; n];
    plan.write_samples(spectrum, &mut output);
    output
}

/// An invalid length or slice shape for a real inverse FFT.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IrfftError {
    /// The length cannot be represented by the required buffers.
    Length,
    /// The input bins or output samples disagree with the plan's length.
    Shape,
}

/// Plans repeated inverse real FFTs of a fixed sample length.
///
/// Retains the even-length unpacking factors and uses the shared FFT and
/// scratch pools. Each call to [`Self::apply_into`] writes to a caller-supplied
/// output slice with the normalisation and bin conventions of [`irfft`].
#[derive(Debug, Clone)]
pub struct IrfftPlan {
    n: usize,
    twiddles: Option<Arc<[Complex64]>>,
}

impl IrfftPlan {
    /// Prepares an inverse transform whose output has `n` real samples.
    ///
    /// # Errors
    /// Returns [`IrfftError::Length`] if the output or FFT workspace size
    /// exceeds the maximum size of a slice.
    pub fn new(n: usize) -> Result<Self, IrfftError> {
        let even = n.is_multiple_of(2);
        let fft_len = if even { n / 2 } else { n };
        if Layout::array::<f64>(n).is_err() || Layout::array::<Complex64>(fft_len).is_err() {
            return Err(IrfftError::Length);
        }
        let twiddles = (even && n >= 4).then(|| real_twiddles(fft_len));
        Ok(Self { n, twiddles })
    }

    /// Real sample count of the inverse transform.
    #[must_use]
    pub const fn output_len(&self) -> usize {
        self.n
    }

    /// Number of nonredundant bins, or zero for an empty transform.
    #[must_use]
    pub const fn spectrum_len(&self) -> usize {
        if self.n == 0 {
            0
        } else {
            self.n / 2 + 1
        }
    }

    /// Writes one inverse real FFT to the supplied sample slice.
    ///
    /// # Errors
    /// Returns [`IrfftError::Shape`] and leaves the output untouched unless
    /// both slices match the plan.
    pub fn apply_into(
        &self,
        spectrum: &[Complex64],
        output: &mut [f64],
    ) -> Result<(), IrfftError> {
        if spectrum.len() != self.spectrum_len() || output.len() != self.n {
            return Err(IrfftError::Shape);
        }
        self.write_samples(spectrum, output);
        Ok(())
    }

    fn write_samples(
        &self,
        spectrum: &[Complex64],
        output: &mut [f64],
    ) {
        let n = self.n;
        match n {
            | 0 => return,
            | 1 => {
                output[0] = spectrum[0].re;
                return;
            },
            | 2 => {
                let (dc, nyquist) = (spectrum[0].re, spectrum[1].re);
                output[0] = f64::midpoint(dc, nyquist);
                output[1] = f64::midpoint(dc, -nyquist);
                return;
            },
            | _ => {},
        }
        if let Some(input_scale) = inverse_real_scale(spectrum, n) {
            let (dc, nyquist, interior) = split_real_spectrum(spectrum, n);
            // A constant real spectrum is an exact impulse. Avoid rounding
            // its leading sample past f64::MAX when restoring the scale.
            if dc.to_bits() == nyquist.unwrap_or(dc).to_bits()
                && interior
                    .iter()
                    .all(|z| z.re.to_bits() == dc.to_bits() && z.im == 0.0)
            {
                output.fill(0.0);
                output[0] = dc;
                return;
            }
            self.write_transform::<true>(spectrum, output, input_scale);
        } else {
            self.write_transform::<false>(spectrum, output, 1.0);
        }
    }

    fn write_transform<const RESCALE: bool>(
        &self,
        spectrum: &[Complex64],
        output: &mut [f64],
        input_scale: f64,
    ) {
        let n = self.n;
        let (input_scale, output_scale) = if RESCALE {
            (input_scale, input_scale.recip())
        } else {
            (1.0, 1.0)
        };
        if let Some(twiddles) = &self.twiddles {
            let m = n / 2;
            let mut packed = fft_scratch_for(m);
            let half_scale = 0.5 * input_scale;
            let dc = spectrum[0].re * half_scale;
            let nyquist = spectrum[m].re * half_scale;
            packed.push(Complex64::new(dc + nyquist, dc - nyquist));
            for (k, twiddle) in (1..m).zip(twiddles.iter()) {
                let a = spectrum[k] * half_scale;
                let b = spectrum[m - k].conj() * half_scale;
                let even = a + b;
                let odd = (a - b) * twiddle.conj();
                packed.push(even + Complex64::new(0.0, 1.0) * odd);
            }
            transform_in_place(&mut packed, true, Norm::Backward);
            for (pair, value) in output.as_chunks_mut::<2>().0.iter_mut().zip(packed.iter()) {
                *pair = [value.re * output_scale, value.im * output_scale];
            }
        } else {
            let half = n / 2;
            let mut full = fft_scratch_for(n);
            full.resize(n, Complex64::ZERO);
            full[0] = Complex64::new(spectrum[0].re * input_scale, 0.0);
            for k in 1..=half {
                let value = spectrum[k] * input_scale;
                full[k] = value;
                full[n - k] = value.conj();
            }
            transform_in_place(&mut full, true, Norm::Backward);
            for (sample, value) in output.iter_mut().zip(full.iter()) {
                *sample = value.re * output_scale;
            }
        }
    }
}

fn split_real_spectrum(
    spectrum: &[Complex64],
    n: usize,
) -> (f64, Option<f64>, &[Complex64]) {
    let dc = spectrum[0].re;
    if n.is_multiple_of(2) {
        let end = spectrum.len() - 1;
        (dc, Some(spectrum[end].re), &spectrum[1..end])
    } else {
        (dc, None, &spectrum[1..])
    }
}

#[allow(clippy::needless_bitwise_bool)] // Eager boolean reductions allow vectorisation.
fn inverse_real_scale(
    spectrum: &[Complex64],
    n: usize,
) -> Option<f64> {
    // The imaginary parts of self-conjugate bins never enter the real FFT.
    let (dc, nyquist, interior) = split_real_spectrum(spectrum, n);
    let edge_a = dc.abs();
    let edge_b = nyquist.unwrap_or(0.0).abs();
    let length = n as f64;
    // Leave headroom for packing and for FFT algorithms that use an
    // unnormalised convolution internally.
    let upper = f64::MAX / (16.0 * length * length);
    let lower = 4.0 * length * f64::MIN_POSITIVE;
    let mut large = (edge_a > upper) | (edge_b > upper);
    let mut small = (edge_a < lower) & (edge_b < lower);
    // Independent boolean reductions keep this screening pass vectorisable.
    for bin in interior {
        let re = bin.re.abs();
        let im = bin.im.abs();
        large |= (re > upper) | (im > upper);
        small &= (re < lower) & (im < lower);
    }
    if large {
        // Absolute float bit patterns are ordered by magnitude. Non-finite
        // bins retain the ordinary FFT's propagation behavior.
        let mut peak = edge_a.to_bits().max(edge_b.to_bits());
        for bin in interior {
            peak = peak.max(bin.re.abs().to_bits()).max(bin.im.abs().to_bits());
        }
        if !f64::from_bits(peak).is_finite() {
            return None;
        }
        let shift = (peak >> 52) - (upper.to_bits() >> 52) + 1;
        Some(f64::from_bits((1023 - shift) << 52))
    } else if small
        && (edge_a > 0.0 || edge_b > 0.0 || interior.iter().any(|z| z.re != 0.0 || z.im != 0.0))
    {
        // Even the smallest subnormal becomes normal at this working scale.
        Some(f64::from_bits((1023 + 256) << 52))
    } else {
        None
    }
}

/// Forward twiddles for the interior bins of an even real FFT.
fn real_twiddles(half_len: usize) -> Arc<[Complex64]> {
    static TWIDDLES: PlanCache<usize, [Complex64]> = PlanCache::new(16);
    let build = || {
        let scale = -PI / half_len as f64;
        (1..half_len)
            .map(|bin| Complex64::from_polar(1.0, scale * bin as f64))
            .collect()
    };
    if half_len <= 32_768 {
        TWIDDLES.get_or_insert_with(half_len, build)
    } else {
        build()
    }
}

/// Computes the Forward FFT for a slice.
///
/// **Use this for efficient, parallel chunk processing (e.g., with rayon).**
/// Works at any length and operates in place: a power-of-two length uses the
/// scalar radix-2 kernel and other lengths use the shared FFT planner.
/// Bit-identical to [`fft`].
///
/// # Arguments
/// * `data` - A mutable slice of `Complex<f64>` representing the input sequence.
///
/// # Example
/// ```rust
/// use num_complex::Complex;
/// use rssn::kernels::transforms::fft_slice;
///
/// let mut data = vec![
///     Complex::new(1.0, 0.0),
///     Complex::new(1.0, 0.0),
///     Complex::new(0.0, 0.0),
///     Complex::new(0.0, 0.0),
/// ];
///
/// fft_slice(&mut data);
///
/// assert!((data[0].re - 2.0).abs() < 1e-9);
/// ```
pub fn fft_slice(data: &mut [Complex<f64>]) {
    fft(data);
}

/// Computes the Inverse FFT for a slice.
///
/// **Use this for efficient, parallel chunk processing (e.g., with rayon).**
/// Works at any length and operates in place; the `1/n` normalisation is
/// included. Bit-identical to [`ifft`].
///
/// # Arguments
/// * `data` - A mutable slice of `Complex<f64>` representing the input frequency-domain sequence.
///
/// # Example
/// ```rust
/// use num_complex::Complex;
/// use rssn::kernels::transforms::fft_slice;
/// use rssn::kernels::transforms::ifft_slice;
///
/// let mut data = vec![Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)];
///
/// let original = data.clone();
///
/// fft_slice(&mut data);
///
/// ifft_slice(&mut data);
///
/// assert!((data[0].re - original[0].re).abs() < 1e-9);
/// ```
pub fn ifft_slice(data: &mut [Complex<f64>]) {
    ifft(data);
}

/// Forward FFT of every row, concurrently, in row order.
///
/// Each row runs the serial [`fft`] kernel and only different rows run
/// concurrently, so the result is bit-identical to transforming the rows one at
/// a time and therefore independent of the thread count. Fewer than `grain`
/// rows (or `grain < 2`) run sequentially, so a small batch does not pay for the
/// task spawn. Use the batch entry point instead of a hand-written `rayon` map
/// when the rows are independent (the usual row/column pass of a 2-D transform).
#[must_use]
pub fn fft_many_par(
    rows: &[Vec<Complex<f64>>],
    grain: usize,
) -> Vec<Vec<Complex<f64>>> {
    crate::kernels::par::map_slice(rows, grain, |row| {
        let mut buffer = row.clone();
        fft(&mut buffer);
        buffer
    })
}

/// [`fft_many_par`] for the inverse transform.
#[must_use]
pub fn ifft_many_par(
    rows: &[Vec<Complex<f64>>],
    grain: usize,
) -> Vec<Vec<Complex<f64>>> {
    crate::kernels::par::map_slice(rows, grain, |row| {
        let mut buffer = row.clone();
        ifft(&mut buffer);
        buffer
    })
}
