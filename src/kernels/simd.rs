//! Runtime-dispatched SIMD BLAS-1 and complex-vector microkernels.
//!
//! Uses [`crate::execution::host_environment`] (probed once at startup and
//! cached in a `OnceLock`) to select AVX2+FMA kernels on `x86_64`, NEON kernels
//! on `aarch64`, and falls back to 4-way unrolled scalar loops on other targets.
//! When compiled with `-C target-feature=+avx2,+fma` or `+neon`, dispatch folds
//! to a compile-time constant with zero runtime overhead.
//! Every reduction uses a fixed lane tree so results are deterministic across
//! runs on the same host.

#![allow(unsafe_code)]

use num_complex::Complex64;

#[inline]
#[cfg(target_arch = "x86_64")]
fn has_avx2_fma() -> bool {
    cfg!(all(target_feature = "avx2", target_feature = "fma")) || {
        let simd = crate::execution::host_environment().simd();
        simd.avx2 && simd.fma
    }
}

#[inline]
#[cfg(target_arch = "aarch64")]
fn has_neon() -> bool {
    cfg!(target_feature = "neon") || crate::execution::host_environment().simd().neon
}

/// Computes the dot product `sum_i a[i] * b[i]` over `0..a.len().min(b.len())`.
#[must_use]
pub fn dot_f64(
    a: &[f64],
    b: &[f64],
) -> f64 {
    let n = a.len().min(b.len());
    if n == 0 {
        return 0.0;
    }
    #[cfg(target_arch = "x86_64")]
    if has_avx2_fma() && n >= 8 {
        // SAFETY: `has_avx2_fma()` verified AVX2 and FMA support at runtime,
        // and `n <= a.len()` and `n <= b.len()`.
        return unsafe { avx2::dot_f64(a.as_ptr(), b.as_ptr(), n) };
    }
    #[cfg(target_arch = "aarch64")]
    if has_neon() && n >= 4 {
        // SAFETY: `has_neon()` verified NEON support, and `n <= a.len().min(b.len())`.
        return unsafe { neon::dot_f64(a.as_ptr(), b.as_ptr(), n) };
    }
    scalar::dot_f64(&a[..n], &b[..n])
}

/// Computes the sum `sum_i x[i]`.
#[must_use]
pub fn sum_f64(x: &[f64]) -> f64 {
    if x.is_empty() {
        return 0.0;
    }
    #[cfg(target_arch = "x86_64")]
    if has_avx2_fma() && x.len() >= 8 {
        // SAFETY: `has_avx2_fma()` verified AVX2 support at runtime.
        return unsafe { avx2::sum_f64(x.as_ptr(), x.len()) };
    }
    #[cfg(target_arch = "aarch64")]
    if has_neon() && x.len() >= 4 {
        // SAFETY: `has_neon()` verified NEON support.
        return unsafe { neon::sum_f64(x.as_ptr(), x.len()) };
    }
    scalar::sum_f64(x)
}

/// Computes the squared Euclidean norm `sum_i x[i] * x[i]`.
#[must_use]
pub fn norm2_sq_f64(x: &[f64]) -> f64 {
    if x.is_empty() {
        return 0.0;
    }
    #[cfg(target_arch = "x86_64")]
    if has_avx2_fma() && x.len() >= 8 {
        // SAFETY: `has_avx2_fma()` verified AVX2 and FMA support at runtime.
        return unsafe { avx2::norm2_sq_f64(x.as_ptr(), x.len()) };
    }
    #[cfg(target_arch = "aarch64")]
    if has_neon() && x.len() >= 4 {
        // SAFETY: `has_neon()` verified NEON support.
        return unsafe { neon::norm2_sq_f64(x.as_ptr(), x.len()) };
    }
    scalar::norm2_sq_f64(x)
}

/// Computes the L1 norm `sum_i |x[i]|`.
#[must_use]
pub fn l1_norm_f64(x: &[f64]) -> f64 {
    if x.is_empty() {
        return 0.0;
    }
    #[cfg(target_arch = "x86_64")]
    if has_avx2_fma() && x.len() >= 8 {
        // SAFETY: `has_avx2_fma()` verified AVX2 support at runtime.
        return unsafe { avx2::l1_norm_f64(x.as_ptr(), x.len()) };
    }
    #[cfg(target_arch = "aarch64")]
    if has_neon() && x.len() >= 4 {
        // SAFETY: `has_neon()` verified NEON support.
        return unsafe { neon::l1_norm_f64(x.as_ptr(), x.len()) };
    }
    scalar::l1_norm_f64(x)
}

/// Computes the L-infinity norm `max_i |x[i]|` (propagates NaN).
#[must_use]
pub fn linf_norm_f64(x: &[f64]) -> f64 {
    if x.is_empty() {
        return 0.0;
    }
    #[cfg(target_arch = "x86_64")]
    if has_avx2_fma() && x.len() >= 8 {
        // SAFETY: `has_avx2_fma()` verified AVX2 support at runtime.
        return unsafe { avx2::linf_norm_f64(x.as_ptr(), x.len()) };
    }
    #[cfg(target_arch = "aarch64")]
    if has_neon() && x.len() >= 4 {
        // SAFETY: `has_neon()` verified NEON support.
        return unsafe { neon::linf_norm_f64(x.as_ptr(), x.len()) };
    }
    scalar::linf_norm_f64(x)
}

/// Computes `y[i] += alpha * x[i]` for `i in 0..x.len().min(y.len())`.
pub fn axpy_f64(
    alpha: f64,
    x: &[f64],
    y: &mut [f64],
) {
    let n = x.len().min(y.len());
    if n == 0 || alpha == 0.0 {
        return;
    }
    #[cfg(target_arch = "x86_64")]
    if has_avx2_fma() && n >= 8 {
        // SAFETY: `has_avx2_fma()` verified AVX2 and FMA support at runtime,
        // and `n <= x.len()` and `n <= y.len()`.
        unsafe { avx2::axpy_f64(alpha, x.as_ptr(), y.as_mut_ptr(), n) };
        return;
    }
    #[cfg(target_arch = "aarch64")]
    if has_neon() && n >= 4 {
        // SAFETY: `has_neon()` verified NEON support, and `n <= x.len().min(y.len())`.
        unsafe { neon::axpy_f64(alpha, x.as_ptr(), y.as_mut_ptr(), n) };
        return;
    }
    scalar::axpy_f64(alpha, &x[..n], &mut y[..n]);
}

/// Computes `y[i] = alpha * x[i] + beta * y[i]` for `i in 0..x.len().min(y.len())`.
pub fn axpby_f64(
    alpha: f64,
    x: &[f64],
    beta: f64,
    y: &mut [f64],
) {
    let n = x.len().min(y.len());
    if n == 0 {
        return;
    }
    #[cfg(target_arch = "x86_64")]
    if has_avx2_fma() && n >= 8 {
        // SAFETY: `has_avx2_fma()` verified AVX2 and FMA support at runtime,
        // and `n <= x.len()` and `n <= y.len()`.
        unsafe { avx2::axpby_f64(alpha, x.as_ptr(), beta, y.as_mut_ptr(), n) };
        return;
    }
    #[cfg(target_arch = "aarch64")]
    if has_neon() && n >= 4 {
        // SAFETY: `has_neon()` verified NEON support, and `n <= x.len().min(y.len())`.
        unsafe { neon::axpby_f64(alpha, x.as_ptr(), beta, y.as_mut_ptr(), n) };
        return;
    }
    scalar::axpby_f64(alpha, &x[..n], beta, &mut y[..n]);
}

/// Scales `x[i] *= alpha` in place.
pub fn scale_f64(
    alpha: f64,
    x: &mut [f64],
) {
    if x.is_empty() {
        return;
    }
    #[cfg(target_arch = "x86_64")]
    if has_avx2_fma() && x.len() >= 8 {
        // SAFETY: `has_avx2_fma()` verified AVX2 support at runtime.
        unsafe { avx2::scale_f64(alpha, x.as_mut_ptr(), x.len()) };
        return;
    }
    #[cfg(target_arch = "aarch64")]
    if has_neon() && x.len() >= 4 {
        // SAFETY: `has_neon()` verified NEON support.
        unsafe { neon::scale_f64(alpha, x.as_mut_ptr(), x.len()) };
        return;
    }
    scalar::scale_f64(alpha, x);
}

/// Computes `out[i] = a[i] + b[i]` over `0..min(a.len(), b.len(), out.len())`.
pub fn add_into_f64(
    a: &[f64],
    b: &[f64],
    out: &mut [f64],
) {
    let n = a.len().min(b.len()).min(out.len());
    if n == 0 {
        return;
    }
    #[cfg(target_arch = "x86_64")]
    if has_avx2_fma() && n >= 8 {
        // SAFETY: `has_avx2_fma()` verified AVX2 support at runtime,
        // and `n` does not exceed any slice length.
        unsafe { avx2::add_into_f64(a.as_ptr(), b.as_ptr(), out.as_mut_ptr(), n) };
        return;
    }
    #[cfg(target_arch = "aarch64")]
    if has_neon() && n >= 4 {
        // SAFETY: `has_neon()` verified NEON support, and `n` does not exceed any slice length.
        unsafe { neon::add_into_f64(a.as_ptr(), b.as_ptr(), out.as_mut_ptr(), n) };
        return;
    }
    scalar::add_into_f64(&a[..n], &b[..n], &mut out[..n]);
}

/// Computes `out[i] = a[i] - b[i]` over `0..min(a.len(), b.len(), out.len())`.
pub fn sub_into_f64(
    a: &[f64],
    b: &[f64],
    out: &mut [f64],
) {
    let n = a.len().min(b.len()).min(out.len());
    if n == 0 {
        return;
    }
    #[cfg(target_arch = "x86_64")]
    if has_avx2_fma() && n >= 8 {
        // SAFETY: `has_avx2_fma()` verified AVX2 support at runtime,
        // and `n` does not exceed any slice length.
        unsafe { avx2::sub_into_f64(a.as_ptr(), b.as_ptr(), out.as_mut_ptr(), n) };
        return;
    }
    #[cfg(target_arch = "aarch64")]
    if has_neon() && n >= 4 {
        // SAFETY: `has_neon()` verified NEON support, and `n` does not exceed any slice length.
        unsafe { neon::sub_into_f64(a.as_ptr(), b.as_ptr(), out.as_mut_ptr(), n) };
        return;
    }
    scalar::sub_into_f64(&a[..n], &b[..n], &mut out[..n]);
}

/// Computes `out[i] = a[i] * b[i]` over `0..min(a.len(), b.len(), out.len())`.
pub fn hadamard_mul_into_f64(
    a: &[f64],
    b: &[f64],
    out: &mut [f64],
) {
    let n = a.len().min(b.len()).min(out.len());
    if n == 0 {
        return;
    }
    #[cfg(target_arch = "x86_64")]
    if has_avx2_fma() && n >= 8 {
        // SAFETY: `has_avx2_fma()` verified AVX2 support at runtime,
        // and `n` does not exceed any slice length.
        unsafe { avx2::hadamard_mul_into_f64(a.as_ptr(), b.as_ptr(), out.as_mut_ptr(), n) };
        return;
    }
    #[cfg(target_arch = "aarch64")]
    if has_neon() && n >= 4 {
        // SAFETY: `has_neon()` verified NEON support, and `n` does not exceed any slice length.
        unsafe { neon::hadamard_mul_into_f64(a.as_ptr(), b.as_ptr(), out.as_mut_ptr(), n) };
        return;
    }
    scalar::hadamard_mul_into_f64(&a[..n], &b[..n], &mut out[..n]);
}

/// Computes `out[i] = a[i] * b[i]` for complex slices.
pub fn complex_mul_into(
    a: &[Complex64],
    b: &[Complex64],
    out: &mut [Complex64],
) {
    let n = a.len().min(b.len()).min(out.len());
    if n == 0 {
        return;
    }
    #[cfg(target_arch = "x86_64")]
    if has_avx2_fma() && n >= 4 {
        // SAFETY: `Complex64` is `#[repr(C)] { re: f64, im: f64 }`, so a slice of
        // `n` `Complex64` values has identical layout to `2 * n` contiguous `f64`s.
        unsafe {
            avx2::complex_mul_into(
                a.as_ptr().cast::<f64>(),
                b.as_ptr().cast::<f64>(),
                out.as_mut_ptr().cast::<f64>(),
                n,
            );
        }
        return;
    }
    #[cfg(target_arch = "aarch64")]
    if has_neon() && n >= 2 {
        // SAFETY: `Complex64` has layout `[f64; 2]` and `n` does not exceed any slice length.
        unsafe {
            neon::complex_mul_into(
                a.as_ptr().cast::<f64>(),
                b.as_ptr().cast::<f64>(),
                out.as_mut_ptr().cast::<f64>(),
                n,
            );
        }
        return;
    }
    scalar::complex_mul_into(&a[..n], &b[..n], &mut out[..n]);
}

/// Computes `out[i] = a[i] * b[i].conj()` for complex slices.
pub fn complex_mul_conj_into(
    a: &[Complex64],
    b: &[Complex64],
    out: &mut [Complex64],
) {
    let n = a.len().min(b.len()).min(out.len());
    if n == 0 {
        return;
    }
    #[cfg(target_arch = "x86_64")]
    if has_avx2_fma() && n >= 4 {
        // SAFETY: `Complex64` is `#[repr(C)] { re: f64, im: f64 }`, matching
        // `2 * n` contiguous `f64` values.
        unsafe {
            avx2::complex_mul_conj_into(
                a.as_ptr().cast::<f64>(),
                b.as_ptr().cast::<f64>(),
                out.as_mut_ptr().cast::<f64>(),
                n,
            );
        }
        return;
    }
    #[cfg(target_arch = "aarch64")]
    if has_neon() && n >= 2 {
        // SAFETY: `Complex64` has layout `[f64; 2]` and `n` does not exceed any slice length.
        unsafe {
            neon::complex_mul_conj_into(
                a.as_ptr().cast::<f64>(),
                b.as_ptr().cast::<f64>(),
                out.as_mut_ptr().cast::<f64>(),
                n,
            );
        }
        return;
    }
    scalar::complex_mul_conj_into(&a[..n], &b[..n], &mut out[..n]);
}

/// Multiplies every element of `data` by real scalar `scale` in place.
pub fn complex_scale_in_place(
    scale: f64,
    data: &mut [Complex64],
) {
    if data.is_empty() {
        return;
    }
    let floats_len = data.len() * 2;
    // SAFETY: `Complex64` is `#[repr(C)] { re: f64, im: f64 }` with alignment
    // equal to `f64`, so reinterpreting `&mut [Complex64]` as `&mut [f64]` of
    // length `2 * len` is valid and preserves lifetime/aliasing rules.
    let floats =
        unsafe { std::slice::from_raw_parts_mut(data.as_mut_ptr().cast::<f64>(), floats_len) };
    scale_f64(scale, floats);
}

/// Multiplies `a[i] *= b[i]` in place over `0..a.len().min(b.len())`.
pub fn complex_mul_in_place(
    a: &mut [Complex64],
    b: &[Complex64],
) {
    let n = a.len().min(b.len());
    if n == 0 {
        return;
    }
    #[cfg(target_arch = "x86_64")]
    if has_avx2_fma() && n >= 4 {
        // SAFETY: `Complex64` has layout `[f64; 2]`; `a` and `b` are disjoint
        // borrows of length at least `n`.
        unsafe {
            avx2::complex_mul_in_place(a.as_mut_ptr().cast::<f64>(), b.as_ptr().cast::<f64>(), n);
        }
        return;
    }
    #[cfg(target_arch = "aarch64")]
    if has_neon() && n >= 2 {
        // SAFETY: `Complex64` has layout `[f64; 2]`; `a` and `b` are disjoint borrows.
        unsafe {
            neon::complex_mul_in_place(a.as_mut_ptr().cast::<f64>(), b.as_ptr().cast::<f64>(), n);
        }
        return;
    }
    scalar::complex_mul_in_place(&mut a[..n], &b[..n]);
}

/// Computes the Hermitian inner product `sum_i a[i].conj() * b[i]` over
/// `0..a.len().min(b.len())`.
#[must_use]
pub fn complex_dot_conj(
    a: &[Complex64],
    b: &[Complex64],
) -> Complex64 {
    let n = a.len().min(b.len());
    if n == 0 {
        return Complex64::new(0.0, 0.0);
    }
    #[cfg(target_arch = "x86_64")]
    if has_avx2_fma() && n >= 4 {
        // SAFETY: `Complex64` has layout `[f64; 2]`, and `n <= a.len().min(b.len())`.
        return unsafe {
            avx2::complex_dot_conj(a.as_ptr().cast::<f64>(), b.as_ptr().cast::<f64>(), n)
        };
    }
    #[cfg(target_arch = "aarch64")]
    if has_neon() && n >= 2 {
        // SAFETY: `Complex64` has layout `[f64; 2]`, and `n <= a.len().min(b.len())`.
        return unsafe {
            neon::complex_dot_conj(a.as_ptr().cast::<f64>(), b.as_ptr().cast::<f64>(), n)
        };
    }
    scalar::complex_dot_conj(&a[..n], &b[..n])
}

/// Computes the squared Euclidean norm `sum_i |x[i]|^2` for a complex slice.
#[must_use]
pub fn complex_norm2_sq(x: &[Complex64]) -> f64 {
    if x.is_empty() {
        return 0.0;
    }
    let floats_len = x.len() * 2;
    // SAFETY: `Complex64` is `#[repr(C)] { re: f64, im: f64 }` with alignment
    // equal to `f64`, so `sum_i (re_i^2 + im_i^2)` equals `norm2_sq_f64` over
    // the reinterpreted `2 * len` slice.
    let floats = unsafe { std::slice::from_raw_parts(x.as_ptr().cast::<f64>(), floats_len) };
    norm2_sq_f64(floats)
}

/// Computes `(min, max)` over `x` in a single pass, returning `None` if empty.
/// Propagates `NaN` if any element is `NaN`.
#[must_use]
pub fn min_max_f64(x: &[f64]) -> Option<(f64, f64)> {
    if x.is_empty() {
        return None;
    }
    #[cfg(target_arch = "x86_64")]
    if has_avx2_fma() && x.len() >= 8 {
        // SAFETY: `has_avx2_fma()` verified AVX2 support at runtime and `x` is non-empty.
        return Some(unsafe { avx2::min_max_f64(x.as_ptr(), x.len()) });
    }
    #[cfg(target_arch = "aarch64")]
    if has_neon() && x.len() >= 4 {
        // SAFETY: `has_neon()` verified NEON support and `x.len() >= 4`.
        return Some(unsafe { neon::min_max_f64(x.as_ptr(), x.len()) });
    }
    Some(scalar::min_max_f64(x))
}

/// Computes `(sum (a - mean_a)*(b - mean_b), sum (a - mean_a)^2, sum (b - mean_b)^2)`
/// in a single fused pass over `0..a.len().min(b.len())`.
#[must_use]
pub fn centered_stats_f64(
    a: &[f64],
    b: &[f64],
    mean_a: f64,
    mean_b: f64,
) -> (f64, f64, f64) {
    let n = a.len().min(b.len());
    if n == 0 {
        return (0.0, 0.0, 0.0);
    }
    #[cfg(target_arch = "x86_64")]
    if has_avx2_fma() && n >= 8 {
        // SAFETY: `has_avx2_fma()` verified AVX2+FMA support and `n <= a.len().min(b.len())`.
        return unsafe { avx2::centered_stats_f64(a.as_ptr(), b.as_ptr(), n, mean_a, mean_b) };
    }
    #[cfg(target_arch = "aarch64")]
    if has_neon() && n >= 4 {
        // SAFETY: `has_neon()` verified NEON support and `n <= a.len().min(b.len())`.
        return unsafe { neon::centered_stats_f64(a.as_ptr(), b.as_ptr(), n, mean_a, mean_b) };
    }
    scalar::centered_stats_f64(&a[..n], &b[..n], mean_a, mean_b)
}

pub(crate) mod scalar {
    use num_complex::Complex64;

    pub fn dot_f64(
        a: &[f64],
        b: &[f64],
    ) -> f64 {
        let mut s0 = 0.0_f64;
        let mut s1 = 0.0_f64;
        let mut s2 = 0.0_f64;
        let mut s3 = 0.0_f64;
        let (a_chunks, a_tail) = a.as_chunks::<4>();
        let (b_chunks, b_tail) = b.as_chunks::<4>();
        for (xa, xb) in a_chunks.iter().zip(b_chunks) {
            s0 = xa[0].mul_add(xb[0], s0);
            s1 = xa[1].mul_add(xb[1], s1);
            s2 = xa[2].mul_add(xb[2], s2);
            s3 = xa[3].mul_add(xb[3], s3);
        }
        let mut total = (s0 + s1) + (s2 + s3);
        for (&xa, &xb) in a_tail.iter().zip(b_tail) {
            total = xa.mul_add(xb, total);
        }
        total
    }

    pub fn sum_f64(x: &[f64]) -> f64 {
        let mut s0 = 0.0_f64;
        let mut s1 = 0.0_f64;
        let mut s2 = 0.0_f64;
        let mut s3 = 0.0_f64;
        let (chunks, chunk_tail) = x.as_chunks::<4>();
        for c in chunks {
            s0 += c[0];
            s1 += c[1];
            s2 += c[2];
            s3 += c[3];
        }
        let mut total = (s0 + s1) + (s2 + s3);
        for &v in chunk_tail {
            total += v;
        }
        total
    }

    pub fn norm2_sq_f64(x: &[f64]) -> f64 {
        let mut s0 = 0.0_f64;
        let mut s1 = 0.0_f64;
        let mut s2 = 0.0_f64;
        let mut s3 = 0.0_f64;
        let (chunks, chunk_tail) = x.as_chunks::<4>();
        for c in chunks {
            s0 = c[0].mul_add(c[0], s0);
            s1 = c[1].mul_add(c[1], s1);
            s2 = c[2].mul_add(c[2], s2);
            s3 = c[3].mul_add(c[3], s3);
        }
        let mut total = (s0 + s1) + (s2 + s3);
        for &v in chunk_tail {
            total = v.mul_add(v, total);
        }
        total
    }

    pub fn l1_norm_f64(x: &[f64]) -> f64 {
        let mut s0 = 0.0_f64;
        let mut s1 = 0.0_f64;
        let mut s2 = 0.0_f64;
        let mut s3 = 0.0_f64;
        let (chunks, chunk_tail) = x.as_chunks::<4>();
        for c in chunks {
            s0 += c[0].abs();
            s1 += c[1].abs();
            s2 += c[2].abs();
            s3 += c[3].abs();
        }
        let mut total = (s0 + s1) + (s2 + s3);
        for &v in chunk_tail {
            total += v.abs();
        }
        total
    }

    pub fn linf_norm_f64(x: &[f64]) -> f64 {
        let mut best = 0.0_f64;
        for &v in x {
            let a = v.abs();
            if a.is_nan() {
                return f64::NAN;
            }
            if a > best {
                best = a;
            }
        }
        best
    }

    pub fn axpy_f64(
        alpha: f64,
        x: &[f64],
        y: &mut [f64],
    ) {
        for (&xi, yi) in x.iter().zip(y.iter_mut()) {
            *yi = alpha.mul_add(xi, *yi);
        }
    }

    pub fn axpby_f64(
        alpha: f64,
        x: &[f64],
        beta: f64,
        y: &mut [f64],
    ) {
        for (&xi, yi) in x.iter().zip(y.iter_mut()) {
            *yi = alpha.mul_add(xi, beta * *yi);
        }
    }

    pub fn scale_f64(
        alpha: f64,
        x: &mut [f64],
    ) {
        for xi in x.iter_mut() {
            *xi *= alpha;
        }
    }

    pub fn add_into_f64(
        a: &[f64],
        b: &[f64],
        out: &mut [f64],
    ) {
        for ((&ai, &bi), oi) in a.iter().zip(b.iter()).zip(out.iter_mut()) {
            *oi = ai + bi;
        }
    }

    pub fn sub_into_f64(
        a: &[f64],
        b: &[f64],
        out: &mut [f64],
    ) {
        for ((&ai, &bi), oi) in a.iter().zip(b.iter()).zip(out.iter_mut()) {
            *oi = ai - bi;
        }
    }

    pub fn hadamard_mul_into_f64(
        a: &[f64],
        b: &[f64],
        out: &mut [f64],
    ) {
        for ((&ai, &bi), oi) in a.iter().zip(b.iter()).zip(out.iter_mut()) {
            *oi = ai * bi;
        }
    }

    pub fn complex_mul_into(
        a: &[Complex64],
        b: &[Complex64],
        out: &mut [Complex64],
    ) {
        for ((&ai, &bi), oi) in a.iter().zip(b.iter()).zip(out.iter_mut()) {
            *oi = ai * bi;
        }
    }

    pub fn complex_mul_conj_into(
        a: &[Complex64],
        b: &[Complex64],
        out: &mut [Complex64],
    ) {
        for ((&ai, &bi), oi) in a.iter().zip(b.iter()).zip(out.iter_mut()) {
            *oi = ai * bi.conj();
        }
    }

    pub fn complex_mul_in_place(
        a: &mut [Complex64],
        b: &[Complex64],
    ) {
        for (ai, &bi) in a.iter_mut().zip(b.iter()) {
            *ai *= bi;
        }
    }

    pub fn complex_dot_conj(
        a: &[Complex64],
        b: &[Complex64],
    ) -> Complex64 {
        let mut re = 0.0_f64;
        let mut im = 0.0_f64;
        for (&ai, &bi) in a.iter().zip(b.iter()) {
            re = ai.re.mul_add(bi.re, ai.im.mul_add(bi.im, re));
            im = ai.re.mul_add(bi.im, (-ai.im).mul_add(bi.re, im));
        }
        Complex64::new(re, im)
    }

    pub fn min_max_f64(x: &[f64]) -> (f64, f64) {
        let mut lo = x[0];
        let mut hi = x[0];
        if lo.is_nan() {
            return (f64::NAN, f64::NAN);
        }
        for &v in &x[1..] {
            if v.is_nan() {
                return (f64::NAN, f64::NAN);
            }
            if v < lo {
                lo = v;
            }
            if v > hi {
                hi = v;
            }
        }
        (lo, hi)
    }

    pub fn centered_stats_f64(
        a: &[f64],
        b: &[f64],
        mean_a: f64,
        mean_b: f64,
    ) -> (f64, f64, f64) {
        let mut cov = 0.0_f64;
        let mut va = 0.0_f64;
        let mut vb = 0.0_f64;
        for (&xa, &xb) in a.iter().zip(b.iter()) {
            let da = xa - mean_a;
            let db = xb - mean_b;
            cov = da.mul_add(db, cov);
            va = da.mul_add(da, va);
            vb = db.mul_add(db, vb);
        }
        (cov, va, vb)
    }
}

#[cfg(target_arch = "x86_64")]
mod avx2 {
    use num_complex::Complex64;
    use std::arch::x86_64::__m256d;
    use std::arch::x86_64::_CMP_UNORD_Q;
    use std::arch::x86_64::_mm256_add_pd;
    use std::arch::x86_64::_mm256_addsub_pd;
    use std::arch::x86_64::_mm256_andnot_pd;
    use std::arch::x86_64::_mm256_cmp_pd;
    use std::arch::x86_64::_mm256_fmadd_pd;
    use std::arch::x86_64::_mm256_loadu_pd;
    use std::arch::x86_64::_mm256_max_pd;
    use std::arch::x86_64::_mm256_min_pd;
    use std::arch::x86_64::_mm256_movedup_pd;
    use std::arch::x86_64::_mm256_movemask_pd;
    use std::arch::x86_64::_mm256_mul_pd;
    use std::arch::x86_64::_mm256_or_pd;
    use std::arch::x86_64::_mm256_permute_pd;
    use std::arch::x86_64::_mm256_set_pd;
    use std::arch::x86_64::_mm256_set1_pd;
    use std::arch::x86_64::_mm256_setzero_pd;
    use std::arch::x86_64::_mm256_storeu_pd;
    use std::arch::x86_64::_mm256_sub_pd;
    use std::arch::x86_64::_mm256_unpackhi_pd;
    use std::arch::x86_64::_mm256_xor_pd;

    #[inline]
    #[target_feature(enable = "avx2")]
    unsafe fn hsum256(v: __m256d) -> f64 {
        let mut lanes = [0.0_f64; 4];
        // SAFETY: `lanes` has room for 4 `f64` values.
        unsafe { _mm256_storeu_pd(lanes.as_mut_ptr(), v) };
        (lanes[0] + lanes[1]) + (lanes[2] + lanes[3])
    }

    #[target_feature(enable = "avx2,fma")]
    pub(super) unsafe fn dot_f64(
        a: *const f64,
        b: *const f64,
        n: usize,
    ) -> f64 {
        let mut acc0 = _mm256_setzero_pd();
        let mut acc1 = _mm256_setzero_pd();
        let mut i = 0_usize;
        while i + 8 <= n {
            // SAFETY: `i + 8 <= n`, so loads at `i` and `i + 4` are in bounds.
            unsafe {
                let va0 = _mm256_loadu_pd(a.add(i));
                let vb0 = _mm256_loadu_pd(b.add(i));
                let va1 = _mm256_loadu_pd(a.add(i + 4));
                let vb1 = _mm256_loadu_pd(b.add(i + 4));
                acc0 = _mm256_fmadd_pd(va0, vb0, acc0);
                acc1 = _mm256_fmadd_pd(va1, vb1, acc1);
            }
            i += 8;
        }
        if i + 4 <= n {
            // SAFETY: `i + 4 <= n`.
            unsafe {
                let va = _mm256_loadu_pd(a.add(i));
                let vb = _mm256_loadu_pd(b.add(i));
                acc0 = _mm256_fmadd_pd(va, vb, acc0);
            }
            i += 4;
        }
        // SAFETY: AVX2 is enabled by this function's target_feature.
        let mut total = unsafe { hsum256(_mm256_add_pd(acc0, acc1)) };
        while i < n {
            // SAFETY: `i < n`.
            unsafe {
                total = (*a.add(i)).mul_add(*b.add(i), total);
            }
            i += 1;
        }
        total
    }

    #[target_feature(enable = "avx2,fma")]
    pub(super) unsafe fn sum_f64(
        x: *const f64,
        n: usize,
    ) -> f64 {
        let mut acc0 = _mm256_setzero_pd();
        let mut acc1 = _mm256_setzero_pd();
        let mut i = 0_usize;
        while i + 8 <= n {
            // SAFETY: `i + 8 <= n`.
            unsafe {
                acc0 = _mm256_add_pd(acc0, _mm256_loadu_pd(x.add(i)));
                acc1 = _mm256_add_pd(acc1, _mm256_loadu_pd(x.add(i + 4)));
            }
            i += 8;
        }
        if i + 4 <= n {
            // SAFETY: `i + 4 <= n`.
            unsafe {
                acc0 = _mm256_add_pd(acc0, _mm256_loadu_pd(x.add(i)));
            }
            i += 4;
        }
        // SAFETY: AVX2 is enabled.
        let mut total = unsafe { hsum256(_mm256_add_pd(acc0, acc1)) };
        while i < n {
            // SAFETY: `i < n`.
            unsafe {
                total += *x.add(i);
            }
            i += 1;
        }
        total
    }

    #[target_feature(enable = "avx2,fma")]
    pub(super) unsafe fn norm2_sq_f64(
        x: *const f64,
        n: usize,
    ) -> f64 {
        let mut acc0 = _mm256_setzero_pd();
        let mut acc1 = _mm256_setzero_pd();
        let mut i = 0_usize;
        while i + 8 <= n {
            // SAFETY: `i + 8 <= n`.
            unsafe {
                let v0 = _mm256_loadu_pd(x.add(i));
                let v1 = _mm256_loadu_pd(x.add(i + 4));
                acc0 = _mm256_fmadd_pd(v0, v0, acc0);
                acc1 = _mm256_fmadd_pd(v1, v1, acc1);
            }
            i += 8;
        }
        if i + 4 <= n {
            // SAFETY: `i + 4 <= n`.
            unsafe {
                let v = _mm256_loadu_pd(x.add(i));
                acc0 = _mm256_fmadd_pd(v, v, acc0);
            }
            i += 4;
        }
        // SAFETY: AVX2 is enabled.
        let mut total = unsafe { hsum256(_mm256_add_pd(acc0, acc1)) };
        while i < n {
            // SAFETY: `i < n`.
            unsafe {
                let v = *x.add(i);
                total = v.mul_add(v, total);
            }
            i += 1;
        }
        total
    }

    #[target_feature(enable = "avx2,fma")]
    pub(super) unsafe fn l1_norm_f64(
        x: *const f64,
        n: usize,
    ) -> f64 {
        let sign_mask = _mm256_set1_pd(-0.0_f64);
        let mut acc0 = _mm256_setzero_pd();
        let mut acc1 = _mm256_setzero_pd();
        let mut i = 0_usize;
        while i + 8 <= n {
            // SAFETY: `i + 8 <= n`.
            unsafe {
                let v0 = _mm256_andnot_pd(sign_mask, _mm256_loadu_pd(x.add(i)));
                let v1 = _mm256_andnot_pd(sign_mask, _mm256_loadu_pd(x.add(i + 4)));
                acc0 = _mm256_add_pd(acc0, v0);
                acc1 = _mm256_add_pd(acc1, v1);
            }
            i += 8;
        }
        if i + 4 <= n {
            // SAFETY: `i + 4 <= n`.
            unsafe {
                let v = _mm256_andnot_pd(sign_mask, _mm256_loadu_pd(x.add(i)));
                acc0 = _mm256_add_pd(acc0, v);
            }
            i += 4;
        }
        // SAFETY: AVX2 is enabled.
        let mut total = unsafe { hsum256(_mm256_add_pd(acc0, acc1)) };
        while i < n {
            // SAFETY: `i < n`.
            unsafe {
                total += (*x.add(i)).abs();
            }
            i += 1;
        }
        total
    }

    #[target_feature(enable = "avx2,fma")]
    pub(super) unsafe fn linf_norm_f64(
        x: *const f64,
        n: usize,
    ) -> f64 {
        let sign_mask = _mm256_set1_pd(-0.0_f64);
        let mut vmax0 = _mm256_setzero_pd();
        let mut vmax1 = _mm256_setzero_pd();
        let mut nan_mask = _mm256_setzero_pd();
        let mut i = 0_usize;
        while i + 8 <= n {
            // SAFETY: `i + 8 <= n`.
            unsafe {
                let raw0 = _mm256_loadu_pd(x.add(i));
                let raw1 = _mm256_loadu_pd(x.add(i + 4));
                nan_mask = _mm256_or_pd(nan_mask, _mm256_cmp_pd::<_CMP_UNORD_Q>(raw0, raw0));
                nan_mask = _mm256_or_pd(nan_mask, _mm256_cmp_pd::<_CMP_UNORD_Q>(raw1, raw1));
                let v0 = _mm256_andnot_pd(sign_mask, raw0);
                let v1 = _mm256_andnot_pd(sign_mask, raw1);
                vmax0 = _mm256_max_pd(vmax0, v0);
                vmax1 = _mm256_max_pd(vmax1, v1);
            }
            i += 8;
        }
        if i + 4 <= n {
            // SAFETY: `i + 4 <= n`.
            unsafe {
                let raw = _mm256_loadu_pd(x.add(i));
                nan_mask = _mm256_or_pd(nan_mask, _mm256_cmp_pd::<_CMP_UNORD_Q>(raw, raw));
                let v = _mm256_andnot_pd(sign_mask, raw);
                vmax0 = _mm256_max_pd(vmax0, v);
            }
            i += 4;
        }
        if _mm256_movemask_pd(nan_mask) != 0 {
            return f64::NAN;
        }
        let vmax = _mm256_max_pd(vmax0, vmax1);
        let mut lanes = [0.0_f64; 4];
        // SAFETY: `lanes` holds 4 `f64` values.
        unsafe { _mm256_storeu_pd(lanes.as_mut_ptr(), vmax) };
        let mut best = lanes[0].max(lanes[1]).max(lanes[2].max(lanes[3]));
        while i < n {
            // SAFETY: `i < n`.
            unsafe {
                let a = (*x.add(i)).abs();
                if a.is_nan() {
                    return f64::NAN;
                }
                if a > best {
                    best = a;
                }
            }
            i += 1;
        }
        best
    }

    #[target_feature(enable = "avx2,fma")]
    pub(super) unsafe fn axpy_f64(
        alpha: f64,
        x: *const f64,
        y: *mut f64,
        n: usize,
    ) {
        let va = _mm256_set1_pd(alpha);
        let mut i = 0_usize;
        while i + 8 <= n {
            // SAFETY: `i + 8 <= n`.
            unsafe {
                let vx0 = _mm256_loadu_pd(x.add(i));
                let vy0 = _mm256_loadu_pd(y.add(i));
                let vx1 = _mm256_loadu_pd(x.add(i + 4));
                let vy1 = _mm256_loadu_pd(y.add(i + 4));
                _mm256_storeu_pd(y.add(i), _mm256_fmadd_pd(va, vx0, vy0));
                _mm256_storeu_pd(y.add(i + 4), _mm256_fmadd_pd(va, vx1, vy1));
            }
            i += 8;
        }
        if i + 4 <= n {
            // SAFETY: `i + 4 <= n`.
            unsafe {
                let vx = _mm256_loadu_pd(x.add(i));
                let vy = _mm256_loadu_pd(y.add(i));
                _mm256_storeu_pd(y.add(i), _mm256_fmadd_pd(va, vx, vy));
            }
            i += 4;
        }
        while i < n {
            // SAFETY: `i < n`.
            unsafe {
                *y.add(i) = alpha.mul_add(*x.add(i), *y.add(i));
            }
            i += 1;
        }
    }

    #[target_feature(enable = "avx2,fma")]
    pub(super) unsafe fn axpby_f64(
        alpha: f64,
        x: *const f64,
        beta: f64,
        y: *mut f64,
        n: usize,
    ) {
        let va = _mm256_set1_pd(alpha);
        let vb = _mm256_set1_pd(beta);
        let mut i = 0_usize;
        while i + 8 <= n {
            // SAFETY: `i + 8 <= n`.
            unsafe {
                let vx0 = _mm256_loadu_pd(x.add(i));
                let vy0 = _mm256_loadu_pd(y.add(i));
                let vx1 = _mm256_loadu_pd(x.add(i + 4));
                let vy1 = _mm256_loadu_pd(y.add(i + 4));
                let by0 = _mm256_mul_pd(vb, vy0);
                let by1 = _mm256_mul_pd(vb, vy1);
                _mm256_storeu_pd(y.add(i), _mm256_fmadd_pd(va, vx0, by0));
                _mm256_storeu_pd(y.add(i + 4), _mm256_fmadd_pd(va, vx1, by1));
            }
            i += 8;
        }
        if i + 4 <= n {
            // SAFETY: `i + 4 <= n`.
            unsafe {
                let vx = _mm256_loadu_pd(x.add(i));
                let vy = _mm256_loadu_pd(y.add(i));
                let by = _mm256_mul_pd(vb, vy);
                _mm256_storeu_pd(y.add(i), _mm256_fmadd_pd(va, vx, by));
            }
            i += 4;
        }
        while i < n {
            // SAFETY: `i < n`.
            unsafe {
                *y.add(i) = alpha.mul_add(*x.add(i), beta * *y.add(i));
            }
            i += 1;
        }
    }

    #[target_feature(enable = "avx2,fma")]
    pub(super) unsafe fn scale_f64(
        alpha: f64,
        x: *mut f64,
        n: usize,
    ) {
        let va = _mm256_set1_pd(alpha);
        let mut i = 0_usize;
        while i + 8 <= n {
            // SAFETY: `i + 8 <= n`.
            unsafe {
                let v0 = _mm256_loadu_pd(x.add(i));
                let v1 = _mm256_loadu_pd(x.add(i + 4));
                _mm256_storeu_pd(x.add(i), _mm256_mul_pd(va, v0));
                _mm256_storeu_pd(x.add(i + 4), _mm256_mul_pd(va, v1));
            }
            i += 8;
        }
        if i + 4 <= n {
            // SAFETY: `i + 4 <= n`.
            unsafe {
                let v = _mm256_loadu_pd(x.add(i));
                _mm256_storeu_pd(x.add(i), _mm256_mul_pd(va, v));
            }
            i += 4;
        }
        while i < n {
            // SAFETY: `i < n`.
            unsafe {
                *x.add(i) *= alpha;
            }
            i += 1;
        }
    }

    #[target_feature(enable = "avx2,fma")]
    pub(super) unsafe fn add_into_f64(
        a: *const f64,
        b: *const f64,
        out: *mut f64,
        n: usize,
    ) {
        let mut i = 0_usize;
        while i + 8 <= n {
            // SAFETY: `i + 8 <= n`.
            unsafe {
                let r0 = _mm256_add_pd(_mm256_loadu_pd(a.add(i)), _mm256_loadu_pd(b.add(i)));
                let r1 =
                    _mm256_add_pd(_mm256_loadu_pd(a.add(i + 4)), _mm256_loadu_pd(b.add(i + 4)));
                _mm256_storeu_pd(out.add(i), r0);
                _mm256_storeu_pd(out.add(i + 4), r1);
            }
            i += 8;
        }
        if i + 4 <= n {
            // SAFETY: `i + 4 <= n`.
            unsafe {
                let r = _mm256_add_pd(_mm256_loadu_pd(a.add(i)), _mm256_loadu_pd(b.add(i)));
                _mm256_storeu_pd(out.add(i), r);
            }
            i += 4;
        }
        while i < n {
            // SAFETY: `i < n`.
            unsafe {
                *out.add(i) = *a.add(i) + *b.add(i);
            }
            i += 1;
        }
    }

    #[target_feature(enable = "avx2,fma")]
    pub(super) unsafe fn sub_into_f64(
        a: *const f64,
        b: *const f64,
        out: *mut f64,
        n: usize,
    ) {
        let mut i = 0_usize;
        while i + 8 <= n {
            // SAFETY: `i + 8 <= n`.
            unsafe {
                let r0 = _mm256_sub_pd(_mm256_loadu_pd(a.add(i)), _mm256_loadu_pd(b.add(i)));
                let r1 =
                    _mm256_sub_pd(_mm256_loadu_pd(a.add(i + 4)), _mm256_loadu_pd(b.add(i + 4)));
                _mm256_storeu_pd(out.add(i), r0);
                _mm256_storeu_pd(out.add(i + 4), r1);
            }
            i += 8;
        }
        if i + 4 <= n {
            // SAFETY: `i + 4 <= n`.
            unsafe {
                let r = _mm256_sub_pd(_mm256_loadu_pd(a.add(i)), _mm256_loadu_pd(b.add(i)));
                _mm256_storeu_pd(out.add(i), r);
            }
            i += 4;
        }
        while i < n {
            // SAFETY: `i < n`.
            unsafe {
                *out.add(i) = *a.add(i) - *b.add(i);
            }
            i += 1;
        }
    }

    #[target_feature(enable = "avx2,fma")]
    pub(super) unsafe fn hadamard_mul_into_f64(
        a: *const f64,
        b: *const f64,
        out: *mut f64,
        n: usize,
    ) {
        let mut i = 0_usize;
        while i + 8 <= n {
            // SAFETY: `i + 8 <= n`.
            unsafe {
                let r0 = _mm256_mul_pd(_mm256_loadu_pd(a.add(i)), _mm256_loadu_pd(b.add(i)));
                let r1 =
                    _mm256_mul_pd(_mm256_loadu_pd(a.add(i + 4)), _mm256_loadu_pd(b.add(i + 4)));
                _mm256_storeu_pd(out.add(i), r0);
                _mm256_storeu_pd(out.add(i + 4), r1);
            }
            i += 8;
        }
        if i + 4 <= n {
            // SAFETY: `i + 4 <= n`.
            unsafe {
                let r = _mm256_mul_pd(_mm256_loadu_pd(a.add(i)), _mm256_loadu_pd(b.add(i)));
                _mm256_storeu_pd(out.add(i), r);
            }
            i += 4;
        }
        while i < n {
            // SAFETY: `i < n`.
            unsafe {
                *out.add(i) = *a.add(i) * *b.add(i);
            }
            i += 1;
        }
    }

    /// Processes 2 complex numbers (`[re0, im0, re1, im1]`) per AVX2 register.
    #[target_feature(enable = "avx2,fma")]
    pub(super) unsafe fn complex_mul_into(
        a: *const f64,
        b: *const f64,
        out: *mut f64,
        n: usize,
    ) {
        let total_f64 = n * 2;
        let mut i = 0_usize;
        while i + 4 <= total_f64 {
            // SAFETY: `i + 4 <= 2 * n`.
            unsafe {
                let va = _mm256_loadu_pd(a.add(i));
                let vb = _mm256_loadu_pd(b.add(i));
                let a_re = _mm256_movedup_pd(va);
                let a_im = _mm256_unpackhi_pd(va, va);
                let b_swap = _mm256_permute_pd::<0b0101>(vb);
                let prod1 = _mm256_mul_pd(a_re, vb);
                let prod2 = _mm256_mul_pd(a_im, b_swap);
                let res = _mm256_addsub_pd(prod1, prod2);
                _mm256_storeu_pd(out.add(i), res);
            }
            i += 4;
        }
        while i < total_f64 {
            // SAFETY: `i + 1 < 2 * n`.
            unsafe {
                let ar = *a.add(i);
                let ai = *a.add(i + 1);
                let br = *b.add(i);
                let bi = *b.add(i + 1);
                *out.add(i) = ar * br - ai * bi;
                *out.add(i + 1) = ar * bi + ai * br;
            }
            i += 2;
        }
    }

    /// Computes `a * conj(b)` for 2 complex numbers per AVX2 register.
    #[target_feature(enable = "avx2,fma")]
    pub(super) unsafe fn complex_mul_conj_into(
        a: *const f64,
        b: *const f64,
        out: *mut f64,
        n: usize,
    ) {
        let conj_mask = _mm256_set_pd(-0.0, 0.0, -0.0, 0.0);
        let total_f64 = n * 2;
        let mut i = 0_usize;
        while i + 4 <= total_f64 {
            // SAFETY: `i + 4 <= 2 * n`.
            unsafe {
                let va = _mm256_loadu_pd(a.add(i));
                let vb = _mm256_loadu_pd(b.add(i));
                let vb_conj = _mm256_xor_pd(vb, conj_mask);
                let a_re = _mm256_movedup_pd(va);
                let a_im = _mm256_unpackhi_pd(va, va);
                let b_conj_swap = _mm256_permute_pd::<0b0101>(vb_conj);
                let p1 = _mm256_mul_pd(a_re, vb_conj);
                let p2 = _mm256_mul_pd(a_im, b_conj_swap);
                _mm256_storeu_pd(out.add(i), _mm256_addsub_pd(p1, p2));
            }
            i += 4;
        }
        while i < total_f64 {
            // SAFETY: `i + 1 < 2 * n`.
            unsafe {
                let ar = *a.add(i);
                let ai = *a.add(i + 1);
                let br = *b.add(i);
                let bi = *b.add(i + 1);
                *out.add(i) = ar * br + ai * bi;
                *out.add(i + 1) = ai * br - ar * bi;
            }
            i += 2;
        }
    }

    /// Multiplies `a[k] *= b[k]` in place, 2 complex numbers per AVX2 register.
    #[target_feature(enable = "avx2,fma")]
    pub(super) unsafe fn complex_mul_in_place(
        a: *mut f64,
        b: *const f64,
        n: usize,
    ) {
        let total_f64 = n * 2;
        let mut i = 0_usize;
        while i + 4 <= total_f64 {
            // SAFETY: `i + 4 <= 2 * n`.
            unsafe {
                let va = _mm256_loadu_pd(a.add(i));
                let vb = _mm256_loadu_pd(b.add(i));
                let a_re = _mm256_movedup_pd(va);
                let a_im = _mm256_unpackhi_pd(va, va);
                let b_swap = _mm256_permute_pd::<0b0101>(vb);
                let prod1 = _mm256_mul_pd(a_re, vb);
                let prod2 = _mm256_mul_pd(a_im, b_swap);
                _mm256_storeu_pd(a.add(i), _mm256_addsub_pd(prod1, prod2));
            }
            i += 4;
        }
        while i < total_f64 {
            // SAFETY: `i + 1 < 2 * n`.
            unsafe {
                let ar = *a.add(i);
                let ai = *a.add(i + 1);
                let br = *b.add(i);
                let bi = *b.add(i + 1);
                *a.add(i) = ar * br - ai * bi;
                *a.add(i + 1) = ar * bi + ai * br;
            }
            i += 2;
        }
    }

    /// Computes Hermitian inner product `sum_k conj(a[k]) * b[k]`.
    #[target_feature(enable = "avx2,fma")]
    pub(super) unsafe fn complex_dot_conj(
        a: *const f64,
        b: *const f64,
        n: usize,
    ) -> Complex64 {
        let total_f64 = n * 2;
        let mut acc_re = _mm256_setzero_pd();
        let mut acc_im = _mm256_setzero_pd();
        // Sign flip on odd lanes (imag positions) so add-horizontally yields
        // `ar * bi - ai * br` after swapping `b`.
        let neg_odd = _mm256_set_pd(-0.0, 0.0, -0.0, 0.0);
        let mut i = 0_usize;
        while i + 4 <= total_f64 {
            // SAFETY: `i + 4 <= 2 * n`.
            unsafe {
                let va = _mm256_loadu_pd(a.add(i));
                let vb = _mm256_loadu_pd(b.add(i));
                // Real part: `ar * br + ai * bi` across all 4 lanes.
                acc_re = _mm256_fmadd_pd(va, vb, acc_re);
                // Imag part: `va_signed = [ar0, -ai0, ar1, -ai1]`,
                // `vb_swap = [bi0, br0, bi1, br1]` -> lane products are
                // `[ar0*bi0, -ai0*br0, ar1*bi1, -ai1*br1]`.
                let va_conj = _mm256_xor_pd(va, neg_odd);
                let vb_swap = _mm256_permute_pd::<0b0101>(vb);
                acc_im = _mm256_fmadd_pd(va_conj, vb_swap, acc_im);
            }
            i += 4;
        }
        // SAFETY: AVX2 is enabled.
        let mut re = unsafe { hsum256(acc_re) };
        let mut im = unsafe { hsum256(acc_im) };
        while i < total_f64 {
            // SAFETY: `i + 1 < 2 * n`.
            unsafe {
                let ar = *a.add(i);
                let ai = *a.add(i + 1);
                let br = *b.add(i);
                let bi = *b.add(i + 1);
                re = ar.mul_add(br, ai.mul_add(bi, re));
                im = ar.mul_add(bi, (-ai).mul_add(br, im));
            }
            i += 2;
        }
        Complex64::new(re, im)
    }

    #[target_feature(enable = "avx2,fma")]
    pub(super) unsafe fn min_max_f64(
        x: *const f64,
        n: usize,
    ) -> (f64, f64) {
        // SAFETY: caller guarantees `n >= 8`.
        let first = unsafe { _mm256_loadu_pd(x) };
        let mut vmin0 = first;
        let mut vmax0 = first;
        let mut vmin1 = first;
        let mut vmax1 = first;
        let mut nan_mask = _mm256_cmp_pd::<_CMP_UNORD_Q>(first, first);
        let mut i = 4_usize;
        while i + 8 <= n {
            // SAFETY: `i + 8 <= n`.
            unsafe {
                let v0 = _mm256_loadu_pd(x.add(i));
                let v1 = _mm256_loadu_pd(x.add(i + 4));
                nan_mask = _mm256_or_pd(nan_mask, _mm256_cmp_pd::<_CMP_UNORD_Q>(v0, v0));
                nan_mask = _mm256_or_pd(nan_mask, _mm256_cmp_pd::<_CMP_UNORD_Q>(v1, v1));
                vmin0 = _mm256_min_pd(vmin0, v0);
                vmax0 = _mm256_max_pd(vmax0, v0);
                vmin1 = _mm256_min_pd(vmin1, v1);
                vmax1 = _mm256_max_pd(vmax1, v1);
            }
            i += 8;
        }
        if i + 4 <= n {
            // SAFETY: `i + 4 <= n`.
            unsafe {
                let v = _mm256_loadu_pd(x.add(i));
                nan_mask = _mm256_or_pd(nan_mask, _mm256_cmp_pd::<_CMP_UNORD_Q>(v, v));
                vmin0 = _mm256_min_pd(vmin0, v);
                vmax0 = _mm256_max_pd(vmax0, v);
            }
            i += 4;
        }
        if _mm256_movemask_pd(nan_mask) != 0 {
            return (f64::NAN, f64::NAN);
        }
        let vmin = _mm256_min_pd(vmin0, vmin1);
        let vmax = _mm256_max_pd(vmax0, vmax1);
        let mut lo_lanes = [0.0_f64; 4];
        let mut hi_lanes = [0.0_f64; 4];
        // SAFETY: buffers hold 4 `f64` values each.
        unsafe {
            _mm256_storeu_pd(lo_lanes.as_mut_ptr(), vmin);
            _mm256_storeu_pd(hi_lanes.as_mut_ptr(), vmax);
        }
        let mut lo = lo_lanes[0]
            .min(lo_lanes[1])
            .min(lo_lanes[2].min(lo_lanes[3]));
        let mut hi = hi_lanes[0]
            .max(hi_lanes[1])
            .max(hi_lanes[2].max(hi_lanes[3]));
        while i < n {
            // SAFETY: `i < n`.
            unsafe {
                let v = *x.add(i);
                if v.is_nan() {
                    return (f64::NAN, f64::NAN);
                }
                if v < lo {
                    lo = v;
                }
                if v > hi {
                    hi = v;
                }
            }
            i += 1;
        }
        (lo, hi)
    }

    #[target_feature(enable = "avx2,fma")]
    pub(super) unsafe fn centered_stats_f64(
        a: *const f64,
        b: *const f64,
        n: usize,
        mean_a: f64,
        mean_b: f64,
    ) -> (f64, f64, f64) {
        let vma = _mm256_set1_pd(mean_a);
        let vmb = _mm256_set1_pd(mean_b);
        let mut cov0 = _mm256_setzero_pd();
        let mut va0 = _mm256_setzero_pd();
        let mut vb0 = _mm256_setzero_pd();
        let mut cov1 = _mm256_setzero_pd();
        let mut va1 = _mm256_setzero_pd();
        let mut vb1 = _mm256_setzero_pd();
        let mut i = 0_usize;
        while i + 8 <= n {
            // SAFETY: `i + 8 <= n`.
            unsafe {
                let da0 = _mm256_sub_pd(_mm256_loadu_pd(a.add(i)), vma);
                let db0 = _mm256_sub_pd(_mm256_loadu_pd(b.add(i)), vmb);
                let da1 = _mm256_sub_pd(_mm256_loadu_pd(a.add(i + 4)), vma);
                let db1 = _mm256_sub_pd(_mm256_loadu_pd(b.add(i + 4)), vmb);
                cov0 = _mm256_fmadd_pd(da0, db0, cov0);
                va0 = _mm256_fmadd_pd(da0, da0, va0);
                vb0 = _mm256_fmadd_pd(db0, db0, vb0);
                cov1 = _mm256_fmadd_pd(da1, db1, cov1);
                va1 = _mm256_fmadd_pd(da1, da1, va1);
                vb1 = _mm256_fmadd_pd(db1, db1, vb1);
            }
            i += 8;
        }
        if i + 4 <= n {
            // SAFETY: `i + 4 <= n`.
            unsafe {
                let da = _mm256_sub_pd(_mm256_loadu_pd(a.add(i)), vma);
                let db = _mm256_sub_pd(_mm256_loadu_pd(b.add(i)), vmb);
                cov0 = _mm256_fmadd_pd(da, db, cov0);
                va0 = _mm256_fmadd_pd(da, da, va0);
                vb0 = _mm256_fmadd_pd(db, db, vb0);
            }
            i += 4;
        }
        // SAFETY: AVX2 is enabled.
        let mut cov = unsafe { hsum256(_mm256_add_pd(cov0, cov1)) };
        let mut va = unsafe { hsum256(_mm256_add_pd(va0, va1)) };
        let mut vb = unsafe { hsum256(_mm256_add_pd(vb0, vb1)) };
        while i < n {
            // SAFETY: `i < n`.
            unsafe {
                let da = *a.add(i) - mean_a;
                let db = *b.add(i) - mean_b;
                cov = da.mul_add(db, cov);
                va = da.mul_add(da, va);
                vb = db.mul_add(db, vb);
            }
            i += 1;
        }
        (cov, va, vb)
    }
}

#[cfg(target_arch = "aarch64")]
mod neon {
    use num_complex::Complex64;
    use std::arch::aarch64::float64x2x2_t;
    use std::arch::aarch64::vabsq_f64;
    use std::arch::aarch64::vaddq_f64;
    use std::arch::aarch64::vaddvq_f64;
    use std::arch::aarch64::vdupq_n_f64;
    use std::arch::aarch64::vfmaq_f64;
    use std::arch::aarch64::vfmsq_f64;
    use std::arch::aarch64::vld1q_f64;
    use std::arch::aarch64::vld2q_f64;
    use std::arch::aarch64::vmaxq_f64;
    use std::arch::aarch64::vmaxvq_f64;
    use std::arch::aarch64::vminq_f64;
    use std::arch::aarch64::vminvq_f64;
    use std::arch::aarch64::vmulq_f64;
    use std::arch::aarch64::vst1q_f64;
    use std::arch::aarch64::vst2q_f64;
    use std::arch::aarch64::vsubq_f64;

    #[target_feature(enable = "neon")]
    pub(super) unsafe fn dot_f64(
        a: *const f64,
        b: *const f64,
        n: usize,
    ) -> f64 {
        let mut acc0 = vdupq_n_f64(0.0);
        let mut acc1 = vdupq_n_f64(0.0);
        let mut i = 0_usize;
        while i + 4 <= n {
            // SAFETY: `i + 4 <= n`.
            unsafe {
                let va0 = vld1q_f64(a.add(i));
                let vb0 = vld1q_f64(b.add(i));
                let va1 = vld1q_f64(a.add(i + 2));
                let vb1 = vld1q_f64(b.add(i + 2));
                acc0 = vfmaq_f64(acc0, va0, vb0);
                acc1 = vfmaq_f64(acc1, va1, vb1);
            }
            i += 4;
        }
        if i + 2 <= n {
            // SAFETY: `i + 2 <= n`.
            unsafe {
                acc0 = vfmaq_f64(acc0, vld1q_f64(a.add(i)), vld1q_f64(b.add(i)));
            }
            i += 2;
        }
        let mut total = vaddvq_f64(vaddq_f64(acc0, acc1));
        while i < n {
            // SAFETY: `i < n`.
            unsafe {
                total = (*a.add(i)).mul_add(*b.add(i), total);
            }
            i += 1;
        }
        total
    }

    #[target_feature(enable = "neon")]
    pub(super) unsafe fn sum_f64(
        x: *const f64,
        n: usize,
    ) -> f64 {
        let mut acc0 = vdupq_n_f64(0.0);
        let mut acc1 = vdupq_n_f64(0.0);
        let mut i = 0_usize;
        while i + 4 <= n {
            // SAFETY: `i + 4 <= n`.
            unsafe {
                acc0 = vaddq_f64(acc0, vld1q_f64(x.add(i)));
                acc1 = vaddq_f64(acc1, vld1q_f64(x.add(i + 2)));
            }
            i += 4;
        }
        if i + 2 <= n {
            // SAFETY: `i + 2 <= n`.
            unsafe {
                acc0 = vaddq_f64(acc0, vld1q_f64(x.add(i)));
            }
            i += 2;
        }
        let mut total = vaddvq_f64(vaddq_f64(acc0, acc1));
        while i < n {
            // SAFETY: `i < n`.
            unsafe {
                total += *x.add(i);
            }
            i += 1;
        }
        total
    }

    #[target_feature(enable = "neon")]
    pub(super) unsafe fn norm2_sq_f64(
        x: *const f64,
        n: usize,
    ) -> f64 {
        let mut acc0 = vdupq_n_f64(0.0);
        let mut acc1 = vdupq_n_f64(0.0);
        let mut i = 0_usize;
        while i + 4 <= n {
            // SAFETY: `i + 4 <= n`.
            unsafe {
                let v0 = vld1q_f64(x.add(i));
                let v1 = vld1q_f64(x.add(i + 2));
                acc0 = vfmaq_f64(acc0, v0, v0);
                acc1 = vfmaq_f64(acc1, v1, v1);
            }
            i += 4;
        }
        if i + 2 <= n {
            // SAFETY: `i + 2 <= n`.
            unsafe {
                let v = vld1q_f64(x.add(i));
                acc0 = vfmaq_f64(acc0, v, v);
            }
            i += 2;
        }
        let mut total = vaddvq_f64(vaddq_f64(acc0, acc1));
        while i < n {
            // SAFETY: `i < n`.
            unsafe {
                let v = *x.add(i);
                total = v.mul_add(v, total);
            }
            i += 1;
        }
        total
    }

    #[target_feature(enable = "neon")]
    pub(super) unsafe fn l1_norm_f64(
        x: *const f64,
        n: usize,
    ) -> f64 {
        let mut acc0 = vdupq_n_f64(0.0);
        let mut acc1 = vdupq_n_f64(0.0);
        let mut i = 0_usize;
        while i + 4 <= n {
            // SAFETY: `i + 4 <= n`.
            unsafe {
                acc0 = vaddq_f64(acc0, vabsq_f64(vld1q_f64(x.add(i))));
                acc1 = vaddq_f64(acc1, vabsq_f64(vld1q_f64(x.add(i + 2))));
            }
            i += 4;
        }
        if i + 2 <= n {
            // SAFETY: `i + 2 <= n`.
            unsafe {
                acc0 = vaddq_f64(acc0, vabsq_f64(vld1q_f64(x.add(i))));
            }
            i += 2;
        }
        let mut total = vaddvq_f64(vaddq_f64(acc0, acc1));
        while i < n {
            // SAFETY: `i < n`.
            unsafe {
                total += (*x.add(i)).abs();
            }
            i += 1;
        }
        total
    }

    #[target_feature(enable = "neon")]
    pub(super) unsafe fn linf_norm_f64(
        x: *const f64,
        n: usize,
    ) -> f64 {
        // ARMv8 `vmaxq_f64` / `vmaxvq_f64` (FMAX) propagate NaN directly.
        let mut vmax0 = vdupq_n_f64(0.0);
        let mut vmax1 = vdupq_n_f64(0.0);
        let mut i = 0_usize;
        while i + 4 <= n {
            // SAFETY: `i + 4 <= n`.
            unsafe {
                vmax0 = vmaxq_f64(vmax0, vabsq_f64(vld1q_f64(x.add(i))));
                vmax1 = vmaxq_f64(vmax1, vabsq_f64(vld1q_f64(x.add(i + 2))));
            }
            i += 4;
        }
        if i + 2 <= n {
            // SAFETY: `i + 2 <= n`.
            unsafe {
                vmax0 = vmaxq_f64(vmax0, vabsq_f64(vld1q_f64(x.add(i))));
            }
            i += 2;
        }
        let mut best = vmaxvq_f64(vmaxq_f64(vmax0, vmax1));
        if best.is_nan() {
            return f64::NAN;
        }
        while i < n {
            // SAFETY: `i < n`.
            unsafe {
                let a = (*x.add(i)).abs();
                if a.is_nan() {
                    return f64::NAN;
                }
                if a > best {
                    best = a;
                }
            }
            i += 1;
        }
        best
    }

    #[target_feature(enable = "neon")]
    pub(super) unsafe fn axpy_f64(
        alpha: f64,
        x: *const f64,
        y: *mut f64,
        n: usize,
    ) {
        let va = vdupq_n_f64(alpha);
        let mut i = 0_usize;
        while i + 4 <= n {
            // SAFETY: `i + 4 <= n`.
            unsafe {
                let r0 = vfmaq_f64(vld1q_f64(y.add(i)), va, vld1q_f64(x.add(i)));
                let r1 = vfmaq_f64(vld1q_f64(y.add(i + 2)), va, vld1q_f64(x.add(i + 2)));
                vst1q_f64(y.add(i), r0);
                vst1q_f64(y.add(i + 2), r1);
            }
            i += 4;
        }
        if i + 2 <= n {
            // SAFETY: `i + 2 <= n`.
            unsafe {
                let r = vfmaq_f64(vld1q_f64(y.add(i)), va, vld1q_f64(x.add(i)));
                vst1q_f64(y.add(i), r);
            }
            i += 2;
        }
        while i < n {
            // SAFETY: `i < n`.
            unsafe {
                *y.add(i) = alpha.mul_add(*x.add(i), *y.add(i));
            }
            i += 1;
        }
    }

    #[target_feature(enable = "neon")]
    pub(super) unsafe fn axpby_f64(
        alpha: f64,
        x: *const f64,
        beta: f64,
        y: *mut f64,
        n: usize,
    ) {
        let va = vdupq_n_f64(alpha);
        let vb = vdupq_n_f64(beta);
        let mut i = 0_usize;
        while i + 4 <= n {
            // SAFETY: `i + 4 <= n`.
            unsafe {
                let by0 = vmulq_f64(vb, vld1q_f64(y.add(i)));
                let by1 = vmulq_f64(vb, vld1q_f64(y.add(i + 2)));
                vst1q_f64(y.add(i), vfmaq_f64(by0, va, vld1q_f64(x.add(i))));
                vst1q_f64(y.add(i + 2), vfmaq_f64(by1, va, vld1q_f64(x.add(i + 2))));
            }
            i += 4;
        }
        if i + 2 <= n {
            // SAFETY: `i + 2 <= n`.
            unsafe {
                let by = vmulq_f64(vb, vld1q_f64(y.add(i)));
                vst1q_f64(y.add(i), vfmaq_f64(by, va, vld1q_f64(x.add(i))));
            }
            i += 2;
        }
        while i < n {
            // SAFETY: `i < n`.
            unsafe {
                *y.add(i) = alpha.mul_add(*x.add(i), beta * *y.add(i));
            }
            i += 1;
        }
    }

    #[target_feature(enable = "neon")]
    pub(super) unsafe fn scale_f64(
        alpha: f64,
        x: *mut f64,
        n: usize,
    ) {
        let va = vdupq_n_f64(alpha);
        let mut i = 0_usize;
        while i + 4 <= n {
            // SAFETY: `i + 4 <= n`.
            unsafe {
                vst1q_f64(x.add(i), vmulq_f64(va, vld1q_f64(x.add(i))));
                vst1q_f64(x.add(i + 2), vmulq_f64(va, vld1q_f64(x.add(i + 2))));
            }
            i += 4;
        }
        if i + 2 <= n {
            // SAFETY: `i + 2 <= n`.
            unsafe {
                vst1q_f64(x.add(i), vmulq_f64(va, vld1q_f64(x.add(i))));
            }
            i += 2;
        }
        while i < n {
            // SAFETY: `i < n`.
            unsafe {
                *x.add(i) *= alpha;
            }
            i += 1;
        }
    }

    #[target_feature(enable = "neon")]
    pub(super) unsafe fn add_into_f64(
        a: *const f64,
        b: *const f64,
        out: *mut f64,
        n: usize,
    ) {
        let mut i = 0_usize;
        while i + 4 <= n {
            // SAFETY: `i + 4 <= n`.
            unsafe {
                vst1q_f64(
                    out.add(i),
                    vaddq_f64(vld1q_f64(a.add(i)), vld1q_f64(b.add(i))),
                );
                vst1q_f64(
                    out.add(i + 2),
                    vaddq_f64(vld1q_f64(a.add(i + 2)), vld1q_f64(b.add(i + 2))),
                );
            }
            i += 4;
        }
        if i + 2 <= n {
            // SAFETY: `i + 2 <= n`.
            unsafe {
                vst1q_f64(
                    out.add(i),
                    vaddq_f64(vld1q_f64(a.add(i)), vld1q_f64(b.add(i))),
                );
            }
            i += 2;
        }
        while i < n {
            // SAFETY: `i < n`.
            unsafe {
                *out.add(i) = *a.add(i) + *b.add(i);
            }
            i += 1;
        }
    }

    #[target_feature(enable = "neon")]
    pub(super) unsafe fn sub_into_f64(
        a: *const f64,
        b: *const f64,
        out: *mut f64,
        n: usize,
    ) {
        let mut i = 0_usize;
        while i + 4 <= n {
            // SAFETY: `i + 4 <= n`.
            unsafe {
                vst1q_f64(
                    out.add(i),
                    vsubq_f64(vld1q_f64(a.add(i)), vld1q_f64(b.add(i))),
                );
                vst1q_f64(
                    out.add(i + 2),
                    vsubq_f64(vld1q_f64(a.add(i + 2)), vld1q_f64(b.add(i + 2))),
                );
            }
            i += 4;
        }
        if i + 2 <= n {
            // SAFETY: `i + 2 <= n`.
            unsafe {
                vst1q_f64(
                    out.add(i),
                    vsubq_f64(vld1q_f64(a.add(i)), vld1q_f64(b.add(i))),
                );
            }
            i += 2;
        }
        while i < n {
            // SAFETY: `i < n`.
            unsafe {
                *out.add(i) = *a.add(i) - *b.add(i);
            }
            i += 1;
        }
    }

    #[target_feature(enable = "neon")]
    pub(super) unsafe fn hadamard_mul_into_f64(
        a: *const f64,
        b: *const f64,
        out: *mut f64,
        n: usize,
    ) {
        let mut i = 0_usize;
        while i + 4 <= n {
            // SAFETY: `i + 4 <= n`.
            unsafe {
                vst1q_f64(
                    out.add(i),
                    vmulq_f64(vld1q_f64(a.add(i)), vld1q_f64(b.add(i))),
                );
                vst1q_f64(
                    out.add(i + 2),
                    vmulq_f64(vld1q_f64(a.add(i + 2)), vld1q_f64(b.add(i + 2))),
                );
            }
            i += 4;
        }
        if i + 2 <= n {
            // SAFETY: `i + 2 <= n`.
            unsafe {
                vst1q_f64(
                    out.add(i),
                    vmulq_f64(vld1q_f64(a.add(i)), vld1q_f64(b.add(i))),
                );
            }
            i += 2;
        }
        while i < n {
            // SAFETY: `i < n`.
            unsafe {
                *out.add(i) = *a.add(i) * *b.add(i);
            }
            i += 1;
        }
    }

    #[target_feature(enable = "neon")]
    pub(super) unsafe fn complex_mul_into(
        a: *const f64,
        b: *const f64,
        out: *mut f64,
        n: usize,
    ) {
        let total_f64 = n * 2;
        let mut i = 0_usize;
        while i + 4 <= total_f64 {
            // SAFETY: `i + 4 <= 2 * n`.
            unsafe {
                let va = vld2q_f64(a.add(i));
                let vb = vld2q_f64(b.add(i));
                let re = vfmsq_f64(vmulq_f64(va.0, vb.0), va.1, vb.1);
                let im = vfmaq_f64(vmulq_f64(va.0, vb.1), va.1, vb.0);
                vst2q_f64(out.add(i), float64x2x2_t(re, im));
            }
            i += 4;
        }
        while i < total_f64 {
            // SAFETY: `i + 1 < 2 * n`.
            unsafe {
                let ar = *a.add(i);
                let ai = *a.add(i + 1);
                let br = *b.add(i);
                let bi = *b.add(i + 1);
                *out.add(i) = ar * br - ai * bi;
                *out.add(i + 1) = ar * bi + ai * br;
            }
            i += 2;
        }
    }

    #[target_feature(enable = "neon")]
    pub(super) unsafe fn complex_mul_conj_into(
        a: *const f64,
        b: *const f64,
        out: *mut f64,
        n: usize,
    ) {
        let total_f64 = n * 2;
        let mut i = 0_usize;
        while i + 4 <= total_f64 {
            // SAFETY: `i + 4 <= 2 * n`.
            unsafe {
                let va = vld2q_f64(a.add(i));
                let vb = vld2q_f64(b.add(i));
                let re = vfmaq_f64(vmulq_f64(va.0, vb.0), va.1, vb.1);
                let im = vfmsq_f64(vmulq_f64(va.1, vb.0), va.0, vb.1);
                vst2q_f64(out.add(i), float64x2x2_t(re, im));
            }
            i += 4;
        }
        while i < total_f64 {
            // SAFETY: `i + 1 < 2 * n`.
            unsafe {
                let ar = *a.add(i);
                let ai = *a.add(i + 1);
                let br = *b.add(i);
                let bi = *b.add(i + 1);
                *out.add(i) = ar * br + ai * bi;
                *out.add(i + 1) = ai * br - ar * bi;
            }
            i += 2;
        }
    }

    #[target_feature(enable = "neon")]
    pub(super) unsafe fn complex_mul_in_place(
        a: *mut f64,
        b: *const f64,
        n: usize,
    ) {
        let total_f64 = n * 2;
        let mut i = 0_usize;
        while i + 4 <= total_f64 {
            // SAFETY: `i + 4 <= 2 * n`.
            unsafe {
                let va = vld2q_f64(a.add(i));
                let vb = vld2q_f64(b.add(i));
                let re = vfmsq_f64(vmulq_f64(va.0, vb.0), va.1, vb.1);
                let im = vfmaq_f64(vmulq_f64(va.0, vb.1), va.1, vb.0);
                vst2q_f64(a.add(i), float64x2x2_t(re, im));
            }
            i += 4;
        }
        while i < total_f64 {
            // SAFETY: `i + 1 < 2 * n`.
            unsafe {
                let ar = *a.add(i);
                let ai = *a.add(i + 1);
                let br = *b.add(i);
                let bi = *b.add(i + 1);
                *a.add(i) = ar * br - ai * bi;
                *a.add(i + 1) = ar * bi + ai * br;
            }
            i += 2;
        }
    }

    #[target_feature(enable = "neon")]
    pub(super) unsafe fn complex_dot_conj(
        a: *const f64,
        b: *const f64,
        n: usize,
    ) -> Complex64 {
        let total_f64 = n * 2;
        let mut acc_re = vdupq_n_f64(0.0);
        let mut acc_im = vdupq_n_f64(0.0);
        let mut i = 0_usize;
        while i + 4 <= total_f64 {
            // SAFETY: `i + 4 <= 2 * n`.
            unsafe {
                let va = vld2q_f64(a.add(i));
                let vb = vld2q_f64(b.add(i));
                acc_re = vfmaq_f64(vfmaq_f64(acc_re, va.0, vb.0), va.1, vb.1);
                acc_im = vfmsq_f64(vfmaq_f64(acc_im, va.0, vb.1), va.1, vb.0);
            }
            i += 4;
        }
        let mut re = vaddvq_f64(acc_re);
        let mut im = vaddvq_f64(acc_im);
        while i < total_f64 {
            // SAFETY: `i + 1 < 2 * n`.
            unsafe {
                let ar = *a.add(i);
                let ai = *a.add(i + 1);
                let br = *b.add(i);
                let bi = *b.add(i + 1);
                re = ar.mul_add(br, ai.mul_add(bi, re));
                im = ar.mul_add(bi, (-ai).mul_add(br, im));
            }
            i += 2;
        }
        Complex64::new(re, im)
    }

    #[target_feature(enable = "neon")]
    pub(super) unsafe fn min_max_f64(
        x: *const f64,
        n: usize,
    ) -> (f64, f64) {
        // SAFETY: caller guarantees `n >= 4`.
        let first = unsafe { vld1q_f64(x) };
        let mut vmin0 = first;
        let mut vmax0 = first;
        let mut vmin1 = first;
        let mut vmax1 = first;
        let mut i = 2_usize;
        while i + 4 <= n {
            // SAFETY: `i + 4 <= n`.
            unsafe {
                let v0 = vld1q_f64(x.add(i));
                let v1 = vld1q_f64(x.add(i + 2));
                vmin0 = vminq_f64(vmin0, v0);
                vmax0 = vmaxq_f64(vmax0, v0);
                vmin1 = vminq_f64(vmin1, v1);
                vmax1 = vmaxq_f64(vmax1, v1);
            }
            i += 4;
        }
        if i + 2 <= n {
            // SAFETY: `i + 2 <= n`.
            unsafe {
                let v = vld1q_f64(x.add(i));
                vmin0 = vminq_f64(vmin0, v);
                vmax0 = vmaxq_f64(vmax0, v);
            }
            i += 2;
        }
        let mut lo = vminvq_f64(vminq_f64(vmin0, vmin1));
        let mut hi = vmaxvq_f64(vmaxq_f64(vmax0, vmax1));
        if lo.is_nan() || hi.is_nan() {
            return (f64::NAN, f64::NAN);
        }
        while i < n {
            // SAFETY: `i < n`.
            unsafe {
                let v = *x.add(i);
                if v.is_nan() {
                    return (f64::NAN, f64::NAN);
                }
                if v < lo {
                    lo = v;
                }
                if v > hi {
                    hi = v;
                }
            }
            i += 1;
        }
        (lo, hi)
    }

    #[target_feature(enable = "neon")]
    pub(super) unsafe fn centered_stats_f64(
        a: *const f64,
        b: *const f64,
        n: usize,
        mean_a: f64,
        mean_b: f64,
    ) -> (f64, f64, f64) {
        let vma = vdupq_n_f64(mean_a);
        let vmb = vdupq_n_f64(mean_b);
        let mut cov0 = vdupq_n_f64(0.0);
        let mut va0 = vdupq_n_f64(0.0);
        let mut vb0 = vdupq_n_f64(0.0);
        let mut cov1 = vdupq_n_f64(0.0);
        let mut va1 = vdupq_n_f64(0.0);
        let mut vb1 = vdupq_n_f64(0.0);
        let mut i = 0_usize;
        while i + 4 <= n {
            // SAFETY: `i + 4 <= n`.
            unsafe {
                let da0 = vsubq_f64(vld1q_f64(a.add(i)), vma);
                let db0 = vsubq_f64(vld1q_f64(b.add(i)), vmb);
                let da1 = vsubq_f64(vld1q_f64(a.add(i + 2)), vma);
                let db1 = vsubq_f64(vld1q_f64(b.add(i + 2)), vmb);
                cov0 = vfmaq_f64(cov0, da0, db0);
                va0 = vfmaq_f64(va0, da0, da0);
                vb0 = vfmaq_f64(vb0, db0, db0);
                cov1 = vfmaq_f64(cov1, da1, db1);
                va1 = vfmaq_f64(va1, da1, da1);
                vb1 = vfmaq_f64(vb1, db1, db1);
            }
            i += 4;
        }
        if i + 2 <= n {
            // SAFETY: `i + 2 <= n`.
            unsafe {
                let da = vsubq_f64(vld1q_f64(a.add(i)), vma);
                let db = vsubq_f64(vld1q_f64(b.add(i)), vmb);
                cov0 = vfmaq_f64(cov0, da, db);
                va0 = vfmaq_f64(va0, da, da);
                vb0 = vfmaq_f64(vb0, db, db);
            }
            i += 2;
        }
        let mut cov = vaddvq_f64(vaddq_f64(cov0, cov1));
        let mut va = vaddvq_f64(vaddq_f64(va0, va1));
        let mut vb = vaddvq_f64(vaddq_f64(vb0, vb1));
        while i < n {
            // SAFETY: `i < n`.
            unsafe {
                let da = *a.add(i) - mean_a;
                let db = *b.add(i) - mean_b;
                cov = da.mul_add(db, cov);
                va = da.mul_add(da, va);
                vb = db.mul_add(db, vb);
            }
            i += 1;
        }
        (cov, va, vb)
    }
}

#[cfg(test)]
mod tests {
    use num_complex::Complex64;

    use super::*;

    fn sample_vec(
        n: usize,
        offset: f64,
    ) -> Vec<f64> {
        (0..n)
            .map(|i| ((i as f64 + offset) * 0.37).sin() * 3.25 - 0.5)
            .collect()
    }

    #[test]
    fn simd_matches_scalar_across_lengths() {
        for &n in &[
            0_usize, 1, 2, 3, 4, 5, 7, 8, 9, 15, 16, 17, 63, 64, 129, 513,
        ] {
            let a = sample_vec(n, 0.1);
            let b = sample_vec(n, 1.7);

            let tol = 1e-12 * (n as f64).max(1.0);
            assert!((dot_f64(&a, &b) - scalar::dot_f64(&a, &b)).abs() <= tol);
            assert!((sum_f64(&a) - scalar::sum_f64(&a)).abs() <= tol);
            assert!((norm2_sq_f64(&a) - scalar::norm2_sq_f64(&a)).abs() <= tol);
            assert!((l1_norm_f64(&a) - scalar::l1_norm_f64(&a)).abs() <= tol);
            assert!((linf_norm_f64(&a) - scalar::linf_norm_f64(&a)).abs() <= 1e-15);

            let mut y_simd = b.clone();
            let mut y_ref = b.clone();
            axpy_f64(-2.25, &a, &mut y_simd);
            scalar::axpy_f64(-2.25, &a, &mut y_ref);
            for (s, r) in y_simd.iter().zip(&y_ref) {
                assert!((s - r).abs() <= 1e-14);
            }

            let mut y2_simd = b.clone();
            let mut y2_ref = b.clone();
            axpby_f64(1.5, &a, -0.75, &mut y2_simd);
            scalar::axpby_f64(1.5, &a, -0.75, &mut y2_ref);
            for (s, r) in y2_simd.iter().zip(&y2_ref) {
                assert!((s - r).abs() <= 1e-14);
            }

            let mut sc_simd = a.clone();
            let mut sc_ref = a.clone();
            scale_f64(0.125, &mut sc_simd);
            scalar::scale_f64(0.125, &mut sc_ref);
            assert_eq!(sc_simd, sc_ref);

            let mut out_simd = vec![0.0; n];
            let mut out_ref = vec![0.0; n];
            add_into_f64(&a, &b, &mut out_simd);
            scalar::add_into_f64(&a, &b, &mut out_ref);
            assert_eq!(out_simd, out_ref);

            sub_into_f64(&a, &b, &mut out_simd);
            scalar::sub_into_f64(&a, &b, &mut out_ref);
            assert_eq!(out_simd, out_ref);

            hadamard_mul_into_f64(&a, &b, &mut out_simd);
            scalar::hadamard_mul_into_f64(&a, &b, &mut out_ref);
            assert_eq!(out_simd, out_ref);
        }
    }

    #[test]
    fn complex_simd_matches_scalar_across_lengths() {
        for &n in &[0_usize, 1, 2, 3, 4, 5, 7, 8, 17, 65] {
            let a: Vec<Complex64> = (0..n)
                .map(|i| Complex64::new((i as f64 * 0.3).sin(), (i as f64 * 0.7).cos()))
                .collect();
            let b: Vec<Complex64> = (0..n)
                .map(|i| Complex64::new((i as f64 * 0.5).cos(), (i as f64 * 0.2).sin()))
                .collect();

            let mut out_simd = vec![Complex64::new(0.0, 0.0); n];
            let mut out_ref = vec![Complex64::new(0.0, 0.0); n];
            complex_mul_into(&a, &b, &mut out_simd);
            scalar::complex_mul_into(&a, &b, &mut out_ref);
            for (s, r) in out_simd.iter().zip(&out_ref) {
                assert!((s - r).norm() <= 1e-14);
            }

            complex_mul_conj_into(&a, &b, &mut out_simd);
            scalar::complex_mul_conj_into(&a, &b, &mut out_ref);
            for (s, r) in out_simd.iter().zip(&out_ref) {
                assert!((s - r).norm() <= 1e-14);
            }

            let mut sc = a.clone();
            complex_scale_in_place(0.25, &mut sc);
            for (got, orig) in sc.iter().zip(&a) {
                assert!((*got - *orig * 0.25).norm() <= 1e-15);
            }
        }
    }

    #[test]
    fn linf_propagates_nan() {
        let mut x = sample_vec(16, 0.2);
        x[11] = f64::NAN;
        assert!(linf_norm_f64(&x).is_nan());
        let mut tail_nan = sample_vec(19, 0.2);
        tail_nan[18] = f64::NAN;
        assert!(linf_norm_f64(&tail_nan).is_nan());
    }

    #[test]
    fn extended_simd_primitives_match_scalar() {
        for &n in &[0_usize, 1, 2, 3, 4, 5, 7, 8, 9, 15, 16, 17, 65] {
            let a = sample_vec(n, 0.3);
            let b = sample_vec(n, 1.1);
            let tol = 1e-12 * (n as f64).max(1.0);

            if n == 0 {
                assert_eq!(min_max_f64(&a), None);
            } else {
                assert_eq!(min_max_f64(&a), Some(scalar::min_max_f64(&a)));
            }

            let ma = if n == 0 {
                0.0
            } else {
                sum_f64(&a) / n as f64
            };
            let mb = if n == 0 {
                0.0
            } else {
                sum_f64(&b) / n as f64
            };
            let (c1, va1, vb1) = centered_stats_f64(&a, &b, ma, mb);
            let (c2, va2, vb2) = scalar::centered_stats_f64(&a, &b, ma, mb);
            assert!((c1 - c2).abs() <= tol);
            assert!((va1 - va2).abs() <= tol);
            assert!((vb1 - vb2).abs() <= tol);

            let ca: Vec<Complex64> = (0..n)
                .map(|i| Complex64::new((i as f64 * 0.4).sin(), (i as f64 * 0.6).cos()))
                .collect();
            let cb: Vec<Complex64> = (0..n)
                .map(|i| Complex64::new((i as f64 * 0.7).cos(), (i as f64 * 0.2).sin()))
                .collect();

            let mut ip_simd = ca.clone();
            let mut ip_ref = ca.clone();
            complex_mul_in_place(&mut ip_simd, &cb);
            scalar::complex_mul_in_place(&mut ip_ref, &cb);
            for (s, r) in ip_simd.iter().zip(&ip_ref) {
                assert!((s - r).norm() <= 1e-14);
            }

            let d_simd = complex_dot_conj(&ca, &cb);
            let d_ref = scalar::complex_dot_conj(&ca, &cb);
            assert!((d_simd - d_ref).norm() <= tol);

            let n2_simd = complex_norm2_sq(&ca);
            let n2_ref: f64 = ca.iter().map(Complex64::norm_sqr).sum();
            assert!((n2_simd - n2_ref).abs() <= tol);
        }

        let mut with_nan = sample_vec(19, 0.5);
        with_nan[10] = f64::NAN;
        let (lo, hi) = min_max_f64(&with_nan).unwrap_or((0.0, 0.0));
        assert!(lo.is_nan() && hi.is_nan());
    }
}
