//! # FFT and fast convolution
//!
//! Complex FFT for arbitrary lengths (radix-2 for powers of two,
//! Bluestein's chirp-z algorithm otherwise), real convolution and
//! correlation through the FFT, and exact-ish polynomial multiplication.
#![allow(
    clippy::manual_midpoint,
    clippy::missing_const_for_fn,
    clippy::struct_field_names,
    clippy::or_fun_call,
    clippy::manual_swap,
    clippy::if_not_else,
    clippy::unnecessary_sort_by,
    clippy::while_float,
    clippy::too_long_first_doc_paragraph,
    clippy::cast_sign_loss,
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::needless_pass_by_value,
    clippy::manual_map,
    clippy::unnecessary_map_or,
    clippy::suboptimal_flops,
    clippy::similar_names,
    clippy::unreadable_literal,
    clippy::excessive_precision,
    clippy::needless_range_loop,
    clippy::float_cmp,
    clippy::too_many_lines,
    clippy::cognitive_complexity,
    clippy::option_if_let_else,
    clippy::many_single_char_names,
    clippy::type_complexity,
    clippy::too_many_arguments,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::cast_precision_loss,
    clippy::map_unwrap_or,
    clippy::cast_lossless
)]

use std::f64::consts::PI;

/// Complex number as `(re, im)`.
pub type Complex = (f64, f64);

fn cm(a: Complex, b: Complex) -> Complex {
    (a.0 * b.0 - a.1 * b.1, a.0 * b.1 + a.1 * b.0)
}

fn fft_pow2(a: &mut [Complex], inverse: bool) {
    let n = a.len();
    let mut j = 0;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if i < j {
            a.swap(i, j);
        }
    }
    let mut len = 2;
    while len <= n {
        let ang = 2.0 * PI / len as f64 * if inverse { 1.0 } else { -1.0 };
        for i in (0..n).step_by(len) {
            for k in 0..len / 2 {
                let w = ((ang * k as f64).cos(), (ang * k as f64).sin());
                let u = a[i + k];
                let v = cm(a[i + k + len / 2], w);
                a[i + k] = (u.0 + v.0, u.1 + v.1);
                a[i + k + len / 2] = (u.0 - v.0, u.1 - v.1);
            }
        }
        len <<= 1;
    }
}

/// Discrete Fourier transform of any length (`inverse` includes the `1/n` scaling).
#[must_use]
pub fn fft(x: &[Complex], inverse: bool) -> Vec<Complex> {
    let n = x.len();
    if n <= 1 {
        return x.to_vec();
    }
    if n.is_power_of_two() {
        let mut a = x.to_vec();
        fft_pow2(&mut a, inverse);
        if inverse {
            for v in &mut a {
                *v = (v.0 / n as f64, v.1 / n as f64);
            }
        }
        return a;
    }
    // Bluestein: X_k = conj(w_k) sum (x_j conj(w_j)) w_{k-j}, w_j = e^{i pi j^2 / n}
    let sign = if inverse { 1.0 } else { -1.0 };
    let chirp: Vec<Complex> = (0..n)
        .map(|j| {
            let jj = (j * j) % (2 * n);
            let ang = sign * PI * jj as f64 / n as f64;
            (ang.cos(), ang.sin())
        })
        .collect();
    let m = (2 * n - 1).next_power_of_two();
    let mut a = vec![(0.0, 0.0); m];
    let mut b = vec![(0.0, 0.0); m];
    for j in 0..n {
        a[j] = cm(x[j], chirp[j]);
    }
    b[0] = (chirp[0].0, -chirp[0].1);
    for j in 1..n {
        let c = (chirp[j].0, -chirp[j].1);
        b[j] = c;
        b[m - j] = c;
    }
    fft_pow2(&mut a, false);
    fft_pow2(&mut b, false);
    for i in 0..m {
        a[i] = cm(a[i], b[i]);
    }
    fft_pow2(&mut a, true);
    (0..n)
        .map(|k| {
            let v = cm((a[k].0 / m as f64, a[k].1 / m as f64), chirp[k]);
            if inverse { (v.0 / n as f64, v.1 / n as f64) } else { v }
        })
        .collect()
}

/// Linear convolution of two real sequences (`len = a.len() + b.len() - 1`).
#[must_use]
pub fn convolve(a: &[f64], b: &[f64]) -> Vec<f64> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let n = a.len() + b.len() - 1;
    if a.len().min(b.len()) <= 16 {
        let mut out = vec![0.0; n];
        for (i, &x) in a.iter().enumerate() {
            for (j, &y) in b.iter().enumerate() {
                out[i + j] += x * y;
            }
        }
        return out;
    }
    let m = n.next_power_of_two();
    let mut fa: Vec<Complex> = (0..m).map(|i| (a.get(i).copied().unwrap_or(0.0), 0.0)).collect();
    let mut fb: Vec<Complex> = (0..m).map(|i| (b.get(i).copied().unwrap_or(0.0), 0.0)).collect();
    fft_pow2(&mut fa, false);
    fft_pow2(&mut fb, false);
    for i in 0..m {
        fa[i] = cm(fa[i], fb[i]);
    }
    fft_pow2(&mut fa, true);
    fa.iter().take(n).map(|v| v.0 / m as f64).collect()
}

/// Circular convolution of two equal-length real sequences.
#[must_use]
pub fn circular_convolve(a: &[f64], b: &[f64]) -> Vec<f64> {
    let n = a.len();
    let fa = fft(&a.iter().map(|&v| (v, 0.0)).collect::<Vec<_>>(), false);
    let fb = fft(&b.iter().map(|&v| (v, 0.0)).collect::<Vec<_>>(), false);
    let prod: Vec<Complex> = (0..n).map(|i| cm(fa[i], fb[i])).collect();
    fft(&prod, true).iter().map(|v| v.0).collect()
}

/// Cross-correlation `c[k] = sum_i a[i + k - (len(b) - 1)] * b[i]`
/// (full mode, length `a.len() + b.len() - 1`).
#[must_use]
pub fn correlate(a: &[f64], b: &[f64]) -> Vec<f64> {
    let rb: Vec<f64> = b.iter().rev().copied().collect();
    convolve(a, &rb)
}

/// Multiplies two integer polynomials (ascending coefficients) exactly by
/// splitting into 15-bit limbs so the floating-point FFT stays exact for
/// moderate sizes.
#[must_use]
pub fn poly_mul_i64(a: &[i64], b: &[i64]) -> Vec<i64> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let split = |v: &[i64]| -> (Vec<f64>, Vec<f64>) {
        (
            v.iter().map(|&x| (x & 0x7fff) as f64).collect(),
            v.iter().map(|&x| (x >> 15) as f64).collect(),
        )
    };
    let (alo, ahi) = split(a);
    let (blo, bhi) = split(b);
    let ll = convolve(&alo, &blo);
    let lh = convolve(&alo, &bhi);
    let hl = convolve(&ahi, &blo);
    let hh = convolve(&ahi, &bhi);
    (0..ll.len())
        .map(|i| {
            (ll[i].round() as i64)
                + (((lh[i] + hl[i]).round() as i64) << 15)
                + ((hh[i].round() as i64) << 30)
        })
        .collect()
}
