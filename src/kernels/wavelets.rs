//! # Discrete wavelet transforms and the discrete fractional Fourier transform
//!
//! * [`Wavelet`]: Haar and Daubechies `db2`, `db3`, `db4` orthogonal filter banks
//!   ([`Wavelet::low_pass`], [`Wavelet::high_pass`]);
//! * [`dwt`] / [`idwt`]: one level of the periodised orthogonal transform (a signal of even
//!   length `n` maps to `n/2` approximation and `n/2` detail coefficients, and back, exactly);
//! * [`wavedec`] / [`waverec`]: the multilevel decomposition and its inverse;
//! * [`frft_matrix`] / [`frft`]: the **discrete fractional Fourier transform** of Candan, Kutay
//!   and Ozaktas / Pei: the eigenvectors of the real symmetric matrix `S` that commutes with the
//!   unitary DFT are ordered like the Hermite functions and raised to the fractional power
//!   `exp(-i pi a k / 2)`; the result is unitary, `a = 1` is the unitary DFT, `a = 0` the
//!   identity, `a = 2` the reversal `x[n] -> x[-n mod N]`, and the transforms form a group in
//!   `a`.
//!
//! Limitations: the wavelet family is limited to `Haar` and `db2`..`db4` (no biorthogonal or
//! symmetric families, no 2-D transform, periodic extension only); the fractional Fourier matrix
//! costs an `O(N^3)` symmetric eigensolve (use it for lengths up to a few hundred; the continuous
//! transform of a sampled chirp-like function needs the sampling rate chosen by the caller).

use num_complex::Complex64;

use crate::kernels::dense::Mat;
use crate::kernels::dense::eigen_symmetric;

/// An orthogonal wavelet family.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Wavelet {
    /// Haar (`db1`), 1 vanishing moment.
    Haar,
    /// Daubechies with 2 vanishing moments (4 taps).
    Db2,
    /// Daubechies with 3 vanishing moments (6 taps).
    Db3,
    /// Daubechies with 4 vanishing moments (8 taps).
    Db4,
}

impl Wavelet {
    /// The scaling (low-pass) filter `h`, normalised to `sum h = sqrt 2`, `sum h^2 = 1`.
    #[must_use]
    pub fn low_pass(self) -> Vec<f64> {
        match self {
            | Self::Haar => vec![std::f64::consts::FRAC_1_SQRT_2; 2],
            | Self::Db2 => {
                let s3 = 3.0_f64.sqrt();
                let d = 4.0 * std::f64::consts::SQRT_2;
                vec![
                    (1.0 + s3) / d,
                    (3.0 + s3) / d,
                    (3.0 - s3) / d,
                    (1.0 - s3) / d,
                ]
            },
            | Self::Db3 => {
                vec![
                    0.332_670_552_950_082_5,
                    0.806_891_509_311_092_4,
                    0.459_877_502_118_491_4,
                    -0.135_011_020_010_254_6,
                    -0.085_441_273_882_026_7,
                    0.035_226_291_885_709_5,
                ]
            },
            | Self::Db4 => {
                vec![
                    0.230_377_813_308_896_4,
                    0.714_846_570_552_915_4,
                    0.630_880_767_929_858_7,
                    -0.027_983_769_416_859_9,
                    -0.187_034_811_719_093_1,
                    0.030_841_381_835_560_7,
                    0.032_883_011_666_885_2,
                    -0.010_597_401_785_069_0,
                ]
            },
        }
    }

    /// The wavelet (high-pass) filter `g[k] = (-1)^k h[L - 1 - k]`.
    #[must_use]
    pub fn high_pass(self) -> Vec<f64> {
        let h = self.low_pass();
        let l = h.len();
        (0..l)
            .map(|k| {
                if k % 2 == 0 {
                    h[l - 1 - k]
                } else {
                    -h[l - 1 - k]
                }
            })
            .collect()
    }

    /// Number of filter taps.
    #[must_use]
    pub fn taps(self) -> usize {
        self.low_pass().len()
    }
}

/// Errors of the wavelet transforms.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum WaveletError {
    /// The signal length is odd or zero.
    OddLength,
    /// More levels were asked for than the length allows.
    TooManyLevels,
    /// The coefficient lists do not fit together.
    Shape,
}

/// One level of the periodised discrete wavelet transform: returns `(approximation, detail)`.
///
/// `a[k] = sum_m h[m] x[(2k + m) mod n]`, `d[k] = sum_m g[m] x[(2k + m) mod n]`.
///
/// # Errors
/// [`WaveletError::OddLength`] for an empty or odd-length signal.
pub fn dwt(
    x: &[f64],
    w: Wavelet,
) -> Result<(Vec<f64>, Vec<f64>), WaveletError> {
    let n = x.len();
    if n == 0 || !n.is_multiple_of(2) {
        return Err(WaveletError::OddLength);
    }
    let (h, g) = (w.low_pass(), w.high_pass());
    let half = n / 2;
    let mut a = vec![0.0; half];
    let mut d = vec![0.0; half];
    for k in 0..half {
        for (m, (hm, gm)) in h.iter().zip(&g).enumerate() {
            let v = x[(2 * k + m) % n];
            a[k] += hm * v;
            d[k] += gm * v;
        }
    }
    Ok((a, d))
}

/// Inverse of [`dwt`] (exact up to rounding).
///
/// # Errors
/// [`WaveletError::Shape`] when the two coefficient lists differ in length or are empty.
pub fn idwt(
    a: &[f64],
    d: &[f64],
    w: Wavelet,
) -> Result<Vec<f64>, WaveletError> {
    if a.len() != d.len() || a.is_empty() {
        return Err(WaveletError::Shape);
    }
    let (h, g) = (w.low_pass(), w.high_pass());
    let n = 2 * a.len();
    let mut x = vec![0.0; n];
    // the analysis matrix is orthogonal, so synthesis is its transpose
    for k in 0..a.len() {
        for (m, (hm, gm)) in h.iter().zip(&g).enumerate() {
            x[(2 * k + m) % n] += hm * a[k] + gm * d[k];
        }
    }
    Ok(x)
}

/// Multilevel decomposition: `[approximation_L, detail_L, ..., detail_1]`.
///
/// # Errors
/// [`WaveletError::TooManyLevels`] when `2^levels` does not divide the length, or the length is
/// shorter than the filter at some level.
pub fn wavedec(
    x: &[f64],
    w: Wavelet,
    levels: usize,
) -> Result<Vec<Vec<f64>>, WaveletError> {
    if x.is_empty() || levels == 0 || !x.len().is_multiple_of(1_usize << levels.min(60)) {
        return Err(WaveletError::TooManyLevels);
    }
    let mut details = Vec::new();
    let mut current = x.to_vec();
    for _ in 0..levels {
        if current.len() < w.taps() {
            return Err(WaveletError::TooManyLevels);
        }
        let (a, d) = dwt(&current, w)?;
        details.push(d);
        current = a;
    }
    let mut out = vec![current];
    out.extend(details.into_iter().rev());
    Ok(out)
}

/// Inverse of [`wavedec`].
///
/// # Errors
/// [`WaveletError::Shape`] for an inconsistent list.
pub fn waverec(
    coeffs: &[Vec<f64>],
    w: Wavelet,
) -> Result<Vec<f64>, WaveletError> {
    let (first, rest) = coeffs.split_first().ok_or(WaveletError::Shape)?;
    let mut current = first.clone();
    for d in rest {
        current = idwt(&current, d, w)?;
    }
    Ok(current)
}

/// The unitary DFT matrix `F[j][k] = exp(-2 pi i j k / N) / sqrt N`.
fn dft_apply(
    n: usize,
    v: &[f64],
) -> Vec<Complex64> {
    let s = 1.0 / (n as f64).sqrt();
    (0..n)
        .map(|j| {
            let mut acc = Complex64::new(0.0, 0.0);
            for (k, x) in v.iter().enumerate() {
                let ang = -2.0 * std::f64::consts::PI * ((j * k) % n) as f64 / n as f64;
                acc += Complex64::from_polar(*x, ang);
            }
            acc * s
        })
        .collect()
}

/// The `n x n` matrix of the discrete fractional Fourier transform of order `a` (row-major).
///
/// `a = 1` is the unitary DFT; the matrix is unitary and `frft_matrix(n, a) * frft_matrix(n, b)
/// = frft_matrix(n, a + b)`.
///
/// # Errors
/// `Err` for `n < 2` or when the eigensolver fails.
pub fn frft_matrix(
    n: usize,
    a: f64,
) -> Result<Vec<Complex64>, String> {
    if n < 2 {
        return Err("the fractional Fourier transform needs n >= 2".into());
    }
    // S commutes with the DFT; its eigenvectors approximate the Hermite-Gauss functions
    let mut s = Mat::zeros(n, n);
    for i in 0..n {
        s.set(
            i,
            i,
            2.0 * (2.0 * std::f64::consts::PI * i as f64 / n as f64).cos() - 4.0,
        );
        let (p, q) = ((i + 1) % n, (i + n - 1) % n);
        s.set(i, p, s.at(i, p) + 1.0);
        s.set(i, q, s.at(i, q) + 1.0);
    }
    let (theta, vecs) = eigen_symmetric(&s).map_err(|e| e.to_string())?;
    // Hermite order = descending eigenvalue of S
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&p, &q| theta[q].total_cmp(&theta[p]));
    // Eigenvectors of S with (numerically) equal eigenvalues may mix DFT eigenspaces; each such
    // cluster is rotated to pure DFT eigenvectors: the restriction of the DFT to the cluster is a
    // unitary *symmetric* matrix M = A + iB with commuting real symmetric A, B, whose common
    // real eigenvectors are the pure combinations.
    let scale = theta.iter().fold(1.0_f64, |m, t| m.max(t.abs()));
    let mut vectors: Vec<(Vec<f64>, usize, usize)> = Vec::with_capacity(n); // (vector, class, rank)
    let mut start = 0;
    while start < n {
        let mut end = start + 1;
        while end < n && (theta[order[end - 1]] - theta[order[end]]).abs() < 1e-8 * scale {
            end += 1;
        }
        let cols: Vec<Vec<f64>> = order[start..end]
            .iter()
            .map(|&c| (0..n).map(|i| vecs.at(i, c)).collect())
            .collect();
        let m = cols.len();
        let pure: Vec<Vec<f64>> = if m == 1 {
            cols
        } else {
            let fcols: Vec<Vec<Complex64>> = cols.iter().map(|u| dft_apply(n, u)).collect();
            let mut r = Mat::zeros(m, m);
            for i in 0..m {
                for j in 0..m {
                    let mij: Complex64 = fcols[j].iter().zip(&cols[i]).map(|(f, x)| f * x).sum();
                    r.set(i, j, mij.re + 0.618_033_988_75 * mij.im);
                }
            }
            for i in 0..m {
                for j in 0..i {
                    let v = f64::midpoint(r.at(i, j), r.at(j, i));
                    r.set(i, j, v);
                    r.set(j, i, v);
                }
            }
            let (_, w) = eigen_symmetric(&r).map_err(|e| e.to_string())?;
            (0..m)
                .map(|c| {
                    (0..n)
                        .map(|i| (0..m).map(|j| w.at(j, c) * cols[j][i]).sum::<f64>())
                        .collect()
                })
                .collect()
        };
        for (offset, u) in pure.into_iter().enumerate() {
            let fu = dft_apply(n, &u);
            let lam: Complex64 = fu.iter().zip(&u).map(|(f, x)| f * x).sum();
            let class = [
                Complex64::new(1.0, 0.0),
                Complex64::new(0.0, -1.0),
                Complex64::new(-1.0, 0.0),
                Complex64::new(0.0, 1.0),
            ]
            .iter()
            .enumerate()
            .max_by(|(_, p), (_, q)| (lam * p.conj()).re.total_cmp(&(lam * q.conj()).re))
            .map_or(0, |(c, _)| c);
            vectors.push((u, class, start + offset));
        }
        start = end;
    }
    // Hermite index of each vector: a free index with the right class modulo 4, nearest to its
    // rank (the missing index of an even length is skipped this way)
    let mut used = vec![false; n + 1];
    vectors.sort_by_key(|v| v.2);
    let mut mat = vec![Complex64::new(0.0, 0.0); n * n];
    for (u, class, rank) in &vectors {
        let k = (0..=n)
            .filter(|k| k % 4 == *class && !used[*k])
            .min_by_key(|k| k.abs_diff(*rank))
            .unwrap_or(*rank);
        used[k.min(n)] = true;
        let phase = Complex64::from_polar(1.0, -std::f64::consts::FRAC_PI_2 * a * k as f64);
        for i in 0..n {
            for j in 0..n {
                mat[i * n + j] += phase * (u[i] * u[j]);
            }
        }
    }
    Ok(mat)
}

/// Discrete fractional Fourier transform of order `a` of the complex signal `x`.
///
/// # Errors
/// As [`frft_matrix`].
pub fn frft(
    x: &[Complex64],
    a: f64,
) -> Result<Vec<Complex64>, String> {
    let n = x.len();
    let m = frft_matrix(n, a)?;
    Ok((0..n)
        .map(|i| (0..n).map(|j| m[i * n + j] * x[j]).sum())
        .collect())
}
