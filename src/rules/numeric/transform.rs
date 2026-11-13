//! Multilevel and two-dimensional wavelet transforms and the fractional
//! Fourier transform.
//!
//! | operator | value |
//! |---|---|
//! | `wavedec(x, wavelet, levels)` | `list(a_L, d_L, ..., d_1)`: the final approximation, then the details from the coarsest level to the finest (periodised; `2^levels` must divide the length) |
//! | `waverec(coeffs, wavelet)` | the signal from the list `wavedec` returns |
//! | `dwt2(m, wavelet)` | `list(LL, LH, HL, HH)`, the four sub-bands of one level of the separable transform of a matrix with an even number of rows and columns; `LH` is low along the rows' index and high along the columns' |
//! | `idwt2(LL, LH, HL, HH, wavelet)` | the matrix the four sub-bands came from |
//! | `frft_re(re, im, a)`, `frft_im(re, im, a)` | the real and imaginary parts of the discrete fractional Fourier transform of order `a` of the signal `re + i im` (`a = 1` is the unitary DFT; the signal has 2 to 128 samples) |
//!
//! The wavelet is a symbol, as for `dwt_approx`. There is one operator for
//! the matrix transform whatever the wavelet and no `_2d` suffix: a vector
//! goes to `dwt_approx`/`dwt_detail`, a matrix to `dwt2`.

use num_complex::Complex64;

use super::bridge::Data;
use super::bridge::Entry;
use super::bridge::MAX_LEN;
use super::bridge::count;
use super::bridge::matrix;
use super::bridge::rows;
use super::bridge::scalar;
use super::bridge::vec_out;
use super::bridge::vector;
use super::bridge::wavelet;
use crate::kernels::wavelets;
use crate::kernels::wavelets::Boundary;

/// Most elements of a matrix transformed by `dwt2`.
const MAX_MATRIX: usize = 1_000_000;
/// Most samples of a fractional Fourier transform (it diagonalises an `n x n` matrix).
const MAX_FRFT: usize = 128;

fn wavedec(a: &[Data]) -> Option<Data> {
    let (x, w) = (vector(a, 0)?, wavelet(a, 1)?);
    let levels = count(a, 2, 30)?;
    (levels >= 1 && !x.is_empty()).then_some(())?;
    let coefficients = wavelets::wavedec(x, w, levels).ok()?;
    Some(Data::Tree(
        coefficients.into_iter().map(Data::vector).collect(),
    ))
}

fn waverec(a: &[Data]) -> Option<Data> {
    let coefficients = rows(a.first()?)?;
    let w = wavelet(a, 1)?;
    (coefficients.len() >= 2 && coefficients.iter().all(|c| c.len() <= MAX_LEN)).then_some(())?;
    vec_out(wavelets::waverec(&coefficients, w).ok()?)
}

/// A matrix as a `rows x cols` array.
fn matrix_out(
    rows: usize,
    cols: usize,
    data: Vec<f64>,
) -> Data {
    Data::Array {
        shape: vec![rows, cols],
        data,
    }
}

/// Rows `r0..r1` and columns `c0..c1` of the row-major `m` (with `width` columns).
fn block(
    m: &[f64],
    width: usize,
    (r0, r1): (usize, usize),
    (c0, c1): (usize, usize),
) -> Vec<f64> {
    (r0..r1)
        .flat_map(|r| m[r * width + c0..r * width + c1].iter().copied())
        .collect()
}

fn dwt2(a: &[Data]) -> Option<Data> {
    let (rows, cols, x) = matrix(a, 0)?;
    let w = wavelet(a, 1)?;
    (rows % 2 == 0 && cols % 2 == 0 && rows * cols <= MAX_MATRIX).then_some(())?;
    let packed = wavelets::dwt2(x, rows, cols, w, Boundary::Periodic).ok()?;
    let (ra, ca) = (rows / 2, cols / 2);
    let width = 2 * ca;
    (packed.len() == 4 * ra * ca).then_some(())?;
    let band =
        |r: (usize, usize), c: (usize, usize)| matrix_out(ra, ca, block(&packed, width, r, c));
    Some(Data::Tree(vec![
        band((0, ra), (0, ca)),
        band((0, ra), (ca, width)),
        band((ra, 2 * ra), (0, ca)),
        band((ra, 2 * ra), (ca, width)),
    ]))
}

fn idwt2(a: &[Data]) -> Option<Data> {
    let w = wavelet(a, 4)?;
    let (ra, ca, ll) = matrix(a, 0)?;
    let mut bands = vec![ll];
    for at in 1..4 {
        let (r, c, band) = matrix(a, at)?;
        (r == ra && c == ca).then_some(())?;
        bands.push(band);
    }
    (4 * ra * ca <= MAX_MATRIX).then_some(())?;
    let (rows, cols) = (2 * ra, 2 * ca);
    let mut packed = vec![0.0; rows * cols];
    for (k, band) in bands.iter().enumerate() {
        let (r0, c0) = ((k / 2) * ra, (k % 2) * ca);
        for r in 0..ra {
            let dst = (r0 + r) * cols + c0;
            packed[dst..dst + ca].copy_from_slice(&band[r * ca..(r + 1) * ca]);
        }
    }
    let out = wavelets::idwt2(&packed, rows, cols, w, Boundary::Periodic).ok()?;
    Some(matrix_out(rows, cols, out))
}

/// The transform of the signal `re + i im` of order `a`.
fn fractional(a: &[Data]) -> Option<Vec<Complex64>> {
    let (re, im) = (vector(a, 0)?, vector(a, 1)?);
    let order = scalar(a, 2)?;
    (re.len() == im.len() && (2..=MAX_FRFT).contains(&re.len()) && order.abs() <= 1.0e6)
        .then_some(())?;
    let signal: Vec<Complex64> = re
        .iter()
        .zip(im)
        .map(|(&r, &i)| Complex64::new(r, i))
        .collect();
    wavelets::frft(&signal, order).ok()
}

fn frft_re(a: &[Data]) -> Option<Data> {
    vec_out(fractional(a)?.iter().map(|c| c.re).collect())
}

fn frft_im(a: &[Data]) -> Option<Data> {
    vec_out(fractional(a)?.iter().map(|c| c.im).collect())
}

pub(super) const ENTRIES: &[Entry] = &[
    Entry {
        name: "wavedec",
        params: "vws",
        eval: wavedec,
    },
    Entry {
        name: "waverec",
        params: "nw",
        eval: waverec,
    },
    Entry {
        name: "dwt2",
        params: "nw",
        eval: dwt2,
    },
    Entry {
        name: "idwt2",
        params: "nnnnw",
        eval: idwt2,
    },
    Entry {
        name: "frft_re",
        params: "vvs",
        eval: frft_re,
    },
    Entry {
        name: "frft_im",
        params: "vvs",
        eval: frft_im,
    },
];

#[cfg(test)]
mod tests {
    use crate::rules::numeric::numeric;
    use crate::rules::testing::reduce_with;
    use crate::rules::testing::simplify;

    fn run(src: &str) -> String {
        simplify(&[numeric()], src)
    }

    fn flat(text: &str) -> Vec<f64> {
        text.replace(['(', ')', ','], " ")
            .split_whitespace()
            .filter(|word| *word != "list")
            .map(|word| word.parse().unwrap_or(f64::NAN))
            .collect()
    }

    fn close(
        got: &[f64],
        want: &[f64],
        tol: f64,
    ) {
        assert_eq!(got.len(), want.len(), "{got:?} vs {want:?}");
        for (g, w) in got.iter().zip(want) {
            assert!((g - w).abs() <= tol, "{got:?} vs {want:?}");
        }
    }

    #[test]
    fn multilevel_decomposition_has_the_pywt_layout_and_inverts() {
        // Two Haar levels of a constant: the approximation is 2 x constant (sqrt2^2),
        // every detail is zero.
        let out = run("wavedec(list(1, 1, 1, 1, 1, 1, 1, 1), haar, 2)");
        assert_eq!(out.matches("list(").count(), 4, "{out}");
        close(
            &flat(&out),
            &[2.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            1e-12,
        );
        // The first level agrees with dwt_approx/dwt_detail, a second implementation.
        let x = "list(3, 1, 4, 1, 5, 9, 2, 6)";
        let level_one = flat(&run(&format!("wavedec({x}, db2, 1)")));
        let a = flat(&run(&format!("dwt_approx({x}, db2)")));
        let d = flat(&run(&format!("dwt_detail({x}, db2)")));
        close(&level_one, &[a, d].concat(), 1e-12);
        // A filter longer than the current signal stops the cascade, so the
        // longer wavelets get one level of this eight-sample signal.
        for (name, deepest) in [("haar", 3), ("db2", 2), ("sym4", 1), ("bior2_2", 1)] {
            for levels in 1..=deepest {
                let round = run(&format!("waverec(wavedec({x}, {name}, {levels}), {name})"));
                close(&flat(&round), &flat(x), 1e-9);
            }
        }
        // Parseval for an orthogonal wavelet: the energy is kept.
        let energy: f64 = flat(&run(&format!("wavedec({x}, db2, 2)")))
            .iter()
            .map(|v| v * v)
            .sum();
        assert!((energy - 173.0).abs() < 1e-9, "{energy}");
    }

    #[test]
    fn matrix_transform_splits_into_four_bands_and_inverts() {
        let m =
            "list(list(1, 2, 3, 4), list(5, 6, 7, 8), list(9, 10, 11, 12), list(13, 14, 15, 16))";
        let bands = run(&format!("dwt2({m}, haar)"));
        assert_eq!(flat(&bands).len(), 16, "{bands}");
        // The Haar LL band is the 2x2 block sum halved: (1 + 2 + 5 + 6) / 2 = 7.
        assert!((flat(&bands)[0] - 7.0).abs() < 1e-12, "{bands}");
        // Energy is kept by the orthogonal transform.
        let energy: f64 = flat(&bands).iter().map(|v| v * v).sum();
        let input: f64 = (1..=16).map(|v| f64::from(v * v)).sum();
        assert!((energy - input).abs() < 1e-9);
        for name in ["haar", "db2"] {
            let numbers = flat(&run(&format!("dwt2({m}, {name})")));
            let band = |k: usize| {
                let v = &numbers[4 * k..4 * k + 4];
                format!("list(list({}, {}), list({}, {}))", v[0], v[1], v[2], v[3])
            };
            let round = run(&format!(
                "idwt2({}, {}, {}, {}, {name})",
                band(0),
                band(1),
                band(2),
                band(3)
            ));
            close(&flat(&round), &flat(m), 1e-9);
        }
    }

    #[test]
    fn fractional_fourier_transform_has_the_dft_as_order_one() {
        // Order 0 is the identity and order 1 is the unitary DFT, i.e. fft / sqrt(n).
        let x = "list(1, 2, 0, -1, 3, 0.5, -2, 1)";
        let zero = flat(&run(&format!(
            "frft_re({x}, list(0, 0, 0, 0, 0, 0, 0, 0), 0)"
        )));
        close(&zero, &flat(x), 1e-9);
        let n = 8.0_f64.sqrt();
        let want_re: Vec<f64> = flat(&run(&format!("fft_re({x})")))
            .iter()
            .map(|v| v / n)
            .collect();
        let want_im: Vec<f64> = flat(&run(&format!("fft_im({x})")))
            .iter()
            .map(|v| v / n)
            .collect();
        let zeros = "list(0, 0, 0, 0, 0, 0, 0, 0)";
        close(
            &flat(&run(&format!("frft_re({x}, {zeros}, 1)"))),
            &want_re,
            1e-8,
        );
        close(
            &flat(&run(&format!("frft_im({x}, {zeros}, 1)"))),
            &want_im,
            1e-8,
        );
        // Unitarity: the norm is kept at a fractional order.
        let re = flat(&run(&format!("frft_re({x}, {zeros}, 0.37)")));
        let im = flat(&run(&format!("frft_im({x}, {zeros}, 0.37)")));
        let energy: f64 = re.iter().zip(&im).map(|(r, i)| r * r + i * i).sum();
        let input: f64 = flat(x).iter().map(|v| v * v).sum();
        assert!((energy - input).abs() < 1e-8, "{energy} vs {input}");
    }

    #[test]
    fn invalid_transforms_are_left_as_requests() {
        for src in [
            "wavedec(list(1, 2, 3), haar, 1)",
            "wavedec(list(1, 2, 3, 4), haar, 3)",
            "wavedec(list(1, 2, 3, 4), haar, 0)",
            "wavedec(list(1, 2, 3, 4), nonesuch, 1)",
            "waverec(list(list(1), list(1, 2, 3)), haar)",
            "waverec(list(1, 2, 3, 4), haar)",
            "dwt2(list(list(1, 2, 3), list(4, 5, 6)), haar)",
            "dwt2(list(1, 2, 3, 4), haar)",
            "frft_re(list(1), list(0), 1)",
            "frft_re(list(1, 2), list(0), 1)",
            "frft_re(list(1, a), list(0, 0), 1)",
        ] {
            let (text, reduced) = reduce_with(&[numeric()], src, &[]);
            assert!(!reduced, "{text}");
        }
    }
}
