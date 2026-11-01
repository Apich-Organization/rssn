//! # Discrete wavelet transforms and the discrete fractional Fourier transform
//!
//! * [`Wavelet`]: fifteen filter-bank families — the orthogonal `Haar` and Daubechies
//!   `db2`..`db4`, the least-asymmetric `sym4`..`sym6`, and the symmetric biorthogonal CDF
//!   splines `bior1.1`..`bior4.4`. A family carries two low-pass banks, the analysis one
//!   ([`Wavelet::analysis_low`]) and the synthesis one ([`Wavelet::synthesis_low`]); the
//!   high-pass banks are their alternating mirror. For the orthogonal families the two pairs
//!   coincide, for the biorthogonal ones they are distinct and every bank is symmetric or
//!   antisymmetric, hence linear phase.
//! * [`dwt`] / [`idwt`]: one level of the periodised transform — an even-length signal of `n`
//!   samples maps to `n/2` approximation and `n/2` detail coefficients, and back, exactly.
//! * [`dwt_with`] / [`idwt_with`]: one level with an explicit [`Boundary`] extension mode, for
//!   odd as well as even lengths.
//! * [`wavedec`] / [`waverec`] and [`wavedec_with`] / [`waverec_with`]: the multilevel
//!   decomposition and its inverse.
//! * [`dwt2`] / [`idwt2`]: the separable two-dimensional transform on a rectangular grid.
//! * [`frft_matrix`] / [`frft`]: the **discrete fractional Fourier transform** of Candan, Kutay
//!   and Ozaktas / Pei: the eigenvectors of the real symmetric matrix `S` that commutes with the
//!   unitary DFT are ordered like the Hermite functions and raised to the fractional power
//!   `exp(-i pi a k / 2)`; the result is unitary, `a = 1` is the unitary DFT, `a = 0` the
//!   identity, `a = 2` the reversal `x[n] -> x[-n mod N]`, and the transforms form a group in
//!   `a`.
//!
//! ## Conventions
//!
//! **Filters.** `h~` ([`Wavelet::analysis_low`]) is the analysis low-pass and `h`
//! ([`Wavelet::synthesis_low`]) the synthesis low-pass, both of length `L`. The high-pass banks
//! are the alternating mirror of the *dual* low-pass,
//! `g~[n] = (-1)^n h[L-1-n]` ([`Wavelet::analysis_high`]) and `g[n] = (-1)^n h~[L-1-n]`
//! ([`Wavelet::synthesis_high`]) — the CDF relation `g~[n] = (-1)^n h[1-n]` written for an array
//! indexed from zero. The tables are `PyWavelets`' `filter_bank` for each family: `h~` is `dec_lo`
//! reversed and `h` is `rec_lo`, so that the transforms below are correlations. Orthogonal
//! families satisfy `sum_n h[n] h[n+2k] = delta_k` and `h = h~`, biorthogonal ones satisfy the
//! biorthogonality `sum_n h~[n] h[n+2k] = delta_k`; either way `idwt(dwt(x)) = x`.
//!
//! **Periodic mode.** [`dwt`] and [`Boundary::Periodic`] use the circular form
//! `a[k] = sum_m h~[m] x[(2k + m) mod n]` on `ceil(n/2)` coefficients. For an odd `n` the signal
//! is first extended by one sample (`x[n] = x[0]`), transformed at the even length `n+1`, and the
//! inverse truncates the result back to `n`.
//!
//! **Extension modes.** `Zero`, `Symmetric` and `Reflect` evaluate the analysis on the extension
//! `x~` of `x` (see [`extend`]) at `a[j] = sum_m h~[m] x~[2j + m]`, and reconstruct with the
//! transposed filter bank, `x[i] += sum_j h[m] a[j] + g[m] d[j]` over the pairs with `2j + m = i`.
//! The synthesis banks decide which `j` can reach a sample of `x`, so the coefficient count is
//! `j_hi - j_lo + 1` with `j_lo = -floor(m_hi/2)` and `j_hi = floor((n - 1 - m_lo)/2)`, where
//! `[m_lo, m_hi]` is the support of the two synthesis banks. That support is not centred on zero
//! — `PyWavelets`' zero padding fixes the phase, and for `bior2.2` it is `[1, 5]` — so the two
//! sides are extended by different amounts, and the total is about `n/2 + (m_hi - m_lo)/2`
//! coefficients. The analysis low-pass bank is the longer of the two in the biorthogonal
//! families, so the analysis window reaches further outside the signal than the synthesis window;
//! only the synthesis support, however, decides how many coefficients are kept.
//!
//! **2-D layout.** [`dwt2`] transforms the rows and then the columns of a `rows x cols` array and
//! packs the four sub-bands into one `(2*ra) x (2*ca)` array, where `ra`/`ca` are the row/column
//! coefficient counts: `LL` in rows `0..ra`, columns `0..ca`; `LH` in rows `0..ra`, columns
//! `ca..2ca`; `HL` in rows `ra..2ra`, columns `0..ca`; `HH` in rows `ra..2ra`, columns `ca..2ca`.
//! The first letter names the band along the row index (the vertical direction), the second along
//! the column index.
//!
//! ## Limitations
//!
//! The families stop where the tables stop: no `db5` or higher, no `sym7` or higher, no coiflets,
//! discrete Meyer or other orthogonal families, and of the biorthogonal splines only the eight
//! listed (`bior5.5`, `bior6.8` are absent). There are no complex wavelets, no wavelet packets,
//! no lifting scheme and no non-separable or 3-D transforms. The 2-D transform is one level only:
//! the coefficient counts of the extension modes are not dyadic, so a multilevel 2-D scheme needs
//! shape bookkeeping that [`dwt2`] does not provide. Only the four [`Boundary`] modes exist, and
//! one mode applies to both ends. The published symlet tables are orthonormal only to about
//! `1e-12` (and `bior4.4` biorthogonal to about `1e-12`), so those families round-trip to
//! `1e-12 * max|x|` rather than to machine epsilon. Because the extension modes return about
//! `n/2 + (m_hi - m_lo)/2` coefficients, a multilevel decomposition stops shrinking once the
//! length is comparable to the filter length. The fractional Fourier matrix costs an `O(N^3)` symmetric
//! eigensolve (use it for lengths up to a few hundred; the continuous transform of a sampled
//! chirp-like function needs the sampling rate chosen by the caller).

use num_complex::Complex64;

use crate::kernels::dense::Mat;
use crate::kernels::dense::eigen_symmetric;

/// A wavelet family: orthogonal (`Haar`, `db2`..`db4`, `sym4`..`sym6`) or biorthogonal (`bior*`).
///
/// The biorthogonal variants are the CDF spline families of Cohen, Daubechies and Feauveau; the
/// digits are the `PyWavelets` name `biorNr.Nd`, with `Nr` vanishing moments in the analysis
/// wavelet and `Nd` in the synthesis wavelet ([`Wavelet::vanishing_moments`]).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Wavelet {
    /// Haar (`db1`): orthogonal, 1 vanishing moment, 2 taps.
    Haar,
    /// Daubechies with 2 vanishing moments (4 taps).
    Db2,
    /// Daubechies with 3 vanishing moments (6 taps).
    Db3,
    /// Daubechies with 4 vanishing moments (8 taps).
    Db4,
    /// Symlet with 4 vanishing moments (8 taps), least asymmetric.
    Sym4,
    /// Symlet with 5 vanishing moments (10 taps), least asymmetric.
    Sym5,
    /// Symlet with 6 vanishing moments (12 taps), least asymmetric.
    Sym6,
    /// Biorthogonal spline `bior1.1`: 1 analysis and 1 synthesis vanishing moment, 2 taps.
    Bior1_1,
    /// Biorthogonal spline `bior1.3`: 1 analysis and 3 synthesis vanishing moments, 6 taps.
    Bior1_3,
    /// Biorthogonal spline `bior1.5`: 1 analysis and 5 synthesis vanishing moments, 10 taps.
    Bior1_5,
    /// Biorthogonal spline `bior2.2`: 2 analysis and 2 synthesis vanishing moments, 6 taps.
    Bior2_2,
    /// Biorthogonal spline `bior2.4`: 2 analysis and 4 synthesis vanishing moments, 10 taps.
    Bior2_4,
    /// Biorthogonal spline `bior3.3`: 3 analysis and 3 synthesis vanishing moments, 8 taps.
    Bior3_3,
    /// Biorthogonal spline `bior3.5`: 3 analysis and 5 synthesis vanishing moments, 12 taps.
    Bior3_5,
    /// Biorthogonal spline `bior4.4` (the `9/7` pair): 4 analysis and 4 synthesis moments, 10 taps.
    Bior4_4,
}

impl Wavelet {
    /// Every family, orthogonal first.
    pub const ALL: [Self; 15] = [
        Self::Haar,
        Self::Db2,
        Self::Db3,
        Self::Db4,
        Self::Sym4,
        Self::Sym5,
        Self::Sym6,
        Self::Bior1_1,
        Self::Bior1_3,
        Self::Bior1_5,
        Self::Bior2_2,
        Self::Bior2_4,
        Self::Bior3_3,
        Self::Bior3_5,
        Self::Bior4_4,
    ];

    /// `true` when one bank serves both analysis and synthesis (orthonormal filters).
    #[must_use]
    pub const fn is_orthogonal(self) -> bool {
        !matches!(
            self,
            Self::Bior1_1
                | Self::Bior1_3
                | Self::Bior1_5
                | Self::Bior2_2
                | Self::Bior2_4
                | Self::Bior3_3
                | Self::Bior3_5
                | Self::Bior4_4
        )
    }

    /// The analysis (decomposition) low-pass filter `h~`, normalised to `sum h~ = sqrt 2`.
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
            | Self::Sym4 => SYM4.to_vec(),
            | Self::Sym5 => SYM5.to_vec(),
            | Self::Sym6 => SYM6.to_vec(),
            | Self::Bior1_1 => BIOR11_LOW.to_vec(),
            | Self::Bior1_3 => BIOR13_LOW.to_vec(),
            | Self::Bior1_5 => BIOR15_LOW.to_vec(),
            | Self::Bior2_2 => BIOR22_LOW.to_vec(),
            | Self::Bior2_4 => BIOR24_LOW.to_vec(),
            | Self::Bior3_3 => BIOR33_LOW.to_vec(),
            | Self::Bior3_5 => BIOR35_LOW.to_vec(),
            | Self::Bior4_4 => BIOR44_LOW.to_vec(),
        }
    }

    /// The analysis (decomposition) high-pass filter `g~[n] = (-1)^n h[L-1-n]`.
    #[must_use]
    pub fn high_pass(self) -> Vec<f64> {
        self.analysis_high()
    }

    /// The analysis high-pass filter `g~`, the alternating mirror of the synthesis low-pass.
    #[must_use]
    pub fn analysis_high(self) -> Vec<f64> {
        alternating_mirror(&self.synthesis_low())
    }

    /// The analysis low-pass filter `h~` (the same as [`Wavelet::low_pass`]).
    #[must_use]
    pub fn analysis_low(self) -> Vec<f64> {
        self.low_pass()
    }

    /// The synthesis (reconstruction) low-pass filter `h`; equal to [`Wavelet::low_pass`] for the
    /// orthogonal families.
    #[must_use]
    pub fn synthesis_low(self) -> Vec<f64> {
        match self {
            | Self::Haar
            | Self::Db2
            | Self::Db3
            | Self::Db4
            | Self::Sym4
            | Self::Sym5
            | Self::Sym6 => self.low_pass(),
            | Self::Bior1_1 => BIOR11_SYN.to_vec(),
            | Self::Bior1_3 => BIOR13_SYN.to_vec(),
            | Self::Bior1_5 => BIOR15_SYN.to_vec(),
            | Self::Bior2_2 => BIOR22_SYN.to_vec(),
            | Self::Bior2_4 => BIOR24_SYN.to_vec(),
            | Self::Bior3_3 => BIOR33_SYN.to_vec(),
            | Self::Bior3_5 => BIOR35_SYN.to_vec(),
            | Self::Bior4_4 => BIOR44_SYN.to_vec(),
        }
    }

    /// The synthesis high-pass filter `g[n] = (-1)^n h~[L-1-n]`.
    #[must_use]
    pub fn synthesis_high(self) -> Vec<f64> {
        alternating_mirror(&self.low_pass())
    }

    /// Number of filter taps (the length of the banks).
    #[must_use]
    pub fn taps(self) -> usize {
        self.low_pass().len()
    }

    /// The vanishing moments of the analysis and of the synthesis wavelet, `(Nr, Nd)`.
    ///
    /// Counted from the filters: the analysis wavelet annihilates polynomials of degree below
    /// `Nr`, the synthesis wavelet those below `Nd`. For the orthogonal families the two agree
    /// and equal the Daubechies order.
    #[must_use]
    pub const fn vanishing_moments(self) -> (usize, usize) {
        match self {
            | Self::Haar | Self::Bior1_1 => (1, 1),
            | Self::Db2 | Self::Bior2_2 => (2, 2),
            | Self::Db3 | Self::Bior3_3 => (3, 3),
            | Self::Db4 | Self::Sym4 | Self::Bior4_4 => (4, 4),
            | Self::Sym5 => (5, 5),
            | Self::Sym6 => (6, 6),
            | Self::Bior1_3 => (1, 3),
            | Self::Bior1_5 => (1, 5),
            | Self::Bior2_4 => (2, 4),
            | Self::Bior3_5 => (3, 5),
        }
    }

    /// The number of coefficients [`dwt_with`] returns for a signal of length `n`.
    ///
    /// `ceil(n/2)` in [`Boundary::Periodic`] mode, and `j_hi - j_lo + 1` — about
    /// `n/2 + (m_hi - m_lo)/2` — in the extension modes, where `[m_lo, m_hi]` is the support of
    /// the synthesis banks (see the module documentation).
    #[must_use]
    pub fn coeff_len(
        self,
        n: usize,
        b: Boundary,
    ) -> usize {
        match b {
            | Boundary::Periodic => n.div_ceil(2),
            | _ => {
                let (_, _, h, g) = self.banks();
                let (_, count) = ext_layout(&h, &g, n);
                count
            },
        }
    }

    /// The four banks `(h~, g~, h, g)` in one call.
    fn banks(self) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
        let low = self.low_pass();
        let synthesis = self.synthesis_low();
        let analysis_high = alternating_mirror(&synthesis);
        let synthesis_high = alternating_mirror(&low);
        (low, analysis_high, synthesis, synthesis_high)
    }
}

/// `g[k] = (-1)^k h[L-1-k]`, the high-pass bank that pairs with the low-pass `h`.
fn alternating_mirror(h: &[f64]) -> Vec<f64> {
    let l = h.len();
    (0..l)
        .map(|k| {
            let v = h.get(l - 1 - k).copied().unwrap_or(0.0);
            if k % 2 == 0 { v } else { -v }
        })
        .collect()
}

// `PyWavelets` `Wavelet(family).filter_bank`: `*_LOW` is `dec_lo` reversed (the analysis low-pass),
// `*_HIGH` is `rec_lo` (the synthesis low-pass). Every bank keeps `PyWavelets`' zero padding, which
// fixes the relative alignment of the two filters.
const SYM4: [f64; 8] = [
    0.032_223_100_604_042_7,
    -0.012_603_967_262_037_833,
    -0.099_219_543_576_847_22,
    0.297_857_795_605_277_36,
    0.803_738_751_805_916_1,
    0.497_618_667_632_015_45,
    -0.029_635_527_645_998_51,
    -0.075_765_714_789_273_33,
];
const SYM5: [f64; 10] = [
    0.019_538_882_735_286_728,
    -0.021_101_834_024_758_855,
    -0.175_328_089_908_450_47,
    0.016_602_105_764_522_32,
    0.633_978_963_458_211_9,
    0.723_407_690_402_420_6,
    0.199_397_533_977_393_6,
    -0.039_134_249_302_383_094,
    0.029_519_490_925_774_643,
    0.027_333_068_345_077_982,
];
const SYM6: [f64; 12] = [
    -0.007_800_708_325_034_148,
    0.001_767_711_864_242_803_6,
    0.044_724_901_770_665_78,
    -0.021_060_292_512_300_564,
    -0.072_637_522_786_462_52,
    0.337_929_421_727_621_8,
    0.787_641_141_030_194,
    0.491_055_941_926_746_6,
    -0.048_311_742_585_633,
    -0.117_990_111_148_190_57,
    0.003_490_712_084_217_470_2,
    0.015_404_109_327_027_373,
];
const BIOR11_LOW: [f64; 2] = [std::f64::consts::FRAC_1_SQRT_2, std::f64::consts::FRAC_1_SQRT_2];
const BIOR11_SYN: [f64; 2] = [std::f64::consts::FRAC_1_SQRT_2, std::f64::consts::FRAC_1_SQRT_2];
const BIOR13_LOW: [f64; 6] = [
    -0.088_388_347_648_318_45,
    0.088_388_347_648_318_45,
    std::f64::consts::FRAC_1_SQRT_2,
    std::f64::consts::FRAC_1_SQRT_2,
    0.088_388_347_648_318_45,
    -0.088_388_347_648_318_45,
];
const BIOR13_SYN: [f64; 6] = [
    0.0,
    0.0,
    std::f64::consts::FRAC_1_SQRT_2,
    std::f64::consts::FRAC_1_SQRT_2,
    0.0,
    0.0,
];
const BIOR15_LOW: [f64; 10] = [
    0.016_572_815_184_059_706,
    -0.016_572_815_184_059_706,
    -0.121_533_978_016_437_85,
    0.121_533_978_016_437_85,
    std::f64::consts::FRAC_1_SQRT_2,
    std::f64::consts::FRAC_1_SQRT_2,
    0.121_533_978_016_437_85,
    -0.121_533_978_016_437_85,
    -0.016_572_815_184_059_706,
    0.016_572_815_184_059_706,
];
const BIOR15_SYN: [f64; 10] = [
    0.0,
    0.0,
    0.0,
    0.0,
    std::f64::consts::FRAC_1_SQRT_2,
    std::f64::consts::FRAC_1_SQRT_2,
    0.0,
    0.0,
    0.0,
    0.0,
];
const BIOR22_LOW: [f64; 6] = [
    -0.176_776_695_296_636_9,
    0.353_553_390_593_273_8,
    1.060_660_171_779_821_2,
    0.353_553_390_593_273_8,
    -0.176_776_695_296_636_9,
    0.0,
];
const BIOR22_SYN: [f64; 6] = [
    0.0,
    0.353_553_390_593_273_8,
    std::f64::consts::FRAC_1_SQRT_2,
    0.353_553_390_593_273_8,
    0.0,
    0.0,
];
const BIOR24_LOW: [f64; 10] = [
    0.033_145_630_368_119_41,
    -0.066_291_260_736_238_82,
    -0.176_776_695_296_636_9,
    0.419_844_651_329_512_6,
    0.994_368_911_043_582_5,
    0.419_844_651_329_512_6,
    -0.176_776_695_296_636_9,
    -0.066_291_260_736_238_82,
    0.033_145_630_368_119_41,
    0.0,
];
const BIOR24_SYN: [f64; 10] = [
    0.0,
    0.0,
    0.0,
    0.353_553_390_593_273_8,
    std::f64::consts::FRAC_1_SQRT_2,
    0.353_553_390_593_273_8,
    0.0,
    0.0,
    0.0,
    0.0,
];
const BIOR33_LOW: [f64; 8] = [
    0.066_291_260_736_238_82,
    -0.198_873_782_208_716_5,
    -0.154_679_608_384_557_27,
    0.994_368_911_043_582_5,
    0.994_368_911_043_582_5,
    -0.154_679_608_384_557_27,
    -0.198_873_782_208_716_5,
    0.066_291_260_736_238_82,
];
const BIOR33_SYN: [f64; 8] = [
    0.0,
    0.0,
    0.176_776_695_296_636_9,
    0.530_330_085_889_910_6,
    0.530_330_085_889_910_6,
    0.176_776_695_296_636_9,
    0.0,
    0.0,
];
const BIOR35_LOW: [f64; 12] = [
    -0.013_810_679_320_049_757,
    0.041_432_037_960_149_27,
    0.052_480_581_416_189_075,
    -0.267_927_178_808_965_27,
    -0.071_815_532_464_258_73,
    0.966_747_552_403_483,
    0.966_747_552_403_483,
    -0.071_815_532_464_258_73,
    -0.267_927_178_808_965_27,
    0.052_480_581_416_189_075,
    0.041_432_037_960_149_27,
    -0.013_810_679_320_049_757,
];
const BIOR35_SYN: [f64; 12] = [
    0.0,
    0.0,
    0.0,
    0.0,
    0.176_776_695_296_636_9,
    0.530_330_085_889_910_6,
    0.530_330_085_889_910_6,
    0.176_776_695_296_636_9,
    0.0,
    0.0,
    0.0,
    0.0,
];
const BIOR44_LOW: [f64; 10] = [
    0.037_828_455_507_264_04,
    -0.023_849_465_019_556_843,
    -0.110_624_404_418_437_18,
    0.377_402_855_612_830_66,
    0.852_698_679_008_893_8,
    0.377_402_855_612_830_66,
    -0.110_624_404_418_437_18,
    -0.023_849_465_019_556_843,
    0.037_828_455_507_264_04,
    0.0,
];
const BIOR44_SYN: [f64; 10] = [
    0.0,
    -0.064_538_882_628_697_06,
    -0.040_689_417_609_164_06,
    0.418_092_273_221_617_24,
    0.788_485_616_405_582_9,
    0.418_092_273_221_617_24,
    -0.040_689_417_609_164_06,
    -0.064_538_882_628_697_06,
    0.0,
    0.0,
];

/// How a signal is extended past its ends when a filter reaches outside it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Boundary {
    /// Periodic: `x[i] = x[i mod n]`. This is the circular (periodised) transform of [`dwt`].
    Periodic,
    /// Zero padding: samples outside `[0, n)` are zero.
    Zero,
    /// Half-sample symmetric: `.. x[1], x[0] | x[0], x[1], .. x[n-1] | x[n-1], x[n-2] ..`, i.e.
    /// the period of the extension is `2n` and the edge samples are repeated once.
    Symmetric,
    /// Whole-sample symmetric: `.. x[2], x[1] | x[0], x[1], .. x[n-1] | x[n-2], x[n-3] ..`, i.e.
    /// the period of the extension is `2(n-1)` and the edge samples are not repeated.
    Reflect,
}

/// Sample `i` of the boundary extension of `x` in mode `b`.
///
/// `Periodic` has period `n`, `Symmetric` period `2n` (mirror with the edge sample repeated) and
/// `Reflect` period `2(n-1)` (mirror without repeating it); `Zero` returns `0.0` outside
/// `[0, n)`. For `n = 1` every mode returns `x[0]`, and an empty signal has no samples.
#[must_use]
pub fn extend(
    x: &[f64],
    i: isize,
    b: Boundary,
) -> f64 {
    let n = x.len();
    let Ok(len) = isize::try_from(n) else {
        return 0.0;
    };
    if len == 0 {
        return 0.0;
    }
    let idx = match b {
        | Boundary::Zero => {
            if i < 0 || i >= len {
                return 0.0;
            }
            i
        },
        | Boundary::Periodic => i.rem_euclid(len),
        | Boundary::Symmetric => {
            let r = i.rem_euclid(2 * len);
            if r < len { r } else { 2 * len - 1 - r }
        },
        | Boundary::Reflect => {
            if len == 1 {
                0
            } else {
                let p = 2 * (len - 1);
                let r = i.rem_euclid(p);
                if r < len { r } else { p - r }
            }
        },
    };
    usize::try_from(idx)
        .ok()
        .and_then(|k| x.get(k))
        .copied()
        .unwrap_or(0.0)
}

/// The first coefficient index `j_lo` and the coefficient count of an extension-mode transform.
///
/// `j_lo = -floor(m_hi/2)` and `j_hi = floor((n-1-m_lo)/2)`, with `[m_lo, m_hi]` the support of
/// the two synthesis banks: those are exactly the `j` whose synthesis window `2j + m` reaches a
/// sample of `x`.
fn ext_layout(
    low: &[f64],
    high: &[f64],
    n: usize,
) -> (isize, usize) {
    let mut lo = usize::MAX;
    let mut hi = 0;
    for (m, v) in low.iter().enumerate().chain(high.iter().enumerate()) {
        if *v != 0.0 {
            lo = lo.min(m);
            hi = m;
        }
    }
    if lo == usize::MAX {
        return (0, 0);
    }
    let Ok(len) = isize::try_from(n) else {
        return (0, 0);
    };
    // j_lo = ceil(-m_hi/2) = -floor(m_hi/2); j_hi = floor((n - 1 - m_lo)/2)
    let j_lo = -(isize::try_from(hi).unwrap_or(0) / 2);
    let j_hi = (len - 1 - isize::try_from(lo).unwrap_or(0)).div_euclid(2);
    if j_hi < j_lo {
        return (j_lo, 0);
    }
    (j_lo, usize::try_from(j_hi - j_lo + 1).unwrap_or(0))
}

/// Errors of the wavelet transforms.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum WaveletError {
    /// The signal length is odd or zero, where an even length is required.
    OddLength,
    /// More levels were asked for than the length allows.
    TooManyLevels,
    /// The coefficient lists do not fit together, or a shape does not match the data.
    Shape,
}

/// One level of the periodised discrete wavelet transform: returns `(approximation, detail)`.
///
/// `a[k] = sum_m h~[m] x[(2k + m) mod n]`, `d[k] = sum_m g~[m] x[(2k + m) mod n]` with the
/// analysis banks. An even signal of length `n` maps to `n/2` coefficients of each kind. For odd
/// lengths or another boundary rule use [`dwt_with`].
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
    dwt_with(x, w, Boundary::Periodic)
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
    idwt_with(a, d, w, Boundary::Periodic, 2 * a.len())
}

/// One level of the discrete wavelet transform with an explicit [`Boundary`] mode.
///
/// Returns `(approximation, detail)`, each of the length that [`Wavelet::coeff_len`] reports for
/// `x.len()` and `b`. In [`Boundary::Periodic`] mode this is the circular transform of [`dwt`],
/// extended by one sample for an odd length; in the other modes the analysis is evaluated on the
/// boundary extension of the signal (see the module documentation for the exact index rule).
///
/// # Errors
/// [`WaveletError::OddLength`] for an empty signal.
pub fn dwt_with(
    x: &[f64],
    w: Wavelet,
    b: Boundary,
) -> Result<(Vec<f64>, Vec<f64>), WaveletError> {
    let n = x.len();
    if n == 0 {
        return Err(WaveletError::OddLength);
    }
    let (ha, ga, hs, gs) = w.banks();
    match b {
        | Boundary::Periodic => {
            // an odd signal is extended by one sample so that the circular transform acts on an
            // even length; the inverse truncates the extra sample away again
            let mut y = x.to_vec();
            if !n.is_multiple_of(2) {
                y.extend(x.first().copied());
            }
            let len = y.len();
            let half = len / 2;
            let a = (0..half)
                .map(|k| {
                    ha.iter()
                        .enumerate()
                        .map(|(m, hm)| hm * y.get((2 * k + m) % len).copied().unwrap_or(0.0))
                        .sum()
                })
                .collect();
            let d = (0..half)
                .map(|k| {
                    ga.iter()
                        .enumerate()
                        .map(|(m, gm)| gm * y.get((2 * k + m) % len).copied().unwrap_or(0.0))
                        .sum()
                })
                .collect();
            Ok((a, d))
        },
        | _ => {
            let (j_lo, count) = ext_layout(&hs, &gs, n);
            let a = (0..count)
                .map(|q| {
                    let j = j_lo + isize::try_from(q).unwrap_or(0);
                    ha.iter()
                        .enumerate()
                        .map(|(m, hm)| hm * extend(x, 2 * j + isize::try_from(m).unwrap_or(0), b))
                        .sum()
                })
                .collect();
            let d = (0..count)
                .map(|q| {
                    let j = j_lo + isize::try_from(q).unwrap_or(0);
                    ga.iter()
                        .enumerate()
                        .map(|(m, gm)| gm * extend(x, 2 * j + isize::try_from(m).unwrap_or(0), b))
                        .sum()
                })
                .collect();
            Ok((a, d))
        },
    }
}

/// Inverse of [`dwt_with`]; `n` is the length of the signal that produced the coefficients.
///
/// # Errors
/// [`WaveletError::Shape`] when the lists differ in length, are empty, or do not match
/// `w.coeff_len(n, b)`.
pub fn idwt_with(
    a: &[f64],
    d: &[f64],
    w: Wavelet,
    b: Boundary,
    n: usize,
) -> Result<Vec<f64>, WaveletError> {
    if a.len() != d.len() || a.is_empty() || a.len() != w.coeff_len(n, b) {
        return Err(WaveletError::Shape);
    }
    let (_, _, hs, gs) = w.banks();
    match b {
        | Boundary::Periodic => {
            // the odd case reconstructs the padded even signal and drops the extra sample
            let len = if n.is_multiple_of(2) { n } else { n + 1 };
            let mut x = vec![0.0; len];
            for (k, (ak, dk)) in a.iter().zip(d).enumerate() {
                for (m, (hm, gm)) in hs.iter().zip(&gs).enumerate() {
                    if let Some(slot) = x.get_mut((2 * k + m) % len) {
                        *slot += hm * ak + gm * dk;
                    }
                }
            }
            x.truncate(n);
            Ok(x)
        },
        | _ => {
            let (j_lo, _) = ext_layout(&hs, &gs, n);
            let mut x = vec![0.0; n];
            for (q, (aq, dq)) in a.iter().zip(d).enumerate() {
                let j = j_lo + isize::try_from(q).unwrap_or(0);
                for (m, (hm, gm)) in hs.iter().zip(&gs).enumerate() {
                    let i = 2 * j + isize::try_from(m).unwrap_or(0);
                    let slot = usize::try_from(i).ok().and_then(|k| x.get_mut(k));
                    if let Some(slot) = slot {
                        *slot += hm * aq + gm * dq;
                    }
                }
            }
            Ok(x)
        },
    }
}

/// Multilevel decomposition: `[approximation_L, detail_L, ..., detail_1]`.
///
/// Requires `2^levels` to divide the length, as before; use [`wavedec_with`] for the general case.
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

/// Multilevel decomposition with a [`Boundary`] mode, for any non-empty signal.
///
/// Every level keeps the number of coefficients that [`Wavelet::coeff_len`] reports for its
/// length and `b`, so the levels stop
/// shrinking once the length is comparable to the filter length (see the module documentation).
///
/// # Errors
/// [`WaveletError::TooManyLevels`] when `levels` is zero or `2^(levels-1) > x.len()`.
pub fn wavedec_with(
    x: &[f64],
    w: Wavelet,
    levels: usize,
    b: Boundary,
) -> Result<Vec<Vec<f64>>, WaveletError> {
    let n = x.len();
    if n == 0 || levels == 0 || levels > 60 || (n >> (levels - 1)) == 0 {
        return Err(WaveletError::TooManyLevels);
    }
    let mut details = Vec::with_capacity(levels);
    let mut current = x.to_vec();
    for _ in 0..levels {
        if current.is_empty() {
            return Err(WaveletError::TooManyLevels);
        }
        let (a, d) = dwt_with(&current, w, b)?;
        details.push(d);
        current = a;
    }
    let mut out = vec![current];
    out.extend(details.into_iter().rev());
    Ok(out)
}

/// Inverse of [`wavedec_with`]; `n` is the length of the decomposed signal.
///
/// # Errors
/// [`WaveletError::Shape`] when the list or one of its lengths does not fit `n` and the mode.
pub fn waverec_with(
    coeffs: &[Vec<f64>],
    w: Wavelet,
    b: Boundary,
    n: usize,
) -> Result<Vec<f64>, WaveletError> {
    let (first, rest) = coeffs.split_first().ok_or(WaveletError::Shape)?;
    let levels = rest.len();
    let mut lens = Vec::with_capacity(levels + 1);
    let mut current_len = n;
    lens.push(current_len);
    for _ in 0..levels {
        current_len = w.coeff_len(current_len, b);
        lens.push(current_len);
    }
    if first.len() != lens.get(levels).copied().unwrap_or(0) {
        return Err(WaveletError::Shape);
    }
    for (k, d) in rest.iter().enumerate() {
        if d.len() != lens.get(levels - k).copied().unwrap_or(0) {
            return Err(WaveletError::Shape);
        }
    }
    let mut current = first.clone();
    for (k, d) in rest.iter().enumerate() {
        let out_len = lens.get(levels - k - 1).copied().unwrap_or(0);
        current = idwt_with(&current, d, w, b, out_len)?;
    }
    Ok(current)
}

/// One level of the separable two-dimensional transform, `x` row-major with `rows x cols` entries.
///
/// The rows are transformed first, then the columns, and the four sub-bands are packed into a
/// `(2*ra) x (2*ca)` row-major array with `ra = w.coeff_len(rows, b)` and
/// `ca = w.coeff_len(cols, b)`: `LL` at rows `0..ra`, columns `0..ca`; `LH` at rows `0..ra`,
/// columns `ca..2ca`; `HL` at rows `ra..2ra`, columns `0..ca`; `HH` at rows `ra..2ra`, columns
/// `ca..2ca`. `LH` is the low band along the row index and the high band along the column index.
///
/// # Errors
/// [`WaveletError::Shape`] when `rows * cols` does not match `x.len()`, or when a row or column
/// transform fails.
pub fn dwt2(
    x: &[f64],
    rows: usize,
    cols: usize,
    w: Wavelet,
    b: Boundary,
) -> Result<Vec<f64>, WaveletError> {
    let total = rows.checked_mul(cols).ok_or(WaveletError::Shape)?;
    if rows == 0 || cols == 0 || x.len() != total {
        return Err(WaveletError::Shape);
    }
    let ca = w.coeff_len(cols, b);
    let ra = w.coeff_len(rows, b);
    let band = rows.checked_mul(ca).ok_or(WaveletError::Shape)?;
    let mut low = vec![0.0; band];
    let mut high = vec![0.0; band];
    for r in 0..rows {
        let start = r * cols;
        let row: Vec<f64> = x
            .get(start..start + cols)
            .map_or_else(Vec::new, <[f64]>::to_vec);
        let (a, d) = dwt_with(&row, w, b)?;
        for (c, v) in a.iter().enumerate() {
            if let Some(slot) = low.get_mut(r * ca + c) {
                *slot = *v;
            }
        }
        for (c, v) in d.iter().enumerate() {
            if let Some(slot) = high.get_mut(r * ca + c) {
                *slot = *v;
            }
        }
    }
    let out_rows = 2 * ra;
    let out_cols = 2 * ca;
    let mut out = vec![0.0; out_rows.checked_mul(out_cols).ok_or(WaveletError::Shape)?];
    for c in 0..ca {
        for (src, offset) in [(&low, 0_usize), (&high, ra)] {
            let column: Vec<f64> = (0..rows)
                .map(|r| src.get(r * ca + c).copied().unwrap_or(0.0))
                .collect();
            let (a, d) = dwt_with(&column, w, b)?;
            for (i, v) in a.iter().enumerate() {
                if let Some(slot) = out.get_mut((offset + i) * out_cols + c) {
                    *slot = *v;
                }
            }
            for (i, v) in d.iter().enumerate() {
                if let Some(slot) = out.get_mut((offset + i) * out_cols + ca + c) {
                    *slot = *v;
                }
            }
        }
    }
    Ok(out)
}

/// Inverse of [`dwt2`]: `rows` and `cols` are the shape of the original array.
///
/// # Errors
/// [`WaveletError::Shape`] when the packed array does not have the expected
/// `(2*ra) x (2*ca)` size, or when a column or row transform fails.
pub fn idwt2(
    c: &[f64],
    rows: usize,
    cols: usize,
    w: Wavelet,
    b: Boundary,
) -> Result<Vec<f64>, WaveletError> {
    if rows == 0 || cols == 0 {
        return Err(WaveletError::Shape);
    }
    let ca = w.coeff_len(cols, b);
    let ra = w.coeff_len(rows, b);
    let out_rows = 2 * ra;
    let out_cols = 2 * ca;
    if c.len() != out_rows.checked_mul(out_cols).ok_or(WaveletError::Shape)? {
        return Err(WaveletError::Shape);
    }
    let band = rows.checked_mul(ca).ok_or(WaveletError::Shape)?;
    let mut low = vec![0.0; band];
    let mut high = vec![0.0; band];
    for col in 0..ca {
        for (dst, offset) in [(&mut low, 0_usize), (&mut high, ra)] {
            let ll: Vec<f64> = (0..ra)
                .map(|i| c.get((offset + i) * out_cols + col).copied().unwrap_or(0.0))
                .collect();
            let hh: Vec<f64> = (0..ra)
                .map(|i| {
                    c.get((offset + i) * out_cols + ca + col)
                        .copied()
                        .unwrap_or(0.0)
                })
                .collect();
            let column = idwt_with(&ll, &hh, w, b, rows)?;
            for (r, v) in column.iter().enumerate() {
                if let Some(slot) = dst.get_mut(r * ca + col) {
                    *slot = *v;
                }
            }
        }
    }
    let mut out = vec![0.0; rows.checked_mul(cols).ok_or(WaveletError::Shape)?];
    for r in 0..rows {
        let row_low: Vec<f64> = (0..ca)
            .map(|k| low.get(r * ca + k).copied().unwrap_or(0.0))
            .collect();
        let row_high: Vec<f64> = (0..ca)
            .map(|k| high.get(r * ca + k).copied().unwrap_or(0.0))
            .collect();
        let row = idwt_with(&row_low, &row_high, w, b, cols)?;
        for (k, v) in row.iter().enumerate() {
            if let Some(slot) = out.get_mut(r * cols + k) {
                *slot = *v;
            }
        }
    }
    Ok(out)
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

#[cfg(test)]
mod tests {
    //! Unit tests for the parts the integration tests cannot reach: the private bank assembly,
    //! the coefficient-count arithmetic and the self-duality of the orthogonal families.

    #![allow(clippy::indexing_slicing)]

    use super::Boundary;
    use super::Wavelet;
    use super::alternating_mirror;
    use super::dwt;
    use super::dwt_with;
    use super::ext_layout;
    use super::idwt;
    use super::idwt_with;

    const MODES: [Boundary; 4] = [
        Boundary::Periodic,
        Boundary::Zero,
        Boundary::Symmetric,
        Boundary::Reflect,
    ];

    #[test]
    fn the_alternating_mirror_reverses_and_flips_every_second_tap() {
        assert_eq!(alternating_mirror(&[1.0, 2.0, 3.0, 4.0]), vec![4.0, -3.0, 2.0, -1.0]);
        assert_eq!(alternating_mirror(&[2.0]), vec![2.0]);
    }

    #[test]
    fn the_coefficient_count_is_what_the_transform_returns() {
        for w in Wavelet::ALL {
            for b in MODES {
                for n in [1_usize, 2, 3, 7, 8, 9, 16, 17, 100] {
                    let x: Vec<f64> = (0..n).map(|i| (i as f64).sin()).collect();
                    let (a, d) = dwt_with(&x, w, b).unwrap();
                    assert_eq!(a.len(), w.coeff_len(n, b), "{w:?} {b:?} {n}");
                    assert_eq!(d.len(), a.len(), "{w:?} {b:?} {n}");
                    assert!(!a.is_empty(), "{w:?} {b:?} {n}");
                }
            }
        }
    }

    #[test]
    fn the_extension_layout_covers_the_synthesis_support() {
        // bior1.3: the synthesis banks span [0, 5], so the first coefficient of a length-7 signal
        // sits at j = -2 and six coefficients reach every sample of the signal
        let high = alternating_mirror(&Wavelet::Bior1_3.analysis_low());
        let (j_lo, count) = ext_layout(&Wavelet::Bior1_3.synthesis_low(), &high, 7);
        assert_eq!((j_lo, count), (-2, 6));
        // Haar: support [0, 1], so the layout is the plain ceil(n/2) of the periodised transform
        let haar = Wavelet::Haar;
        let (j_lo, count) = ext_layout(
            &haar.synthesis_low(),
            &alternating_mirror(&haar.analysis_low()),
            9,
        );
        assert_eq!((j_lo, count), (0, 5));
    }

    #[test]
    fn the_orthogonal_families_are_self_dual() {
        for w in Wavelet::ALL.iter().filter(|w| w.is_orthogonal()) {
            assert_eq!(w.analysis_low(), w.synthesis_low(), "{w:?}");
            assert_eq!(w.analysis_high(), w.synthesis_high(), "{w:?}");
            assert_eq!(w.high_pass(), w.analysis_high(), "{w:?}");
            assert_eq!(w.low_pass(), w.analysis_low(), "{w:?}");
            assert_eq!(w.taps(), w.analysis_low().len(), "{w:?}");
        }
        // the biorthogonal pairs are not self-dual, except bior1.1 which is the Haar pair
        for w in Wavelet::ALL.iter().filter(|w| !w.is_orthogonal()) {
            if *w == Wavelet::Bior1_1 {
                assert_eq!(w.analysis_low(), w.synthesis_low());
            } else {
                assert_ne!(w.analysis_low(), w.synthesis_low(), "{w:?}");
            }
        }
    }

    #[test]
    fn the_periodised_transform_of_an_odd_signal_is_the_padded_even_one() {
        // a signal of length 7 is transformed at length 8 with x[7] = x[0]; the coefficients
        // therefore agree with those of the padded even signal
        let x: Vec<f64> = (0..7).map(|i| f64::from(i) * 0.25 - 0.5).collect();
        let mut padded = x.clone();
        padded.push(x[0]);
        for w in Wavelet::ALL {
            let (a, d) = dwt_with(&x, w, Boundary::Periodic).unwrap();
            let (a2, d2) = dwt(&padded, w).unwrap();
            assert_eq!(a, a2, "{w:?}");
            assert_eq!(d, d2, "{w:?}");
            let back = idwt_with(&a, &d, w, Boundary::Periodic, 7).unwrap();
            assert_eq!(back.len(), 7);
            assert!(back.iter().zip(&x).all(|(p, q)| (p - q).abs() < 1e-9));
            assert_eq!(
                idwt(&a, &d, w).map(|v| v.len()),
                Ok(8),
                "{w:?}: the even inverse keeps the padded sample"
            );
        }
    }
}
