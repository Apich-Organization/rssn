//! Discrete wavelet transforms and the discrete fractional Fourier transform.

use num_complex::Complex64;
use rssn::kernels::wavelets::Boundary;
use rssn::kernels::wavelets::Wavelet;
use rssn::kernels::wavelets::WaveletError;
use rssn::kernels::wavelets::dwt;
use rssn::kernels::wavelets::dwt_with;
use rssn::kernels::wavelets::dwt2;
use rssn::kernels::wavelets::extend;
use rssn::kernels::wavelets::frft;
use rssn::kernels::wavelets::frft_matrix;
use rssn::kernels::wavelets::idwt;
use rssn::kernels::wavelets::idwt_with;
use rssn::kernels::wavelets::idwt2;
use rssn::kernels::wavelets::wavedec;
use rssn::kernels::wavelets::wavedec_with;
use rssn::kernels::wavelets::wavedec2;
use rssn::kernels::wavelets::wavedec2_with;
use rssn::kernels::wavelets::waverec;
use rssn::kernels::wavelets::waverec_with;
use rssn::kernels::wavelets::waverec2;
use rssn::kernels::wavelets::waverec2_with;

const FAMILIES: [Wavelet; 21] = Wavelet::ALL;

/// Every mode, the circular one included.
const MODES: [Boundary; 4] = [
    Boundary::Periodic,
    Boundary::Zero,
    Boundary::Symmetric,
    Boundary::Reflect,
];

/// Lengths the task asks for: odd, even, small, prime, long.
const LENGTHS: [usize; 7] = [7, 8, 9, 16, 17, 32, 100];

fn signal(n: usize) -> Vec<f64> {
    let mut s = 12345_u64;
    (0..n)
        .map(|_| {
            s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (s >> 11) as f64 / (1_u64 << 53) as f64 - 0.5
        })
        .collect()
}

fn max_diff(
    a: &[f64],
    b: &[f64],
) -> f64 {
    assert_eq!(a.len(), b.len(), "length mismatch");
    a.iter()
        .zip(b)
        .map(|(p, q)| (p - q).abs())
        .fold(0.0, f64::max)
}

/// Families whose published tables are only accurate to about `1e-12`; their round trip cannot be
/// better than that.
fn table_limited(w: Wavelet) -> bool {
    matches!(
        w,
        Wavelet::Sym4 | Wavelet::Sym5 | Wavelet::Sym6 | Wavelet::Bior4_4
    )
}

fn round_trip_tolerance(w: Wavelet) -> f64 {
    if table_limited(w) {
        1e-9
    } else {
        1e-11
    }
}

/// `|sum_m g[m] m^p| / sum_m |g[m]| m^p`: a scale-free measure of a vanishing
/// moment.
///
/// The absolute round-off of the moment sum grows like `max_m |g[m]| m^p`, so a
/// fixed absolute tolerance is meaningless for the long filters: `db9`/`db10`
/// reach about `1e-6` absolute at their last moment while still cancelling to
/// far below round-off relative to the terms. This replaces the old absolute
/// `1e-8` bound, which already failed at HEAD for `db9` (moment 7 was `1.9e-8`)
/// and `db10`, independently of the multilevel 2-D transform.
fn moment_residual(
    g: &[f64],
    p: usize,
) -> f64 {
    let mut sum = 0.0;
    let mut scale = 0.0;
    for (m, v) in g.iter().enumerate() {
        let w = (m as f64).powi(i32::try_from(p).unwrap_or(0));
        sum += v * w;
        scale += v.abs() * w;
    }
    if scale == 0.0 {
        0.0
    } else {
        sum.abs() / scale
    }
}

#[test]
fn filters_are_orthonormal_or_biorthogonal_with_the_advertised_vanishing_moments() {
    for w in FAMILIES {
        let h = w.analysis_low();
        let h_syn = w.synthesis_low();
        let (na, ns) = w.vanishing_moments();
        assert_eq!(h.len(), h_syn.len(), "{w:?}: bank lengths");
        let sum: f64 = h.iter().sum();
        assert!((sum - 2.0_f64.sqrt()).abs() < 1e-12, "{w:?}: sum h = {sum}");
        let sum_syn: f64 = h_syn.iter().sum();
        assert!(
            (sum_syn - 2.0_f64.sqrt()).abs() < 1e-12,
            "{w:?}: sum h_syn = {sum_syn}"
        );
        // biorthogonality (orthonormality for the self-dual families)
        let tol = if table_limited(w) {
            1e-11
        } else {
            1e-12
        };
        for shift in 0..=h.len() / 2 {
            let dot: f64 = (0..h.len())
                .filter(|m| m + 2 * shift < h.len())
                .map(|m| h[m] * h_syn[m + 2 * shift])
                .sum();
            let want = if shift == 0 { 1.0 } else { 0.0 };
            assert!(
                (dot - want).abs() < tol,
                "{w:?}: <h, h_syn(. + 2*{shift})> = {dot}"
            );
        }
        if w.is_orthogonal() {
            for (name, g) in [
                ("analysis", w.analysis_high()),
                ("synthesis", w.synthesis_high()),
            ] {
                for shift in 0..=g.len() / 2 {
                    let dot: f64 = (0..g.len())
                        .filter(|m| m + 2 * shift < g.len())
                        .map(|m| g[m] * g[m + 2 * shift])
                        .sum();
                    let want = if shift == 0 { 1.0 } else { 0.0 };
                    assert!(
                        (dot - want).abs() < tol,
                        "{w:?} {name}: <g, g(. + 2*{shift})> = {dot}"
                    );
                }
            }
            // the two banks are orthogonal to each other
            let g = w.analysis_high();
            for shift in (0..h.len()).step_by(2) {
                let cross: f64 = (0..h.len())
                    .filter(|m| m + shift < h.len())
                    .map(|m| g[m] * h[m + shift])
                    .sum();
                assert!(cross.abs() < tol, "{w:?}: <g, h(. + {shift})> = {cross}");
            }
        }
        // the vanishing moments, counted rather than assumed: the first `na` moments of the
        // analysis high-pass vanish and the next one does not
        let ga = w.analysis_high();
        let gs = w.synthesis_high();
        for p in 0..na {
            let rel = moment_residual(&ga, p);
            assert!(
                rel < 1e-10,
                "{w:?}: analysis moment {p} is {rel:.3e} of its term scale"
            );
        }
        assert!(
            moment_residual(&ga, na) > 1e-8,
            "{w:?}: analysis moment {na} should not vanish (relative {:.3e})",
            moment_residual(&ga, na)
        );
        for p in 0..ns {
            let rel = moment_residual(&gs, p);
            assert!(
                rel < 1e-10,
                "{w:?}: synthesis moment {p} is {rel:.3e} of its term scale"
            );
        }
        assert!(
            moment_residual(&gs, ns) > 1e-8,
            "{w:?}: synthesis moment {ns} should not vanish (relative {:.3e})",
            moment_residual(&gs, ns)
        );
    }
}

#[test]
fn the_high_pass_banks_are_the_alternating_mirror_of_the_dual_low_pass() {
    for w in FAMILIES {
        let (h, g, h_syn, g_syn) = (
            w.analysis_low(),
            w.analysis_high(),
            w.synthesis_low(),
            w.synthesis_high(),
        );
        let l = h.len();
        for k in 0..l {
            let want_a = if k % 2 == 0 {
                h_syn[l - 1 - k]
            } else {
                -h_syn[l - 1 - k]
            };
            let want_s = if k % 2 == 0 {
                h[l - 1 - k]
            } else {
                -h[l - 1 - k]
            };
            assert!((g[k] - want_a).abs() < 1e-16, "{w:?}: g~[{k}]");
            assert!((g_syn[k] - want_s).abs() < 1e-16, "{w:?}: g[{k}]");
        }
    }
}

/// Support `[lo, hi]` of a filter, and whether it is symmetric (`+1`), antisymmetric (`-1`) or
/// neither (`0`) about the middle of that support.
fn symmetry(g: &[f64]) -> (usize, usize, i8) {
    let lo = g.iter().position(|v| *v != 0.0).unwrap_or(0);
    let hi = g.iter().rposition(|v| *v != 0.0).unwrap_or(0);
    let sym = (0..=hi - lo).all(|t| (g[lo + t] - g[hi - t]).abs() < 1e-15);
    let anti = (0..=hi - lo).all(|t| (g[lo + t] + g[hi - t]).abs() < 1e-15);
    let kind = if sym {
        1
    } else if anti {
        -1
    } else {
        0
    };
    (lo, hi, kind)
}

/// Spread of the group delay over the open frequency interval, in samples.
fn group_delay_spread(g: &[f64]) -> f64 {
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for step in 0..64 {
        let w = std::f64::consts::PI * (step as f64 + 0.5) / 64.0;
        let mut num = Complex64::new(0.0, 0.0);
        let mut den = Complex64::new(0.0, 0.0);
        for (n, v) in g.iter().enumerate() {
            let e = Complex64::from_polar(1.0, -w * n as f64);
            num += e * (n as f64) * v;
            den += e * v;
        }
        let tau = (num / den).re;
        lo = lo.min(tau);
        hi = hi.max(tau);
    }
    hi - lo
}

#[test]
fn the_biorthogonal_filters_are_symmetric_and_the_symlets_nearly_linear_phase() {
    for w in FAMILIES.iter().filter(|w| !w.is_orthogonal()) {
        for (name, bank) in [
            ("analysis low", w.analysis_low()),
            ("analysis high", w.analysis_high()),
            ("synthesis low", w.synthesis_low()),
            ("synthesis high", w.synthesis_high()),
        ] {
            // exact (anti)symmetry about the middle of the support: that is linear phase, so the
            // group delay is flat by construction (it is not measured here because the ratio is
            // ill-conditioned near the zeros of the response)
            let (lo, hi, kind) = symmetry(&bank);
            assert_ne!(kind, 0, "{w:?} {name}: not (anti)symmetric on [{lo}, {hi}]");
        }
    }
    // symlets are only approximately symmetric, but their group delay is much flatter than the
    // Daubechies filters of the same length
    for w in [Wavelet::Sym4, Wavelet::Sym5, Wavelet::Sym6] {
        let spread = group_delay_spread(&w.analysis_low());
        assert!(spread < 1.0, "{w:?}: group delay spread {spread}");
    }
    assert!(
        group_delay_spread(&Wavelet::Sym4.analysis_low())
            < 0.5 * group_delay_spread(&Wavelet::Db4.analysis_low())
    );
}

#[test]
fn one_level_is_an_orthogonal_transform_with_exact_reconstruction() {
    for w in FAMILIES.iter().filter(|w| w.is_orthogonal()) {
        let x = signal(64);
        let (a, d) = dwt(&x, *w).unwrap();
        assert_eq!((a.len(), d.len()), (32, 32));
        let energy: f64 = x.iter().map(|v| v * v).sum();
        let coeff: f64 = a.iter().chain(&d).map(|v| v * v).sum();
        assert!((energy - coeff).abs() < 1e-10 * energy, "{w:?}: Parseval");
        let back = idwt(&a, &d, *w).unwrap();
        assert!(
            max_diff(&x, &back) < round_trip_tolerance(*w),
            "{w:?}: reconstruction"
        );
    }
    // a constant has no detail; an interior linear ramp has none for db2 and up
    let (_, d) = dwt(&[3.0; 16], Wavelet::Db2).unwrap();
    assert!(d.iter().all(|v| v.abs() < 1e-12));
    let ramp: Vec<f64> = (0..32).map(f64::from).collect();
    let (_, d) = dwt(&ramp, Wavelet::Db2).unwrap();
    assert!(d[..d.len() - 2].iter().all(|v| v.abs() < 1e-10), "{d:?}");
}

#[test]
fn every_family_and_boundary_mode_reconstructs_one_level() {
    for w in FAMILIES {
        for b in MODES {
            for n in LENGTHS {
                let x = signal(n);
                let (a, d) = dwt_with(&x, w, b).unwrap();
                assert_eq!(a.len(), d.len(), "{w:?} {b:?} {n}");
                assert_eq!(a.len(), w.coeff_len(n, b), "{w:?} {b:?} {n}: count");
                let back = idwt_with(&a, &d, w, b, n).unwrap();
                let err = max_diff(&x, &back);
                assert!(
                    err < round_trip_tolerance(w),
                    "{w:?} {b:?} n = {n}: reconstruction error {err}"
                );
            }
        }
    }
}

#[test]
fn every_family_and_boundary_mode_reconstructs_multilevel() {
    for w in FAMILIES {
        for b in MODES {
            for n in LENGTHS {
                let x = signal(n);
                let levels = if n >= 16 { 3 } else { 2 };
                let c = wavedec_with(&x, w, levels, b).unwrap();
                assert_eq!(c.len(), levels + 1, "{w:?} {b:?} {n}");
                let back = waverec_with(&c, w, b, n).unwrap();
                let err = max_diff(&x, &back);
                assert!(
                    err < 10.0 * round_trip_tolerance(w),
                    "{w:?} {b:?} n = {n} levels = {levels}: error {err}"
                );
            }
        }
    }
    // the periodised transform halves every level; the extension modes shrink by about n/2 once
    // the signal is long compared to the filter
    for w in FAMILIES {
        let c = wavedec_with(&signal(100), w, 3, Boundary::Symmetric).unwrap();
        assert!(c[0].len() < c[3].len(), "{w:?}: multilevel does not shrink");
        let p = wavedec_with(&signal(100), w, 3, Boundary::Periodic).unwrap();
        assert_eq!(p[0].len(), 13, "{w:?}");
    }
}

#[test]
fn the_periodic_mode_is_the_legacy_transform() {
    for w in FAMILIES {
        for n in [8_usize, 16, 32, 100] {
            let x = signal(n);
            let (a, d) = dwt(&x, w).unwrap();
            let (a2, d2) = dwt_with(&x, w, Boundary::Periodic).unwrap();
            assert_eq!(a, a2, "{w:?} {n}");
            assert_eq!(d, d2, "{w:?} {n}");
        }
        // an odd length is supported by dwt_with only
        assert_eq!(dwt(&signal(7), w), Err(WaveletError::OddLength));
        let (a, d) = dwt_with(&signal(7), w, Boundary::Periodic).unwrap();
        assert_eq!(a.len(), 4);
        assert!(
            max_diff(
                &signal(7),
                &idwt_with(&a, &d, w, Boundary::Periodic, 7).unwrap()
            ) < 1e-9
        );
    }
}

#[test]
fn the_boundary_extensions_are_the_documented_mirrors() {
    let x = [1.0, 2.0, 3.0, 4.0, 5.0];
    let periodic: Vec<f64> = (-7..=7)
        .map(|i| extend(&x, i, Boundary::Periodic))
        .collect();
    let zero: Vec<f64> = (-7..=7).map(|i| extend(&x, i, Boundary::Zero)).collect();
    let symmetric: Vec<f64> = (-7..=7)
        .map(|i| extend(&x, i, Boundary::Symmetric))
        .collect();
    let reflect: Vec<f64> = (-7..=7).map(|i| extend(&x, i, Boundary::Reflect)).collect();
    // period 5
    assert_eq!(
        periodic,
        [
            4.0, 5.0, 1.0, 2.0, 3.0, 4.0, 5.0, 1.0, 2.0, 3.0, 4.0, 5.0, 1.0, 2.0, 3.0
        ]
    );
    // nothing outside the signal
    assert_eq!(
        zero,
        [
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 0.0, 0.0, 0.0
        ]
    );
    // half-sample: the edge sample is repeated, period 10
    assert_eq!(
        symmetric,
        [
            4.0, 5.0, 5.0, 4.0, 3.0, 2.0, 1.0, 1.0, 2.0, 3.0, 4.0, 5.0, 5.0, 4.0, 3.0
        ]
    );
    // whole-sample: no repeated edge sample, period 8
    assert_eq!(
        reflect,
        [
            2.0, 3.0, 4.0, 5.0, 4.0, 3.0, 2.0, 1.0, 2.0, 3.0, 4.0, 5.0, 4.0, 3.0, 2.0
        ]
    );
    // a single sample extends to itself in every mirroring mode
    for b in [Boundary::Periodic, Boundary::Symmetric, Boundary::Reflect] {
        for i in -3..=3 {
            assert_eq!(extend(&[7.0], i, b), 7.0, "{b:?} at {i}");
        }
    }
    assert_eq!(extend(&[7.0], 0, Boundary::Zero), 7.0);
    assert_eq!(extend(&[7.0], 1, Boundary::Zero), 0.0);
    assert_eq!(extend(&[], 0, Boundary::Periodic), 0.0);
}

#[test]
fn the_two_dimensional_transform_reconstructs_and_matches_the_separable_reference() {
    let w = Wavelet::Bior3_3;
    let b = Boundary::Symmetric;
    let (rows, cols) = (5_usize, 8_usize);
    let x = signal(rows * cols);
    let packed = dwt2(&x, rows, cols, w, b).unwrap();
    let ca = w.coeff_len(cols, b);
    let ra = w.coeff_len(rows, b);
    assert_eq!(packed.len(), 4 * ra * ca);
    // reference: rows first, then columns, written out independently of dwt2
    let mut row_low = vec![vec![0.0; ca]; rows];
    let mut row_high = vec![vec![0.0; ca]; rows];
    for r in 0..rows {
        let row: Vec<f64> = (0..cols).map(|c| x[r * cols + c]).collect();
        let (a, d) = dwt_with(&row, w, b).unwrap();
        row_low[r].copy_from_slice(&a);
        row_high[r].copy_from_slice(&d);
    }
    let column = |m: &Vec<Vec<f64>>, c: usize| -> Vec<f64> { (0..rows).map(|r| m[r][c]).collect() };
    let mut want = vec![vec![0.0; 2 * ca]; 2 * ra];
    for c in 0..ca {
        let (ll, lh) = dwt_with(&column(&row_low, c), w, b).unwrap();
        let (hl, hh) = dwt_with(&column(&row_high, c), w, b).unwrap();
        for i in 0..ra {
            want[i][c] = ll[i];
            want[i][ca + c] = lh[i];
            want[ra + i][c] = hl[i];
            want[ra + i][ca + c] = hh[i];
        }
    }
    for (r, line) in want.iter().enumerate() {
        for (c, v) in line.iter().enumerate() {
            assert!(
                (packed[r * 2 * ca + c] - v).abs() < 1e-14,
                "sub-band layout at ({r}, {c})"
            );
        }
    }
    let back = idwt2(&packed, rows, cols, w, b).unwrap();
    assert!(max_diff(&x, &back) < 1e-10, "2-D reconstruction");
    // rectangular grids of every family and mode, including non-dyadic ones
    for w in FAMILIES {
        for b in MODES {
            for (rows, cols) in [(3_usize, 5_usize), (4, 4), (7, 6), (1, 9), (9, 1)] {
                let x = signal(rows * cols);
                let packed = dwt2(&x, rows, cols, w, b).unwrap();
                let back = idwt2(&packed, rows, cols, w, b).unwrap();
                let err = max_diff(&x, &back);
                assert!(
                    err < 10.0 * round_trip_tolerance(w),
                    "{w:?} {b:?} {rows}x{cols}: 2-D error {err}"
                );
            }
        }
    }
}

#[test]
fn the_biorthogonal_coefficients_are_the_hand_computed_cdf_values() {
    let s2 = std::f64::consts::SQRT_2;
    // bior2.2: h~ = [-1, 2, 6, 2, -1] / (4 sqrt 2), zero padded to six taps
    let h = Wavelet::Bior2_2.analysis_low();
    let want = [
        -1.0 / (4.0 * s2),
        2.0 / (4.0 * s2),
        6.0 / (4.0 * s2),
        2.0 / (4.0 * s2),
        -1.0 / (4.0 * s2),
        0.0,
    ];
    for (k, v) in want.iter().enumerate() {
        assert!((h[k] - v).abs() < 1e-15, "bior2.2 h~[{k}] = {}", h[k]);
    }
    // bior2.2 synthesis low-pass is the 3-tap [1, 2, 1] / (2 sqrt 2) spline
    let h_syn = Wavelet::Bior2_2.synthesis_low();
    let want_syn = [
        0.0,
        1.0 / (2.0 * s2),
        2.0 / (2.0 * s2),
        1.0 / (2.0 * s2),
        0.0,
        0.0,
    ];
    for (k, v) in want_syn.iter().enumerate() {
        assert!((h_syn[k] - v).abs() < 1e-15, "bior2.2 h[{k}]");
    }
    // bior1.3: h~ = [-1, 1, 8, 8, 1, -1] / (8 sqrt 2), the synthesis bank is Haar
    let h13 = Wavelet::Bior1_3.analysis_low();
    let want13 = [-1.0, 1.0, 8.0, 8.0, 1.0, -1.0];
    for (k, v) in want13.iter().enumerate() {
        assert!((h13[k] - v / (8.0 * s2)).abs() < 1e-15, "bior1.3 h~[{k}]");
    }
    // a hand-computed transform value: for x = e_0 the circular bior2.2 approximation is the
    // filter itself read at 2k mod n, so a[0] = h~[0] and a[3] = h~[2]
    let n = 8;
    let mut x = vec![0.0; n];
    x[0] = 1.0;
    let (a, d) = dwt(&x, Wavelet::Bior2_2).unwrap();
    assert!((a[0] + 1.0 / (4.0 * s2)).abs() < 1e-15, "a[0] = {}", a[0]);
    assert!((a[3] - 6.0 / (4.0 * s2)).abs() < 1e-15, "a[3] = {}", a[3]);
    assert!((d[0]).abs() < 1e-15, "d[0] = {}", d[0]);
    assert!(max_diff(&x, &idwt(&a, &d, Wavelet::Bior2_2).unwrap()) < 1e-15);
}

#[test]
fn invalid_inputs_are_errors() {
    assert_eq!(
        dwt(&[1.0, 2.0, 3.0], Wavelet::Haar),
        Err(WaveletError::OddLength)
    );
    assert_eq!(dwt(&[], Wavelet::Haar), Err(WaveletError::OddLength));
    assert_eq!(
        dwt_with(&[], Wavelet::Haar, Boundary::Symmetric),
        Err(WaveletError::OddLength)
    );
    assert_eq!(
        idwt(&[1.0], &[1.0, 2.0], Wavelet::Haar),
        Err(WaveletError::Shape)
    );
    assert_eq!(
        idwt_with(&[1.0], &[1.0], Wavelet::Haar, Boundary::Zero, 7),
        Err(WaveletError::Shape)
    );
    assert_eq!(
        wavedec(&signal(12), Wavelet::Haar, 3),
        Err(WaveletError::TooManyLevels)
    );
    assert_eq!(
        wavedec(&signal(16), Wavelet::Db4, 3),
        Err(WaveletError::TooManyLevels)
    );
    assert_eq!(
        wavedec_with(&signal(8), Wavelet::Haar, 0, Boundary::Zero),
        Err(WaveletError::TooManyLevels)
    );
    assert_eq!(
        wavedec_with(&signal(3), Wavelet::Haar, 3, Boundary::Zero),
        Err(WaveletError::TooManyLevels)
    );
    assert!(waverec(&[], Wavelet::Haar).is_err());
    assert!(waverec_with(&[], Wavelet::Haar, Boundary::Zero, 8).is_err());
    let c = wavedec_with(&signal(16), Wavelet::Haar, 2, Boundary::Zero).unwrap();
    assert!(waverec_with(&c, Wavelet::Haar, Boundary::Zero, 100).is_err());
    assert_eq!(
        dwt2(&signal(15), 4, 4, Wavelet::Haar, Boundary::Zero),
        Err(WaveletError::Shape)
    );
    assert_eq!(
        dwt2(&signal(16), 0, 4, Wavelet::Haar, Boundary::Zero),
        Err(WaveletError::Shape)
    );
    assert_eq!(
        idwt2(&signal(16), 5, 5, Wavelet::Haar, Boundary::Zero),
        Err(WaveletError::Shape)
    );
}

fn cplx(n: usize) -> Vec<Complex64> {
    let (a, b) = (signal(n), signal(n + 7));
    a.iter()
        .zip(&b[7..])
        .map(|(p, q)| Complex64::new(*p, *q))
        .collect()
}

fn max_diff_complex(
    a: &[Complex64],
    b: &[Complex64],
) -> f64 {
    a.iter()
        .zip(b)
        .map(|(p, q)| (p - q).norm())
        .fold(0.0, f64::max)
}

#[test]
fn the_fractional_fourier_transform_is_a_unitary_group_through_dft_and_reversal() {
    for n in [8, 15, 16, 24] {
        let x = cplx(n);
        // a = 0: identity; a = 1: unitary DFT; a = 2: reversal; a = 4: identity again
        assert!(
            max_diff_complex(&frft(&x, 0.0).unwrap(), &x) < 1e-9,
            "n = {n}: a = 0"
        );
        let dft: Vec<Complex64> = (0..n)
            .map(|j| {
                (0..n)
                    .map(|k| {
                        Complex64::from_polar(
                            1.0,
                            -2.0 * std::f64::consts::PI * (j * k) as f64 / n as f64,
                        ) * x[k]
                    })
                    .sum::<Complex64>()
                    / (n as f64).sqrt()
            })
            .collect();
        assert!(
            max_diff_complex(&frft(&x, 1.0).unwrap(), &dft) < 1e-9,
            "n = {n}: a = 1"
        );
        let rev: Vec<Complex64> = (0..n).map(|j| x[(n - j) % n]).collect();
        assert!(
            max_diff_complex(&frft(&x, 2.0).unwrap(), &rev) < 1e-9,
            "n = {n}: a = 2"
        );
        assert!(
            max_diff_complex(&frft(&x, 4.0).unwrap(), &x) < 1e-9,
            "n = {n}: a = 4"
        );
        // index additivity and unitarity at fractional orders
        let step = frft(&frft(&x, 0.3).unwrap(), 0.45).unwrap();
        assert!(
            max_diff_complex(&step, &frft(&x, 0.75).unwrap()) < 1e-9,
            "n = {n}: additivity"
        );
        let y = frft(&x, 0.37).unwrap();
        let (e0, e1): (f64, f64) = (
            x.iter().map(|v| v.norm_sqr()).sum(),
            y.iter().map(|v| v.norm_sqr()).sum(),
        );
        assert!((e0 - e1).abs() < 1e-9 * e0, "n = {n}: norm");
    }
}

#[test]
fn the_fractional_matrix_is_unitary_and_validates_the_size() {
    let n = 12;
    let m = frft_matrix(n, 0.6).unwrap();
    for i in 0..n {
        for j in 0..n {
            let dot: Complex64 = (0..n).map(|k| m[i * n + k] * m[j * n + k].conj()).sum();
            let want = if i == j { 1.0 } else { 0.0 };
            assert!((dot - want).norm() < 1e-9, "({i}, {j}): {dot}");
        }
    }
    assert!(frft_matrix(1, 0.5).is_err());
}

// -------------------------------------------------------------------------
// Multilevel 2-D transform
// -------------------------------------------------------------------------

/// The `LL` band (top-left `ra x ca` block) of a [`dwt2`] packing.
fn ll_band(
    packed: &[f64],
    ra: usize,
    ca: usize,
) -> Vec<f64> {
    let out_cols = 2 * ca;
    let mut ll = Vec::with_capacity(ra * ca);
    for r in 0..ra {
        for c in 0..ca {
            ll.push(packed[r * out_cols + c]);
        }
    }
    ll
}

fn array_energy(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum()
}

#[test]
fn multilevel_two_dimensional_round_trip_on_dyadic_shapes() {
    let families = [Wavelet::Haar, Wavelet::Db2, Wavelet::Db4, Wavelet::Bior2_2];
    for (rows, cols, levels) in [(8usize, 8usize, 2usize), (16, 16, 3)] {
        let x = signal(rows * cols);
        for w in families {
            let c = wavedec2(&x, rows, cols, w, levels).unwrap();
            assert_eq!(c.len(), levels + 1, "{w:?}: coefficient count");
            // Shape bookkeeping: LL_L is rows/2^L x cols/2^L, and packed_j has
            // the size of the array level j was applied to.
            assert_eq!(
                c[0].len(),
                (rows >> levels) * (cols >> levels),
                "{w:?}: LL shape"
            );
            for j in 1..=levels {
                let scale = 1usize << (j - 1);
                assert_eq!(
                    c[levels - j + 1].len(),
                    (rows / scale) * (cols / scale),
                    "{w:?}: packed level {j}"
                );
            }
            let back = waverec2(&c, rows, cols, w).unwrap();
            let err = max_diff(&x, &back);
            assert!(
                err < round_trip_tolerance(w),
                "{w:?} {rows}x{cols} levels = {levels}: error {err}"
            );
            // The explicit periodic variant must agree with wavedec2.
            let c2 = wavedec2_with(&x, rows, cols, w, levels, Boundary::Periodic).unwrap();
            assert_eq!(c, c2, "{w:?}: periodic variant");
        }
    }
}

#[test]
fn multilevel_two_dimensional_layout_is_the_recursive_one_level_transform() {
    let (rows, cols, levels) = (16usize, 16usize, 3usize);
    let x = signal(rows * cols);
    for w in [Wavelet::Haar, Wavelet::Db2, Wavelet::Db4, Wavelet::Bior3_3] {
        let c = wavedec2(&x, rows, cols, w, levels).unwrap();
        let mut current = x.clone();
        let (mut rr, mut cc) = (rows, cols);
        let mut reference = Vec::new();
        for _ in 0..levels {
            let packed = dwt2(&current, rr, cc, w, Boundary::Periodic).unwrap();
            let ra = w.coeff_len(rr, Boundary::Periodic);
            let ca = w.coeff_len(cc, Boundary::Periodic);
            reference.push(packed.clone());
            current = ll_band(&packed, ra, ca);
            rr = ra;
            cc = ca;
        }
        reference.reverse();
        assert_eq!(c[0], current, "{w:?}: final approximation");
        for (k, want) in reference.iter().enumerate() {
            assert_eq!(&c[k + 1], want, "{w:?}: packed level {}", levels - k);
        }
    }
}

#[test]
fn multilevel_two_dimensional_parseval_per_level() {
    let (rows, cols, levels) = (16usize, 16usize, 3usize);
    let x = signal(rows * cols);
    for w in [Wavelet::Haar, Wavelet::Db2, Wavelet::Db4] {
        let c = wavedec2(&x, rows, cols, w, levels).unwrap();
        let mut input_energy = array_energy(&x);
        let (mut rr, mut cc) = (rows, cols);
        for j in 1..=levels {
            let idx = levels - j + 1;
            let packed = &c[idx];
            let packed_energy = array_energy(packed);
            assert!(
                (packed_energy - input_energy).abs() <= 1e-10 * input_energy.max(1.0),
                "{w:?}: Parseval at level {j}: {packed_energy} vs {input_energy}"
            );
            let ra = w.coeff_len(rr, Boundary::Periodic);
            let ca = w.coeff_len(cc, Boundary::Periodic);
            let ll = ll_band(packed, ra, ca);
            input_energy = array_energy(&ll);
            rr = ra;
            cc = ca;
        }
        assert_eq!(c[0].len(), rr * cc, "{w:?}: final LL shape");
        assert!(
            (array_energy(&c[0]) - input_energy).abs() < 1e-10 * input_energy.max(1.0),
            "{w:?}: final LL energy"
        );
    }
}

#[test]
fn multilevel_two_dimensional_extension_modes_round_trip_and_stop() {
    for (w, b) in [
        (Wavelet::Db4, Boundary::Symmetric),
        (Wavelet::Bior3_3, Boundary::Reflect),
        (Wavelet::Haar, Boundary::Zero),
    ] {
        let (rows, cols) = (32usize, 32usize);
        let x = signal(rows * cols);
        let c = wavedec2_with(&x, rows, cols, w, 2, b).unwrap();
        assert_eq!(c.len(), 3, "{w:?} {b:?}: coefficient count");
        let back = waverec2_with(&c, rows, cols, w, b).unwrap();
        let err = max_diff(&x, &back);
        assert!(
            err < 10.0 * round_trip_tolerance(w),
            "{w:?} {b:?}: error {err}"
        );
        assert!(c[0].len() < c[2].len(), "{w:?} {b:?}: shrinking");
    }
    // db4 has 8 taps: its symmetric coefficient count at 8x8 is 7, so a second
    // level would not shrink the band and the recursion reports it.
    assert_eq!(
        wavedec2_with(&signal(64), 8, 8, Wavelet::Db4, 2, Boundary::Symmetric),
        Err(WaveletError::TooManyLevels)
    );
}

#[test]
fn multilevel_two_dimensional_invalid_inputs() {
    assert_eq!(
        wavedec2(&signal(15), 4, 4, Wavelet::Haar, 2),
        Err(WaveletError::Shape)
    );
    assert_eq!(
        wavedec2(&[], 0, 4, Wavelet::Haar, 2),
        Err(WaveletError::Shape)
    );
    assert_eq!(
        wavedec2(&signal(64), 8, 8, Wavelet::Haar, 0),
        Err(WaveletError::TooManyLevels)
    );
    // beyond the dyadic depth
    assert_eq!(
        wavedec2(&signal(64), 8, 8, Wavelet::Haar, 4),
        Err(WaveletError::TooManyLevels)
    );
    // non-dyadic shape
    assert_eq!(
        wavedec2(&signal(36), 6, 6, Wavelet::Haar, 2),
        Err(WaveletError::TooManyLevels)
    );
    let c = wavedec2(&signal(64), 8, 8, Wavelet::Haar, 2).unwrap();
    assert_eq!(waverec2(&c, 16, 8, Wavelet::Haar), Err(WaveletError::Shape));
    assert_eq!(waverec2(&[], 8, 8, Wavelet::Haar), Err(WaveletError::Shape));
    assert_eq!(
        waverec2(&c[..2], 8, 8, Wavelet::Haar),
        Err(WaveletError::Shape)
    );
    let mut bad_ll = c.clone();
    bad_ll[0].push(0.0);
    assert_eq!(
        waverec2(&bad_ll, 8, 8, Wavelet::Haar),
        Err(WaveletError::Shape)
    );
    let mut bad_packed = c.clone();
    bad_packed[1].pop();
    assert_eq!(
        waverec2(&bad_packed, 8, 8, Wavelet::Haar),
        Err(WaveletError::Shape)
    );
}
