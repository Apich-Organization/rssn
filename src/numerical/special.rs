//! # Numerical Special Functions
//!
//! This module provides numerical implementations of mathematical special functions.
//! It includes gamma functions, beta functions, error functions, Bessel functions,
//! orthogonal polynomials, and other commonly used special functions.

use statrs::function::beta::beta;
use statrs::function::beta::ln_beta;
use statrs::function::erf::erf;
use statrs::function::erf::erfc;
use statrs::function::gamma::digamma;
use statrs::function::gamma::gamma;
use statrs::function::gamma::ln_gamma;

// ============================================================================
// Gamma Functions
// ============================================================================

/// Computes the gamma function, Γ(x).
/// Γ(n) = (n-1)! for positive integers.
#[must_use]

pub fn gamma_numerical(x: f64) -> f64 {

    gamma(x)
}

/// Computes the natural logarithm of the gamma function, ln(Γ(x)).
#[must_use]

pub fn ln_gamma_numerical(
    x: f64
) -> f64 {

    ln_gamma(x)
}

/// Computes the digamma function, ψ(x) = d/dx ln(Γ(x)) = Γ'(x)/Γ(x).
#[must_use]

pub fn digamma_numerical(
    x: f64
) -> f64 {

    digamma(x)
}

/// Computes the lower incomplete gamma function, γ(s, x) = ∫₀ˣ t^(s-1) e^(-t) dt.
/// Uses series expansion for computation.
#[must_use]

pub fn lower_incomplete_gamma(
    s: f64,
    x: f64,
) -> f64 {

    if x < 0.0 || s <= 0.0 {

        return f64::NAN;
    }

    if x == 0.0 {

        return 0.0;
    }

    // Use series expansion: γ(s,x) = x^s * e^(-x) * Σ(x^n / Γ(s+n+1))
    // Or regularized: P(s,x) = γ(s,x) / Γ(s)
    // Using statrs regularized gamma if available, otherwise series
    let mut sum = 0.0;

    let mut term = 1.0 / s;

    sum += term;

    for n in 1 .. 200 {

        term *= x / (s + f64::from(n));

        sum += term;

        if term.abs()
            < 1e-15 * sum.abs()
        {

            break;
        }
    }

    x.powf(s) * (-x).exp() * sum
}

/// Computes the upper incomplete gamma function, Γ(s, x) = ∫ₓ^∞ t^(s-1) e^(-t) dt.
#[must_use]

pub fn upper_incomplete_gamma(
    s: f64,
    x: f64,
) -> f64 {

    gamma(s)
        - lower_incomplete_gamma(s, x)
}

/// Computes the regularized lower incomplete gamma function, P(s, x) = γ(s, x) / Γ(s).
#[must_use]

pub fn regularized_lower_gamma(
    s: f64,
    x: f64,
) -> f64 {

    lower_incomplete_gamma(s, x)
        / gamma(s)
}

/// Computes the regularized upper incomplete gamma function, Q(s, x) = Γ(s, x) / Γ(s).
#[must_use]

pub fn regularized_upper_gamma(
    s: f64,
    x: f64,
) -> f64 {

    1.0 - regularized_lower_gamma(s, x)
}

// ============================================================================
// Beta Functions
// ============================================================================

/// Computes the beta function, B(a, b) = Γ(a)Γ(b)/Γ(a+b).
#[must_use]

pub fn beta_numerical(
    a: f64,
    b: f64,
) -> f64 {

    beta(a, b)
}

/// Computes the natural logarithm of the beta function, ln(B(a, b)).
#[must_use]

pub fn ln_beta_numerical(
    a: f64,
    b: f64,
) -> f64 {

    ln_beta(a, b)
}

/// Computes the incomplete beta function, B(x; a, b) = ∫₀ˣ t^(a-1) (1-t)^(b-1) dt.
#[must_use]

pub fn incomplete_beta(
    x: f64,
    a: f64,
    b: f64,
) -> f64 {

    if !(0.0 ..= 1.0).contains(&x)
        || a <= 0.0
        || b <= 0.0
    {

        return f64::NAN;
    }

    if x == 0.0 {

        return 0.0;
    }

    if x == 1.0 {

        return beta(a, b);
    }

    // Use series expansion
    regularized_beta(x, a, b)
        * beta(a, b)
}

/// Computes the regularized incomplete beta function, `I_x(a`, b) = B(x; a, b) / B(a, b).
/// Uses continued fraction expansion for better convergence.
#[must_use]

pub fn regularized_beta(
    x: f64,
    a: f64,
    b: f64,
) -> f64 {

    if !(0.0 ..= 1.0).contains(&x)
        || a <= 0.0
        || b <= 0.0
    {

        return f64::NAN;
    }

    if x == 0.0 {

        return 0.0;
    }

    if x == 1.0 {

        return 1.0;
    }

    // Use symmetry for better convergence
    if x > (a + 1.0) / (a + b + 2.0) {

        return 1.0
            - regularized_beta(
                1.0 - x,
                b,
                a,
            );
    }

    // Continued fraction expansion
    let bt = if x == 0.0 || x == 1.0 {

        0.0
    } else {

        (a.mul_add(
            x.ln(),
            ln_gamma(a + b)
                - ln_gamma(a)
                - ln_gamma(b),
        ) + b * (1.0 - x).ln())
        .exp()
    };

    let mut am = 1.0;

    let mut bm = 1.0;

    let mut az = 1.0;

    let qab = a + b;

    let qap = a + 1.0;

    let qam = a - 1.0;

    let mut bz = 1.0 - qab * x / qap;

    for m in 1 .. 200 {

        let em = f64::from(m);

        let tem = em + em;

        let d = em * (b - em) * x
            / ((qam + tem) * (a + tem));

        let ap = az + d * am;

        let bp = bz + d * bm;

        let d = -(a + em)
            * (qab + em)
            * x
            / ((a + tem) * (qap + tem));

        let app = ap + d * az;

        let bpp = bp + d * bz;

        let aold = az;

        am = ap / bpp;

        bm = bp / bpp;

        az = app / bpp;

        bz = 1.0;

        if (az - aold).abs()
            < 1e-14 * az.abs()
        {

            return bt * az / a;
        }
    }

    bt * az / a
}

// ============================================================================
// Error Functions
// ============================================================================

/// Computes the error function, erf(x) = (2/√π) ∫₀ˣ e^(-t²) dt.
#[must_use]

pub fn erf_numerical(x: f64) -> f64 {

    erf(x)
}

/// Computes the complementary error function, erfc(x) = 1 - erf(x).
#[must_use]

pub fn erfc_numerical(x: f64) -> f64 {

    erfc(x)
}

/// Computes the inverse error function, erf⁻¹(x).
/// Uses a rational approximation for the initial guess and Newton-Raphson refinement.
#[must_use]

pub fn inverse_erf_numerical(
    x: f64
) -> f64 {

    if x <= -1.0 {

        return f64::NEG_INFINITY;
    }

    if x >= 1.0 {

        return f64::INFINITY;
    }

    if x.abs() < 1e-15 {

        return 0.0;
    }

    // Use an approximation formula
    let sign = if x < 0.0 {

        -1.0
    } else {

        1.0
    };

    let x = x.abs();

    // Rational approximation for |x| < 0.7
    let result = if x < 0.7 {

        let x2 = x * x;

        let num = x
            * (1.0
                + x2 * x2.mul_add(
                    0.014_000_2,
                    -0.140_543_331,
                ));

        let den = 1.0
            + x2 * x2.mul_add(
                0.049_988,
                -0.453_004_011,
            );

        num / den
    } else {

        // For larger x, use a different approximation
        let y =
            (-(1.0 - x).ln()).sqrt();

        let num = y
            * (1.0
                + y * (-0.094_1
                    + y * 0.003_27));

        let den = 1.0
            + y * (-0.188 + y * 0.0329);

        num / den * 0.886_226_899 // √(π/2)
    };

    // Refine with Newton-Raphson iterations
    let two_over_sqrt_pi = 2.0
        / std::f64::consts::PI.sqrt();

    let mut y = result;

    for _ in 0 .. 3 {

        let err = erf(y) - x;

        let deriv = two_over_sqrt_pi
            * (-y * y).exp();

        y -= err / deriv;
    }

    sign * y
}

// ============================================================================
// Bessel Functions
// ============================================================================

/// Computes the Bessel function of the first kind, J₀(x).
#[must_use]

pub fn bessel_j0(x: f64) -> f64 {

    if x == 0.0 {

        return 1.0;
    }

    let ax = x.abs();

    if ax < 8.0 {

        let y = x * x;

        let ans1 = 57_568_490_574.0
            + y * (-13_362_590_354.0
                + y * (651_619_640.7 + y * (-11_214_424.18 + y * y.mul_add(-184.905_245_6, 77_392.330_17))));

        let ans2 = 57_568_490_411.0
            + y * (1_029_532_985.0
                + y * (9_494_680.718
                    + y * y.mul_add(
                        267.853_271_2
                            + y,
                        59_272.648_53,
                    )));

        ans1 / ans2
    } else {

        let z = 8.0 / ax;

        let y = z * z;

        let xx =
            ax - 0.785_398_163_397_448;

        let ans1 = 1.0
            + y * (-0.109_862_862_7e-2
                + y * (0.273_451_040_7e-4 + y * (-0.207_337_063_9e-5 + y * 0.209_388_721_1e-6)));

        let ans2 = -0.156_249_999_5e-1
            + y * (0.143_048_876_5e-3
                + y * (-0.691_114_765_1e-5 + y * (0.762_109_516_1e-6 - y * 0.934_945_152e-7)));

        (0.636_619_772 / ax).sqrt()
            * xx.cos().mul_add(
                ans1,
                -(z * xx.sin() * ans2),
            )
    }
}

/// Computes the Bessel function of the first kind, J₁(x).
#[must_use]

pub fn bessel_j1(x: f64) -> f64 {

    if x == 0.0 {

        return 0.0;
    }

    let ax = x.abs();

    if ax < 8.0 {

        let y = x * x;

        let ans1 = x
            * (723_626_142_320.0
                + y * (-789_505_923_500.0
                    + y * (242_396_853.1 + y * (-2_972_611.439 + y * y.mul_add(-30.160_366_06, 15_704.482_60)))));

        let ans2 = 144_725_228_442.0
            + y * (230_053_517_800.0
                + y * (185_833_047.4
                    + y * y.mul_add(
                        376.999_139_7
                            + y,
                        994_474.339_4,
                    )));

        ans1 / ans2
    } else {

        let z = 8.0 / ax;

        let y = z * z;

        let xx = ax - 2.356_194_491;

        let ans1 = 1.0
            + y * (0.183_105e-2
                + y * (-0.351_639_649_6e-4 + y * (0.245_752_017_4e-5 + y * (-0.240_337_019e-6))));

        let ans2 = 0.046_874_999_95
            + y * (-0.200_269_087_3e-3
                + y * (0.844_919_909_6e-5 + y * (-0.882_289_87e-6 + y * 0.105_787_412e-6)));

        let ans = (0.636_619_772 / ax)
            .sqrt()
            * xx.cos().mul_add(
                ans1,
                -(z * xx.sin() * ans2),
            );

        if x < 0.0 {

            -ans
        } else {

            ans
        }
    }
}

/// Computes the Bessel function of the second kind, Y₀(x).
#[must_use]

pub fn bessel_y0(x: f64) -> f64 {

    if x < 0.0 {

        return f64::NAN;
    }

    if x < 8.0 {

        let y = x * x;

        let ans1 = -2_957_821_389.0
            + y * (7_062_834_065.0
                + y * (-512_359_803.6 + y * (10_879_881.29 + y * y.mul_add(228.462_273_3, -86_327.927_57))));

        let ans2 = 4_007_654_426.9
            + y * (745_249_964.8
                + y * (7_189_466.438
                    + y * y.mul_add(
                        226.103_024_4
                            + y,
                        474_472.647_0,
                    )));

        (0.636_619_772 * bessel_j0(x))
            .mul_add(
                x.ln(),
                ans1 / ans2,
            )
    } else {

        let z = 8.0 / x;

        let y = z * z;

        let xx =
            x - 0.785_398_163_397_448;

        let ans1 = 1.0
            + y * (-0.109_862_862_7e-2
                + y * (0.273_451_040_7e-4 + y * (-0.207_337_063_9e-5 + y * 0.209_388_721_1e-6)));

        let ans2 = -0.156_249_999_5e-1
            + y * (0.143_048_876_5e-3
                + y * (-0.691_114_765_1e-5 + y * (0.762_109_516_1e-6 + y * (-0.934_945_152e-7))));

        (0.636_619_772 / x).sqrt()
            * xx.sin().mul_add(
                ans1,
                z * xx.cos() * ans2,
            )
    }
}

/// Computes the Bessel function of the second kind, Y₁(x).
#[must_use]

pub fn bessel_y1(x: f64) -> f64 {

    if x < 0.0 {

        return f64::NAN;
    }

    if x < 8.0 {

        let y = x * x;

        let ans1 = x
            * (-49_006_049_430_000_000.0
                + y * (12_752_743_900_000_000.0
                    + y * (-5_153_438_139_000.0
                        + y * (7_349_264_551.0 + y * y.mul_add(85_119_379.35, -42_379_227.26)))));

        let ans2 = 24_995_805_700_000_000.0
            + y * (42_444_196_640_000.0
                + y * (37_336_503_670.0
                    + y * (22_459_040.0 + y * y.mul_add(354_963.288_5 + y, 102_042.605))));

        ans1 / ans2
            + 0.636_619_772
                * bessel_j1(x).mul_add(
                    x.ln(),
                    -(1.0 / x),
                )
    } else {

        let z = 8.0 / x;

        let y = z * z;

        let xx = x - 2.356_194_491;

        let ans1 = 1.0
            + y * (0.183_105e-2
                + y * (-0.351_639_649_6e-4 + y * (0.245_752_017_4e-5 + y * (-0.240_337_019e-6))));

        let ans2 = 0.046_874_999_95
            + y * (-0.200_269_087_3e-3
                + y * (0.844_919_909_6e-5 + y * (-0.882_289_87e-6 + y * 0.105_787_412e-6)));

        (0.636_619_772 / x).sqrt()
            * xx.sin().mul_add(
                ans1,
                z * xx.cos() * ans2,
            )
    }
}

/// Computes the modified Bessel function of the first kind, I₀(x).
#[must_use]

pub fn bessel_i0(x: f64) -> f64 {

    if x == 0.0 {

        return 1.0;
    }

    let ax = x.abs();

    if ax < 3.75 {

        let y = (x / 3.75).powi(2);

        1.0 + y
            * (3.515_622_9
                + y * (3.089_942_4
                    + y * (1.206_749_2 + y * (0.265_973_2 + y * (0.036_076_8 + y * 0.004_581_3)))))
    } else {

        let y = 3.75 / ax;

        (ax.exp() / ax.sqrt())
            * (0.398_942_28
                + y * (0.013_285_92
                    + y * (0.002_253_19
                        + y * (-0.001_575_65
                            + y * (0.009_162_81
                                + y * (-0.020_577_06
                                    + y * (0.026_355_37 + y * (-0.016_476_33 + y * 0.003_923_77))))))))
    }
}

/// Computes the modified Bessel function of the first kind, I₁(x).
#[must_use]

pub fn bessel_i1(x: f64) -> f64 {

    let ax = x.abs();

    let ans = if ax < 3.75 {

        let y = (x / 3.75).powi(2);

        ax * (0.5
            + y * (0.878_905_94
                + y * (0.514_988_69
                    + y * (0.150_849_34 + y * (0.026_587_33 + y * (0.003_015_32 + y * 0.000_324_11))))))
    } else {

        let y = 3.75 / ax;

        let ans = 0.398_942_28
            + y * (-0.039_880_24
                + y * (-0.003_620_18
                    + y * (0.001_638_01
                        + y * (-0.010_315_55
                            + y * (0.022_829_67
                                + y * (-0.028_953_12 + y * (0.017_876_54 - y * 0.004_200_59)))))));

        (ax.exp() / ax.sqrt()) * ans
    };

    if x < 0.0 {

        -ans
    } else {

        ans
    }
}

// ============================================================================
// Orthogonal Polynomials
// ============================================================================

/// Computes the Legendre polynomial `P_n(x)` using recurrence relation.
#[must_use]

pub fn legendre_p(
    n: u32,
    x: f64,
) -> f64 {

    if n == 0 {

        return 1.0;
    }

    if n == 1 {

        return x;
    }

    let mut p_prev = 1.0;

    let mut p_curr = x;

    for k in 2 ..= n {

        let p_next =
            (f64::from(2 * k - 1)
                * x
                * p_curr
                - f64::from(k - 1)
                    * p_prev)
                / f64::from(k);

        p_prev = p_curr;

        p_curr = p_next;
    }

    p_curr
}

/// Computes the Chebyshev polynomial of the first kind `T_n(x)`.
#[must_use]

pub fn chebyshev_t(
    n: u32,
    x: f64,
) -> f64 {

    if n == 0 {

        return 1.0;
    }

    if n == 1 {

        return x;
    }

    let mut t_prev = 1.0;

    let mut t_curr = x;

    for _ in 2 ..= n {

        let t_next = (2.0 * x)
            .mul_add(t_curr, -t_prev);

        t_prev = t_curr;

        t_curr = t_next;
    }

    t_curr
}

/// Computes the Chebyshev polynomial of the second kind `U_n(x)`.
#[must_use]

pub fn chebyshev_u(
    n: u32,
    x: f64,
) -> f64 {

    if n == 0 {

        return 1.0;
    }

    if n == 1 {

        return 2.0 * x;
    }

    let mut u_prev = 1.0;

    let mut u_curr = 2.0 * x;

    for _ in 2 ..= n {

        let u_next = (2.0 * x)
            .mul_add(u_curr, -u_prev);

        u_prev = u_curr;

        u_curr = u_next;
    }

    u_curr
}

/// Computes the (physicists') Hermite polynomial `H_n(x)`.
#[must_use]

pub fn hermite_h(
    n: u32,
    x: f64,
) -> f64 {

    if n == 0 {

        return 1.0;
    }

    if n == 1 {

        return 2.0 * x;
    }

    let mut h_prev = 1.0;

    let mut h_curr = 2.0 * x;

    for k in 2 ..= n {

        let h_next = 2.0 * x * h_curr
            - 2.0
                * f64::from(k - 1)
                * h_prev;

        h_prev = h_curr;

        h_curr = h_next;
    }

    h_curr
}

/// Computes the Laguerre polynomial `L_n(x)`.
#[must_use]

pub fn laguerre_l(
    n: u32,
    x: f64,
) -> f64 {

    if n == 0 {

        return 1.0;
    }

    if n == 1 {

        return 1.0 - x;
    }

    let mut l_prev = 1.0;

    let mut l_curr = 1.0 - x;

    for k in 2 ..= n {

        let l_next =
            (f64::from(2 * k - 1) - x)
                * l_curr
                / f64::from(k)
                - f64::from(k - 1)
                    * l_prev
                    / f64::from(k);

        l_prev = l_curr;

        l_curr = l_next;
    }

    l_curr
}

// ============================================================================
// Other Special Functions
// ============================================================================

/// Computes the factorial n!
#[must_use]

pub fn factorial(n: u64) -> f64 {

    if n <= 1 {

        return 1.0;
    }

    gamma((n + 1) as f64)
}

/// Computes the double factorial n!!
/// n!! = n * (n-2) * (n-4) * ... * (1 or 2)
#[must_use]

pub fn double_factorial(n: u64) -> f64 {

    if n <= 1 {

        return 1.0;
    }

    let mut result = 1.0;

    let mut k = n;

    while k > 1 {

        result *= k as f64;

        k -= 2;
    }

    result
}

/// Computes the binomial coefficient C(n, k) = n! / (k! * (n-k)!)
#[must_use]

pub fn binomial(
    n: u64,
    k: u64,
) -> f64 {

    if k > n {

        return 0.0;
    }

    if k == 0 || k == n {

        return 1.0;
    }

    // Use gamma for large values to avoid overflow
    gamma((n + 1) as f64)
        / (gamma((k + 1) as f64)
            * gamma((n - k + 1) as f64))
}

/// Computes the Riemann zeta function ζ(s) for real s > 1.
/// Uses the Euler product formula for computation.
#[must_use]

pub fn riemann_zeta(s: f64) -> f64 {

    if s <= 1.0 {

        if s == 1.0 {

            return f64::INFINITY;
        }

        // For s < 1, use reflection formula (simplified)
        return f64::NAN; // TODO: implement for s < 1
    }

    // Direct summation for s > 1
    let mut sum = 0.0;

    for n in 1 .. 10000 {

        let term =
            1.0 / f64::from(n).powf(s);

        sum += term;

        if term < 1e-15 * sum.abs() {

            break;
        }
    }

    sum
}

/// Computes the sinc function sinc(x) = sin(πx) / (πx).
#[must_use]

pub fn sinc(x: f64) -> f64 {

    if x.abs() < 1e-10 {

        return 1.0;
    }

    let px = std::f64::consts::PI * x;

    px.sin() / px
}

/// Computes the logit function logit(p) = ln(p / (1-p)).
#[must_use]

pub fn logit(p: f64) -> f64 {

    if p <= 0.0 || p >= 1.0 {

        return f64::NAN;
    }

    (p / (1.0 - p)).ln()
}

/// Computes the logistic (sigmoid) function σ(x) = 1 / (1 + e^(-x)).
#[must_use]

pub fn sigmoid(x: f64) -> f64 {

    1.0 / (1.0 + (-x).exp())
}

/// Computes the softplus function softplus(x) = ln(1 + e^x).
#[must_use]

pub fn softplus(x: f64) -> f64 {

    if x > 20.0 {

        x // Avoid overflow
    } else if x < -20.0 {

        x.exp()
    } else {

        x.exp().ln_1p()
    }
}
