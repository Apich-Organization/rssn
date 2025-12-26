//! # Numerical Special Functions
//!
//! This module provides numerical implementations of mathematical special functions.
//! It includes gamma functions, beta functions, error functions, Bessel functions,
//! orthogonal polynomials, and other commonly used special functions.

use statrs::function::beta::{beta, ln_beta};
use statrs::function::erf::{erf, erfc};
use statrs::function::gamma::{digamma, gamma, ln_gamma};

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

pub fn ln_gamma_numerical(x: f64) -> f64 {

    ln_gamma(x)
}

/// Computes the digamma function, ψ(x) = d/dx ln(Γ(x)) = Γ'(x)/Γ(x).
#[must_use]

pub fn digamma_numerical(x: f64) -> f64 {

    digamma(x)
}

/// Computes the lower incomplete gamma function, γ(s, x) = ∫₀ˣ t^(s-1) e^(-t) dt.
/// Uses series expansion for computation.
#[must_use]

pub fn lower_incomplete_gamma(s: f64, x: f64) -> f64 {

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

    for n in 1..200 {

        term *= x / (s + n as f64);

        sum += term;

        if term.abs() < 1e-15 * sum.abs() {

            break;
        }
    }

    x.powf(s) * (-x).exp() * sum
}

/// Computes the upper incomplete gamma function, Γ(s, x) = ∫ₓ^∞ t^(s-1) e^(-t) dt.
#[must_use]

pub fn upper_incomplete_gamma(s: f64, x: f64) -> f64 {

    gamma(s) - lower_incomplete_gamma(s, x)
}

/// Computes the regularized lower incomplete gamma function, P(s, x) = γ(s, x) / Γ(s).
#[must_use]

pub fn regularized_lower_gamma(s: f64, x: f64) -> f64 {

    lower_incomplete_gamma(s, x) / gamma(s)
}

/// Computes the regularized upper incomplete gamma function, Q(s, x) = Γ(s, x) / Γ(s).
#[must_use]

pub fn regularized_upper_gamma(s: f64, x: f64) -> f64 {

    1.0 - regularized_lower_gamma(s, x)
}

// ============================================================================
// Beta Functions
// ============================================================================

/// Computes the beta function, B(a, b) = Γ(a)Γ(b)/Γ(a+b).
#[must_use]

pub fn beta_numerical(a: f64, b: f64) -> f64 {

    beta(a, b)
}

/// Computes the natural logarithm of the beta function, ln(B(a, b)).
#[must_use]

pub fn ln_beta_numerical(a: f64, b: f64) -> f64 {

    ln_beta(a, b)
}

/// Computes the incomplete beta function, B(x; a, b) = ∫₀ˣ t^(a-1) (1-t)^(b-1) dt.
#[must_use]

pub fn incomplete_beta(x: f64, a: f64, b: f64) -> f64 {

    if x < 0.0 || x > 1.0 || a <= 0.0 || b <= 0.0 {

        return f64::NAN;
    }

    if x == 0.0 {

        return 0.0;
    }

    if x == 1.0 {

        return beta(a, b);
    }

    // Use series expansion
    regularized_beta(x, a, b) * beta(a, b)
}

/// Computes the regularized incomplete beta function, I_x(a, b) = B(x; a, b) / B(a, b).
/// Uses continued fraction expansion for better convergence.
#[must_use]

pub fn regularized_beta(x: f64, a: f64, b: f64) -> f64 {

    if x < 0.0 || x > 1.0 || a <= 0.0 || b <= 0.0 {

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

        return 1.0 - regularized_beta(1.0 - x, b, a);
    }

    // Continued fraction expansion
    let bt = if x == 0.0 || x == 1.0 {

        0.0
    } else {

        (ln_gamma(a + b) - ln_gamma(a) - ln_gamma(b) + a * x.ln() + b * (1.0 - x).ln()).exp()
    };

    let mut am = 1.0;

    let mut bm = 1.0;

    let mut az = 1.0;

    let qab = a + b;

    let qap = a + 1.0;

    let qam = a - 1.0;

    let mut bz = 1.0 - qab * x / qap;

    for m in 1..200 {

        let em = m as f64;

        let tem = em + em;

        let d = em * (b - em) * x / ((qam + tem) * (a + tem));

        let ap = az + d * am;

        let bp = bz + d * bm;

        let d = -(a + em) * (qab + em) * x / ((a + tem) * (qap + tem));

        let app = ap + d * az;

        let bpp = bp + d * bz;

        let aold = az;

        am = ap / bpp;

        bm = bp / bpp;

        az = app / bpp;

        bz = 1.0;

        if (az - aold).abs() < 1e-14 * az.abs() {

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

pub fn inverse_erf_numerical(x: f64) -> f64 {

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
    let sign = if x < 0.0 { -1.0 } else { 1.0 };

    let x = x.abs();

    // Rational approximation for |x| < 0.7
    let result = if x < 0.7 {

        let x2 = x * x;

        let num = x * (1.0 + x2 * (-0.140543331 + x2 * 0.0140002));

        let den = 1.0 + x2 * (-0.453004011 + x2 * 0.049988);

        num / den
    } else {

        // For larger x, use a different approximation
        let y = (-(1.0 - x).ln()).sqrt();

        let num = y * (1.0 + y * (-0.0941 + y * 0.00327));

        let den = 1.0 + y * (-0.188 + y * 0.0329);

        num / den * 0.886226899 // √(π/2)
    };

    // Refine with Newton-Raphson iterations
    let two_over_sqrt_pi = 2.0 / std::f64::consts::PI.sqrt();

    let mut y = result;

    for _ in 0..3 {

        let err = erf(y) - x;

        let deriv = two_over_sqrt_pi * (-y * y).exp();

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

        let ans1 = 57568490574.0
            + y * (-13362590354.0
                + y * (651619640.7 + y * (-11214424.18 + y * (77392.33017 + y * (-184.9052456)))));

        let ans2 = 57568490411.0
            + y * (1029532985.0 + y * (9494680.718 + y * (59272.64853 + y * (267.8532712 + y))));

        ans1 / ans2
    } else {

        let z = 8.0 / ax;

        let y = z * z;

        let xx = ax - 0.785398163397448;

        let ans1 = 1.0
            + y * (-0.1098628627e-2
                + y * (0.2734510407e-4 + y * (-0.2073370639e-5 + y * 0.2093887211e-6)));

        let ans2 = -0.1562499995e-1
            + y * (0.1430488765e-3
                + y * (-0.6911147651e-5 + y * (0.7621095161e-6 - y * 0.934945152e-7)));

        (0.636619772 / ax).sqrt() * (xx.cos() * ans1 - z * xx.sin() * ans2)
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
            * (72362614232.0
                + y * (-7895059235.0
                    + y * (242396853.1
                        + y * (-2972611.439 + y * (15704.48260 + y * (-30.16036606))))));

        let ans2 = 144725228442.0
            + y * (2300535178.0 + y * (18583304.74 + y * (99447.43394 + y * (376.9991397 + y))));

        ans1 / ans2
    } else {

        let z = 8.0 / ax;

        let y = z * z;

        let xx = ax - 2.356194491;

        let ans1 = 1.0
            + y * (0.183105e-2
                + y * (-0.3516396496e-4 + y * (0.2457520174e-5 + y * (-0.240337019e-6))));

        let ans2 = 0.04687499995
            + y * (-0.2002690873e-3
                + y * (0.8449199096e-5 + y * (-0.88228987e-6 + y * 0.105787412e-6)));

        let ans = (0.636619772 / ax).sqrt() * (xx.cos() * ans1 - z * xx.sin() * ans2);

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

        let ans1 = -2957821389.0
            + y * (7062834065.0
                + y * (-512359803.6 + y * (10879881.29 + y * (-86327.92757 + y * 228.4622733))));

        let ans2 = 40076544269.0
            + y * (745249964.8 + y * (7189466.438 + y * (47447.26470 + y * (226.1030244 + y))));

        ans1 / ans2 + 0.636619772 * bessel_j0(x) * x.ln()
    } else {

        let z = 8.0 / x;

        let y = z * z;

        let xx = x - 0.785398163397448;

        let ans1 = 1.0
            + y * (-0.1098628627e-2
                + y * (0.2734510407e-4 + y * (-0.2073370639e-5 + y * 0.2093887211e-6)));

        let ans2 = -0.1562499995e-1
            + y * (0.1430488765e-3
                + y * (-0.6911147651e-5 + y * (0.7621095161e-6 + y * (-0.934945152e-7))));

        (0.636619772 / x).sqrt() * (xx.sin() * ans1 + z * xx.cos() * ans2)
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
            * (-0.4900604943e13
                + y * (0.1275274390e13
                    + y * (-0.5153438139e11
                        + y * (0.7349264551e9 + y * (-0.4237922726e7 + y * 0.8511937935e4)))));

        let ans2 = 0.2499580570e14
            + y * (0.4244419664e12
                + y * (0.3733650367e10
                    + y * (0.2245904002e8 + y * (0.1020426050e6 + y * (0.3549632885e3 + y)))));

        ans1 / ans2 + 0.636619772 * (bessel_j1(x) * x.ln() - 1.0 / x)
    } else {

        let z = 8.0 / x;

        let y = z * z;

        let xx = x - 2.356194491;

        let ans1 = 1.0
            + y * (0.183105e-2
                + y * (-0.3516396496e-4 + y * (0.2457520174e-5 + y * (-0.240337019e-6))));

        let ans2 = 0.04687499995
            + y * (-0.2002690873e-3
                + y * (0.8449199096e-5 + y * (-0.88228987e-6 + y * 0.105787412e-6)));

        (0.636619772 / x).sqrt() * (xx.sin() * ans1 + z * xx.cos() * ans2)
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
            * (3.5156229
                + y * (3.0899424
                    + y * (1.2067492 + y * (0.2659732 + y * (0.0360768 + y * 0.0045813)))))
    } else {

        let y = 3.75 / ax;

        (ax.exp() / ax.sqrt())
            * (0.39894228
                + y * (0.01328592
                    + y * (0.00225319
                        + y * (-0.00157565
                            + y * (0.00916281
                                + y * (-0.02057706
                                    + y * (0.02635537 + y * (-0.01647633 + y * 0.00392377))))))))
    }
}

/// Computes the modified Bessel function of the first kind, I₁(x).
#[must_use]

pub fn bessel_i1(x: f64) -> f64 {

    let ax = x.abs();

    let ans = if ax < 3.75 {

        let y = (x / 3.75).powi(2);

        ax * (0.5
            + y * (0.87890594
                + y * (0.51498869
                    + y * (0.15084934 + y * (0.02658733 + y * (0.00301532 + y * 0.00032411))))))
    } else {

        let y = 3.75 / ax;

        0.02282967 + y * (-0.02895312 + y * (0.01787654 - y * 0.00420059));

        let ans = 0.39894228
            + y * (-0.03988024
                + y * (-0.00362018
                    + y * (0.00163801
                        + y * (-0.01031555
                            + y * (0.02282967
                                + y * (-0.02895312 + y * (0.01787654 - y * 0.00420059)))))));

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

/// Computes the Legendre polynomial P_n(x) using recurrence relation.
#[must_use]

pub fn legendre_p(n: u32, x: f64) -> f64 {

    if n == 0 {

        return 1.0;
    }

    if n == 1 {

        return x;
    }

    let mut p_prev = 1.0;

    let mut p_curr = x;

    for k in 2..=n {

        let p_next = ((2 * k - 1) as f64 * x * p_curr - (k - 1) as f64 * p_prev) / k as f64;

        p_prev = p_curr;

        p_curr = p_next;
    }

    p_curr
}

/// Computes the Chebyshev polynomial of the first kind T_n(x).
#[must_use]

pub fn chebyshev_t(n: u32, x: f64) -> f64 {

    if n == 0 {

        return 1.0;
    }

    if n == 1 {

        return x;
    }

    let mut t_prev = 1.0;

    let mut t_curr = x;

    for _ in 2..=n {

        let t_next = 2.0 * x * t_curr - t_prev;

        t_prev = t_curr;

        t_curr = t_next;
    }

    t_curr
}

/// Computes the Chebyshev polynomial of the second kind U_n(x).
#[must_use]

pub fn chebyshev_u(n: u32, x: f64) -> f64 {

    if n == 0 {

        return 1.0;
    }

    if n == 1 {

        return 2.0 * x;
    }

    let mut u_prev = 1.0;

    let mut u_curr = 2.0 * x;

    for _ in 2..=n {

        let u_next = 2.0 * x * u_curr - u_prev;

        u_prev = u_curr;

        u_curr = u_next;
    }

    u_curr
}

/// Computes the (physicists') Hermite polynomial H_n(x).
#[must_use]

pub fn hermite_h(n: u32, x: f64) -> f64 {

    if n == 0 {

        return 1.0;
    }

    if n == 1 {

        return 2.0 * x;
    }

    let mut h_prev = 1.0;

    let mut h_curr = 2.0 * x;

    for k in 2..=n {

        let h_next = 2.0 * x * h_curr - 2.0 * (k - 1) as f64 * h_prev;

        h_prev = h_curr;

        h_curr = h_next;
    }

    h_curr
}

/// Computes the Laguerre polynomial L_n(x).
#[must_use]

pub fn laguerre_l(n: u32, x: f64) -> f64 {

    if n == 0 {

        return 1.0;
    }

    if n == 1 {

        return 1.0 - x;
    }

    let mut l_prev = 1.0;

    let mut l_curr = 1.0 - x;

    for k in 2..=n {

        let l_next =
            ((2 * k - 1) as f64 - x) * l_curr / k as f64 - (k - 1) as f64 * l_prev / k as f64;

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

pub fn binomial(n: u64, k: u64) -> f64 {

    if k > n {

        return 0.0;
    }

    if k == 0 || k == n {

        return 1.0;
    }

    // Use gamma for large values to avoid overflow
    gamma((n + 1) as f64) / (gamma((k + 1) as f64) * gamma((n - k + 1) as f64))
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

    for n in 1..10000 {

        let term = 1.0 / (n as f64).powf(s);

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

        (1.0 + x.exp()).ln()
    }
}
