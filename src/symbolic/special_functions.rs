//! # Symbolic Special Functions
//!
//! This module provides symbolic representations and "smart" constructors for various
//! special functions, such as Gamma, Beta, Error Function (erf), Riemann Zeta, Bessel,
//! Legendre, Laguerre, and Hermite polynomials. These constructors include simplification
//! logic for common arguments and identities.

use crate::symbolic::calculus::differentiate;
use crate::symbolic::calculus::factorial;
use crate::symbolic::core::Expr;
use crate::symbolic::simplify::{is_zero, simplify};
//use num_traits::ToPrimitive;

// =====================================================================================
// region: "Smart" Constructors with Simplification Logic
// =====================================================================================

pub fn gamma(arg: Expr) -> Expr {
    /// Symbolic representation and smart constructor for the Gamma function, `Γ(z)`.
    ///
    /// The Gamma function is an extension of the factorial function to complex numbers.
    /// This constructor includes simplification rules for integer and half-integer arguments,
    /// as well as the recurrence relation `Γ(z+1) = zΓ(z)`.
    ///
    /// # Arguments
    /// * `arg` - The argument `z` of the Gamma function.
    ///
    /// # Returns
    /// An `Expr` representing `Γ(z)`.
    let s_arg = simplify(arg);
    if let Some(n) = s_arg.to_f64() {
        if n > 0.0 && n.fract() == 0.0 {
            return Expr::Constant(factorial((n - 1.0) as usize));
        }
        if (n - 0.5).abs() < 1e-9 {
            return Expr::Sqrt(Box::new(Expr::Pi));
        }
    }
    if let Expr::Constant(1.0) = s_arg {
        return Expr::Constant(1.0);
    }
    if let Expr::Add(a, b) = &s_arg {
        if let Expr::Constant(1.0) = **b {
            return simplify(Expr::Mul(a.clone(), Box::new(gamma(*a.clone()))));
        }
        if let Expr::Constant(1.0) = **a {
            return simplify(Expr::Mul(b.clone(), Box::new(gamma(*b.clone()))));
        }
    }
    Expr::Gamma(Box::new(s_arg))
}

/// Symbolic representation and smart constructor for the Beta function, `B(x, y)`.
///
/// The Beta function is closely related to the Gamma function by the identity:
/// `B(x, y) = Γ(x)Γ(y) / Γ(x+y)`.
///
/// # Arguments
/// * `a` - The first argument `x`.
/// * `b` - The second argument `y`.
///
/// # Returns
/// An `Expr` representing `B(x, y)`.
pub fn beta(a: Expr, b: Expr) -> Expr {
    let gamma_a = gamma(a.clone());
    let gamma_b = gamma(b.clone());
    let gamma_a_plus_b = gamma(simplify(Expr::Add(Box::new(a), Box::new(b))));
    simplify(Expr::Div(
        Box::new(Expr::Mul(Box::new(gamma_a), Box::new(gamma_b))),
        Box::new(gamma_a_plus_b),
    ))
}

/// Symbolic representation and smart constructor for the Error Function, `erf(z)`.
///
/// The error function is a special function of sigmoid shape that arises in probability,
/// statistics, and partial differential equations. This constructor includes simplification
/// rules for `erf(0)` and `erf(-z)`.
///
/// # Arguments
/// * `arg` - The argument `z` of the error function.
///
/// # Returns
/// An `Expr` representing `erf(z)`.
pub fn erf(arg: Expr) -> Expr {
    let s_arg = simplify(arg);
    if is_zero(&s_arg) {
        return Expr::Constant(0.0);
    }
    if let Expr::Neg(inner) = s_arg {
        return Expr::Neg(Box::new(erf(*inner)));
    }
    Expr::Erf(Box::new(s_arg))
}

/// Symbolic representation and smart constructor for the Complementary Error Function, `erfc(z)`.
///
/// Defined as `erfc(z) = 1 - erf(z)`.
///
/// # Arguments
/// * `arg` - The argument `z`.
///
/// # Returns
/// An `Expr` representing `erfc(z)`.
pub fn erfc(arg: Expr) -> Expr {
    simplify(Expr::Sub(Box::new(Expr::Constant(1.0)), Box::new(erf(arg))))
}

/// Symbolic representation and smart constructor for the Imaginary Error Function, `erfi(z)`.
///
/// Defined as `erfi(z) = -i * erf(iz)`.
///
/// # Arguments
/// * `arg` - The argument `z`.
///
/// # Returns
/// An `Expr` representing `erfi(z)`.
pub fn erfi(arg: Expr) -> Expr {
    let i = Expr::Complex(Box::new(Expr::Constant(0.0)), Box::new(Expr::Constant(1.0)));
    simplify(Expr::Mul(
        Box::new(Expr::Neg(Box::new(i.clone()))),
        Box::new(erf(Expr::Mul(Box::new(i), Box::new(arg)))),
    ))
}

/// Symbolic representation and smart constructor for the Riemann Zeta function, `ζ(s)`.
///
/// The Riemann Zeta function is a central object in number theory, with connections
/// to the distribution of prime numbers. This constructor includes simplification
/// rules for specific integer arguments (e.g., `ζ(0)`, `ζ(1)`, `ζ(2)`) and trivial zeros.
///
/// # Arguments
/// * `arg` - The argument `s` of the Zeta function.
///
/// # Returns
/// An `Expr` representing `ζ(s)`.
pub fn zeta(arg: Expr) -> Expr {
    let s_arg = simplify(arg);
    if let Some(n) = s_arg.to_f64() {
        if n.fract() == 0.0 {
            let n_int = n as i32;
            if n_int == 0 {
                return Expr::Constant(-0.5);
            }
            if n_int == 1 {
                return Expr::Infinity;
            } // Pole
            if n_int == 2 {
                return simplify(Expr::Div(
                    Box::new(Expr::Power(
                        Box::new(Expr::Pi),
                        Box::new(Expr::Constant(2.0)),
                    )),
                    Box::new(Expr::Constant(6.0)),
                ));
            }
            if n_int < 0 && n_int % 2 == 0 {
                return Expr::Constant(0.0);
            } // Trivial zeros
        }
    }
    Expr::Zeta(Box::new(s_arg))
}

/// Symbolic representation and smart constructor for the Digamma function, `ψ(z)`.
///
/// The Digamma function is the logarithmic derivative of the Gamma function: `ψ(z) = Γ'(z)/Γ(z)`.
/// This constructor includes simplification rules for `ψ(1)` (related to Euler-Mascheroni constant)
/// and the recurrence relation `ψ(z+1) = ψ(z) + 1/z`.
///
/// # Arguments
/// * `arg` - The argument `z` of the Digamma function.
///
/// # Returns
/// An `Expr` representing `ψ(z)`.
pub fn digamma(arg: Expr) -> Expr {
    let s_arg = simplify(arg);
    if let Some(n) = s_arg.to_f64() {
        if (n - 1.0).abs() < 1e-9 {
            // Euler-Mascheroni constant, approx -0.57721
            return Expr::Variable("-gamma".to_string());
        }
    }
    if let Expr::Add(a, b) = &s_arg {
        if let Expr::Constant(1.0) = **b {
            return simplify(Expr::Add(
                Box::new(digamma(*a.clone())),
                Box::new(Expr::Div(Box::new(Expr::Constant(1.0)), a.clone())),
            ));
        }
        if let Expr::Constant(1.0) = **a {
            return simplify(Expr::Add(
                Box::new(digamma(*b.clone())),
                Box::new(Expr::Div(Box::new(Expr::Constant(1.0)), b.clone())),
            ));
        }
    }
    Expr::Digamma(Box::new(s_arg))
}

/// Symbolic representation and smart constructor for the Bessel function of the first kind, `J_n(x)`.
///
/// Bessel functions are solutions to Bessel's differential equation and arise in many
/// problems of wave propagation and potential theory. This constructor includes simplification
/// rules for `J_n(0)` and `J_{-n}(x)`.
///
/// # Arguments
/// * `order` - The order `n` of the Bessel function.
/// * `arg` - The argument `x` of the Bessel function.
///
/// # Returns
/// An `Expr` representing `J_n(x)`.
pub fn bessel_j(order: Expr, arg: Expr) -> Expr {
    let s_order = simplify(order);
    let s_arg = simplify(arg);
    if is_zero(&s_arg) {
        if let Some(n) = s_order.to_f64() {
            if n == 0.0 {
                return Expr::Constant(1.0);
            }
            if n > 0.0 {
                return Expr::Constant(0.0);
            }
        }
    }
    if let Expr::Neg(inner_order) = &s_order {
        if let Some(n) = inner_order.to_f64() {
            if n.fract() == 0.0 {
                let factor =
                    Expr::Power(Box::new(Expr::Constant(-1.0)), Box::new(Expr::Constant(n)));
                return simplify(Expr::Mul(
                    Box::new(factor),
                    Box::new(bessel_j(*inner_order.clone(), s_arg)),
                ));
            }
        }
    }
    Expr::BesselJ(Box::new(s_order), Box::new(s_arg))
}

/// Symbolic representation and smart constructor for the Bessel function of the second kind, `Y_n(x)`.
///
/// Bessel functions of the second kind are also solutions to Bessel's differential equation,
/// linearly independent from `J_n(x)`. They typically have a singularity at the origin.
///
/// # Arguments
/// * `order` - The order `n` of the Bessel function.
/// * `arg` - The argument `x` of the Bessel function.
///
/// # Returns
/// An `Expr` representing `Y_n(x)`.
pub fn bessel_y(order: Expr, arg: Expr) -> Expr {
    let s_order = simplify(order);
    let s_arg = simplify(arg);
    if is_zero(&s_arg) {
        return Expr::NegativeInfinity;
    }
    Expr::BesselY(Box::new(s_order), Box::new(s_arg))
}

/// Symbolic representation and smart constructor for the Legendre Polynomials, `P_n(x)`.
///
/// Legendre polynomials are solutions to Legendre's differential equation and form a complete
/// orthogonal set over the interval `[-1, 1]`. They are widely used in physics and engineering.
/// This constructor includes simplification rules for `P_0(x)`, `P_1(x)`, and the recurrence relation.
///
/// # Arguments
/// * `degree` - The degree `n` of the Legendre polynomial.
/// * `arg` - The argument `x` of the polynomial.
///
/// # Returns
/// An `Expr` representing `P_n(x)`.
pub fn legendre_p(degree: Expr, arg: Expr) -> Expr {
    let s_degree = simplify(degree);
    if let Some(n) = s_degree.to_f64() {
        let n_int = n as i32;
        if n >= 0.0 && n.fract() == 0.0 {
            if n_int == 0 {
                return Expr::Constant(1.0);
            }
            if n_int == 1 {
                return arg;
            }
            let p_n = legendre_p(Expr::Constant(n - 1.0), arg.clone());
            let p_n_minus_1 = legendre_p(Expr::Constant(n - 2.0), arg.clone());
            let term1 = Expr::Mul(
                Box::new(Expr::Constant(2.0 * n - 1.0)),
                Box::new(Expr::Mul(Box::new(arg), Box::new(p_n))),
            );
            let term2 = Expr::Mul(Box::new(Expr::Constant(n - 1.0)), Box::new(p_n_minus_1));
            return simplify(Expr::Div(
                Box::new(Expr::Sub(Box::new(term1), Box::new(term2))),
                Box::new(Expr::Constant(n)),
            ));
        }
    }
    Expr::LegendreP(Box::new(s_degree), Box::new(arg))
}

/// Symbolic representation and smart constructor for the Laguerre Polynomials, `L_n(x)`.
///
/// Laguerre polynomials are solutions to Laguerre's differential equation and are important
/// in quantum mechanics (e.g., hydrogen atom wave functions) and other areas.
/// This constructor includes simplification rules for `L_0(x)`, `L_1(x)`, and the recurrence relation.
///
/// # Arguments
/// * `degree` - The degree `n` of the Laguerre polynomial.
/// * `arg` - The argument `x` of the polynomial.
///
/// # Returns
/// An `Expr` representing `L_n(x)`.
pub fn laguerre_l(degree: Expr, arg: Expr) -> Expr {
    let s_degree = simplify(degree);
    if let Some(n) = s_degree.to_f64() {
        let n_int = n as i32;
        if n >= 0.0 && n.fract() == 0.0 {
            if n_int == 0 {
                return Expr::Constant(1.0);
            }
            if n_int == 1 {
                return simplify(Expr::Sub(Box::new(Expr::Constant(1.0)), Box::new(arg)));
            }
            let l_n = laguerre_l(Expr::Constant(n - 1.0), arg.clone());
            let l_n_minus_1 = laguerre_l(Expr::Constant(n - 2.0), arg.clone());
            let term1_factor = simplify(Expr::Sub(
                Box::new(Expr::Constant(2.0 * n - 1.0)),
                Box::new(arg),
            ));
            let term1 = Expr::Mul(Box::new(term1_factor), Box::new(l_n));
            let term2 = Expr::Mul(Box::new(Expr::Constant(n - 1.0)), Box::new(l_n_minus_1));
            return simplify(Expr::Div(
                Box::new(Expr::Sub(Box::new(term1), Box::new(term2))),
                Box::new(Expr::Constant(n)),
            ));
        }
    }
    Expr::LaguerreL(Box::new(s_degree), Box::new(arg))
}

/// Symbolic representation and smart constructor for the Hermite Polynomials, `H_n(x)`.
///
/// Hermite polynomials are solutions to Hermite's differential equation and are crucial
/// in quantum mechanics (e.g., harmonic oscillator) and probability theory.
/// This constructor includes simplification rules for `H_0(x)`, `H_1(x)`, and the recurrence relation.
///
/// # Arguments
/// * `degree` - The degree `n` of the Hermite polynomial.
/// * `arg` - The argument `x` of the polynomial.
///
/// # Returns
/// An `Expr` representing `H_n(x)`.
pub fn hermite_h(degree: Expr, arg: Expr) -> Expr {
    let s_degree = simplify(degree);
    if let Some(n) = s_degree.to_f64() {
        let n_int = n as i32;
        if n >= 0.0 && n.fract() == 0.0 {
            if n_int == 0 {
                return Expr::Constant(1.0);
            }
            if n_int == 1 {
                return simplify(Expr::Mul(Box::new(Expr::Constant(2.0)), Box::new(arg)));
            }
            let h_n = hermite_h(Expr::Constant(n - 1.0), arg.clone());
            let h_n_minus_1 = hermite_h(Expr::Constant(n - 2.0), arg.clone());
            let term1 = Expr::Mul(
                Box::new(Expr::Constant(2.0)),
                Box::new(Expr::Mul(Box::new(arg), Box::new(h_n))),
            );
            let term2 = Expr::Mul(
                Box::new(Expr::Constant(2.0 * (n - 1.0))),
                Box::new(h_n_minus_1),
            );
            return simplify(Expr::Sub(Box::new(term1), Box::new(term2)));
        }
    }
    Expr::HermiteH(Box::new(s_degree), Box::new(arg))
}

// =====================================================================================
// endregion: "Smart" Constructors
// =====================================================================================

// =====================================================================================
// region: Differential Equations and Formulas (Preserved from original)
// =====================================================================================

/// Represents Bessel's differential equation: `x^2 y'' + x y' + (x^2 - n^2) y = 0`.
///
/// This second-order linear ordinary differential equation is fundamental in many areas
/// of physics and engineering, particularly in problems involving cylindrical symmetry.
/// Its solutions are known as Bessel functions.
///
/// # Arguments
/// * `y` - The unknown function `y(x)`.
/// * `x` - The independent variable `x`.
/// * `n` - The order `n` of the Bessel equation.
///
/// # Returns
/// An `Expr::Eq` representing Bessel's differential equation.
pub fn bessel_differential_equation(y: &Expr, x: &Expr, n: &Expr) -> Expr {
    let y_prime = differentiate(y, "x");
    let y_double_prime = differentiate(&y_prime, "x");
    let term1 = Expr::Mul(
        Box::new(Expr::Power(
            Box::new(x.clone()),
            Box::new(Expr::Constant(2.0)),
        )),
        Box::new(y_double_prime),
    );
    let term2 = Expr::Mul(Box::new(x.clone()), Box::new(y_prime));
    let term3 = Expr::Mul(
        Box::new(Expr::Sub(
            Box::new(Expr::Power(
                Box::new(x.clone()),
                Box::new(Expr::Constant(2.0)),
            )),
            Box::new(Expr::Power(
                Box::new(n.clone()),
                Box::new(Expr::Constant(2.0)),
            )),
        )),
        Box::new(y.clone()),
    );
    Expr::Eq(
        Box::new(Expr::Add(
            Box::new(term1),
            Box::new(Expr::Add(Box::new(term2), Box::new(term3))),
        )),
        Box::new(Expr::Constant(0.0)),
    )
}

/// Represents Legendre's differential equation: `(1-x^2)y'' - 2xy' + n(n+1)y = 0`.
///
/// This second-order linear ordinary differential equation arises in the solution of
/// Laplace's equation in spherical coordinates. Its solutions are Legendre polynomials.
///
/// # Arguments
/// * `y` - The unknown function `y(x)`.
/// * `x` - The independent variable `x`.
/// * `n` - The degree `n` of the Legendre equation.
///
/// # Returns
/// An `Expr::Eq` representing Legendre's differential equation.
pub fn legendre_differential_equation(y: &Expr, x: &Expr, n: &Expr) -> Expr {
    let y_prime = differentiate(y, "x");
    let y_double_prime = differentiate(&y_prime, "x");
    let term1 = Expr::Mul(
        Box::new(Expr::Sub(
            Box::new(Expr::Constant(1.0)),
            Box::new(Expr::Power(
                Box::new(x.clone()),
                Box::new(Expr::Constant(2.0)),
            )),
        )),
        Box::new(y_double_prime),
    );
    let term2 = Expr::Mul(
        Box::new(Expr::Constant(-2.0)),
        Box::new(Expr::Mul(Box::new(x.clone()), Box::new(y_prime))),
    );
    let term3 = Expr::Mul(
        Box::new(Expr::Mul(
            Box::new(n.clone()),
            Box::new(Expr::Add(
                Box::new(n.clone()),
                Box::new(Expr::Constant(1.0)),
            )),
        )),
        Box::new(y.clone()),
    );
    Expr::Eq(
        Box::new(Expr::Add(
            Box::new(term1),
            Box::new(Expr::Add(Box::new(term2), Box::new(term3))),
        )),
        Box::new(Expr::Constant(0.0)),
    )
}

/// Represents Rodrigues' Formula for Legendre Polynomials: `P_n(x) = (1/(2^n n!)) * d^n/dx^n [(x^2 - 1)^n]`.
///
/// This formula provides a compact way to define Legendre polynomials and is useful
/// for deriving their properties.
///
/// # Arguments
/// * `n` - The degree `n` of the Legendre polynomial.
/// * `x` - The independent variable `x`.
///
/// # Returns
/// An `Expr::Eq` representing Rodrigues' formula.
pub fn legendre_rodrigues_formula(n: &Expr, x: &Expr) -> Expr {
    let n_f64 = if let Expr::Constant(val) = n {
        *val
    } else {
        return Expr::LegendreP(Box::new(n.clone()), Box::new(x.clone()));
    };
    let n_factorial = Expr::Constant(factorial(n_f64 as usize));
    Expr::Eq(
        Box::new(legendre_p(n.clone(), x.clone())),
        Box::new(Expr::Mul(
            Box::new(Expr::Div(
                Box::new(Expr::Constant(1.0)),
                Box::new(Expr::Mul(
                    Box::new(Expr::Power(
                        Box::new(Expr::Constant(2.0)),
                        Box::new(n.clone()),
                    )),
                    Box::new(n_factorial),
                )),
            )),
            Box::new(Expr::DerivativeN(
                Box::new(Expr::Power(
                    Box::new(Expr::Sub(
                        Box::new(Expr::Power(
                            Box::new(x.clone()),
                            Box::new(Expr::Constant(2.0)),
                        )),
                        Box::new(Expr::Constant(1.0)),
                    )),
                    Box::new(n.clone()),
                )),
                "x".to_string(),
                Box::new(n.clone()),
            )),
        )),
    )
}

/// Represents Laguerre's differential equation: `xy'' + (1-x)y' + ny = 0`.
///
/// This second-order linear ordinary differential equation arises in the quantum mechanical
/// treatment of the hydrogen atom. Its solutions are Laguerre polynomials.
///
/// # Arguments
/// * `y` - The unknown function `y(x)`.
/// * `x` - The independent variable `x`.
/// * `n` - The degree `n` of the Laguerre equation.
///
/// # Returns
/// An `Expr::Eq` representing Laguerre's differential equation.
pub fn laguerre_differential_equation(y: &Expr, x: &Expr, n: &Expr) -> Expr {
    let y_prime = differentiate(y, "x");
    let y_double_prime = differentiate(&y_prime, "x");
    let term1 = Expr::Mul(Box::new(x.clone()), Box::new(y_double_prime));
    let term2 = Expr::Mul(
        Box::new(Expr::Sub(
            Box::new(Expr::Constant(1.0)),
            Box::new(x.clone()),
        )),
        Box::new(y_prime),
    );
    let term3 = Expr::Mul(Box::new(n.clone()), Box::new(y.clone()));
    Expr::Eq(
        Box::new(Expr::Add(
            Box::new(term1),
            Box::new(Expr::Add(Box::new(term2), Box::new(term3))),
        )),
        Box::new(Expr::Constant(0.0)),
    )
}

/// Represents Hermite's differential equation: `y'' - 2xy' + 2ny = 0`.
///
/// This second-order linear ordinary differential equation is important in quantum mechanics
/// (e.g., for the quantum harmonic oscillator) and probability theory. Its solutions are Hermite polynomials.
///
/// # Arguments
/// * `y` - The unknown function `y(x)`.
/// * `x` - The independent variable `x`.
/// * `n` - The degree `n` of the Hermite equation.
///
/// # Returns
/// An `Expr::Eq` representing Hermite's differential equation.
pub fn hermite_differential_equation(y: &Expr, x: &Expr, n: &Expr) -> Expr {
    let y_prime = differentiate(y, "x");
    let y_double_prime = differentiate(&y_prime, "x");
    let term1 = y_double_prime;
    let term2 = Expr::Mul(
        Box::new(Expr::Constant(-2.0)),
        Box::new(Expr::Mul(Box::new(x.clone()), Box::new(y_prime))),
    );
    let term3 = Expr::Mul(
        Box::new(Expr::Mul(
            Box::new(Expr::Constant(2.0)),
            Box::new(n.clone()),
        )),
        Box::new(y.clone()),
    );
    Expr::Eq(
        Box::new(Expr::Add(
            Box::new(term1),
            Box::new(Expr::Add(Box::new(term2), Box::new(term3))),
        )),
        Box::new(Expr::Constant(0.0)),
    )
}

/// Represents Rodrigues' Formula for Hermite Polynomials: `H_n(x) = (-1)^n e^(x^2) d^n/dx^n [e^(-x^2)]`.
///
/// This formula provides a compact way to define Hermite polynomials and is useful
/// for deriving their properties.
///
/// # Arguments
/// * `n` - The degree `n` of the Hermite polynomial.
/// * `x` - The independent variable `x`.
///
/// # Returns
/// An `Expr::Eq` representing Rodrigues' formula.
pub fn hermite_rodrigues_formula(n: &Expr, x: &Expr) -> Expr {
    Expr::Eq(
        Box::new(hermite_h(n.clone(), x.clone())),
        Box::new(Expr::Mul(
            Box::new(Expr::Power(
                Box::new(Expr::Constant(-1.0)),
                Box::new(n.clone()),
            )),
            Box::new(Expr::Mul(
                Box::new(Expr::Exp(Box::new(Expr::Power(
                    Box::new(x.clone()),
                    Box::new(Expr::Constant(2.0)),
                )))),
                Box::new(Expr::DerivativeN(
                    Box::new(Expr::Exp(Box::new(Expr::Neg(Box::new(Expr::Power(
                        Box::new(x.clone()),
                        Box::new(Expr::Constant(2.0)),
                    )))))),
                    "x".to_string(),
                    Box::new(n.clone()),
                )),
            )),
        )),
    )
}

// =====================================================================================
// endregion: Differential Equations and Formulas
// =====================================================================================
