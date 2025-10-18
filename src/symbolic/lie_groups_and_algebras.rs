//! # Lie Groups and Lie Algebras
//!
//! This module provides structures and functions for working with Lie groups and Lie algebras.
//! It includes representations for Lie algebra elements and Lie algebras themselves,
//! along with fundamental operations such as the Lie bracket, the exponential map,
//! and adjoint representations. Specific examples like so(3) and su(2) are also provided.
use crate::symbolic::core::Expr;
use crate::symbolic::matrix;
use num_bigint::BigInt;
use num_traits::One;
use std::sync::Arc;
/// Represents an element of a Lie algebra, which is typically a matrix.
#[derive(Debug, Clone, PartialEq)]
pub struct LieAlgebraElement(pub Expr);
/// Represents a Lie algebra, defined by its name and basis elements.
#[derive(Debug, Clone)]
pub struct LieAlgebra {
    pub name: String,
    pub basis: Vec<LieAlgebraElement>,
    pub dimension: usize,
}
/// Computes the Lie bracket `[X, Y] = XY - YX` for matrix Lie algebras.
///
/// The Lie bracket is the fundamental operation defining a Lie algebra.
/// For matrix Lie algebras, it is defined as the commutator of the matrices.
///
/// # Arguments
/// * `x` - The first Lie algebra element (matrix).
/// * `y` - The second Lie algebra element (matrix).
///
/// # Returns
/// A `Result` containing an `Expr` representing the Lie bracket,
/// or an error string if operands are not valid matrices.
pub fn lie_bracket(x: &Expr, y: &Expr) -> Result<Expr, String> {
    let xy = matrix::mul_matrices(x, y);
    let yx = matrix::mul_matrices(y, x);
    if !matches!(xy, Expr::Matrix(_)) || !matches!(yx, Expr::Matrix(_)) {
        return Err("Operands for lie_bracket must be valid matrices.".to_string());
    }
    Ok(matrix::sub_matrices(&xy, &yx))
}
/// Computes the exponential map `e^X` for a Lie algebra element `X` using a Taylor series expansion.
///
/// This connects the Lie algebra to the corresponding Lie group. The exponential map
/// is defined by the matrix exponential series: `e^X = I + X + X^2/2! + X^3/3! + ...`.
///
/// # Arguments
/// * `x` - The Lie algebra element (matrix).
/// * `order` - The number of terms to include in the series expansion.
///
/// # Returns
/// A `Result` containing an `Expr` representing the Lie group element,
/// or an error string if the input is not a square matrix.
pub fn exponential_map(x: &Expr, order: usize) -> Result<Expr, String> {
    let (rows, cols) =
        matrix::get_matrix_dims(x).ok_or_else(|| "Input must be a valid matrix.".to_string())?;
    if rows != cols {
        return Err("Matrix must be square for exponential map.".to_string());
    }
    let n = rows;
    let mut result = matrix::identity_matrix(n);
    let mut x_power = x.clone();
    let mut factorial = BigInt::one();
    for i in 1..=order {
        factorial *= i;
        let factor = Expr::new_div(Expr::BigInt(BigInt::one()), Expr::BigInt(factorial.clone()));
        let term = matrix::scalar_mul_matrix(&factor, &x_power);
        result = matrix::add_matrices(&result, &term);
        x_power = matrix::mul_matrices(&x_power, x);
    }
    Ok(result)
}
/// Computes the adjoint representation of a Lie group element `g` on a Lie algebra element `X`.
///
/// The adjoint representation `Ad_g(X)` is defined as `g * X * g^-1`.
/// It describes how Lie algebra elements transform under conjugation by Lie group elements.
///
/// # Arguments
/// * `g` - The Lie group element (matrix).
/// * `x` - The Lie algebra element (matrix).
///
/// # Returns
/// A `Result` containing an `Expr` representing `Ad_g(X)`,
/// or an error string if `g` is not invertible.
pub fn adjoint_representation_group(g: &Expr, x: &Expr) -> Result<Expr, String> {
    let g_inv = matrix::inverse_matrix(g);
    if let Expr::Variable(s) = &g_inv {
        if s.starts_with("Error:") {
            return Err(format!("Failed to invert group element g: {}", s));
        }
    }
    let gx = matrix::mul_matrices(g, x);
    let gxg_inv = matrix::mul_matrices(&gx, &g_inv);
    Ok(gxg_inv)
}
/// Computes the adjoint representation of a Lie algebra element `X` on another Lie algebra element `Y`.
///
/// The adjoint representation `ad_X(Y)` is defined as the Lie bracket `[X, Y]`.
///
/// # Arguments
/// * `x` - The first Lie algebra element (matrix).
/// * `y` - The second Lie algebra element (matrix).
///
/// # Returns
/// A `Result` containing an `Expr` representing `ad_X(Y)`.
pub fn adjoint_representation_algebra(x: &Expr, y: &Expr) -> Result<Expr, String> {
    lie_bracket(x, y)
}
/// Returns the basis generators for the `so(3)` Lie algebra (infinitesimal rotations).
///
/// These are the angular momentum operators in 3D quantum mechanics (up to a factor).
/// The basis consists of three matrices representing rotations around the x, y, and z axes.
///
/// # Returns
/// A `Vec<LieAlgebraElement>` containing the three `so(3)` generators.
pub fn so3_generators() -> Vec<LieAlgebraElement> {
    let lz = Expr::Matrix(vec![
        vec![
            Expr::Constant(0.0),
            Expr::Constant(-1.0),
            Expr::Constant(0.0),
        ],
        vec![
            Expr::Constant(1.0),
            Expr::Constant(0.0),
            Expr::Constant(0.0),
        ],
        vec![
            Expr::Constant(0.0),
            Expr::Constant(0.0),
            Expr::Constant(0.0),
        ],
    ]);
    let ly = Expr::Matrix(vec![
        vec![
            Expr::Constant(0.0),
            Expr::Constant(0.0),
            Expr::Constant(1.0),
        ],
        vec![
            Expr::Constant(0.0),
            Expr::Constant(0.0),
            Expr::Constant(0.0),
        ],
        vec![
            Expr::Constant(-1.0),
            Expr::Constant(0.0),
            Expr::Constant(0.0),
        ],
    ]);
    let lx = Expr::Matrix(vec![
        vec![
            Expr::Constant(0.0),
            Expr::Constant(0.0),
            Expr::Constant(0.0),
        ],
        vec![
            Expr::Constant(0.0),
            Expr::Constant(0.0),
            Expr::Constant(-1.0),
        ],
        vec![
            Expr::Constant(0.0),
            Expr::Constant(1.0),
            Expr::Constant(0.0),
        ],
    ]);
    vec![
        LieAlgebraElement(lx),
        LieAlgebraElement(ly),
        LieAlgebraElement(lz),
    ]
}
/// Creates the `so(3)` Lie algebra.
///
/// # Returns
/// A `LieAlgebra` struct representing `so(3)`.
pub fn so3() -> LieAlgebra {
    let basis = so3_generators();
    LieAlgebra {
        name: "so(3)".to_string(),
        dimension: basis.len(),
        basis,
    }
}
/// Returns the basis generators for the `su(2)` Lie algebra.
///
/// These are proportional to the Pauli matrices and describe spin-1/2 systems.
/// The basis is `{i*sigma_x/2, i*sigma_y/2, i*sigma_z/2}`.
///
/// # Returns
/// A `Vec<LieAlgebraElement>` containing the three `su(2)` generators.
pub fn su2_generators() -> Vec<LieAlgebraElement> {
    let i = Expr::Variable("i".to_string());
    let half = Expr::new_div(Expr::BigInt(One::one()), Expr::BigInt(BigInt::from(2)));
    let i_half = Expr::new_mul(i.clone(), half.clone());
    let sx = matrix::scalar_mul_matrix(
        &i_half,
        &Expr::Matrix(vec![
            vec![Expr::Constant(0.0), Expr::Constant(1.0)],
            vec![Expr::Constant(1.0), Expr::Constant(0.0)],
        ]),
    );
    let sy = matrix::scalar_mul_matrix(
        &i_half,
        &Expr::Matrix(vec![
            vec![
                Expr::Constant(0.0),
                Expr::Mul(Arc::new(Expr::Constant(-1.0)), Arc::new(i.clone())),
            ],
            vec![i.clone(), Expr::Constant(0.0)],
        ]),
    );
    let sz = matrix::scalar_mul_matrix(
        &i_half,
        &Expr::Matrix(vec![
            vec![Expr::Constant(1.0), Expr::Constant(0.0)],
            vec![Expr::Constant(0.0), Expr::Constant(-1.0)],
        ]),
    );
    vec![
        LieAlgebraElement(sx),
        LieAlgebraElement(sy),
        LieAlgebraElement(sz),
    ]
}
/// Creates the `su(2)` Lie algebra.
///
/// # Returns
/// A `LieAlgebra` struct representing `su(2)`.
pub fn su2() -> LieAlgebra {
    let basis = su2_generators();
    LieAlgebra {
        name: "su(2)".to_string(),
        dimension: basis.len(),
        basis,
    }
}
