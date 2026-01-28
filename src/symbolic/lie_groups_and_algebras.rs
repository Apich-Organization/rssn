//! # Lie Groups and Lie Algebras
//!
//! This module provides structures and functions for working with Lie groups and Lie algebras.
//! It includes representations for Lie algebra elements and Lie algebras themselves,
//! along with fundamental operations such as the Lie bracket, the exponential map,
//! and adjoint representations. Specific examples like so(3) and su(2) are also provided.

use std::sync::Arc;

use num_bigint::BigInt;
use num_traits::One;
use serde::Deserialize;
use serde::Serialize;

use crate::symbolic::core::Expr;
use crate::symbolic::matrix;

/// Represents an element of a Lie algebra, which is typically a matrix.
#[derive(
    Debug,
    Clone,
    PartialEq,
    Eq,
    Serialize,
    Deserialize,
)]

pub struct LieAlgebraElement(pub Expr);

/// Represents a Lie algebra, defined by its name and basis elements.
#[derive(
    Debug, Clone, Serialize, Deserialize,
)]

pub struct LieAlgebra {
    /// The name of the Lie algebra (e.g., "so(3)", "su(2)").
    pub name: String,
    /// The set of basis elements (generators) for the Lie algebra.
    pub basis: Vec<LieAlgebraElement>,
    /// The dimension of the Lie algebra (number of basis elements).
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
/// A `Result` containing an `Expr` representing the Lie bracket.
///
/// # Errors
///
/// This function will return an error if the operands `x` or `y` are not valid matrices
/// (i.e., `matrix::mul_matrices` returns a non-matrix `Expr`).

pub fn lie_bracket(
    x: &Expr,
    y: &Expr,
) -> Result<Expr, String> {

    let xy = matrix::mul_matrices(x, y);

    let yx = matrix::mul_matrices(y, x);

    if !matches!(xy, Expr::Matrix(_))
        || !matches!(
            yx,
            Expr::Matrix(_)
        )
    {

        return Err("Operands for \
                    lie_bracket \
                    must be valid \
                    matrices."
            .to_string());
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
/// A `Result` containing an `Expr` representing the Lie group element.
///
/// # Errors
///
/// This function will return an error if the input `x` is not a valid matrix,
/// or if it is not a square matrix.

pub fn exponential_map(
    x: &Expr,
    order: usize,
) -> Result<Expr, String> {

    let (rows, cols) =
        matrix::get_matrix_dims(x)
            .ok_or_else(|| {

                "Input must be a valid \
                 matrix."
                    .to_string()
            })?;

    if rows != cols {

        return Err("Matrix must be \
                    square for \
                    exponential map."
            .to_string());
    }

    let n = rows;

    let mut result =
        matrix::identity_matrix(n);

    let mut x_power = x.clone();

    let mut factorial = BigInt::one();

    for i in 1 ..= order {

        factorial *= i;

        let factor = Expr::new_div(
            Expr::new_bigint(BigInt::one()),
            Expr::new_bigint(
                factorial.clone(),
            ),
        );

        let term =
            matrix::scalar_mul_matrix(
                &factor,
                &x_power,
            );

        result = matrix::add_matrices(
            &result,
            &term,
        );

        x_power = matrix::mul_matrices(
            &x_power,
            x,
        );
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
/// A `Result` containing an `Expr` representing `Ad_g(X)`.
///
/// # Errors
///
/// This function will return an error if the Lie group element `g` is not invertible.

pub fn adjoint_representation_group(
    g: &Expr,
    x: &Expr,
) -> Result<Expr, String> {

    let g_inv =
        matrix::inverse_matrix(g);

    if let Expr::Variable(s) = &g_inv {

        if s.starts_with("Error:") {

            return Err(format!(
                "Failed to invert \
                 group element g: {s}"
            ));
        }
    }

    let gx = matrix::mul_matrices(g, x);

    let gxg_inv = matrix::mul_matrices(
        &gx,
        &g_inv,
    );

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
///
/// # Errors
///
/// This function will return an error if `lie_bracket` fails, which occurs
/// if `x` or `y` are not valid matrices.

pub fn adjoint_representation_algebra(
    x: &Expr,
    y: &Expr,
) -> Result<Expr, String> {

    lie_bracket(x, y)
}

/// Computes the commutator table for a Lie algebra.
///
/// # Arguments
/// * `algebra` - The Lie algebra.
///
/// # Returns
/// A `Result` containing a vector of vectors (matrix) where the element at `[i][j]`
/// is the Lie bracket `[basis[i], basis[j]]`.
///
/// # Errors
///
/// This function will return an error if `lie_bracket` fails for any pair of basis elements,
/// which occurs if the basis elements are not valid matrices.

pub fn commutator_table(
    algebra: &LieAlgebra
) -> Result<Vec<Vec<Expr>>, String> {

    let n = algebra.dimension;

    let mut table =
        Vec::with_capacity(n);

    for i in 0 .. n {

        let mut row =
            Vec::with_capacity(n);

        for j in 0 .. n {

            let bracket = lie_bracket(
                &algebra.basis[i].0,
                &algebra.basis[j].0,
            )?;

            row.push(bracket);
        }

        table.push(row);
    }

    Ok(table)
}

/// Checks if the basis of a Lie algebra satisfies the Jacobi identity.
///
/// The Jacobi identity states that `[x, [y, z]] + [y, [z, x]] + [z, [x, y]] = 0`
/// for all `x, y, z` in the algebra.
///
/// # Arguments
/// * `algebra` - The Lie algebra.
///
/// # Returns
/// `Ok(true)` if the Jacobi identity holds for all basis elements (symbolically simplifies to zero matrix),
/// `Ok(false)` otherwise.
///
/// # Errors
///
/// This function will return an error if `lie_bracket` fails during the computation of any
/// nested Lie bracket (e.g., due to invalid matrix operands).
///
/// Note: This check relies on symbolic simplification. If simplification is incomplete,
/// it might return false negatives (but shouldn't return false positives if simplification is correct).

pub fn check_jacobi_identity(
    algebra: &LieAlgebra
) -> Result<bool, String> {

    let n = algebra.dimension;

    let basis = &algebra.basis;

    for i in 0 .. n {

        for j in 0 .. n {

            for k in 0 .. n {

                let x = &basis[i].0;

                let y = &basis[j].0;

                let z = &basis[k].0;

                // [y, z]
                let yz =
                    lie_bracket(y, z)?;

                // [x, [y, z]]
                let term1 =
                    lie_bracket(
                        x, &yz,
                    )?;

                // [z, x]
                let zx =
                    lie_bracket(z, x)?;

                // [y, [z, x]]
                let term2 =
                    lie_bracket(
                        y, &zx,
                    )?;

                // [x, y]
                let xy =
                    lie_bracket(x, y)?;

                // [z, [x, y]]
                let term3 =
                    lie_bracket(
                        z, &xy,
                    )?;

                // Sum = term1 + term2 + term3
                let sum = matrix::add_matrices(
                    &matrix::add_matrices(&term1, &term2),
                    &term3,
                );

                // Check if sum is zero matrix
                if !matrix::is_zero_matrix(&sum) {

                    // Try to simplify?
                    // For now, we rely on matrix operations doing some simplification.
                    // If it's not explicitly zero, we might want to return false or try harder.
                    // Let's assume is_zero_matrix handles basic zero checks.
                    return Ok(false);
                }
            }
        }
    }

    Ok(true)
}

/// Returns the basis generators for the `so(3)` Lie algebra (infinitesimal rotations).
///
/// These are the angular momentum operators in 3D quantum mechanics (up to a factor).
/// The basis consists of three matrices representing rotations around the x, y, and z axes.
///
/// # Returns
/// A `Vec<LieAlgebraElement>` containing the three `so(3)` generators.
#[must_use]

pub fn so3_generators()
-> Vec<LieAlgebraElement> {

    let lz = Expr::Matrix(vec![
        vec![
            Expr::new_constant(0.0),
            Expr::new_constant(-1.0),
            Expr::new_constant(0.0),
        ],
        vec![
            Expr::new_constant(1.0),
            Expr::new_constant(0.0),
            Expr::new_constant(0.0),
        ],
        vec![
            Expr::new_constant(0.0),
            Expr::new_constant(0.0),
            Expr::new_constant(0.0),
        ],
    ]);

    let ly = Expr::Matrix(vec![
        vec![
            Expr::new_constant(0.0),
            Expr::new_constant(0.0),
            Expr::new_constant(1.0),
        ],
        vec![
            Expr::new_constant(0.0),
            Expr::new_constant(0.0),
            Expr::new_constant(0.0),
        ],
        vec![
            Expr::new_constant(-1.0),
            Expr::new_constant(0.0),
            Expr::new_constant(0.0),
        ],
    ]);

    let lx = Expr::Matrix(vec![
        vec![
            Expr::new_constant(0.0),
            Expr::new_constant(0.0),
            Expr::new_constant(0.0),
        ],
        vec![
            Expr::new_constant(0.0),
            Expr::new_constant(0.0),
            Expr::new_constant(-1.0),
        ],
        vec![
            Expr::new_constant(0.0),
            Expr::new_constant(1.0),
            Expr::new_constant(0.0),
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
#[must_use]

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
#[must_use]

pub fn su2_generators()
-> Vec<LieAlgebraElement> {

    let i =
        Expr::Variable("i".to_string());

    let half = Expr::new_div(
        Expr::new_bigint(One::one()),
        Expr::new_bigint(BigInt::from(2)),
    );

    let i_half =
        Expr::new_mul(i.clone(), half);

    let sx = matrix::scalar_mul_matrix(
        &i_half,
        &Expr::Matrix(vec![
            vec![
                Expr::new_constant(0.0),
                Expr::new_constant(1.0),
            ],
            vec![
                Expr::new_constant(1.0),
                Expr::new_constant(0.0),
            ],
        ]),
    );

    let sy = matrix::scalar_mul_matrix(
        &i_half,
        &Expr::Matrix(vec![
            vec![
                Expr::new_constant(0.0),
                Expr::Mul(
                    Arc::new(
                        Expr::new_constant(
                            -1.0,
                        ),
                    ),
                    Arc::new(i.clone()),
                ),
            ],
            vec![
                i,
                Expr::new_constant(0.0),
            ],
        ]),
    );

    let sz = matrix::scalar_mul_matrix(
        &i_half,
        &Expr::Matrix(vec![
            vec![
                Expr::new_constant(1.0),
                Expr::new_constant(0.0),
            ],
            vec![
                Expr::new_constant(0.0),
                Expr::new_constant(-1.0),
            ],
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
#[must_use]

pub fn su2() -> LieAlgebra {

    let basis = su2_generators();

    LieAlgebra {
        name: "su(2)".to_string(),
        dimension: basis.len(),
        basis,
    }
}
