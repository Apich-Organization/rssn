//! # Symbolic Matrix Operations
//!
//! This module provides functions for symbolic matrix operations, where matrix elements
//! are represented by `Expr` (symbolic expressions). It includes basic matrix arithmetic
//! (addition, subtraction, multiplication), transposition, determinant calculation,
//! inversion, RREF (Reduced Row Echelon Form), null space computation, and eigenvalue
//! decomposition.

use num_bigint::BigInt;
use num_traits::One;
use num_traits::Zero;

use crate::symbolic::core::Expr;
use crate::symbolic::simplify::is_zero;
use crate::symbolic::simplify_dag::simplify;
use crate::symbolic::solve::solve;

/// Helper to get dimensions of a matrix `Expr`.
///
/// # Arguments
/// * `matrix` - The matrix expression.
///
/// # Returns
/// An `Option<(usize, usize)>` containing `(rows, cols)` if the expression is a valid matrix,
/// `None` otherwise.
#[must_use]

pub fn get_matrix_dims(
    matrix: &Expr
) -> Option<(usize, usize)> {

    if let Expr::Matrix(rows) = matrix {

        if rows.is_empty() {

            return Some((0, 0));
        }

        let num_rows = rows.len();

        let num_cols = rows[0].len();

        if rows
            .iter()
            .all(|row| {

                row.len() == num_cols
            })
        {

            Some((num_rows, num_cols))
        } else {

            None
        }
    } else {

        None
    }
}

#[must_use]

/// Creates a matrix of the specified dimensions filled with symbolic zeros.
///
/// # Arguments
/// * `rows` - The number of rows.
/// * `cols` - The number of columns.
///
/// # Returns
/// A 2D vector of `Expr` representing the zero matrix.

pub fn create_empty_matrix(
    rows: usize,
    cols: usize,
) -> Vec<Vec<Expr>> {

    vec![
        vec![
            Expr::BigInt(BigInt::zero());
            cols
        ];
        rows
    ]
}

/// Creates an identity matrix of a given size.
///
/// An identity matrix is a square matrix with ones on the main diagonal
/// and zeros elsewhere. It acts as the multiplicative identity in matrix multiplication.
///
/// # Arguments
/// * `size` - The dimension of the square identity matrix.
///
/// # Returns
/// An `Expr::Matrix` representing the identity matrix.
#[must_use]

pub fn identity_matrix(
    size: usize
) -> Expr {

    let mut rows =
        create_empty_matrix(size, size);

    for (i, row) in rows.iter_mut().enumerate() {

        row[i] =
            Expr::BigInt(BigInt::one());
    }

    Expr::Matrix(rows)
}

/// Adds two matrices element-wise.
///
/// # Arguments
/// * `m1` - The first matrix as an `Expr::Matrix`.
/// * `m2` - The second matrix as an `Expr::Matrix`.
///
/// # Returns
/// An `Expr::Matrix` representing the sum, or an unevaluated `Expr::Add` if dimensions are incompatible.
#[must_use]

pub fn add_matrices(
    m1: &Expr,
    m2: &Expr,
) -> Expr {

    let dims1 = get_matrix_dims(m1);

    let dims2 = get_matrix_dims(m2);

    if let (
        Some((r1, c1)),
        Some((r2, c2)),
    ) = (dims1, dims2)
    {

        if r1 != r2 || c1 != c2 {

            return Expr::new_add(
                m1.clone(),
                m2.clone(),
            );
        }

        let Expr::Matrix(rows1) = m1
        else {

            unreachable!()
        };

        let Expr::Matrix(rows2) = m2
        else {

            unreachable!()
        };

        let mut result_rows =
            create_empty_matrix(r1, c1);

        for i in 0 .. r1 {

            for j in 0 .. c1 {

                result_rows[i][j] =
                    simplify(
                        &Expr::new_add(
                            rows1[i][j]
                                .clone(
                                ),
                            rows2[i][j]
                                .clone(
                                ),
                        ),
                    );
            }
        }

        Expr::Matrix(result_rows)
    } else {

        Expr::new_add(
            m1.clone(),
            m2.clone(),
        )
    }
}

/// Subtracts one matrix from another element-wise.
///
/// # Arguments
/// * `m1` - The first matrix as an `Expr::Matrix`.
/// * `m2` - The second matrix as an `Expr::Matrix`.
///
/// # Returns
/// An `Expr::Matrix` representing the difference, or an unevaluated `Expr::Sub` if dimensions are incompatible.
#[must_use]

pub fn sub_matrices(
    m1: &Expr,
    m2: &Expr,
) -> Expr {

    let dims1 = get_matrix_dims(m1);

    let dims2 = get_matrix_dims(m2);

    if let (
        Some((r1, c1)),
        Some((r2, c2)),
    ) = (dims1, dims2)
    {

        if r1 != r2 || c1 != c2 {

            return Expr::new_sub(
                m1.clone(),
                m2.clone(),
            );
        }

        let Expr::Matrix(rows1) = m1
        else {

            unreachable!()
        };

        let Expr::Matrix(rows2) = m2
        else {

            unreachable!()
        };

        let mut result_rows =
            create_empty_matrix(r1, c1);

        for i in 0 .. r1 {

            for j in 0 .. c1 {

                result_rows[i][j] =
                    simplify(
                        &Expr::new_sub(
                            rows1[i][j]
                                .clone(
                                ),
                            rows2[i][j]
                                .clone(
                                ),
                        ),
                    );
            }
        }

        Expr::Matrix(result_rows)
    } else {

        Expr::new_sub(
            m1.clone(),
            m2.clone(),
        )
    }
}

/// Multiplies two matrices.
///
/// # Arguments
/// * `m1` - The first matrix as an `Expr::Matrix`.
/// * `m2` - The second matrix as an `Expr::Matrix`.
///
/// # Returns
/// An `Expr::Matrix` representing the product, or an unevaluated `Expr::Mul` if dimensions are incompatible.
#[must_use]

pub fn mul_matrices(
    m1: &Expr,
    m2: &Expr,
) -> Expr {

    let dims1 = get_matrix_dims(m1);

    let dims2 = get_matrix_dims(m2);

    if let (
        Some((r1, c1)),
        Some((r2, c2)),
    ) = (dims1, dims2)
    {

        if c1 != r2 {

            return Expr::new_mul(
                m1.clone(),
                m2.clone(),
            );
        }

        let Expr::Matrix(rows1) = m1
        else {

            unreachable!()
        };

        let Expr::Matrix(rows2) = m2
        else {

            unreachable!()
        };

        let mut result_rows =
            create_empty_matrix(r1, c2);

        for i in 0 .. r1 {

            for j in 0 .. c2 {

                let mut sum_term =
                    Expr::BigInt(
                        BigInt::zero(),
                    );

                for k in 0 .. c1 {

                    sum_term = simplify(&Expr::new_add(
                        sum_term,
                        simplify(&Expr::new_mul(
                            rows1[i][k].clone(),
                            rows2[k][j].clone(),
                        )),
                    ));
                }

                result_rows[i][j] =
                    sum_term;
            }
        }

        Expr::Matrix(result_rows)
    } else {

        Expr::new_mul(
            m1.clone(),
            m2.clone(),
        )
    }
}

/// Multiplies a matrix by a scalar expression.
///
/// # Arguments
/// * `scalar` - The scalar expression.
/// * `matrix` - The matrix as an `Expr::Matrix`.
///
/// # Returns
/// An `Expr::Matrix` representing the scaled matrix, or an unevaluated `Expr::Mul` if the matrix is invalid.
#[must_use]

pub fn scalar_mul_matrix(
    scalar: &Expr,
    matrix: &Expr,
) -> Expr {

    if let Some((r, c)) =
        get_matrix_dims(matrix)
    {

        let Expr::Matrix(rows) = matrix
        else {

            unreachable!()
        };

        let mut result_rows =
            create_empty_matrix(r, c);

        for i in 0 .. r {

            for j in 0 .. c {

                result_rows[i][j] =
                    simplify(
                        &Expr::new_mul(
                            scalar
                                .clone(
                                ),
                            rows[i][j]
                                .clone(
                                ),
                        ),
                    );
            }
        }

        Expr::Matrix(result_rows)
    } else {

        Expr::new_mul(
            scalar.clone(),
            matrix.clone(),
        )
    }
}

/// Computes the transpose of a matrix.
///
/// # Arguments
/// * `matrix` - The matrix as an `Expr::Matrix`.
///
/// # Returns
/// An `Expr::Matrix` representing the transposed matrix, or an unevaluated `Expr::Power` if the matrix is invalid.
#[must_use]

pub fn transpose_matrix(
    matrix: &Expr
) -> Expr {

    if let Some((r, c)) =
        get_matrix_dims(matrix)
    {

        let Expr::Matrix(rows) = matrix
        else {

            unreachable!()
        };

        let mut result_rows =
            create_empty_matrix(c, r);

        for i in 0 .. r {

            for j in 0 .. c {

                result_rows[j][i] =
                    rows[i][j].clone();
            }
        }

        Expr::Matrix(result_rows)
    } else {

        Expr::new_pow(
            matrix.clone(),
            Expr::Variable(
                "T".to_string(),
            ),
        )
    }
}

/// Computes the determinant of a square matrix.
///
/// This function uses Laplace expansion (cofactor expansion) along the first row.
/// It is computationally expensive for large matrices.
///
/// # Arguments
/// * `matrix` - The square matrix as an `Expr::Matrix`.
///
/// # Returns
/// An `Expr` representing the determinant, or an error message if the matrix is not square.
#[must_use]

pub fn determinant(
    matrix: &Expr
) -> Expr {

    if let Some((r, c)) =
        get_matrix_dims(matrix)
    {

        if r != c {

            return Expr::Variable(
                "Error: Matrix must \
                 be square"
                    .to_string(),
            );
        }

        if r == 0 {

            return Expr::BigInt(
                BigInt::one(),
            );
        }

        if r == 1 {

            let Expr::Matrix(rows) =
                matrix
            else {

                unreachable!()
            };

            return rows[0][0].clone();
        }

        if r == 2 {

            let Expr::Matrix(rows) =
                matrix
            else {

                unreachable!()
            };

            let a = &rows[0][0];

            let b = &rows[0][1];

            let c = &rows[1][0];

            let d = &rows[1][1];

            return simplify(
                &Expr::new_sub(
                    Expr::new_mul(
                        a.clone(),
                        d.clone(),
                    ),
                    Expr::new_mul(
                        b.clone(),
                        c.clone(),
                    ),
                ),
            );
        }

        let Expr::Matrix(rows) = matrix
        else {

            unreachable!()
        };

        let mut det = Expr::BigInt(
            BigInt::zero(),
        );

        for (j, row0j) in rows[0].iter().enumerate() {

            let minor =
                get_minor(matrix, 0, j);

            let sign = if j % 2 == 0 {

                Expr::BigInt(
                    BigInt::one(),
                )
            } else {

                Expr::BigInt(
                    BigInt::from(-1),
                )
            };

            let term = simplify(
                &Expr::new_mul(
                    row0j.clone(),
                    determinant(&minor),
                ),
            );

            det = simplify(
                &Expr::new_add(
                    det,
                    Expr::new_mul(
                        sign, term,
                    ),
                ),
            );
        }

        det
    } else {

        Expr::Variable(
            "Error: Not a valid matrix"
                .to_string(),
        )
    }
}

pub(crate) fn get_minor(
    matrix: &Expr,
    row_to_remove: usize,
    col_to_remove: usize,
) -> Expr {

    if let Some((_r, _c)) =
        get_matrix_dims(matrix)
    {

        let Expr::Matrix(rows) = matrix
        else {

            unreachable!()
        };

        let mut minor_rows = Vec::new();

        for (i, row) in rows.iter().enumerate() {

            if i == row_to_remove {

                continue;
            }

            let mut new_row =
                Vec::new();

            for (j, item) in row.iter().enumerate() {

                if j == col_to_remove {

                    continue;
                }

                new_row.push(
                    item.clone(),
                );
            }

            minor_rows.push(new_row);
        }

        Expr::Matrix(minor_rows)
    } else {

        matrix.clone()
    }
}

/// Computes the inverse of a square matrix using the adjugate matrix method.
///
/// The inverse of a matrix `A` is given by `A^-1 = (1/det(A)) * adj(A)`,
/// where `adj(A)` is the adjugate matrix (transpose of the cofactor matrix).
///
/// # Arguments
/// * `matrix` - The square matrix as an `Expr::Matrix`.
///
/// # Returns
/// An `Expr::Matrix` representing the inverse, or an error message if the matrix is singular or not square.
#[must_use]

pub fn inverse_matrix(
    matrix: &Expr
) -> Expr {

    let det = determinant(matrix);

    if let Expr::Variable(_) = det {

        return det;
    }

    if is_zero(&det) {

        return Expr::Variable(
            "Error: Matrix is \
             singular and cannot be \
             inverted"
                .to_string(),
        );
    }

    if let Some((r, c)) =
        get_matrix_dims(matrix)
    {

        if r != c {

            return Expr::Variable(
                "Error: Matrix must \
                 be square to have an \
                 inverse"
                    .to_string(),
            );
        }

        let mut adj_rows =
            create_empty_matrix(r, c);

        for (i, adj_row_i) in adj_rows.iter_mut().enumerate() {

            for (j, adj_item) in adj_row_i.iter_mut().enumerate() {

                let minor = get_minor(
                    matrix,
                    j,
                    i,
                );

                let sign = if (i + j)
                    % 2
                    == 0
                {

                    Expr::BigInt(
                        BigInt::one(),
                    )
                } else {

                    Expr::BigInt(
                        BigInt::from(
                            -1,
                        ),
                    )
                };

                let cofactor = simplify(
                    &Expr::new_mul(
                        sign,
                        determinant(
                            &minor,
                        ),
                    ),
                );

                *adj_item =
                    cofactor;
            }
        }

        let adj_matrix =
            Expr::Matrix(adj_rows);

        scalar_mul_matrix(
            &simplify(&Expr::new_div(
                Expr::BigInt(
                    BigInt::one(),
                ),
                det,
            )),
            &adj_matrix,
        )
    } else {

        Expr::Variable(
            "Error: Not a valid matrix"
                .to_string(),
        )
    }
}

/// Solves a system of linear equations `Ax = b` for any `M x N` matrix `A`.
///
/// This function constructs an augmented matrix `[A | b]`, computes its Reduced Row Echelon Form (RREF),
/// and then analyzes the RREF to determine the nature of the solution:
/// - **Unique Solution**: Returns a column vector `x`.
/// - **Infinite Solutions**: Returns a parametric solution (particular solution + null space basis).
/// - **No Solution**: Returns `Expr::NoSolution`.
///
/// # Arguments
/// * `a` - An `Expr::Matrix` representing the coefficient matrix `A`.
/// * `b` - An `Expr::Matrix` representing the constant vector `b` (must be a column vector).
///
/// # Returns
/// A `Result` containing an `Expr` representing the solution (matrix, system, or no solution).
///
/// # Errors
///
/// This function will return an error if:
/// - `A` or `b` are not valid matrices.
/// - The row dimensions of `A` and `b` are incompatible.
/// - `b` is not a column vector.
/// - The `rref` computation fails.
/// - The `null_space` computation fails.

pub fn solve_linear_system(
    a: &Expr,
    b: &Expr,
) -> Result<Expr, String> {

    let (a_rows, a_cols) =
        get_matrix_dims(a).ok_or_else(
            || {

                "A is not a valid \
                 matrix"
                    .to_string()
            },
        )?;

    let (b_rows, b_cols) =
        get_matrix_dims(b).ok_or_else(
            || {

                "b is not a valid \
                 matrix"
                    .to_string()
            },
        )?;

    if a_rows != b_rows {

        return Err("Matrix A and \
                    vector b have \
                    incompatible \
                    row dimensions"
            .to_string());
    }

    if b_cols != 1 {

        return Err("b must be a \
                    column vector"
            .to_string());
    }

    let Expr::Matrix(a_mat) = a else {

        unreachable!()
    };

    let Expr::Matrix(b_mat) = b else {

        unreachable!()
    };

    let mut augmented_mat =
        a_mat.clone();

    for i in 0 .. a_rows {

        augmented_mat[i]
            .push(b_mat[i][0].clone());
    }

    let rref_expr = rref(
        &Expr::Matrix(augmented_mat),
    )?;

    let Expr::Matrix(rref_mat) =
        rref_expr
    else {

        unreachable!()
    };

    for row in rref_mat.iter().take(a_rows) {

        let is_lhs_zero = row
            [0 .. a_cols]
            .iter()
            .all(is_zero);

        if is_lhs_zero
            && !is_zero(
                &row[a_cols],
            )
        {

            return Ok(
                Expr::NoSolution,
            );
        }
    }

    let mut pivot_cols = Vec::new();

    let mut lead = 0;

    for row in rref_mat.iter().take(a_rows) {

        if lead >= a_cols {

            break;
        }

        let mut i = lead;

        while i < a_cols
            && is_zero(&row[i])
        {

            i += 1;
        }

        if i < a_cols {

            pivot_cols.push(i);

            lead = i + 1;
        }
    }

    let free_cols: Vec<usize> = (0
        .. a_cols)
        .filter(|c| {

            !pivot_cols.contains(c)
        })
        .collect();

    if free_cols.is_empty() {

        let mut solution =
            create_empty_matrix(
                a_cols,
                1,
            );

        for (i, &p_col) in pivot_cols
            .iter()
            .enumerate()
        {

            solution[p_col][0] =
                rref_mat[i][a_cols]
                    .clone();
        }

        Ok(Expr::Matrix(
            solution,
        ))
    } else {

        let particular_solution = {

            let mut sol =
                create_empty_matrix(
                    a_cols,
                    1,
                );

            for (i, &p_col) in
                pivot_cols
                    .iter()
                    .enumerate()
            {

                sol[p_col][0] =
                    rref_mat[i][a_cols]
                        .clone();
            }

            sol
        };

        let null_space_basis =
            null_space(a)?;

        Ok(Expr::System(vec![
            Expr::Matrix(
                particular_solution,
            ),
            null_space_basis,
        ]))
    }
}

/// Computes the trace of a square matrix.
/// Computes the trace of a square matrix.
///
/// The trace of a square matrix is the sum of the elements on its main diagonal.
///
/// # Arguments
/// * `matrix` - The square matrix as an `Expr::Matrix`.
///
/// # Returns
/// A `Result` containing an `Expr` representing the trace.
///
/// # Errors
///
/// This function will return an error if the input `matrix` is not a valid matrix
/// or if it is not a square matrix.

pub fn trace(
    matrix: &Expr
) -> Result<Expr, String> {

    let (rows, cols) =
        get_matrix_dims(matrix)
            .ok_or("Invalid matrix")?;

    if rows != cols {

        return Err("Matrix must be \
                    square"
            .to_string());
    }

    let Expr::Matrix(mat) = matrix
    else {

        unreachable!()
    };

    let mut tr =
        Expr::BigInt(BigInt::zero());

    for (i, row) in mat.iter().enumerate() {

        tr = simplify(&Expr::new_add(
            tr,
            row[i].clone(),
        ));
    }

    Ok(tr)
}

/// Computes the characteristic polynomial of a square matrix.
/// Computes the characteristic polynomial of a square matrix.
///
/// The characteristic polynomial `p(λ) = det(A - λI)` is a polynomial whose roots
/// are the eigenvalues of the matrix `A`.
///
/// # Arguments
/// * `matrix` - The square matrix as an `Expr::Matrix`.
/// * `lambda_var` - The name of the variable representing the eigenvalue (e.g., "lambda").
///
/// # Returns
/// A `Result` containing an `Expr` representing the characteristic polynomial.
///
/// # Errors
///
/// This function will return an error if the input `matrix` is not a valid matrix
/// or if it is not a square matrix.

pub fn characteristic_polynomial(
    matrix: &Expr,
    lambda_var: &str,
) -> Result<Expr, String> {

    let (rows, cols) =
        get_matrix_dims(matrix)
            .ok_or("Invalid matrix")?;

    if rows != cols {

        return Err("Matrix must be \
                    square"
            .to_string());
    }

    let lambda = Expr::Variable(
        lambda_var.to_string(),
    );

    let lambda_i = scalar_mul_matrix(
        &lambda,
        &identity_matrix(rows),
    );

    let a_minus_lambda_i =
        sub_matrices(matrix, &lambda_i);

    Ok(determinant(
        &a_minus_lambda_i,
    ))
}

/// Performs LU decomposition of a matrix. A = LU.
/// Performs LU decomposition of a square matrix `A` into a lower triangular matrix `L`
/// and an upper triangular matrix `U`, such that `A = LU`.
///
/// # Arguments
/// * `matrix` - The square matrix as an `Expr::Matrix`.
///
/// # Returns
/// A `Result` containing a tuple `(L, U)` as `Expr::Matrix`.
///
/// # Errors
///
/// This function will return an error if the input `matrix` is not a valid matrix,
/// is not a square matrix, or if the matrix is singular (i.e., a pivot element is zero),
/// preventing decomposition.

pub fn lu_decomposition(
    matrix: &Expr
) -> Result<(Expr, Expr), String> {

    let (rows, cols) =
        get_matrix_dims(matrix)
            .ok_or("Invalid matrix")?;

    if rows != cols {

        return Err("Matrix must be \
                    square for LU \
                    decomposition"
            .to_string());
    }

    let Expr::Matrix(a) = matrix else {

        unreachable!()
    };

    let n = rows;

    let mut l =
        create_empty_matrix(n, n);

    let mut u =
        create_empty_matrix(n, n);

    for (i, row) in l.iter_mut().enumerate() {

        row[i] =
            Expr::BigInt(BigInt::one());
    }

    for j in 0 .. n {

        for i in 0 ..= j {

            let mut sum = Expr::BigInt(
                BigInt::zero(),
            );

            for (k, _item) in u
                .iter()
                .enumerate()
                .take(i)
            {

                sum = simplify(
                    &Expr::new_add(
                        sum,
                        Expr::new_mul(
                            l[i][k]
                                .clone(
                                ),
                            u[k][j]
                                .clone(
                                ),
                        ),
                    ),
                );
            }

            u[i][j] = simplify(
                &Expr::new_sub(
                    a[i][j].clone(),
                    sum,
                ),
            );
        }

        for i in (j + 1) .. n {

            let mut sum = Expr::BigInt(
                BigInt::zero(),
            );

            for (k, _item) in u
                .iter()
                .enumerate()
                .take(j)
            {

                sum = simplify(
                    &Expr::new_add(
                        sum,
                        Expr::new_mul(
                            l[i][k]
                                .clone(
                                ),
                            u[k][j]
                                .clone(
                                ),
                        ),
                    ),
                );
            }

            if is_zero(&u[j][j]) {

                return Err("Matrix is singular and cannot be decomposed.".to_string());
            }

            l[i][j] = simplify(
                &Expr::new_div(
                    simplify(
                        &Expr::new_sub(
                            a[i][j]
                                .clone(
                                ),
                            sum,
                        ),
                    ),
                    u[j][j].clone(),
                ),
            );
        }
    }

    Ok((
        Expr::Matrix(l),
        Expr::Matrix(u),
    ))
}

/// Performs QR decomposition of a matrix using Gram-Schmidt process. A = QR.
///
/// Performs QR decomposition of a matrix `A` into an orthogonal matrix `Q`
/// and an upper triangular matrix `R`, such that `A = QR`.
/// This implementation uses the Gram-Schmidt process.
///
/// # Arguments
/// * `matrix` - The matrix as an `Expr::Matrix`.
///
/// # Returns
/// A `Result` containing a tuple `(Q, R)` as `Expr::Matrix`.
///
/// # Errors
///
/// This function will return an error if the input `matrix` is not a valid matrix.
/// Potential future errors might include issues with normalization (division by zero)
/// if columns are linearly dependent, but current symbolic implementation may defer
/// such simplification issues.

pub fn qr_decomposition(
    matrix: &Expr
) -> Result<(Expr, Expr), String> {

    let (rows, cols) =
        get_matrix_dims(matrix)
            .ok_or("Invalid matrix")?;

    let Expr::Matrix(a) = matrix else {

        unreachable!()
    };

    let mut q_cols = Vec::new();

    let mut r =
        create_empty_matrix(cols, cols);

    for j in 0 .. cols {

        let mut u_j = a
            .iter()
            .map(|row| row[j].clone())
            .collect::<Vec<_>>();

        for i in 0 .. j {

            let q_i: &Vec<Expr> =
                &q_cols[i];

            let mut dot_a_q =
                Expr::BigInt(
                    BigInt::zero(),
                );

            for k in 0 .. rows {

                dot_a_q = simplify(
                    &Expr::new_add(
                        dot_a_q,
                        Expr::new_mul(
                            a[k][j]
                                .clone(
                                ),
                            q_i[k]
                                .clone(
                                ),
                        ),
                    ),
                );
            }

            r[i][j] = dot_a_q;

            for k in 0 .. rows {

                let proj_term =
                    simplify(
                        &Expr::new_mul(
                            r[i][j]
                                .clone(
                                ),
                            q_i[k]
                                .clone(
                                ),
                        ),
                    );

                u_j[k] = simplify(
                    &Expr::new_sub(
                        u_j[k].clone(),
                        proj_term,
                    ),
                );
            }
        }

        let mut norm_u_j_sq =
            Expr::BigInt(BigInt::zero());

        for item in &u_j {

            norm_u_j_sq = simplify(&Expr::new_add(
                norm_u_j_sq,
                Expr::new_pow(
                    item.clone(),
                    Expr::BigInt(BigInt::from(2)),
                ),
            ));
        }

        let norm_u_j =
            simplify(&Expr::new_sqrt(
                norm_u_j_sq,
            ));

        r[j][j] = norm_u_j.clone();

        let mut q_j = Vec::new();

        for item in &u_j {

            q_j.push(simplify(
                &Expr::new_div(
                    item.clone(),
                    norm_u_j.clone(),
                ),
            ));
        }

        q_cols.push(q_j);
    }

    let mut q_rows =
        create_empty_matrix(rows, cols);

    for i in 0 .. rows {

        for (j, _item) in q_cols
            .iter()
            .enumerate()
            .take(cols)
        {

            q_rows[i][j] =
                q_cols[j][i].clone();
        }
    }

    Ok((
        Expr::Matrix(q_rows),
        Expr::Matrix(r),
    ))
}

/// Computes the reduced row echelon form (RREF) of a matrix.
///
/// This function applies Gaussian elimination to transform the matrix into its RREF.
/// It is used for solving linear systems, finding matrix inverses, and determining rank.
///
/// # Arguments
/// * `matrix` - The matrix as an `Expr::Matrix`.
///
/// # Returns
/// A `Result` containing an `Expr::Matrix` representing the RREF.
///
/// # Errors
///
/// This function will return an error if the input `matrix` is not a valid matrix.

pub fn rref(
    matrix: &Expr
) -> Result<Expr, String> {

    let (rows, cols) = get_matrix_dims(
        matrix,
    )
    .ok_or("Invalid matrix for RREF")?;

    let Expr::Matrix(mut mat) =
        matrix.clone()
    else {

        return Err("Input must be a \
                    matrix"
            .to_string());
    };

    let mut lead = 0;

    for r in 0 .. rows {

        if lead >= cols {

            break;
        }

        let mut i = r;

        while is_zero(&mat[i][lead]) {

            i += 1;

            if i == rows {

                i = r;

                lead += 1;

                if lead == cols {

                    return Ok(
                        Expr::Matrix(
                            mat,
                        ),
                    );
                }
            }
        }

        // Swap rows i and r
        mat.swap(i, r);

        // Divide row r by mat[r][lead]
        let val = mat[r][lead].clone();


        for item in &mut mat[r] {

            *item = simplify(
                &Expr::new_div(
                    item.clone(),
                    val.clone(),
                ),
            );
        }

        // Subtract row r from other rows
        let pivot_row = mat[r].clone();
        for (i, row) in mat.iter_mut().enumerate() {

            if i != r {

                let val = row[lead]
                    .clone();

                for (mij, mrj) in row.iter_mut().zip(&pivot_row) {

                    let term = simplify(
                        &Expr::new_mul(
                            val.clone(),
                            mrj.clone(),
                        ),
                    );

                    *mij = simplify(&Expr::new_sub(
                        mij.clone(),
                        term,
                    ));
                }
            }
        }

        lead += 1;
    }

    Ok(Expr::Matrix(mat))
}

/// Computes the null space (kernel) of a matrix.
///
/// The null space of a matrix `A` is the set of all vectors `x` such that `Ax = 0`.
/// This method finds the null space by first computing the RREF of the matrix,
/// identifying pivot and free variables, and then constructing basis vectors.
///
/// # Arguments
/// * `matrix` - The matrix as an `Expr::Matrix`.
///
/// # Returns
/// A `Result` containing an `Expr::Matrix` whose columns form a basis for the null space.
///
/// # Errors
///
/// This function will return an error if the input `matrix` is not a valid matrix
/// or if the `rref` computation fails.

pub fn null_space(
    matrix: &Expr
) -> Result<Expr, String> {

    let (_rows, cols) =
        get_matrix_dims(matrix).ok_or(
            "Invalid matrix for null \
             space",
        )?;

    let rref_matrix = rref(matrix)?;

    let Expr::Matrix(rref_mat) =
        rref_matrix
    else {

        unreachable!()
    };

    let mut pivot_cols = Vec::new();

    let mut lead = 0;

    for row in &rref_mat {

        if lead >= cols {

            break;
        }

        let mut i = lead;

        while i < cols
            && is_zero(&row[i])
        {

            i += 1;
        }

        if i < cols {

            pivot_cols.push(i);

            lead = i + 1;
        }
    }

    let free_cols: Vec<usize> = (0
        .. cols)
        .filter(|c| {

            !pivot_cols.contains(c)
        })
        .collect();

    let mut basis_vectors = Vec::new();

    for &free_col in &free_cols {

        let mut vec =
            create_empty_matrix(
                cols, 1,
            );

        vec[free_col][0] =
            Expr::BigInt(BigInt::one());

        for (i, &pivot_col) in
            pivot_cols
                .iter()
                .enumerate()
        {

            if !is_zero(
                &rref_mat[i][free_col],
            ) {

                vec[pivot_col][0] = simplify(&Expr::new_neg(
                    rref_mat[i][free_col].clone(),
                ));
            }
        }

        basis_vectors.push(vec);
    }

    if basis_vectors.is_empty() {

        return Ok(Expr::Matrix(
            create_empty_matrix(
                cols, 0,
            ),
        ));
    }

    let num_basis_vectors =
        basis_vectors.len();

    let mut result_matrix =
        create_empty_matrix(
            cols,
            num_basis_vectors,
        );

    for (j, vec) in basis_vectors
        .iter()
        .enumerate()
    {

        for i in 0 .. cols {

            result_matrix[i][j] =
                vec[i][0].clone();
        }
    }

    Ok(Expr::Matrix(
        result_matrix,
    ))
}

/// Performs eigenvalue decomposition of a square matrix.
///
/// This function computes the eigenvalues and corresponding eigenvectors of a square matrix.
/// It first finds the characteristic polynomial, solves for its roots (eigenvalues),
/// and then for each eigenvalue, finds the basis for the null space of `(A - λI)`
/// to determine the eigenvectors.
///
/// # Arguments
/// * `matrix` - The square matrix as an `Expr::Matrix`.
///
/// # Returns
/// A `Result` containing a tuple `(eigenvalues, eigenvectors)`.
/// `eigenvalues` is a column vector of eigenvalues.
/// `eigenvectors` is a matrix where each column is an eigenvector.
/// Returns an error string if the matrix is not square or eigenvalues cannot be found.
///
/// # Errors
///
/// This function will return an error if:
/// - The input `matrix` is not a valid matrix or is not square.
/// - `characteristic_polynomial` fails.
/// - `solve` fails to find any eigenvalues.
/// - `null_space` fails during the eigenvector computation.

pub fn eigen_decomposition(
    matrix: &Expr
) -> Result<(Expr, Expr), String> {

    let (rows, cols) =
        get_matrix_dims(matrix)
            .ok_or("Invalid matrix")?;

    if rows != cols {

        return Err("Matrix must be \
                    square for \
                    eigenvalue \
                    decomposition"
            .to_string());
    }

    let n = rows;

    let lambda_var = "lambda";

    let char_poly =
        characteristic_polynomial(
            matrix,
            lambda_var,
        )?;

    let eigenvalues = solve(
        &char_poly,
        lambda_var,
    );

    if eigenvalues.is_empty() {

        return Err("Could not find \
                    eigenvalues."
            .to_string());
    }

    let unique_eigenvalues : Vec<Expr> = eigenvalues
        .into_iter()
        .collect::<std::collections::HashSet<_>>()
        .into_iter()
        .collect();

    let mut all_eigenvectors_matrix =
        create_empty_matrix(n, n);

    let mut current_col = 0;

    for lambda in &unique_eigenvalues {

        let lambda_i =
            scalar_mul_matrix(
                lambda,
                &identity_matrix(n),
            );

        let a_minus_lambda_i =
            sub_matrices(
                matrix,
                &lambda_i,
            );

        let basis = null_space(
            &a_minus_lambda_i,
        )?;

        if let Expr::Matrix(
            basis_vectors,
        ) = basis
        {

            let (_, num_vectors) =
                get_matrix_dims(
                    &Expr::Matrix(
                        basis_vectors
                            .clone(),
                    ),
                )
                .unwrap_or((0, 0));

            for j in 0 .. num_vectors {

                if current_col < n {

                    for (row_dest, row_src) in
                        all_eigenvectors_matrix
                            .iter_mut()
                            .zip(
                                basis_vectors
                                    .iter(),
                            )
                    {

                        row_dest[current_col] =
                            row_src[j].clone();
                    }

                    current_col += 1;
                }
            }
        }
    }

    while current_col < n {

        for row in &mut all_eigenvectors_matrix {

            row[current_col] =
                Expr::Variable("Not_enough_eigenvectors".to_string());
        }

        current_col += 1;
    }

    let eigenvalues_col_vec =
        Expr::Matrix(
            unique_eigenvalues
                .into_iter()
                .map(|e| vec![e])
                .collect(),
        );

    Ok((
        eigenvalues_col_vec,
        Expr::Matrix(
            all_eigenvectors_matrix,
        ),
    ))
}

/// Performs Singular Value Decomposition (SVD) of a matrix `A`.
///
/// SVD decomposes a matrix `A` into three matrices: `U`, `Σ` (Sigma), and `V^T`,
/// such that `A = UΣV^T`.
/// - `U` is an orthogonal matrix whose columns are the left singular vectors.
/// - `Σ` is a diagonal matrix containing the singular values.
/// - `V^T` is the transpose of an orthogonal matrix `V`, whose columns are the right singular vectors.
///
/// # Arguments
/// * `matrix` - The matrix as an `Expr::Matrix`.
///
/// # Returns
/// A `Result` containing a tuple `(U, Σ, V_transpose)` as `Expr::Matrix`.
///
/// # Errors
///
/// This function will return an error if:
/// - The input `matrix` is not a valid matrix.
/// - `eigen_decomposition` fails during computation of eigenvalues/eigenvectors for `A^T A`.
/// - Eigenvalues are not returned as a matrix or fail to be computed numerically.

pub fn svd_decomposition(
    matrix: &Expr
) -> Result<(Expr, Expr, Expr), String>
{

    let (rows, cols) =
        get_matrix_dims(matrix)
            .ok_or("Invalid matrix")?;

    let a_t = transpose_matrix(matrix);

    let a_t_a =
        mul_matrices(&a_t, matrix);

    let (eigenvalues_sq, v_matrix) =
        eigen_decomposition(&a_t_a)?;

    let singular_values_vec =
        if let Expr::Matrix(
            eig_vals_mat,
        ) = &eigenvalues_sq
        {

            let mut singular_values =
                Vec::new();

            for r in eig_vals_mat {

                singular_values.push(
                    simplify(
                        &Expr::new_sqrt(
                            r[0].clone(
                            ),
                        ),
                    ),
                );
            }

            singular_values
        } else {

            return Err(
                "Failed to compute \
                 eigenvalues for SVD"
                    .to_string(),
            );
        };

    let mut sigma_mat =
        create_empty_matrix(rows, cols);

    for (i, sv) in singular_values_vec
        .iter()
        .enumerate()
    {

        if i < rows && i < cols {

            sigma_mat[i][i] =
                sv.clone();
        }
    }

    let sigma = Expr::Matrix(sigma_mat);

    let v_t =
        transpose_matrix(&v_matrix);

    let mut u_cols = Vec::new();

    if let Expr::Matrix(v_mat) =
        &v_matrix
    {

        for (i, _item) in
            singular_values_vec
                .iter()
                .enumerate()
                .take(cols)
        {

            let v_i = Expr::Matrix(
                (0 .. cols)
                    .map(|r| {

                        vec![v_mat[r]
                            [i]
                            .clone()]
                    })
                    .collect(),
            );

            let a_v_i = mul_matrices(
                matrix,
                &v_i,
            );

            let sigma_i =
                &singular_values_vec[i];

            let u_i = if is_zero(
                sigma_i,
            ) {

                Expr::Matrix(
                    create_empty_matrix(
                        rows, 1,
                    ),
                )
            } else {

                scalar_mul_matrix(
                    &simplify(&Expr::new_div(
                        Expr::BigInt(BigInt::one()),
                        sigma_i.clone(),
                    )),
                    &a_v_i,
                )
            };

            if let Expr::Matrix(
                u_i_mat,
            ) = u_i
            {

                u_cols.push(
                    u_i_mat
                        .into_iter()
                        .map(|r| r[0].clone())
                        .collect::<Vec<_>>(),
                );
            }
        }
    }

    let mut u_mat =
        create_empty_matrix(rows, rows);

    for i in 0 .. rows {

        for (j, _item) in u_cols
            .iter()
            .enumerate()
        {

            u_mat[i][j] =
                u_cols[j][i].clone();
        }
    }

    Ok((
        Expr::Matrix(u_mat),
        sigma,
        v_t,
    ))
}

/// Computes the rank of a matrix.
///
/// The rank of a matrix is the number of linearly independent rows or columns.
/// This function computes the rank by converting the matrix to RREF and counting the number of non-zero rows.
///
/// # Arguments
/// * `matrix` - The matrix as an `Expr::Matrix`.
///
/// # Returns
/// A `Result` containing the rank as a `usize`.
///
/// # Errors
///
/// This function will return an error if the input `matrix` is invalid or if `rref` computation fails
/// or does not return a matrix.

pub fn rank(
    matrix: &Expr
) -> Result<usize, String> {

    let rref_matrix = rref(matrix)?;

    if let Expr::Matrix(rows) =
        rref_matrix
    {

        let rank = rows
            .iter()
            .filter(|row| {

                !row.iter()
                    .all(is_zero)
            })
            .count();

        Ok(rank)
    } else {

        Err(
            "RREF did not return a \
             matrix"
                .to_string(),
        )
    }
}

/// Performs Gaussian elimination on a matrix to produce its row echelon form.
///
/// # Arguments
/// * `matrix` - The matrix as an `Expr::Matrix`.
///
/// # Returns
/// A `Result` containing an `Expr::Matrix` in row echelon form.
///
/// # Errors
///
/// This function will return an error if the input `matrix` is invalid.

pub fn gaussian_elimination(
    matrix: &Expr
) -> Result<Expr, String> {

    let (rows, cols) =
        get_matrix_dims(matrix).ok_or(
            "Invalid matrix for \
             Gaussian elimination",
        )?;

    let Expr::Matrix(mut mat) =
        matrix.clone()
    else {

        return Err("Input must be a \
                    matrix"
            .to_string());
    };

    let mut pivot_row = 0;

    for j in 0 .. cols {

        if pivot_row >= rows {

            break;
        }

        let mut i = pivot_row;

        while is_zero(&mat[i][j]) {

            i += 1;

            if i >= rows {

                i = pivot_row;

                break;
            }
        }

        if i < rows
            && !is_zero(&mat[i][j])
        {

            mat.swap(i, pivot_row);

            let pivot_row_vec = mat[pivot_row].clone();

            for row in
                mat.iter_mut().skip(pivot_row + 1)
            {

                let factor = simplify(
                    &Expr::new_div(
                        row[j]
                            .clone(),
                        pivot_row_vec
                            [j]
                            .clone(),
                    ),
                );

                for (mipk, mprk) in row[j..]
                    .iter_mut()
                    .zip(&pivot_row_vec[j..])
                {

                    let term = simplify(&Expr::new_mul(
                        factor.clone(),
                        mprk.clone(),
                    ));

                    *mipk = simplify(&Expr::new_sub(
                        mipk.clone(),
                        term,
                    ));
                }
            }

            pivot_row += 1;
        }
    }

    Ok(Expr::Matrix(mat))
}

/// Checks if a matrix is a zero matrix (all elements are zero).
///
/// # Arguments
/// * `matrix` - The matrix as an `Expr::Matrix`.
///
/// # Returns
/// `true` if all elements are zero, `false` otherwise.
#[must_use]

pub fn is_zero_matrix(
    matrix: &Expr
) -> bool {

    if let Some((_rows, _cols)) =
        get_matrix_dims(matrix)
    {

        let Expr::Matrix(mat) = matrix
        else {

            return false;
        };

        mat.iter()
            .flatten()
            .all(is_zero)
    } else {

        false
    }
}
