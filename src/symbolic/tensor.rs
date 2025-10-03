//! # Symbolic Tensor Calculus
//!
//! This module provides a `Tensor` struct for symbolic tensor manipulation and calculus.
//! It supports fundamental tensor operations like outer product and contraction, as well as
//! more advanced concepts from differential geometry, including metric tensors, Christoffel
//! symbols, the Riemann curvature tensor, and covariant derivatives.

use crate::symbolic::calculus::differentiate;
use crate::symbolic::core::Expr;
use crate::symbolic::matrix::inverse_matrix;
use crate::symbolic::simplify::simplify;
use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{One, Zero};

#[derive(Clone, Debug, PartialEq)]
pub struct Tensor {
    pub components: Vec<Expr>,
    pub shape: Vec<usize>,
}

impl Tensor {
    /// Creates a new `Tensor` with the given components and shape.
    ///
    /// # Arguments
    /// * `components` - A `Vec<Expr>` containing the tensor's components in row-major order.
    /// * `shape` - A `Vec<usize>` defining the dimensions of the tensor (e.g., `[dim1, dim2, ...]`).
    ///
    /// # Returns
    /// A `Result` containing the new `Tensor` or an error string if the number of components
    /// does not match the product of the shape dimensions.
    pub fn new(components: Vec<Expr>, shape: Vec<usize>) -> Result<Self, String> {
        let expected_len: usize = shape.iter().product();
        if components.len() != expected_len {
            return Err(format!(
                "Number of components ({}) does not match shape ({:?})",
                components.len(),
                shape
            ));
        }
        Ok(Tensor { components, shape })
    }

    /// Returns the rank (order) of the tensor.
    ///
    /// The rank is the number of indices required to uniquely specify each component
    /// of the tensor, equivalent to the length of its shape vector.
    ///
    /// # Returns
    /// The rank as a `usize`.
    pub fn rank(&self) -> usize {
        self.shape.len()
    }

    /// Returns an immutable reference to the component at the specified indices.
    ///
    /// # Arguments
    /// * `indices` - A slice of `usize` representing the indices of the component.
    ///
    /// # Returns
    /// A `Result` containing an immutable reference to the `Expr` component,
    /// or an error string if the number of indices is incorrect or an index is out of bounds.
    pub fn get(&self, indices: &[usize]) -> Result<&Expr, String> {
        if indices.len() != self.rank() {
            return Err("Incorrect number of indices for tensor rank".to_string());
        }
        let mut flat_index = 0;
        let mut stride = 1;
        for (i, &dim) in self.shape.iter().enumerate().rev() {
            if indices[i] >= dim {
                return Err(format!(
                    "Index {} out of bounds for dimension {}",
                    indices[i], i
                ));
            }
            flat_index += indices[i] * stride;
            stride *= dim;
        }
        Ok(&self.components[flat_index])
    }

    pub(crate) fn get_mut(&mut self, indices: &[usize]) -> Result<&mut Expr, String> {
        if indices.len() != self.rank() {
            return Err("Incorrect number of indices for tensor rank".to_string());
        }
        let mut flat_index = 0;
        let mut stride = 1;
        for (i, &dim) in self.shape.iter().enumerate().rev() {
            if indices[i] >= dim {
                return Err(format!(
                    "Index {} out of bounds for dimension {}",
                    indices[i], i
                ));
            }
            flat_index += indices[i] * stride;
            stride *= dim;
        }
        Ok(&mut self.components[flat_index])
    }

    /// Performs tensor addition with another tensor.
    ///
    /// # Arguments
    /// * `other` - The other `Tensor` to add.
    ///
    /// # Returns
    /// A `Result` containing a new `Tensor` representing the sum,
    /// or an error string if the tensors have incompatible shapes.
    pub fn add(&self, other: &Tensor) -> Result<Tensor, String> {
        if self.shape != other.shape {
            return Err("Tensors must have the same shape for addition".to_string());
        }
        let new_components = self
            .components
            .iter()
            .zip(other.components.iter())
            .map(|(a, b)| simplify(Expr::Add(Box::new(a.clone()), Box::new(b.clone()))))
            .collect();
        Tensor::new(new_components, self.shape.clone())
    }

    /// Performs tensor subtraction with another tensor.
    ///
    /// # Arguments
    /// * `other` - The other `Tensor` to subtract.
    ///
    /// # Returns
    /// A `Result` containing a new `Tensor` representing the difference,
    /// or an error string if the tensors have incompatible shapes.
    pub fn sub(&self, other: &Tensor) -> Result<Tensor, String> {
        if self.shape != other.shape {
            return Err("Tensors must have the same shape for subtraction".to_string());
        }
        let new_components = self
            .components
            .iter()
            .zip(other.components.iter())
            .map(|(a, b)| simplify(Expr::Sub(Box::new(a.clone()), Box::new(b.clone()))))
            .collect();
        Tensor::new(new_components, self.shape.clone())
    }

    /// Multiplies the tensor by a scalar expression.
    ///
    /// Each component of the tensor is multiplied by the given scalar.
    ///
    /// # Arguments
    /// * `scalar` - The `Expr` to multiply the tensor by.
    ///
    /// # Returns
    /// A new `Tensor` representing the result of the scalar multiplication.
    pub fn scalar_mul(&self, scalar: &Expr) -> Tensor {
        let new_components = self
            .components
            .iter()
            .map(|c| simplify(Expr::Mul(Box::new(scalar.clone()), Box::new(c.clone()))))
            .collect();
        Tensor::new(new_components, self.shape.clone()).unwrap()
    }

    /// Computes the outer product of this tensor with another tensor.
    ///
    /// The outer product of two tensors `A` (rank `r`) and `B` (rank `s`)
    /// results in a new tensor `C` of rank `r + s`. Each component of `C`
    /// is the product of a component from `A` and a component from `B`.
    ///
    /// # Arguments
    /// * `other` - The other `Tensor` to compute the outer product with.
    ///
    /// # Returns
    /// A new `Tensor` representing the outer product.
    pub fn outer_product(&self, other: &Tensor) -> Tensor {
        let new_shape: Vec<usize> = self
            .shape
            .iter()
            .chain(other.shape.iter())
            .cloned()
            .collect();
        let mut new_components = Vec::with_capacity(self.components.len() * other.components.len());
        for c1 in &self.components {
            for c2 in &other.components {
                new_components.push(simplify(Expr::Mul(
                    Box::new(c1.clone()),
                    Box::new(c2.clone()),
                )));
            }
        }
        Tensor::new(new_components, new_shape).unwrap()
    }

    /// Contracts two specified axes of the tensor.
    ///
    /// Tensor contraction (also known as trace) is an operation that reduces the rank
    /// of a tensor by summing over components where two indices are equal.
    /// The dimensions of the contracted axes must be equal.
    ///
    /// # Arguments
    /// * `axis1` - The index of the first axis to contract.
    /// * `axis2` - The index of the second axis to contract.
    ///
    /// # Returns
    /// A `Result` containing a new `Tensor` with reduced rank,
    /// or an error string if axes are out of bounds or have unequal dimensions.
    pub fn contract(&self, axis1: usize, axis2: usize) -> Result<Tensor, String> {
        if axis1 >= self.rank() || axis2 >= self.rank() {
            return Err("Axis out of bounds".to_string());
        }
        if self.shape[axis1] != self.shape[axis2] {
            return Err("Dimensions of contracted axes must be equal".to_string());
        }
        if axis1 == axis2 {
            return Err("Cannot contract an axis with itself".to_string());
        }

        let mut new_shape = self.shape.clone();
        let dim = self.shape[axis1];
        new_shape.remove(axis1.max(axis2));
        new_shape.remove(axis1.min(axis2));

        let new_len: usize = if new_shape.is_empty() {
            1
        } else {
            new_shape.iter().product()
        };
        let new_components = vec![Expr::BigInt(BigInt::zero()); new_len];
        let mut new_tensor = Tensor::new(new_components, new_shape.clone())?;

        let mut current_indices = vec![0; self.rank()];
        loop {
            let mut sum_val = Expr::BigInt(BigInt::zero());
            for i in 0..dim {
                current_indices[axis1] = i;
                current_indices[axis2] = i;
                sum_val = simplify(Expr::Add(
                    Box::new(sum_val),
                    Box::new(self.get(&current_indices)?.clone()),
                ));
            }

            let new_indices: Vec<usize> = current_indices
                .iter()
                .enumerate()
                .filter(|(idx, _)| *idx != axis1 && *idx != axis2)
                .map(|(_, &val)| val)
                .collect();

            if new_tensor.rank() > 0 {
                *new_tensor.get_mut(&new_indices)? = sum_val;
            } else {
                new_tensor.components[0] = sum_val;
            }

            // Increment indices
            let mut carry = self.rank() - 1;
            while carry > 0 {
                current_indices[carry] += 1;
                if current_indices[carry] < self.shape[carry] {
                    break;
                }
                current_indices[carry] = 0;
                carry -= 1;
            }
            if carry == 0 && current_indices[0] >= self.shape[0] {
                break;
            }
        }

        Ok(new_tensor)
    }

    /// Converts a rank-2 tensor into an `Expr::Matrix`.
    ///
    /// This is useful for interoperability with matrix operations defined elsewhere.
    ///
    /// # Returns
    /// A `Result` containing an `Expr::Matrix` or an error string if the tensor is not rank-2.
    pub fn to_matrix_expr(&self) -> Result<Expr, String> {
        if self.rank() != 2 {
            return Err("Can only convert a rank-2 tensor to a matrix expression.".to_string());
        }
        Ok(Expr::Matrix(
            self.components
                .chunks(self.shape[1])
                .map(|c| c.to_vec())
                .collect(),
        ))
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct MetricTensor {
    pub g: Tensor,
    pub g_inv: Tensor,
}

impl MetricTensor {
    /// Creates a new `MetricTensor` from a rank-2 tensor representing the metric `g`.
    ///
    /// This constructor also computes and stores the inverse metric `g_inv`,
    /// which is essential for raising and lowering indices.
    ///
    /// # Arguments
    /// * `g` - A rank-2 `Tensor` representing the metric tensor.
    ///
    /// # Returns
    /// A `Result` containing the new `MetricTensor` or an error string if `g` is not
    /// a square rank-2 tensor or cannot be inverted.
    pub fn new(g: Tensor) -> Result<Self, String> {
        if g.rank() != 2 || g.shape[0] != g.shape[1] {
            return Err("Metric tensor must be a square rank-2 tensor.".to_string());
        }
        let g_matrix = Expr::Matrix(
            g.components
                .chunks(g.shape[1])
                .map(|c| c.to_vec())
                .collect(),
        );
        let g_inv_matrix = inverse_matrix(&g_matrix);
        let g_inv = if let Expr::Matrix(rows) = g_inv_matrix {
            Tensor::new(rows.into_iter().flatten().collect(), g.shape.clone())?
        } else {
            return Err("Failed to invert metric tensor".to_string());
        };
        Ok(MetricTensor { g, g_inv })
    }

    /// Raises an index of a covector (rank-1 tensor with lower index) to a vector (upper index).
    ///
    /// This operation uses the inverse metric tensor `g^ij` to transform a covector `v_j`
    /// into a vector `v^i = g^ij * v_j`.
    ///
    /// # Arguments
    /// * `covector` - A rank-1 `Tensor` representing the covector.
    ///
    /// # Returns
    /// A `Result` containing a new `Tensor` representing the vector with the raised index,
    /// or an error string if the input is not a rank-1 tensor.
    pub fn raise_index(&self, covector: &Tensor) -> Result<Tensor, String> {
        if covector.rank() != 1 {
            return Err("Can only raise index of a rank-1 tensor (covector).".to_string());
        }
        let product = self.g_inv.outer_product(covector);
        product.contract(1, 2)
    }

    /// Lowers an index of a vector (rank-1 tensor with upper index) to a covector (lower index).
    ///
    /// This operation uses the metric tensor `g_ij` to transform a vector `v^j`
    /// into a covector `v_i = g_ij * v^j`.
    ///
    /// # Arguments
    /// * `vector` - A rank-1 `Tensor` representing the vector.
    ///
    /// # Returns
    /// A `Result` containing a new `Tensor` representing the covector with the lowered index,
    /// or an error string if the input is not a rank-1 tensor.
    pub fn lower_index(&self, vector: &Tensor) -> Result<Tensor, String> {
        if vector.rank() != 1 {
            return Err("Can only lower index of a rank-1 tensor (vector).".to_string());
        }
        let product = self.g.outer_product(vector);
        product.contract(1, 2)
    }
}

/// Computes the Christoffel symbols of the first kind `Γ_{ijk}`.
///
/// The Christoffel symbols describe the connection coefficients of a metric tensor.
/// They are used to define covariant derivatives and curvature.
/// The first kind symbols are defined as:
/// `Γ_{ijk} = 1/2 * (∂_j g_{ik} + ∂_i g_{jk} - ∂_k g_{ij})`.
///
/// # Arguments
/// * `metric` - The `MetricTensor` of the manifold.
/// * `vars` - A slice of string slices representing the coordinate variables.
///
/// # Returns
/// A `Result` containing a rank-3 `Tensor` representing the Christoffel symbols of the first kind,
/// or an error string if dimensions mismatch.
pub fn christoffel_symbols_first_kind(
    metric: &MetricTensor,
    vars: &[&str],
) -> Result<Tensor, String> {
    let dim = metric.g.shape[0];
    if vars.len() != dim {
        return Err("Number of variables must match metric dimension".to_string());
    }
    let mut components = Vec::new();
    for i in 0..dim {
        for j in 0..dim {
            for k in 0..dim {
                let g_ik = metric.g.get(&[i, k])?;
                let g_jk = metric.g.get(&[j, k])?;
                let g_ij = metric.g.get(&[i, j])?;

                let d_g_ik_dj = differentiate(g_ik, vars[j]);
                let d_g_jk_di = differentiate(g_jk, vars[i]);
                let d_g_ij_dk = differentiate(g_ij, vars[k]);

                let term1 = simplify(Expr::Add(Box::new(d_g_ik_dj), Box::new(d_g_jk_di)));
                let term2 = simplify(Expr::Sub(Box::new(term1), Box::new(d_g_ij_dk)));
                let christoffel = simplify(Expr::Mul(
                    Box::new(Expr::Rational(BigRational::new(
                        BigInt::one(),
                        BigInt::from(2),
                    ))),
                    Box::new(term2),
                ));
                components.push(christoffel);
            }
        }
    }
    Tensor::new(components, vec![dim, dim, dim])
}

/// Computes the Christoffel symbols of the second kind `Γ^i_{jk}`.
///
/// These symbols are obtained by raising the first index of the Christoffel symbols
/// of the first kind using the inverse metric tensor `g^il`:
/// `Γ^i_{jk} = g^il * Γ_{ljk}`.
///
/// # Arguments
/// * `metric` - The `MetricTensor` of the manifold.
/// * `vars` - A slice of string slices representing the coordinate variables.
///
/// # Returns
/// A `Result` containing a rank-3 `Tensor` representing the Christoffel symbols of the second kind,
/// or an error string if dimensions mismatch or computation fails.
pub fn christoffel_symbols_second_kind(
    metric: &MetricTensor,
    vars: &[&str],
) -> Result<Tensor, String> {
    let christoffel_1st = christoffel_symbols_first_kind(metric, vars)?;
    let product = metric.g_inv.outer_product(&christoffel_1st);
    // Contract g^{il} with Γ_{ljk} on l
    product.contract(1, 2)
}

/// Computes the Riemann curvature tensor `R^i_{jkl}`.
///
/// The Riemann curvature tensor is a fundamental object in differential geometry
/// that describes the curvature of Riemannian manifolds. It quantifies the failure
/// of parallel transport to preserve the direction of a vector.
/// It is defined in terms of Christoffel symbols and their derivatives.
///
/// # Arguments
/// * `metric` - The `MetricTensor` of the manifold.
/// * `vars` - A slice of string slices representing the coordinate variables.
///
/// # Returns
/// A `Result` containing a rank-4 `Tensor` representing the Riemann curvature tensor,
/// or an error string if computation fails.
pub fn riemann_curvature_tensor(metric: &MetricTensor, vars: &[&str]) -> Result<Tensor, String> {
    let dim = metric.g.shape[0];
    let christoffel_2nd = christoffel_symbols_second_kind(metric, vars)?;
    let mut components = Vec::new();

    for i in 0..dim {
        for j in 0..dim {
            for k in 0..dim {
                for l in 0..dim {
                    // ∂_k Γ^i_{jl}
                    let term1 = differentiate(christoffel_2nd.get(&[i, j, l])?, vars[k]);
                    // ∂_l Γ^i_{jk}
                    let term2 = differentiate(christoffel_2nd.get(&[i, j, k])?, vars[l]);

                    let mut term3 = Expr::BigInt(BigInt::zero());
                    for m in 0..dim {
                        let g_mjl = christoffel_2nd.get(&[m, j, l])?;
                        let g_imk = christoffel_2nd.get(&[i, m, k])?;
                        term3 = simplify(Expr::Add(
                            Box::new(term3),
                            Box::new(Expr::Mul(Box::new(g_mjl.clone()), Box::new(g_imk.clone()))),
                        ));
                    }

                    let mut term4 = Expr::BigInt(BigInt::zero());
                    for m in 0..dim {
                        let g_mjk = christoffel_2nd.get(&[m, j, k])?;
                        let g_iml = christoffel_2nd.get(&[i, m, l])?;
                        term4 = simplify(Expr::Add(
                            Box::new(term4),
                            Box::new(Expr::Mul(Box::new(g_mjk.clone()), Box::new(g_iml.clone()))),
                        ));
                    }

                    let r_ijkl = simplify(Expr::Sub(
                        Box::new(simplify(Expr::Add(Box::new(term1), Box::new(term3)))),
                        Box::new(simplify(Expr::Add(Box::new(term2), Box::new(term4)))),
                    ));
                    components.push(r_ijkl);
                }
            }
        }
    }
    Tensor::new(components, vec![dim, dim, dim, dim])
}

/// Computes the covariant derivative of a vector field `V^i` with respect to a coordinate `x^k`.
///
/// The covariant derivative `∇_k V^i` accounts for the change in the vector field itself
/// and the change in the basis vectors due to the curvature of the manifold.
/// It is defined as `∇_k V^i = ∂_k V^i + Γ^i_{jk} V^j`.
///
/// # Arguments
/// * `vector_field` - A rank-1 `Tensor` representing the vector field.
/// * `metric` - The `MetricTensor` of the manifold.
/// * `vars` - A slice of string slices representing the coordinate variables.
///
/// # Returns
/// A `Result` containing a rank-2 `Tensor` representing the covariant derivative,
/// or an error string if dimensions mismatch or computation fails.
pub fn covariant_derivative_vector(
    vector_field: &Tensor,
    metric: &MetricTensor,
    vars: &[&str],
) -> Result<Tensor, String> {
    if vector_field.rank() != 1 {
        return Err("Input must be a vector field (rank-1 tensor)".to_string());
    }
    let dim = vector_field.shape[0];
    let christoffel_2nd = christoffel_symbols_second_kind(metric, vars)?;
    let mut components = Vec::new();

    for i in 0..dim {
        //for k in 0..dim {
        for (k, _item) in vars.iter().enumerate().take(dim) {
            // ∂_k V^i
            let partial_deriv = differentiate(vector_field.get(&[i])?, vars[k]);

            let mut christoffel_term = Expr::BigInt(BigInt::zero());
            for j in 0..dim {
                let g_ijk = christoffel_2nd.get(&[i, j, k])?;
                let v_j = vector_field.get(&[j])?;
                christoffel_term = simplify(Expr::Add(
                    Box::new(christoffel_term),
                    Box::new(Expr::Mul(Box::new(g_ijk.clone()), Box::new(v_j.clone()))),
                ));
            }

            let nabla_v = simplify(Expr::Add(
                Box::new(partial_deriv),
                Box::new(christoffel_term),
            ));
            components.push(nabla_v);
        }
    }
    Tensor::new(components, vec![dim, dim])
}
