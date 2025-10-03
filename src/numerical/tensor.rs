//! # Numerical Tensor Operations
//!
//! This module provides numerical tensor operations, primarily using `ndarray`
//! for efficient multi-dimensional array manipulation. It includes functions
//! for tensor contraction (tensordot), outer product, and Einstein summation (`einsum`).

use ndarray::{ArrayD, IxDyn};

/// Performs tensor contraction between two N-dimensional arrays (tensordot).
///
/// # Arguments
/// * `a`, `b` - The two tensors (`ndarray::ArrayD<f64>`) to contract.
/// * `axes_a`, `axes_b` - The axes to contract for tensor `a` and `b` respectively.
///
/// # Returns
/// The resulting contracted tensor as an `ndarray::ArrayD<f64>`.
pub fn tensordot(
    a: &ArrayD<f64>,
    b: &ArrayD<f64>,
    axes_a: &[usize],
    axes_b: &[usize],
) -> Result<ArrayD<f64>, String> {
    if axes_a.len() != axes_b.len() {
        return Err("Contracted axes must have the same length.".to_string());
    }

    // Validate that contracted dimensions match
    for (&ax_a, &ax_b) in axes_a.iter().zip(axes_b.iter()) {
        if a.shape()[ax_a] != b.shape()[ax_b] {
            return Err(format!(
                "Dimension mismatch on contracted axes: {} != {}",
                a.shape()[ax_a],
                b.shape()[ax_b]
            ));
        }
    }

    // Identify free and contracted axes
    let free_axes_a: Vec<_> = (0..a.ndim()).filter(|i| !axes_a.contains(i)).collect();
    let free_axes_b: Vec<_> = (0..b.ndim()).filter(|i| !axes_b.contains(i)).collect();

    // Permute axes to bring contracted axes to the end for `a` and beginning for `b`
    let perm_a: Vec<_> = free_axes_a.iter().chain(axes_a.iter()).cloned().collect();
    let perm_b: Vec<_> = axes_b.iter().chain(free_axes_b.iter()).cloned().collect();
    let a_perm = a.clone().permuted_axes(perm_a);
    let b_perm = b.clone().permuted_axes(perm_b);

    // Reshape into 2D matrices for multiplication
    let free_dim_a = free_axes_a.iter().map(|&i| a.shape()[i]).product::<usize>();
    let free_dim_b = free_axes_b.iter().map(|&i| b.shape()[i]).product::<usize>();
    let contracted_dim = axes_a.iter().map(|&i| a.shape()[i]).product::<usize>();

    let a_mat = a_perm
        .to_shape((free_dim_a, contracted_dim))
        .unwrap()
        .to_owned();
    let b_mat = b_perm
        .to_shape((contracted_dim, free_dim_b))
        .unwrap()
        .to_owned();

    // Perform matrix multiplication
    let result_mat = a_mat.dot(&b_mat);

    // Reshape result to final tensor shape
    let mut final_shape_dims = Vec::new();
    final_shape_dims.extend(free_axes_a.iter().map(|&i| a.shape()[i]));
    final_shape_dims.extend(free_axes_b.iter().map(|&i| b.shape()[i]));

    Ok(result_mat
        .to_shape(IxDyn(&final_shape_dims))
        .unwrap()
        .to_owned())
}

/// Computes the outer product of two tensors.
///
/// The outer product of two tensors `A` (rank `r`) and `B` (rank `s`)
/// results in a new tensor `C` of rank `r + s`. Each component of `C`
/// is the product of a component from `A` and a component from `B`.
///
/// # Arguments
/// * `a` - The first tensor (`ndarray::ArrayD<f64>`).
/// * `b` - The second tensor (`ndarray::ArrayD<f64>`).
///
/// # Returns
/// The resulting outer product tensor as an `ndarray::ArrayD<f64>`.
pub fn outer_product(a: &ArrayD<f64>, b: &ArrayD<f64>) -> ArrayD<f64> {
    let mut new_shape = a.shape().to_vec();
    new_shape.extend_from_slice(b.shape());

    let a_flat = a.as_slice().unwrap();
    let b_flat = b.as_slice().unwrap();

    let mut result_data = Vec::with_capacity(a.len() * b.len());
    for val_a in a_flat {
        for val_b in b_flat {
            result_data.push(val_a * val_b);
        }
    }

    ArrayD::from_shape_vec(IxDyn(&new_shape), result_data).unwrap()
}

/// Performs tensor operations using Einstein summation convention.
///
/// This is a simplified version that handles two tensors and specific operation strings.
/// Example: `"ij,jk->ik"` for matrix multiplication.
///
/// # Arguments
/// * `op_str` - The Einstein summation string (e.g., "ij,jk->ik").
/// * `tensors` - A slice of references to `ndarray::ArrayD<f64>` tensors.
///
/// # Returns
/// A `Result` containing the resulting `ndarray::ArrayD<f64>` tensor, or an error string
/// if the operation string is not supported or dimensions mismatch.
pub fn einsum(op_str: &str, tensors: &[&ArrayD<f64>]) -> Result<ArrayD<f64>, String> {
    // This is a complex parsing task. The implementation below is a placeholder
    // demonstrating the concept for matrix multiplication.
    if op_str == "ij,jk->ik" && tensors.len() == 2 {
        let a = tensors[0];
        let b = tensors[1];
        if a.ndim() != 2 || b.ndim() != 2 {
            return Err("Inputs must be 2D matrices for 'ij,jk->ik' operation.".to_string());
        }
        // This is matrix multiplication, which ndarray handles with `dot`.
        let a_2d = a.view().into_dimensionality::<ndarray::Ix2>().unwrap();
        let b_2d = b.view().into_dimensionality::<ndarray::Ix2>().unwrap();
        let result_2d = a_2d.dot(&b_2d);
        return Ok(result_2d.into_dyn());
    }
    Err(format!(
        "Einsum operation '{}' is not supported in this version.",
        op_str
    ))
}
