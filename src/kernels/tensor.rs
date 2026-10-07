//! # Numerical Tensor Operations
//!
//! This module provides numerical tensor operations, primarily using `ndarray`
//! for efficient multi-dimensional array manipulation. It includes functions
//! for tensor contraction (tensordot), outer product, and Einstein summation (`einsum`).

use ndarray::ArrayD;
use ndarray::IxDyn;

/// Performs tensor contraction between two N-dimensional arrays (tensordot).
///
/// # Arguments
/// * `a`, `b` - The two tensors (`ndarray::ArrayD<f64>`) to contract.
/// * `axes_a`, `axes_b` - The axes to contract for tensor `a` and `b` respectively.
///
/// # Returns
/// The resulting contracted tensor as an `ndarray::ArrayD<f64>`.
///
/// # Errors
/// Returns an error if the number of axes to contract mismatch, or if the dimensions along these axes do not match.
pub fn tensordot(
    a: &ArrayD<f64>,
    b: &ArrayD<f64>,
    axes_a: &[usize],
    axes_b: &[usize],
) -> Result<ArrayD<f64>, String> {
    if axes_a.len() != axes_b.len() {
        return Err("Contracted axes \
                    must have the \
                    same length."
            .to_string());
    }

    for (&ax_a, &ax_b) in axes_a.iter().zip(axes_b.iter()) {
        if a.shape()[ax_a] != b.shape()[ax_b] {
            return Err(format!(
                "Dimension mismatch \
                 on contracted axes: \
                 {} != {}",
                a.shape()[ax_a],
                b.shape()[ax_b]
            ));
        }
    }

    let free_axes_a: Vec<_> = (0..a.ndim()).filter(|i| !axes_a.contains(i)).collect();

    let free_axes_b: Vec<_> = (0..b.ndim()).filter(|i| !axes_b.contains(i)).collect();

    let perm_a: Vec<_> = free_axes_a.iter().chain(axes_a.iter()).copied().collect();

    let perm_b: Vec<_> = axes_b.iter().chain(free_axes_b.iter()).copied().collect();

    let a_perm = a.clone().permuted_axes(perm_a);

    let b_perm = b.clone().permuted_axes(perm_b);

    let free_dim_a = free_axes_a.iter().map(|&i| a.shape()[i]).product::<usize>();

    let free_dim_b = free_axes_b.iter().map(|&i| b.shape()[i]).product::<usize>();

    let contracted_dim = axes_a.iter().map(|&i| a.shape()[i]).product::<usize>();

    let a_mat = a_perm
        .to_shape((free_dim_a, contracted_dim))
        .map_err(|e| e.to_string())?
        .to_owned();

    let b_mat = b_perm
        .to_shape((contracted_dim, free_dim_b))
        .map_err(|e| e.to_string())?
        .to_owned();

    let result_mat = a_mat.dot(&b_mat);

    let mut final_shape_dims = Vec::new();

    final_shape_dims.extend(free_axes_a.iter().map(|&i| a.shape()[i]));

    final_shape_dims.extend(free_axes_b.iter().map(|&i| b.shape()[i]));

    Ok(result_mat
        .to_shape(IxDyn(&final_shape_dims))
        .map_err(|e| e.to_string())?
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
///
/// # Errors
/// Returns an error if input tensors are not contiguous.
pub fn outer_product(
    a: &ArrayD<f64>,
    b: &ArrayD<f64>,
) -> Result<ArrayD<f64>, String> {
    let mut new_shape = a.shape().to_vec();

    new_shape.extend_from_slice(b.shape());

    let a_flat = a.as_slice().ok_or_else(|| {
        "Input tensor 'a' is not \
             contiguous"
            .to_string()
    })?;

    let b_flat = b.as_slice().ok_or_else(|| {
        "Input tensor 'b' is not \
             contiguous"
            .to_string()
    })?;

    let mut result_data = Vec::with_capacity(a.len() * b.len());

    for val_a in a_flat {
        for val_b in b_flat {
            result_data.push(val_a * val_b);
        }
    }

    ArrayD::from_shape_vec(IxDyn(&new_shape), result_data).map_err(|e| e.to_string())
}

/// Performs tensor-vector multiplication.
///
/// # Errors
/// Returns an error if the tensor has zero dimensions or if the last dimension mismatches the vector size.
pub fn tensor_vec_mul(
    tensor: &ArrayD<f64>,
    vector: &[f64],
) -> Result<ArrayD<f64>, String> {
    if tensor.ndim() < 1 {
        return Err("Tensor must \
                    have at least \
                    one dimension."
            .to_string());
    }

    let last_dim = tensor.shape()[tensor.ndim() - 1];

    if last_dim != vector.len() {
        return Err(format!(
            "Dimension mismatch: last \
             tensor dim {} != vector \
             length {}",
            last_dim,
            vector.len()
        ));
    }

    let vec_arr = ndarray::Array1::from_vec(vector.to_vec());

    let res = tensordot(tensor, &vec_arr.into_dyn(), &[tensor.ndim() - 1], &[0])?;

    Ok(res)
}

/// Computes the inner product of two tensors of the same shape.
///
/// # Errors
/// Returns an error if the shapes mismatch or if tensors are not contiguous.
pub fn inner_product(
    a: &ArrayD<f64>,
    b: &ArrayD<f64>,
) -> Result<f64, String> {
    if a.shape() != b.shape() {
        return Err("Tensors must \
                    have the same \
                    shape for inner \
                    product."
            .to_string());
    }

    let a_flat = a.as_slice().ok_or("Tensor 'a' is not contiguous")?;

    let b_flat = b.as_slice().ok_or("Tensor 'b' is not contiguous")?;

    Ok(a_flat.iter().zip(b_flat.iter()).map(|(x, y)| x * y).sum())
}

/// Contracts a single tensor along two specified axes (a generalised
/// trace): the result has rank `ndim - 2` and entries
/// `sum_t a[..., t, ..., t, ...]`, the remaining axes in their original
/// order.
///
/// # Errors
/// Returns an error if the axes are equal, out of range, or of different
/// lengths.
pub fn contract(
    a: &ArrayD<f64>,
    axis1: usize,
    axis2: usize,
) -> Result<ArrayD<f64>, String> {
    if axis1 == axis2 {
        return Err("Axes must be different for contraction.".to_string());
    }
    let (Some(&n), Some(&m)) = (a.shape().get(axis1), a.shape().get(axis2)) else {
        return Err("Contraction axis out of range.".to_string());
    };
    if n != m {
        return Err("Dimensions along contraction axes must be equal.".to_string());
    }
    // Remove the higher axis first so the lower index stays valid.
    let (low, high) = if axis1 < axis2 { (axis1, axis2) } else { (axis2, axis1) };
    let shape: Vec<usize> =
        a.shape().iter().enumerate().filter(|&(i, _)| i != low && i != high).map(|(_, &d)| d).collect();
    let mut out = ArrayD::<f64>::zeros(IxDyn(&shape));
    for t in 0..n {
        let diagonal = a.index_axis(ndarray::Axis(high), t);
        let diagonal = diagonal.index_axis(ndarray::Axis(low), t);
        out += &diagonal;
    }
    Ok(out)
}

/// Computes the Frobenius norm of a tensor.
#[must_use]
pub fn norm(a: &ArrayD<f64>) -> f64 {
    a.iter().map(|x| x * x).sum::<f64>().sqrt()
}

use serde::Deserialize;
use serde::Serialize;

/// A serializable representation of an N-dimensional tensor.
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct TensorData {
    /// Dimensions of the tensor.
    pub shape: Vec<usize>,
    /// Flat vector of tensor data.
    pub data: Vec<f64>,
}

impl From<&ArrayD<f64>> for TensorData {
    fn from(arr: &ArrayD<f64>) -> Self {
        Self {
            shape: arr.shape().to_vec(),
            data: arr.clone().into_raw_vec_and_offset().0,
        }
    }
}

impl TensorData {
    /// Converts back to an `ndarray::ArrayD`.
    ///
    /// # Errors
    /// Returns an error if the shape and data are inconsistent.
    pub fn to_arrayd(&self) -> Result<ArrayD<f64>, String> {
        ArrayD::from_shape_vec(IxDyn(&self.shape), self.data.clone()).map_err(|e| e.to_string())
    }
}

#[cfg(test)]
mod tests {

    use ndarray::array;

    use super::*;

    #[test]
    fn test_tensordot() {
        let a = array![[1.0, 2.0], [3.0, 4.0]].into_dyn();

        let b = array![[5.0, 6.0], [7.0, 8.0]].into_dyn();

        let res = tensordot(&a, &b, &[1], &[0]).unwrap();

        // Standard matrix multiplication
        assert_eq!(res.shape(), &[2, 2]);

        assert!((res[[0, 0]] - (1.0 * 5.0 + 2.0 * 7.0)).abs() < 1e-9);
    }

    #[test]
    fn test_outer_product() {
        let a = array![1.0, 2.0].into_dyn();

        let b = array![3.0, 4.0].into_dyn();

        let res = outer_product(&a, &b).unwrap();

        assert_eq!(res.shape(), &[2, 2]);

        assert!((res[[0, 0]] - 3.0).abs() < 1e-9);

        assert!((res[[1, 1]] - 8.0).abs() < 1e-9);
    }
}
