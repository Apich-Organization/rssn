//! # Numerical Vector Operations
//!
//! This module provides numerical implementations for N-dimensional vector operations.
//! It includes basic vector arithmetic (addition, subtraction, scalar multiplication),
//! dot and cross products, magnitude (norm), and distance and angle calculations.

/// Adds two vectors.
///
/// # Arguments
/// * `v1` - The first vector.
/// * `v2` - The second vector.
///
/// # Returns
/// A `Result` containing a new `Vec<f64>` representing the sum, or an error string if dimensions mismatch.
pub fn vec_add(v1: &[f64], v2: &[f64]) -> Result<Vec<f64>, String> {
    if v1.len() != v2.len() {
        return Err("Vectors must have the same dimension.".to_string());
    }
    Ok(v1.iter().zip(v2.iter()).map(|(a, b)| a + b).collect())
}

/// Subtracts a vector from another.
///
/// # Arguments
/// * `v1` - The first vector.
/// * `v2` - The second vector.
///
/// # Returns
/// A `Result` containing a new `Vec<f64>` representing the difference, or an error string if dimensions mismatch.
pub fn vec_sub(v1: &[f64], v2: &[f64]) -> Result<Vec<f64>, String> {
    if v1.len() != v2.len() {
        return Err("Vectors must have the same dimension.".to_string());
    }
    Ok(v1.iter().zip(v2.iter()).map(|(a, b)| a - b).collect())
}

/// Multiplies a vector by a scalar.
///
/// # Arguments
/// * `v` - The vector.
/// * `s` - The scalar.
///
/// # Returns
/// A new `Vec<f64>` representing the scaled vector.
pub fn scalar_mul(v: &[f64], s: f64) -> Vec<f64> {
    v.iter().map(|&a| a * s).collect()
}

/// Computes the dot product of two vectors.
///
/// # Arguments
/// * `v1` - The first vector.
/// * `v2` - The second vector.
///
/// # Returns
/// A `Result` containing the dot product as an `f64`, or an error string if dimensions mismatch.
pub fn dot_product(v1: &[f64], v2: &[f64]) -> Result<f64, String> {
    if v1.len() != v2.len() {
        return Err("Vectors must have the same dimension.".to_string());
    }
    Ok(v1.iter().zip(v2.iter()).map(|(a, b)| a * b).sum())
}

/// Computes the magnitude (L2 norm) of a vector.
///
/// # Arguments
/// * `v` - The vector.
///
/// # Returns
/// The magnitude of the vector as an `f64`.
pub fn norm(v: &[f64]) -> f64 {
    v.iter().map(|&a| a.powi(2)).sum::<f64>().sqrt()
}

/// Computes the cross product of two 3D vectors.
///
/// # Arguments
/// * `v1` - The first 3D vector.
/// * `v2` - The second 3D vector.
///
/// # Returns
/// A `Result` containing a new `Vec<f64>` representing the cross product, or an error string if vectors are not 3D.
pub fn cross_product(v1: &[f64], v2: &[f64]) -> Result<Vec<f64>, String> {
    if v1.len() != 3 || v2.len() != 3 {
        return Err("Cross product is only defined for 3D vectors.".to_string());
    }
    Ok(vec![
        v1[1] * v2[2] - v1[2] * v2[1],
        v1[2] * v2[0] - v1[0] * v2[2],
        v1[0] * v2[1] - v1[1] * v2[0],
    ])
}

/// Computes the Euclidean distance between two vectors.
///
/// # Arguments
/// * `v1` - The first vector.
/// * `v2` - The second vector.
///
/// # Returns
/// A `Result` containing the Euclidean distance as an `f64`, or an error string if dimensions mismatch.
pub fn distance(v1: &[f64], v2: &[f64]) -> Result<f64, String> {
    let diff = vec_sub(v1, v2)?;
    Ok(norm(&diff))
}

/// Computes the angle between two vectors in radians.
///
/// # Arguments
/// * `v1` - The first vector.
/// * `v2` - The second vector.
///
/// # Returns
/// A `Result` containing the angle in radians as an `f64`, or an error string if dimensions mismatch.
pub fn angle(v1: &[f64], v2: &[f64]) -> Result<f64, String> {
    let dot = dot_product(v1, v2)?;
    let norm1 = norm(v1);
    let norm2 = norm(v2);
    if norm1 == 0.0 || norm2 == 0.0 {
        return Ok(0.0);
    }
    Ok((dot / (norm1 * norm2)).acos())
}
