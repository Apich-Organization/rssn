//! # Numerical Vector Operations
//!
//! This module provides numerical implementations for N-dimensional vector operations.
//! It includes basic vector arithmetic, products, norms, and spatial relationship calculations.
//! All operations are designed to be efficient and work with standard Rust slices and vectors.
//!
//! ## Key Formulas
//!
//! - **Dot Product**: $\mathbf{a} \cdot \mathbf{b} = \sum `a_i` `b_i`$
//! - **Euclidean Norm ($`L_2`$)**: $||\mathbf{v}||_2 = \sqrt{\sum `v_i^2`}$
//! - **Cross Product (3D)**: $\mathbf{a} \times \mathbf{b} = (`a_2` `b_3` - `a_3` `b_2`, `a_3` `b_1` - `a_1` `b_3`, `a_1` `b_2` - `a_2` `b_1`)$
//! - **Projection**: $\text{proj}_{\mathbf{b}}(\mathbf{a}) = \frac{\mathbf{a} \cdot \mathbf{b}}{||\mathbf{b}||^2} \mathbf{b}$

/// Adds two vectors element-wise.
///
/// Formula: $\mathbf{r}_i = \mathbf{v1}_i + \mathbf{v2}_i$
///
/// # Arguments
/// * `v1` - The first vector.
/// * `v2` - The second vector.
///
/// # Returns
/// * `Ok(Vec<f64>)` - The sum of the two vectors.
/// * `Err(String)` - If the vectors have different lengths.

pub fn vec_add(
    v1 : &[f64],
    v2 : &[f64],
) -> Result<Vec<f64>, String> {

    if v1.len() != v2.len() {

        return Err(format!(
            "Dimension mismatch: v1.len() = {}, v2.len() = {}",
            v1.len(),
            v2.len()
        ));
    }

    Ok(v1
        .iter()
        .zip(v2.iter())
        .map(|(a, b)| a + b)
        .collect())
}

/// Subtracts one vector from another element-wise.
///
/// Formula: $\mathbf{r}_i = \mathbf{v1}_i - \mathbf{v2}_i$
///
/// # Arguments
/// * `v1` - The first vector (minuend).
/// * `v2` - The second vector (subtrahend).
///
/// # Returns
/// * `Ok(Vec<f64>)` - The difference of the two vectors.
/// * `Err(String)` - If the vectors have different lengths.

pub fn vec_sub(
    v1 : &[f64],
    v2 : &[f64],
) -> Result<Vec<f64>, String> {

    if v1.len() != v2.len() {

        return Err(format!(
            "Dimension mismatch: v1.len() = {}, v2.len() = {}",
            v1.len(),
            v2.len()
        ));
    }

    Ok(v1
        .iter()
        .zip(v2.iter())
        .map(|(a, b)| a - b)
        .collect())
}

/// Multiplies a vector by a scalar.
///
/// Formula: $\mathbf{r}_i = s \cdot \mathbf{v}_i$
///
/// # Arguments
/// * `v` - The vector to scale.
/// * `s` - The scalar factor.
#[must_use]

pub fn scalar_mul(
    v : &[f64],
    s : f64,
) -> Vec<f64> {

    v.iter()
        .map(|&a| a * s)
        .collect()
}

/// Computes the dot product of two vectors.
///
/// Formula: $\mathbf{a} \cdot \mathbf{b} = \sum_{i=1}^n `a_i` `b_i`$
///
/// # Arguments
/// * `v1` - The first vector.
/// * `v2` - The second vector.
///
/// # Returns
/// * `Ok(f64)` - The dot product.
/// * `Err(String)` - If the vectors have different lengths.

pub fn dot_product(
    v1 : &[f64],
    v2 : &[f64],
) -> Result<f64, String> {

    if v1.len() != v2.len() {

        return Err(format!(
            "Dimension mismatch: v1.len() = {}, v2.len() = {}",
            v1.len(),
            v2.len()
        ));
    }

    Ok(v1
        .iter()
        .zip(v2.iter())
        .map(|(a, b)| a * b)
        .sum())
}

/// Computes the Euclidean norm ($`L_2`$ norm) of a vector.
///
/// Formula: $||\mathbf{v}||_2 = \sqrt{\sum_{i=1}^n `v_i^2`}$
#[must_use]

pub fn norm(v : &[f64]) -> f64 {

    v.iter()
        .map(|&a| a * a)
        .sum::<f64>()
        .sqrt()
}

/// Computes the Manhattan norm ($`L_1`$ norm) of a vector.
///
/// Formula: $||\mathbf{v}||_1 = \sum_{i=1}^n |`v_i`|$
#[must_use]

pub fn l1_norm(v : &[f64]) -> f64 {

    v.iter()
        .map(|&a| a.abs())
        .sum()
}

/// Computes the Infinity norm ($L_\infty$ norm) of a vector.
///
/// Formula: $||\mathbf{v}||_\infty = \max_{i} |`v_i`|$
#[must_use]

pub fn linf_norm(v : &[f64]) -> f64 {

    v.iter()
        .map(|&a| a.abs())
        .fold(0.0, |max, val| {
            if val > max {

                val
            } else {

                max
            }
        })
}

/// Computes the $`L_p`$ norm of a vector.
///
/// Formula: $||\mathbf{v}||_p = (\sum_{i=1}^n |`v_i|^p)^{1/p`}$
#[must_use]

pub fn lp_norm(
    v : &[f64],
    p : f64,
) -> f64 {

    if p <= 0.0 {

        return f64::NAN;
    }

    if p == 1.0 {

        return l1_norm(v);
    }

    if p.is_infinite() {

        return linf_norm(v);
    }

    v.iter()
        .map(|&a| a.abs().powf(p))
        .sum::<f64>()
        .powf(1.0 / p)
}

/// Normalizes a vector to have a magnitude of 1.
///
/// # Returns
/// * `Ok(Vec<f64>)` - The unit vector.
/// * `Err(String)` - If the vector is a zero vector and cannot be normalized.

pub fn normalize(v : &[f64]) -> Result<Vec<f64>, String> {

    let n = norm(v);

    if n == 0.0 {

        return Err("Cannot normalize a zero vector.".to_string());
    }

    Ok(scalar_mul(
        v,
        1.0 / n,
    ))
}

/// Computes the cross product of two 3D vectors.
///
/// # Arguments
/// * `v1` - First 3D vector.
/// * `v2` - Second 3D vector.
///
/// # Returns
/// * `Ok(Vec<f64>)` - The 3D cross product vector.
/// * `Err(String)` - If the input vectors are not 3D.

pub fn cross_product(
    v1 : &[f64],
    v2 : &[f64],
) -> Result<Vec<f64>, String> {

    if v1.len() != 3 || v2.len() != 3 {

        return Err("Cross product is only defined for 3D vectors.".to_string());
    }

    Ok(vec![
        v1[1].mul_add(
            v2[2],
            -(v1[2] * v2[1]),
        ),
        v1[2].mul_add(
            v2[0],
            -(v1[0] * v2[2]),
        ),
        v1[0].mul_add(
            v2[1],
            -(v1[1] * v2[0]),
        ),
    ])
}

/// Computes the Euclidean distance between two vectors.

pub fn distance(
    v1 : &[f64],
    v2 : &[f64],
) -> Result<f64, String> {

    let diff = vec_sub(v1, v2)?;

    Ok(norm(&diff))
}

/// Computes the angle between two vectors in radians.
///
/// Formula: $\theta = \arccos\left(\frac{\mathbf{a} \cdot \mathbf{b}}{||\mathbf{a}|| ||\mathbf{b}||}\right)$

pub fn angle(
    v1 : &[f64],
    v2 : &[f64],
) -> Result<f64, String> {

    let n1 = norm(v1);

    let n2 = norm(v2);

    if n1 == 0.0 || n2 == 0.0 {

        return Ok(0.0);
    }

    let dot = dot_product(v1, v2)?;

    let cos_theta = (dot / (n1 * n2)).clamp(-1.0, 1.0);

    Ok(cos_theta.acos())
}

/// Projects vector `v1` onto vector `v2`.
///
/// Formula: $\text{proj}_{\mathbf{v2}}(\mathbf{v1}) = \frac{\mathbf{v1} \cdot \mathbf{v2}}{||\mathbf{v2}||^2} \mathbf{v2}$

pub fn project(
    v1 : &[f64],
    v2 : &[f64],
) -> Result<Vec<f64>, String> {

    if v1.len() != v2.len() {

        return Err("Vectors must have the same dimension.".to_string());
    }

    let dot = dot_product(v1, v2)?;

    let dot2 = dot_product(v2, v2)?;

    if dot2 == 0.0 {

        return Ok(vec![0.0; v1.len()]);
    }

    Ok(scalar_mul(
        v2,
        dot / dot2,
    ))
}

/// Reflects vector `v` about a normal vector `n`.
///
/// Formula: $\mathbf{r} = \mathbf{v} - 2(\mathbf{v} \cdot \mathbf{n})\mathbf{n}$
/// Note: The normal vector `n` must be normalized.

pub fn reflect(
    v : &[f64],
    n : &[f64],
) -> Result<Vec<f64>, String> {

    if v.len() != n.len() {

        return Err("Vectors must have the same dimension.".to_string());
    }

    let dot = dot_product(v, n)?;

    let scaled_n = scalar_mul(n, 2.0 * dot);

    vec_sub(v, &scaled_n)
}

/// Linearly interpolates between two vectors.
///
/// Formula: $\mathbf{r} = (1 - t)\mathbf{v1} + t\mathbf{v2}$
///
/// # Arguments
/// * `t` - Interpolation factor, typically in $[0, 1]$.

pub fn lerp(
    v1 : &[f64],
    v2 : &[f64],
    t : f64,
) -> Result<Vec<f64>, String> {

    if v1.len() != v2.len() {

        return Err("Vectors must have the same dimension.".to_string());
    }

    Ok(v1
        .iter()
        .zip(v2.iter())
        .map(|(&a, &b)| (1.0 - t).mul_add(a, t * b))
        .collect())
}

/// Checks if two vectors are orthogonal.

pub fn is_orthogonal(
    v1 : &[f64],
    v2 : &[f64],
    epsilon : f64,
) -> Result<bool, String> {

    let dot = dot_product(v1, v2)?;

    Ok(dot.abs() < epsilon)
}

/// Checks if two vectors are parallel.

pub fn is_parallel(
    v1 : &[f64],
    v2 : &[f64],
    epsilon : f64,
) -> Result<bool, String> {

    if v1.len() != v2.len() {

        return Err("Vectors must have the same dimension.".to_string());
    }

    let n1 = norm(v1);

    let n2 = norm(v2);

    if n1 == 0.0 || n2 == 0.0 {

        return Ok(true);
    }

    let dot = dot_product(v1, v2)?;

    let cos_theta = (dot / (n1 * n2)).abs();

    Ok((1.0 - cos_theta).abs() < epsilon)
}

/// Computes the cosine similarity between two vectors.

pub fn cosine_similarity(
    v1 : &[f64],
    v2 : &[f64],
) -> Result<f64, String> {

    let n1 = norm(v1);

    let n2 = norm(v2);

    if n1 == 0.0 || n2 == 0.0 {

        return Err("Cosine similarity is undefined for zero vectors.".to_string());
    }

    let dot = dot_product(v1, v2)?;

    Ok(dot / (n1 * n2))
}
