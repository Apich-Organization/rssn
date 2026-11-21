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
//!
//! # Limitations and coverage
//!
//! * Vectors are plain `&[f64]` slices. [`normalize_into`] accepts caller-owned
//!   output and scalar distance uses no intermediate vector. The other
//!   vector-valued results allocate `Vec`s. Only [`cross_product`] is fixed to
//!   three dimensions.
//! * There is no matrix, quaternion or geometric-algebra type here (see
//!   `matrix`, `computer_graphics` and `geometric_algebra`) and no operator
//!   overloading; the API is a set of functions returning `Err` on a length
//!   mismatch.
//! * Degenerate cases return a value, not an error: [`angle`] is `0` when a
//!   vector is zero, [`project`] is the zero vector for a zero direction,
//!   [`is_parallel`] is `true` for a zero vector, and [`reflect`] assumes the
//!   normal is already normalised.
//! * [`is_orthogonal`] / [`is_parallel`] take the caller's `epsilon` and test
//!   the raw or normalised dot product, so the former is an absolute test and
//!   the latter is scale-free; neither is exact for nearly degenerate inputs.

fn same_length(
    v1: &[f64],
    v2: &[f64],
) -> Result<(), String> {
    if v1.len() == v2.len() {
        return Ok(());
    }
    Err(format!(
        "Dimension mismatch: v1.len() = {}, v2.len() = {}",
        v1.len(),
        v2.len()
    ))
}

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
///
/// # Errors
/// Returns an error if the vectors have different lengths.
pub fn vec_add(
    v1: &[f64],
    v2: &[f64],
) -> Result<Vec<f64>, String> {
    same_length(v1, v2)?;

    Ok(v1.iter().zip(v2.iter()).map(|(a, b)| a + b).collect())
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
///
/// # Errors
/// Returns an error if the vectors have different lengths.
pub fn vec_sub(
    v1: &[f64],
    v2: &[f64],
) -> Result<Vec<f64>, String> {
    same_length(v1, v2)?;

    Ok(v1.iter().zip(v2.iter()).map(|(a, b)| a - b).collect())
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
    v: &[f64],
    s: f64,
) -> Vec<f64> {
    v.iter().map(|&a| a * s).collect()
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
///
/// # Errors
/// Returns an error if the vectors have different lengths.
pub fn dot_product(
    v1: &[f64],
    v2: &[f64],
) -> Result<f64, String> {
    same_length(v1, v2)?;

    Ok(v1.iter().zip(v2.iter()).map(|(a, b)| a * b).sum())
}

fn norm_parts(values: impl Iterator<Item = f64> + Clone) -> (Option<f64>, f64) {
    let squares = values.clone().map(|x| x * x).sum::<f64>();
    if squares.is_normal() || squares.is_nan() {
        return (None, squares.sqrt());
    }
    let scale = values.clone().map(f64::abs).fold(0.0, f64::max);
    if scale == 0.0 || scale.is_infinite() {
        return (None, scale);
    }
    let squares = values
        .map(|x| {
            let x = x / scale;
            x * x
        })
        .sum::<f64>();
    (Some(scale), squares.sqrt())
}

/// Computes the Euclidean norm ($`L_2`$ norm) of a vector.
///
/// Formula: $||\mathbf{v}||_2 = \sqrt{\sum_{i=1}^n `v_i^2`}$
#[must_use]
pub fn norm(v: &[f64]) -> f64 {
    let (scale, length) = norm_parts(v.iter().copied());
    scale.unwrap_or(1.0) * length
}

/// Computes the Manhattan norm ($`L_1`$ norm) of a vector.
///
/// Formula: $||\mathbf{v}||_1 = \sum_{i=1}^n |`v_i`|$
#[must_use]
pub fn l1_norm(v: &[f64]) -> f64 {
    v.iter().map(|&a| a.abs()).sum()
}

/// Computes the Infinity norm ($L_\infty$ norm) of a vector.
///
/// Formula: $||\mathbf{v}||_\infty = \max_{i} |`v_i`|$
#[must_use]
pub fn linf_norm(v: &[f64]) -> f64 {
    v.iter()
        .map(|&a| a.abs())
        .fold(0.0, |max, val| if val > max { val } else { max })
}

/// Computes the $`L_p`$ norm of a vector.
///
/// Formula: $||\mathbf{v}||_p = (\sum_{i=1}^n |`v_i|^p)^{1/p`}$
#[must_use]
pub fn lp_norm(
    v: &[f64],
    p: f64,
) -> f64 {
    if p <= 0.0 {
        return f64::NAN;
    }

    if (p - 1.0).abs() < f64::EPSILON {
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
/// # Errors
/// Returns an error if the vector is a zero vector.
pub fn normalize(v: &[f64]) -> Result<Vec<f64>, String> {
    Ok(normalized(v)?.collect())
}

/// Writes the normalized vector into an existing slice of the same length.
/// Finite nonzero vectors can be normalized even when their norm overflows.
/// # Errors
/// Returns an error for a zero vector or an output with a different length.
/// Leaves `out` unchanged on error.
pub fn normalize_into(
    v: &[f64],
    out: &mut [f64],
) -> Result<(), String> {
    if v.len() != out.len() {
        return Err("Normalization output must have the same length as the input.".into());
    }
    for (dest, value) in out.iter_mut().zip(normalized(v)?) {
        *dest = value;
    }
    Ok(())
}

fn normalized(v: &[f64]) -> Result<impl Iterator<Item = f64> + '_, String> {
    let (scale, length) = norm_parts(v.iter().copied());
    if length == 0.0 {
        return Err("Cannot normalize a zero vector.".into());
    }
    let inverse = length.recip();
    Ok(v.iter().map(move |&value| {
        let value = scale.map_or(value, |scale| value / scale);
        value * inverse
    }))
}

/// Computes the cross product of two 3D vectors.
///
/// # Arguments
/// * `v1` - First 3D vector.
/// * `v2` - Second 3D vector.
///
/// # Returns
/// * `Ok(Vec<f64>)` - The 3D cross product vector.
///
/// # Errors
/// Returns an error if the input vectors are not 3D.
pub fn cross_product(
    v1: &[f64],
    v2: &[f64],
) -> Result<Vec<f64>, String> {
    if v1.len() != 3 || v2.len() != 3 {
        return Err("Cross product \
                    is only defined \
                    for 3D vectors."
            .to_string());
    }

    Ok(vec![
        v1[1].mul_add(v2[2], -(v1[2] * v2[1])),
        v1[2].mul_add(v2[0], -(v1[0] * v2[2])),
        v1[0].mul_add(v2[1], -(v1[1] * v2[0])),
    ])
}

/// Computes the Euclidean distance between two vectors without allocating an
/// intermediate vector.
///
/// # Errors
/// Returns an error if the vectors have different lengths.
pub fn distance(
    v1: &[f64],
    v2: &[f64],
) -> Result<f64, String> {
    same_length(v1, v2)?;
    let (scale, length) = norm_parts(v1.iter().zip(v2).map(|(a, b)| a - b));
    Ok(scale.unwrap_or(1.0) * length)
}

fn direction_cosine(
    v1: &[f64],
    v2: &[f64],
) -> Option<f64> {
    let (scale1, length1) = norm_parts(v1.iter().copied());
    let (scale2, length2) = norm_parts(v2.iter().copied());
    if length1 == 0.0 || length2 == 0.0 {
        return None;
    }

    let pairs = v1.iter().zip(v2);
    let denominator = length1 * length2;
    if scale1.is_none()
        && scale2.is_none()
        && (denominator.is_normal() || !length1.is_finite() || !length2.is_finite())
    {
        return Some(pairs.map(|(a, b)| a * b).sum::<f64>() / denominator);
    }

    let unit1 = length1.recip();
    let unit2 = length2.recip();
    Some(
        pairs
            .map(|(&a, &b)| {
                let a = scale1.map_or(a, |s| a / s) * unit1;
                let b = scale2.map_or(b, |s| b / s) * unit2;
                a * b
            })
            .sum(),
    )
}

/// Computes the angle between two vectors in radians.
///
/// Formula: $\theta = \arccos\left(\frac{\mathbf{a} \cdot \mathbf{b}}{||\mathbf{a}|| ||\mathbf{b}||}\right)$
///
/// # Errors
/// Returns an error if the vectors have different lengths.
pub fn angle(
    v1: &[f64],
    v2: &[f64],
) -> Result<f64, String> {
    same_length(v1, v2)?;
    Ok(direction_cosine(v1, v2).map_or(0.0, |cosine| cosine.clamp(-1.0, 1.0).acos()))
}

/// Projects vector `v1` onto vector `v2`.
///
/// Formula: $\text{proj}_{\mathbf{v2}}(\mathbf{v1}) = \frac{\mathbf{v1} \cdot \mathbf{v2}}{||\mathbf{v2}||^2} \mathbf{v2}$
///
/// # Errors
/// Returns an error if the vectors have different lengths.
pub fn project(
    v1: &[f64],
    v2: &[f64],
) -> Result<Vec<f64>, String> {
    if v1.len() != v2.len() {
        return Err("Vectors must \
                    have the same \
                    dimension."
            .to_string());
    }

    let dot = dot_product(v1, v2)?;

    let dot2 = dot_product(v2, v2)?;

    if dot2 == 0.0 {
        return Ok(vec![0.0; v1.len()]);
    }

    Ok(scalar_mul(v2, dot / dot2))
}

/// Reflects vector `v` about a normal vector `n`.
///
/// Formula: $\mathbf{r} = \mathbf{v} - 2(\mathbf{v} \cdot \mathbf{n})\mathbf{n}$
/// Note: The normal vector `n` must be normalized.
/// The result can remain finite even when the intermediate correction overflows.
///
/// # Errors
/// Returns an error if the vectors have different lengths.
pub fn reflect(
    v: &[f64],
    n: &[f64],
) -> Result<Vec<f64>, String> {
    if v.len() != n.len() {
        return Err("Vectors must \
                    have the same \
                    dimension."
            .to_string());
    }

    let dot = dot_product(v, n)?;
    let factor = 2.0 * dot;
    let components = v.iter().zip(n);
    if factor.is_finite() || n.iter().any(|x| x.abs() > 1.0 || x.is_nan()) {
        return Ok(components.map(|(x, normal)| x - normal * factor).collect());
    }

    // The product can exceed the Float range while adding the original
    // coordinate brings the result back into range.
    if dot.is_finite() {
        return Ok(components
            .map(|(&x, &normal)| (-2.0 * normal).mul_add(dot, x))
            .collect());
    }

    if v.iter().any(|x| !x.is_finite()) {
        return Ok(components.map(|(x, normal)| x - normal * factor).collect());
    }
    let scale = linf_norm(v);
    let (sum, correction) =
        components
            .clone()
            .fold((0.0_f64, 0.0), |(sum, correction), (x, normal)| {
                let product = (x / scale) * normal;
                let next = sum + product;
                let residual = if sum.abs() >= product.abs() {
                    (sum - next) + product
                } else {
                    (product - next) + sum
                };
                (next, correction + residual)
            });
    let dot = sum + correction;
    let factor = 2.0 * dot;
    if factor == 0.0 {
        return Ok(v.to_vec());
    }

    Ok(components
        .map(|(&x, &normal)| {
            let signed = if factor.is_sign_negative() {
                -normal
            } else {
                normal
            };
            match crate::kernels::numeric_common::scale_product(signed, &[factor.abs(), scale], &[])
            {
                | Ok(correction) => x - correction,
                | Err(_) => (-2.0 * normal).mul_add(dot, x / scale) * scale,
            }
        })
        .collect())
}

/// Linearly interpolates between two vectors.
///
/// Formula: $\mathbf{r} = (1 - t)\mathbf{v1} + t\mathbf{v2}$
///
/// # Arguments
/// * `t` - Interpolation factor, typically in $[0, 1]$.
///
/// # Errors
/// Returns an error if the vectors have different lengths.
pub fn lerp(
    v1: &[f64],
    v2: &[f64],
    t: f64,
) -> Result<Vec<f64>, String> {
    if v1.len() != v2.len() {
        return Err("Vectors must \
                    have the same \
                    dimension."
            .to_string());
    }

    Ok(v1
        .iter()
        .zip(v2.iter())
        .map(|(&a, &b)| (1.0 - t).mul_add(a, t * b))
        .collect())
}

/// Checks if two vectors are orthogonal.
///
/// # Errors
/// Returns an error if the vectors have different lengths.
pub fn is_orthogonal(
    v1: &[f64],
    v2: &[f64],
    epsilon: f64,
) -> Result<bool, String> {
    let dot = dot_product(v1, v2)?;

    Ok(dot.abs() < epsilon)
}

/// Checks if two vectors are parallel.
///
/// # Errors
/// Returns an error if the vectors have different lengths.
pub fn is_parallel(
    v1: &[f64],
    v2: &[f64],
    epsilon: f64,
) -> Result<bool, String> {
    if v1.len() != v2.len() {
        return Err("Vectors must \
                    have the same \
                    dimension."
            .to_string());
    }

    Ok(direction_cosine(v1, v2).is_none_or(|cosine| (1.0 - cosine.abs()).abs() < epsilon))
}

/// Computes the cosine similarity between two vectors.
///
/// # Errors
/// Returns an error if the vectors have different lengths or if any vector is zero.
pub fn cosine_similarity(
    v1: &[f64],
    v2: &[f64],
) -> Result<f64, String> {
    same_length(v1, v2)?;
    direction_cosine(v1, v2)
        .ok_or_else(|| "Cosine similarity is undefined for zero vectors.".into())
}
