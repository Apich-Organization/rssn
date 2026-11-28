//! Common barycentric interpolation weights for real, distinct nodes.

use std::sync::Arc;

use super::reuse::PlanCache;

static NODE_WEIGHTS: PlanCache<Vec<u64>, [f64]> = PlanCache::new_lru(16);
const MAX_CACHED_NODES: usize = 512;

fn normalize(mut weights: Vec<f64>) -> Option<Vec<f64>> {
    let peak = weights.iter().copied().map(f64::abs).fold(0.0, f64::max);
    if !peak.is_normal() {
        return None;
    }
    for weight in &mut weights {
        *weight /= peak;
        if *weight == 0.0 || !weight.is_finite() {
            return None;
        }
    }
    Some(weights)
}

fn differences(
    nodes: &[f64],
    mut use_difference: impl FnMut(usize, usize, f64, f64) -> Result<(), String>,
) -> Result<(), String> {
    for (i, &left) in nodes.iter().enumerate() {
        for (j, &right) in nodes[..i].iter().enumerate() {
            if super::exact_eq(left, right) {
                return Err("barycentric interpolation nodes must be distinct".into());
            }
            use_difference(i, j, left, right)?;
        }
    }
    Ok(())
}

fn logarithmic_weights(nodes: &[f64]) -> Result<Vec<f64>, String> {
    let mut exponents = vec![0.0_f64; nodes.len()];
    let mut signs = vec![1.0_f64; nodes.len()];
    differences(nodes, |i, j, left, right| {
        let difference = left - right;
        let magnitude = if difference.is_finite() {
            difference.abs().log2()
        } else {
            (0.5 * left - 0.5 * right).abs().log2() + 1.0
        };
        exponents[i] -= magnitude;
        exponents[j] -= magnitude;
        signs[i] *= difference.signum();
        signs[j] *= -difference.signum();
        Ok(())
    })?;
    let peak = exponents.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let weights: Vec<_> = exponents
        .into_iter()
        .zip(signs)
        .map(|(exponent, sign)| sign * (exponent - peak).exp2())
        .collect();
    if weights
        .iter()
        .any(|&weight| weight == 0.0 || !weight.is_finite())
    {
        return Err("barycentric weight ratios exceed Float range".into());
    }
    Ok(weights)
}

fn build(nodes: &[f64]) -> Result<Arc<[f64]>, String> {
    if nodes.iter().any(|node| !node.is_finite()) {
        return Err("barycentric interpolation nodes must be finite".into());
    }
    let maximum = nodes.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let minimum = nodes.iter().copied().fold(f64::INFINITY, f64::min);
    let span = maximum - minimum;
    if span.is_normal() {
        // A pair contributes opposite differences to the two products.
        // Evaluate it once, in a common scale chosen for the entire grid.
        let scale = 4.0 / span;
        let mut products = vec![1.0; nodes.len()];
        differences(nodes, |i, j, left, right| {
            let difference = (left - right) * scale;
            products[i] *= difference;
            products[j] *= -difference;
            Ok(())
        })?;
        if products.iter().all(|value| value.is_normal())
            && let Some(weights) = normalize(products.into_iter().map(f64::recip).collect())
        {
            return Ok(weights.into());
        }
    }
    logarithmic_weights(nodes).map(Into::into)
}

/// Relative weights with a common bounded scale. Small grids reuse weights
/// by their exact abscissas, independently of ordinates or query positions.
pub(super) fn weights(nodes: &[f64]) -> Result<Arc<[f64]>, String> {
    if nodes.len() > MAX_CACHED_NODES {
        return build(nodes);
    }
    let key = nodes.iter().map(|&value| value.to_bits()).collect();
    NODE_WEIGHTS.try_get_or_insert_with(key, || build(nodes))
}
