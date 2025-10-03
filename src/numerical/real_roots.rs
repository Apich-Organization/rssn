use crate::numerical::polynomial::Polynomial;

/// Generates the Sturm sequence for a numerical polynomial.
///
/// The Sturm sequence is a sequence of polynomials derived from the original polynomial
/// and its derivative. It is used to determine the number of distinct real roots
/// of a polynomial in a given interval.
///
/// # Arguments
/// * `poly` - The input polynomial as a `Polynomial`.
///
/// # Returns
/// A `Vec<Polynomial>` representing the Sturm sequence.
pub fn sturm_sequence(poly: &Polynomial) -> Vec<Polynomial> {
    let mut seq = Vec::new();
    if poly.coeffs.is_empty() || (poly.coeffs.len() == 1 && poly.coeffs[0] == 0.0) {
        return seq;
    }

    // For now, we assume the polynomial is square-free.
    // A full implementation would compute p / gcd(p, p').
    seq.push(poly.clone());
    let p1 = poly.derivative();
    if p1.coeffs.is_empty() || (p1.coeffs.len() == 1 && p1.coeffs[0] == 0.0) {
        return seq;
    }
    seq.push(p1);

    let mut i = 1;
    while seq[i].coeffs.len() > 1 || seq[i].coeffs[0] != 0.0 {
        let p_prev = &seq[i - 1];
        let p_curr = &seq[i];
        let (_, mut remainder) = p_prev.clone().long_division(p_curr.clone());

        if remainder.coeffs.is_empty()
            || (remainder.coeffs.len() == 1 && remainder.coeffs[0] == 0.0)
        {
            break;
        }
        // Negate the remainder
        for c in &mut remainder.coeffs {
            *c = -*c;
        }
        seq.push(remainder);
        i += 1;
    }
    seq
}

/// Counts sign changes of the Sturm sequence at a point.
pub(crate) fn count_sign_changes(sequence: &[Polynomial], point: f64) -> usize {
    let mut changes = 0;
    let mut last_sign: Option<i8> = None;

    for poly in sequence {
        let val = poly.eval(point);
        let sign = if val > 1e-9 {
            Some(1)
        } else if val < -1e-9 {
            Some(-1)
        } else {
            None
        };

        if let Some(s) = sign {
            if let Some(ls) = last_sign {
                if s != ls {
                    changes += 1;
                }
            }
            last_sign = Some(s);
        }
    }
    changes
}

/// Isolates real roots of a numerical polynomial.
///
/// This function uses a bisection method combined with Sturm's theorem to recursively
/// narrow down intervals until each interval contains exactly one distinct real root.
///
/// # Arguments
/// * `poly` - The input polynomial as a `Polynomial`.
/// * `precision` - The desired maximum width of the isolating intervals.
///
/// # Returns
/// A `Result` containing a `Vec<(f64, f64)>` of tuples, where each tuple `(a, b)`
/// represents an interval `[a, b]` containing exactly one root. Returns an error string
/// if root bounding or counting fails.
pub fn isolate_real_roots(poly: &Polynomial, precision: f64) -> Result<Vec<(f64, f64)>, String> {
    let seq = sturm_sequence(poly);
    let bound = root_bound(poly)?;
    let mut roots = Vec::new();
    let mut stack = vec![(-bound, bound)];

    while let Some((a, b)) = stack.pop() {
        if b - a < precision {
            continue;
        }

        let changes_a = count_sign_changes(&seq, a);
        let changes_b = count_sign_changes(&seq, b);
        let num_roots = changes_a.saturating_sub(changes_b);

        if num_roots == 1 {
            roots.push((a, b));
        } else if num_roots > 1 {
            let mid = (a + b) / 2.0;
            stack.push((a, mid));
            stack.push((mid, b));
        }
    }
    roots.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    Ok(roots)
}

/// Computes an upper bound for the absolute value of the real roots (Cauchy's bound).
pub(crate) fn root_bound(poly: &Polynomial) -> Result<f64, String> {
    if poly.coeffs.is_empty() {
        return Ok(1.0);
    }
    let lc = poly.coeffs[0];
    if lc == 0.0 {
        return Err("Leading coefficient cannot be zero.".to_string());
    }
    let max_coeff = poly
        .coeffs
        .iter()
        .skip(1)
        .map(|c| c.abs())
        .fold(0.0, f64::max);
    Ok(1.0 + max_coeff / lc.abs())
}
