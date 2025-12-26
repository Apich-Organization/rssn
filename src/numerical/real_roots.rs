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
#[must_use]

pub fn sturm_sequence(poly: &Polynomial) -> Vec<Polynomial> {

    let mut seq = Vec::new();

    if poly.coeffs.is_empty() || (poly.coeffs.len() == 1 && poly.coeffs[0] == 0.0) {

        return seq;
    }

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

        let (_, mut remainder) = p_prev.clone().long_division(&p_curr.clone());

        if remainder.coeffs.is_empty()
            || (remainder.coeffs.len() == 1 && remainder.coeffs[0] == 0.0)
        {

            break;
        }

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

            let changes_a = count_sign_changes(&seq, a);

            let changes_b = count_sign_changes(&seq, b);

            let num_roots = changes_a.saturating_sub(changes_b);

            if num_roots >= 1 {

                // Interval is small but contains roots. We can't separate further within precision.
                // Treat as one interval.
                roots.push((a, b));
            }

            continue;
        }

        let changes_a = count_sign_changes(&seq, a);

        let changes_b = count_sign_changes(&seq, b);

        let num_roots = changes_a.saturating_sub(changes_b);

        if num_roots == 0 {

            continue;
        } else if num_roots == 1 {

            roots.push((a, b));
        } else {

            let mid = f64::midpoint(a, b);

            stack.push((a, mid));

            stack.push((mid, b));
        }
    }

    roots.sort_by(|a, b| {

        a.0.partial_cmp(&b.0).unwrap_or_else(|| {
            if a.0.is_nan() && !b.0.is_nan() {

                std::cmp::Ordering::Greater
            } else if !a.0.is_nan() && b.0.is_nan() {

                std::cmp::Ordering::Less
            } else {

                std::cmp::Ordering::Equal
            }
        })
    });

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

/// Refines a root within an isolated interval using the bisection method.
///
/// # Arguments
/// * `poly` - The polynomial.
/// * `interval` - A tuple `(min, max)` known to contain exactly one root.
/// * `tolerance` - The convergence tolerance.
///
/// # Returns
/// The approximated root.

pub fn refine_root_bisection(poly: &Polynomial, interval: (f64, f64), tolerance: f64) -> f64 {

    let (mut a, mut b) = interval;

    let mut mid = f64::midpoint(a, b);

    // Check endpoints first
    if poly.eval(a).abs() < tolerance {

        return a;
    }

    if poly.eval(b).abs() < tolerance {

        return b;
    }

    while (b - a) > tolerance {

        mid = f64::midpoint(a, b);

        if poly.eval(mid) == 0.0 {

            return mid;
        }

        if poly.eval(a).signum() * poly.eval(mid).signum() < 0.0 {

            b = mid;
        } else {

            a = mid;
        }
    }

    mid
}

/// Finds all distinct real roots of a polynomial.
///
/// This function first isolates roots using Sturm sequences and then refines them
/// using the bisection method.
///
/// # Arguments
/// * `poly` - The polynomial.
/// * `tolerance` - The tolerance for root refinement.
///
/// # Returns
/// A `Result` containing a sorted `Vec<f64>` of roots.

pub fn find_roots(poly: &Polynomial, tolerance: f64) -> Result<Vec<f64>, String> {

    // Stage 1: Isolate roots. Use a coarser precision for isolation to save time,
    // but ensuring it's fine enough to separate close roots.
    // Sturm guarantees separation, so precision here just needs to be reasonable
    // to stop the isolation recursion. 0.1 is often too coarse if roots are close,
    // but the isolation logic splits until count is 1. The 'precision' arg in 'isolate_real_roots'
    // is a stop condition for when to give up if the interval is too small but count > 1?
    // We use the requested tolerance for isolation to ensure we separate close roots.
    // If roots are closer than 'tolerance', they might be treated as a single root (or cluster),
    // which is practically acceptable for numerical finding.
    let isolation_precision = tolerance;

    let intervals = isolate_real_roots(poly, isolation_precision).or_else(|e| {

        // If isolation fails (e.g., zero polynomial), just return error or empty
        Err(e)
    })?;

    let mut roots = Vec::with_capacity(intervals.len());

    for interval in intervals {

        let root = refine_root_bisection(poly, interval, tolerance);

        roots.push(root);
    }

    // Sort and dedup just in case, though Sturm guarantees distinct intervals
    roots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    Ok(roots)
}
