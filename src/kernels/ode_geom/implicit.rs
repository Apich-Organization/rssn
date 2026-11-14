//! Variable-step BDF1/BDF2 for residual-form DAEs with real local-error control.
use super::DaeOptions;
use super::GeomSolution;
use super::OdeGeomError;
use super::solve_dense;
use std::cell::Cell;

fn checked<R: FnMut(f64, &[f64], &[f64]) -> Vec<f64> + ?Sized>(
    f: &mut R,
    t: f64,
    y: &[f64],
    yp: &[f64],
) -> Result<Vec<f64>, OdeGeomError> {
    let r = f(t, y, yp);
    if r.len() != y.len() {
        return Err(OdeGeomError::ShapeMismatch(
            "implicit DAE residual length must equal state length".into(),
        ));
    }
    if r.iter().any(|x| !x.is_finite()) {
        return Err(OdeGeomError::Domain(
            "implicit DAE residual is nonfinite".into(),
        ));
    }
    Ok(r)
}
fn norm(v: &[f64]) -> f64 {
    v.iter().fold(0.0_f64, |m, x| m.max(x.abs()))
}
#[allow(clippy::too_many_arguments)] // One BDF history state, residual and Newton controls.
fn step<R: FnMut(f64, &[f64], &[f64]) -> Vec<f64> + ?Sized>(
    f: &mut R,
    t: f64,
    y: &[f64],
    yp: &[f64],
    previous: Option<(f64, &[f64])>,
    h: f64,
    options: DaeOptions,
) -> Result<(Vec<f64>, Vec<f64>), OdeGeomError> {
    let n = y.len();
    let mut next: Vec<f64> = y.iter().zip(yp).map(|(a, b)| a + h * b).collect();
    let derivative = |next: &[f64]| -> Vec<f64> {
        if let Some((tp, yprev)) = previous {
            let ratio = h / (t - tp);
            let a0 = (1.0 + 2.0 * ratio) / (1.0 + ratio);
            let a1 = -(1.0 + ratio);
            let a2 = ratio * ratio / (1.0 + ratio);
            (0..n)
                .map(|i| (a0 * next[i] + a1 * y[i] + a2 * yprev[i]) / h)
                .collect()
        } else {
            (0..n).map(|i| (next[i] - y[i]) / h).collect()
        }
    };
    let time = t + h;
    for _ in 0..options.max_newton_iter {
        let dy = derivative(&next);
        let residual = checked(f, time, &next, &dy)?;
        if norm(&residual) <= options.newton_tol {
            return Ok((next, dy));
        }
        let mut jac = vec![vec![0.0; n]; n];
        for j in 0..n {
            let epsilon = f64::EPSILON.sqrt() * next[j].abs().max(1.0);
            let mut shifted = next.clone();
            shifted[j] += epsilon;
            let actual = shifted[j] - next[j];
            if actual == 0.0 {
                return Err(OdeGeomError::NotConverged);
            }
            let r = checked(f, time, &shifted, &derivative(&shifted))?;
            for i in 0..n {
                jac[i][j] = (r[i] - residual[i]) / actual;
            }
        }
        let delta = solve_dense(jac, residual.iter().map(|x| -x).collect())?;
        let mut accepted = false;
        let mut damping = 1.0;
        for _ in 0..16 {
            let candidate: Vec<f64> = next
                .iter()
                .zip(&delta)
                .map(|(x, d)| x + damping * d)
                .collect();
            if candidate.iter().all(|x| x.is_finite()) {
                let r = checked(f, time, &candidate, &derivative(&candidate))?;
                if norm(&r) < norm(&residual) {
                    next = candidate;
                    accepted = true;
                    break;
                }
            }
            damping *= 0.5;
        }
        if !accepted {
            return Err(OdeGeomError::NotConverged);
        }
    }
    Err(OdeGeomError::NotConverged)
}

pub(super) fn solve<R: FnMut(f64, &[f64], &[f64]) -> Vec<f64> + ?Sized>(
    residual: &mut R,
    y0: &[f64],
    yp0: &[f64],
    span: (f64, f64),
    options: DaeOptions,
) -> Result<GeomSolution, OdeGeomError> {
    super::validate_dae_options(options)?;
    let (start, end) = span;
    let n = y0.len();
    if n == 0 {
        return Err(OdeGeomError::Empty);
    }
    if yp0.len() != n {
        return Err(OdeGeomError::ShapeMismatch(
            "initial derivative length must equal state length".into(),
        ));
    }
    if !start.is_finite()
        || !end.is_finite()
        || end < start
        || !(end - start).is_finite()
        || y0.iter().chain(yp0).any(|x| !x.is_finite())
    {
        return Err(OdeGeomError::Domain(
            "invalid implicit DAE initial data".into(),
        ));
    }
    let used = Cell::new(0usize);
    let exhausted = Cell::new(false);
    let mut f = |t: f64, y: &[f64], yp: &[f64]| {
        if exhausted.get() || used.get() >= options.max_evaluations {
            exhausted.set(true);
            return vec![f64::NAN; y.len()];
        }
        used.set(used.get() + 1);
        residual(t, y, yp)
    };
    let outcome = (|| {
        // Algebraic initial states must already be consistent. Differential
        // derivatives may be corrected when dF/dy' is nonsingular.
        let mut initial_rate = yp0.to_vec();
        for _ in 0..options.max_newton_iter {
            let r = checked(&mut f, start, y0, &initial_rate)?;
            if norm(&r) <= options.newton_tol {
                break;
            }
            let mut jac = vec![vec![0.0; n]; n];
            for j in 0..n {
                let epsilon = f64::EPSILON.sqrt() * initial_rate[j].abs().max(1.0);
                let mut yp = initial_rate.clone();
                yp[j] += epsilon;
                let r1 = checked(&mut f, start, y0, &yp)?;
                for i in 0..n {
                    jac[i][j] = (r1[i] - r[i]) / epsilon;
                }
            }
            let d = solve_dense(jac, r.iter().map(|x| -x).collect()).map_err(|_| {
                OdeGeomError::Domain(
                    "inconsistent initial DAE state; provide a consistent y0 and derivative guess"
                        .into(),
                )
            })?;
            for i in 0..n {
                initial_rate[i] += d[i];
            }
        }
        if norm(&checked(&mut f, start, y0, &initial_rate)?) > options.newton_tol {
            return Err(OdeGeomError::NotConverged);
        }
        let mut result = GeomSolution {
            t: vec![start],
            y: vec![y0.to_vec()],
            rates: vec![initial_rate.clone()],
            nfev: 0,
            rejected: 0,
        };
        let mut t = start;
        let mut y = y0.to_vec();
        let mut rate = initial_rate;
        let mut h = if options.h0 > 0.0 {
            options.h0
        } else {
            (end - start) / 100.0
        }
        .min(options.h_max)
        .max(options.h_min);
        for _ in 0..options.max_steps {
            if t >= end {
                result.nfev = used.get();
                return Ok(result);
            }
            let remaining = end - t;
            h = h.min(remaining).min(options.h_max);
            if !h.is_finite()
                || crate::kernels::exact_eq(t + h, t)
                || (h < options.h_min && remaining > options.h_min)
            {
                return Err(OdeGeomError::NotConverged);
            }
            let next = if h >= remaining { end } else { t + h };
            let middle = f64::midpoint(t, next);
            let previous = if result.t.len() > 1 {
                let i = result.t.len() - 2;
                Some((result.t[i], result.y[i].as_slice()))
            } else {
                None
            };
            let trial = (|| {
                let (coarse, _) = step(&mut f, t, &y, &rate, previous, next - t, options)?;
                let (half, rhalf) = step(&mut f, t, &y, &rate, previous, middle - t, options)?;
                let (fine, rfine) = step(
                    &mut f,
                    middle,
                    &half,
                    &rhalf,
                    Some((t, &y)),
                    next - middle,
                    options,
                )?;
                Ok::<_, OdeGeomError>((coarse, half, rhalf, fine, rfine))
            })();
            if exhausted.get() {
                return Err(OdeGeomError::NotConverged);
            }
            let (coarse, half, rhalf, fine, rfine) = match trial {
                | Ok(x) => x,
                | Err(OdeGeomError::NotConverged | OdeGeomError::Singular) => {
                    result.rejected += 1;
                    h *= 0.5;
                    continue;
                },
                | Err(e) => return Err(e),
            };
            let divisor = if previous.is_some() {
                3.0
            } else {
                1.0
            };
            let mut error = 0.0_f64;
            for i in 0..n {
                let scale = options.atol + options.rtol * y[i].abs().max(fine[i].abs());
                error = error.max(((fine[i] - coarse[i]) / divisor / scale).abs());
            }
            if !error.is_finite() {
                return Err(OdeGeomError::NotConverged);
            }
            let power = if previous.is_some() {
                -1.0 / 3.0
            } else {
                -0.5
            };
            if error <= 1.0 {
                result.t.push(middle);
                result.y.push(half);
                result.rates.push(rhalf);
                result.t.push(next);
                result.y.push(fine.clone());
                result.rates.push(rfine.clone());
                t = next;
                y = fine;
                rate = rfine;
                h *= if super::super::exact_eq(error, 0.0) {
                    2.0
                } else {
                    (0.85 * error.powf(power)).clamp(0.3, 2.0)
                };
            } else {
                result.rejected += 1;
                h *= (0.85 * error.powf(power)).clamp(0.1, 0.8);
            }
        }
        if t >= end {
            result.nfev = used.get();
            Ok(result)
        } else {
            Err(OdeGeomError::NotConverged)
        }
    })();
    if exhausted.get() {
        Err(OdeGeomError::NotConverged)
    } else {
        outcome
    }
}
