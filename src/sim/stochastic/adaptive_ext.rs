//! Rejection-consistent adaptive stochastic integration on bounded noise trees.
//!
//! Independent finest-cell increments fix one noise realization before stepping.
//! Every coarser increment is a sum of those same cells (iterated Brownian
//! integrals use Chen's relation). No rejected proposal samples replacement noise.
use super::Adapters;
use super::Cell;
use super::DEFAULT_CELLS;
use super::DriftImplicitScheme;
use super::JumpDistribution;
use super::JumpPath;
use super::JumpSpec;
use super::Method;
use super::Options;
use super::Path;
use super::RefCell;
use super::Rng;
use super::Schedule;
use super::SdeError;
use super::Stream;
use super::default_terms;
use super::drift_implicit_core;
use super::ensemble_seed;
use super::kernel;
use super::levy_increments;
use super::map_sde_ext_err;
use super::sde_order_for;
use super::simulate_terminal_with_increments;
use super::uses_levy_areas;
use super::wiener;

const MAX_NOISE_VALUES: usize = 8_388_608;
const MAX_CELLS: usize = 1_048_576;
/// The rejection at the finest noise level: only a finer tree can help.
const FINEST_RESOLUTION: &str =
    "adaptive error exceeds tolerance at the finest coupled noise resolution; increase cells";

#[derive(Clone, Default)]
struct Block {
    dw: Vec<f64>,
    areas: Option<Vec<f64>>,
}
struct NoiseTree {
    cells: usize,
    blocks: Vec<Block>,
}
fn checked_cells(
    requested: usize,
    noise: usize,
    areas: bool,
) -> Result<usize, SdeError> {
    let cells = if requested == 0 {
        DEFAULT_CELLS
    } else {
        requested
    };
    let cells = cells
        .max(2)
        .checked_next_power_of_two()
        .ok_or_else(|| SdeError::invalid("adaptive noise resolution overflows"))?;
    let width = if areas {
        noise.checked_mul(noise).and_then(|n| n.checked_add(noise))
    } else {
        Some(noise)
    };
    if cells > MAX_CELLS
        || width
            .and_then(|n| cells.checked_mul(2)?.checked_mul(n))
            .is_none_or(|n| n > MAX_NOISE_VALUES)
    {
        return Err(SdeError::invalid(
            "adaptive noise tree exceeds the finite allocation budget (1048576 cells, 8388608 coefficients)",
        ));
    }
    Ok(cells)
}
impl NoiseTree {
    fn new(
        increments: Vec<Vec<f64>>,
        areas: Option<Vec<Vec<f64>>>,
    ) -> Result<Self, SdeError> {
        if increments.is_empty()
            || increments[0].is_empty()
            || increments
                .iter()
                .any(|row| row.len() != increments[0].len() || row.iter().any(|v| !v.is_finite()))
            || areas.as_ref().is_some_and(|a| {
                a.len() != increments.len()
                    || a.iter().any(|row| row.iter().any(|v| !v.is_finite()))
            })
        {
            return Err(SdeError::NonFinite(
                "adaptive noise sample has invalid/non-finite storage".into(),
            ));
        }
        let cells = increments.len();
        let noise = increments[0].len();
        let mut blocks = vec![Block::default(); 2 * cells];
        let mut areas = areas.map(Vec::into_iter);
        for (i, dw) in increments.into_iter().enumerate() {
            blocks[cells + i] = Block {
                dw,
                areas: areas.as_mut().and_then(Iterator::next),
            };
        }
        for i in (1..cells).rev() {
            let (left, right) = (&blocks[2 * i], &blocks[2 * i + 1]);
            let dw = left.dw.iter().zip(&right.dw).map(|(a, b)| a + b).collect();
            let areas = left.areas.as_ref().zip(right.areas.as_ref()).map(|(a, b)| {
                let mut result = vec![0.0; noise * noise];
                for j in 0..noise {
                    for k in 0..noise {
                        // Chen's relation uses two distinct noise channels j,k.
                        let left_j = left.dw[j];
                        let right_k = right.dw[k];
                        let cross = left_j * right_k;
                        result[j * noise + k] = a[j * noise + k] + b[j * noise + k] + cross;
                    }
                }
                result
            });
            blocks[i] = Block { dw, areas };
        }
        Ok(Self { cells, blocks })
    }

    fn brownian(
        seed: u64,
        noise: usize,
        span: f64,
        cells: usize,
        areas: bool,
        terms: usize,
    ) -> Result<Self, SdeError> {
        let cells = checked_cells(cells, noise, areas)?;
        let h = span / cells as f64;
        if !h.is_finite() || h <= 0.0 {
            return Err(SdeError::StepLimit(
                "finest Brownian cell is not representable".into(),
            ));
        }
        if areas {
            let terms = if terms == 0 {
                default_terms(h)
            } else {
                terms
            };
            if terms > 4096 {
                return Err(SdeError::invalid(
                    "adaptive Levy-area Fourier budget must not exceed 4096 terms",
                ));
            }
            let work = cells
                .checked_mul(noise)
                .and_then(|n| n.checked_mul(noise))
                .and_then(|n| n.checked_mul(terms));
            if work.is_none_or(|n| n > 268_435_456) {
                return Err(SdeError::invalid(
                    "adaptive Levy-area work exceeds 268435456 coefficient updates",
                ));
            }
            let mut rng = Rng::new(seed);
            let sample = kernel::sample_noise(&mut rng, noise, cells, h, terms);
            Self::new(sample.dws, Some(sample.areas))
        } else {
            Self::new(wiener(seed, noise, cells, h), None)
        }
    }
}

struct Control {
    t0: f64,
    t1: f64,
    rtol: f64,
    atol: f64,
    max_steps: usize,
    trials: usize,
    rejected: usize,
    moment: f64,
    order: f64,
    reference: Vec<f64>,
}
impl Control {
    fn new(
        options: &Options,
        initial: &[f64],
        moment: f64,
    ) -> Result<Self, SdeError> {
        let Schedule::Adaptive {
            t1,
            rtol,
            atol,
            max_steps,
            ..
        } = options.schedule
        else {
            return Err(SdeError::invalid(
                "adaptive controller requires an adaptive schedule",
            ));
        };
        Ok(Self {
            t0: options.t0,
            t1,
            rtol,
            atol,
            max_steps,
            trials: 0,
            rejected: 0,
            moment,
            order: sde_order_for(options.method, options.structure).strong,
            reference: initial.iter().map(|v| v.abs()).collect(),
        })
    }

    fn difference(
        &self,
        a: &[f64],
        b: &[f64],
        initial: &[f64],
    ) -> f64 {
        a.iter().zip(b).zip(initial).zip(&self.reference).fold(
            0.0_f64,
            |error, (((a, b), y), reference)| {
                error.max(
                    (a - b).abs()
                        / (self.atol
                            + self.rtol * a.abs().max(b.abs()).max(y.abs()).max(*reference))
                        .max(f64::MIN_POSITIVE),
                )
            },
        )
    }
}

// Both the stochastic and drift-only two-half-step defects are controlled.
// The stochastic allocation uses a finite q-moment: q=2 for Brownian noise and
// q<alpha for stable noise. The drift allocation scales linearly with time.
// These are local coupled-error estimates, not a certified global error bound.
#[allow(clippy::too_many_arguments)] // explicit state, noise and schedule inputs stay in one numerical step contract
fn segment<S, D>(
    control: &mut Control,
    tree: &NoiseTree,
    start: f64,
    end: f64,
    initial: &[f64],
    times: &mut Vec<f64>,
    states: &mut Vec<Vec<f64>>,
    mut step: S,
    mut drift_step: D,
) -> Result<(), SdeError>
where
    S: FnMut(f64, f64, &[f64], &Block) -> Result<Vec<f64>, SdeError>,
    D: FnMut(f64, f64, &[f64]) -> Result<Vec<f64>, SdeError>,
{
    let mut y = initial.to_vec();
    let mut stack = vec![(1usize, 0usize, tree.cells)];
    let time = |index: usize| {
        if index == tree.cells {
            end
        } else {
            start + (end - start) * (index as f64 / tree.cells as f64)
        }
    };
    while let Some((node, first, count)) = stack.pop() {
        control.trials += 1;
        if control.trials > control.max_steps {
            return Err(SdeError::StepLimit(
                "adaptive proposal budget exhausted".into(),
            ));
        }
        let (t, middle, next) = (time(first), time(first + count / 2), time(first + count));
        if middle <= t || next <= middle {
            return Err(SdeError::StepLimit(
                "adaptive time subdivision collapsed".into(),
            ));
        }
        let trial = (|| {
            let coarse = step(t, next - t, &y, &tree.blocks[node])?;
            let left = step(t, middle - t, &y, &tree.blocks[2 * node])?;
            let fine = step(middle, next - middle, &left, &tree.blocks[2 * node + 1])?;
            let drift_coarse = drift_step(t, next - t, &y)?;
            let drift_left = drift_step(t, middle - t, &y)?;
            let drift_fine = drift_step(middle, next - middle, &drift_left)?;
            let fraction = (next - t) / (control.t1 - control.t0);
            let stochastic = control.difference(&fine, &coarse, &y)
                / ((2.0_f64.powf(control.order) - 1.0) * fraction.powf(1.0 / control.moment));
            let bias = control.difference(&drift_fine, &drift_coarse, &y) / fraction;
            Ok::<_, SdeError>((left, fine, stochastic.max(bias)))
        })();
        match trial {
            | Ok((left, fine, error)) if error.is_finite() && error <= 1.0 => {
                if times.len().saturating_add(2) > control.max_steps.saturating_add(1)
                    || times
                        .len()
                        .saturating_add(2)
                        .checked_mul(initial.len())
                        .is_none_or(|n| n > MAX_NOISE_VALUES)
                {
                    return Err(SdeError::StepLimit(
                        "adaptive accepted-path storage/step budget exhausted".into(),
                    ));
                }
                times.push(middle);
                states.push(left);
                times.push(next);
                states.push(fine.clone());
                y = fine;
            },
            | Err(
                error @ (SdeError::Invalid(_) | SdeError::Shape(_) | SdeError::Unsupported { .. }),
            ) => return Err(error),
            | _ => {
                control.rejected += 1;
                if count <= 2 {
                    return Err(SdeError::StepLimit(FINEST_RESOLUTION.into()));
                }
                stack.push((2 * node + 1, first + count / 2, count / 2));
                stack.push((2 * node, first, count / 2));
            },
        }
    }
    Ok(())
}

const fn method_scheme(method: Method) -> Option<DriftImplicitScheme> {
    match method {
        | Method::DriftImplicitEuler => Some(DriftImplicitScheme::BackwardEuler),
        | Method::ThetaMaruyama => Some(DriftImplicitScheme::Trapezoidal),
        | Method::DriftImplicitMilstein => Some(DriftImplicitScheme::DriftImplicitMilstein),
        | _ => None,
    }
}
#[allow(clippy::too_many_arguments)] // explicit state, noise and schedule inputs stay in one numerical step contract
fn brownian_step<D, F>(
    drift: &D,
    diffusion: &F,
    noise: usize,
    method: Method,
    scheme: Option<DriftImplicitScheme>,
    t: f64,
    h: f64,
    y: &[f64],
    block: &Block,
) -> Result<Vec<f64>, SdeError>
where
    D: Fn(f64, &[f64], &mut [f64]),
    F: Fn(f64, &[f64], &mut [f64]),
{
    if let Some(scheme) = scheme.or_else(|| method_scheme(method)) {
        if let Some(areas) = &block.areas {
            return crate::kernels::sde_ext::drift_implicit_step_with_areas(
                drift, diffusion, noise, t, h, y, &block.dw, areas, scheme,
            )
            .map_err(map_sde_ext_err);
        }
        return drift_implicit_core(
            drift,
            diffusion,
            noise,
            t,
            h,
            y,
            std::slice::from_ref(&block.dw),
            scheme,
        )?
        .pop()
        .ok_or_else(|| SdeError::invalid("empty implicit step"));
    }
    if let Some(areas) = &block.areas {
        let sample = kernel::NoiseSample {
            dws: vec![block.dw.clone()],
            areas: vec![areas.clone()],
        };
        let states = match method {
            | Method::MilsteinStratonovich | Method::MilsteinNoncommutativeStratonovich => {
                kernel::milstein_noncommutative_stratonovich(
                    drift, diffusion, noise, t, h, y, &sample,
                )
            },
            | Method::SrkPlaten => {
                kernel::srk_noncommutative(drift, diffusion, noise, t, h, y, &sample)
            },
            | _ => kernel::milstein_noncommutative(drift, diffusion, noise, t, h, y, &sample),
        }
        .map_err(SdeError::from)?;
        return states
            .last()
            .cloned()
            .ok_or_else(|| SdeError::invalid("empty stochastic step"));
    }
    simulate_terminal_with_increments(
        drift,
        diffusion,
        noise,
        t,
        h,
        y,
        std::slice::from_ref(&block.dw),
        method,
    )
}

pub(super) fn brownian<A, G>(
    adapters: &Adapters<A, G>,
    y0: &[f64],
    options: &Options,
    scheme: Option<DriftImplicitScheme>,
) -> Result<Path, SdeError>
where
    A: FnMut(f64, &[f64]) -> Result<Vec<f64>, SdeError>,
    G: FnMut(f64, &[f64]) -> Result<Vec<f64>, SdeError>,
{
    // `cells == 0` asks the solver to choose the noise resolution. A
    // strong-order-1/2 scheme's normalised defect shrinks only like the
    // square root of the cell width, so at a tight tolerance the default
    // tree can be too coarse; a finer one (a new, equally valid realisation
    // for the seed) is tried before giving up. An explicit `cells` is
    // honoured as given.
    let Schedule::Adaptive { cells: 0, .. } = options.schedule else {
        return brownian_on(adapters, y0, options, scheme);
    };
    let first = brownian_on(
        adapters,
        y0,
        &with_cells(options, super::DEFAULT_CELLS),
        scheme,
    );
    match first {
        | Err(SdeError::StepLimit(message)) if message == FINEST_RESOLUTION => {
            brownian_on_finer(adapters, y0, options, scheme, super::DEFAULT_CELLS * 4)
                .unwrap_or(Err(SdeError::StepLimit(message)))
        },
        | other => other,
    }
}

/// `options` with the adaptive schedule's finest cells set to `cells`.
fn with_cells(
    options: &Options,
    cells: usize,
) -> Options {
    let mut attempt = options.clone();
    if let Schedule::Adaptive { cells: slot, .. } = &mut attempt.schedule {
        *slot = cells;
    }
    attempt
}

/// [`brownian`] on trees of `cells`, `4 cells`, ... finest cells, while the
/// finest level is the limit and the tree fits the allocation budget;
/// `None` when every tree that fits still hits the finest level.
fn brownian_on_finer<A, G>(
    adapters: &Adapters<A, G>,
    y0: &[f64],
    options: &Options,
    scheme: Option<DriftImplicitScheme>,
    mut cells: usize,
) -> Option<Result<Path, SdeError>>
where
    A: FnMut(f64, &[f64]) -> Result<Vec<f64>, SdeError>,
    G: FnMut(f64, &[f64]) -> Result<Vec<f64>, SdeError>,
{
    loop {
        if cells > MAX_CELLS {
            return None;
        }
        match brownian_on(adapters, y0, &with_cells(options, cells), scheme) {
            | Err(SdeError::StepLimit(message)) if message == FINEST_RESOLUTION => {
                cells = cells.saturating_mul(4);
            },
            // The finer tree does not fit the budget: the coarser failure
            // is the meaningful one.
            | Err(SdeError::Invalid(_)) => return None,
            | other => return Some(other),
        }
    }
}

fn brownian_on<A, G>(
    adapters: &Adapters<A, G>,
    y0: &[f64],
    options: &Options,
    scheme: Option<DriftImplicitScheme>,
) -> Result<Path, SdeError>
where
    A: FnMut(f64, &[f64]) -> Result<Vec<f64>, SdeError>,
    G: FnMut(f64, &[f64]) -> Result<Vec<f64>, SdeError>,
{
    let mut actual = options.clone();
    if let Some(scheme) = scheme {
        actual.method = match scheme {
            | DriftImplicitScheme::Trapezoidal | DriftImplicitScheme::Theta(_) => {
                Method::ThetaMaruyama
            },
            | DriftImplicitScheme::DriftImplicitMilstein => Method::DriftImplicitMilstein,
            | DriftImplicitScheme::BackwardEuler | DriftImplicitScheme::SplitStepBackwardEuler => {
                Method::DriftImplicitEuler
            },
        };
    }
    let options = &actual;
    let noise = options.structure.noise(y0.len());
    let Schedule::Adaptive { t1, cells, .. } = options.schedule else {
        return Err(SdeError::invalid("adaptive schedule required"));
    };
    let areas = uses_levy_areas(options)
        || (options.method == Method::DriftImplicitMilstein && !options.structure.is_commutative());
    let tree = NoiseTree::brownian(
        options.seed,
        noise,
        t1 - options.t0,
        cells,
        areas,
        options.terms,
    )?;
    let drift = |t: f64, y: &[f64], out: &mut [f64]| adapters.drift_at(t, y, out);
    let diffusion = |t: f64, y: &[f64], out: &mut [f64]| adapters.diffusion_at(t, y, out);
    let zero = |_: f64, _: &[f64], out: &mut [f64]| out.fill(0.0);
    let mut control = Control::new(options, y0, 2.0)?;
    let mut t = vec![options.t0];
    let mut y = vec![y0.to_vec()];
    let run = segment(
        &mut control,
        &tree,
        options.t0,
        t1,
        y0,
        &mut t,
        &mut y,
        |t, h, y, block| {
            let result = brownian_step(
                &drift,
                &diffusion,
                noise,
                options.method,
                scheme,
                t,
                h,
                y,
                block,
            );
            if adapters.pending_failure().is_some() {
                return Err(SdeError::invalid("coefficient callback failed"));
            }
            result
        },
        |t, h, y| {
            let block = Block {
                dw: vec![0.0; noise],
                areas: None,
            };
            let result = brownian_step(
                &drift,
                &zero,
                noise,
                options.method,
                scheme,
                t,
                h,
                y,
                &block,
            );
            if adapters.pending_failure().is_some() {
                return Err(SdeError::invalid("coefficient callback failed"));
            }
            result
        },
    );
    if let Some(error) = adapters.failure() {
        return Err(error);
    }
    run?;
    Ok(Path {
        steps: t.len() - 1,
        t,
        y,
        rejected: control.rejected,
        nfev: adapters.nfev(),
        seed: options.seed,
        noise,
        method: options.method,
    })
}

pub(super) fn jump<A, G, J>(
    mut drift: A,
    diffusion: G,
    jump_coeff: J,
    y0: &[f64],
    spec: &JumpSpec,
    options: &Options,
) -> Result<JumpPath, SdeError>
where
    A: FnMut(f64, &[f64]) -> Result<Vec<f64>, SdeError>,
    G: FnMut(f64, &[f64]) -> Result<Vec<f64>, SdeError>,
    J: FnMut(f64, &[f64], &[f64]) -> Result<Vec<f64>, SdeError>,
{
    let mut control = Control::new(options, y0, 2.0)?;
    let Schedule::Adaptive { cells, .. } = options.schedule else {
        return Err(SdeError::invalid("adaptive jump schedule required"));
    };
    let n = y0.len();
    let noise = options.structure.noise(n);
    let areas = uses_levy_areas(options);
    let cells = checked_cells(cells, noise, areas)?;
    let constant = if spec.compensated {
        match &spec.distribution {
            | JumpDistribution::Constant(mark) => Some(mark.clone()),
            | _ => {
                return Err(SdeError::invalid(
                    "nonconstant compensated marks need integrate_jump_diffusion_compensated",
                ));
            },
        }
    } else {
        None
    };
    let mut arrivals = Vec::new();
    let mut stream = Stream::new(ensemble_seed(options.seed, 0x4A55_4D50));
    let mut clock = options.t0;
    let mut mark_values = 0usize;
    if spec.intensity > 0.0 {
        loop {
            let next = clock - stream.next_unit().max(f64::MIN_POSITIVE).ln() / spec.intensity;
            if next >= control.t1 {
                break;
            }
            if !next.is_finite() || next <= clock {
                return Err(SdeError::StepLimit(
                    "Poisson arrival clock collapsed".into(),
                ));
            }
            if arrivals.len() >= control.max_steps {
                return Err(SdeError::StepLimit(
                    "Poisson arrival/mark budget exhausted".into(),
                ));
            }
            let mark = spec.distribution.sample(&mut stream)?;
            mark_values = mark_values
                .checked_add(mark.len())
                .ok_or_else(|| SdeError::invalid("jump mark allocation overflows"))?;
            if mark_values > MAX_NOISE_VALUES || mark.iter().any(|v| !v.is_finite()) {
                return Err(SdeError::invalid(
                    "jump marks exceed the finite coefficient allocation budget",
                ));
            }
            arrivals.push((next, mark));
            clock = next;
        }
    }
    let total_jumps = arrivals.len();
    let jump = RefCell::new(jump_coeff);
    let jump_calls = Cell::new(0usize);
    let adjusted = |t: f64, y: &[f64]| {
        let mut value = drift(t, y)?;
        if value.len() != n {
            return Err(SdeError::shape("invalid jump drift shape"));
        }
        if let Some(mark) = &constant {
            let expected = jump.borrow_mut()(t, y, mark)?;
            jump_calls.set(jump_calls.get() + 1);
            if expected.len() != n || expected.iter().any(|v| !v.is_finite()) {
                return Err(SdeError::shape("invalid compensator output"));
            }
            for (a, c) in value.iter_mut().zip(expected) {
                *a -= spec.intensity * c;
            }
        }
        Ok(value)
    };
    let adapters = Adapters::new(adjusted, diffusion, n, options.structure);
    let drift = |t: f64, y: &[f64], out: &mut [f64]| adapters.drift_at(t, y, out);
    let diffusion = |t: f64, y: &[f64], out: &mut [f64]| adapters.diffusion_at(t, y, out);
    let zero = |_: f64, _: &[f64], out: &mut [f64]| out.fill(0.0);
    let mut times = vec![options.t0];
    let mut states = vec![y0.to_vec()];
    let mut counts = Vec::new();
    let mut start = options.t0;
    for index in 0..=total_jumps {
        let end = arrivals.get(index).map_or(control.t1, |(t, _)| *t);
        let segment_cells =
            ((cells as f64 * (end - start) / (control.t1 - control.t0)).ceil() as usize).max(2);
        let tree = NoiseTree::brownian(
            ensemble_seed(options.seed, index),
            noise,
            end - start,
            segment_cells,
            areas,
            options.terms,
        )?;
        let initial = states
            .last()
            .cloned()
            .ok_or_else(|| SdeError::invalid("empty jump trajectory"))?;
        let run = segment(
            &mut control,
            &tree,
            start,
            end,
            &initial,
            &mut times,
            &mut states,
            |t, h, y, block| {
                let result = brownian_step(
                    &drift,
                    &diffusion,
                    noise,
                    options.method,
                    None,
                    t,
                    h,
                    y,
                    block,
                );
                if adapters.pending_failure().is_some() {
                    return Err(SdeError::invalid("coefficient callback failed"));
                }
                result
            },
            |t, h, y| {
                let block = Block {
                    dw: vec![0.0; noise],
                    areas: None,
                };
                let result =
                    brownian_step(&drift, &zero, noise, options.method, None, t, h, y, &block);
                if adapters.pending_failure().is_some() {
                    return Err(SdeError::invalid("coefficient callback failed"));
                }
                result
            },
        );
        if let Some(e) = adapters.failure() {
            return Err(e);
        }
        run?;
        counts.resize(times.len() - 1, 0);
        if let Some((time, mark)) = arrivals.get(index) {
            let y = states
                .last_mut()
                .ok_or_else(|| SdeError::invalid("empty pre-jump state"))?;
            let delta = jump.borrow_mut()(*time, y, mark)?;
            jump_calls.set(jump_calls.get() + 1);
            if delta.len() != n {
                return Err(SdeError::shape("invalid jump coefficient shape"));
            }
            for (value, change) in y.iter_mut().zip(delta) {
                *value += change;
            }
            if y.iter().any(|v| !v.is_finite()) {
                return Err(SdeError::NonFinite(
                    "post-jump state became non-finite".into(),
                ));
            }
            if let Some(last) = counts.last_mut() {
                *last = 1;
            }
        }
        start = end;
    }
    Ok(JumpPath {
        path: Path {
            steps: times.len() - 1,
            t: times,
            y: states,
            rejected: control.rejected,
            nfev: adapters.nfev() + jump_calls.get(),
            seed: options.seed,
            noise,
            method: options.method,
        },
        total_jumps,
        jump_counts: counts,
    })
}

pub(super) fn levy<A, G>(
    drift: A,
    diffusion: G,
    y0: &[f64],
    alpha: f64,
    beta: f64,
    options: &Options,
) -> Result<Path, SdeError>
where
    A: FnMut(f64, &[f64]) -> Result<Vec<f64>, SdeError>,
    G: FnMut(f64, &[f64]) -> Result<Vec<f64>, SdeError>,
{
    let noise = options.structure.noise(y0.len());
    let Schedule::Adaptive { t1, cells, .. } = options.schedule else {
        return Err(SdeError::invalid("adaptive Levy schedule required"));
    };
    let cells = checked_cells(cells, noise, false)?;
    let increments = levy_increments(
        options.seed,
        noise,
        cells,
        (t1 - options.t0) / cells as f64,
        alpha,
        beta,
    )?;
    let tree = NoiseTree::new(increments, None)?;
    let adapters = Adapters::new(drift, diffusion, y0.len(), options.structure);
    let drift = |t: f64, y: &[f64], out: &mut [f64]| adapters.drift_at(t, y, out);
    let diffusion = |t: f64, y: &[f64], out: &mut [f64]| adapters.diffusion_at(t, y, out);
    let zero = |_: f64, _: &[f64], out: &mut [f64]| out.fill(0.0);
    // Stable increments have no second moment for alpha<2. Use q=0.9 alpha,
    // strictly inside the finite-moment range, rather than a Brownian RMS law.
    let mut control = Control::new(options, y0, 0.9 * alpha)?;
    control.order = (1.0 / alpha).min(1.0);
    let mut times = vec![options.t0];
    let mut states = vec![y0.to_vec()];
    let run = segment(
        &mut control,
        &tree,
        options.t0,
        t1,
        y0,
        &mut times,
        &mut states,
        |t, h, y, block| {
            let result = brownian_step(
                &drift,
                &diffusion,
                noise,
                options.method,
                None,
                t,
                h,
                y,
                block,
            );
            if adapters.pending_failure().is_some() {
                return Err(SdeError::invalid("coefficient callback failed"));
            }
            result
        },
        |t, h, y| {
            let block = Block {
                dw: vec![0.0; noise],
                areas: None,
            };
            let result = brownian_step(&drift, &zero, noise, options.method, None, t, h, y, &block);
            if adapters.pending_failure().is_some() {
                return Err(SdeError::invalid("coefficient callback failed"));
            }
            result
        },
    );
    if let Some(e) = adapters.failure() {
        return Err(e);
    }
    run?;
    Ok(Path {
        steps: times.len() - 1,
        t: times,
        y: states,
        rejected: control.rejected,
        nfev: adapters.nfev(),
        seed: options.seed,
        noise,
        method: options.method,
    })
}

#[allow(clippy::too_many_arguments)] // explicit state, noise and schedule inputs stay in one numerical step contract
pub(super) fn prescribed_step<D, F>(
    drift: &D,
    diffusion: &F,
    noise: usize,
    method: Method,
    t: f64,
    h: f64,
    y: &[f64],
    dw: &[f64],
    areas: Option<&[f64]>,
) -> Result<Vec<f64>, SdeError>
where
    D: Fn(f64, &[f64], &mut [f64]),
    F: Fn(f64, &[f64], &mut [f64]),
{
    brownian_step(
        drift,
        diffusion,
        noise,
        method,
        None,
        t,
        h,
        y,
        &Block {
            dw: dw.to_vec(),
            areas: areas.map(<[f64]>::to_vec),
        },
    )
}
