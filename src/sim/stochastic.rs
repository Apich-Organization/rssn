//! Stochastic differential equations (§8.9.1, new in this round).
//!
//! The differential-equation chapter now has a second surface next to
//! [`crate::kernels::ode`]: `dX = a(t, X) dt + G(t, X) dW` with a typed drift and a
//! diffusion, integrated by the schemes of `crate::kernels::sde`. This crate owns
//! the *interface*: method/calculus/structure selection, the state as a flat
//! `Vec<f64>` (the runtime flattens and rebuilds records and tensors exactly as
//! for the ODE path, §8.9), and one deterministic RNG contract.
//!
//! # What is implemented
//!
//! The **method name selects the scheme family and the calculus; the structure
//! selects the commutative or the non-commutative kernel**:
//!
//! | method | calculus | commutative noise (`Diagonal`, `General{commutative: true}`) | non-commutative noise (`General{commutative: false}`) | order strong/weak |
//! |---|---|---|---|---|
//! | `EulerMaruyama` | Itô | Euler–Maruyama | Euler–Maruyama | 0.5 / 1 |
//! | `Milstein` | Itô | Milstein | Milstein with Lévy areas | 1 / 1 |
//! | `MilsteinStratonovich` | Stratonovich | Milstein (Strat.) | Milstein (Strat.) with Lévy areas | 1 / 1 |
//! | `HeunStratonovich` | Stratonovich | Heun | Heun | 1 / 1 (commutative); **0.5 / 1** (non-commutative) |
//! | `SrkPlaten` | Itô | Platen's derivative-free SRK | Kloeden–Platen SRK (11.1.7) with Lévy areas | 1 / 1 |
//! | `MilsteinNoncommutative` | Itô | Milstein with Lévy areas (explicit alias) | same | 1 / 1 |
//! | `MilsteinNoncommutativeStratonovich` | Stratonovich | Milstein (Strat.) with Lévy areas (alias) | same | 1 / 1 |
//!
//! The calculus is part of the method name, not an orthogonal option: nearly
//! every Itô/Stratonovich mix-up in the wild is a scheme whose name does not say
//! which integral it converges to. [`sde_order_for`] reports the order for a
//! `(method, structure)` pair; [`sde_order`] is the commutative-noise order.
//!
//! Step policy is fixed (`dt`/`steps`, or an explicit **uniform** grid) or
//! adaptive (`rtol`/`atol`, Brownian-bridge rejection in the kernel). Adaptive
//! integration works for every method:
//!
//! * Itô, commutative noise (`EulerMaruyama`, `Milstein`, `SrkPlaten`): the
//!   kernel's embedded EM/Milstein controller, which returns the *Milstein*
//!   solution, so these requests give the same controlled path (the requested
//!   name is kept in [`Path::method`]);
//! * Stratonovich, commutative noise (`MilsteinStratonovich`,
//!   `HeunStratonovich`): **exact conversion to Itô form** — the drift becomes
//!   `a~ = a + 1/2 sum_j (D_j g)_{.j}` and the Itô controller runs on `(a~, g)`.
//!   The returned path is the controlled Itô–Milstein solution of the equivalent
//!   Itô SDE, which converges to the Stratonovich solution;
//! * non-commutative noise (and the `*Noncommutative*` methods): the fine-cell
//!   controller with Chen's relation (`cells`, default [`DEFAULT_CELLS`]), the
//!   Itô one, or the Itô-converted one for Stratonovich methods.
//!
//! # Determinism and the two RNG streams
//!
//! * Fixed-step paths on the **language stream** ([`wiener`], [`Stream`]): the
//!   stream is the interpreter's `rng(seed)` LCG (state → `seed * 6364136223846793005 + 1`,
//!   `(state >> 11) / 2^53`), with standard normals by Box–Muller (two uniforms
//!   per normal, cosine branch). [`wiener`] and a fixed-step diagonal/general path
//!   therefore consume *the same numbers* — a test asserts the path equals a
//!   hand-driven kernel EM run on `wiener`'s increments. Same build ⇒
//!   bit-identical; cross-platform bit identity is not promised.
//! * Fixed-step schemes that use Lévy areas (every method on non-commutative
//!   noise except Euler–Maruyama and Heun, and the `*Noncommutative*` methods)
//!   and every adaptive path need the
//!   kernel's internal xoshiro stream (`crate::kernels::random::Rng`) because the
//!   kernel samples Brownian bridges and Lévy areas as it goes. They are
//!   seeded from [`Options::seed`] and deterministic for a given build; they do
//!   not share the language LCG stream. This is a documented two-stream
//!   reality: the language stream is the one `wiener`/`rng` expose, the kernel
//!   stream is what bridge splitting consumes.
//!
//! # Limitations (deliberate, not silent)
//!
//! * Fixed-step kernels assume a uniform step; an explicit non-uniform grid is
//!   rejected with `E0802` (use the adaptive schedule for variable steps).
//! * A `General { commutative: true }` structure is a **claim by the caller**:
//!   the Milstein/SRK corrections without Lévy areas are only valid for
//!   commutative noise. Declare `commutative: false` and the same method names
//!   dispatch to the Lévy-area kernels instead (no error).
//! * Lévy areas are simulated with a finite Fourier expansion
//!   ([`Options::terms`]) plus the Gaussian tail correction; accuracy improves
//!   with the mode count. `terms == 0` means "default": [`default_terms`]
//!   `= clamp(ceil(1/sqrt(h)), 16, 4096)` for the step (or the Brownian cell) `h`.
//!   Non-commutative schemes are strong order one at best.
//! * `HeunStratonovich` on non-commutative noise is accepted but is only strong
//!   order 1/2 there (it uses no Lévy areas); use `MilsteinStratonovich` or a
//!   `*Noncommutative*` method for order one.
//! * No events/root-finding and no dense output in this round (`SdeError::Unsupported`,
//!   capability `sde.events`, is reserved for them; no current combination produces it).
//! * Ensembles run paths sequentially: an Aevia interpreter closure is
//!   single-threaded (`!Send`), so the kernel's rayon ensemble cannot be used
//!   without an `unsafe` callback. The per-path seed is derived from the pair
//!   `(seed, index)` only, so the result is independent of execution order.
//! * Record/tensor states are handled by the caller: this module sees the flat
//!   slice produced by the same flatten/rebuild contract as [`crate::kernels::ode`].
//!
//! The statistical order/convergence tests live in the kernel's own suite; the
//! tests here are deterministic (pinned seeds, path equality, error codes).

use crate::kernels::random::Rng;
use crate::kernels::sde as kernel;
use std::cell::Cell;
use std::cell::RefCell;

// ---------------------------------------------------------------------------
// Errors
// ---------------------------------------------------------------------------

/// Why an SDE integration refused or failed.
///
/// The code is the stable Aevia diagnostic (appendix B): `E0802` for a shape or
/// parameter error, `E0806` for a value outside the domain (a non-finite state),
/// `E0807` for a step collapse, and `E1001` for a method/structure/step-policy
/// combination this build does not implement — naming the capability.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SdeError {
    /// Inconsistent dimensions or parameters (`E0802`).
    Invalid(String),
    /// A shape (`E0802`); kept separate from [`Self::Invalid`] for messages.
    Shape(String),
    /// The state became non-finite (`E0806`), with the kernel's reason.
    NonFinite(String),
    /// The adaptive step collapsed or exhausted its budget (`E0807`).
    StepLimit(String),
    /// An unsupported combination, naming the capability (`E1001`).
    Unsupported {
        /// The capability a build would need (`core.sde`, `sde.events`, …). The
        /// former `sde.noncommutative` and `sde.interpretation` gaps are closed.
        capability: &'static str,
        /// What exactly is not implemented, for the note.
        detail: String,
    },
}

impl SdeError {
    /// The stable Aevia diagnostic code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            | Self::Invalid(_) | Self::Shape(_) => "E0802",
            | Self::NonFinite(_) => "E0806",
            | Self::StepLimit(_) => "E0807",
            | Self::Unsupported { .. } => "E1001",
        }
    }

    /// The capability an `E1001` names, when it is one.
    #[must_use]
    pub const fn capability(&self) -> Option<&'static str> {
        match self {
            | Self::Unsupported { capability, .. } => Some(capability),
            | Self::Invalid(_) | Self::Shape(_) | Self::NonFinite(_) | Self::StepLimit(_) => None,
        }
    }

    fn shape(message: impl Into<String>) -> Self {
        Self::Shape(message.into())
    }

    fn invalid(message: impl Into<String>) -> Self {
        Self::Invalid(message.into())
    }

    // Reserved for capabilities that are genuinely missing (`sde.events`); no current
    // method/structure/step-policy combination produces it.
    #[allow(dead_code)]
    fn unsupported(
        capability: &'static str,
        detail: impl Into<String>,
    ) -> Self {
        Self::Unsupported {
            capability,
            detail: detail.into(),
        }
    }
}

impl std::fmt::Display for SdeError {
    fn fmt(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        match self {
            | Self::Invalid(message) | Self::Shape(message) => f.write_str(message),
            | Self::NonFinite(reason) => write!(f, "the state became non-finite: {reason}"),
            | Self::StepLimit(reason) => f.write_str(reason),
            | Self::Unsupported { capability, detail } => {
                write!(f, "this build does not provide `{capability}`: {detail}")
            },
        }
    }
}

impl std::error::Error for SdeError {}

impl From<kernel::SdeError> for SdeError {
    fn from(error: kernel::SdeError) -> Self {
        match error {
            | kernel::SdeError::Invalid { what } => {
                Self::invalid(format!("the kernel rejected the SDE parameters: {what}"))
            },
            | kernel::SdeError::NonFinite { step } => {
                Self::NonFinite(format!("the state left the reals at step {step}"))
            },
            | kernel::SdeError::StepLimit { limit } => {
                Self::StepLimit(format!("the step budget of {limit} was exhausted"))
            },
            | kernel::SdeError::StepUnderflow { t, h } => {
                Self::StepLimit(format!("the step size underflowed at t = {t} (h = {h:e})"))
            },
            | kernel::SdeError::StepBelowCell { cells } => {
                Self::StepLimit(format!(
                    "an adaptive step cannot be shorter than one of the {cells} Brownian cells; \
                 use more cells"
                ))
            },
        }
    }
}

// ---------------------------------------------------------------------------
// Methods, structures, options
// ---------------------------------------------------------------------------

/// A stochastic integration scheme. The calculus is part of the name.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Method {
    /// Euler–Maruyama (Itô), strong order 1/2.
    EulerMaruyama,
    /// Milstein (Itô), strong order 1 for commutative noise.
    Milstein,
    /// Milstein (Stratonovich), strong order 1 for commutative noise.
    MilsteinStratonovich,
    /// Heun predictor–corrector (Stratonovich), strong order 1 for commutative
    /// noise.
    HeunStratonovich,
    /// SRK–Platen (Itô), strong order 1 for commutative noise.
    SrkPlaten,
    /// Milstein (Itô) with simulated Lévy areas, whatever the structure (an
    /// explicit alias; `Milstein` already dispatches here on non-commutative
    /// noise).
    MilsteinNoncommutative,
    /// Milstein (Stratonovich) with simulated Lévy areas, whatever the structure
    /// (an explicit alias of `MilsteinStratonovich` on non-commutative noise).
    MilsteinNoncommutativeStratonovich,
    /// Drift-implicit Euler–Maruyama (backward Euler on drift, explicit on diffusion;
    /// L-stable in the drift for stiff SDEs), strong order 1/2.
    DriftImplicitEuler,
    /// Symmetric drift-implicit trapezoidal / theta-Maruyama (`theta = 0.5`, Itô),
    /// A-stable in the drift, strong order 1/2.
    ThetaMaruyama,
    /// Drift-implicit Milstein (backward Euler on drift + Milstein diffusion correction,
    /// Itô), L-stable in the drift and strong order 1 for commutative noise.
    DriftImplicitMilstein,
}

/// Which stochastic integral a scheme converges to.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Calculus {
    /// The Itô integral.
    Ito,
    /// The Stratonovich integral.
    Stratonovich,
}

impl Method {
    /// The spelling used by the language surface and diagnostics.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            | Self::EulerMaruyama => "euler_maruyama",
            | Self::Milstein => "milstein",
            | Self::MilsteinStratonovich => "milstein_stratonovich",
            | Self::HeunStratonovich => "heun_stratonovich",
            | Self::SrkPlaten => "srk_platen",
            | Self::MilsteinNoncommutative => "milstein_noncommutative",
            | Self::MilsteinNoncommutativeStratonovich => "milstein_noncommutative_stratonovich",
            | Self::DriftImplicitEuler => "drift_implicit_euler",
            | Self::ThetaMaruyama => "theta_maruyama",
            | Self::DriftImplicitMilstein => "drift_implicit_milstein",
        }
    }

    /// The method named `name`, if it is one this build knows.
    #[must_use]
    pub fn from_name(name: &str) -> Option<Self> {
        Some(match name {
            | "euler_maruyama" | "em" => Self::EulerMaruyama,
            | "milstein" => Self::Milstein,
            | "milstein_stratonovich" => Self::MilsteinStratonovich,
            | "heun_stratonovich" | "heun" => Self::HeunStratonovich,
            | "srk_platen" | "srk" => Self::SrkPlaten,
            | "milstein_noncommutative" => Self::MilsteinNoncommutative,
            | "milstein_noncommutative_stratonovich" => Self::MilsteinNoncommutativeStratonovich,
            | "drift_implicit_euler" | "implicit_euler" => Self::DriftImplicitEuler,
            | "theta_maruyama" | "theta" => Self::ThetaMaruyama,
            | "drift_implicit_milstein" | "implicit_milstein" => Self::DriftImplicitMilstein,
            | _ => return None,
        })
    }

    /// The integral the scheme converges to.
    #[must_use]
    pub const fn calculus(self) -> Calculus {
        match self {
            | Self::EulerMaruyama
            | Self::Milstein
            | Self::SrkPlaten
            | Self::MilsteinNoncommutative
            | Self::DriftImplicitEuler
            | Self::ThetaMaruyama
            | Self::DriftImplicitMilstein => Calculus::Ito,
            | Self::MilsteinStratonovich
            | Self::HeunStratonovich
            | Self::MilsteinNoncommutativeStratonovich => Calculus::Stratonovich,
        }
    }
}

/// The order of a scheme, as documented by the kernel.
///
/// The kernel states that *no* scheme here is better than strong order one and
/// that its Milstein variants are weak order one (not two), so weak order two
/// is never claimed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Order {
    /// Strong (pathwise) order.
    pub strong: f64,
    /// Weak (moment) order.
    pub weak: f64,
}

/// The strong/weak order of `method` on **commutative** noise (see
/// [`sde_order_for`] for a given structure), as the kernel's documentation states it.
#[must_use]
pub const fn sde_order(method: Method) -> Order {
    match method {
        // Euler–Maruyama and drift-implicit Euler / theta-Maruyama are strong 1/2, weak 1.
        | Method::EulerMaruyama | Method::DriftImplicitEuler | Method::ThetaMaruyama => {
            Order {
                strong: 0.5,
                weak: 1.0,
            }
        },
        // Every other implemented scheme is strong order one. The kernel is
        // explicit that they are weak order one too (not two).
        | Method::Milstein
        | Method::MilsteinStratonovich
        | Method::HeunStratonovich
        | Method::SrkPlaten
        | Method::MilsteinNoncommutative
        | Method::MilsteinNoncommutativeStratonovich
        | Method::DriftImplicitMilstein => {
            Order {
                strong: 1.0,
                weak: 1.0,
            }
        },
    }
}

/// The strong/weak order of `method` on noise of the given `structure`.
///
/// As [`sde_order`], except that `HeunStratonovich` on non-commutative noise is
/// only strong order 1/2 (it has no Lévy-area term), and that the schemes which
/// dispatch to the Lévy-area kernels on non-commutative noise keep strong order
/// one.
#[must_use]
pub fn sde_order_for(
    method: Method,
    structure: Structure,
) -> Order {
    if method == Method::HeunStratonovich && !structure.is_commutative() {
        return Order {
            strong: 0.5,
            weak: 1.0,
        };
    }
    sde_order(method)
}

/// How the diffusion enters the equation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Structure {
    /// One noise channel per state component: `g` is `n × n` and diagonal, and
    /// the diffusion callback returns its `n` diagonal entries.
    Diagonal,
    /// `m` noise channels; the diffusion callback returns the `n × m` matrix in
    /// row-major order. `commutative` is the caller's claim: the derivative
    /// corrections used by Milstein/SRK/Heun are only valid for commutative
    /// noise.
    General {
        /// Number of independent Wiener processes (`m`).
        noise: usize,
        /// Whether the caller guarantees the noise is commutative.
        commutative: bool,
    },
}

impl Structure {
    /// The number of noise channels (`m`) for a state of `n` components.
    #[must_use]
    pub const fn noise(
        self,
        state: usize,
    ) -> usize {
        match self {
            | Self::Diagonal => state,
            | Self::General { noise, .. } => noise,
        }
    }

    /// Whether the structure is (declared) commutative.
    #[must_use]
    pub const fn is_commutative(self) -> bool {
        match self {
            | Self::Diagonal => true,
            | Self::General { commutative, .. } => commutative,
        }
    }

    /// How many values the diffusion callback returns.
    #[must_use]
    pub const fn coefficients(
        self,
        state: usize,
    ) -> usize {
        match self {
            | Self::Diagonal => state,
            | Self::General { noise, .. } => state * noise,
        }
    }
}

/// The step policy.
#[derive(Clone, Debug, PartialEq)]
pub enum Schedule {
    /// Fixed steps of `dt`, `steps` of them: `t_k = t0 + k dt`.
    Fixed {
        /// Step size.
        dt: f64,
        /// Number of steps.
        steps: usize,
    },
    /// An explicit, **uniform** time grid (the first entry is `t0`).
    Grid {
        /// The times, strictly increasing with equal spacing.
        times: Vec<f64>,
    },
    /// Adaptive steps on `[t0, t1]`.
    Adaptive {
        /// End time.
        t1: f64,
        /// Relative tolerance of the local error.
        rtol: f64,
        /// Absolute tolerance of the local error.
        atol: f64,
        /// Step budget; exhaustion is `E0807`.
        max_steps: usize,
        /// Brownian cells for the non-commutative adaptive solver on `[t0,t1]`.
        /// Ignored by the commutative solvers; `0` selects [`DEFAULT_CELLS`] when
        /// the non-commutative solver is used.
        cells: usize,
    },
}

/// Options of one integration.
///
/// There is no `Default`: a seed and a step are not guessable, so they are
/// constructor arguments (the spec requires an explicit seed for reproducibility
/// and an explicit step policy). [`Options::with_terms`] is optional: `0` selects
/// [`default_terms`].
#[derive(Clone, Debug, PartialEq)]
pub struct Options {
    /// The scheme.
    pub method: Method,
    /// The noise structure.
    pub structure: Structure,
    /// The seed of the Brownian path (and of the language stream for fixed
    /// steps).
    pub seed: u64,
    /// Initial time.
    pub t0: f64,
    /// The step policy.
    pub schedule: Schedule,
    /// Fourier modes of the simulated Lévy areas; `0` (the default) selects
    /// [`default_terms`] for the step (or Brownian cell) size.
    pub terms: usize,
}

impl Options {
    /// `steps` fixed steps of `dt` from `t0`.
    #[must_use]
    pub const fn fixed(
        method: Method,
        structure: Structure,
        seed: u64,
        t0: f64,
        dt: f64,
        steps: usize,
    ) -> Self {
        Self {
            method,
            structure,
            seed,
            t0,
            schedule: Schedule::Fixed { dt, steps },
            terms: 0,
        }
    }

    /// Fixed steps covering `[t0, t1]`: the step count is rounded, and `dt` is
    /// adjusted to `(t1 - t0) / steps` so the end time is reached exactly.
    #[must_use]
    pub fn fixed_span(
        method: Method,
        structure: Structure,
        seed: u64,
        t0: f64,
        t1: f64,
        dt: f64,
    ) -> Self {
        let steps = if dt > 0.0 && t1 > t0 {
            ((t1 - t0) / dt).round().max(1.0) as usize
        } else {
            1
        };
        let adjusted = (t1 - t0) / steps as f64;
        Self::fixed(method, structure, seed, t0, adjusted, steps)
    }

    /// An explicit uniform grid.
    #[must_use]
    pub fn grid(
        method: Method,
        structure: Structure,
        seed: u64,
        times: Vec<f64>,
    ) -> Self {
        let t0 = times.first().copied().unwrap_or(0.0);
        Self {
            method,
            structure,
            seed,
            t0,
            schedule: Schedule::Grid { times },
            terms: 0,
        }
    }

    /// Adaptive integration on `[t0, t1]` (`cells` is used only by the
    /// non-commutative method; pass `0` otherwise).
    #[must_use]
    #[allow(clippy::too_many_arguments)] // the option record is the language's, not a builder
    pub const fn adaptive(
        method: Method,
        structure: Structure,
        seed: u64,
        t0: f64,
        t1: f64,
        rtol: f64,
        atol: f64,
        max_steps: usize,
        cells: usize,
    ) -> Self {
        Self {
            method,
            structure,
            seed,
            t0,
            schedule: Schedule::Adaptive {
                t1,
                rtol,
                atol,
                max_steps,
                cells,
            },
            terms: 0,
        }
    }

    /// Sets the number of Fourier modes for the simulated Lévy areas.
    #[must_use]
    pub const fn with_terms(
        mut self,
        terms: usize,
    ) -> Self {
        self.terms = terms;
        self
    }
}

/// Brownian cells of the non-commutative adaptive solver when `cells == 0`.
pub const DEFAULT_CELLS: usize = 4096;

/// The default number of Fourier modes of the simulated Lévy areas for a step
/// (or Brownian cell) of size `h`: `clamp(ceil(1 / sqrt(h)), 16, 4096)`.
///
/// The truncation error of the area decays like `h p^(-3/2)` with the kernel's
/// Gaussian tail correction, so this rule keeps it far below the `O(h)` strong
/// error of the scheme at an `O(m^2 / sqrt(h))` sampling cost per step.
#[must_use]
pub fn default_terms(h: f64) -> usize {
    let h = h.abs();
    if !(h.is_finite() && h > 0.0) {
        return 16;
    }
    let p = (1.0 / h.sqrt()).ceil();
    if p >= 4096.0 {
        4096
    } else {
        (p as usize).max(16)
    }
}

// ---------------------------------------------------------------------------
// The language RNG stream
// ---------------------------------------------------------------------------

/// The interpreter's `rng(seed)` stream, re-implemented here so the standard
/// library can produce the same numbers the language does.
///
/// The formula matches the Aevia interpreter's `next_random` contract: a
/// 64-bit LCG with multiplier `6364136223846793005` and increment `1`, taking the
/// top 53 bits. [`Stream::new`] mirrors the interpreter's default seed only if
/// given the same value; a language program's `rng(seed)` and
/// `sdeint(..., seed)` share the sequence when they share the seed.
///
/// Standard normals are Box–Muller, cosine branch, two uniforms per normal. A
/// zero uniform is replaced by the smallest positive double so `ln` stays
/// finite; the probability of that is below `2^-53`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Stream {
    state: u64,
}

impl Stream {
    /// A stream seeded with `seed` (the same state `rng(seed)` starts from).
    #[must_use]
    pub const fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    /// The next uniform in `[0, 1)`.
    pub fn next_unit(&mut self) -> f64 {
        self.state = self
            .state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1);
        ((self.state >> 11) as f64) / 9_007_199_254_740_992.0
    }

    /// The next standard normal, Box–Muller (cosine branch).
    pub fn next_normal(&mut self) -> f64 {
        let u1 = self.next_unit().max(f64::MIN_POSITIVE);
        let u2 = self.next_unit();
        (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
    }
}

/// Brownian increments for `steps` steps of `dt`, `m` components each, on the
/// language stream seeded with `seed`.
///
/// The argument order is the kernel's `sde_wiener(seed, m, steps, h)` so the
/// language surface, the manifest operator and this helper cannot drift:
/// `wiener(seed, m, steps, dt)[k][j]` is the increment of the `j`-th Wiener
/// process over `[t0 + k dt, t0 + (k+1) dt]`, drawn step-major (component order
/// inside a step), which is exactly the layout the kernel's fixed-step solvers
/// consume. A fixed-step integration with the same seed and structure therefore
/// consumes these numbers.
#[must_use]
pub fn wiener(
    seed: u64,
    m: usize,
    steps: usize,
    dt: f64,
) -> Vec<Vec<f64>> {
    let mut stream = Stream::new(seed);
    wiener_from(&mut stream, m, dt, steps)
}

fn wiener_from(
    stream: &mut Stream,
    m: usize,
    dt: f64,
    steps: usize,
) -> Vec<Vec<f64>> {
    let scale = dt.abs().sqrt();
    (0..steps)
        .map(|_| {
            (0..m)
                .map(|_| scale * stream.next_normal())
                .collect::<Vec<f64>>()
        })
        .collect()
}

/// Seed of path `index` of an ensemble with base `seed`.
///
/// A function of the pair only (the kernel's derivation), so path `i` is the
/// same whichever order or concurrency runs the ensemble; this build runs the
/// paths sequentially because an interpreter closure is `!Send`.
#[must_use]
pub const fn ensemble_seed(
    seed: u64,
    index: usize,
) -> u64 {
    kernel::ensemble_seed(seed, index)
}

// ---------------------------------------------------------------------------
// Results
// ---------------------------------------------------------------------------

/// One integrated path.
#[derive(Clone, Debug, PartialEq)]
pub struct Path {
    /// Accepted time points (`t` and `y` are parallel).
    pub t: Vec<f64>,
    /// States at the accepted times, in time order.
    pub y: Vec<Vec<f64>>,
    /// Accepted steps.
    pub steps: usize,
    /// Steps rejected by the adaptive controller (fixed-step paths: 0).
    pub rejected: usize,
    /// Number of drift + diffusion callbacks made.
    pub nfev: usize,
    /// The seed actually used.
    pub seed: u64,
    /// The number of noise channels (`m`).
    pub noise: usize,
    /// The method actually used (adaptive `EulerMaruyama` runs the kernel's
    /// controlled Milstein solution, but keeps the requested name).
    pub method: Method,
}

impl Path {
    /// The final state, when the path has one.
    #[must_use]
    pub fn last(&self) -> Option<&Vec<f64>> {
        self.y.last()
    }
}

// ---------------------------------------------------------------------------
// Callback adapters
// ---------------------------------------------------------------------------

/// Adapts the Aevia-side `FnMut` callbacks to the kernel's `Fn` callbacks.
///
/// The kernel takes `Fn` (it may call a closure inside a rayon ensemble, where
/// the closure must be `Sync`); an Aevia interpreter closure is `!Send` and
/// cannot cross. This adapter owns the closures in a `RefCell`, writes into the
/// kernel's output slice, and records the first failure — because the kernel's
/// callbacks return `()`, an error has to travel through the adapter, not the
/// return value. No `unsafe`: the interior mutability is the whole mechanism.
struct Adapters<A, G> {
    drift: RefCell<A>,
    diffusion: RefCell<G>,
    drift_error: RefCell<Option<SdeError>>,
    diffusion_error: RefCell<Option<SdeError>>,
    state: usize,
    coefficients: usize,
    diagonal: bool,
    calls: Cell<usize>,
}

impl<A, G> Adapters<A, G>
where
    A: FnMut(f64, &[f64]) -> Result<Vec<f64>, SdeError>,
    G: FnMut(f64, &[f64]) -> Result<Vec<f64>, SdeError>,
{
    fn new(
        drift: A,
        diffusion: G,
        state: usize,
        structure: Structure,
    ) -> Self {
        Self {
            drift: RefCell::new(drift),
            diffusion: RefCell::new(diffusion),
            drift_error: RefCell::new(None),
            diffusion_error: RefCell::new(None),
            state,
            coefficients: structure.coefficients(state),
            diagonal: structure == Structure::Diagonal,
            calls: Cell::new(0),
        }
    }

    /// Writes the drift into the kernel's output slice.
    fn drift_at(
        &self,
        t: f64,
        y: &[f64],
        out: &mut [f64],
    ) {
        if self.drift_error.borrow().is_some() {
            out.fill(0.0);
            return;
        }
        self.calls.set(self.calls.get() + 1);
        let result = {
            let mut drift = self.drift.borrow_mut();
            (*drift)(t, y)
        };
        match result {
            | Ok(values) if values.len() == self.state => out.copy_from_slice(&values),
            | Ok(values) => {
                *self.drift_error.borrow_mut() = Some(SdeError::shape(format!(
                    "the drift returned {} value(s) for a state of {}",
                    values.len(),
                    self.state
                )));
                out.fill(0.0);
            },
            | Err(error) => {
                *self.drift_error.borrow_mut() = Some(error);
                out.fill(0.0);
            },
        }
    }

    /// Writes the diffusion into the kernel's `n × m` output slice.
    fn diffusion_at(
        &self,
        t: f64,
        y: &[f64],
        out: &mut [f64],
    ) {
        if self.diffusion_error.borrow().is_some() {
            out.fill(0.0);
            return;
        }
        self.calls.set(self.calls.get() + 1);
        let result = {
            let mut diffusion = self.diffusion.borrow_mut();
            (*diffusion)(t, y)
        };
        match result {
            | Ok(values) if values.len() == self.coefficients => {
                if self.diagonal {
                    // Expand the diagonal to the kernel's dense n×n layout.
                    out.fill(0.0);
                    for index in 0..self.state {
                        out[index * self.state + index] = values[index];
                    }
                } else {
                    out.copy_from_slice(&values);
                }
            },
            | Ok(values) => {
                *self.diffusion_error.borrow_mut() = Some(SdeError::shape(format!(
                    "the diffusion returned {} coefficient(s); this structure needs {}",
                    values.len(),
                    self.coefficients
                )));
                out.fill(0.0);
            },
            | Err(error) => {
                *self.diffusion_error.borrow_mut() = Some(error);
                out.fill(0.0);
            },
        }
    }

    /// The first callback failure, if any.
    fn failure(&self) -> Option<SdeError> {
        self.drift_error
            .borrow_mut()
            .take()
            .or_else(|| self.diffusion_error.borrow_mut().take())
    }

    const fn nfev(&self) -> usize {
        self.calls.get()
    }
}

// ---------------------------------------------------------------------------
// Integration
// ---------------------------------------------------------------------------

/// Integrates `dX = a(t, X) dt + G(t, X) dW` from `y0`.
///
/// `drift` returns `n` values; `diffusion` returns `n` values for
/// [`Structure::Diagonal`] (the diagonal entries) or `n × m` values in
/// row-major order for [`Structure::General`].
///
/// The step policy comes from [`Options::schedule`]. Fixed steps consume
/// [`wiener`]'s language-stream increments; the non-commutative method and every
/// adaptive path consume the kernel's internal stream (seeded from
/// [`Options::seed`]); see the module docs.
///
/// The adaptive controller integrates with the kernel's Milstein-controlled
/// step size: the returned states are the kernel's Milstein solution for the
/// `EulerMaruyama`, `Milstein` and `SrkPlaten` requests (Itô). For the
/// Stratonovich methods the drift is converted exactly to its Itô equivalent
/// (`a~ = a + 1/2 sum_j (D_j g)_{.j}`) and the same controller runs on the
/// converted SDE: the returned path is the controlled Itô–Milstein solution of
/// the equivalent Itô SDE, which converges to the Stratonovich solution.
/// Non-commutative noise uses the fine-cell controller (`cells`, default
/// [`DEFAULT_CELLS`]). This is the kernel's documented contract, not a hidden
/// substitution of the requested scheme's control law.
///
/// # Errors
///
/// `E0802` for shapes/parameters, `E0806` for a non-finite state, `E0807` for a
/// step collapse. No method/structure/step-policy combination is refused with
/// `E1001` any more.
pub fn integrate<A, G>(
    drift: A,
    diffusion: G,
    y0: &[f64],
    options: &Options,
) -> Result<Path, SdeError>
where
    A: FnMut(f64, &[f64]) -> Result<Vec<f64>, SdeError>,
    G: FnMut(f64, &[f64]) -> Result<Vec<f64>, SdeError>,
{
    validate_options(options, y0)?;
    if matches!(
        options.method,
        Method::DriftImplicitEuler | Method::ThetaMaruyama | Method::DriftImplicitMilstein
    ) && matches!(options.schedule, Schedule::Adaptive { .. })
    {
        return Err(SdeError::invalid(
            "drift-implicit schemes require a fixed or uniform-grid schedule",
        ));
    }
    let state = y0.len();
    let noise = options.structure.noise(state);
    let adapters = Adapters::new(drift, diffusion, state, options.structure);
    let run = match &options.schedule {
        | Schedule::Fixed { dt, steps } => fixed(&adapters, y0, options, *dt, *steps, noise),
        | Schedule::Grid { times } => {
            let (dt, steps) = uniform_grid(times)?;
            if (times[0] - options.t0).abs() > 1e-12 * options.t0.abs().max(1.0) {
                return Err(SdeError::invalid(
                    "the explicit grid does not start at the initial time",
                ));
            }
            fixed(&adapters, y0, options, dt, steps, noise)
        },
        | Schedule::Adaptive {
            t1,
            rtol,
            atol,
            max_steps,
            cells,
        } => {
            adaptive(
                &adapters, y0, options, *t1, *rtol, *atol, *max_steps, *cells,
            )
        },
    };
    // A callback failure is recorded inside the adapter (the kernel callbacks
    // return `()`); it takes precedence over the kernel's own error code because
    // it is the cause, not the symptom.
    if let Some(failure) = adapters.failure() {
        return Err(failure);
    }
    let mut path = run?;
    path.nfev = adapters.nfev();
    path.seed = options.seed;
    path.noise = noise;
    path.method = options.method;
    Ok(path)
}

/// Integrates one path per ensemble member, sequentially and order-preserving.
///
/// Member `i` uses [`ensemble_seed`] with `(seed, i)` as its own seed, so the result
/// does not depend on execution order or concurrency. The paths run one after
/// another because the Aevia interpreter closures are `!Send`; the kernel's
/// rayon ensemble would need a `Send + Sync` callback, which this build cannot
/// produce. The returned paths are in index order.
///
/// # Errors
///
/// The first member error, in index order.
pub fn ensemble<A, G>(
    mut drift: A,
    mut diffusion: G,
    y0: &[f64],
    options: &Options,
    paths: usize,
) -> Result<Vec<Path>, SdeError>
where
    A: FnMut(f64, &[f64]) -> Result<Vec<f64>, SdeError>,
    G: FnMut(f64, &[f64]) -> Result<Vec<f64>, SdeError>,
{
    if paths == 0 {
        return Err(SdeError::invalid("an ensemble needs at least one path"));
    }
    let mut out = Vec::with_capacity(paths);
    for index in 0..paths {
        let mut member = options.clone();
        member.seed = ensemble_seed(options.seed, index);
        out.push(integrate(&mut drift, &mut diffusion, y0, &member)?);
    }
    Ok(out)
}

/// Checks the options against a state length before any work is done.
fn validate_options(
    options: &Options,
    y0: &[f64],
) -> Result<(), SdeError> {
    if y0.is_empty() {
        return Err(SdeError::shape("the initial state is empty"));
    }
    if !y0.iter().all(|value| value.is_finite()) {
        return Err(SdeError::NonFinite(
            "the initial state is not finite".to_string(),
        ));
    }
    if !options.t0.is_finite() {
        return Err(SdeError::invalid("the initial time is not finite"));
    }
    if options.structure.noise(y0.len()) == 0 {
        return Err(SdeError::invalid("the noise dimension is zero"));
    }
    match &options.schedule {
        | Schedule::Fixed { dt, steps } => {
            if !dt.is_finite() || *dt <= 0.0 {
                return Err(SdeError::invalid(
                    "the step size must be finite and positive",
                ));
            }
            if *steps == 0 {
                return Err(SdeError::invalid("the step count must be positive"));
            }
        },
        | Schedule::Grid { times } => {
            uniform_grid(times)?;
        },
        | Schedule::Adaptive {
            t1,
            rtol,
            atol,
            max_steps,
            cells,
        } => {
            if !t1.is_finite() || *t1 <= options.t0 {
                return Err(SdeError::invalid(
                    "the end time must be after the initial time",
                ));
            }
            if !rtol.is_finite() || *rtol <= 0.0 {
                return Err(SdeError::invalid("`rtol` must be finite and positive"));
            }
            if !atol.is_finite() || *atol < 0.0 {
                return Err(SdeError::invalid("`atol` must be finite and non-negative"));
            }
            if *max_steps == 0 {
                return Err(SdeError::invalid("the step budget must be positive"));
            }
            let _ = cells;
        },
    }
    Ok(())
}

/// The step size and count of a uniform explicit grid.
fn uniform_grid(times: &[f64]) -> Result<(f64, usize), SdeError> {
    if times.len() < 2 {
        return Err(SdeError::invalid("the time grid needs at least two points"));
    }
    if !times.iter().all(|time| time.is_finite()) {
        return Err(SdeError::invalid(
            "the time grid contains a non-finite time",
        ));
    }
    let dt = times[1] - times[0];
    if dt <= 0.0 {
        return Err(SdeError::invalid(
            "the time grid must be strictly increasing",
        ));
    }
    let tolerance = 1e-9 * dt.abs().max(1.0);
    for window in times.windows(2) {
        if (window[1] - window[0] - dt).abs() > tolerance {
            return Err(SdeError::invalid(
                "a fixed-step path needs a uniform grid; use the adaptive \
                 schedule for variable steps",
            ));
        }
    }
    Ok((dt, times.len() - 1))
}

/// Whether a fixed-step run needs the kernel's Lévy-area noise: the explicit
/// `*Noncommutative*` methods always, and `Milstein`, `MilsteinStratonovich` and
/// `SrkPlaten` on non-commutative noise. Euler–Maruyama and Heun use plain
/// increments (and Heun is then only strong order 1/2).
const fn uses_levy_areas(options: &Options) -> bool {
    match options.method {
        | Method::MilsteinNoncommutative | Method::MilsteinNoncommutativeStratonovich => true,
        | Method::Milstein | Method::MilsteinStratonovich | Method::SrkPlaten => {
            !options.structure.is_commutative()
        },
        | Method::EulerMaruyama
        | Method::HeunStratonovich
        | Method::DriftImplicitEuler
        | Method::ThetaMaruyama
        | Method::DriftImplicitMilstein => false,
    }
}

/// Whether an adaptive run uses the non-commutative (fine-cell) controller.
const fn adaptive_noncommutative(options: &Options) -> bool {
    !options.structure.is_commutative()
        || matches!(
            options.method,
            Method::MilsteinNoncommutative | Method::MilsteinNoncommutativeStratonovich
        )
}

/// A fixed-step integration.
fn fixed<A, G>(
    adapters: &Adapters<A, G>,
    y0: &[f64],
    options: &Options,
    dt: f64,
    steps: usize,
    noise: usize,
) -> Result<Path, SdeError>
where
    A: FnMut(f64, &[f64]) -> Result<Vec<f64>, SdeError>,
    G: FnMut(f64, &[f64]) -> Result<Vec<f64>, SdeError>,
{
    let t0 = options.t0;
    let drift = |t: f64, y: &[f64], out: &mut [f64]| adapters.drift_at(t, y, out);
    let diffusion = |t: f64, y: &[f64], out: &mut [f64]| adapters.diffusion_at(t, y, out);
    let states = if uses_levy_areas(options) {
        // Lévy areas are sampled by the kernel's stream (it needs the same Rng
        // that draws the bridge/fourier coefficients).
        let mut rng = Rng::new(options.seed);
        let terms = if options.terms == 0 {
            default_terms(dt)
        } else {
            options.terms
        };
        let sample = kernel::sample_noise(&mut rng, noise, steps, dt, terms);
        match options.method {
            | Method::MilsteinStratonovich | Method::MilsteinNoncommutativeStratonovich => {
                kernel::milstein_noncommutative_stratonovich(
                    drift, diffusion, noise, t0, dt, y0, &sample,
                )
            },
            | Method::SrkPlaten => {
                kernel::srk_noncommutative(drift, diffusion, noise, t0, dt, y0, &sample)
            },
            | _ => kernel::milstein_noncommutative(drift, diffusion, noise, t0, dt, y0, &sample),
        }
        .map_err(SdeError::from)?
    } else {
        let increments = wiener(options.seed, noise, steps, dt);
        match options.method {
            | Method::EulerMaruyama => {
                kernel::euler_maruyama(drift, diffusion, noise, t0, dt, y0, &increments)
            },
            | Method::Milstein => {
                kernel::milstein(drift, diffusion, noise, t0, dt, y0, &increments)
            },
            | Method::MilsteinStratonovich => {
                kernel::milstein_stratonovich(drift, diffusion, noise, t0, dt, y0, &increments)
            },
            | Method::HeunStratonovich => {
                kernel::heun_stratonovich(drift, diffusion, noise, t0, dt, y0, &increments)
            },
            | Method::SrkPlaten => {
                kernel::srk_platen(drift, diffusion, noise, t0, dt, y0, &increments)
            },
            | Method::DriftImplicitEuler
            | Method::ThetaMaruyama
            | Method::DriftImplicitMilstein => {
                let scheme = match options.method {
                    | Method::DriftImplicitEuler => DriftImplicitScheme::BackwardEuler,
                    | Method::ThetaMaruyama => DriftImplicitScheme::Trapezoidal,
                    | Method::DriftImplicitMilstein => DriftImplicitScheme::DriftImplicitMilstein,
                    | _ => unreachable!(),
                };
                Ok(drift_implicit_core(
                    &drift,
                    &diffusion,
                    noise,
                    t0,
                    dt,
                    y0,
                    &increments,
                    scheme,
                )?)
            },
            // The explicit Lévy-area methods are routed by `uses_levy_areas`.
            | Method::MilsteinNoncommutative | Method::MilsteinNoncommutativeStratonovich => {
                return Err(SdeError::invalid(
                    "internal: a Lévy-area method took the plain path",
                ));
            },
        }
        .map_err(SdeError::from)?
    };
    if let Some(failure) = adapters.failure() {
        return Err(failure);
    }
    let t: Vec<f64> = (0..states.len())
        .map(|step| t0 + step as f64 * dt)
        .collect();
    Ok(Path {
        t,
        steps: states.len().saturating_sub(1),
        y: states,
        rejected: 0,
        nfev: adapters.nfev(),
        seed: options.seed,
        noise,
        method: options.method,
    })
}

/// An adaptive-step integration.
#[allow(clippy::too_many_arguments)] // mirrors the kernel's adaptive entry points
fn adaptive<A, G>(
    adapters: &Adapters<A, G>,
    y0: &[f64],
    options: &Options,
    t1: f64,
    rtol: f64,
    atol: f64,
    max_steps: usize,
    cells: usize,
) -> Result<Path, SdeError>
where
    A: FnMut(f64, &[f64]) -> Result<Vec<f64>, SdeError>,
    G: FnMut(f64, &[f64]) -> Result<Vec<f64>, SdeError>,
{
    let noise = options.structure.noise(y0.len());
    let drift = |t: f64, y: &[f64], out: &mut [f64]| adapters.drift_at(t, y, out);
    let diffusion = |t: f64, y: &[f64], out: &mut [f64]| adapters.diffusion_at(t, y, out);
    let stratonovich = options.method.calculus() == Calculus::Stratonovich;
    let base = kernel::SdeAdaptiveOptions {
        rtol,
        atol,
        seed: options.seed,
        max_steps,
    };
    let path = if adaptive_noncommutative(options) {
        let cell_count = if cells == 0 {
            DEFAULT_CELLS
        } else {
            cells
        };
        let cell_width = (t1 - options.t0) / cell_count as f64;
        let terms = if options.terms == 0 {
            default_terms(cell_width)
        } else {
            options.terms
        };
        let nc = kernel::SdeNoncommutativeOptions {
            base,
            cells: cell_count,
            terms,
        };
        if stratonovich {
            kernel::sde_adaptive_noncommutative_stratonovich_with(
                drift, diffusion, noise, options.t0, t1, y0, &nc,
            )
        } else {
            kernel::sde_adaptive_noncommutative_with(
                drift, diffusion, noise, options.t0, t1, y0, &nc,
            )
        }
    } else if stratonovich {
        kernel::sde_adaptive_stratonovich_with(drift, diffusion, noise, options.t0, t1, y0, &base)
    } else {
        kernel::sde_adaptive_with(drift, diffusion, noise, options.t0, t1, y0, &base)
    }
    .map_err(SdeError::from)?;
    if let Some(failure) = adapters.failure() {
        return Err(failure);
    }
    Ok(Path {
        steps: path.y.len().saturating_sub(1),
        t: path.t,
        y: path.y,
        rejected: path.rejected,
        nfev: adapters.nfev(),
        seed: options.seed,
        noise,
        method: options.method,
    })
}

// ---------------------------------------------------------------------------
// Continuous output and events on paths (§8.9.1)
// ---------------------------------------------------------------------------

impl Path {
    /// The state at time `t`, by linear interpolation between the accepted
    /// nodes; `None` outside `[t_first, t_last]`.
    ///
    /// A numerical SDE path is only known at its nodes. Between them the
    /// interpolant is the chord, which is continuous and deterministic. It is
    /// as accurate as the scheme itself (strong order at most 1, and a Brownian
    /// path is only Hölder-1/2 between nodes anyway); it does not add the
    /// within-step Brownian-bridge variance, so statistics of `eval` between
    /// nodes understate the within-step variation by `O(dt)`. Use the node
    /// values (`y`) when exact distributional statistics at the nodes are
    /// wanted.
    #[must_use]
    pub fn eval(
        &self,
        t: f64,
    ) -> Option<Vec<f64>> {
        let (&first, &last) = (self.t.first()?, self.t.last()?);
        if !(first..=last).contains(&t) {
            return None;
        }
        let upper = self.t.partition_point(|&node| node < t);
        if upper == 0 {
            return self.y.first().cloned();
        }
        if crate::kernels::exact_eq(self.t[upper], t) {
            return Some(self.y[upper].clone());
        }
        let lower = upper - 1;
        let weight = (t - self.t[lower]) / (self.t[upper] - self.t[lower]);
        Some(
            self.y[lower]
                .iter()
                .zip(&self.y[upper])
                .map(|(a, b)| a + weight * (b - a))
                .collect(),
        )
    }

    /// The states at the requested times (`saveat`): [`Path::eval`] at each.
    ///
    /// # Errors
    ///
    /// `E0802` when a time lies outside the path's span.
    pub fn resample(
        &self,
        times: &[f64],
    ) -> Result<Vec<Vec<f64>>, SdeError> {
        times
            .iter()
            .map(|&t| {
                self.eval(t).ok_or_else(|| {
                    SdeError::invalid("a requested time lies outside the path's span")
                })
            })
            .collect()
    }
}

/// How one event function of a path is watched.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EventSpec {
    /// The path is cut at the first such event.
    pub terminal: bool,
    /// `1` rising only, `-1` falling only, `0` both.
    pub direction: i8,
}

/// An event on a path.
#[derive(Clone, Debug, PartialEq)]
pub struct Hit {
    /// Index of the event function.
    pub index: usize,
    /// The crossing time.
    pub t: f64,
    /// The interpolated state at the crossing.
    pub y: Vec<f64>,
}

/// Detects sign changes of `g_i(t, y)` along a path.
///
/// Between two consecutive nodes the path is the chord of [`Path::eval`]; a
/// crossing is a sign change of `g` at the two nodes in an allowed direction
/// (falling `g0 > 0 >= g1`, rising `g0 < 0 <= g1`, the kernel's convention
/// for ODE events), and its time is the root of `g` along the chord, found by
/// bisection. A node where `g` is exactly zero counts once (as the end of the
/// step that reaches it). Only the discretised path is watched: the
/// continuous process may cross and return between two nodes, which no
/// discrete scheme can see.
///
/// Hits come in time order. If an event is `terminal`, the returned path ends
/// at the first terminal hit (the hit is appended as its last node) and later
/// hits are dropped. The returned path is the input when nothing terminates.
pub fn events(
    path: &Path,
    specs: &[EventSpec],
    g: &mut dyn FnMut(usize, f64, &[f64]) -> f64,
) -> (Path, Vec<Hit>) {
    let mut hits: Vec<Hit> = Vec::new();
    let mut previous: Vec<f64> = match (path.t.first(), path.y.first()) {
        | (Some(&t), Some(y)) => (0..specs.len()).map(|i| g(i, t, y)).collect(),
        | _ => return (path.clone(), hits),
    };
    for k in 1..path.t.len() {
        let (t0, t1) = (path.t[k - 1], path.t[k]);
        let mut step_hits: Vec<Hit> = Vec::new();
        for (index, spec) in specs.iter().enumerate() {
            let g1 = g(index, t1, &path.y[k]);
            let g0 = previous[index];
            let crossing = (g0 < 0.0 && g1 >= 0.0 && spec.direction >= 0)
                || (g0 > 0.0 && g1 <= 0.0 && spec.direction <= 0);
            if crossing {
                let at = |s: f64| -> Vec<f64> {
                    path.y[k - 1]
                        .iter()
                        .zip(&path.y[k])
                        .map(|(a, b)| a + s * (b - a))
                        .collect()
                };
                let (mut lo, mut hi) = (0.0f64, 1.0f64);
                for _ in 0..64 {
                    let mid = f64::midpoint(lo, hi);
                    let value = g(index, t0 + mid * (t1 - t0), &at(mid));
                    if value == 0.0 {
                        hi = mid;
                        break;
                    }
                    if (value < 0.0) == (g0 < 0.0) {
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }
                let s = hi;
                step_hits.push(Hit {
                    index,
                    t: t0 + s * (t1 - t0),
                    y: at(s),
                });
            }
            previous[index] = g1;
        }
        step_hits.sort_by(|a, b| a.t.total_cmp(&b.t));
        if let Some(position) = step_hits.iter().position(|hit| specs[hit.index].terminal) {
            step_hits.truncate(position + 1);
            let last = step_hits[position].clone();
            hits.extend(step_hits);
            let mut cut = path.clone();
            cut.t.truncate(k);
            cut.y.truncate(k);
            cut.t.push(last.t);
            cut.y.push(last.y);
            cut.steps = cut.t.len() - 1;
            return (cut, hits);
        }
        hits.extend(step_hits);
    }
    (path.clone(), hits)
}

// ---------------------------------------------------------------------------
// Drift-Implicit Schemes for Stiff SDEs
// ---------------------------------------------------------------------------

/// Drift-implicit integration schemes for stiff Itô SDEs
/// `dX_t = a(t, X_t) dt + G(t, X_t) dW_t`.
pub use crate::kernels::sde_ext::DriftImplicitScheme;

fn map_sde_ext_err(err: crate::kernels::sde_ext::SdeExtError) -> SdeError {
    match err {
        | crate::kernels::sde_ext::SdeExtError::Invalid(s) => SdeError::invalid(s),
        | crate::kernels::sde_ext::SdeExtError::NonFinite(s) => SdeError::NonFinite(s),
        | crate::kernels::sde_ext::SdeExtError::StepLimit(s) => SdeError::StepLimit(s),
    }
}

#[allow(clippy::too_many_arguments)]
fn drift_implicit_core<D, F>(
    drift: &D,
    diffusion: &F,
    noise: usize,
    t0: f64,
    dt: f64,
    y0: &[f64],
    increments: &[Vec<f64>],
    scheme: DriftImplicitScheme,
) -> Result<Vec<Vec<f64>>, SdeError>
where
    D: Fn(f64, &[f64], &mut [f64]),
    F: Fn(f64, &[f64], &mut [f64]),
{
    crate::kernels::sde_ext::drift_implicit_core(
        drift, diffusion, noise, t0, dt, y0, increments, scheme,
    )
    .map_err(map_sde_ext_err)
}

/// Integrates a (potentially stiff) SDE `dX = a(t, X) dt + G(t, X) dW` using the
/// specified [`DriftImplicitScheme`].
///
/// Supports `BackwardEuler` (`theta = 1`), `Trapezoidal` (`theta = 0.5`), general
/// `Theta(theta)`, `DriftImplicitMilstein` (strong order 1.0), and
/// `SplitStepBackwardEuler` (for nonlinear stiff superlinear drifts).
///
/// # Errors
///
/// Returns [`SdeError`] on invalid options, mismatched callback dimensions,
/// non-finite values, or Newton non-convergence.
pub fn integrate_drift_implicit<A, G>(
    drift: A,
    diffusion: G,
    y0: &[f64],
    scheme: DriftImplicitScheme,
    options: &Options,
) -> Result<Path, SdeError>
where
    A: FnMut(f64, &[f64]) -> Result<Vec<f64>, SdeError>,
    G: FnMut(f64, &[f64]) -> Result<Vec<f64>, SdeError>,
{
    validate_options(options, y0)?;
    let noise = options.structure.noise(y0.len());
    let (dt, steps) = match &options.schedule {
        | Schedule::Fixed { dt, steps } => (*dt, *steps),
        | Schedule::Grid { times } => uniform_grid(times)?,
        | Schedule::Adaptive { .. } => {
            return Err(SdeError::invalid(
                "this solver requires a fixed or uniform-grid schedule",
            ));
        },
    };
    let adapters = Adapters::new(drift, diffusion, y0.len(), options.structure);
    let drift_fn = |t: f64, y: &[f64], out: &mut [f64]| adapters.drift_at(t, y, out);
    let diff_fn = |t: f64, y: &[f64], out: &mut [f64]| adapters.diffusion_at(t, y, out);
    let increments = wiener(options.seed, noise, steps, dt);
    let states = drift_implicit_core(
        &drift_fn,
        &diff_fn,
        noise,
        options.t0,
        dt,
        y0,
        &increments,
        scheme,
    )?;
    if let Some(failure) = adapters.failure() {
        return Err(failure);
    }
    let t: Vec<f64> = (0..states.len())
        .map(|k| options.t0 + k as f64 * dt)
        .collect();
    Ok(Path {
        steps: states.len().saturating_sub(1),
        t,
        y: states,
        rejected: 0,
        nfev: adapters.nfev(),
        seed: options.seed,
        noise,
        method: match scheme {
            | DriftImplicitScheme::Trapezoidal | DriftImplicitScheme::Theta(_) => {
                Method::ThetaMaruyama
            },
            | DriftImplicitScheme::DriftImplicitMilstein => Method::DriftImplicitMilstein,
            | DriftImplicitScheme::BackwardEuler | DriftImplicitScheme::SplitStepBackwardEuler => {
                Method::DriftImplicitEuler
            },
        },
    })
}

// ---------------------------------------------------------------------------
// Jump-Diffusion and Lévy-Driven SDEs
// ---------------------------------------------------------------------------

/// Distribution of jump marks `Z` in a compound Poisson jump-diffusion SDE
/// `dX_t = a(t, X_{t-}) dt + G(t, X_{t-}) dW_t + \sum_{q=1}^{dN_t} c(t, X_{t-}, Z_q)`.
#[derive(Clone, Debug, PartialEq)]
pub enum JumpDistribution {
    /// Deterministic jump mark vector `z`.
    Constant(Vec<f64>),
    /// Multivariate independent Normal jump marks `Z_j ~ N(mean[j], std[j]^2)`
    /// (Merton 1976 jump-diffusion model).
    Normal {
        /// Mean of each jump mark component.
        mean: Vec<f64>,
        /// Standard deviation of each jump mark component (`>= 0`).
        std: Vec<f64>,
    },
    /// Kou (2002) asymmetric double-exponential jump marks: with probability `p_up`
    /// draws `+Exp(eta_up)`, and with probability `1 - p_up` draws `-Exp(eta_down)`.
    DoubleExponential {
        /// Probability of an upward jump (`p_up in [0, 1]`).
        p_up: f64,
        /// Rate parameter of upward jumps (`eta_up > 0`, mean upward jump `1/eta_up`).
        eta_up: f64,
        /// Rate parameter of downward jumps (`eta_down > 0`, mean downward jump `1/eta_down`).
        eta_down: f64,
        /// Dimension of the jump mark vector.
        dim: usize,
    },
    /// Independent uniform jump marks `Z_j ~ Uniform(low[j], high[j])`.
    Uniform {
        /// Lower bounds.
        low: Vec<f64>,
        /// Upper bounds.
        high: Vec<f64>,
    },
    /// General `\alpha`-stable jump marks `Z_j ~ S_\alpha(scale, beta, location)`
    /// simulated via the Chambers–Mallows–Stuck algorithm.
    AlphaStable {
        /// Stability index `\alpha in (0, 2]`.
        alpha: f64,
        /// Skewness parameter `\beta in [-1, 1]`.
        beta: f64,
        /// Scale parameter `\gamma > 0`.
        scale: f64,
        /// Location parameter `\delta`.
        location: f64,
        /// Dimension of the jump mark vector.
        dim: usize,
    },
}

impl JumpDistribution {
    /// Samples one jump mark vector `Z` from the distribution using `stream`.
    ///
    /// # Errors
    /// Invalid shapes or parameters, non-finite inputs, or numerical solver failure.
    pub fn sample(
        &self,
        stream: &mut Stream,
    ) -> Result<Vec<f64>, SdeError> {
        match self {
            | Self::Constant(z) => {
                if z.is_empty() || !z.iter().all(|v| v.is_finite()) {
                    return Err(SdeError::invalid(
                        "constant jump mark must be non-empty and finite",
                    ));
                }
                Ok(z.clone())
            },
            | Self::Normal { mean, std } => {
                if mean.is_empty()
                    || mean.len() != std.len()
                    || !mean.iter().chain(std).all(|v| v.is_finite())
                    || std.iter().any(|&s| s < 0.0)
                {
                    return Err(SdeError::invalid("invalid Normal jump mark parameters"));
                }
                Ok(mean
                    .iter()
                    .zip(std)
                    .map(|(&m, &s)| m + s * stream.next_normal())
                    .collect())
            },
            | Self::DoubleExponential {
                p_up,
                eta_up,
                eta_down,
                dim,
            } => {
                if *dim == 0
                    || !(0.0..=1.0).contains(p_up)
                    || *eta_up <= 0.0
                    || *eta_down <= 0.0
                    || !p_up.is_finite()
                    || !eta_up.is_finite()
                    || !eta_down.is_finite()
                {
                    return Err(SdeError::invalid(
                        "invalid DoubleExponential (Kou) jump parameters",
                    ));
                }
                let mut out = Vec::with_capacity(*dim);
                for _ in 0..*dim {
                    let u_dir = stream.next_unit();
                    let u_exp = stream.next_unit().max(f64::MIN_POSITIVE);
                    if u_dir < *p_up {
                        out.push(-u_exp.ln() / *eta_up);
                    } else {
                        out.push(u_exp.ln() / *eta_down);
                    }
                }
                Ok(out)
            },
            | Self::Uniform { low, high } => {
                if low.is_empty()
                    || low.len() != high.len()
                    || !low.iter().chain(high).all(|v| v.is_finite())
                    || low.iter().zip(high).any(|(l, h)| l > h)
                {
                    return Err(SdeError::invalid("invalid Uniform jump parameters"));
                }
                Ok(low
                    .iter()
                    .zip(high)
                    .map(|(&l, &h)| l + (h - l) * stream.next_unit())
                    .collect())
            },
            | Self::AlphaStable {
                alpha,
                beta,
                scale,
                location,
                dim,
            } => {
                if *dim == 0 {
                    return Err(SdeError::invalid("AlphaStable jump dimension is zero"));
                }
                let mut out = Vec::with_capacity(*dim);
                for _ in 0..*dim {
                    out.push(sample_alpha_stable(
                        stream, *alpha, *beta, *scale, *location,
                    )?);
                }
                Ok(out)
            },
        }
    }

    /// Analytical mean mark `E[Z]` when finite and known (used for compensated Poisson jumps).
    #[must_use]
    pub fn mean_mark(&self) -> Option<Vec<f64>> {
        match self {
            | Self::Constant(z) => Some(z.clone()),
            | Self::Normal { mean, .. } => Some(mean.clone()),
            | Self::DoubleExponential {
                p_up,
                eta_up,
                eta_down,
                dim,
            } => {
                let m = *p_up / *eta_up - (1.0 - *p_up) / *eta_down;
                Some(vec![m; *dim])
            },
            | Self::Uniform { low, high } => {
                Some(
                    low.iter()
                        .zip(high)
                        .map(|(&l, &h)| f64::midpoint(l, h))
                        .collect(),
                )
            },
            | Self::AlphaStable {
                alpha, location, dim, ..
            } => {
                if *alpha > 1.0 {
                    Some(vec![*location; *dim])
                } else {
                    None
                }
            },
        }
    }
}

/// Specification of a compound Poisson jump process.
#[derive(Clone, Debug, PartialEq)]
pub struct JumpSpec {
    /// Poisson arrival intensity `\lambda >= 0` (expected number of jumps per unit time).
    pub intensity: f64,
    /// Distribution of the jump marks `Z`.
    pub distribution: JumpDistribution,
    /// Subtract the exact compensator for constant marks. Nonconstant marks
    /// are rejected because a general nonlinear coefficient needs E[c(t, X, Z)].
    pub compensated: bool,
}

/// Simulated trajectory of a jump-diffusion SDE.
#[derive(Clone, Debug, PartialEq)]
pub struct JumpPath {
    /// The underlying state trajectory (`t`, `y`, `steps`, `method`, etc.).
    pub path: Path,
    /// Total number of Poisson jumps that occurred over the trajectory.
    pub total_jumps: usize,
    /// Number of Poisson jumps that occurred in each time step `0..steps`.
    pub jump_counts: Vec<usize>,
}

/// Samples a `Poisson(mean)` random variable from `stream` using Knuth's exact
/// multiplication algorithm for `mean <= 30` and exact transformed rejection
/// for `mean > 30`.
#[must_use]
pub fn sample_poisson(
    stream: &mut Stream,
    mean: f64,
) -> usize {
    let mut next_u = || stream.next_unit();
    let mut unused_normal = || 0.0;
    crate::kernels::sde_ext::sample_poisson_with(&mut next_u, &mut unused_normal, mean)
}

/// Samples a general `\alpha`-stable random variable `S_\alpha(scale, beta, location)`
/// using the exact Chambers–Mallows–Stuck (1976) algorithm.
///
/// # Errors
/// Invalid shapes or parameters, non-finite inputs, or numerical solver failure.
pub fn sample_alpha_stable(
    stream: &mut Stream,
    alpha: f64,
    beta: f64,
    scale: f64,
    location: f64,
) -> Result<f64, SdeError> {
    if (alpha - 2.0).abs() < 1e-12 && alpha > 0.0 && (-1.0..=1.0).contains(&beta) && scale > 0.0 {
        let mut dummy_u = || 0.5;
        let mut next_n = || stream.next_normal();
        return crate::kernels::sde_ext::sample_alpha_stable_with(
            &mut dummy_u,
            &mut next_n,
            alpha,
            beta,
            scale,
            location,
        )
        .map_err(map_sde_ext_err);
    }
    let mut next_u = || stream.next_unit();
    let mut dummy_n = || 0.0;
    crate::kernels::sde_ext::sample_alpha_stable_with(
        &mut next_u,
        &mut dummy_n,
        alpha,
        beta,
        scale,
        location,
    )
    .map_err(map_sde_ext_err)
}

/// Generates `steps` vectors of `dim`-dimensional `\alpha`-stable Lévy increments
/// `\Delta L_k ~ S_\alpha(dt^{1/\alpha}, \beta, 0)` from the deterministic language stream.
///
/// # Errors
///
/// Returns [`SdeError::Invalid`] if `dim == 0`, `dt <= 0`, or `\alpha, \beta` are out of range.
pub fn levy_increments(
    seed: u64,
    dim: usize,
    steps: usize,
    dt: f64,
    alpha: f64,
    beta: f64,
) -> Result<Vec<Vec<f64>>, SdeError> {
    if dim == 0 || !dt.is_finite() || dt <= 0.0 {
        return Err(SdeError::invalid(
            "levy_increments requires dim > 0 and finite dt > 0",
        ));
    }
    let mut stream = Stream::new(seed);
    let scale = dt.powf(1.0 / alpha);
    let mut out = Vec::with_capacity(steps);
    for _ in 0..steps {
        let mut step_vec = Vec::with_capacity(dim);
        for _ in 0..dim {
            step_vec.push(sample_alpha_stable(&mut stream, alpha, beta, scale, 0.0)?);
        }
        out.push(step_vec);
    }
    Ok(out)
}

/// Integrates a jump-diffusion SDE
/// `dX_t = a(t, X_{t-}) dt + G(t, X_{t-}) dW_t + \sum_{q=1}^{dN_t} c(t, X_{t-}, Z_q)`
/// over a uniform step grid specified in `options`.
///
/// * `drift(t, y)` returns the drift vector `a(t, y) in R^n`.
/// * `diffusion(t, y)` returns the diffusion coefficients (`n` for `Diagonal`, `n * m` for `General`).
/// * `jump_coeff(t, y, z)` returns the state jump increment `c(t, y, z) in R^n` produced by a single
///   jump mark `z` (for example, in Merton's geometric jump-diffusion `dS = \mu S dt + \sigma S dW + S (e^Z - 1) dN`,
///   `c(t, S, z) = [S_0 * (z_0.exp() - 1.0)]`).
///
/// Supports `EulerMaruyama`, `Milstein`, `HeunStratonovich`, `DriftImplicitEuler`,
/// `ThetaMaruyama`, and `DriftImplicitMilstein` for the continuous diffusion part.
///
/// # Errors
///
/// Returns [`SdeError`] if parameters or shapes are invalid or the state becomes non-finite.
pub fn integrate_jump_diffusion<A, G, J>(
    mut drift: A,
    mut diffusion: G,
    mut jump_coeff: J,
    y0: &[f64],
    jump_spec: &JumpSpec,
    options: &Options,
) -> Result<JumpPath, SdeError>
where
    A: FnMut(f64, &[f64]) -> Result<Vec<f64>, SdeError>,
    G: FnMut(f64, &[f64]) -> Result<Vec<f64>, SdeError>,
    J: FnMut(f64, &[f64], &[f64]) -> Result<Vec<f64>, SdeError>,
{
    validate_options(options, y0)?;
    if !jump_spec.intensity.is_finite() || jump_spec.intensity < 0.0 {
        return Err(SdeError::invalid(
            "jump intensity lambda must be finite and non-negative",
        ));
    }
    let n = y0.len();
    let noise = options.structure.noise(n);

    let (dt, steps) = match &options.schedule {
        | Schedule::Fixed { dt, steps } => (*dt, *steps),
        | Schedule::Grid { times } => uniform_grid(times)?,
        | Schedule::Adaptive { .. } => {
            return Err(SdeError::invalid(
                "this solver requires a fixed or uniform-grid schedule",
            ));
        },
    };

    if !options.structure.is_commutative()
        || matches!(
            options.method,
            Method::MilsteinNoncommutative | Method::MilsteinNoncommutativeStratonovich
        )
    {
        return Err(SdeError::invalid(
            "jump diffusion currently requires commutative noise",
        ));
    }
    // c(E[Z]) is not E[c(Z)] for a nonlinear jump coefficient. Only constant
    // marks have an exact compensator without an additional expectation callback.
    let mean_mark = if jump_spec.compensated {
        match &jump_spec.distribution {
            | JumpDistribution::Constant(mark) => Some(mark.clone()),
            | _ => {
                return Err(SdeError::invalid(
                    "compensation with nonconstant marks requires an explicit expected jump coefficient",
                ));
            },
        }
    } else {
        None
    };

    let increments = wiener(options.seed, noise, steps, dt);
    // Use a distinct deterministic seed stream for the Poisson jump counts and marks
    let mut jump_stream = Stream::new(ensemble_seed(options.seed, 0x4A55_4D50));

    let mut ts = Vec::with_capacity(steps + 1);
    let mut ys = Vec::with_capacity(steps + 1);
    let mut jump_counts = Vec::with_capacity(steps);
    let mut total_jumps = 0usize;
    let mut nfev = 0usize;

    let mut y = y0.to_vec();
    ts.push(options.t0);
    ys.push(y.clone());

    for (step_idx, dw) in increments.iter().enumerate() {
        let t_curr = options.t0 + step_idx as f64 * dt;
        let t_next = options.t0 + (step_idx + 1) as f64 * dt;

        let mut y_cont = {
            let adjusted_drift = |t: f64, state: &[f64]| {
                let mut value = drift(t, state)?;
                if value.len() != n {
                    return Err(SdeError::shape("invalid jump-diffusion drift shape"));
                }
                if let Some(mark) = &mean_mark {
                    let comp = jump_coeff(t, state, mark)?;
                    if comp.len() != n {
                        return Err(SdeError::shape("invalid compensator shape"));
                    }
                    for i in 0..n {
                        value[i] -= jump_spec.intensity * comp[i];
                    }
                }
                Ok(value)
            };
            let adapters = Adapters::new(adjusted_drift, &mut diffusion, n, options.structure);
            let drift_fn =
                |t: f64, state: &[f64], out: &mut [f64]| adapters.drift_at(t, state, out);
            let diff_fn =
                |t: f64, state: &[f64], out: &mut [f64]| adapters.diffusion_at(t, state, out);
            let next = simulate_terminal_with_increments(
                &drift_fn,
                &diff_fn,
                noise,
                t_curr,
                dt,
                &y,
                std::slice::from_ref(dw),
                options.method,
            );
            if let Some(error) = adapters.failure() {
                return Err(error);
            }
            nfev += adapters.nfev();
            next?
        };

        // Compound Poisson jumps on [t_curr, t_next]
        let n_jumps = sample_poisson(&mut jump_stream, jump_spec.intensity * dt);
        jump_counts.push(n_jumps);
        total_jumps += n_jumps;

        for _ in 0..n_jumps {
            let z = jump_spec.distribution.sample(&mut jump_stream)?;
            let delta_jump = jump_coeff(t_curr, &y_cont, &z)?;
            if delta_jump.len() != n {
                return Err(SdeError::Shape(format!(
                    "jump_coeff returned {} values for state of {n}",
                    delta_jump.len()
                )));
            }
            if !delta_jump.iter().all(|v| v.is_finite()) {
                return Err(SdeError::NonFinite(
                    "jump_coeff returned non-finite value".into(),
                ));
            }
            for i in 0..n {
                y_cont[i] += delta_jump[i];
            }
        }

        if !y_cont.iter().all(|v| v.is_finite()) {
            return Err(SdeError::NonFinite(
                "jump-diffusion state became non-finite".into(),
            ));
        }

        y = y_cont;
        ts.push(t_next);
        ys.push(y.clone());
    }

    Ok(JumpPath {
        path: Path {
            t: ts,
            y: ys,
            steps,
            rejected: 0,
            nfev,
            seed: options.seed,
            noise,
            method: options.method,
        },
        total_jumps,
        jump_counts,
    })
}

/// Integrate a jump diffusion with an explicit compensator.
///
/// Integrate a compensated jump diffusion with an explicit expected jump coefficient.
/// `expected_jump(t, y)` must return E[c(t, y, Z)], not c(t, y, E[Z]).
/// This supports nonlinear coefficients and any mark distribution with a finite expectation.
///
/// # Errors
/// Invalid integration inputs, callback shapes, or non-finite expected coefficients.
pub fn integrate_jump_diffusion_compensated<A, G, J, C>(
    mut drift: A,
    diffusion: G,
    jump_coeff: J,
    mut expected_jump: C,
    y0: &[f64],
    jump_spec: &JumpSpec,
    options: &Options,
) -> Result<JumpPath, SdeError>
where
    A: FnMut(f64, &[f64]) -> Result<Vec<f64>, SdeError>,
    G: FnMut(f64, &[f64]) -> Result<Vec<f64>, SdeError>,
    J: FnMut(f64, &[f64], &[f64]) -> Result<Vec<f64>, SdeError>,
    C: FnMut(f64, &[f64]) -> Result<Vec<f64>, SdeError>,
{
    let adjusted = |t: f64, y: &[f64]| {
        let mut value = drift(t, y)?;
        let expected = expected_jump(t, y)?;
        if value.len() != y.len() || expected.len() != y.len() {
            return Err(SdeError::shape("invalid expected jump coefficient shape"));
        }
        if expected.iter().any(|v| !v.is_finite()) {
            return Err(SdeError::NonFinite(
                "non-finite expected jump coefficient".into(),
            ));
        }
        for (v, e) in value.iter_mut().zip(expected) {
            *v -= jump_spec.intensity * e;
        }
        Ok(value)
    };
    let mut uncompensated = jump_spec.clone();
    uncompensated.compensated = false;
    integrate_jump_diffusion(adjusted, diffusion, jump_coeff, y0, &uncompensated, options)
}

/// Integrates an `\alpha`-stable Lévy-driven SDE
/// `dX_t = a(t, X_{t-}) dt + G(t, X_{t-}) dL_t^{(\alpha, \beta)}`
/// where `dL_t` has independent `S_\alpha(dt^{1/\alpha}, \beta, 0)` increments.
///
/// # Errors
///
/// Returns [`SdeError`] for invalid `\alpha \notin (0, 2]`, `\beta \notin [-1, 1]`,
/// mismatched callback shapes, non-finite states, or a method other than `EulerMaruyama`.
pub fn integrate_levy<A, G>(
    mut drift: A,
    mut levy_diffusion: G,
    y0: &[f64],
    alpha: f64,
    beta: f64,
    options: &Options,
) -> Result<Path, SdeError>
where
    A: FnMut(f64, &[f64]) -> Result<Vec<f64>, SdeError>,
    G: FnMut(f64, &[f64]) -> Result<Vec<f64>, SdeError>,
{
    validate_options(options, y0)?;
    if options.method != Method::EulerMaruyama {
        return Err(SdeError::invalid(
            "alpha-stable Levy integration currently requires EulerMaruyama",
        ));
    }
    let n = y0.len();
    let noise = options.structure.noise(n);
    let expected_diff_len = options.structure.coefficients(n);
    let (dt, steps) = match &options.schedule {
        | Schedule::Fixed { dt, steps } => (*dt, *steps),
        | Schedule::Grid { times } => uniform_grid(times)?,
        | Schedule::Adaptive { .. } => {
            return Err(SdeError::invalid(
                "this solver requires a fixed or uniform-grid schedule",
            ));
        },
    };
    let dl_increments = levy_increments(options.seed, noise, steps, dt, alpha, beta)?;
    let mut ts = Vec::with_capacity(steps + 1);
    let mut ys = Vec::with_capacity(steps + 1);
    let mut y = y0.to_vec();
    let mut nfev = 0usize;

    ts.push(options.t0);
    ys.push(y.clone());

    for (step_idx, dl) in dl_increments.iter().enumerate() {
        let t_curr = options.t0 + step_idx as f64 * dt;
        let t_next = options.t0 + (step_idx + 1) as f64 * dt;

        let a_val = drift(t_curr, &y)?;
        let g_val = levy_diffusion(t_curr, &y)?;
        nfev += 1;
        if a_val.len() != n || g_val.len() != expected_diff_len {
            return Err(SdeError::Shape(
                "mismatched drift or diffusion shape in integrate_levy".into(),
            ));
        }
        if !a_val.iter().chain(&g_val).all(|v| v.is_finite()) {
            return Err(SdeError::NonFinite(
                "non-finite rate in integrate_levy".into(),
            ));
        }
        for i in 0..n {
            let diff_i = match options.structure {
                | Structure::Diagonal => g_val[i] * dl[i],
                | Structure::General { .. } => {
                    (0..noise).map(|j| g_val[i * noise + j] * dl[j]).sum()
                },
            };
            y[i] += a_val[i] * dt + diff_i;
        }
        if !y.iter().all(|v| v.is_finite()) {
            return Err(SdeError::NonFinite(
                "state left the reals in integrate_levy".into(),
            ));
        }
        ts.push(t_next);
        ys.push(y.clone());
    }

    Ok(Path {
        t: ts,
        y: ys,
        steps,
        rejected: 0,
        nfev,
        seed: options.seed,
        noise,
        method: options.method,
    })
}

// ---------------------------------------------------------------------------
// Multilevel Monte Carlo (MLMC)
// ---------------------------------------------------------------------------

/// Configuration for Multilevel Monte Carlo (MLMC) estimation (Giles 2008).
#[derive(Clone, Copy, Debug)]
pub struct MlmcOptions {
    /// Time-stepping scheme used on each level (`EulerMaruyama`, `Milstein`,
    /// `DriftImplicitEuler`, `ThetaMaruyama`, `DriftImplicitMilstein`, `HeunStratonovich`,
    /// or `SrkPlaten`).
    pub method: Method,
    /// Noise structure (`Diagonal` or `General`).
    pub structure: Structure,
    /// Base RNG seed.
    pub seed: u64,
    /// Initial time `t0`.
    pub t0: f64,
    /// Terminal time `t1 > t0`.
    pub t1: f64,
    /// Number of steps on the coarsest level `l = 0` (`steps_l = base_steps * 2^l`).
    pub base_steps: usize,
    /// Minimum hierarchy depth `L >= min_level` before testing weak bias convergence.
    pub min_level: usize,
    /// Maximum hierarchy depth `L <= max_level`.
    pub max_level: usize,
    /// Initial warmup samples per level used to estimate level variances `V_l`.
    pub warmup_samples: usize,
    /// Maximum samples per level cap (guards against runaway allocations in tests).
    pub max_samples_per_level: usize,
    /// Target root-mean-square error `\varepsilon > 0`.
    pub target_rmse: f64,
}

impl Default for MlmcOptions {
    fn default() -> Self {
        Self {
            method: Method::Milstein,
            structure: Structure::Diagonal,
            seed: 42,
            t0: 0.0,
            t1: 1.0,
            base_steps: 4,
            min_level: 2,
            max_level: 6,
            warmup_samples: 40,
            max_samples_per_level: 20_000,
            target_rmse: 1e-2,
        }
    }
}

/// Per-level statistics of a Multilevel Monte Carlo estimation.
#[derive(Clone, Debug, PartialEq)]
pub struct MlmcLevelStats {
    /// Level index `l` (`0..=L`).
    pub level: usize,
    /// Fine time steps `base_steps * 2^l` on this level.
    pub steps: usize,
    /// Number of coupled paths `N_l` simulated at this level.
    pub samples: usize,
    /// Sample mean of `P_f - P_c` (`P_0` on level `0`).
    pub mean_diff: f64,
    /// Sample variance `V_l = Var(P_f - P_c)` (`Var(P_0)` on level `0`).
    pub variance: f64,
}

/// Result of a Multilevel Monte Carlo (MLMC) estimation of `E[P(X_{t1})]`.
#[derive(Clone, Debug, PartialEq)]
pub struct MlmcResult {
    /// Multilevel telescoping estimate `\hat{Y} = \sum_{l=0}^L \hat{Y}_l`.
    pub estimate: f64,
    /// Total estimator variance `\sum_{l=0}^L V_l / N_l`.
    pub variance: f64,
    /// Standard error `\sqrt{\text{variance}}`.
    pub std_error: f64,
    /// Statistics for each level `l = 0..=L`.
    pub levels: Vec<MlmcLevelStats>,
    /// Total number of SDE time steps evaluated across all levels and samples.
    pub total_cost: usize,
}

/// Evaluates a terminal path state using `method` on pre-sampled Brownian increments.
// Preserve the established solver API with explicit grid, rates, and callbacks.
#[allow(clippy::too_many_arguments)]
fn simulate_terminal_with_increments<D, F>(
    drift: &D,
    diffusion: &F,
    noise: usize,
    t0: f64,
    dt: f64,
    y0: &[f64],
    increments: &[Vec<f64>],
    method: Method,
) -> Result<Vec<f64>, SdeError>
where
    D: Fn(f64, &[f64], &mut [f64]),
    F: Fn(f64, &[f64], &mut [f64]),
{
    let states = match method {
        | Method::EulerMaruyama => {
            kernel::euler_maruyama(drift, diffusion, noise, t0, dt, y0, increments)
                .map_err(SdeError::from)?
        },
        | Method::Milstein | Method::MilsteinNoncommutative => {
            kernel::milstein(drift, diffusion, noise, t0, dt, y0, increments)
                .map_err(SdeError::from)?
        },
        | Method::MilsteinStratonovich | Method::MilsteinNoncommutativeStratonovich => {
            kernel::milstein_stratonovich(drift, diffusion, noise, t0, dt, y0, increments)
                .map_err(SdeError::from)?
        },
        | Method::HeunStratonovich => {
            kernel::heun_stratonovich(drift, diffusion, noise, t0, dt, y0, increments)
                .map_err(SdeError::from)?
        },
        | Method::SrkPlaten => {
            kernel::srk_platen(drift, diffusion, noise, t0, dt, y0, increments)
                .map_err(SdeError::from)?
        },
        | Method::DriftImplicitEuler => {
            drift_implicit_core(
                drift,
                diffusion,
                noise,
                t0,
                dt,
                y0,
                increments,
                DriftImplicitScheme::BackwardEuler,
            )?
        },
        | Method::ThetaMaruyama => {
            drift_implicit_core(
                drift,
                diffusion,
                noise,
                t0,
                dt,
                y0,
                increments,
                DriftImplicitScheme::Trapezoidal,
            )?
        },
        | Method::DriftImplicitMilstein => {
            drift_implicit_core(
                drift,
                diffusion,
                noise,
                t0,
                dt,
                y0,
                increments,
                DriftImplicitScheme::DriftImplicitMilstein,
            )?
        },
    };
    states
        .last()
        .cloned()
        .ok_or_else(|| SdeError::invalid("empty trajectory in MLMC"))
}

/// Estimate a terminal SDE expectation with multilevel Monte Carlo.
///
/// Estimates the expected value `E[payoff(X_{t1})]` of a terminal functional of an
/// SDE `dX = a(t, X) dt + G(t, X) dW` using Giles' Multilevel Monte Carlo (MLMC)
/// method with coupled fine/coarse Brownian paths (`\Delta W^{(c)}_k = \Delta W^{(f)}_{2k} + \Delta W^{(f)}_{2k+1}`)
/// and optimal sample allocation `N_l \propto \sqrt{V_l / C_l}`.
///
/// # Errors
///
/// Returns [`SdeError`] if options or callback dimensions are invalid or if any
/// simulated trajectory becomes non-finite.
pub fn mlmc<A, G, P>(
    drift: A,
    diffusion: G,
    mut payoff: P,
    y0: &[f64],
    options: &MlmcOptions,
) -> Result<MlmcResult, SdeError>
where
    A: FnMut(f64, &[f64]) -> Result<Vec<f64>, SdeError>,
    G: FnMut(f64, &[f64]) -> Result<Vec<f64>, SdeError>,
    P: FnMut(&[f64]) -> Result<f64, SdeError>,
{
    if y0.is_empty() || !y0.iter().all(|v| v.is_finite()) {
        return Err(SdeError::invalid(
            "MLMC requires a non-empty finite initial state",
        ));
    }
    if !options.t0.is_finite()
        || !options.t1.is_finite()
        || options.t1 <= options.t0
        || options.base_steps == 0
        || options.target_rmse <= 0.0
        || !options.target_rmse.is_finite()
        || options.max_level < options.min_level
    {
        return Err(SdeError::invalid("invalid MlmcOptions parameters"));
    }
    let n = y0.len();
    let noise = options.structure.noise(n);
    if noise == 0 {
        return Err(SdeError::invalid("MLMC noise dimension is zero"));
    }

    if !options.structure.is_commutative()
        || matches!(
            options.method,
            Method::MilsteinNoncommutative | Method::MilsteinNoncommutativeStratonovich
        )
    {
        return Err(SdeError::invalid(
            "MLMC coupling currently requires commutative noise",
        ));
    }
    if options.max_level >= usize::BITS as usize
        || options
            .base_steps
            .checked_shl(options.max_level as u32)
            .is_none()
        || options.base_steps > (usize::MAX >> options.max_level)
    {
        return Err(SdeError::invalid("MLMC level step count overflows"));
    }
    let adapters = Adapters::new(drift, diffusion, n, options.structure);
    let drift_fn = |t: f64, y: &[f64], out: &mut [f64]| adapters.drift_at(t, y, out);
    let diff_fn = |t: f64, y: &[f64], out: &mut [f64]| adapters.diffusion_at(t, y, out);

    let span = options.t1 - options.t0;
    let max_lv = options.max_level;
    let mut n_samples = vec![0usize; max_lv + 1];
    let mut sum_y = vec![0.0f64; max_lv + 1];
    let mut sum_y2 = vec![0.0f64; max_lv + 1];
    let mut total_cost = 0usize;

    let sample_level_batch = |level: usize,
                              count: usize,
                              start_idx: usize,
                              payoff: &mut P|
     -> Result<(f64, f64, usize), SdeError> {
        let fine_steps = options.base_steps * (1usize << level);
        let dt_fine = span / fine_steps as f64;
        let mut s1 = 0.0f64;
        let mut s2 = 0.0f64;
        let mut cost = 0usize;

        for idx in start_idx..start_idx + count {
            let path_seed = ensemble_seed(
                options
                    .seed
                    .wrapping_add((level as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15)),
                idx,
            );
            let fine_inc = wiener(path_seed, noise, fine_steps, dt_fine);
            let y_fine = simulate_terminal_with_increments(
                &drift_fn,
                &diff_fn,
                noise,
                options.t0,
                dt_fine,
                y0,
                &fine_inc,
                options.method,
            )?;
            if let Some(err) = adapters.failure() {
                return Err(err);
            }
            let p_fine = payoff(&y_fine)?;
            if !p_fine.is_finite() {
                return Err(SdeError::NonFinite(
                    "MLMC payoff returned non-finite".into(),
                ));
            }

            let diff_val = if level == 0 {
                cost += fine_steps;
                p_fine
            } else {
                let coarse_steps = fine_steps / 2;
                let dt_coarse = 2.0 * dt_fine;
                // Exact Brownian coupling: dW_coarse[k] = dW_fine[2k] + dW_fine[2k + 1]
                let mut coarse_inc = Vec::with_capacity(coarse_steps);
                for k in 0..coarse_steps {
                    let mut cw = Vec::with_capacity(noise);
                    for j in 0..noise {
                        cw.push(fine_inc[2 * k][j] + fine_inc[2 * k + 1][j]);
                    }
                    coarse_inc.push(cw);
                }
                let y_coarse = simulate_terminal_with_increments(
                    &drift_fn,
                    &diff_fn,
                    noise,
                    options.t0,
                    dt_coarse,
                    y0,
                    &coarse_inc,
                    options.method,
                )?;
                if let Some(err) = adapters.failure() {
                    return Err(err);
                }
                let p_coarse = payoff(&y_coarse)?;
                if !p_coarse.is_finite() {
                    return Err(SdeError::NonFinite(
                        "MLMC payoff returned non-finite".into(),
                    ));
                }
                cost += fine_steps + coarse_steps;
                p_fine - p_coarse
            };
            s1 += diff_val;
            s2 += diff_val * diff_val;
        }
        Ok((s1, s2, cost))
    };

    let warmup = options.warmup_samples.max(10);
    let mut active_max_l = options.min_level.min(max_lv);

    for l in 0..=active_max_l {
        let (s1, s2, c) = sample_level_batch(l, warmup, 0, &mut payoff)?;
        n_samples[l] += warmup;
        sum_y[l] += s1;
        sum_y2[l] += s2;
        total_cost += c;
    }

    loop {
        // Compute variances V_l and costs C_l for levels 0..=active_max_l
        let mut variances = vec![0.0f64; active_max_l + 1];
        let mut costs = vec![0.0f64; active_max_l + 1];
        let mut sum_sqrt_vc = 0.0f64;

        for l in 0..=active_max_l {
            let nl = n_samples[l] as f64;
            let mean = sum_y[l] / nl;
            let var = ((sum_y2[l] / nl) - mean * mean).max(0.0);
            variances[l] = var;
            let fine_steps = (options.base_steps * (1usize << l)) as f64;
            costs[l] = if l == 0 {
                fine_steps
            } else {
                1.5 * fine_steps
            };
            sum_sqrt_vc += (var * costs[l]).sqrt();
        }

        // Giles optimal sample counts: N_l = ceil((2 / eps^2) * sqrt(V_l / C_l) * sum_k sqrt(V_k C_k))
        let eps2 = options.target_rmse * options.target_rmse;
        for l in 0..=active_max_l {
            let target_nl =
                ((2.0 / eps2) * (variances[l] / costs[l]).sqrt() * sum_sqrt_vc).ceil() as usize;
            let capped_nl = target_nl.clamp(warmup, options.max_samples_per_level.max(warmup));
            if capped_nl > n_samples[l] {
                let extra = capped_nl - n_samples[l];
                let (s1, s2, c) = sample_level_batch(l, extra, n_samples[l], &mut payoff)?;
                n_samples[l] += extra;
                sum_y[l] += s1;
                sum_y2[l] += s2;
                total_cost += c;
            }
        }

        // Test weak convergence bias at active_max_l
        let mean_l = (sum_y[active_max_l] / n_samples[active_max_l] as f64).abs();
        let mean_prev = if active_max_l >= 1 {
            (sum_y[active_max_l - 1] / n_samples[active_max_l - 1] as f64).abs()
        } else {
            0.0
        };
        let rem_bias = mean_l.max(0.5 * mean_prev);
        if rem_bias <= options.target_rmse / 2.0f64.sqrt() || active_max_l >= max_lv {
            break;
        }
        active_max_l += 1;
        let (s1, s2, c) = sample_level_batch(active_max_l, warmup, 0, &mut payoff)?;
        n_samples[active_max_l] += warmup;
        sum_y[active_max_l] += s1;
        sum_y2[active_max_l] += s2;
        total_cost += c;
    }

    let mut estimate = 0.0f64;
    let mut total_var = 0.0f64;
    let mut levels = Vec::with_capacity(active_max_l + 1);

    for l in 0..=active_max_l {
        let nl = n_samples[l];
        let nl_f = nl as f64;
        let mean_diff = sum_y[l] / nl_f;
        let variance = if nl > 1 {
            ((sum_y2[l] - nl_f * mean_diff * mean_diff) / (nl_f - 1.0)).max(0.0)
        } else {
            0.0
        };
        estimate += mean_diff;
        total_var += variance / nl_f;
        levels.push(MlmcLevelStats {
            level: l,
            steps: options.base_steps * (1usize << l),
            samples: nl,
            mean_diff,
            variance,
        });
    }

    Ok(MlmcResult {
        estimate,
        variance: total_var,
        std_error: total_var.sqrt(),
        levels,
        total_cost,
    })
}
