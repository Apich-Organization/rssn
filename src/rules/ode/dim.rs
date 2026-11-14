//! Semantic dimension rules for the differential-equation operators.
//!
//! The numeric spectral operators (`bessel_zero`, `sl_eigen`, `hill_mode`, the
//! Coulomb and spheroidal functions, ...) are pure-number maps: their kernels
//! receive and return `f64`, so they are annotated with
//! `DimRule::Dimensionless` and the physical dimension is carried by the
//! symbolic solver around them. The operators whose value is an equation, a
//! list of symmetries or a conservation-law report are `NonMetric`.
//!
//! The remaining operators share one dispatcher, [`rule`], because their
//! result dimension is a formula in their arguments: a solution carries the
//! dimension of the dependent variable, a soliton amplitude scales with the
//! spectral parameter, and a Green's function scales with the domain radius.

// The dispatcher must be `pub(crate)` for the `pde` rule set to install it,
// while `rules::ode::dim` is only reachable inside the crate; `pub(crate)` is
// already the widest visibility, so the lint's `pub(super)` suggestion would
// hide the function from `pde`.
#![allow(clippy::redundant_pub_crate)]

use crate::graph::dim::DimCtx;
use crate::graph::dim::DimError;
use crate::graph::dim::DimTerm;

/// The dimension rule of an `ode`/`pde` operator whose result is an expression
/// rather than an equation or a label.
///
/// `args` are the dimensions of the operator's arguments, so the rule only has
/// to describe how the value is built from them.
// The `Result` is part of the `DimRule::Custom` signature; no arm of the
// dispatch can fail, so every rule returns `Ok`.
#[allow(clippy::unnecessary_wraps)]
pub(crate) fn rule(
    ctx: &DimCtx<'_>,
    args: &[DimTerm],
) -> Result<DimTerm, DimError> {
    let op = ctx.graph.op(ctx.node);
    let name = &*ctx.graph.ops().get(op).name;
    let dimension = match name {
        // A solution, separated product or Green's function is an expression
        // in the dependent variable `u(...)`/`y(...)` of argument 1.
        | "ode_asymptotic"
        | "ode_borel"
        | "pde_separate"
        | "pde_greens_function"
        | "conservation_laws" => dependent(args, 1),
        // The Noether current density carries the dimension of the Lagrangian
        // density; the first hodograph solution component has that of `F(r)`.
        | "noether_current" | "pde_hodograph_solution" => dependent(args, 0),
        // `u = -2 ∂² ln det(I + A)` for a soliton: the determinant is
        // dimensionless and every derivative of `ln` scales as the spectral
        // parameter `κ`, so `u` scales as `κ²`.
        | "pde_kdv_soliton" => square(&dependent(args, 0)),
        | "pde_kdv_reflectionless" => square(&dependent(args, 1)),
        // Unit-normalised or phase-valued results of purely symbolic origin:
        // the PDE order (a count), a point value of the numeric Gelfand-Levitan
        // solver, and the sine-Gordon kink (a phase).
        | "pde_order" | "pde_kdv_scattering" | "pde_sine_gordon_soliton" => {
            DimTerm::dimensionless()
        },
        // Dirichlet/Neumann Green's functions: `1/(4π|x - s|)` in 3D (so
        // `1/R`) and a logarithm in 2D.
        | "pde_greens_ball" | "pde_greens_ball_neumann" => {
            if child_index(ctx, 1) == Some(3) {
                dependent(args, 0).inverse()
            } else {
                DimTerm::dimensionless()
            }
        },
        // Half-space Robin: `G = G₀ + G₀* - 2h ∫ e^{-hw} G₀ dw`. With
        // `G₀ ~ 1/|x|` in 3D the integral contributes `1/h`, so `G` has the
        // dimension of `h`; in 2D `G₀` is a logarithm and `G` is
        // dimensionless.
        | "pde_greens_robin" => {
            if child_index(ctx, 0) == Some(3) {
                dependent(args, 1)
            } else {
                DimTerm::dimensionless()
            }
        },
        // Interior ball Robin: `1/(4π|x - s|) + Σ A_l (r s)^l P_l(cos γ)`,
        // whose every term scales as `1/R`.
        | "pde_greens_ball_robin" => dependent(args, 0).inverse(),
        // Exterior Robin: a logarithm in 2D and `1/R` in 3D.
        | "pde_greens_exterior_robin" => {
            if child_index(ctx, 2) == Some(3) {
                dependent(args, 0).inverse()
            } else {
                DimTerm::dimensionless()
            }
        },
        // Only operators registered with this rule reach the dispatcher.
        | _ => DimTerm::dimensionless(),
    };
    Ok(dimension)
}

/// The dimension of argument `index`, or `1` when that argument is not a
/// scalar (an equation, a condition list, a structural container).
fn dependent(
    args: &[DimTerm],
    index: usize,
) -> DimTerm {
    args.get(index)
        .filter(|term| !term.is_non_metric())
        .cloned()
        .unwrap_or_else(DimTerm::dimensionless)
}

/// The square of a dimension.
fn square(dimension: &DimTerm) -> DimTerm {
    dimension.mul(dimension)
}

/// The literal integer in the `index`-th child, when it is one (the spatial
/// dimension argument of the Green's-function operators).
fn child_index(
    ctx: &DimCtx<'_>,
    index: usize,
) -> Option<i64> {
    ctx.graph
        .children(ctx.node)
        .get(index)
        .and_then(|&node| ctx.graph.number_of(node))
        .and_then(crate::graph::Number::to_i64)
}
