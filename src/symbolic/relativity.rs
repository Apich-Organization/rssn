//! # Relativity Module
//!
//! This module provides symbolic tools for Special and General Relativity.
//!
//! ## Key Concepts:
//! - **Special Relativity**: Lorentz transformations, time dilation, length contraction, and relativistic dynamics.
//! - **General Relativity**: Einstein field equations, Schwarzschild metric, and geodesic equations.
//! - **Four-Vectors**: Unified representation of space and time components.

use std::sync::Arc;

use crate::symbolic::core::Expr;

/// Calculates the Lorentz factor: $\gamma = \frac{1}{\sqrt{1 - v^2/c^2}}$.
///
/// This factor appears in all relativistic calculations for time dilation,
/// length contraction, and relativistic mass.
#[must_use]

pub fn lorentz_factor(velocity : &Expr) -> Expr {

    let c = Expr::new_variable("c");

    let v_sq = Expr::new_pow(
        velocity.clone(),
        Expr::Constant(2.0),
    );

    let c_sq = Expr::new_pow(
        c,
        Expr::Constant(2.0),
    );

    let beta_sq = Expr::new_div(v_sq, c_sq);

    Expr::new_pow(
        Expr::new_sub(
            Expr::Constant(1.0),
            beta_sq,
        ),
        Expr::Constant(-0.5),
    )
}

/// Performs a Lorentz transformation in the x-direction.
///
/// Formulas:
/// $x' = \gamma (x - vt)$
/// $t' = \gamma (t - vx/c^2)$
#[must_use]

pub fn lorentz_transformation_x(
    x : &Expr,
    t : &Expr,
    v : &Expr,
) -> (Expr, Expr) {

    let gamma = lorentz_factor(v);

    let c = Expr::new_variable("c");

    let c_sq = Expr::new_pow(
        c,
        Expr::Constant(2.0),
    );

    let x_prime = Expr::new_mul(
        gamma.clone(),
        Expr::new_sub(
            x.clone(),
            Expr::new_mul(v.clone(), t.clone()),
        ),
    );

    let t_prime = Expr::new_mul(
        gamma,
        Expr::new_sub(
            t.clone(),
            Expr::new_div(
                Expr::new_mul(v.clone(), x.clone()),
                c_sq,
            ),
        ),
    );

    (x_prime, t_prime)
}

/// Calculates relativistic velocity addition.
///
/// Formula: $u = \frac{v + u'}{1 + v u'/c^2}$
#[must_use]

pub fn velocity_addition(
    v : &Expr,
    u_prime : &Expr,
) -> Expr {

    let c = Expr::new_variable("c");

    let c_sq = Expr::new_pow(
        c,
        Expr::Constant(2.0),
    );

    let numerator = Expr::new_add(
        v.clone(),
        u_prime.clone(),
    );

    let denominator = Expr::new_add(
        Expr::Constant(1.0),
        Expr::new_div(
            Expr::new_mul(
                v.clone(),
                u_prime.clone(),
            ),
            c_sq,
        ),
    );

    Expr::new_div(
        numerator,
        denominator,
    )
}

/// Calculates mass-energy equivalence: $E = mc^2$.
#[must_use]

pub fn mass_energy_equivalence(mass : &Expr) -> Expr {

    let c = Expr::new_variable("c");

    let c_sq = Expr::new_pow(
        c,
        Expr::Constant(2.0),
    );

    Expr::new_mul(mass.clone(), c_sq)
}

/// Calculates relativistic momentum: $p = \gamma m v$.
#[must_use]

pub fn relativistic_momentum(
    mass : &Expr,
    velocity : &Expr,
) -> Expr {

    let gamma = lorentz_factor(velocity);

    Expr::new_mul(
        gamma,
        Expr::new_mul(
            mass.clone(),
            velocity.clone(),
        ),
    )
}

/// Calculates the Relativistic Doppler Effect for source and observer moving apart.
///
/// Formula: $f_{obs} = f_{src} \sqrt{\frac{1 - \beta}{1 + \beta}}$ where $\beta = v/c$.
#[must_use]

pub fn doppler_effect(
    f_src : &Expr,
    v : &Expr,
) -> Expr {

    let c = Expr::new_variable("c");

    let beta = Expr::new_div(v.clone(), c);

    let ratio = Expr::new_div(
        Expr::new_sub(
            Expr::Constant(1.0),
            beta.clone(),
        ),
        Expr::new_add(
            Expr::Constant(1.0),
            beta,
        ),
    );

    Expr::new_mul(
        f_src.clone(),
        Expr::new_pow(
            ratio,
            Expr::Constant(0.5),
        ),
    )
}

/// Calculates the Schwarzschild Radius: $`r_s` = \frac{2GM}{c^2}$.
#[must_use]

pub fn schwarzschild_radius(mass : &Expr) -> Expr {

    let g = Expr::new_variable("G");

    let c = Expr::new_variable("c");

    let c_sq = Expr::new_pow(
        c,
        Expr::Constant(2.0),
    );

    Expr::new_div(
        Expr::new_mul(
            Expr::Constant(2.0),
            Expr::new_mul(g, mass.clone()),
        ),
        c_sq,
    )
}

/// Calculates gravitational time dilation in the Schwarzschild metric.
///
/// Formula: $t_{proper} = t_{far} \sqrt{1 - \`frac{r_s}{r`}}$
#[must_use]

pub fn gravitational_time_dilation(
    t_far : &Expr,
    r : &Expr,
    mass : &Expr,
) -> Expr {

    let r_s = schwarzschild_radius(mass);

    let ratio = Expr::new_div(r_s, r.clone());

    Expr::new_mul(
        t_far.clone(),
        Expr::new_pow(
            Expr::new_sub(
                Expr::Constant(1.0),
                ratio,
            ),
            Expr::Constant(0.5),
        ),
    )
}

/// Represents the simplified Einstein Field Equations (LHS).
///
/// Formula: $G_{\mu\nu} = R_{\mu\nu} - \frac{1}{2} R g_{\mu\nu}$
#[must_use]

pub fn einstein_tensor(
    ricci_tensor : &Expr,
    scalar_curvature : &Expr,
    metric : &Expr,
) -> Expr {

    Expr::new_sub(
        ricci_tensor.clone(),
        Expr::new_mul(
            Expr::Constant(0.5),
            Expr::new_mul(
                scalar_curvature.clone(),
                metric.clone(),
            ),
        ),
    )
}

/// Represents the Geodesic Equation term: $\frac{d^2 x^\mu}{d\tau^2} = -\Gamma^\mu_{\alpha\beta} \frac{dx^\alpha}{d\tau} \frac{dx^\beta}{d\tau}$
///
/// This function returns the right-hand side of the geodesic equation.
#[must_use]

pub fn geodesic_acceleration(
    christoffel : &Expr,
    u_alpha : &Expr,
    u_beta : &Expr,
) -> Expr {

    Expr::new_neg(Arc::new(
        Expr::new_mul(
            christoffel.clone(),
            Expr::new_mul(
                u_alpha.clone(),
                u_beta.clone(),
            ),
        ),
    ))
}

// --- Backward Compatibility Aliases ---

/// Alias for `lorentz_transformation_x`.
#[must_use]

pub fn lorentz_transformation(
    x : &Expr,
    t : &Expr,
    v : &Expr,
) -> (Expr, Expr) {

    lorentz_transformation_x(x, t, v)
}

/// Placeholder for `einstein_field_equations` (now use `einstein_tensor`).
#[must_use]

pub fn einstein_field_equations(
    ricci : &Expr,
    scalar : &Expr,
    metric : &Expr,
    _stress : &Expr,
) -> Expr {

    einstein_tensor(
        ricci,
        scalar,
        metric,
    )
}

/// Placeholder for `geodesic_equation` (now use `geodesic_acceleration`).
#[must_use]

pub fn geodesic_equation(
    christoffel : &Expr,
    u_alpha : &Expr,
    u_beta : &Expr,
    _tau : &str,
) -> Expr {

    geodesic_acceleration(
        christoffel,
        u_alpha,
        u_beta,
    )
}
