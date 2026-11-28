// Finite and bounded checked adapters for the dimensionless map z -> z² + c.
struct QuadraticEscape {
    steps: usize,
    escaped: bool,
    potential: Option<f64>,
}
record!(QuadraticEscape { steps: usize, escaped: bool, potential: Option<f64> });

struct QuadraticPeriodicPoints {
    points: Vec<num_complex::Complex64>,
    residuals: Vec<f64>,
    period: usize,
}
record!(QuadraticPeriodicPoints { points: Vec<num_complex::Complex64>, residuals: Vec<f64>, period: usize });

struct QuadraticMultiplier {
    multiplier: num_complex::Complex64,
    kind: String,
    period: usize,
    residual: f64,
}
record!(QuadraticMultiplier {
    multiplier: num_complex::Complex64,
    kind: String,
    period: usize,
    residual: f64
});

struct QuadraticRay {
    point: num_complex::Complex64,
    residual: f64,
}
record!(QuadraticRay {
    point: num_complex::Complex64,
    residual: f64
});

const QUADRATIC_MAX_STEPS: usize = 16_384;

fn quadratic_escape_budget(
    max_steps: usize,
    radius: f64,
) -> Result<(), OpError> {
    if max_steps > QUADRATIC_MAX_STEPS {
        return Err(OpError::BudgetExceeded);
    }
    if max_steps == 0 || radius < 2.0 {
        return Err(OpError::Domain(
            "escape iteration needs a positive budget and radius at least 2".into(),
        ));
    }
    crate::cancellation::checkpoint()?;
    Ok(())
}

fn quadratic_escape_result(
    steps: usize,
    potential: f64,
    max_steps: usize,
) -> Result<QuadraticEscape, OpError> {
    crate::cancellation::checkpoint()?;
    if !potential.is_finite() {
        return Err(OpError::Domain(
            "quadratic escape iterate is not finite".into(),
        ));
    }
    let escaped = steps < max_steps;
    Ok(QuadraticEscape {
        steps,
        escaped,
        potential: escaped.then_some(potential),
    })
}

fn quadratic_step(
    z: num_complex::Complex64,
    c: num_complex::Complex64,
) -> num_complex::Complex64 {
    if crate::cancellation::is_cancelled() {
        num_complex::Complex64::new(f64::NAN, f64::NAN)
    } else {
        z * z + c
    }
}

fn quadratic_dynamics_error(message: String) -> OpError {
    match crate::cancellation::checkpoint() {
        | Ok(()) => OpError::Domain(message),
        | Err(error) => error.into(),
    }
}

binding!(
    complex_dynamics_mandelbrot_escape,
    "num_complex_dynamics_mandelbrot_escape",
    "numerical_special",
    "`num_complex_dynamics_mandelbrot_escape(c, max_steps, radius)` checks iterates preceding the budget for a radius crossing. Returns the crossing step and a potential estimate, or the exhausted budget and no estimate.",
    "crate::kernels::complex_dynamics::mandelbrot_escape_time",
    DimensionPolicy::Dimensionless,
    (c: num_complex::Complex64, max_steps: usize, radius: f64) -> Result<super::QuadraticEscape, OpError>,
    |c: num_complex::Complex64, max_steps: usize, radius: f64| {
        super::quadratic_escape_budget(max_steps, radius)?;
        let (steps, potential) = crate::kernels::complex_dynamics::mandelbrot_escape_time(c, max_steps, radius);
        super::quadratic_escape_result(steps, potential, max_steps)
    }
);

binding!(
    complex_dynamics_julia_escape,
    "num_complex_dynamics_julia_escape",
    "numerical_special",
    "`num_complex_dynamics_julia_escape(initial, c, max_steps, radius)` checks iterates preceding the budget for a radius crossing. Returns the crossing step and a potential estimate, or the exhausted budget and no estimate.",
    "crate::kernels::complex_dynamics::julia_escape_time",
    DimensionPolicy::Dimensionless,
    (initial: num_complex::Complex64, c: num_complex::Complex64, max_steps: usize, radius: f64) -> Result<super::QuadraticEscape, OpError>,
    |initial: num_complex::Complex64, c: num_complex::Complex64, max_steps: usize, radius: f64| {
        super::quadratic_escape_budget(max_steps, radius)?;
        let (steps, potential) = crate::kernels::complex_dynamics::julia_escape_time(initial, c, max_steps, radius);
        super::quadratic_escape_result(steps, potential, max_steps)
    }
);

binding!(
    complex_dynamics_in_main_components,
    "num_complex_dynamics_in_main_components",
    "numerical_special",
    "`num_complex_dynamics_in_main_components(c)` tests membership in the Mandelbrot main cardioid or period-two bulb. Other bounded parameters may return false.",
    "crate::kernels::complex_dynamics::in_main_cardioid_or_period2_bulb",
    DimensionPolicy::Dimensionless,
    (c: num_complex::Complex64) -> bool,
    crate::kernels::complex_dynamics::in_main_cardioid_or_period2_bulb
);

binding!(
    complex_dynamics_mandelbrot_distance,
    "num_complex_dynamics_mandelbrot_distance",
    "numerical_special",
    "`num_complex_dynamics_mandelbrot_distance(c, max_steps)` returns an exterior distance estimate when an iterate crosses the native threshold, or None when no crossing was observed within the budget.",
    "crate::kernels::complex_dynamics::mandelbrot_distance_estimator",
    DimensionPolicy::Dimensionless,
    (c: num_complex::Complex64, max_steps: usize) -> Result<Option<f64>, OpError>,
    |c: num_complex::Complex64, max_steps: usize| {
        super::quadratic_escape_budget(max_steps, 2.0)?;
        let estimate = crate::kernels::complex_dynamics::mandelbrot_distance_estimator(c, max_steps);
        crate::cancellation::checkpoint()?;
        if !estimate.is_finite() {
            return Err(OpError::Domain("Mandelbrot distance estimate is outside binary64 range".into()));
        }
        Ok((estimate > 0.0).then_some(estimate))
    }
);

binding!(
    complex_dynamics_quadratic_periodic_points,
    "num_complex_dynamics_quadratic_periodic_points",
    "numerical_special",
    "`num_complex_dynamics_quadratic_periodic_points(c, period, seeds, max_steps, tolerance)` uses bounded Newton searches from the supplied seeds and returns deduplicated candidates for the requested primitive period with their residuals. Empty results report that the supplied searches found no candidate.",
    "crate::kernels::complex_dynamics::find_periodic_points_quadratic",
    DimensionPolicy::Dimensionless,
    (c: num_complex::Complex64, period: usize, seeds: Vec<num_complex::Complex64>, max_steps: usize, tolerance: f64) -> Result<super::QuadraticPeriodicPoints, OpError>,
    |c: num_complex::Complex64, period: usize, seeds: Vec<num_complex::Complex64>, max_steps: usize, tolerance: f64| {
        if period == 0 || max_steps == 0 || tolerance <= 0.0 || !(tolerance * 1e3).is_finite() {
            return Err(OpError::Domain("period and iterations must be positive, and the root tolerance and deduplication scale must be finite and positive".into()));
        }
        let work_per_seed = period.saturating_mul(max_steps.saturating_add(period));
        if seeds.len() > 2048 || period > 1024 || max_steps > super::QUADRATIC_MAX_STEPS
            || seeds.len().saturating_mul(work_per_seed.saturating_add(seeds.len())) > (1 << 23)
        {
            return Err(OpError::BudgetExceeded);
        }
        crate::cancellation::checkpoint()?;
        let points = crate::kernels::complex_dynamics::find_periodic_points_quadratic(c, period, &seeds, max_steps, tolerance);
        crate::cancellation::checkpoint()?;
        let residuals = points.iter().map(|&z| ((0..period).fold(z, |w, _| super::quadratic_step(w,c)) - z).norm()).collect();
        crate::cancellation::checkpoint()?;
        Ok(super::QuadraticPeriodicPoints { points, residuals, period })
    }
);

binding!(
    complex_dynamics_bottcher,
    "num_complex_dynamics_bottcher",
    "numerical_special",
    "`num_complex_dynamics_bottcher(point, c, max_steps)` returns a truncated principal Böttcher coordinate estimate in the exterior quadratic basin. Requires a positive iteration budget and an initial point outside the conservative radius sqrt(2) sqrt(max(abs(re(c)), abs(im(c)))) + 1, with a minimum radius of 2.",
    "crate::kernels::complex_dynamics::bottcher_coordinate",
    DimensionPolicy::Dimensionless,
    (point: num_complex::Complex64, c: num_complex::Complex64, max_steps: usize) -> Result<num_complex::Complex64, OpError>,
    |point: num_complex::Complex64, c: num_complex::Complex64, max_steps: usize| {
        if max_steps > 1024 { return Err(OpError::BudgetExceeded); }
        let exterior = crate::kernels::complex_dynamics::quadratic_exterior_radius(c);
        if max_steps == 0 || point.norm() <= exterior {
            return Err(OpError::Domain("Bottcher coordinate needs a positive budget and an initial point in the exterior basin".into()));
        }
        crate::kernels::complex_dynamics::bottcher_coordinate(point,c,max_steps).map_err(super::quadratic_dynamics_error)
    }
);

binding!(
    complex_dynamics_external_ray,
    "num_complex_dynamics_external_ray",
    "numerical_special",
    "`num_complex_dynamics_external_ray(c, angle, potential, max_steps)` estimates a quadratic external-ray point and reports a relative Böttcher-coordinate residual. Angles are turns in [0, 1); the potential is positive and at most 325.",
    "crate::kernels::complex_dynamics::external_ray_point",
    DimensionPolicy::Dimensionless,
    (c: num_complex::Complex64, angle: f64, potential: f64, max_steps: usize) -> Result<super::QuadraticRay, OpError>,
    |c: num_complex::Complex64, angle: f64, potential: f64, max_steps: usize| {
        if max_steps > 30 { return Err(OpError::BudgetExceeded); }
        if !(0.0..1.0).contains(&angle) { return Err(OpError::Domain("external ray angle must lie in [0, 1) turns".into())); }
        crate::cancellation::checkpoint()?;
        let point = crate::kernels::complex_dynamics::external_ray_point(c,angle,potential,max_steps).map_err(super::quadratic_dynamics_error)?;
        let at_point = crate::kernels::complex_dynamics::bottcher_coordinate(point,c,24).map_err(super::quadratic_dynamics_error)?;
        let modulus = potential.exp();
        let target = num_complex::Complex64::from_polar(modulus,std::f64::consts::TAU*angle);
        let residual = ((at_point-target)/modulus).norm();
        crate::cancellation::checkpoint()?;
        Ok(super::QuadraticRay { point, residual })
    }
);

binding!(
    complex_dynamics_quadratic_orbit,
    "num_complex_dynamics_quadratic_orbit",
    "numerical_special",
    "`num_complex_dynamics_quadratic_orbit(initial, c, steps)` returns the dimensionless forward orbit of z -> z²+c, including the initial point. A nonfinite iterate produces an error.",
    "crate::kernels::complex_dynamics::orbit",
    DimensionPolicy::Dimensionless,
    (initial: num_complex::Complex64, c: num_complex::Complex64, steps: usize) -> Result<Vec<num_complex::Complex64>, OpError>,
    |initial: num_complex::Complex64, c: num_complex::Complex64, steps: usize| {
        if steps > super::QUADRATIC_MAX_STEPS { return Err(OpError::BudgetExceeded); }
        crate::cancellation::checkpoint()?;
        let orbit = crate::kernels::complex_dynamics::orbit(|z|super::quadratic_step(z,c),initial,steps);
        crate::cancellation::checkpoint()?;
        Ok(orbit)
    }
);

binding!(
    complex_dynamics_quadratic_multiplier,
    "num_complex_dynamics_quadratic_multiplier",
    "numerical_special",
    "`num_complex_dynamics_quadratic_multiplier(point, c, period, tolerance)` returns the multiplier, closure residual and numerical class for the requested quadratic iterate, or an error when it fails to close within tolerance.",
    "crate::kernels::complex_dynamics::classify_periodic_point",
    DimensionPolicy::Dimensionless,
    (point: num_complex::Complex64, c: num_complex::Complex64, period: usize, tolerance: f64) -> Result<super::QuadraticMultiplier, OpError>,
    |point: num_complex::Complex64, c: num_complex::Complex64, period: usize, tolerance: f64| {
        use crate::kernels::complex_dynamics::FixedPointKind;
        if period > 1024 { return Err(OpError::BudgetExceeded); }
        if period == 0 || tolerance <= 0.0 { return Err(OpError::Domain("quadratic multiplier needs a positive period and closure tolerance".into())); }
        let last = std::cell::Cell::new(point);
        let (multiplier, kind) = crate::kernels::complex_dynamics::classify_periodic_point(|z| { let w=super::quadratic_step(z,c); last.set(w); w }, |z| z*2.0,point,period);
        crate::cancellation::checkpoint()?;
        let residual = (last.get() - point).norm();
        if !residual.is_finite() || residual > tolerance {
            return Err(OpError::Domain("quadratic orbit does not close within the requested tolerance".into()));
        }
        let kind = match kind {
            | FixedPointKind::Superattracting => "superattracting",
            | FixedPointKind::Attracting => "attracting",
            | FixedPointKind::Neutral => "neutral",
            | FixedPointKind::Repelling => "repelling",
        };
        Ok(super::QuadraticMultiplier { multiplier, kind: kind.into(), period, residual })
    }
);

binding!(
    complex_dynamics_quadratic_lyapunov,
    "num_complex_dynamics_quadratic_lyapunov",
    "numerical_special",
    "`num_complex_dynamics_quadratic_lyapunov(initial, c, steps, burn_in)` estimates the mean log derivative magnitude along the dimensionless quadratic orbit after the requested burn-in. A zero or nonfinite derivative is an error.",
    "crate::kernels::complex_dynamics::lyapunov_exponent",
    DimensionPolicy::Dimensionless,
    (initial: num_complex::Complex64, c: num_complex::Complex64, steps: usize, burn_in: usize) -> Result<f64, OpError>,
    |initial: num_complex::Complex64, c: num_complex::Complex64, steps: usize, burn_in: usize| {
        if steps.saturating_add(burn_in) > super::QUADRATIC_MAX_STEPS { return Err(OpError::BudgetExceeded); }
        crate::cancellation::checkpoint()?;
        let result = crate::kernels::complex_dynamics::lyapunov_exponent(|z|super::quadratic_step(z,c),|z|z*2.0,initial,steps,burn_in).map_err(super::quadratic_dynamics_error)?;
        crate::cancellation::checkpoint()?;
        Ok(result)
    }
);
