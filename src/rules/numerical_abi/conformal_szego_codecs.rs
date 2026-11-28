// Portable warm-start state for Szegő collocation.
use crate::kernels::conformal::SzegoKerzmanSteinSolver;

struct SzegoWire {
    points: Vec<Complex64>,
    weights: Vec<Complex64>,
    center: Complex64,
}
record!(SzegoWire { points: Vec<Complex64>, weights: Vec<Complex64>, center: Complex64 });
struct SzegoMapValue {
    mapped: Complex64,
    szego: Complex64,
}
record!(SzegoMapValue {
    mapped: Complex64,
    szego: Complex64
});

impl Arg for SzegoKerzmanSteinSolver {
    fn ty() -> TypeExpr {
        SzegoWire::ty()
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        const FIELDS: &[&str] = &["points", "weights", "center"];
        let field = |name| model_field(value, name, FIELDS);
        let points = Vec::<Complex64>::decode(field("points")?, index)?;
        let weights = Vec::<Complex64>::decode(field("weights")?, index)?;
        let center = Complex64::decode(field("center")?, index)?;
        Self::from_weights(points, weights, center).map_err(OpError::Domain)
    }
}

impl Output for SzegoKerzmanSteinSolver {
    fn ty() -> TypeExpr {
        SzegoWire::ty()
    }

    fn encode(self) -> Result<Value, OpError> {
        SzegoWire {
            points: self.points,
            weights: self.weights,
            center: self.center,
        }
        .encode()
    }
}

pub(in crate::rules) fn szego_model_type() -> TypeExpr {
    SzegoWire::ty()
}

binding!(
    conformal_szego_prepare_sampled,
    "num_conformal_szego_prepare_sampled",
    "numerical_special",
    "Prepares reusable Szegő quadrature from 16 to 192 equally spaced counterclockwise boundary samples and their complex parameter derivatives, in dimensionless normalized coordinates.",
    "crate::kernels::conformal::SzegoKerzmanSteinSolver::solve",
    DimensionPolicy::Dimensionless,
    (points: Vec<num_complex::Complex64>, derivatives: Vec<num_complex::Complex64>, center: num_complex::Complex64) -> Result<super::SzegoKerzmanSteinSolver, OpError>,
    |points: Vec<num_complex::Complex64>, derivatives: Vec<num_complex::Complex64>, center: num_complex::Complex64| {
        let result=super::SzegoKerzmanSteinSolver::prepare(points,&derivatives,center);
        crate::cancellation::checkpoint()?;
        result.map_err(OpError::Domain)
    }
);
binding!(
    conformal_szego_eval,
    "num_conformal_szego_evaluate",
    "numerical_special",
    "Evaluates an interior point against a portable prepared Szegő quadrature.",
    "crate::kernels::conformal::SzegoKerzmanSteinSolver::eval_szego",
    DimensionPolicy::Dimensionless,
    (model: super::SzegoKerzmanSteinSolver, point: num_complex::Complex64) -> Result<num_complex::Complex64, OpError>,
    |model: super::SzegoKerzmanSteinSolver, point: num_complex::Complex64| {
        model.checked_point(point).map_err(OpError::Domain)?;
        Ok(model.eval_szego(point))
    }
);
binding!(
    conformal_szego_eval_map,
    "num_conformal_szego_evaluate_map",
    "numerical_special",
    "Evaluates a Riemann map and Szegő kernel from prepared quadrature. The straight segment from the center to the query must stay interior.",
    "crate::kernels::conformal::SzegoKerzmanSteinSolver::eval_map",
    DimensionPolicy::Dimensionless,
    (model: super::SzegoKerzmanSteinSolver, point: num_complex::Complex64) -> Result<super::SzegoMapValue, OpError>,
    |model: super::SzegoKerzmanSteinSolver, point: num_complex::Complex64| {
        let result=model.eval_map_checked(point);
        crate::cancellation::checkpoint()?;
        let (mapped,szego)=result.map_err(OpError::Domain)?;
        Ok(super::SzegoMapValue {mapped,szego})
    }
);
binding!(
    conformal_szego_eval_bergman,
    "num_conformal_szego_evaluate_bergman",
    "numerical_special",
    "Evaluates the Bergman kernel via the prepared Szegő quadrature at an interior point.",
    "crate::kernels::conformal::SzegoKerzmanSteinSolver::eval_bergman",
    DimensionPolicy::Dimensionless,
    (model: super::SzegoKerzmanSteinSolver, point: num_complex::Complex64) -> Result<num_complex::Complex64, OpError>,
    |model: super::SzegoKerzmanSteinSolver, point: num_complex::Complex64| {
        model.checked_point(point).map_err(OpError::Domain)?;
        Ok(model.eval_bergman(point))
    }
);
