// Shared measured normalization and status conversion for the LP interfaces.
use super::native_object_api::scale_without_underflow as scale_lp_value;
use crate::kernels::convex::lp as bounded_lp;
use crate::kernels::optim as tableau_lp;

#[derive(Clone, Copy)]
struct LpFrame {
    coordinates: f64,
    objective: f64,
    constraints: f64,
}

struct NormalizedLp {
    objective: Vec<f64>,
    rows: Vec<Vec<f64>>,
    rhs: Vec<f64>,
}

impl LpFrame {
    fn normalize(
        self,
        mut objective: Vec<f64>,
        mut rows: Vec<Vec<f64>>,
        mut rhs: Vec<f64>,
        relations: &[tableau_lp::Relation],
    ) -> Result<NormalizedLp, OpError> {
        analysis_optimization_shape(
            objective.len(),
            &[&rows],
            &[&rhs],
            &[self.coordinates, self.objective, self.constraints],
        )?;
        if relations.len() != rhs.len() {
            return Err(OpError::Domain(
                "LP relation count must match the constraint rows".into(),
            ));
        }
        for value in &mut objective {
            *value = scale_lp_value(*value, &[self.coordinates], &[self.objective])?;
        }
        for value in rows.iter_mut().flatten() {
            *value = scale_lp_value(*value, &[self.coordinates], &[self.constraints])?;
        }
        for value in &mut rhs {
            *value = scale_lp_value(*value, &[], &[self.constraints])?;
        }
        Ok(NormalizedLp {
            objective,
            rows,
            rhs,
        })
    }

    fn output(
        self,
        result: tableau_lp::LpStatus,
        iterations: Option<usize>,
        bland_pivots: Option<usize>,
    ) -> Result<AnalysisLpOutput, OpError> {
        let (status, x, value) = match result {
            tableau_lp::LpStatus::Optimal { mut x, value } => {
                for coordinate in &mut x {
                    *coordinate = scale_lp_value(*coordinate, &[self.coordinates], &[])?;
                }
                (
                    "Optimal",
                    x,
                    Some(scale_lp_value(value, &[self.objective], &[])?),
                )
            }
            tableau_lp::LpStatus::Infeasible => ("Infeasible", Vec::new(), None),
            tableau_lp::LpStatus::Unbounded => ("Unbounded", Vec::new(), None),
        };
        Ok(AnalysisLpOutput {
            status: status.into(),
            x,
            value,
            iterations,
            bland_pivots,
        })
    }

    fn bounds(
        self,
        n: usize,
        mut lower: Vec<f64>,
        upper: Vec<f64>,
        upper_present: Vec<bool>,
    ) -> Result<(Vec<f64>, Vec<Option<f64>>), OpError> {
        if lower.len() != n || upper.len() != n || upper_present.len() != n {
            return Err(OpError::Domain(
                "LP bound vectors must match the objective length".into(),
            ));
        }
        for value in &mut lower {
            *value = scale_lp_value(*value, &[], &[self.coordinates])?;
        }
        let upper = upper
            .into_iter()
            .zip(upper_present)
            .map(|(value, supplied)| {
                supplied
                    .then(|| scale_lp_value(value, &[], &[self.coordinates]))
                    .transpose()
            })
            .collect::<Result<_, _>>()?;
        Ok((lower, upper))
    }
}

fn require_lp_budget(max_pivots: usize) -> Result<(), OpError> {
    if max_pivots > 100_000 {
        return Err(OpError::BudgetExceeded);
    }
    crate::cancellation::checkpoint()?;
    Ok(())
}

fn lp_solver_error(error: String) -> OpError {
    if let Err(cancelled) = crate::cancellation::checkpoint() {
        return cancelled.into();
    }
    if error == "iteration limit exceeded" {
        OpError::BudgetExceeded
    } else {
        OpError::Domain(error)
    }
}

binding!(
    analysis_simplex,
    "num_analysis_simplex_scaled",
    "numerical_vector",
    "Two-phase nonnegative simplex with measured coordinate, objective and constraint references.",
    "crate::kernels::optim::simplex",
    DimensionPolicy::AnalysisSimplex,
    (objective: Vec<f64>, constraints: Vec<Vec<f64>>, relations: Vec<super::tableau_lp::Relation>, values: Vec<f64>,
     coordinate_reference: f64, objective_reference: f64, constraint_reference: f64) -> Result<super::AnalysisLpOutput, OpError>,
    |objective, constraints, relations: Vec<super::tableau_lp::Relation>, values, coordinates, objective_scale, rhs_scale| {
        let frame = super::LpFrame { coordinates, objective: objective_scale, constraints: rhs_scale };
        let input = frame.normalize(objective, constraints, values, &relations)?;
        let result = super::tableau_lp::simplex(&input.objective, &input.rows, &relations, &input.rhs);
        crate::cancellation::checkpoint()?;
        frame.output(result, None, None)
    }
);

binding!(
    analysis_simplex_with,
    "num_analysis_simplex_with_scaled",
    "numerical_vector",
    "Two-phase nonnegative simplex with measured reference scales, selectable pricing and pivot statistics.",
    "crate::kernels::optim::simplex_with",
    DimensionPolicy::AnalysisSimplexWith,
    (objective: Vec<f64>, constraints: Vec<Vec<f64>>, relations: Vec<super::tableau_lp::Relation>, values: Vec<f64>,
     pricing: super::tableau_lp::Pricing, coordinate_reference: f64, objective_reference: f64,
     constraint_reference: f64) -> Result<super::AnalysisLpOutput, OpError>,
    |objective, constraints, relations: Vec<super::tableau_lp::Relation>, values, pricing, coordinates, objective_scale, rhs_scale| {
        let frame = super::LpFrame { coordinates, objective: objective_scale, constraints: rhs_scale };
        let input = frame.normalize(objective, constraints, values, &relations)?;
        let result = super::tableau_lp::simplex_with(&input.objective, &input.rows, &relations, &input.rhs, pricing);
        crate::cancellation::checkpoint()?;
        frame.output(result.status, Some(result.iterations), Some(result.bland_pivots))
    }
);

binding!(
    analysis_simplex_budgeted,
    "num_analysis_simplex_budgeted",
    "numerical_vector",
    "Two-phase nonnegative simplex with selectable Dantzig, SteepestEdge or Bland pricing and an exact total pivot budget.",
    "crate::kernels::optim::simplex_with_max_iter",
    DimensionPolicy::AnalysisSimplexWith,
    (objective: Vec<f64>, constraints: Vec<Vec<f64>>, relations: Vec<super::tableau_lp::Relation>, values: Vec<f64>,
     pricing: super::tableau_lp::Pricing, coordinate_reference: f64, objective_reference: f64,
     constraint_reference: f64, max_pivots: usize) -> Result<super::AnalysisLpOutput, OpError>,
    |objective, constraints, relations: Vec<super::tableau_lp::Relation>, values, pricing, coordinates, objective_scale, rhs_scale, max_pivots| {
        super::require_lp_budget(max_pivots)?;
        let frame = super::LpFrame { coordinates, objective: objective_scale, constraints: rhs_scale };
        let input = frame.normalize(objective, constraints, values, &relations)?;
        let result = super::tableau_lp::simplex_with_max_iter(&input.objective, &input.rows, &relations, &input.rhs, pricing, max_pivots)
            .map_err(super::lp_solver_error)?;
        frame.output(result.status, Some(result.iterations), Some(result.bland_pivots))
    }
);

binding!(
    analysis_bounded_lp,
    "num_convex_solve_lp",
    "numerical_vector",
    "Bounded dense LP with physical reference scales, finite lower bounds, masked upper bounds and a total simplex pivot budget.",
    "crate::kernels::convex::lp::solve_lp",
    DimensionPolicy::AnalysisBoundedLp,
    (objective: Vec<f64>, constraints: Vec<Vec<f64>>, relations: Vec<super::tableau_lp::Relation>, values: Vec<f64>,
     coordinate_reference: f64, objective_reference: f64, constraint_reference: f64, lower: Vec<f64>,
     upper: Vec<f64>, upper_present: Vec<bool>, max_pivots: usize) -> Result<super::AnalysisLpOutput, OpError>,
    |objective, constraints, relations: Vec<super::tableau_lp::Relation>, values, coordinates, objective_scale, rhs_scale, lower, upper, upper_present, max_pivots| {
        super::require_lp_budget(max_pivots)?;
        let frame = super::LpFrame { coordinates, objective: objective_scale, constraints: rhs_scale };
        let input = frame.normalize(objective, constraints, values, &relations)?;
        let (lower, upper) = frame.bounds(input.objective.len(), lower, upper, upper_present)?;
        let senses = relations.into_iter().map(|relation| match relation {
            | super::tableau_lp::Relation::Le => super::bounded_lp::Relation::Le,
            | super::tableau_lp::Relation::Ge => super::bounded_lp::Relation::Ge,
            | super::tableau_lp::Relation::Eq => super::bounded_lp::Relation::Eq,
        }).collect();
        let result = super::bounded_lp::solve_lp(&super::bounded_lp::LpProblem {
            c: input.objective, a: input.rows, rel: senses, b: input.rhs, lower, upper,
        }, max_pivots).map_err(super::lp_solver_error)?;
        let outcome = match result.status {
            | super::bounded_lp::LpStatus::Optimal => super::tableau_lp::LpStatus::Optimal { x: result.x, value: result.objective },
            | super::bounded_lp::LpStatus::Infeasible => super::tableau_lp::LpStatus::Infeasible,
            | super::bounded_lp::LpStatus::Unbounded => super::tableau_lp::LpStatus::Unbounded,
        };
        frame.output(outcome, Some(result.iterations), None)
    }
);
