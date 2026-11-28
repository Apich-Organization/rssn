// Quadratic and linear optimization with independent physical reference domains.
analysis_enum_arg!(crate::kernels::optim::Relation;Le,Ge,Eq);
analysis_enum_arg!(crate::kernels::optim::Pricing;Dantzig,SteepestEdge,Bland);
impl Output for crate::kernels::optim::Relation {
    fn ty() -> TypeExpr {
        TypeExpr::String
    }
    fn encode(self) -> Result<Value, OpError> {
        format!("{self:?}").encode()
    }
}
native_vector_codec!(crate::kernels::optim::Relation);
record!(crate::kernels::optim::QpSolution{d:Vec<f64>,lambda:Vec<f64>,mu:Vec<f64>});
struct AnalysisLpOutput {
    status: String,
    x: Vec<f64>,
    value: Option<f64>,
    iterations: Option<usize>,
    bland_pivots: Option<usize>,
}
record!(AnalysisLpOutput{status:String,x:Vec<f64>,value:Option<f64>,iterations:Option<usize>,bland_pivots:Option<usize>});
fn analysis_lp_output(
    status: crate::kernels::optim::LpStatus,
    xs: f64,
    ys: f64,
    stats: Option<(usize, usize)>,
) -> AnalysisLpOutput {
    let (status, x, value) = match status {
        | crate::kernels::optim::LpStatus::Optimal { x, value } => (
            "Optimal",
            x.into_iter().map(|x| x * xs).collect(),
            Some(value * ys),
        ),
        | crate::kernels::optim::LpStatus::Infeasible => ("Infeasible", vec![], None),
        | crate::kernels::optim::LpStatus::Unbounded => ("Unbounded", vec![], None),
    };
    AnalysisLpOutput {
        status: status.into(),
        x,
        value,
        iterations: stats.map(|s| s.0),
        bland_pivots: stats.map(|s| s.1),
    }
}
fn analysis_optimization_shape(
    n: usize,
    rows: &[&Vec<Vec<f64>>],
    rhs: &[&Vec<f64>],
    references: &[f64],
) -> Result<(), OpError> {
    if n == 0 || n > 128 || references.iter().any(|s| !s.is_finite() || *s <= 0.0) {
        return Err(OpError::Domain(
            "optimization needs1..128 variables and positive finite reference scales".into(),
        ));
    }
    for (a, b) in rows.iter().zip(rhs) {
        if a.len() != b.len() || a.len() > 256 || a.iter().any(|r| r.len() != n) {
            return Err(OpError::Domain(
                "constraint shape or256row budget exceeded".into(),
            ));
        }
    }
    Ok(())
}
binding!(analysis_qp,"num_analysis_qp_active_set_scaled","numerical_vector","Native active-set QP with positive coordinate,objective,equality and inequality references. H has E/D^2,g E/D,Aeq C/D,beq C,Aineq G/D,bineq G; result step D and multipliers E/C,E/G.","crate::kernels::optim::qp_active_set",DimensionPolicy::AnalysisQp,(hess:dense::Mat,gradient:Vec<f64>,equalities:Vec<Vec<f64>>,equality_values:Vec<f64>,inequalities:Vec<Vec<f64>>,inequality_values:Vec<f64>,coordinate_reference:f64,objective_reference:f64,equality_reference:f64,inequality_reference:f64)->Result<Option<crate::kernels::optim::QpSolution>,OpError>,|mut hess:dense::Mat,gradient:Vec<f64>,equalities:Vec<Vec<f64>>,equality_values:Vec<f64>,inequalities:Vec<Vec<f64>>,inequality_values:Vec<f64>,coordinate_reference:f64,objective_reference:f64,equality_reference:f64,inequality_reference:f64|{
 super::analysis_optimization_shape(gradient.len(),&[&equalities,&inequalities],&[&equality_values,&inequality_values],&[coordinate_reference,objective_reference,equality_reference,inequality_reference])?;if hess.rows!=gradient.len()||hess.cols!=gradient.len(){return Err(OpError::Domain("QP Hessian shape mismatch".into()));}
 for h in&mut hess.data{*h=(*h*coordinate_reference/objective_reference)*coordinate_reference;}
 let g=gradient.into_iter().map(|v|v*(coordinate_reference/objective_reference)).collect::<Vec<_>>();let ae=equalities.into_iter().map(|r|r.into_iter().map(|v|v*(coordinate_reference/equality_reference)).collect()).collect::<Vec<_>>();let be=equality_values.into_iter().map(|v|v/equality_reference).collect::<Vec<_>>();let ai=inequalities.into_iter().map(|r|r.into_iter().map(|v|v*(coordinate_reference/inequality_reference)).collect()).collect::<Vec<_>>();let bi=inequality_values.into_iter().map(|v|v/inequality_reference).collect::<Vec<_>>();
 Ok(crate::kernels::optim::qp_active_set(&hess,&g,&ae,&be,&ai,&bi).map(|mut r|{for v in&mut r.d{*v*=coordinate_reference;}for v in&mut r.lambda{*v*=objective_reference/equality_reference;}for v in&mut r.mu{*v*=objective_reference/inequality_reference;}r}))
});
binding!(analysis_simplex,"num_analysis_simplex_scaled","numerical_vector","Native simplex LP with normalized coefficient conversion using positive coordinate/objective/constraint references; preserves Optimal,Infeasible,Unbounded without fabricated coordinates or objective.","crate::kernels::optim::simplex",DimensionPolicy::AnalysisSimplex,(objective:Vec<f64>,constraints:Vec<Vec<f64>>,relations:Vec<crate::kernels::optim::Relation>,values:Vec<f64>,coordinate_reference:f64,objective_reference:f64,constraint_reference:f64)->Result<super::AnalysisLpOutput,OpError>,|objective:Vec<f64>,constraints:Vec<Vec<f64>>,relations:Vec<crate::kernels::optim::Relation>,values:Vec<f64>,coordinate_reference:f64,objective_reference:f64,constraint_reference:f64|{
 super::analysis_optimization_shape(objective.len(),&[&constraints],&[&values],&[coordinate_reference,objective_reference,constraint_reference])?;if relations.len()!=values.len(){return Err(OpError::Domain("LP relation count mismatch".into()));}
 let c=objective.into_iter().map(|v|v*(coordinate_reference/objective_reference)).collect::<Vec<_>>();let a=constraints.into_iter().map(|r|r.into_iter().map(|v|v*(coordinate_reference/constraint_reference)).collect()).collect::<Vec<_>>();let b=values.into_iter().map(|v|v/constraint_reference).collect::<Vec<_>>();Ok(super::analysis_lp_output(crate::kernels::optim::simplex(&c,&a,&relations,&b),coordinate_reference,objective_reference,None))
});
binding!(analysis_simplex_with,"num_analysis_simplex_with_scaled","numerical_vector","Native simplex LP with actual selectable Dantzig,SteepestEdge,Bland pricing and retained pivot statistics; independent positive coordinate/objective/constraint references.","crate::kernels::optim::simplex_with",DimensionPolicy::AnalysisSimplexWith,(objective:Vec<f64>,constraints:Vec<Vec<f64>>,relations:Vec<crate::kernels::optim::Relation>,values:Vec<f64>,pricing:crate::kernels::optim::Pricing,coordinate_reference:f64,objective_reference:f64,constraint_reference:f64)->Result<super::AnalysisLpOutput,OpError>,|objective:Vec<f64>,constraints:Vec<Vec<f64>>,relations:Vec<crate::kernels::optim::Relation>,values:Vec<f64>,pricing:crate::kernels::optim::Pricing,coordinate_reference:f64,objective_reference:f64,constraint_reference:f64|{
 super::analysis_optimization_shape(objective.len(),&[&constraints],&[&values],&[coordinate_reference,objective_reference,constraint_reference])?;if relations.len()!=values.len(){return Err(OpError::Domain("LP relation count mismatch".into()));}
 let c=objective.into_iter().map(|v|v*(coordinate_reference/objective_reference)).collect::<Vec<_>>();let a=constraints.into_iter().map(|r|r.into_iter().map(|v|v*(coordinate_reference/constraint_reference)).collect()).collect::<Vec<_>>();let b=values.into_iter().map(|v|v/constraint_reference).collect::<Vec<_>>();let r=crate::kernels::optim::simplex_with(&c,&a,&relations,&b,pricing);Ok(super::analysis_lp_output(r.status,coordinate_reference,objective_reference,Some((r.iterations,r.bland_pivots))))
});
binding!(analysis_kkt,"num_analysis_kkt_residual_scaled","numerical_vector","Actual native KKT stationarity residual in the normalized coordinate/objective metric. Physical gradient E/D,constraint Jacobians C/D,G/D,multipliers E/C,E/G are normalized by four positive references; result is dimensionless.","crate::kernels::optim::kkt_residual",DimensionPolicy::AnalysisKkt,(gradient:Vec<f64>,equality_jacobian:Vec<Vec<f64>>,inequality_jacobian:Vec<Vec<f64>>,lambda:Vec<f64>,mu:Vec<f64>,coordinate_reference:f64,objective_reference:f64,equality_reference:f64,inequality_reference:f64)->Result<f64,OpError>,|gradient:Vec<f64>,equality_jacobian:Vec<Vec<f64>>,inequality_jacobian:Vec<Vec<f64>>,lambda:Vec<f64>,mu:Vec<f64>,coordinate_reference:f64,objective_reference:f64,equality_reference:f64,inequality_reference:f64|{
 super::analysis_optimization_shape(gradient.len(),&[&equality_jacobian,&inequality_jacobian],&[&lambda,&mu],&[coordinate_reference,objective_reference,equality_reference,inequality_reference])?;
 let scale=crate::rules::numerical_abi::native_object_api::scale_product;
 let g=gradient.into_iter().map(|v|scale(v,&[coordinate_reference],&[objective_reference])).collect::<Result<Vec<_>,_>>()?;
 let jh=equality_jacobian.into_iter().map(|r|r.into_iter().map(|v|scale(v,&[coordinate_reference],&[equality_reference])).collect::<Result<Vec<_>,_>>()).collect::<Result<Vec<_>,_>>()?;
 let jg=inequality_jacobian.into_iter().map(|r|r.into_iter().map(|v|scale(v,&[coordinate_reference],&[inequality_reference])).collect::<Result<Vec<_>,_>>()).collect::<Result<Vec<_>,_>>()?;
 let lambda=lambda.into_iter().map(|v|scale(v,&[equality_reference],&[objective_reference])).collect::<Result<Vec<_>,_>>()?;let mu=mu.into_iter().map(|v|scale(v,&[inequality_reference],&[objective_reference])).collect::<Result<Vec<_>,_>>()?;
 Ok(crate::kernels::optim::kkt_residual(&g,&jh,&jg,&lambda,&mu))
});
