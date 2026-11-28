// Normalized multivariate polynomial systems retain physical root/residual references.
struct AnalysisHomotopyResiduals(Vec<Option<f64>>);
impl Output for AnalysisHomotopyResiduals {
    fn ty() -> TypeExpr {
        TypeExpr::Tensor {
            rank: NatExpr::Lit(1),
            shape: None,
            elem: Box::new(<Option<f64> as Output>::ty()),
        }
    }

    fn encode(self) -> Result<Value, OpError> {
        Ok(Value::Tensor(TensorValue {
            shape: vec![self.0.len()],
            data: self
                .0
                .into_iter()
                .map(Output::encode)
                .collect::<Result<_, _>>()?,
            elem: <Option<f64> as Output>::ty(),
        }))
    }
}
struct AnalysisHomotopyOutput {
    roots: Vec<num_complex::Complex64>,
    offsets: Vec<usize>,
    finite: Vec<bool>,
    residuals: AnalysisHomotopyResiduals,
}
record!(AnalysisHomotopyOutput{roots:Vec<num_complex::Complex64>,offsets:Vec<usize>,finite:Vec<bool>,residuals:AnalysisHomotopyResiduals});
#[allow(
    clippy::needless_pass_by_value,
    reason = "The typed native ABI transfers decoded arguments by value; this shared dispatcher mirrors its explicit source signature."
)]
fn analysis_polynomial_system(
    coefficients: Vec<f64>,
    exponents: Vec<Vec<u32>>,
    offsets: Vec<usize>,
) -> Result<Vec<Vec<crate::kernels::rootfind::PolyTerm>>, OpError> {
    let variables = offsets
        .len()
        .checked_sub(1)
        .ok_or_else(|| OpError::Domain("polynomial offsets require at least two entries".into()))?;
    if variables == 0
        || variables > 6
        || coefficients.len() > 4096
        || coefficients.len() != exponents.len()
        || offsets.first() != Some(&0)
        || offsets.last() != Some(&coefficients.len())
        || offsets.windows(2).any(|w| w[0] >= w[1])
    {
        return Err(OpError::Domain(
            "invalid square polynomial system descriptor or6variable/4096term budget".into(),
        ));
    }
    let equations = offsets.windows(2).map(|range| {
        let (start, end) = (range[0], range[1]);
        coefficients[start..end]
            .iter()
            .copied()
            .zip(exponents[start..end].iter().map(Vec::as_slice))
    });
    crate::kernels::rootfind::total_degree_plan(variables, equations, 16, 256).map_err(|_| {
        OpError::Domain("invalid homotopy terms, degree above 16, or more than 256 paths".into())
    })?;
    let mut terms = coefficients
        .into_iter()
        .zip(exponents)
        .map(|(coef, exps)| crate::kernels::rootfind::PolyTerm { coef, exps });
    Ok(offsets
        .windows(2)
        .map(|w| terms.by_ref().take(w[1] - w[0]).collect())
        .collect())
}
fn analysis_homotopy_output(
    paths: Vec<crate::kernels::rootfind::PathResult>,
    coordinate_reference: f64,
    value_reference: f64,
) -> Result<AnalysisHomotopyOutput, OpError> {
    let mut roots = Vec::new();
    let mut offsets = vec![0];
    let mut finite = Vec::new();
    let mut residuals = Vec::new();
    for path in paths {
        finite.push(path.finite);
        residuals.push(if path.residual.is_finite() {
            Some(
                crate::rules::numerical_abi::native_object_api::scale_product(
                    path.residual,
                    &[value_reference],
                    &[],
                )?,
            )
        } else {
            None
        });
        if path.finite {
            for (re, im) in path.x {
                roots.push(num_complex::Complex64::new(
                    crate::rules::numerical_abi::native_object_api::scale_product(
                        re,
                        &[coordinate_reference],
                        &[],
                    )?,
                    crate::rules::numerical_abi::native_object_api::scale_product(
                        im,
                        &[coordinate_reference],
                        &[],
                    )?,
                ));
            }
        }
        offsets.push(roots.len());
    }
    Ok(AnalysisHomotopyOutput {
        roots,
        offsets,
        finite,
        residuals: AnalysisHomotopyResiduals(residuals),
    })
}
binding!(analysis_homotopy_system,"num_analysis_homotopy_system_scaled","numerical_vector","Actual native square polynomial homotopy for normalized coefficients/exponent rows/checked polynomial offsets. Positive coordinate/output references restore physical complex roots/residuals. Maximum6variables,4096terms,degree16,and256total-degree paths. Failed paths retain false status and optional residual; their undefined roots occupy empty offset ranges.","crate::kernels::rootfind::homotopy_system",DimensionPolicy::AnalysisHomotopy,(coefficients:Vec<f64>,exponents:Vec<Vec<u32>>,offsets:Vec<usize>,coordinate_reference:f64,value_reference:f64)->Result<super::AnalysisHomotopyOutput,OpError>,|coefficients:Vec<f64>,exponents:Vec<Vec<u32>>,offsets:Vec<usize>,coordinate_reference:f64,value_reference:f64|{if coordinate_reference<=0.0||value_reference<=0.0{return Err(OpError::Domain("homotopy references must be positive".into()));}let system=super::analysis_polynomial_system(coefficients,exponents,offsets)?;let paths=crate::kernels::rootfind::homotopy_system(&system).map_err(|e|OpError::Domain(format!("{e:?}")))?;super::analysis_homotopy_output(paths,coordinate_reference,value_reference)});
binding!(analysis_homotopy_gamma,"num_analysis_homotopy_system_with_gamma_scaled","numerical_vector","Actual gamma-selected native homotopy using normalized polynomial descriptors and independent physical root/residual references; all path statuses and defined coordinates are retained.","crate::kernels::rootfind::homotopy_system_with_gamma",DimensionPolicy::AnalysisHomotopy,(coefficients:Vec<f64>,exponents:Vec<Vec<u32>>,offsets:Vec<usize>,coordinate_reference:f64,value_reference:f64,gamma:num_complex::Complex64)->Result<super::AnalysisHomotopyOutput,OpError>,|coefficients:Vec<f64>,exponents:Vec<Vec<u32>>,offsets:Vec<usize>,coordinate_reference:f64,value_reference:f64,gamma:num_complex::Complex64|{if coordinate_reference<=0.0||value_reference<=0.0{return Err(OpError::Domain("homotopy references must be positive".into()));}let system=super::analysis_polynomial_system(coefficients,exponents,offsets)?;let paths=crate::kernels::rootfind::homotopy_system_with_gamma(&system,(gamma.re,gamma.im)).map_err(|e|OpError::Domain(format!("{e:?}")))?;super::analysis_homotopy_output(paths,coordinate_reference,value_reference)});
