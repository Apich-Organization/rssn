#![allow(clippy::wildcard_imports)] // Fixed-domain native multigrid uses explicit coordinate/solution references.
use super::*;
use crate::rules::numerical_abi::codecs::Arg;
use crate::rules::numerical_abi::codecs::Output;
use crate::sim::physics_mtm as mtm;
fn reference(v: &Input) -> Result<(f64, DimExpr), OpError> {
    api::reference(&v.0)
}
fn result(
    v: Vec<f64>,
    shape: Vec<usize>,
    reference: f64,
    dimension: DimExpr,
    domain: Value,
) -> Result<Dynamic, OpError> {
    let data = v
        .into_iter()
        .map(|x| api::quantity(api::scale_product(x, &[reference], &[])?, dimension.clone()))
        .collect::<Result<Vec<_>, _>>()?;
    Ok(Dynamic(api::record(vec![
        (
            "solution",
            Value::Tensor(TensorValue {
                shape: vec![data.len()],
                data,
                elem: TypeExpr::Float.with_dim(dimension),
            }),
        ),
        ("shape", Output::encode(shape)?),
        ("domain_reference", domain),
    ])))
}
object_binding!(mtm_solve_poisson_1d_multigrid,"num_native_mtm_solve_poisson_1d_multigrid","crate::sim::physics_mtm::solve_poisson_1d_multigrid","Actual native zero-Dirichlet multigrid on its unit interval/cube, mapped to an explicit physical coordinate and solution reference. RHS has solution/coordinate^2 units.1D input contains interior values;2D/3D input includes boundary nodes. Output is flattened in native x-fastest order with explicit shape.",false,(n:usize,rhs:Input,cycles:usize,coordinate_reference:Input,solution_reference:Input)->Dynamic,|_ctx:&CallCtx<'_>,n:usize,rhs:Input,cycles:usize,coordinate_reference:Input,solution_reference:Input|{let(l,ld)=reference(&coordinate_reference)?;let(y,yd)=reference(&solution_reference)?;let expected=DimExpr::Mul(vec![yd.clone(),ld.pow(Exp::Rat(Rational::int(-2)))]);api::same_dimension(&expected,&api::dimension(&rhs.0)?)?;let values=<Vec<f64> as Arg>::decode(&rhs.0,0)?;if n==0||n>4095||!n.checked_add(1).is_some_and(usize::is_power_of_two)||values.len()!=n{return Err(api::error("1D multigrid needsn=2^k-1<=4095andninteriorvalues"));}let nodes=n+2;if cycles==0||cycles>1000||values.len().checked_mul(cycles).is_none_or(|work|work>50_000_000){return Err(api::error("multigrid cyclework exceedsbudget"));}let normalized=values.into_iter().map(|v|api::scale_product(v,&[l,l],&[y])).collect::<Result<Vec<_>,_>>()?;let solution=mtm::solve_poisson_1d_multigrid(n,&normalized,cycles).map_err(api::error)?;result(solution,vec![nodes;1],y,yd,coordinate_reference.0)});
object_binding!(mtm_solve_poisson_2d_multigrid,"num_native_mtm_solve_poisson_2d_multigrid","crate::sim::physics_mtm::solve_poisson_2d_multigrid","Actual native zero-Dirichlet multigrid on its unit interval/cube, mapped to an explicit physical coordinate and solution reference. RHS has solution/coordinate^2 units.1D input contains interior values;2D/3D input includes boundary nodes. Output is flattened in native x-fastest order with explicit shape.",false,(n:usize,rhs:Input,cycles:usize,coordinate_reference:Input,solution_reference:Input)->Dynamic,|_ctx:&CallCtx<'_>,n:usize,rhs:Input,cycles:usize,coordinate_reference:Input,solution_reference:Input|{let(l,ld)=reference(&coordinate_reference)?;let(y,yd)=reference(&solution_reference)?;let expected=DimExpr::Mul(vec![yd.clone(),ld.pow(Exp::Rat(Rational::int(-2)))]);api::same_dimension(&expected,&api::dimension(&rhs.0)?)?;let values=<Vec<f64> as Arg>::decode(&rhs.0,0)?;if !(3..=129).contains(&n)||!(n-1).is_power_of_two()||n.checked_pow(2).is_none_or(|count|count!=values.len()){return Err(api::error("unitcube multigrid needsn=2^k+1andn^Dvalues withinallocation budget"));}let nodes=n;if cycles==0||cycles>1000||values.len().checked_mul(cycles).is_none_or(|work|work>50_000_000){return Err(api::error("multigrid cyclework exceedsbudget"));}let normalized=values.into_iter().map(|v|api::scale_product(v,&[l,l],&[y])).collect::<Result<Vec<_>,_>>()?;let solution=mtm::solve_poisson_2d_multigrid(n,&normalized,cycles).map_err(api::error)?;result(solution,vec![nodes;2],y,yd,coordinate_reference.0)});
object_binding!(mtm_solve_poisson_3d_multigrid,"num_native_mtm_solve_poisson_3d_multigrid","crate::sim::physics_mtm::solve_poisson_3d_multigrid","Actual native zero-Dirichlet multigrid on its unit interval/cube, mapped to an explicit physical coordinate and solution reference. RHS has solution/coordinate^2 units.1D input contains interior values;2D/3D input includes boundary nodes. Output is flattened in native x-fastest order with explicit shape.",false,(n:usize,rhs:Input,cycles:usize,coordinate_reference:Input,solution_reference:Input)->Dynamic,|_ctx:&CallCtx<'_>,n:usize,rhs:Input,cycles:usize,coordinate_reference:Input,solution_reference:Input|{let(l,ld)=reference(&coordinate_reference)?;let(y,yd)=reference(&solution_reference)?;let expected=DimExpr::Mul(vec![yd.clone(),ld.pow(Exp::Rat(Rational::int(-2)))]);api::same_dimension(&expected,&api::dimension(&rhs.0)?)?;let values=<Vec<f64> as Arg>::decode(&rhs.0,0)?;if !(3..=65).contains(&n)||!(n-1).is_power_of_two()||n.checked_pow(3).is_none_or(|count|count!=values.len()){return Err(api::error("unitcube multigrid needsn=2^k+1andn^Dvalues withinallocation budget"));}let nodes=n;if cycles==0||cycles>1000||values.len().checked_mul(cycles).is_none_or(|work|work>50_000_000){return Err(api::error("multigrid cyclework exceedsbudget"));}let normalized=values.into_iter().map(|v|api::scale_product(v,&[l,l],&[y])).collect::<Result<Vec<_>,_>>()?;let solution=mtm::solve_poisson_3d_multigrid(n,&normalized,cycles).map_err(api::error)?;result(solution,vec![nodes;3],y,yd,coordinate_reference.0)});
object_binding!(mtm_simulate_1d_poisson_multigrid_scenario,"num_native_mtm_simulate_1d_poisson_multigrid_scenario","crate::sim::physics_mtm::simulate_1d_poisson_multigrid_scenario","Execute the original native analytic multigrid demonstration in normalized coordinates; explicit positive domain/solution references give a physical realization without reinterpreting its hard-coded coefficients.",false,(coordinate_reference:Input,solution_reference:Input)->Dynamic,|_ctx:&CallCtx<'_>,coordinate_reference:Input,solution_reference:Input|{let(_,_) = reference(&coordinate_reference)?;let(y,d)=reference(&solution_reference)?;let solution=mtm::simulate_1d_poisson_multigrid_scenario().map_err(api::error)?;result(solution,vec![129;1],y,d,coordinate_reference.0)});
object_binding!(mtm_simulate_2d_poisson_multigrid_scenario,"num_native_mtm_simulate_2d_poisson_multigrid_scenario","crate::sim::physics_mtm::simulate_2d_poisson_multigrid_scenario","Execute the original native analytic multigrid demonstration in normalized coordinates; explicit positive domain/solution references give a physical realization without reinterpreting its hard-coded coefficients.",false,(coordinate_reference:Input,solution_reference:Input)->Dynamic,|_ctx:&CallCtx<'_>,coordinate_reference:Input,solution_reference:Input|{let(_,_) = reference(&coordinate_reference)?;let(y,d)=reference(&solution_reference)?;let solution=mtm::simulate_2d_poisson_multigrid_scenario().map_err(api::error)?;result(solution,vec![33;2],y,d,coordinate_reference.0)});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &mtm_solve_poisson_1d_multigrid::BINDING,
            &mtm_solve_poisson_2d_multigrid::BINDING,
            &mtm_solve_poisson_3d_multigrid::BINDING,
            &mtm_simulate_1d_poisson_multigrid_scenario::BINDING,
            &mtm_simulate_2d_poisson_multigrid_scenario::BINDING,
        ],
        "numerical_vector",
    )
}
