#![allow(clippy::wildcard_imports)] // Actual model entry points preserve physical densities and native mutable wave output.
use super::*;
use crate::sim::models_gpe_superfluidity as gpe;
mod schrodinger;

const MAX_GRID_POINTS: usize = 65_536;
fn si_dimension(p: [i32; 4]) -> Result<DimExpr, OpError> {
    super::super::classical::si().dim(p)
}
fn read(
    v: &Input,
    p: [i32; 4],
) -> Result<f64, OpError> {
    api::same_dimension(&api::dimension(&v.0)?, &si_dimension(p)?)?;
    crate::rules::numerical_abi::numeric(&v.0, 0)
}
fn output_density(
    data: impl IntoIterator<Item = f64>,
    shape: Vec<usize>,
) -> Result<Value, OpError> {
    let dimension = DimExpr::base("Length").pow(-2);
    let data = data
        .into_iter()
        .map(|x| api::quantity(x, dimension.clone()))
        .collect::<Result<_, _>>()?;
    Ok(Value::Tensor(TensorValue {
        shape,
        data,
        elem: TypeExpr::Float.with_dim(dimension),
    }))
}
fn bounded(
    nodes: [usize; 2],
    steps: usize,
    minimum_extent: usize,
) -> Result<usize, OpError> {
    let count = nodes[0]
        .checked_mul(nodes[1])
        .filter(|n| *n <= MAX_GRID_POINTS)
        .ok_or_else(|| api::error("quantum grid allocation exceeds 65536 points"))?;
    if nodes.iter().any(|n| *n < minimum_extent)
        || steps == 0
        || steps > 100000
        || count.checked_mul(steps).is_none_or(|n| n > 50_000_000)
    {
        return Err(api::error("quantum grid-step budget exceeded"));
    }
    Ok(count)
}
fn directory(
    path: &str,
    budget: usize,
    required: usize,
) -> Result<(), OpError> {
    if path.is_empty() || path.len() > 4096 || budget < required {
        return Err(api::error(
            "native file scenario needsnonemptydirectory andexplicitgrid-update budget",
        ));
    }
    Ok(())
}
object_binding!(gpe_run_ground_state,"num_native_gpe_run_ground_state","crate::sim::models_gpe_superfluidity::run_gpe_ground_state_finder","Return the normalized two-dimensional Gross–Pitaevskii density in inverse square length on a grid of positive extents. Planck’s reduced constant, mass and lengths use SI. Imaginary step, contact energy times length squared and trap energy per length squared are scaled in the same reference frame. Zero steps or zero imaginary step return the normalized initial density.",false,(nodes:[usize;2],lengths:Input,mass:Input,coordinate_reference:Input,imaginary_dt:Input,steps:usize,coupling:Input,trap_strength:Input)->Dynamic,|_ctx:&CallCtx<'_>,nodes:[usize;2],lengths:Input,mass:Input,coordinate_reference:Input,imaginary_dt:Input,steps:usize,coupling:Input,trap_strength:Input|{let count=bounded(nodes,steps.max(1),1)?;api::same_dimension(&api::dimension(&lengths.0)?,&DimExpr::base("Length"))?;let length=<[f64;2] as Arg>::decode(&lengths.0,0)?;let(l,d)=api::reference(&coordinate_reference.0)?;api::same_dimension(&d,&DimExpr::base("Length"))?;let m=read(&mass,MASS)?;let dt=read(&imaginary_dt,TIME)?;let g=read(&coupling,[1,4,-2,0])?;let trap=read(&trap_strength,[1,0,-2,0])?;if m<=0.0||dt<0.0||trap<0.0||length.iter().any(|x|*x<=0.0){return Err(api::error("GPE mass/time/domain/trap controls invalid"));}let time=api::scale_product(m,&[l,l],&[HBAR])?;let energy=api::scale_product(HBAR,&[],&[time])?;let g=api::scale_product(g,&[],&[energy,l,l])?;let trap=api::scale_product(trap,&[l,l],&[energy])?;let params=gpe::GpeParameters{nx:nodes[0],ny:nodes[1],lx:length[0]/l,ly:length[1]/l,d_tau:api::scale_without_underflow(dt,&[],&[time])?,time_steps:steps,g,trap_strength:trap};let density=gpe::run_gpe_ground_state_finder(&params).map_err(api::error)?;let data=density.into_iter().map(|x|api::scale_product(x,&[],&[l,l])).collect::<Result<Vec<_>,_>>()?;if data.len()!=count{return Err(api::error("native density shape mismatch"));}Ok(Dynamic(output_density(data,vec![nodes[1],nodes[0]])?))});
object_binding!(gpe_vortex_scenario,"num_native_gpe_vortex_scenario","crate::sim::models_gpe_superfluidity::simulate_bose_einstein_vortex_scenario","IO: execute the original reduced-unit128x128,500-step GPE demonstration and write its NPY density to the requested directory. Explicit budget acknowledgement is required; no physical interpretation is inferred for hard-coded model constants.",true,(output_directory:String,grid_update_budget:usize)->bool,|_ctx:&CallCtx<'_>,output_directory:String,grid_update_budget:usize|{directory(&output_directory,grid_update_budget,128*128*500)?;gpe::simulate_bose_einstein_vortex_scenario(std::path::Path::new(&output_directory)).map_err(api::error)?;Ok(true)});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &gpe_run_ground_state::BINDING,
            &gpe_vortex_scenario::BINDING,
        ],
        "numerical_vector",
    )?;
    schrodinger::install(i)
}
