#![allow(clippy::wildcard_imports)] // Actual model entry points preserve physical densities and native mutable wave output.
use super::*;
use crate::sim::models_gpe_superfluidity as gpe;
use crate::sim::models_schrodinger_quantum as quantum;
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
    data: Vec<f64>,
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
) -> Result<usize, OpError> {
    let count = nodes[0]
        .checked_mul(nodes[1])
        .filter(|n| *n <= 65536)
        .ok_or_else(|| api::error("quantum grid allocationbound65536"))?;
    if nodes.iter().any(|n| *n < 2)
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
object_binding!(gpe_run_ground_state,"num_native_gpe_run_ground_state","crate::sim::models_gpe_superfluidity::run_gpe_ground_state_finder","Actual native2D GPE ground-state finder with hbar=physicalSI, explicit coordinate reference and SI mass. Imaginary time, contactEnergy*Length² and trapEnergy/Length² are normalized consistently; returned probability density has Length^-2 and unit integral.",false,(nodes:[usize;2],lengths:Input,mass:Input,coordinate_reference:Input,imaginary_dt:Input,steps:usize,coupling:Input,trap_strength:Input)->Dynamic,|_ctx:&CallCtx<'_>,nodes:[usize;2],lengths:Input,mass:Input,coordinate_reference:Input,imaginary_dt:Input,steps:usize,coupling:Input,trap_strength:Input|{let count=bounded(nodes,steps)?;api::same_dimension(&api::dimension(&lengths.0)?,&DimExpr::base("Length"))?;let length=<[f64;2] as Arg>::decode(&lengths.0,0)?;let(l,d)=api::reference(&coordinate_reference.0)?;api::same_dimension(&d,&DimExpr::base("Length"))?;let m=read(&mass,MASS)?;let dt=read(&imaginary_dt,TIME)?;let g=read(&coupling,[1,4,-2,0])?;let trap=read(&trap_strength,[1,0,-2,0])?;if m<=0.0||dt<=0.0||trap<0.0||length.iter().any(|x|*x<=0.0){return Err(api::error("GPE mass/time/domain/trap controls invalid"));}let time=api::scale_product(m,&[l,l],&[HBAR])?;let energy=api::scale_product(HBAR,&[],&[time])?;let g=api::scale_product(g,&[],&[energy,l,l])?;let trap=api::scale_product(trap,&[l,l],&[energy])?;let params=gpe::GpeParameters{nx:nodes[0],ny:nodes[1],lx:length[0]/l,ly:length[1]/l,d_tau:dt/time,time_steps:steps,g,trap_strength:trap};let density=gpe::run_gpe_ground_state_finder(&params).map_err(api::error)?;let data=density.into_iter().map(|x|api::scale_product(x,&[],&[l,l])).collect::<Result<Vec<_>,_>>()?;if data.len()!=count{return Err(api::error("native density shape mismatch"));}Ok(Dynamic(output_density(data,vec![nodes[1],nodes[0]])?))});
object_binding!(gpe_vortex_scenario,"num_native_gpe_vortex_scenario","crate::sim::models_gpe_superfluidity::simulate_bose_einstein_vortex_scenario","IO: execute the original reduced-unit128x128,500-step GPE demonstration and write its NPY density to the requested directory. Explicit budget acknowledgement is required; no physical interpretation is inferred for hard-coded model constants.",true,(output_directory:String,grid_update_budget:usize)->bool,|_ctx:&CallCtx<'_>,output_directory:String,grid_update_budget:usize|{directory(&output_directory,grid_update_budget,128*128*500)?;gpe::simulate_bose_einstein_vortex_scenario(std::path::Path::new(&output_directory)).map_err(api::error)?;Ok(true)});
object_binding!(quantum_run_schrodinger,"num_native_quantum_run_schrodinger","crate::sim::models_schrodinger_quantum::run_schrodinger_simulation","Actual native2D Schrodinger entry point with SI hbar, mass, lengths, energy potential and Length^-1 wave amplitude. Native mutable final state is returned with Length^-2 density frames and actual snapshot times; all input/output sizes and grid-step work are bounded.",false,(nodes:[usize;2],lengths:Input,dt:Input,steps:usize,hbar:Input,mass:Input,potential:Input,initial_psi:Input)->Dynamic,|_ctx:&CallCtx<'_>,nodes:[usize;2],lengths:Input,dt:Input,steps:usize,hbar:Input,mass:Input,potential:Input,initial_psi:Input|{let count=bounded(nodes,steps)?;if count.checked_mul(steps.div_ceil(10)).is_none_or(|n|n>262144){return Err(api::error("quantum densityhistory output bound exceeded"));}api::same_dimension(&api::dimension(&lengths.0)?,&DimExpr::base("Length"))?;let length=<[f64;2] as Arg>::decode(&lengths.0,0)?;let dt=read(&dt,TIME)?;let hbar=read(&hbar,[1,2,-1,0])?;let mass=read(&mass,MASS)?;if length.iter().any(|l|*l<=0.0)||dt<=0.0||hbar<=0.0||mass<=0.0{return Err(api::error("Schrodinger domain/time/hbar/mass invalid"));}api::same_dimension(&api::dimension(&potential.0)?,&si_dimension(ENERGY)?)?;let potential=<Vec<f64> as Arg>::decode(&potential.0,0)?;let psi_dimension=DimExpr::base("Length").pow(-1);api::same_dimension(&api::dimension(&initial_psi.0)?,&psi_dimension)?;let Value::Tensor(t)=&initial_psi.0 else{return Err(api::error("wavefunction tensor expected"));};
if t.shape!=[count]||potential.len()!=count{return Err(api::error("Schrodinger potential/psi sizes differ fromgrid"));}let mut psi=t.data.iter().map(complex).collect::<Result<Vec<_>,_>>()?;let params=quantum::SchrodingerParameters{nx:nodes[0],ny:nodes[1],lx:length[0],ly:length[1],dt,time_steps:steps,hbar,mass,potential};let frames=quantum::run_schrodinger_simulation(&params,&mut psi).map_err(api::error)?;let n=frames.len();let density=output_density(frames.into_iter().flat_map(IntoIterator::into_iter).collect(),vec![n,nodes[1],nodes[0]])?;let psi=Value::Tensor(TensorValue{shape:vec![count],data:psi.into_iter().map(Output::encode).collect::<Result<_,_>>()?,elem:TypeExpr::Complex.with_dim(psi_dimension)});let times=(0..steps).step_by(10).map(|i|api::quantity((i+1) as f64*dt,DimExpr::base("Time"))).collect::<Result<Vec<_>,_>>()?;Ok(Dynamic(api::record(vec![("density_frames",density),("final_psi",psi),("snapshot_times",Value::Tensor(TensorValue{shape:vec![times.len()],data:times,elem:TypeExpr::Float.with_dim(DimExpr::base("Time"))}))])))});
object_binding!(quantum_double_slit_scenario,"num_native_quantum_double_slit_scenario","crate::sim::models_schrodinger_quantum::simulate_double_slit_scenario","IO: execute the original reduced-unit256x256,300-step double-slit demonstration and write its NPY density. Explicit grid-update acknowledgement is required.",true,(output_directory:String,grid_update_budget:usize)->bool,|_ctx:&CallCtx<'_>,output_directory:String,grid_update_budget:usize|{directory(&output_directory,grid_update_budget,256*256*300)?;quantum::simulate_double_slit_scenario(std::path::Path::new(&output_directory)).map_err(api::error)?;Ok(true)});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &gpe_run_ground_state::BINDING,
            &gpe_vortex_scenario::BINDING,
            &quantum_run_schrodinger::BINDING,
            &quantum_double_slit_scenario::BINDING,
        ],
        "numerical_vector",
    )
}
