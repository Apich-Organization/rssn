//! Native high-level PDE entry points, preserving each public method identity.
use super::CallCtx;
use super::CompiledSnapshot;
use super::DimExpr;
use super::Dynamic;
use super::Exp;
use super::Frame;
use super::Handle;
use super::Input;
use super::Installer;
use super::NativeType;
use super::NewObject;
use super::OpError;
use super::Output;
use super::Rational;
use super::RuleError;
use super::Value;
use super::api;
use super::canonical;
use super::err;
use super::fdm::Bcs;
use super::fdm::Condition;
use super::fdm::NativeBcs;
use super::object_binding;
use super::payload;
use crate::sim::pde_api as pde;
#[derive(Clone)]
struct Fem {
    data: pde::FemSolution,
    frame: Frame,
}
impl NativeType for Fem {
    const KIND: &'static str = "rssn/native/pde-fem-solution/v1";
    const MAX_BYTES: usize = 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        self.frame.validate()?;
        if self.data.values.len() > 1024
            || self.data.values.len() != self.data.dof_coords.len()
            || self.data.values.iter().any(|v| !v.is_finite())
            || self
                .data
                .dof_coords
                .iter()
                .flatten()
                .any(|v| !v.is_finite())
        {
            return Err(err("PDE FEM finite shape/storage bounds exceeded"));
        }
        Ok(())
    }
}
fn degree(p: usize) -> Result<pde::FemDegree, OpError> {
    match p {
        | 1 => Ok(pde::FemDegree::P1),
        | 2 => Ok(pde::FemDegree::P2),
        | _ => Err(err("FEM degree must be one/two")),
    }
}
fn size(
    extents: &[usize],
    p: usize,
) -> Result<(), OpError> {
    if !(1..=3).contains(&extents.len())
        || extents.iter().any(|n| *n == 0 || *n > 512)
        || extents
            .iter()
            .try_fold(1usize, |n, m| {
                n.checked_mul(m.checked_mul(p)?.checked_add(1)?)
            })
            .is_none_or(|n| n > 1024)
    {
        return Err(err("PDE FEM mesh/DOF allocation budget exceeded"));
    }
    Ok(())
}
fn eval(
    c: &CompiledSnapshot,
    f: &Frame,
    x: &[f64],
    out: f64,
) -> f64 {
    c.sample(&x.iter().map(|x| x * f.length).collect::<Vec<_>>()) / out
}
fn source(
    ctx: &CallCtx<'_>,
    f: &Frame,
    program: &Input,
    budget: usize,
    dimension: usize,
    time_power: i64,
) -> Result<CompiledSnapshot, OpError> {
    if budget == 0 || budget > 8_388_608 {
        return Err(err("PDE callback budget invalid"));
    }
    super::structured_ops::checked_snapshot(
        ctx,
        program,
        budget,
        &vec![f.dims[0].clone(); dimension],
        &super::structured_ops::dims_forcing(f, time_power)?,
    )
}
fn coordinates(
    f: &Frame,
    x: Vec<f64>,
    shape: Vec<usize>,
) -> Result<Value, OpError> {
    Ok(api::record(vec![
        ("shape", shape.encode()?),
        ("normalized", x.encode()?),
        ("reference", api::quantity(f.length, f.dims[0].clone())?),
    ]))
}
fn paired(
    f: &Frame,
    x: Vec<f64>,
    u: Vec<f64>,
) -> Result<Dynamic, OpError> {
    let n = x.len();
    Ok(Dynamic(api::record(vec![
        ("coordinates", coordinates(f, x, vec![n])?),
        ("field", payload(f, u, vec![n])?.0),
    ])))
}
fn work(
    n: usize,
    steps: usize,
) -> Result<(), OpError> {
    if n == 0
        || n > 4096
        || steps == 0
        || steps > 4096
        || n.checked_mul(steps).is_none_or(|n| n > 4_194_304)
    {
        return Err(err("PDE evolution step/work budget exceeded"));
    }
    Ok(())
}
fn scheme(s: &str) -> Result<pde::FvmScheme, OpError> {
    match s {
        | "Upwind" => Ok(pde::FvmScheme::Upwind),
        | "LaxFriedrichs" => Ok(pde::FvmScheme::LaxFriedrichs),
        | "MusclMinmod" => Ok(pde::FvmScheme::MusclMinmod),
        | "MusclVanLeer" => Ok(pde::FvmScheme::MusclVanLeer),
        | _ => Err(err("unknown FVM scheme")),
    }
}
fn spectral(
    extents: &[usize],
    lengths: &[f64],
) -> Result<usize, OpError> {
    if !(1..=3).contains(&extents.len())
        || extents.len() != lengths.len()
        || extents.iter().any(|n| *n < 4 || !n.is_power_of_two())
        || lengths.iter().any(|v| *v <= 0.0)
        || extents
            .iter()
            .try_fold(1usize, |n, m| n.checked_mul(*m))
            .is_none_or(|n| n > 4096)
    {
        return Err(err(
            "PDE spectral geometry/power-of-two/allocation bounds invalid",
        ));
    }
    Ok(extents.iter().product())
}
object_binding!(pde_fem1,"num_pde_native_fem_1d","crate::sim::pde_api::solve_fem_1d","Actual public high-level PDE entry point, with explicit Length/Time/field references and original boundary, discretization and numerical method semantics. Arrays and controls are normalized in that frame; callbacks take coherent physical arguments. Native storage, work and callback failures are checked.",false,(bcs:Handle<Bcs>,lower:Vec<f64>,upper:Vec<f64>,cells:Vec<usize>,polynomial_degree:usize,diffusion:f64,reaction:f64,source_program:Input,budget:usize)->NewObject<Fem>,|ctx,bcs,lower,upper,cells,polynomial_degree,diffusion,reaction,source_program,budget|{bcs.with(ctx,|b|{let b_frame_reference_0_m1_1=b.frame.reference([0, -1, 1])?;if lower.len()!=1||upper.len()!=1||cells.len()!=1||lower.iter().zip(&upper).any(|(a,b)|b<=a)||diffusion<=0.0{return Err(err("FEM domain/coefficient shape invalid"));}let degree=degree(polynomial_degree)?;size(&cells,polynomial_degree)?;let c=source(ctx,&b.frame,&source_program,budget,1,1)?;let NativeBcs::D1(bc)=&b.data else{return Err(err("FEM boundary dimension mismatch"));};let out=pde::solve_fem_1d((lower[0],upper[0]),cells[0],degree,diffusion,reaction,|x|eval(&c,&b.frame,&[x],b_frame_reference_0_m1_1),bc).map_err(err);c.finish()?;Ok(NewObject::new(Fem{data:out?,frame:b.frame.clone()},b.frame.dims.to_vec()))})});
object_binding!(pde_fem2,"num_pde_native_fem_2d","crate::sim::pde_api::solve_fem_2d","Actual public high-level PDE entry point, with explicit Length/Time/field references and original boundary, discretization and numerical method semantics. Arrays and controls are normalized in that frame; callbacks take coherent physical arguments. Native storage, work and callback failures are checked.",false,(bcs:Handle<Bcs>,lower:Vec<f64>,upper:Vec<f64>,cells:Vec<usize>,polynomial_degree:usize,diffusion:f64,reaction:f64,source_program:Input,budget:usize)->NewObject<Fem>,|ctx,bcs,lower,upper,cells,polynomial_degree,diffusion,reaction,source_program,budget|{bcs.with(ctx,|b|{let b_frame_reference_0_m1_1=b.frame.reference([0, -1, 1])?;if lower.len()!=2||upper.len()!=2||cells.len()!=2||lower.iter().zip(&upper).any(|(a,b)|b<=a)||diffusion<=0.0{return Err(err("FEM domain/coefficient shape invalid"));}let degree=degree(polynomial_degree)?;size(&cells,polynomial_degree)?;let c=source(ctx,&b.frame,&source_program,budget,2,1)?;let NativeBcs::D2(bc)=&b.data else{return Err(err("FEM boundary dimension mismatch"));};let out=pde::solve_fem_2d((lower[0],upper[0],lower[1],upper[1]),cells[0],cells[1],degree,diffusion,reaction,|x,y|eval(&c,&b.frame,&[x,y],b_frame_reference_0_m1_1),bc).map_err(err);c.finish()?;Ok(NewObject::new(Fem{data:out?,frame:b.frame.clone()},b.frame.dims.to_vec()))})});
object_binding!(pde_fem3,"num_pde_native_fem_3d","crate::sim::pde_api::solve_fem_3d","Actual public high-level PDE entry point, with explicit Length/Time/field references and original boundary, discretization and numerical method semantics. Arrays and controls are normalized in that frame; callbacks take coherent physical arguments. Native storage, work and callback failures are checked.",false,(bcs:Handle<Bcs>,lower:Vec<f64>,upper:Vec<f64>,cells:Vec<usize>,polynomial_degree:usize,diffusion:f64,reaction:f64,source_program:Input,budget:usize)->NewObject<Fem>,|ctx,bcs,lower,upper,cells,polynomial_degree,diffusion,reaction,source_program,budget|{bcs.with(ctx,|b|{let b_frame_reference_0_m1_1=b.frame.reference([0, -1, 1])?;if lower.len()!=3||upper.len()!=3||cells.len()!=3||lower.iter().zip(&upper).any(|(a,b)|b<=a)||diffusion<=0.0{return Err(err("FEM domain/coefficient shape invalid"));}let degree=degree(polynomial_degree)?;size(&cells,polynomial_degree)?;let c=source(ctx,&b.frame,&source_program,budget,3,1)?;let NativeBcs::D3(bc)=&b.data else{return Err(err("FEM boundary dimension mismatch"));};let out=pde::solve_fem_3d(lower.try_into().map_err(|_|err("3D lower"))?,upper.try_into().map_err(|_|err("3D upper"))?,cells.try_into().map_err(|_|err("3D cells"))?,degree,diffusion,reaction,|x,y,z|eval(&c,&b.frame,&[x,y,z],b_frame_reference_0_m1_1),bc).map_err(err);c.finish()?;Ok(NewObject::new(Fem{data:out?,frame:b.frame.clone()},b.frame.dims.to_vec()))})});
object_binding!(pde_fem_adaptive,"num_pde_native_fem_adaptive_2d","crate::sim::pde_api::solve_fem_adaptive_2d","Actual public high-level PDE entry point, with explicit Length/Time/field references and original boundary, discretization and numerical method semantics. Arrays and controls are normalized in that frame; callbacks take coherent physical arguments. Native storage, work and callback failures are checked.",false,(frame:Handle<Frame>,lower:Vec<f64>,upper:Vec<f64>,cells:Vec<usize>,diffusion:f64,reaction:f64,source_program:Input,dirichlet:Input,refinements:usize,theta:f64,budget:usize)->NewObject<Fem>,|ctx,frame,lower,upper,cells,diffusion,reaction,source_program,dirichlet,refinements,theta,budget|{frame.with(ctx,|f|{let f_reference_0_m1_1=f.reference([0, -1, 1])?;if lower.len()!=2||upper.len()!=2||cells.len()!=2||lower.iter().zip(&upper).any(|(a,b)|b<=a)||diffusion<=0.0||refinements>4||theta<=0.0||theta>1.0{return Err(err("adaptive FEM geometry/control domain invalid"));}size(&cells,1)?;if cells.iter().product::<usize>().checked_mul(4usize.pow(refinements as u32)).is_none_or(|n|n>128){return Err(err("adaptive triangle worst-case refinement budget exceeded"));}let c=source(ctx,f,&source_program,budget,2,1)?;let g=source(ctx,f,&dirichlet,budget,2,0)?;let out=pde::solve_fem_adaptive_2d((lower[0],upper[0],lower[1],upper[1]),cells[0],cells[1],diffusion,reaction,|x,y|eval(&c,f,&[x,y],f_reference_0_m1_1),|x,y|eval(&g,f,&[x,y],f.value),refinements,theta).map_err(err);c.finish()?;g.finish()?;Ok(NewObject::new(Fem{data:out?,frame:f.clone()},f.dims.to_vec()))})});
object_binding!(pde_fem_error,"num_pde_native_fem_max_error","crate::sim::pde_api::FemSolution::max_nodal_error","Actual public high-level PDE entry point, with explicit Length/Time/field references and original boundary, discretization and numerical method semantics. Arrays and controls are normalized in that frame; callbacks take coherent physical arguments. Native storage, work and callback failures are checked.",false,(solution:Handle<Fem>,exact:Input,budget:usize)->Dynamic,|ctx,solution,exact,budget|{solution.with(ctx,|s|{let c=source(ctx,&s.frame,&exact,budget,3,0)?;let out=s.data.max_nodal_error(|x|eval(&c,&s.frame,&x,s.frame.value));c.finish()?;Ok(Dynamic(api::quantity(out*s.frame.value,s.frame.dims[2].clone())?))})});
object_binding!(pde_fem_data,"num_pde_native_fem_data","","Actual public high-level PDE entry point, with explicit Length/Time/field references and original boundary, discretization and numerical method semantics. Arrays and controls are normalized in that frame; callbacks take coherent physical arguments. Native storage, work and callback failures are checked.",true,(solution:Handle<Fem>)->Dynamic,|ctx,solution|{solution.with(ctx,|s|Ok(Dynamic(api::record(vec![("coordinates",coordinates(&s.frame,s.data.dof_coords.iter().flatten().copied().collect(),vec![s.data.values.len(),3])?),("field",payload(&s.frame,s.data.values.clone(),vec![s.data.values.len()])?.0)]))))});
object_binding!(pde_dg_advection,"num_pde_native_dg_advection_1d","crate::sim::pde_api::solve_dg_advection_1d","Actual public high-level PDE entry point, with explicit Length/Time/field references and original boundary, discretization and numerical method semantics. Arrays and controls are normalized in that frame; callbacks take coherent physical arguments. Native storage, work and callback failures are checked.",false,(frame:Handle<Frame>,lower:f64,upper:f64,elements:usize,polynomial_degree:usize,velocity:f64,t_end:f64,cfl:f64,initial:Input,budget:usize)->Dynamic,|ctx,frame,lower,upper,elements,polynomial_degree,velocity,t_end,cfl,initial,budget|{frame.with(ctx,|f|{if elements==0||elements>256||polynomial_degree==0||polynomial_degree>8||elements*(polynomial_degree+1)>4096||upper<=lower||t_end<0.0||cfl<=0.0||cfl>1.0||velocity.abs()*t_end*(2*polynomial_degree+1)as f64*elements as f64/((upper-lower)*cfl)>4096.0{return Err(err("DG geometry/CFL/constant-speed step budget invalid"));}let c=source(ctx,f,&initial,budget,1,0)?;let out=pde::solve_dg_advection_1d((lower,upper),elements,polynomial_degree,velocity,t_end,cfl,|x|eval(&c,f,&[x],f.value)).map_err(err);c.finish()?;let(x,u)=out?;paired(f,x,u)})});
object_binding!(pde_solve_fvm_conservation_1d,"num_pde_native_solve_fvm_conservation_1d","crate::sim::pde_api::solve_fvm_conservation_1d","Actual public high-level PDE entry point, with explicit Length/Time/field references and original boundary, discretization and numerical method semantics. Arrays and controls are normalized in that frame; callbacks take coherent physical arguments. Native storage, work and callback failures are checked.",false,(frame:Handle<Frame>,cells:usize,length:f64,dt:f64,steps:usize,periodic:bool,method:String,initial:Input,budget:usize,flux:Input,wave_speed:Input)->Dynamic,|ctx,frame,cells,length,dt,steps,periodic,method,initial,budget,flux,wave_speed|{frame.with(ctx,|f|{let f_reference_1_m1_0=f.reference([1, -1, 0])?;let f_reference_1_m1_1=f.reference([1, -1, 1])?;work(cells,steps)?;if cells<3||length<=0.0||dt<=0.0{return Err(err("FVM grid/time invalid"));}let c=source(ctx,f,&initial,budget,1,0)?;let scheme=scheme(&method)?;let velocitydim=canonical(DimExpr::Mul(vec![f.dims[0].clone(),DimExpr::Inverse(Box::new(f.dims[1].clone()))]))?;let fluxdim=canonical(DimExpr::Mul(vec![f.dims[2].clone(),velocitydim.clone()]))?;let fl=crate::rules::numerical_abi::native_grids::structured_ops::checked_snapshot(ctx,&flux,budget,&[f.dims[2].clone()],&fluxdim)?;let speed=crate::rules::numerical_abi::native_grids::structured_ops::checked_snapshot(ctx,&wave_speed,budget,&[f.dims[2].clone()],&velocitydim)?;let out=pde::solve_fvm_conservation_1d(cells,length,dt,steps,periodic,scheme,|x|eval(&c,f,&[x],f.value),|u|fl.sample(&[u*f.value])/(f_reference_1_m1_1),|u|speed.sample(&[u*f.value])/(f_reference_1_m1_0)).map_err(err);c.finish()?;fl.finish()?;speed.finish()?;let(x,u)=out?;paired(f,x,u)})});
object_binding!(pde_solve_fvm_burgers_1d,"num_pde_native_solve_fvm_burgers_1d","crate::sim::pde_api::solve_fvm_burgers_1d","Actual public high-level PDE entry point, with explicit Length/Time/field references and original boundary, discretization and numerical method semantics. Arrays and controls are normalized in that frame; callbacks take coherent physical arguments. Native storage, work and callback failures are checked.",false,(frame:Handle<Frame>,cells:usize,length:f64,dt:f64,steps:usize,periodic:bool,method:String,initial:Input,budget:usize)->Dynamic,|ctx,frame,cells,length,dt,steps,periodic,method,initial,budget|{frame.with(ctx,|f|{let f_reference_1_m1_0=f.reference([1, -1, 0])?;work(cells,steps)?;if cells<3||length<=0.0||dt<=0.0{return Err(err("FVM grid/time invalid"));}let c=source(ctx,f,&initial,budget,1,0)?;let scheme=scheme(&method)?;let vd=canonical(DimExpr::Mul(vec![f.dims[0].clone(),DimExpr::Inverse(Box::new(f.dims[1].clone()))]))?;api::same_dimension(&f.dims[2],&vd)?;if f.value.to_bits()!=(f_reference_1_m1_0).to_bits(){return Err(err("Burgers field reference must equal L/T"));}let out=pde::solve_fvm_burgers_1d(cells,length,dt,steps,periodic,scheme,|x|eval(&c,f,&[x],f.value)).map_err(err);c.finish()?;let(x,u)=out?;paired(f,x,u)})});
object_binding!(pde_shallow,"num_pde_native_fvm_shallow_water_1d","crate::sim::pde_api::solve_fvm_shallow_water_1d","Actual public high-level PDE entry point, with explicit Length/Time/field references and original boundary, discretization and numerical method semantics. Arrays and controls are normalized in that frame; callbacks take coherent physical arguments. Native storage, work and callback failures are checked.",true,(frame:Handle<Frame>,depth:Vec<f64>,momentum:Vec<f64>,spacing:f64,dt:f64,steps:usize,gravity:f64)->Dynamic,|ctx,frame,depth,momentum,spacing,dt,steps,gravity|{frame.with(ctx,|f|{let f_reference_2_m1_0=f.reference([2, -1, 0])?;work(depth.len(),steps)?;api::same_dimension(&f.dims[2],&f.dims[0])?;if f.value.to_bits()!=f.length.to_bits()||depth.len()!=momentum.len()||spacing<=0.0||dt<=0.0||gravity<=0.0{return Err(err("shallow-water frame/state/grid controls invalid"));}let(h,hu)=pde::solve_fvm_shallow_water_1d(&depth,&momentum,spacing,dt,steps,gravity).map_err(err)?;let momentumdim=canonical(DimExpr::Mul(vec![DimExpr::Pow(Box::new(f.dims[0].clone()),Exp::Rat(Rational::new(2,1))),DimExpr::Inverse(Box::new(f.dims[1].clone()))]))?;Ok(Dynamic(api::record(vec![("depth",payload(f,h,vec![depth.len()])?.0),("momentum_normalized",hu.encode()?),("momentum_reference",api::quantity(f_reference_2_m1_0,momentumdim)?)])))})});
object_binding!(pde_spectral_helm,"num_pde_native_spectral_helmholtz","crate::sim::pde_api::solve_spectral_helmholtz","Actual public high-level PDE entry point, with explicit Length/Time/field references and original boundary, discretization and numerical method semantics. Arrays and controls are normalized in that frame; callbacks take coherent physical arguments. Native storage, work and callback failures are checked.",true,(frame:Handle<Frame>,extents:Vec<usize>,lengths:Vec<f64>,sigma:f64,source:Vec<f64>)->Dynamic,|ctx,frame,extents,lengths,sigma,source|{frame.with(ctx,|f|{let n=spectral(&extents,&lengths)?;if source.len()!=n{return Err(err("spectral source shape mismatch"));}let u=pde::solve_spectral_helmholtz(&extents,&lengths,sigma,&source).map_err(err)?;payload(f,u,extents)})});
object_binding!(pde_spectral_etd,"num_pde_native_spectral_etdrk4","crate::sim::pde_api::solve_spectral_etdrk4","Actual public high-level PDE entry point, with explicit Length/Time/field references and original boundary, discretization and numerical method semantics. Arrays and controls are normalized in that frame; callbacks take coherent physical arguments. Native storage, work and callback failures are checked.",false,(frame:Handle<Frame>,extents:Vec<usize>,lengths:Vec<f64>,diffusion:f64,reaction:f64,t_end:f64,steps:usize,initial:Vec<f64>,nonlinear:Input,budget:usize)->Dynamic,|ctx,frame,extents,lengths,diffusion,reaction,t_end,steps,initial,nonlinear,budget|{frame.with(ctx,|f|{let f_reference_0_m1_1=f.reference([0, -1, 1])?;let n=spectral(&extents,&lengths)?;work(n,steps)?;if n>62||initial.len()!=n||diffusion<0.0||t_end<=0.0||budget==0||budget>8_388_608{return Err(err("spectral coupled callback domain/work budget invalid"));}let Value::Tensor(t)=&nonlinear.0 else{return Err(err("nonlinearity needs one scalar compiled program per physical grid point"));};
if t.shape!=[n]||t.data.len()!=n{return Err(err("nonlinearity callback count mismatch"));}let mut dims=vec![f.dims[1].clone()];dims.extend(vec![f.dims[2].clone();n]);let output=crate::rules::numerical_abi::native_grids::structured_ops::dims_forcing(f,1)?;let c=t.data.iter().map(|p|crate::rules::numerical_abi::native_grids::structured_ops::checked_snapshot(ctx,&Input(p.clone()),budget,&dims,&output)).collect::<Result<Vec<_>,_>>()?;let out=pde::solve_spectral_etdrk4(&extents,&lengths,diffusion,reaction,t_end,steps,&initial,|time,u|{let mut x=vec![time*f.time];x.extend(u.iter().map(|v|v*f.value));c.iter().map(|c|c.sample(&x)/(f_reference_0_m1_1)).collect()}).map_err(err);for c in &c{c.finish()?;}payload(f,out?,extents)})});
object_binding!(pde_cheb_helm,"num_pde_native_chebyshev_helmholtz","crate::sim::pde_api::solve_chebyshev_helmholtz","Actual public high-level PDE entry point, with explicit Length/Time/field references and original boundary, discretization and numerical method semantics. Arrays and controls are normalized in that frame; callbacks take coherent physical arguments. Native storage, work and callback failures are checked.",false,(boundary:Handle<Condition>,degrees:Vec<usize>,lower:Vec<f64>,upper:Vec<f64>,sigma:f64,source_program:Input,budget:usize)->Dynamic,|ctx,boundary,degrees,lower,upper,sigma,source_program,budget|{boundary.with(ctx,|b|{if !(1..=3).contains(&degrees.len())||degrees.iter().any(|n|*n<2||*n>128)||degrees.iter().try_fold(1usize,|n,d|n.checked_mul(d+1)).is_none_or(|n|n>512)||lower.len()!=degrees.len()||upper.len()!=degrees.len()||lower.iter().zip(&upper).any(|(a,b)|b<=a){return Err(err("Chebyshev mesh/dense allocation budget invalid"));}let f=&b.frame;let f_reference_m2_0_1=f.reference([-2, 0, 1])?;let output=canonical(DimExpr::Mul(vec![f.dims[2].clone(),DimExpr::Pow(Box::new(f.dims[0].clone()),Exp::Rat(Rational::new(-2,1)))]))?;let c=crate::rules::numerical_abi::native_grids::structured_ops::checked_snapshot(ctx,&source_program,budget,&vec![f.dims[0].clone();3],&output)?;let out=pde::solve_chebyshev_helmholtz(&degrees,&lower,&upper,sigma,|x|eval(&c,f,&x,f_reference_m2_0_1),&b.data).map_err(err);c.finish()?;let(x,u)=out?;let n=u.len();Ok(Dynamic(api::record(vec![("coordinates",coordinates(f,x.into_iter().flatten().collect(),vec![n,3])?),("field",payload(f,u,vec![n])?.0)])))})});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(i, BINDINGS, "numerical_vector")
}
static BINDINGS: &[&api::Binding] = &[
    &pde_fem1::BINDING,
    &pde_fem2::BINDING,
    &pde_fem3::BINDING,
    &pde_fem_adaptive::BINDING,
    &pde_fem_error::BINDING,
    &pde_fem_data::BINDING,
    &pde_dg_advection::BINDING,
    &pde_solve_fvm_conservation_1d::BINDING,
    &pde_solve_fvm_burgers_1d::BINDING,
    &pde_shallow::BINDING,
    &pde_spectral_helm::BINDING,
    &pde_spectral_etd::BINDING,
    &pde_cheb_helm::BINDING,
];
