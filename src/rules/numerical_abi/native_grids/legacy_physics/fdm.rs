//! Explicit Float specialization of legacy FDM grids and sampled initializer APIs.
use super::DimExpr;
use super::Dynamic;
use super::Faces;
use super::Frame;
use super::Handle;
use super::Input;
use super::Installer;
use super::LIMIT;
use super::NativeType;
use super::NewObject;
use super::OpError;
use super::RuleError;
use super::Value;
use super::api;
use super::array;
use super::canonical;
use super::derivative;
use super::dim;
use super::err;
use super::object_binding;
use super::quantity;
use super::vector;
use crate::sim::physics_fdm as fdm;
struct Fdm {
    data: fdm::FdmGrid<f64>,
    frame: Frame,
}
impl NativeType for Fdm {
    const KIND: &'static str = "rssn/native/legacy-fdm-grid/v1";
    const MAX_BYTES: usize = 1_048_576;

    fn validate(&self) -> Result<(), OpError> {
        self.frame.validate()?;
        if extent(&sizes(self.data.dimensions()))? != self.data.len()
            || self.data.as_slice().iter().any(|v| !v.is_finite())
        {
            return Err(err("FDM native data shape/finite contract violated"));
        }
        Ok(())
    }
}
fn sizes(d: &fdm::Dimensions) -> Vec<usize> {
    match d {
        | fdm::Dimensions::D1(n) => vec![*n],
        | fdm::Dimensions::D2(x, y) => vec![*x, *y],
        | fdm::Dimensions::D3(x, y, z) => vec![*x, *y, *z],
    }
}
fn extent(d: &[usize]) -> Result<usize, OpError> {
    if !(1..=3).contains(&d.len()) || d.iter().any(|n| *n > LIMIT) {
        return Err(err("FDM shape requires1..3 bounded axes"));
    }
    d.iter()
        .try_fold(1usize, |n, d| n.checked_mul(*d))
        .filter(|n| *n <= LIMIT)
        .ok_or_else(|| err("FDM grid allocation budget exceeded"))
}
fn axes(d: &[usize]) -> Result<fdm::Dimensions, OpError> {
    extent(d)?;
    Ok(match d {
        | [n] => fdm::Dimensions::D1(*n),
        | [x, y] => fdm::Dimensions::D2(*x, *y),
        | [x, y, z] => fdm::Dimensions::D3(*x, *y, *z),
        | _ => return Err(err("invalid shape")),
    })
}
fn result(
    data: fdm::FdmGrid<f64>,
    f: &Frame,
) -> NewObject<Fdm> {
    NewObject::new(
        Fdm {
            data,
            frame: f.clone(),
        },
        f.dims.to_vec(),
    )
}
fn measured_data(a: &Fdm) -> Result<Dynamic, OpError> {
    let data = a
        .data
        .as_slice()
        .iter()
        .map(|v| api::scale_product(*v, &[a.frame.value], &[]))
        .collect::<Result<Vec<_>, _>>()?;
    Ok(Dynamic(array(
        data,
        vec![a.data.len()],
        a.frame.dims[2].clone(),
    )?))
}
fn normalized_values(
    f: &Frame,
    v: &Value,
) -> Result<Vec<f64>, OpError> {
    vector(v, &f.dims[2])?
        .into_iter()
        .map(|v| api::scale_product(v, &[], &[f.value]))
        .collect()
}
fn controls(
    a: &Fdm,
    h: &Value,
    dt: &Value,
    steps: usize,
) -> Result<(Vec<usize>, Vec<f64>, f64), OpError> {
    let n = sizes(a.data.dimensions());
    if n.iter().any(|n| *n < 3)
        || steps > 65536
        || steps
            .checked_mul(a.data.len())
            .is_none_or(|v| v > 20_000_000)
    {
        return Err(err("FDM node/step work budget exceeded"));
    }
    let h = vector(h, &a.frame.dims[0])?
        .into_iter()
        .map(|v| api::scale_product(v, &[], &[a.frame.length]))
        .collect::<Result<Vec<_>, _>>()?;
    let dt = quantity(dt, &a.frame.dims[1])? / a.frame.time;
    if h.len() != n.len() || h.iter().any(|v| *v <= 0.0) || !dt.is_finite() || dt <= 0.0 {
        return Err(err("FDM spacing/step controls invalid"));
    }
    Ok((n, h, dt))
}
fn alpha(
    f: &Frame,
    v: &Value,
) -> Result<f64, OpError> {
    let v = quantity(v, &dim([2, -1, 0, 0])?)?;
    let a = api::scale_product(v, &[f.time], &[f.length, f.length])?;
    if !a.is_finite() || a < 0.0 {
        return Err(err("diffusivity must be nonnegative"));
    }
    Ok(a)
}
fn speed(
    f: &Frame,
    v: &Value,
) -> Result<f64, OpError> {
    let v = quantity(v, &dim([1, -1, 0, 0])?)?;
    let c = api::scale_product(v, &[f.time], &[f.length])?;
    if !c.is_finite() || c < 0.0 {
        return Err(err("wave speed must be nonnegative"));
    }
    Ok(c)
}
fn poisson_source(
    source: &Fdm,
    f: &Frame,
) -> Result<fdm::FdmGrid<f64>, OpError> {
    api::same_dimension(&source.frame.dims[2], &derivative(f.dims[2].clone(), -2)?)?;
    let data = source
        .data
        .as_slice()
        .iter()
        .map(|v| api::scale_product(*v, &[source.frame.value, f.length, f.length], &[f.value]))
        .collect::<Result<_, _>>()?;
    Ok(fdm::FdmGrid::from_data(
        data,
        source.data.dimensions().clone(),
    ))
}
object_binding!(new,"num_grid_legacy_fdm_new","crate::sim::physics_fdm::FdmGrid<T>::new","Explicit native Float FDM-grid and sampled-initializer specialization. Measured field/coefficient/time/spacing inputs are normalized by retained physical references. Grid, CFL, work and data bounds precede native work; BC callbacks retain native errors.",false,(frame:Handle<Frame>,shape:Vec<usize>)->NewObject<Fdm>,|ctx,frame,shape|{let _=ctx;frame.with(ctx,|f|Ok(result(fdm::FdmGrid::<f64>::new(axes(&shape)?),f)))});
object_binding!(with_value,"num_grid_legacy_fdm_with_value","crate::sim::physics_fdm::FdmGrid<T>::with_value","Explicit native Float FDM-grid and sampled-initializer specialization. Measured field/coefficient/time/spacing inputs are normalized by retained physical references. Grid, CFL, work and data bounds precede native work; BC callbacks retain native errors.",false,(frame:Handle<Frame>,shape:Vec<usize>,value:Input)->NewObject<Fdm>,|ctx,frame,shape,value|{let _=ctx;frame.with(ctx,|f|{let value=quantity(&value.0,&f.dims[2])?/f.value;Ok(result(fdm::FdmGrid::<f64>::with_value(axes(&shape)?,value),f))})});
object_binding!(from_data,"num_grid_legacy_fdm_from_data","crate::sim::physics_fdm::FdmGrid<T>::from_data","Explicit native Float FDM-grid and sampled-initializer specialization. Measured field/coefficient/time/spacing inputs are normalized by retained physical references. Grid, CFL, work and data bounds precede native work; BC callbacks retain native errors.",false,(frame:Handle<Frame>,shape:Vec<usize>,values:Input)->NewObject<Fdm>,|ctx,frame,shape,values|{let _=ctx;frame.with(ctx,|f|{let n=extent(&shape)?;let values=normalized_values(f,&values.0)?;if n!=values.len(){return Err(err("FDM shape/data length mismatch"));}Ok(result(fdm::FdmGrid::<f64>::from_data(values,axes(&shape)?),f))})});
object_binding!(as_slice,"num_grid_legacy_fdm_as_slice","crate::sim::physics_fdm::FdmGrid<T>::as_slice","Explicit native Float FDM-grid and sampled-initializer specialization. Measured field/coefficient/time/spacing inputs are normalized by retained physical references. Grid, CFL, work and data bounds precede native work; BC callbacks retain native errors.",false,(a:Handle<Fdm>)->Dynamic,|ctx,a|{let _=ctx;a.with(ctx,measured_data)});
object_binding!(dimensions,"num_grid_legacy_fdm_dimensions","crate::sim::physics_fdm::FdmGrid<T>::dimensions","Explicit native Float FDM-grid and sampled-initializer specialization. Measured field/coefficient/time/spacing inputs are normalized by retained physical references. Grid, CFL, work and data bounds precede native work; BC callbacks retain native errors.",false,(a:Handle<Fdm>)->Vec<usize>,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(sizes(a.data.dimensions())))});
object_binding!(len,"num_grid_legacy_fdm_len","crate::sim::physics_fdm::FdmGrid<T>::len","Explicit native Float FDM-grid and sampled-initializer specialization. Measured field/coefficient/time/spacing inputs are normalized by retained physical references. Grid, CFL, work and data bounds precede native work; BC callbacks retain native errors.",false,(a:Handle<Fdm>)->usize,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(a.data.len()))});
object_binding!(is_empty,"num_grid_legacy_fdm_is_empty","crate::sim::physics_fdm::FdmGrid<T>::is_empty","Explicit native Float FDM-grid and sampled-initializer specialization. Measured field/coefficient/time/spacing inputs are normalized by retained physical references. Grid, CFL, work and data bounds precede native work; BC callbacks retain native errors.",false,(a:Handle<Fdm>)->bool,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(a.data.is_empty()))});
object_binding!(as_mut_slice,"num_grid_legacy_fdm_as_mut_slice","crate::sim::physics_fdm::FdmGrid<T>::as_mut_slice","Explicit native Float FDM-grid and sampled-initializer specialization. Measured field/coefficient/time/spacing inputs are normalized by retained physical references. Grid, CFL, work and data bounds precede native work; BC callbacks retain native errors.",false,(a:Handle<Fdm>,values:Input)->bool,|ctx,a,values|{let _=ctx;a.with_mut(ctx,|a|{let values=normalized_values(&a.frame,&values.0)?;if values.len()!=a.data.len(){return Err(err("replacement slice length mismatch"));}a.data.as_mut_slice().copy_from_slice(&values);Ok(true)})});
object_binding!(solve_heat_equation_2d,"num_grid_legacy_fdm_solve_heat_equation_2d","crate::sim::physics_fdm::solve_heat_equation_2d","Explicit native Float FDM-grid and sampled-initializer specialization. Measured field/coefficient/time/spacing inputs are normalized by retained physical references. Grid, CFL, work and data bounds precede native work; BC callbacks retain native errors.",false,(a:Handle<Fdm>,spacing:Input,dt:Input,steps:usize,coefficient:Input)->NewObject<Fdm>,|ctx,a,spacing,dt,steps,coefficient|{let _=ctx;a.with(ctx,|a|{let(n,h,dt)=controls(a,&spacing.0,&dt.0,steps)?;if n.len()!=2{return Err(err("sampled initializer grid rank mismatch"));}let coefficient=alpha(&a.frame,&coefficient.0)?;if coefficient*dt*h.iter().map(|h|1.0/h/h).sum::<f64>()>0.5{return Err(err("explicit FDM stability bound exceeded"));}let data=fdm::solve_heat_equation_2d(&fdm::FdmSolverConfig2D{width:n[0],height:n[1],dx:h[0],dy:h[1],dt,steps},coefficient,|x,y|a.data.as_slice()[x+y*n[0]]);Ok(result(data,&a.frame))})});
object_binding!(solve_wave_equation_2d,"num_grid_legacy_fdm_solve_wave_equation_2d","crate::sim::physics_fdm::solve_wave_equation_2d","Explicit native Float FDM-grid and sampled-initializer specialization. Measured field/coefficient/time/spacing inputs are normalized by retained physical references. Grid, CFL, work and data bounds precede native work; BC callbacks retain native errors.",false,(a:Handle<Fdm>,spacing:Input,dt:Input,steps:usize,coefficient:Input)->NewObject<Fdm>,|ctx,a,spacing,dt,steps,coefficient|{let _=ctx;a.with(ctx,|a|{let(n,h,dt)=controls(a,&spacing.0,&dt.0,steps)?;if n.len()!=2{return Err(err("sampled initializer grid rank mismatch"));}let coefficient=speed(&a.frame,&coefficient.0)?;if coefficient*coefficient*dt*dt*h.iter().map(|h|1.0/h/h).sum::<f64>()>1.0{return Err(err("explicit FDM stability bound exceeded"));}let data=fdm::solve_wave_equation_2d(&fdm::FdmSolverConfig2D{width:n[0],height:n[1],dx:h[0],dy:h[1],dt,steps},coefficient,|x,y|a.data.as_slice()[x+y*n[0]]);Ok(result(data,&a.frame))})});
object_binding!(solve_wave_equation_3d,"num_grid_legacy_fdm_solve_wave_equation_3d","crate::sim::physics_fdm::solve_wave_equation_3d","Explicit native Float FDM-grid and sampled-initializer specialization. Measured field/coefficient/time/spacing inputs are normalized by retained physical references. Grid, CFL, work and data bounds precede native work; BC callbacks retain native errors.",false,(a:Handle<Fdm>,spacing:Input,dt:Input,steps:usize,coefficient:Input)->NewObject<Fdm>,|ctx,a,spacing,dt,steps,coefficient|{let _=ctx;a.with(ctx,|a|{let(n,h,dt)=controls(a,&spacing.0,&dt.0,steps)?;if n.len()!=3{return Err(err("sampled initializer grid rank mismatch"));}let coefficient=speed(&a.frame,&coefficient.0)?;if coefficient*coefficient*dt*dt*h.iter().map(|h|1.0/h/h).sum::<f64>()>1.0{return Err(err("explicit FDM stability bound exceeded"));}let data=fdm::solve_wave_equation_3d(&fdm::FdmSolverConfig3D{width:n[0],height:n[1],depth:n[2],dx:h[0],dy:h[1],dz:h[2],dt,steps},coefficient,|x,y,z|a.data.as_slice()[x+n[0]*(y+n[1]*z)]);Ok(result(data,&a.frame))})});
object_binding!(heat_bc,"num_grid_legacy_fdm_heat_bc","crate::sim::physics_fdm::solve_heat_equation_bc","Explicit native Float FDM-grid and sampled-initializer specialization. Measured field/coefficient/time/spacing inputs are normalized by retained physical references. Grid, CFL, work and data bounds precede native work; BC callbacks retain native errors.",false,(a:Handle<Fdm>,spacing:Input,diffusivity:Input,dt:Input,steps:usize,bcs:Handle<Faces>)->NewObject<Fdm>,|ctx,a,spacing,diffusivity,dt,steps,bcs|{let _=ctx;a.with(ctx,|a|bcs.with(ctx,|bc|{a.frame.same(&bc.frame)?;let(n,h,dt)=controls(a,&spacing.0,&dt.0,steps)?;if steps==0{return Err(err("implicit heat requires at least1step"));}let alpha=alpha(&a.frame,&diffusivity.0)?;let out=fdm::solve_heat_equation_bc(&a.data,&h,alpha,dt,steps,&bc.data);bc.finish_for(n.len())?;Ok(result(out.map_err(err)?,&a.frame))}))});
object_binding!(wave_bc,"num_grid_legacy_fdm_wave_bc","crate::sim::physics_fdm::solve_wave_equation_bc","Explicit native Float FDM-grid and sampled-initializer specialization. Measured field/coefficient/time/spacing inputs are normalized by retained physical references. Grid, CFL, work and data bounds precede native work; BC callbacks retain native errors.",false,(a:Handle<Fdm>,initial_velocity:Handle<Fdm>,spacing:Input,wave_speed:Input,dt:Input,steps:usize,bcs:Handle<Faces>)->NewObject<Fdm>,|ctx,a,initial_velocity,spacing,wave_speed,dt,steps,bcs|{let _=ctx;a.with(ctx,|a|initial_velocity.with(ctx,|v|bcs.with(ctx,|bc|{a.frame.same(&bc.frame)?;let(n,h,dt)=controls(a,&spacing.0,&dt.0,steps)?;if steps==0||sizes(v.data.dimensions())!=n{return Err(err("wave initial shape/steps mismatch"));}api::same_dimension(&v.frame.dims[2],&canonical(DimExpr::Mul(vec![a.frame.dims[2].clone(),a.frame.dims[1].clone().pow(-1)]))?)?;let velocity=v.data.as_slice().iter().map(|v0|api::scale_product(*v0,&[v.frame.value,a.frame.time],&[a.frame.value])).collect::<Result<Vec<_>,_>>()?;let velocity=fdm::FdmGrid::from_data(velocity,a.data.dimensions().clone());let speed=speed(&a.frame,&wave_speed.0)?;if speed*speed*dt*dt*h.iter().map(|h|1.0/h/h).sum::<f64>()>1.0{return Err(err("wave CFL bound exceeded"));}let out=fdm::solve_wave_equation_bc(&a.data,&velocity,&h,speed,dt,steps,&bc.data);bc.finish_for(n.len())?;Ok(result(out.map_err(err)?,&a.frame))})))});
object_binding!(poisson_bc,"num_grid_legacy_fdm_poisson_bc","crate::sim::physics_fdm::solve_poisson_bc","Explicit native Float FDM-grid and sampled-initializer specialization. Measured field/coefficient/time/spacing inputs are normalized by retained physical references. Grid, CFL, work and data bounds precede native work; BC callbacks retain native errors.",false,(source:Handle<Fdm>,spacing:Input,sigma:Input,bcs:Handle<Faces>)->NewObject<Fdm>,|ctx,source,spacing,sigma,bcs|{let _=ctx;source.with(ctx,|s|bcs.with(ctx,|bc|{let f=&bc.frame;let data=poisson_source(s,f)?;let n=sizes(data.dimensions());if n.iter().any(|n|*n<3)||data.len()>1024{return Err(err("Poisson BC solve size budget exceeded"));}let h=vector(&spacing.0,&f.dims[0])?.into_iter().map(|h|h/f.length).collect::<Vec<_>>();if h.len()!=n.len()||h.iter().any(|h|*h<=0.0){return Err(err("spacing shape/positivity mismatch"));}let sigma=api::scale_product(quantity(&sigma.0,&dim([-2,0,0,0])?)?,&[f.length,f.length],&[])?;let out=fdm::solve_poisson_bc(&data,&h,sigma,&bc.data);bc.finish_for(n.len())?;Ok(result(out.map_err(err)?,f))}))});
object_binding!(poisson_2d,"num_grid_legacy_fdm_poisson_2d","crate::sim::physics_fdm::solve_poisson_2d","Explicit native Float FDM-grid and sampled-initializer specialization. Measured field/coefficient/time/spacing inputs are normalized by retained physical references. Grid, CFL, work and data bounds precede native work; BC callbacks retain native errors.",false,(source:Handle<Fdm>,output_frame:Handle<Frame>,spacing:Input,omega:f64,max_iter:usize,tolerance:Input)->NewObject<Fdm>,|ctx,source,output_frame,spacing,omega,max_iter,tolerance|{let _=ctx;source.with(ctx,|s|output_frame.with(ctx,|f|{let data=poisson_source(s,f)?;let n=sizes(data.dimensions());let h=vector(&spacing.0,&f.dims[0])?.into_iter().map(|h|h/f.length).collect::<Vec<_>>();let tol=quantity(&tolerance.0,&f.dims[2])?/f.value;if n.len()!=2||n.iter().any(|n|*n<3)||h.len()!=2||h.iter().any(|h|*h<=0.0)||omega<=0.0||omega>=2.0||tol<=0.0||max_iter==0||max_iter.checked_mul(data.len()).is_none_or(|n|n>20_000_000){return Err(err("Poisson2D controls/work budget invalid"));}let cfg=fdm::PoissonSolverConfig2D{width:n[0],height:n[1],dx:h[0],dy:h[1],omega,max_iter,tolerance:tol};Ok(result(fdm::solve_poisson_2d(&cfg,&data),f))}))});
object_binding!(poisson_nd,"num_grid_legacy_fdm_poisson_nd","crate::sim::physics_fdm::solve_poisson_nd","Explicit native Float FDM-grid and sampled-initializer specialization. Measured field/coefficient/time/spacing inputs are normalized by retained physical references. Grid, CFL, work and data bounds precede native work; BC callbacks retain native errors.",false,(source:Handle<Fdm>,output_frame:Handle<Frame>,spacing:Input,omega:f64,max_iter:usize,tolerance:Input)->Dynamic,|ctx,source,output_frame,spacing,omega,max_iter,tolerance|{let _=ctx;source.with(ctx,|s|output_frame.with(ctx,|f|{let data=poisson_source(s,f)?;let n=sizes(data.dimensions());let h=vector(&spacing.0,&f.dims[0])?.into_iter().map(|h|h/f.length).collect::<Vec<_>>();let tol=quantity(&tolerance.0,&f.dims[2])?/f.value;if n.iter().any(|n|*n<3)||h.len()!=n.len()||h.iter().any(|h|*h<=0.0)||omega<=0.0||omega>=2.0||tol<=0.0||max_iter==0||max_iter.checked_mul(data.len()).is_none_or(|n|n>20_000_000){return Err(err("PoissonND controls/work budget invalid"));}let cfg=fdm::PoissonSolverConfigNd{dims:n.clone(),spacing:h,omega,max_iter,tolerance:tol};let data=fdm::solve_poisson_nd(&cfg,data.as_slice()).map_err(err)?.into_iter().map(|v|api::scale_product(v,&[f.value],&[])).collect::<Result<Vec<_>,_>>()?;Ok(Dynamic(array(data,vec![extent(&n)?],f.dims[2].clone())?))}))});
object_binding!(scenario_heat,"num_grid_legacy_fdm_scenario_heat","crate::sim::physics_fdm::simulate_2d_heat_conduction_scenario","Explicit native Float FDM-grid and sampled-initializer specialization. Measured field/coefficient/time/spacing inputs are normalized by retained physical references. Grid, CFL, work and data bounds precede native work; BC callbacks retain native errors.",false,()->NewObject<Fdm>,|ctx,|{let _=ctx;{let f=Frame{length:1.0,time:1.0,value:1.0,dims:[DimExpr::base("Length"),DimExpr::base("Time"),DimExpr::base("Temperature")]};Ok(result(fdm::simulate_2d_heat_conduction_scenario(),&f))}});
object_binding!(scenario_wave,"num_grid_legacy_fdm_scenario_wave","crate::sim::physics_fdm::simulate_2d_wave_propagation_scenario","Explicit native Float FDM-grid and sampled-initializer specialization. Measured field/coefficient/time/spacing inputs are normalized by retained physical references. Grid, CFL, work and data bounds precede native work; BC callbacks retain native errors.",false,()->NewObject<Fdm>,|ctx,|{let _=ctx;{let f=Frame{length:1.0,time:1.0,value:1.0,dims:[DimExpr::base("Length"),DimExpr::base("Time"),DimExpr::base("Length")]};Ok(result(fdm::simulate_2d_wave_propagation_scenario(),&f))}});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &new::BINDING,
            &with_value::BINDING,
            &from_data::BINDING,
            &as_slice::BINDING,
            &dimensions::BINDING,
            &len::BINDING,
            &is_empty::BINDING,
            &as_mut_slice::BINDING,
            &solve_heat_equation_2d::BINDING,
            &solve_wave_equation_2d::BINDING,
            &solve_wave_equation_3d::BINDING,
            &heat_bc::BINDING,
            &wave_bc::BINDING,
            &poisson_bc::BINDING,
            &poisson_2d::BINDING,
            &poisson_nd::BINDING,
            &scenario_heat::BINDING,
            &scenario_wave::BINDING,
        ],
        "numerical_vector",
    )
}
