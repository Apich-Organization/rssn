//! Native legacy continuum APIs with measured SI carriers and checked f64 specialization.
mod fdm;
mod fvm;
use super::DimExpr;
use super::Dynamic;
use super::Exp;
use super::Faces;
use super::Frame;
use super::Handle;
use super::Input;
use super::Installer;
use super::LIMIT;
use super::NativeType;
use super::NewObject;
use super::OpError;
use super::Rational;
use super::RuleError;
use super::TensorValue;
use super::TypeExpr;
use super::Value;
use super::api;
use super::canonical;
use super::compiled_snapshot;
use super::err;
use super::object_binding;
use crate::kernels::matrix::Matrix;
use crate::rules::numerical_abi::codecs::Arg;
use crate::sim::physics_cfd as cfd;
fn dim(p: [i32; 4]) -> Result<DimExpr, OpError> {
    canonical(DimExpr::Mul(
        ["Length", "Time", "Mass", "Temperature"]
            .into_iter()
            .zip(p)
            .filter(|(_, p)| *p != 0)
            .map(|(d, p)| DimExpr::base(d).pow(Exp::Rat(Rational::int(i64::from(p)))))
            .collect(),
    ))
}
fn quantity(
    v: &Value,
    d: &DimExpr,
) -> Result<f64, OpError> {
    api::same_dimension(&api::dimension(v)?, d)?;
    crate::rules::numerical_abi::numeric(v, 0)
}
fn vector(
    v: &Value,
    d: &DimExpr,
) -> Result<Vec<f64>, OpError> {
    let Value::Tensor(t) = v else {
        return Err(err("expected a measured tensor"));
    };
    if t.shape != [t.data.len()] || t.data.len() > LIMIT {
        return Err(err("continuum field allocation budget exceeded"));
    }
    api::same_dimension(&api::dimension(v)?, d)?;
    t.data
        .iter()
        .map(|v| crate::rules::numerical_abi::numeric(v, 0))
        .collect()
}
fn array(
    data: Vec<f64>,
    shape: Vec<usize>,
    d: DimExpr,
) -> Result<Value, OpError> {
    if shape.iter().try_fold(1usize, |n, v| n.checked_mul(*v)) != Some(data.len())
        || data.len() > LIMIT
    {
        return Err(err("continuum output shape/budget mismatch"));
    }
    let data = data
        .into_iter()
        .map(|x| api::quantity(x, d.clone()))
        .collect::<Result<_, _>>()?;
    Ok(Value::Tensor(TensorValue {
        shape,
        elem: TypeExpr::Float.with_dim(d),
        data,
    }))
}
fn matrix(v: &Value) -> Result<(Matrix<f64>, DimExpr), OpError> {
    let Value::Tensor(t) = v else {
        return Err(err("expected rank2 measured matrix"));
    };
    if t.shape.len() != 2
        || t.shape.contains(&0)
        || t.shape.iter().try_fold(1usize, |n, d| n.checked_mul(*d)) != Some(t.data.len())
        || t.data.len() > LIMIT
    {
        return Err(err("matrix shape/allocation budget exceeded"));
    }
    let d = api::dimension(v)?;
    let data = t
        .data
        .iter()
        .map(|v| crate::rules::numerical_abi::numeric(v, 0))
        .collect::<Result<Vec<_>, _>>()?;
    Ok((Matrix::new(t.shape[0], t.shape[1], data), d))
}
fn matout(
    m: Matrix<f64>,
    d: DimExpr,
) -> Result<Dynamic, OpError> {
    let shape = vec![m.rows(), m.cols()];
    Ok(Dynamic(array(m.into_data(), shape, d)?))
}
fn derivative(
    d: DimExpr,
    p: i32,
) -> Result<DimExpr, OpError> {
    canonical(DimExpr::Mul(vec![
        d,
        DimExpr::base("Length").pow(Exp::Rat(Rational::int(i64::from(p)))),
    ]))
}
fn spacing(
    dx: &Value,
    dy: &Value,
) -> Result<(f64, f64), OpError> {
    let d = DimExpr::base("Length");
    let (x, y) = (quantity(dx, &d)?, quantity(dy, &d)?);
    if x <= 0.0 || y <= 0.0 {
        return Err(err("spacings must be finite and positive"));
    }
    Ok((x, y))
}
fn grid(m: &Matrix<f64>) -> Result<(), OpError> {
    if m.rows() < 3 || m.cols() < 3 {
        return Err(err("finite differences require at least3nodes per axis"));
    }
    Ok(())
}
fn same(
    a: &Matrix<f64>,
    b: &Matrix<f64>,
) -> Result<(), OpError> {
    if a.rows() != b.rows() || a.cols() != b.cols() {
        Err(err("matrix grid extents differ"))
    } else {
        Ok(())
    }
}
fn iterate(
    m: &Matrix<f64>,
    n: usize,
) -> Result<(), OpError> {
    grid(m)?;
    if n == 0
        || n > 65536
        || n.checked_mul(m.rows())
            .and_then(|n| n.checked_mul(m.cols()))
            .is_none_or(|n| n > 20_000_000)
    {
        Err(err("stationary solver work budget exceeded"))
    } else {
        Ok(())
    }
}
struct Fluid(cfd::FluidProperties);
impl NativeType for Fluid {
    const KIND: &'static str = "rssn/native/fluid-properties/v1";
    const MAX_BYTES: usize = 4096;

    fn validate(&self) -> Result<(), OpError> {
        let f = &self.0;
        if [
            f.density,
            f.dynamic_viscosity,
            f.thermal_conductivity,
            f.specific_heat,
        ]
        .iter()
        .any(|x| !x.is_finite())
            || f.density <= 0.0
            || f.specific_heat <= 0.0
            || f.dynamic_viscosity < 0.0
            || f.thermal_conductivity <= 0.0
        {
            return Err(err(
                "fluid density/heatcapacity/conductivity must be positive and viscosity nonnegative",
            ));
        }
        Ok(())
    }
}
object_binding!(fluid_new,"num_grid_legacy_fluid_new","crate::sim::physics_cfd::FluidProperties::new","Measured native continuum adapter. Float matrices and mutable borrows use explicit checked tensor inputs and returned values; SI dimensions, grid shape, work/storage bounds and original native methods are preserved.",false,(density:Input,viscosity:Input,conductivity:Input,heat_capacity:Input)->NewObject<Fluid>,|ctx,density,viscosity,conductivity,heat_capacity|{let _=ctx;Ok(NewObject::new(Fluid(cfd::FluidProperties::new(quantity(&density.0,&dim([-3,0,1,0])?)?,quantity(&viscosity.0,&dim([-1,-1,1,0])?)?,quantity(&conductivity.0,&dim([1,-3,1,-1])?)?,quantity(&heat_capacity.0,&dim([2,-2,0,-1])?)?)),vec![]))});
object_binding!(fluid_air,"num_grid_legacy_fluid_air","crate::sim::physics_cfd::FluidProperties::air","Measured native continuum adapter. Float matrices and mutable borrows use explicit checked tensor inputs and returned values; SI dimensions, grid shape, work/storage bounds and original native methods are preserved.",false,()->NewObject<Fluid>,|ctx,|{let _=ctx;Ok(NewObject::new(Fluid(cfd::FluidProperties::air()),vec![]))});
object_binding!(fluid_water,"num_grid_legacy_fluid_water","crate::sim::physics_cfd::FluidProperties::water","Measured native continuum adapter. Float matrices and mutable borrows use explicit checked tensor inputs and returned values; SI dimensions, grid shape, work/storage bounds and original native methods are preserved.",false,()->NewObject<Fluid>,|ctx,|{let _=ctx;Ok(NewObject::new(Fluid(cfd::FluidProperties::water()),vec![]))});
object_binding!(fluid_kinematic_viscosity,"num_grid_legacy_fluid_kinematic_viscosity","crate::sim::physics_cfd::FluidProperties::kinematic_viscosity","Measured native continuum adapter. Float matrices and mutable borrows use explicit checked tensor inputs and returned values; SI dimensions, grid shape, work/storage bounds and original native methods are preserved.",false,(a:Handle<Fluid>)->Dynamic,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(Dynamic(api::quantity(a.0.kinematic_viscosity(),dim([2,-1,0,0])?)?)))});
object_binding!(fluid_thermal_diffusivity,"num_grid_legacy_fluid_thermal_diffusivity","crate::sim::physics_cfd::FluidProperties::thermal_diffusivity","Measured native continuum adapter. Float matrices and mutable borrows use explicit checked tensor inputs and returned values; SI dimensions, grid shape, work/storage bounds and original native methods are preserved.",false,(a:Handle<Fluid>)->Dynamic,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(Dynamic(api::quantity(a.0.thermal_diffusivity(),dim([2,-1,0,0])?)?)))});
object_binding!(fluid_prandtl,"num_grid_legacy_fluid_prandtl","crate::sim::physics_cfd::FluidProperties::prandtl_number","Measured native continuum adapter. Float matrices and mutable borrows use explicit checked tensor inputs and returned values; SI dimensions, grid shape, work/storage bounds and original native methods are preserved.",false,(a:Handle<Fluid>)->f64,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(a.0.prandtl_number()))});
object_binding!(compute_vorticity,"num_grid_legacy_compute_vorticity","crate::sim::physics_cfd::compute_vorticity","Measured native continuum adapter. Float matrices and mutable borrows use explicit checked tensor inputs and returned values; SI dimensions, grid shape, work/storage bounds and original native methods are preserved.",false,(u:Input,v:Input,dx:Input,dy:Input)->Dynamic,|ctx,u,v,dx,dy|{let _=ctx;{let(u,d)=matrix(&u.0)?;let(v,e)=matrix(&v.0)?;grid(&u)?;same(&u,&v)?;api::same_dimension(&d,&e)?;let(dx,dy)=spacing(&dx.0,&dy.0)?;matout(cfd::compute_vorticity(&u,&v,dx,dy),derivative(d,-1)?)}});
object_binding!(compute_divergence,"num_grid_legacy_compute_divergence","crate::sim::physics_cfd::compute_divergence","Measured native continuum adapter. Float matrices and mutable borrows use explicit checked tensor inputs and returned values; SI dimensions, grid shape, work/storage bounds and original native methods are preserved.",false,(u:Input,v:Input,dx:Input,dy:Input)->Dynamic,|ctx,u,v,dx,dy|{let _=ctx;{let(u,d)=matrix(&u.0)?;let(v,e)=matrix(&v.0)?;grid(&u)?;same(&u,&v)?;api::same_dimension(&d,&e)?;let(dx,dy)=spacing(&dx.0,&dy.0)?;matout(cfd::compute_divergence(&u,&v,dx,dy),derivative(d,-1)?)}});
object_binding!(compute_gradient,"num_grid_legacy_compute_gradient","crate::sim::physics_cfd::compute_gradient","Measured native continuum adapter. Float matrices and mutable borrows use explicit checked tensor inputs and returned values; SI dimensions, grid shape, work/storage bounds and original native methods are preserved.",false,(field:Input,dx:Input,dy:Input)->Dynamic,|ctx,field,dx,dy|{let _=ctx;{let(m,d)=matrix(&field.0)?;grid(&m)?;let(dx,dy)=spacing(&dx.0,&dy.0)?;let(a,b)=cfd::compute_gradient(&m,dx,dy);let d=derivative(d,-1)?;Ok(Dynamic(api::record(vec![("x",matout(a,d.clone())?.0),("y",matout(b,d)?.0)])))}});
object_binding!(velocity_from_stream_function,"num_grid_legacy_velocity_from_stream_function","crate::sim::physics_cfd::velocity_from_stream_function","Measured native continuum adapter. Float matrices and mutable borrows use explicit checked tensor inputs and returned values; SI dimensions, grid shape, work/storage bounds and original native methods are preserved.",false,(field:Input,dx:Input,dy:Input)->Dynamic,|ctx,field,dx,dy|{let _=ctx;{let(m,d)=matrix(&field.0)?;grid(&m)?;let(dx,dy)=spacing(&dx.0,&dy.0)?;let(a,b)=cfd::velocity_from_stream_function(&m,dx,dy);let d=derivative(d,-1)?;Ok(Dynamic(api::record(vec![("x",matout(a,d.clone())?.0),("y",matout(b,d)?.0)])))}});
object_binding!(compute_laplacian,"num_grid_legacy_compute_laplacian","crate::sim::physics_cfd::compute_laplacian","Measured native continuum adapter. Float matrices and mutable borrows use explicit checked tensor inputs and returned values; SI dimensions, grid shape, work/storage bounds and original native methods are preserved.",false,(field:Input,dx:Input,dy:Input)->Dynamic,|ctx,field,dx,dy|{let _=ctx;{let(m,d)=matrix(&field.0)?;grid(&m)?;let(dx,dy)=spacing(&dx.0,&dy.0)?;matout(cfd::compute_laplacian(&m,dx,dy),derivative(d,-2)?)}});
object_binding!(compute_stream_function,"num_grid_legacy_compute_stream_function","crate::sim::physics_cfd::compute_stream_function","Measured native continuum adapter. Float matrices and mutable borrows use explicit checked tensor inputs and returned values; SI dimensions, grid shape, work/storage bounds and original native methods are preserved.",false,(vorticity:Input,dx:Input,dy:Input,max_iter:usize,tolerance:Input)->Dynamic,|ctx,vorticity,dx,dy,max_iter,tolerance|{let _=ctx;{let(m,d)=matrix(&vorticity.0)?;iterate(&m,max_iter)?;let d=derivative(d,2)?;let tol=quantity(&tolerance.0,&d)?;if tol<=0.0{return Err(err("stream function tolerance must be positive"));}let(dx,dy)=spacing(&dx.0,&dy.0)?;matout(cfd::compute_stream_function(&m,dx,dy,max_iter,tol),d)}});
object_binding!(l2_norm,"num_grid_legacy_l2_norm","crate::sim::physics_cfd::l2_norm","Measured native continuum adapter. Float matrices and mutable borrows use explicit checked tensor inputs and returned values; SI dimensions, grid shape, work/storage bounds and original native methods are preserved.",false,(field:Input)->Dynamic,|ctx,field|{let _=ctx;{let(m,d)=matrix(&field.0)?;Ok(Dynamic(api::quantity(cfd::l2_norm(&m),d)?))}});
object_binding!(max_abs,"num_grid_legacy_max_abs","crate::sim::physics_cfd::max_abs","Measured native continuum adapter. Float matrices and mutable borrows use explicit checked tensor inputs and returned values; SI dimensions, grid shape, work/storage bounds and original native methods are preserved.",false,(field:Input)->Dynamic,|ctx,field|{let _=ctx;{let(m,d)=matrix(&field.0)?;Ok(Dynamic(api::quantity(cfd::max_abs(&m),d)?))}});
object_binding!(max_velocity_magnitude,"num_grid_legacy_max_velocity_magnitude","crate::sim::physics_cfd::max_velocity_magnitude","Measured native continuum adapter. Float matrices and mutable borrows use explicit checked tensor inputs and returned values; SI dimensions, grid shape, work/storage bounds and original native methods are preserved.",false,(u:Input,v:Input)->Dynamic,|ctx,u,v|{let _=ctx;{let(u,d)=matrix(&u.0)?;let(v,e)=matrix(&v.0)?;same(&u,&v)?;api::same_dimension(&d,&e)?;Ok(Dynamic(api::quantity(cfd::max_velocity_magnitude(&u,&v),d)?))}});
object_binding!(apply_dirichlet_bc,"num_grid_legacy_apply_dirichlet_bc","crate::sim::physics_cfd::apply_dirichlet_bc","Measured native continuum adapter. Float matrices and mutable borrows use explicit checked tensor inputs and returned values; SI dimensions, grid shape, work/storage bounds and original native methods are preserved.",false,(field:Input,boundary:Input)->Dynamic,|ctx,field,boundary|{let _=ctx;{let(mut m,d)=matrix(&field.0)?;let x=quantity(&boundary.0,&d)?;cfd::apply_dirichlet_bc(&mut m,x);matout(m,d)}});
object_binding!(apply_neumann_bc,"num_grid_legacy_apply_neumann_bc","crate::sim::physics_cfd::apply_neumann_bc","Measured native continuum adapter. Float matrices and mutable borrows use explicit checked tensor inputs and returned values; SI dimensions, grid shape, work/storage bounds and original native methods are preserved.",false,(field:Input)->Dynamic,|ctx,field|{let _=ctx;{let(mut m,d)=matrix(&field.0)?;grid(&m)?;cfd::apply_neumann_bc(&mut m);matout(m,d)}});
object_binding!(poisson_jacobi,"num_grid_legacy_poisson_jacobi","crate::sim::physics_cfd::solve_poisson_2d_jacobi","Measured native continuum adapter. Float matrices and mutable borrows use explicit checked tensor inputs and returned values; SI dimensions, grid shape, work/storage bounds and original native methods are preserved.",false,(source:Input,initial:Input,dx:Input,dy:Input,max_iter:usize,tolerance:Input)->Dynamic,|ctx,source,initial,dx,dy,max_iter,tolerance|{let _=ctx;{let(f,fd)=matrix(&source.0)?;let(u,d)=matrix(&initial.0)?;same(&f,&u)?;iterate(&u,max_iter)?;api::same_dimension(&fd,&derivative(d.clone(),-2)?)?;let(dx,dy)=spacing(&dx.0,&dy.0)?;let tol=quantity(&tolerance.0,&d)?;if tol<=0.0{return Err(err("Poisson tolerance must be positive"));}matout(cfd::solve_poisson_2d_jacobi(&f,&u,dx,dy,max_iter,tol),d)}});
object_binding!(poisson_gauss_seidel,"num_grid_legacy_poisson_gauss_seidel","crate::sim::physics_cfd::solve_poisson_2d_gauss_seidel","Measured native continuum adapter. Float matrices and mutable borrows use explicit checked tensor inputs and returned values; SI dimensions, grid shape, work/storage bounds and original native methods are preserved.",false,(source:Input,initial:Input,dx:Input,dy:Input,max_iter:usize,tolerance:Input)->Dynamic,|ctx,source,initial,dx,dy,max_iter,tolerance|{let _=ctx;{let(f,fd)=matrix(&source.0)?;let(u,d)=matrix(&initial.0)?;same(&f,&u)?;iterate(&u,max_iter)?;api::same_dimension(&fd,&derivative(d.clone(),-2)?)?;let(dx,dy)=spacing(&dx.0,&dy.0)?;let tol=quantity(&tolerance.0,&d)?;if tol<=0.0{return Err(err("Poisson tolerance must be positive"));}matout(cfd::solve_poisson_2d_gauss_seidel(&f,&u,dx,dy,max_iter,tol),d)}});
object_binding!(poisson_sor,"num_grid_legacy_poisson_sor","crate::sim::physics_cfd::solve_poisson_2d_sor","Measured native continuum adapter. Float matrices and mutable borrows use explicit checked tensor inputs and returned values; SI dimensions, grid shape, work/storage bounds and original native methods are preserved.",false,(source:Input,initial:Input,dx:Input,dy:Input,max_iter:usize,tolerance:Input,omega:f64)->Dynamic,|ctx,source,initial,dx,dy,max_iter,tolerance,omega|{let _=ctx;{let(f,fd)=matrix(&source.0)?;let(u,d)=matrix(&initial.0)?;same(&f,&u)?;iterate(&u,max_iter)?;api::same_dimension(&fd,&derivative(d.clone(),-2)?)?;let(dx,dy)=spacing(&dx.0,&dy.0)?;let tol=quantity(&tolerance.0,&d)?;if tol<=0.0{return Err(err("Poisson tolerance must be positive"));}
if !(0.0..2.0).contains(&omega)||omega==0.0{return Err(err("SOR relaxation must be in(0,2)"));}matout(cfd::solve_poisson_2d_sor(&f,&u,dx,dy,omega,max_iter,tol),d)}});
object_binding!(cavity_simple,"num_grid_legacy_cavity_simple","crate::sim::physics_cfd::lid_driven_cavity_simple","Measured native continuum adapter. Float matrices and mutable borrows use explicit checked tensor inputs and returned values; SI dimensions, grid shape, work/storage bounds and original native methods are preserved.",false,(nx:usize,ny:usize,re:f64,lid_velocity:Input,steps:usize,dt:Input)->Dynamic,|ctx,nx,ny,re,lid_velocity,steps,dt|{let _=ctx;{let velocity=quantity(&lid_velocity.0,&dim([1,-1,0,0])?)?;let dt=quantity(&dt.0,&dim([0,1,0,0])?)?;if nx<3||ny<3||nx.checked_mul(ny).and_then(|n|n.checked_mul(steps)).is_none_or(|n|n>100000)||re<=0.0||velocity<=0.0||dt<=0.0{return Err(err("unit-metre cavity controls/work budget invalid"));}let h=1.0/(nx.max(ny)-1) as f64;if dt*(velocity/h+4.0*velocity/re/h/h)>0.5{return Err(err("cavity explicit-step stability bound exceeded"));}let(psi,omega)=cfd::lid_driven_cavity_simple(nx,ny,re,velocity,steps,dt);Ok(Dynamic(api::record(vec![("stream_function",matout(psi,dim([2,-1,0,0])?)?.0),("vorticity",matout(omega,dim([0,-1,0,0])?)?.0)])))}});
object_binding!(cavity_mac,"num_grid_legacy_cavity_mac","crate::sim::physics_cfd::lid_driven_cavity_mac","Measured native continuum adapter. Float matrices and mutable borrows use explicit checked tensor inputs and returned values; SI dimensions, grid shape, work/storage bounds and original native methods are preserved.",false,(references:Handle<crate::rules::numerical_abi::native_grids::mechanics::References>,dimension:usize,n:usize,re:f64,end:Input,dt:Input)->NewObject<crate::rules::numerical_abi::native_grids::navier::Ns>,|ctx,references,dimension,n,re,end,dt|{let _=ctx;references.with(ctx,|r|{let end=quantity(&end.0,&dim([0,1,0,0])?)?/r.scales[1];let dt=quantity(&dt.0,&dim([0,1,0,0])?)?/r.scales[1];if !(2..=3).contains(&dimension)||n<3||re<=0.0||dt<=0.0||end<0.0||n.checked_pow(dimension as u32).is_none_or(|cells|cells>512||((end/dt).round()*cells as f64*cells as f64)>10000000.0){return Err(err("MAC cavity shape/time/work bounds invalid"));}Ok(crate::rules::numerical_abi::native_grids::navier::from_native(cfd::lid_driven_cavity_mac(dimension,n,re,end,dt).map_err(err)?,r.clone()))})});
object_binding!(burgers_1,"num_grid_legacy_burgers_1","crate::sim::physics_cfd::solve_burgers_1d","Measured native continuum adapter. Float matrices and mutable borrows use explicit checked tensor inputs and returned values; SI dimensions, grid shape, work/storage bounds and original native methods are preserved.",false,(initial:Input,dims:Vec<usize>,viscosity:Input,spacing:Input,dt:Input,steps:usize)->Dynamic,|ctx,initial,dims,viscosity,spacing,dt,steps|{let _=ctx;{let u=vector(&initial.0,&dim([1,-1,0,0])?)?;let h=vector(&spacing.0,&dim([1,0,0,0])?)?;let nu=quantity(&viscosity.0,&dim([2,-1,0,0])?)?;let dt=quantity(&dt.0,&dim([0,1,0,0])?)?;if dims.len()!=1||h.len()!=dims.len()||dims.iter().any(|n|*n<3)||dims.iter().try_fold(1usize,|n,d|n.checked_mul(*d))!=Some(u.len())||h.iter().any(|h|*h<=0.0)||dt<=0.0||nu<0.0||(if steps<=65536{Some(steps+1)}else{None}).and_then(|steps|u.len().checked_mul(steps)).is_none_or(|n|n>LIMIT){return Err(err("Burgers physical shape/stability/storage controls invalid"));}let speed=u.iter().copied().map(f64::abs).fold(0.0,f64::max);if dt*h.iter().map(|h|speed/h+2.0*nu/h/h).sum::<f64>()>0.5{return Err(err("Burgers explicit-step stability bound exceeded"));}let result=cfd::solve_burgers_1d(&u,nu,h[0],dt,steps);let rows=result.len();Ok(Dynamic(array(result.into_iter().flatten().collect(),vec![rows,u.len()],dim([1,-1,0,0])?)?))}});
object_binding!(burgers_2,"num_grid_legacy_burgers_2","crate::sim::physics_cfd::solve_burgers_2d","Measured native continuum adapter. Float matrices and mutable borrows use explicit checked tensor inputs and returned values; SI dimensions, grid shape, work/storage bounds and original native methods are preserved.",false,(initial:Input,dims:Vec<usize>,viscosity:Input,spacing:Input,dt:Input,steps:usize)->Dynamic,|ctx,initial,dims,viscosity,spacing,dt,steps|{let _=ctx;{let u=vector(&initial.0,&dim([1,-1,0,0])?)?;let h=vector(&spacing.0,&dim([1,0,0,0])?)?;let nu=quantity(&viscosity.0,&dim([2,-1,0,0])?)?;let dt=quantity(&dt.0,&dim([0,1,0,0])?)?;if dims.len()!=2||h.len()!=dims.len()||dims.iter().any(|n|*n<3)||dims.iter().try_fold(1usize,|n,d|n.checked_mul(*d))!=Some(u.len())||h.iter().any(|h|*h<=0.0)||dt<=0.0||nu<0.0||(if steps<=65536{Some(steps+1)}else{None}).and_then(|steps|u.len().checked_mul(steps)).is_none_or(|n|n>LIMIT){return Err(err("Burgers physical shape/stability/storage controls invalid"));}let speed=u.iter().copied().map(f64::abs).fold(0.0,f64::max);if dt*h.iter().map(|h|speed/h+2.0*nu/h/h).sum::<f64>()>0.5{return Err(err("Burgers explicit-step stability bound exceeded"));}let result=cfd::solve_burgers_2d(&u,dims.try_into().map_err(|_|err("dimension mismatch"))?,nu,h.try_into().map_err(|_|err("spacing mismatch"))?,dt,steps).map_err(err)?;let rows=result.len();Ok(Dynamic(array(result.into_iter().flatten().collect(),vec![rows,u.len()],dim([1,-1,0,0])?)?))}});
object_binding!(burgers_3,"num_grid_legacy_burgers_3","crate::sim::physics_cfd::solve_burgers_3d","Measured native continuum adapter. Float matrices and mutable borrows use explicit checked tensor inputs and returned values; SI dimensions, grid shape, work/storage bounds and original native methods are preserved.",false,(initial:Input,dims:Vec<usize>,viscosity:Input,spacing:Input,dt:Input,steps:usize)->Dynamic,|ctx,initial,dims,viscosity,spacing,dt,steps|{let _=ctx;{let u=vector(&initial.0,&dim([1,-1,0,0])?)?;let h=vector(&spacing.0,&dim([1,0,0,0])?)?;let nu=quantity(&viscosity.0,&dim([2,-1,0,0])?)?;let dt=quantity(&dt.0,&dim([0,1,0,0])?)?;if dims.len()!=3||h.len()!=dims.len()||dims.iter().any(|n|*n<3)||dims.iter().try_fold(1usize,|n,d|n.checked_mul(*d))!=Some(u.len())||h.iter().any(|h|*h<=0.0)||dt<=0.0||nu<0.0||(if steps<=65536{Some(steps+1)}else{None}).and_then(|steps|u.len().checked_mul(steps)).is_none_or(|n|n>LIMIT){return Err(err("Burgers physical shape/stability/storage controls invalid"));}let speed=u.iter().copied().map(f64::abs).fold(0.0,f64::max);if dt*h.iter().map(|h|speed/h+2.0*nu/h/h).sum::<f64>()>0.5{return Err(err("Burgers explicit-step stability bound exceeded"));}let result=cfd::solve_burgers_3d(&u,dims.try_into().map_err(|_|err("dimension mismatch"))?,nu,h.try_into().map_err(|_|err("spacing mismatch"))?,dt,steps).map_err(err)?;let rows=result.len();Ok(Dynamic(array(result.into_iter().flatten().collect(),vec![rows,u.len()],dim([1,-1,0,0])?)?))}});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    fdm::install(i)?;
    fvm::install(i)?;
    api::install(
        i,
        &[
            &fluid_new::BINDING,
            &fluid_air::BINDING,
            &fluid_water::BINDING,
            &fluid_kinematic_viscosity::BINDING,
            &fluid_thermal_diffusivity::BINDING,
            &fluid_prandtl::BINDING,
            &compute_vorticity::BINDING,
            &compute_divergence::BINDING,
            &compute_gradient::BINDING,
            &velocity_from_stream_function::BINDING,
            &compute_laplacian::BINDING,
            &compute_stream_function::BINDING,
            &l2_norm::BINDING,
            &max_abs::BINDING,
            &max_velocity_magnitude::BINDING,
            &apply_dirichlet_bc::BINDING,
            &apply_neumann_bc::BINDING,
            &poisson_jacobi::BINDING,
            &poisson_gauss_seidel::BINDING,
            &poisson_sor::BINDING,
            &cavity_simple::BINDING,
            &cavity_mac::BINDING,
            &burgers_1::BINDING,
            &burgers_2::BINDING,
            &burgers_3::BINDING,
        ],
        "numerical_vector",
    )
}
