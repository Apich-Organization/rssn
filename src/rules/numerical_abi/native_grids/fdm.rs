//! Exact low-level FDM kernel adapters; each named fixed scheme calls its original kernel.
use super::CallCtx;
use super::DimExpr;
use super::Dynamic;
use super::Frame;
use super::Handle;
use super::Input;
use super::Installer;
use super::LIMIT;
use super::NativeArg;
use super::NativeType;
use super::NewObject;
use super::OpError;
use super::Output;
use super::RuleError;
use super::TensorValue;
use super::TypeExpr;
use super::Value;
use super::api;
use super::err;
use super::grid_size;
use super::object_binding;
use super::payload;
use crate::sim::fdm_nd as fd;
#[derive(Clone)]
enum NativeGrid {
    D1(fd::Grid1D),
    D2(fd::Grid),
    D3(fd::Grid3D),
}
#[derive(Clone)]
struct FdGrid {
    data: NativeGrid,
    frame: Frame,
}
impl FdGrid {
    fn extents(&self) -> Vec<usize> {
        match self.data {
            | NativeGrid::D1(g) => vec![g.nx],
            | NativeGrid::D2(g) => vec![g.nx, g.ny],
            | NativeGrid::D3(g) => vec![g.nx, g.ny, g.nz],
        }
    }

    fn len(&self) -> usize {
        self.extents().iter().product()
    }
}
impl NativeType for FdGrid {
    const KIND: &'static str = "rssn/native/fdm-grid/v1";
    const MAX_BYTES: usize = 4096;

    fn validate(&self) -> Result<(), OpError> {
        self.frame.validate()?;
        let ext = self.extents();
        if ext.iter().any(|n| *n < 3 || *n > 4096)
            || ext
                .iter()
                .try_fold(1usize, |n, m| n.checked_mul(*m))
                .is_none_or(|n| n > LIMIT)
        {
            return Err(err("FDM grid extent/allocation budget invalid"));
        }
        let h = match self.data {
            | NativeGrid::D1(g) => g.h,
            | NativeGrid::D2(g) => g.h,
            | NativeGrid::D3(g) => g.h,
        };
        if !h.is_finite() || h <= 0.0 {
            return Err(err("FDM spacing must be positive finite"));
        }
        Ok(())
    }
}
#[derive(Clone)]
pub(super) struct Condition {
    pub(super) data: fd::BoundaryCondition,
    pub(super) frame: Frame,
}
impl NativeType for Condition {
    const KIND: &'static str = "rssn/native/fdm-boundary/v1";
    const MAX_BYTES: usize = 65536;

    fn validate(&self) -> Result<(), OpError> {
        self.frame.validate()?;
        if matches!(&self.data,fd::BoundaryCondition::DirichletProfile(v)if v.len()>4096||v.iter().any(|v|!v.is_finite()))
        {
            return Err(err("boundary profile size/finite budget invalid"));
        }
        Ok(())
    }
}
#[derive(Clone)]
pub(super) enum NativeBcs {
    D1(fd::Boundary1D),
    D2(fd::Boundary2D),
    D3(fd::Boundary3D),
}
#[derive(Clone)]
pub(super) struct Bcs {
    pub(super) data: NativeBcs,
    pub(super) frame: Frame,
}
impl NativeType for Bcs {
    const KIND: &'static str = "rssn/native/fdm-face-bcs/v1";
    const MAX_BYTES: usize = 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        self.frame.validate()
    }
}
#[derive(Clone)]
enum NativeField {
    D1(fd::Field1D),
    D2(fd::Field),
    D3(fd::Field3D),
}
impl NativeField {
    fn values(&self) -> &[f64] {
        match self {
            | Self::D1(f) => &f.values,
            | Self::D2(f) => &f.values,
            | Self::D3(f) => &f.values,
        }
    }
}
#[derive(Clone)]
struct Field {
    data: NativeField,
    frame: Frame,
}
impl NativeType for Field {
    const KIND: &'static str = "rssn/native/fdm-field/v1";
    const MAX_BYTES: usize = LIMIT * 8 + 4096;

    fn validate(&self) -> Result<(), OpError> {
        self.frame.validate()?;
        if self.data.values().len() > LIMIT || self.data.values().iter().any(|v| !v.is_finite()) {
            return Err(err("FDM field finite/allocation budget invalid"));
        }
        Ok(())
    }
}
struct History {
    data: Vec<NativeField>,
    frame: Frame,
    dt: f64,
}
impl NativeType for History {
    const KIND: &'static str = "rssn/native/fdm-history/v1";
    const MAX_BYTES: usize = 8 * 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        self.frame.validate()?;
        if self.data.len() > 8193
            || self
                .data
                .iter()
                .try_fold(0usize, |n, f| n.checked_add(f.values().len()))
                .is_none_or(|n| n > LIMIT)
            || self
                .data
                .iter()
                .flat_map(NativeField::values)
                .any(|v| !v.is_finite())
        {
            return Err(err("FDM history finite/allocation budget invalid"));
        }
        Ok(())
    }
}
fn stepping(name: &str) -> Result<fd::TimeStepping, OpError> {
    Ok(match name {
        | "Explicit" => fd::TimeStepping::Explicit,
        | "BackwardEuler" => fd::TimeStepping::BackwardEuler,
        | "CrankNicolson" => fd::TimeStepping::CrankNicolson,
        | "Adi" => fd::TimeStepping::Adi,
        | _ => return Err(err("unknown native FDM time stepping")),
    })
}
fn field_result(
    data: NativeField,
    f: &Frame,
) -> NewObject<Field> {
    NewObject::new(
        Field {
            data,
            frame: f.clone(),
        },
        f.dims.to_vec(),
    )
}
fn history(
    data: Vec<NativeField>,
    f: &Frame,
    dt: f64,
) -> NewObject<History> {
    NewObject::new(
        History {
            data,
            frame: f.clone(),
            dt,
        },
        f.dims.to_vec(),
    )
}
object_binding!(fdm_fdm_grid,"num_grid_fdm_grid","","Original native low-level FDM method with explicit Length/Time/field references. Coordinates, steps and coefficients are normalized; output fields/histories retain their physical field reference. Native x-fastest layout, selected integrator and stability/domain checks are preserved with allocation/history budgets.",true,(frame:Handle<Frame>,extents:Vec<usize>,spacing:f64)->NewObject<FdGrid>,|ctx,frame,extents,spacing|{let _=ctx;grid_size(&extents)?;let data=match extents.as_slice(){[x]=>NativeGrid::D1(fd::Grid1D{nx:*x,h:spacing}),[x,y]=>NativeGrid::D2(fd::Grid{nx:*x,ny:*y,h:spacing}),[x,y,z]=>NativeGrid::D3(fd::Grid3D{nx:*x,ny:*y,nz:*z,h:spacing}),_=>return Err(err("FDM grid requires one through three axes"))};frame.with(ctx,|f|Ok(NewObject::new(FdGrid{data,frame:f.clone()},f.dims.to_vec())))});
object_binding!(fdm_fdm_zero_dirichlet,"num_grid_fdm_zero_dirichlet","crate::sim::fdm_nd::BoundaryCondition::zero_dirichlet","Original native low-level FDM method with explicit Length/Time/field references. Coordinates, steps and coefficients are normalized; output fields/histories retain their physical field reference. Native x-fastest layout, selected integrator and stability/domain checks are preserved with allocation/history budgets.",true,(frame:Handle<Frame>)->NewObject<Condition>,|ctx,frame|{let _=ctx;frame.with(ctx,|f|Ok(NewObject::new(Condition{data:fd::BoundaryCondition::zero_dirichlet(),frame:f.clone()},f.dims.to_vec())))});
object_binding!(fdm_fdm_zero_neumann,"num_grid_fdm_zero_neumann","crate::sim::fdm_nd::BoundaryCondition::zero_neumann","Original native low-level FDM method with explicit Length/Time/field references. Coordinates, steps and coefficients are normalized; output fields/histories retain their physical field reference. Native x-fastest layout, selected integrator and stability/domain checks are preserved with allocation/history budgets.",true,(frame:Handle<Frame>)->NewObject<Condition>,|ctx,frame|{let _=ctx;frame.with(ctx,|f|Ok(NewObject::new(Condition{data:fd::BoundaryCondition::zero_neumann(),frame:f.clone()},f.dims.to_vec())))});
object_binding!(fdm_fdm_condition,"num_grid_fdm_condition","","Original native low-level FDM method with explicit Length/Time/field references. Coordinates, steps and coefficients are normalized; output fields/histories retain their physical field reference. Native x-fastest layout, selected integrator and stability/domain checks are preserved with allocation/history budgets.",true,(frame:Handle<Frame>,kind:String,value:f64,alpha:f64,beta:f64,profile:Vec<f64>)->NewObject<Condition>,|ctx,frame,kind,value,alpha,beta,profile|{let _=ctx;frame.with(ctx,|f|{let data=match kind.as_str(){"Dirichlet"=>fd::BoundaryCondition::Dirichlet(value),"Neumann"=>fd::BoundaryCondition::Neumann(value),"Robin"=>{if crate::kernels::exact_eq(alpha,0.0)&&crate::kernels::exact_eq(beta,0.0){return Err(err("degenerate Robin condition"));}fd::BoundaryCondition::Robin{alpha,beta,gamma:value}},"Periodic"=>fd::BoundaryCondition::Periodic,"DirichletProfile"=>fd::BoundaryCondition::DirichletProfile(profile),_=>return Err(err("unknown FDM scalar boundary kind"))};Ok(NewObject::new(Condition{data,frame:f.clone()},f.dims.to_vec()))})});
object_binding!(fdm_fdm1_boundary_uniform,"num_grid_fdm1_boundary_uniform","crate::sim::fdm_nd::Boundary1D::uniform","Original native low-level FDM method with explicit Length/Time/field references. Coordinates, steps and coefficients are normalized; output fields/histories retain their physical field reference. Native x-fastest layout, selected integrator and stability/domain checks are preserved with allocation/history budgets.",true,(condition:Handle<Condition>)->NewObject<Bcs>,|ctx,condition|{let _=ctx;condition.with(ctx,|c|Ok(NewObject::new(Bcs{data:NativeBcs::D1(fd::Boundary1D::uniform(c.data.clone())),frame:c.frame.clone()},c.frame.dims.to_vec())))});
object_binding!(fdm_fdm1_boundary_periodic,"num_grid_fdm1_boundary_periodic","crate::sim::fdm_nd::Boundary1D::periodic","Original native low-level FDM method with explicit Length/Time/field references. Coordinates, steps and coefficients are normalized; output fields/histories retain their physical field reference. Native x-fastest layout, selected integrator and stability/domain checks are preserved with allocation/history budgets.",true,(frame:Handle<Frame>)->NewObject<Bcs>,|ctx,frame|{let _=ctx;frame.with(ctx,|f|Ok(NewObject::new(Bcs{data:NativeBcs::D1(fd::Boundary1D::periodic()),frame:f.clone()},f.dims.to_vec())))});
object_binding!(fdm_fdm1_len,"num_grid_fdm1_len","crate::sim::fdm_nd::Grid1D::len","Original native low-level FDM method with explicit Length/Time/field references. Coordinates, steps and coefficients are normalized; output fields/histories retain their physical field reference. Native x-fastest layout, selected integrator and stability/domain checks are preserved with allocation/history budgets.",true,(grid:Handle<FdGrid>)->usize,|ctx,grid|{let _=ctx;grid.with(ctx,|g|{let NativeGrid::D1(n)=g.data else{return Err(err("FDM grid dimension mismatch"));};Ok(n.len())})});
object_binding!(fdm_fdm1_is_empty,"num_grid_fdm1_is_empty","crate::sim::fdm_nd::Grid1D::is_empty","Original native low-level FDM method with explicit Length/Time/field references. Coordinates, steps and coefficients are normalized; output fields/histories retain their physical field reference. Native x-fastest layout, selected integrator and stability/domain checks are preserved with allocation/history budgets.",true,(grid:Handle<FdGrid>)->bool,|ctx,grid|{let _=ctx;grid.with(ctx,|g|{let NativeGrid::D1(n)=g.data else{return Err(err("FDM grid dimension mismatch"));};Ok(n.is_empty())})});
object_binding!(fdm_fdm2_boundary_uniform,"num_grid_fdm2_boundary_uniform","crate::sim::fdm_nd::Boundary2D::uniform","Original native low-level FDM method with explicit Length/Time/field references. Coordinates, steps and coefficients are normalized; output fields/histories retain their physical field reference. Native x-fastest layout, selected integrator and stability/domain checks are preserved with allocation/history budgets.",true,(condition:Handle<Condition>)->NewObject<Bcs>,|ctx,condition|{let _=ctx;condition.with(ctx,|c|Ok(NewObject::new(Bcs{data:NativeBcs::D2(fd::Boundary2D::uniform(c.data.clone())),frame:c.frame.clone()},c.frame.dims.to_vec())))});
object_binding!(fdm_fdm2_boundary_periodic,"num_grid_fdm2_boundary_periodic","crate::sim::fdm_nd::Boundary2D::periodic","Original native low-level FDM method with explicit Length/Time/field references. Coordinates, steps and coefficients are normalized; output fields/histories retain their physical field reference. Native x-fastest layout, selected integrator and stability/domain checks are preserved with allocation/history budgets.",true,(frame:Handle<Frame>)->NewObject<Bcs>,|ctx,frame|{let _=ctx;frame.with(ctx,|f|Ok(NewObject::new(Bcs{data:NativeBcs::D2(fd::Boundary2D::periodic()),frame:f.clone()},f.dims.to_vec())))});
object_binding!(fdm_fdm2_len,"num_grid_fdm2_len","crate::sim::fdm_nd::Grid::len","Original native low-level FDM method with explicit Length/Time/field references. Coordinates, steps and coefficients are normalized; output fields/histories retain their physical field reference. Native x-fastest layout, selected integrator and stability/domain checks are preserved with allocation/history budgets.",true,(grid:Handle<FdGrid>)->usize,|ctx,grid|{let _=ctx;grid.with(ctx,|g|{let NativeGrid::D2(n)=g.data else{return Err(err("FDM grid dimension mismatch"));};Ok(n.len())})});
object_binding!(fdm_fdm2_is_empty,"num_grid_fdm2_is_empty","crate::sim::fdm_nd::Grid::is_empty","Original native low-level FDM method with explicit Length/Time/field references. Coordinates, steps and coefficients are normalized; output fields/histories retain their physical field reference. Native x-fastest layout, selected integrator and stability/domain checks are preserved with allocation/history budgets.",true,(grid:Handle<FdGrid>)->bool,|ctx,grid|{let _=ctx;grid.with(ctx,|g|{let NativeGrid::D2(n)=g.data else{return Err(err("FDM grid dimension mismatch"));};Ok(n.is_empty())})});
object_binding!(fdm_fdm3_boundary_uniform,"num_grid_fdm3_boundary_uniform","crate::sim::fdm_nd::Boundary3D::uniform","Original native low-level FDM method with explicit Length/Time/field references. Coordinates, steps and coefficients are normalized; output fields/histories retain their physical field reference. Native x-fastest layout, selected integrator and stability/domain checks are preserved with allocation/history budgets.",true,(condition:Handle<Condition>)->NewObject<Bcs>,|ctx,condition|{let _=ctx;condition.with(ctx,|c|Ok(NewObject::new(Bcs{data:NativeBcs::D3(fd::Boundary3D::uniform(c.data.clone())),frame:c.frame.clone()},c.frame.dims.to_vec())))});
object_binding!(fdm_fdm3_boundary_periodic,"num_grid_fdm3_boundary_periodic","crate::sim::fdm_nd::Boundary3D::periodic","Original native low-level FDM method with explicit Length/Time/field references. Coordinates, steps and coefficients are normalized; output fields/histories retain their physical field reference. Native x-fastest layout, selected integrator and stability/domain checks are preserved with allocation/history budgets.",true,(frame:Handle<Frame>)->NewObject<Bcs>,|ctx,frame|{let _=ctx;frame.with(ctx,|f|Ok(NewObject::new(Bcs{data:NativeBcs::D3(fd::Boundary3D::periodic()),frame:f.clone()},f.dims.to_vec())))});
object_binding!(fdm_fdm3_len,"num_grid_fdm3_len","crate::sim::fdm_nd::Grid3D::len","Original native low-level FDM method with explicit Length/Time/field references. Coordinates, steps and coefficients are normalized; output fields/histories retain their physical field reference. Native x-fastest layout, selected integrator and stability/domain checks are preserved with allocation/history budgets.",true,(grid:Handle<FdGrid>)->usize,|ctx,grid|{let _=ctx;grid.with(ctx,|g|{let NativeGrid::D3(n)=g.data else{return Err(err("FDM grid dimension mismatch"));};Ok(n.len())})});
object_binding!(fdm_fdm3_is_empty,"num_grid_fdm3_is_empty","crate::sim::fdm_nd::Grid3D::is_empty","Original native low-level FDM method with explicit Length/Time/field references. Coordinates, steps and coefficients are normalized; output fields/histories retain their physical field reference. Native x-fastest layout, selected integrator and stability/domain checks are preserved with allocation/history budgets.",true,(grid:Handle<FdGrid>)->bool,|ctx,grid|{let _=ctx;grid.with(ctx,|g|{let NativeGrid::D3(n)=g.data else{return Err(err("FDM grid dimension mismatch"));};Ok(n.is_empty())})});
object_binding!(fdm_fdm3_index,"num_grid_fdm3_index","crate::sim::fdm_nd::Grid3D::index","Original native low-level FDM method with explicit Length/Time/field references. Coordinates, steps and coefficients are normalized; output fields/histories retain their physical field reference. Native x-fastest layout, selected integrator and stability/domain checks are preserved with allocation/history budgets.",true,(grid:Handle<FdGrid>,x:usize,y:usize,z:usize)->usize,|ctx,grid,x,y,z|{let _=ctx;grid.with(ctx,|g|{let NativeGrid::D3(n)=g.data else{return Err(err("expected three-dimensional FDM grid"));};
if x>=n.nx||y>=n.ny||z>=n.nz{return Err(err("FDM node index out of range"));}Ok(n.index(x,y,z))})});
object_binding!(fdm_fdm3_field_new,"num_grid_fdm3_field_new","","Original native low-level FDM method with explicit Length/Time/field references. Coordinates, steps and coefficients are normalized; output fields/histories retain their physical field reference. Native x-fastest layout, selected integrator and stability/domain checks are preserved with allocation/history budgets.",true,(grid:Handle<FdGrid>,values:Vec<f64>)->NewObject<Field>,|ctx,grid,values|{let _=ctx;grid.with(ctx,|g|{let NativeGrid::D3(n)=g.data else{return Err(err("expected three-dimensional FDM grid"));};
if values.len()!=n.len(){return Err(err("FDM field shape mismatch"));}Ok(field_result(NativeField::D3(fd::Field3D{grid:n,values}),&g.frame))})});
object_binding!(fdm_fdm3_field_at,"num_grid_fdm3_field_at","crate::sim::fdm_nd::Field3D::at","Original native low-level FDM method with explicit Length/Time/field references. Coordinates, steps and coefficients are normalized; output fields/histories retain their physical field reference. Native x-fastest layout, selected integrator and stability/domain checks are preserved with allocation/history budgets.",true,(field:Handle<Field>,x:usize,y:usize,z:usize)->Dynamic,|ctx,field,x,y,z|{let _=ctx;field.with(ctx,|f|{let NativeField::D3(n)=&f.data else{return Err(err("expected three-dimensional FDM field"));};
if x>=n.grid.nx||y>=n.grid.ny||z>=n.grid.nz{return Err(err("FDM field index out of range"));}Ok(Dynamic(api::quantity(n.at(x,y,z)*f.frame.value,f.frame.dims[2].clone())?))})});
object_binding!(fdm_fdm_solve_poisson_1d,"num_grid_fdm_solve_poisson_1d","crate::sim::fdm_nd::solve_poisson_1d","Original native solve_poisson_1d: only its actual controls are exposed. The plain 2D/3D zero-boundary entry points preserve zero boundary semantics; *_with_bc variants use explicit boundary objects. All numeric data use the retained Length/Time/field references, native x-fastest layout and selected method; work/history bounds are enforced.",true,(grid:Handle<FdGrid>,bcs:Handle<Bcs>,source:Vec<f64>)->NewObject<Field>,|ctx,grid,bcs,source|{grid.with(ctx,|g|bcs.with(ctx,|b|{let NativeGrid::D1(n)=&g.data else{return Err(err("FDM grid dimension mismatch"));};g.frame.same(&b.frame)?;let NativeBcs::D1(bc)=&b.data else{return Err(err("FDM boundary dimension mismatch"));};profiles(g,b)?;if source.len()!=g.len(){return Err(err("Poisson source length mismatch"));}let out=fd::solve_poisson_1d(*n,&source,bc).map_err(err)?;Ok(field_result(NativeField::D1(out),&g.frame))}))});
object_binding!(fdm_fdm_solve_heat_1d,"num_grid_fdm_solve_heat_1d","crate::sim::fdm_nd::solve_heat_1d","Original native solve_heat_1d: only its actual controls are exposed. The plain 2D/3D zero-boundary entry points preserve zero boundary semantics; *_with_bc variants use explicit boundary objects. All numeric data use the retained Length/Time/field references, native x-fastest layout and selected method; work/history bounds are enforced.",true,(grid:Handle<FdGrid>,bcs:Handle<Bcs>,initial:Vec<f64>,dt:f64,steps:usize,diffusivity:f64,method:String)->NewObject<History>,|ctx,grid,bcs,initial,dt,steps,diffusivity,method|{grid.with(ctx,|g|bcs.with(ctx,|b|{let NativeGrid::D1(n)=&g.data else{return Err(err("FDM grid dimension mismatch"));};g.frame.same(&b.frame)?;let NativeBcs::D1(bc)=&b.data else{return Err(err("FDM boundary dimension mismatch"));};profiles(g,b)?;if initial.len()!=g.len()||dt<=0.0||steps>8192||g.len().checked_mul(steps+1).is_none_or(|n|n>LIMIT){return Err(err("FDM initial/time/retained-history budget invalid"));}
if diffusivity<0.0{return Err(err("diffusivity must be nonnegative"));}let out=fd::solve_heat_1d(*n,diffusivity,dt,steps,&initial,bc,stepping(&method)?).map_err(err)?;Ok(history(out.into_iter().map(NativeField::D1).collect(),&g.frame,dt))}))});
object_binding!(fdm_fdm_solve_wave_1d,"num_grid_fdm_solve_wave_1d","crate::sim::fdm_nd::solve_wave_1d","Original native solve_wave_1d: only its actual controls are exposed. The plain 2D/3D zero-boundary entry points preserve zero boundary semantics; *_with_bc variants use explicit boundary objects. All numeric data use the retained Length/Time/field references, native x-fastest layout and selected method; work/history bounds are enforced.",true,(grid:Handle<FdGrid>,bcs:Handle<Bcs>,initial:Vec<f64>,dt:f64,steps:usize,speed:f64,initial_velocity:Vec<f64>)->NewObject<History>,|ctx,grid,bcs,initial,dt,steps,speed,initial_velocity|{grid.with(ctx,|g|bcs.with(ctx,|b|{let NativeGrid::D1(n)=&g.data else{return Err(err("FDM grid dimension mismatch"));};g.frame.same(&b.frame)?;let NativeBcs::D1(bc)=&b.data else{return Err(err("FDM boundary dimension mismatch"));};profiles(g,b)?;if initial.len()!=g.len()||dt<=0.0||steps>8192||g.len().checked_mul(steps+1).is_none_or(|n|n>LIMIT){return Err(err("FDM initial/time/retained-history budget invalid"));}
if speed<0.0||initial_velocity.len()!=initial.len(){return Err(err("wave speed/initial velocity invalid"));}let out=fd::solve_wave_1d(*n,speed,dt,steps,&initial,&initial_velocity,bc).map_err(err)?;Ok(history(out.into_iter().map(NativeField::D1).collect(),&g.frame,dt))}))});
object_binding!(fdm_fdm_solve_advection_1d,"num_grid_fdm_solve_advection_1d","crate::sim::fdm_nd::solve_advection_1d","Original native solve_advection_1d: only its actual controls are exposed. The plain 2D/3D zero-boundary entry points preserve zero boundary semantics; *_with_bc variants use explicit boundary objects. All numeric data use the retained Length/Time/field references, native x-fastest layout and selected method; work/history bounds are enforced.",true,(grid:Handle<FdGrid>,bcs:Handle<Bcs>,initial:Vec<f64>,dt:f64,steps:usize,velocity:Vec<f64>)->NewObject<History>,|ctx,grid,bcs,initial,dt,steps,velocity|{grid.with(ctx,|g|bcs.with(ctx,|b|{let NativeGrid::D1(n)=&g.data else{return Err(err("FDM grid dimension mismatch"));};g.frame.same(&b.frame)?;let NativeBcs::D1(bc)=&b.data else{return Err(err("FDM boundary dimension mismatch"));};profiles(g,b)?;if initial.len()!=g.len()||dt<=0.0||steps>8192||g.len().checked_mul(steps+1).is_none_or(|n|n>LIMIT){return Err(err("FDM initial/time/retained-history budget invalid"));}
if velocity.len()!=1{return Err(err("advection velocity dimension mismatch"));}let out=fd::solve_advection_1d(*n,velocity[0],dt,steps,&initial,bc).map_err(err)?;Ok(history(out.into_iter().map(NativeField::D1).collect(),&g.frame,dt))}))});
object_binding!(fdm_fdm_solve_poisson,"num_grid_fdm_solve_poisson","crate::sim::fdm_nd::solve_poisson","Original native solve_poisson: only its actual controls are exposed. The plain 2D/3D zero-boundary entry points preserve zero boundary semantics; *_with_bc variants use explicit boundary objects. All numeric data use the retained Length/Time/field references, native x-fastest layout and selected method; work/history bounds are enforced.",true,(grid:Handle<FdGrid>,source:Vec<f64>,iterations:usize,tolerance:f64)->NewObject<Field>,|ctx,grid,source,iterations,tolerance|{grid.with(ctx,|g|{let NativeGrid::D2(n)=&g.data else{return Err(err("FDM grid dimension mismatch"));};
if source.len()!=g.len(){return Err(err("Poisson source length mismatch"));}
if iterations==0||iterations>8192||tolerance<=0.0||g.len().checked_mul(iterations).is_none_or(|w|w>8_388_608){return Err(err("Poisson iteration/tolerance/work budget exceeded"));}let out=fd::solve_poisson(&fd::PoissonSpec{grid:*n,source,iterations,tolerance}).map_err(err)?;Ok(field_result(NativeField::D2(out),&g.frame))})});
object_binding!(fdm_fdm_solve_poisson_with_bc,"num_grid_fdm_solve_poisson_with_bc","crate::sim::fdm_nd::solve_poisson_with_bc","Original native solve_poisson_with_bc: only its actual controls are exposed. The plain 2D/3D zero-boundary entry points preserve zero boundary semantics; *_with_bc variants use explicit boundary objects. All numeric data use the retained Length/Time/field references, native x-fastest layout and selected method; work/history bounds are enforced.",true,(grid:Handle<FdGrid>,bcs:Handle<Bcs>,source:Vec<f64>,iterations:usize,tolerance:f64)->NewObject<Field>,|ctx,grid,bcs,source,iterations,tolerance|{grid.with(ctx,|g|bcs.with(ctx,|b|{let NativeGrid::D2(n)=&g.data else{return Err(err("FDM grid dimension mismatch"));};g.frame.same(&b.frame)?;let NativeBcs::D2(bc)=&b.data else{return Err(err("FDM boundary dimension mismatch"));};profiles(g,b)?;if source.len()!=g.len(){return Err(err("Poisson source length mismatch"));}
if iterations==0||iterations>8192||tolerance<=0.0||g.len().checked_mul(iterations).is_none_or(|w|w>8_388_608){return Err(err("Poisson iteration/tolerance/work budget exceeded"));}let out=fd::solve_poisson_with_bc(&fd::PoissonSpec{grid:*n,source,iterations,tolerance},bc).map_err(err)?;Ok(field_result(NativeField::D2(out),&g.frame))}))});
object_binding!(fdm_fdm_solve_heat,"num_grid_fdm_solve_heat","crate::sim::fdm_nd::solve_heat","Original native solve_heat: only its actual controls are exposed. The plain 2D/3D zero-boundary entry points preserve zero boundary semantics; *_with_bc variants use explicit boundary objects. All numeric data use the retained Length/Time/field references, native x-fastest layout and selected method; work/history bounds are enforced.",true,(grid:Handle<FdGrid>,initial:Vec<f64>,dt:f64,steps:usize,diffusivity:f64)->NewObject<History>,|ctx,grid,initial,dt,steps,diffusivity|{grid.with(ctx,|g|{let NativeGrid::D2(n)=&g.data else{return Err(err("FDM grid dimension mismatch"));};
if initial.len()!=g.len()||dt<=0.0||steps>8192||g.len().checked_mul(steps+1).is_none_or(|n|n>LIMIT){return Err(err("FDM initial/time/retained-history budget invalid"));}
if diffusivity<0.0{return Err(err("diffusivity must be nonnegative"));}let out=fd::solve_heat(&fd::HeatSpec{grid:*n,alpha:diffusivity,steps,dt,initial}).map_err(err)?;Ok(history(out.into_iter().map(NativeField::D2).collect(),&g.frame,dt))})});
object_binding!(fdm_fdm_solve_heat_with_bc,"num_grid_fdm_solve_heat_with_bc","crate::sim::fdm_nd::solve_heat_with_bc","Original native solve_heat_with_bc: only its actual controls are exposed. The plain 2D/3D zero-boundary entry points preserve zero boundary semantics; *_with_bc variants use explicit boundary objects. All numeric data use the retained Length/Time/field references, native x-fastest layout and selected method; work/history bounds are enforced.",true,(grid:Handle<FdGrid>,bcs:Handle<Bcs>,initial:Vec<f64>,dt:f64,steps:usize,diffusivity:f64,method:String)->NewObject<History>,|ctx,grid,bcs,initial,dt,steps,diffusivity,method|{grid.with(ctx,|g|bcs.with(ctx,|b|{let NativeGrid::D2(n)=&g.data else{return Err(err("FDM grid dimension mismatch"));};g.frame.same(&b.frame)?;let NativeBcs::D2(bc)=&b.data else{return Err(err("FDM boundary dimension mismatch"));};profiles(g,b)?;if initial.len()!=g.len()||dt<=0.0||steps>8192||g.len().checked_mul(steps+1).is_none_or(|n|n>LIMIT){return Err(err("FDM initial/time/retained-history budget invalid"));}
if diffusivity<0.0{return Err(err("diffusivity must be nonnegative"));}let out=fd::solve_heat_with_bc(&fd::HeatSpec{grid:*n,alpha:diffusivity,steps,dt,initial},bc,stepping(&method)?).map_err(err)?;Ok(history(out.into_iter().map(NativeField::D2).collect(),&g.frame,dt))}))});
object_binding!(fdm_fdm_solve_wave,"num_grid_fdm_solve_wave","crate::sim::fdm_nd::solve_wave","Original native solve_wave: only its actual controls are exposed. The plain 2D/3D zero-boundary entry points preserve zero boundary semantics; *_with_bc variants use explicit boundary objects. All numeric data use the retained Length/Time/field references, native x-fastest layout and selected method; work/history bounds are enforced.",true,(grid:Handle<FdGrid>,initial:Vec<f64>,dt:f64,steps:usize,speed:f64,initial_velocity:Vec<f64>)->NewObject<History>,|ctx,grid,initial,dt,steps,speed,initial_velocity|{grid.with(ctx,|g|{let NativeGrid::D2(n)=&g.data else{return Err(err("FDM grid dimension mismatch"));};
if initial.len()!=g.len()||dt<=0.0||steps>8192||g.len().checked_mul(steps+1).is_none_or(|n|n>LIMIT){return Err(err("FDM initial/time/retained-history budget invalid"));}
if speed<0.0||initial_velocity.len()!=initial.len(){return Err(err("wave speed/initial velocity invalid"));}let out=fd::solve_wave(&fd::WaveSpec{grid:*n,c:speed,steps,dt,initial,velocity:initial_velocity}).map_err(err)?;Ok(history(out.into_iter().map(NativeField::D2).collect(),&g.frame,dt))})});
object_binding!(fdm_fdm_solve_wave_with_bc,"num_grid_fdm_solve_wave_with_bc","crate::sim::fdm_nd::solve_wave_with_bc","Original native solve_wave_with_bc: only its actual controls are exposed. The plain 2D/3D zero-boundary entry points preserve zero boundary semantics; *_with_bc variants use explicit boundary objects. All numeric data use the retained Length/Time/field references, native x-fastest layout and selected method; work/history bounds are enforced.",true,(grid:Handle<FdGrid>,bcs:Handle<Bcs>,initial:Vec<f64>,dt:f64,steps:usize,speed:f64,initial_velocity:Vec<f64>)->NewObject<History>,|ctx,grid,bcs,initial,dt,steps,speed,initial_velocity|{grid.with(ctx,|g|bcs.with(ctx,|b|{let NativeGrid::D2(n)=&g.data else{return Err(err("FDM grid dimension mismatch"));};g.frame.same(&b.frame)?;let NativeBcs::D2(bc)=&b.data else{return Err(err("FDM boundary dimension mismatch"));};profiles(g,b)?;if initial.len()!=g.len()||dt<=0.0||steps>8192||g.len().checked_mul(steps+1).is_none_or(|n|n>LIMIT){return Err(err("FDM initial/time/retained-history budget invalid"));}
if speed<0.0||initial_velocity.len()!=initial.len(){return Err(err("wave speed/initial velocity invalid"));}let out=fd::solve_wave_with_bc(&fd::WaveSpec{grid:*n,c:speed,steps,dt,initial,velocity:initial_velocity},bc).map_err(err)?;Ok(history(out.into_iter().map(NativeField::D2).collect(),&g.frame,dt))}))});
object_binding!(fdm_fdm_solve_advection,"num_grid_fdm_solve_advection","crate::sim::fdm_nd::solve_advection","Original native solve_advection: only its actual controls are exposed. The plain 2D/3D zero-boundary entry points preserve zero boundary semantics; *_with_bc variants use explicit boundary objects. All numeric data use the retained Length/Time/field references, native x-fastest layout and selected method; work/history bounds are enforced.",true,(grid:Handle<FdGrid>,initial:Vec<f64>,dt:f64,steps:usize,velocity:Vec<f64>)->NewObject<History>,|ctx,grid,initial,dt,steps,velocity|{grid.with(ctx,|g|{let NativeGrid::D2(n)=&g.data else{return Err(err("FDM grid dimension mismatch"));};
if initial.len()!=g.len()||dt<=0.0||steps>8192||g.len().checked_mul(steps+1).is_none_or(|n|n>LIMIT){return Err(err("FDM initial/time/retained-history budget invalid"));}
if velocity.len()!=2{return Err(err("advection velocity dimension mismatch"));}let out=fd::solve_advection(&fd::AdvectionSpec{grid:*n,vx:velocity[0],vy:velocity[1],steps,dt,initial}).map_err(err)?;Ok(history(out.into_iter().map(NativeField::D2).collect(),&g.frame,dt))})});
object_binding!(fdm_fdm_solve_advection_with_bc,"num_grid_fdm_solve_advection_with_bc","crate::sim::fdm_nd::solve_advection_with_bc","Original native solve_advection_with_bc: only its actual controls are exposed. The plain 2D/3D zero-boundary entry points preserve zero boundary semantics; *_with_bc variants use explicit boundary objects. All numeric data use the retained Length/Time/field references, native x-fastest layout and selected method; work/history bounds are enforced.",true,(grid:Handle<FdGrid>,bcs:Handle<Bcs>,initial:Vec<f64>,dt:f64,steps:usize,velocity:Vec<f64>)->NewObject<History>,|ctx,grid,bcs,initial,dt,steps,velocity|{grid.with(ctx,|g|bcs.with(ctx,|b|{let NativeGrid::D2(n)=&g.data else{return Err(err("FDM grid dimension mismatch"));};g.frame.same(&b.frame)?;let NativeBcs::D2(bc)=&b.data else{return Err(err("FDM boundary dimension mismatch"));};profiles(g,b)?;if initial.len()!=g.len()||dt<=0.0||steps>8192||g.len().checked_mul(steps+1).is_none_or(|n|n>LIMIT){return Err(err("FDM initial/time/retained-history budget invalid"));}
if velocity.len()!=2{return Err(err("advection velocity dimension mismatch"));}let out=fd::solve_advection_with_bc(&fd::AdvectionSpec{grid:*n,vx:velocity[0],vy:velocity[1],steps,dt,initial},bc).map_err(err)?;Ok(history(out.into_iter().map(NativeField::D2).collect(),&g.frame,dt))}))});
object_binding!(fdm_fdm_solve_poisson_3d,"num_grid_fdm_solve_poisson_3d","crate::sim::fdm_nd::solve_poisson_3d","Original native solve_poisson_3d: only its actual controls are exposed. The plain 2D/3D zero-boundary entry points preserve zero boundary semantics; *_with_bc variants use explicit boundary objects. All numeric data use the retained Length/Time/field references, native x-fastest layout and selected method; work/history bounds are enforced.",true,(grid:Handle<FdGrid>,bcs:Handle<Bcs>,source:Vec<f64>,iterations:usize,tolerance:f64)->NewObject<Field>,|ctx,grid,bcs,source,iterations,tolerance|{grid.with(ctx,|g|bcs.with(ctx,|b|{let NativeGrid::D3(n)=&g.data else{return Err(err("FDM grid dimension mismatch"));};g.frame.same(&b.frame)?;let NativeBcs::D3(bc)=&b.data else{return Err(err("FDM boundary dimension mismatch"));};profiles(g,b)?;if source.len()!=g.len(){return Err(err("Poisson source length mismatch"));}
if iterations==0||iterations>8192||tolerance<=0.0||g.len().checked_mul(iterations).is_none_or(|w|w>8_388_608){return Err(err("Poisson iteration/tolerance/work budget exceeded"));}let out=fd::solve_poisson_3d(*n,&source,bc,iterations,tolerance).map_err(err)?;Ok(field_result(NativeField::D3(out),&g.frame))}))});
object_binding!(fdm_fdm_solve_heat_3d,"num_grid_fdm_solve_heat_3d","crate::sim::fdm_nd::solve_heat_3d","Original native solve_heat_3d: only its actual controls are exposed. The plain 2D/3D zero-boundary entry points preserve zero boundary semantics; *_with_bc variants use explicit boundary objects. All numeric data use the retained Length/Time/field references, native x-fastest layout and selected method; work/history bounds are enforced.",true,(grid:Handle<FdGrid>,bcs:Handle<Bcs>,initial:Vec<f64>,dt:f64,steps:usize,diffusivity:f64,method:String)->NewObject<History>,|ctx,grid,bcs,initial,dt,steps,diffusivity,method|{grid.with(ctx,|g|bcs.with(ctx,|b|{let NativeGrid::D3(n)=&g.data else{return Err(err("FDM grid dimension mismatch"));};g.frame.same(&b.frame)?;let NativeBcs::D3(bc)=&b.data else{return Err(err("FDM boundary dimension mismatch"));};profiles(g,b)?;if initial.len()!=g.len()||dt<=0.0||steps>8192||g.len().checked_mul(steps+1).is_none_or(|n|n>LIMIT){return Err(err("FDM initial/time/retained-history budget invalid"));}
if diffusivity<0.0{return Err(err("diffusivity must be nonnegative"));}let out=fd::solve_heat_3d(*n,diffusivity,dt,steps,&initial,bc,stepping(&method)?).map_err(err)?;Ok(history(out.into_iter().map(NativeField::D3).collect(),&g.frame,dt))}))});
object_binding!(fdm_fdm_solve_wave_3d,"num_grid_fdm_solve_wave_3d","crate::sim::fdm_nd::solve_wave_3d","Original native solve_wave_3d: only its actual controls are exposed. The plain 2D/3D zero-boundary entry points preserve zero boundary semantics; *_with_bc variants use explicit boundary objects. All numeric data use the retained Length/Time/field references, native x-fastest layout and selected method; work/history bounds are enforced.",true,(grid:Handle<FdGrid>,initial:Vec<f64>,dt:f64,steps:usize,speed:f64,initial_velocity:Vec<f64>)->NewObject<History>,|ctx,grid,initial,dt,steps,speed,initial_velocity|{grid.with(ctx,|g|{let NativeGrid::D3(n)=&g.data else{return Err(err("FDM grid dimension mismatch"));};
if initial.len()!=g.len()||dt<=0.0||steps>8192||g.len().checked_mul(steps+1).is_none_or(|n|n>LIMIT){return Err(err("FDM initial/time/retained-history budget invalid"));}
if speed<0.0||initial_velocity.len()!=initial.len(){return Err(err("wave speed/initial velocity invalid"));}let out=fd::solve_wave_3d(*n,speed,dt,steps,&initial,&initial_velocity).map_err(err)?;Ok(history(out.into_iter().map(NativeField::D3).collect(),&g.frame,dt))})});
object_binding!(fdm_fdm_solve_wave_3d_with_bc,"num_grid_fdm_solve_wave_3d_with_bc","crate::sim::fdm_nd::solve_wave_3d_with_bc","Original native solve_wave_3d_with_bc: only its actual controls are exposed. The plain 2D/3D zero-boundary entry points preserve zero boundary semantics; *_with_bc variants use explicit boundary objects. All numeric data use the retained Length/Time/field references, native x-fastest layout and selected method; work/history bounds are enforced.",true,(grid:Handle<FdGrid>,bcs:Handle<Bcs>,initial:Vec<f64>,dt:f64,steps:usize,speed:f64,initial_velocity:Vec<f64>)->NewObject<History>,|ctx,grid,bcs,initial,dt,steps,speed,initial_velocity|{grid.with(ctx,|g|bcs.with(ctx,|b|{let NativeGrid::D3(n)=&g.data else{return Err(err("FDM grid dimension mismatch"));};g.frame.same(&b.frame)?;let NativeBcs::D3(bc)=&b.data else{return Err(err("FDM boundary dimension mismatch"));};profiles(g,b)?;if initial.len()!=g.len()||dt<=0.0||steps>8192||g.len().checked_mul(steps+1).is_none_or(|n|n>LIMIT){return Err(err("FDM initial/time/retained-history budget invalid"));}
if speed<0.0||initial_velocity.len()!=initial.len(){return Err(err("wave speed/initial velocity invalid"));}let out=fd::solve_wave_3d_with_bc(*n,speed,dt,steps,&initial,&initial_velocity,bc).map_err(err)?;Ok(history(out.into_iter().map(NativeField::D3).collect(),&g.frame,dt))}))});
object_binding!(fdm_fdm_solve_advection_3d,"num_grid_fdm_solve_advection_3d","crate::sim::fdm_nd::solve_advection_3d","Original native solve_advection_3d: only its actual controls are exposed. The plain 2D/3D zero-boundary entry points preserve zero boundary semantics; *_with_bc variants use explicit boundary objects. All numeric data use the retained Length/Time/field references, native x-fastest layout and selected method; work/history bounds are enforced.",true,(grid:Handle<FdGrid>,bcs:Handle<Bcs>,initial:Vec<f64>,dt:f64,steps:usize,velocity:Vec<f64>)->NewObject<History>,|ctx,grid,bcs,initial,dt,steps,velocity|{grid.with(ctx,|g|bcs.with(ctx,|b|{let NativeGrid::D3(n)=&g.data else{return Err(err("FDM grid dimension mismatch"));};g.frame.same(&b.frame)?;let NativeBcs::D3(bc)=&b.data else{return Err(err("FDM boundary dimension mismatch"));};profiles(g,b)?;if initial.len()!=g.len()||dt<=0.0||steps>8192||g.len().checked_mul(steps+1).is_none_or(|n|n>LIMIT){return Err(err("FDM initial/time/retained-history budget invalid"));}
if velocity.len()!=3{return Err(err("advection velocity dimension mismatch"));}let out=fd::solve_advection_3d(*n,(velocity[0],velocity[1],velocity[2]),dt,steps,&initial,bc).map_err(err)?;Ok(history(out.into_iter().map(NativeField::D3).collect(),&g.frame,dt))}))});
object_binding!(fdm_fdm_history_len,"num_grid_fdm_history_len","","Original native low-level FDM method with explicit Length/Time/field references. Coordinates, steps and coefficients are normalized; output fields/histories retain their physical field reference. Native x-fastest layout, selected integrator and stability/domain checks are preserved with allocation/history budgets.",true,(history:Handle<History>)->usize,|ctx,history|{let _=ctx;history.with(ctx,|h|Ok(h.data.len()))});
object_binding!(fdm_fdm_history_at,"num_grid_fdm_history_at","","Original native low-level FDM method with explicit Length/Time/field references. Coordinates, steps and coefficients are normalized; output fields/histories retain their physical field reference. Native x-fastest layout, selected integrator and stability/domain checks are preserved with allocation/history budgets.",true,(history:Handle<History>,index:usize)->NewObject<Field>,|ctx,history,index|{let _=ctx;history.with(ctx,|h|Ok(field_result(h.data.get(index).ok_or_else(||err("FDM history index out of range"))?.clone(),&h.frame)))});
object_binding!(fdm_fdm_field_data,"num_grid_fdm_field_data","","Original native low-level FDM method with explicit Length/Time/field references. Coordinates, steps and coefficients are normalized; output fields/histories retain their physical field reference. Native x-fastest layout, selected integrator and stability/domain checks are preserved with allocation/history budgets.",true,(field:Handle<Field>)->Dynamic,|ctx,field|{let _=ctx;field.with(ctx,|f|payload(&f.frame,f.data.values().to_vec(),vec![f.data.values().len()]))});
object_binding!(fdm_fdm_history_data,"num_grid_fdm_history_data","","Original native low-level FDM method with explicit Length/Time/field references. Coordinates, steps and coefficients are normalized; output fields/histories retain their physical field reference. Native x-fastest layout, selected integrator and stability/domain checks are preserved with allocation/history budgets.",true,(history:Handle<History>)->Dynamic,|ctx,history|{let _=ctx;history.with(ctx,|h|{let cols=h.data.first().map_or(0,|f|f.values().len());let data=h.data.iter().flat_map(NativeField::values).copied().collect::<Vec<_>>();Ok(Dynamic(api::record(vec![("field",payload(&h.frame,data,vec![h.data.len(),cols])?.0),("dt",api::quantity(h.dt*h.frame.time,h.frame.dims[1].clone())?)])))})});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(i, BINDINGS, "numerical_vector")
}
static BINDINGS: &[&api::Binding] = &[
    &fdm_adaptive_evolve::BINDING,
    &fdm_faces_new::BINDING,
    &fdm_fdm_grid::BINDING,
    &fdm_fdm_zero_dirichlet::BINDING,
    &fdm_fdm_zero_neumann::BINDING,
    &fdm_fdm_condition::BINDING,
    &fdm_fdm1_boundary_uniform::BINDING,
    &fdm_fdm1_boundary_periodic::BINDING,
    &fdm_fdm1_len::BINDING,
    &fdm_fdm1_is_empty::BINDING,
    &fdm_fdm2_boundary_uniform::BINDING,
    &fdm_fdm2_boundary_periodic::BINDING,
    &fdm_fdm2_len::BINDING,
    &fdm_fdm2_is_empty::BINDING,
    &fdm_fdm3_boundary_uniform::BINDING,
    &fdm_fdm3_boundary_periodic::BINDING,
    &fdm_fdm3_len::BINDING,
    &fdm_fdm3_is_empty::BINDING,
    &fdm_fdm3_index::BINDING,
    &fdm_fdm3_field_new::BINDING,
    &fdm_fdm3_field_at::BINDING,
    &fdm_fdm_solve_poisson_1d::BINDING,
    &fdm_fdm_solve_heat_1d::BINDING,
    &fdm_fdm_solve_wave_1d::BINDING,
    &fdm_fdm_solve_advection_1d::BINDING,
    &fdm_fdm_solve_poisson::BINDING,
    &fdm_fdm_solve_poisson_with_bc::BINDING,
    &fdm_fdm_solve_heat::BINDING,
    &fdm_fdm_solve_heat_with_bc::BINDING,
    &fdm_fdm_solve_wave::BINDING,
    &fdm_fdm_solve_wave_with_bc::BINDING,
    &fdm_fdm_solve_advection::BINDING,
    &fdm_fdm_solve_advection_with_bc::BINDING,
    &fdm_fdm_solve_poisson_3d::BINDING,
    &fdm_fdm_solve_heat_3d::BINDING,
    &fdm_fdm_solve_wave_3d::BINDING,
    &fdm_fdm_solve_wave_3d_with_bc::BINDING,
    &fdm_fdm_solve_advection_3d::BINDING,
    &fdm_fdm_history_len::BINDING,
    &fdm_fdm_history_at::BINDING,
    &fdm_fdm_field_data::BINDING,
    &fdm_fdm_history_data::BINDING,
];

fn profile(
    c: &fd::BoundaryCondition,
    n: usize,
) -> Result<(), OpError> {
    if matches!(c,fd::BoundaryCondition::DirichletProfile(v)if v.len()!=n) {
        return Err(err(
            "Dirichlet profile length must exactly match its boundary face",
        ));
    }
    Ok(())
}
fn profiles(
    g: &FdGrid,
    b: &Bcs,
) -> Result<(), OpError> {
    match (&g.data, &b.data) {
        | (NativeGrid::D1(_), NativeBcs::D1(b)) => {
            profile(&b.left, 1)?;
            profile(&b.right, 1)?;
        },
        | (NativeGrid::D2(g), NativeBcs::D2(b)) => {
            profile(&b.left, g.ny)?;
            profile(&b.right, g.ny)?;
            profile(&b.bottom, g.nx)?;
            profile(&b.top, g.nx)?;
        },
        | (NativeGrid::D3(g), NativeBcs::D3(b)) => {
            profile(&b.x_min, g.ny * g.nz)?;
            profile(&b.x_max, g.ny * g.nz)?;
            profile(&b.y_min, g.nx * g.nz)?;
            profile(&b.y_max, g.nx * g.nz)?;
            profile(&b.z_min, g.nx * g.ny)?;
            profile(&b.z_max, g.nx * g.ny)?;
        },
        | _ => return Err(err("FDM grid/boundary dimension mismatch")),
    }
    Ok(())
}
object_binding!(fdm_faces_new,"num_grid_fdm_faces_new","","Native per-face constant/profile FDM boundary construction from exactly2/4/6 scalar condition capabilities. Every condition must use the same physical reference frame; solver calls validate profile face lengths and native paired periodicity.",true,(frame:Handle<Frame>,conditions:Input)->NewObject<Bcs>,|ctx,frame,conditions|{let Value::Tensor(t)=&conditions.0 else{return Err(err("expected boundary condition capability vector"));};
if t.shape!=[t.data.len()]||![2,4,6].contains(&t.data.len()){return Err(err("FDM face list needs exactly2/4/6 conditions"));}frame.with(ctx,|f|{let mut c=t.data.iter().map(|v|api::get::<Condition,_>(ctx,v,Condition::KIND,|c|{f.same(&c.frame)?;Ok(c.data.clone())})).collect::<Result<Vec<_>,_>>()?.into_iter();let data=match t.data.len(){2=>NativeBcs::D1(fd::Boundary1D{left:c.next().ok_or_else(||err("missing face"))?,right:c.next().ok_or_else(||err("missing face"))?}),4=>NativeBcs::D2(fd::Boundary2D{left:c.next().ok_or_else(||err("missing face"))?,right:c.next().ok_or_else(||err("missing face"))?,bottom:c.next().ok_or_else(||err("missing face"))?,top:c.next().ok_or_else(||err("missing face"))?}),6=>NativeBcs::D3(fd::Boundary3D{x_min:c.next().ok_or_else(||err("missing face"))?,x_max:c.next().ok_or_else(||err("missing face"))?,y_min:c.next().ok_or_else(||err("missing face"))?,y_max:c.next().ok_or_else(||err("missing face"))?,z_min:c.next().ok_or_else(||err("missing face"))?,z_max:c.next().ok_or_else(||err("missing face"))?}),_=>return Err(err("invalid face list"))};Ok(NewObject::new(Bcs{data,frame:f.clone()},f.dims.to_vec()))})});
struct AdaptiveControl(fd::adaptive::PdeAdaptiveOptions);
impl NativeArg for AdaptiveControl {
    fn ty() -> TypeExpr {
        let mut fields = ["t0", "t1", "rtol", "atol", "h0", "h_min", "h_max"]
            .into_iter()
            .map(|n| (n.into(), TypeExpr::Float))
            .collect::<Vec<_>>();
        fields.extend(
            ["max_steps", "max_step_calls", "max_rejections"]
                .into_iter()
                .map(|n| (n.into(), TypeExpr::Int)),
        );
        TypeExpr::Record(fields)
    }

    fn decode(
        value: &Value,
        index: usize,
        _: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        api::normalized(value)?;
        let Value::Record(r) = value else {
            return Err(err("adaptive controls must be a record"));
        };
        if r.fields.len() != 10 {
            return Err(err("adaptive controls need exactly ten fields"));
        }
        let f = |n| {
            <f64 as crate::rules::numerical_abi::codecs::Arg>::decode(api::field(value, n)?, index)
        };
        let u = |n| {
            <usize as crate::rules::numerical_abi::codecs::Arg>::decode(
                api::field(value, n)?,
                index,
            )
        };
        let o = fd::adaptive::PdeAdaptiveOptions {
            t0: f("t0")?,
            t1: f("t1")?,
            rtol: f("rtol")?,
            atol: f("atol")?,
            h0: f("h0")?,
            h_min: f("h_min")?,
            h_max: f("h_max")?,
            max_steps: u("max_steps")?,
            max_step_calls: u("max_step_calls")?,
            max_rejections: u("max_rejections")?,
        };
        if o.max_steps == 0
            || o.max_steps > 4096
            || o.max_step_calls == 0
            || o.max_step_calls > 16384
            || o.max_rejections == 0
            || o.max_rejections > 4096
        {
            return Err(err("adaptive PDE control work budget exceeded"));
        }
        Ok(Self(o))
    }
}
object_binding!(fdm_adaptive_evolve,"num_grid_fdm_adaptive_evolve","crate::sim::fdm_nd::adaptive::adaptive_evolve","Actual generic native adaptive controller over caller-provided compiled step/stability/event functions. Each scalar step output takes coherent [time,step,state...] and has field units; stability takes [time,state...] and returns Time, optional event returns dimensionless. Controls/state arrays are normalized in the explicit frame; order is the supplied stepper's actual order. Native step doubling, rejection/work budgets and terminal-event localization are preserved.",true,(frame:Handle<Frame>,initial:Vec<f64>,order:u32,absolute_tolerances:Vec<f64>,controls:AdaptiveControl,step_programs:Input,stability_program:Input,event_program:Input,event_direction:i32,budget:usize)->Dynamic,|ctx,frame,initial,order,absolute_tolerances,controls,step_programs,stability_program,event_program,event_direction,budget|{
 frame.with(ctx,|f|{let n=initial.len();if n==0||n>61||absolute_tolerances.len()!=n||(controls.0.max_steps+1).checked_mul(n).is_none_or(|w|w>LIMIT)||!(-1..=1).contains(&event_direction){return Err(err("adaptive state/retention/direction budget invalid"));}let Value::Tensor(programs)=&step_programs.0 else{return Err(err("step map must be a program vector"));};
if programs.shape!=[n]{return Err(err("one scalar step program required per state coefficient"));}let mut step_domains=vec![f.dims[1].clone();2];step_domains.extend(vec![f.dims[2].clone();n]);let mut state_domains=vec![f.dims[1].clone()];state_domains.extend(vec![f.dims[2].clone();n]);let mut steps=Vec::new();for p in &programs.data{steps.push(crate::rules::numerical_abi::native_grids::structured_ops::checked_snapshot(ctx,&Input(p.clone()),budget,&step_domains,&f.dims[2])?);}let stability=if matches!(&stability_program.0,Value::Variant(v)if v.tag=="None"&&v.payload.is_empty()){None}else{Some(crate::rules::numerical_abi::native_grids::structured_ops::checked_snapshot(ctx,&stability_program,budget,&state_domains,&f.dims[1])?)};let event=if matches!(&event_program.0,Value::Variant(v)if v.tag=="None"&&v.payload.is_empty()){None}else{Some(crate::rules::numerical_abi::native_grids::structured_ops::checked_snapshot(ctx,&event_program,budget,&state_domains,&DimExpr::Dimless)?)};
 let failure=|e:OpError|fd::FdmError::Domain(e.to_string());let step=|t:f64,y:&[f64],h:f64|{let mut input=vec![t*f.time,h*f.time];input.extend(y.iter().map(|v|v*f.value));let out=steps.iter().map(|s|s.sample(&input)/f.value).collect::<Vec<_>>();for s in &steps{s.finish().map_err(failure)?;}Ok(out)};let stable=|t:f64,y:&[f64]|{for s in &steps{s.finish().map_err(failure)?;}let Some(s)=&stability else{return Ok(f64::INFINITY);};let mut input=vec![t*f.time];input.extend(y.iter().map(|v|v*f.value));let out=s.sample(&input)/f.time;s.finish().map_err(failure)?;Ok(out)};let mut event_fn=|t:f64,y:&[f64]|{let mut input=vec![t*f.time];input.extend(y.iter().map(|v|v*f.value));let out=event.as_ref().map_or(0.0,|s|s.sample(&input));if let Some(s)=&event{s.finish().map_err(failure)?;}Ok(out)};let mut terminal=fd::adaptive::PdeTerminalEvent{direction:i8::try_from(event_direction).map_err(|_|err("invalid event direction"))?,evaluate:&mut event_fn};let out=fd::adaptive::adaptive_evolve(&initial,order,&absolute_tolerances,&controls.0,step,stable,if event.is_some(){Some(&mut terminal)}else{None});for s in steps{s.finish()?;}
if let Some(s)=stability{s.finish()?;}
if let Some(s)=event{s.finish()?;}let out=out.map_err(err)?;let times=Value::Tensor(TensorValue{shape:vec![out.times.len()],elem:TypeExpr::Float.with_dim(f.dims[1].clone()),data:out.times.into_iter().map(|t|crate::rules::numerical_abi::output_scalar(t*f.time)).collect::<Result<Vec<_>,_>>()?});let event_time=Value::Variant(crate::kernel::value::VariantValue{tag:if out.event_time.is_some(){"Some"}else{"None"}.into(),payload:out.event_time.map(|t|api::quantity(t*f.time,f.dims[1].clone())).transpose()?.into_iter().collect()});Ok(Dynamic(api::record(vec![("times",times),("states",out.states.encode()?),("value_reference",api::quantity(f.value,f.dims[2].clone())?),("accepted",out.accepted.encode()?),("rejected",out.rejected.encode()?),("step_calls",out.step_calls.encode()?),("error_ratio",out.error_ratio.encode()?),("event_time",event_time)])))})
});
