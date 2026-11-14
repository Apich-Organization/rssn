//! Native simplex finite elements with retained physical source/coefficient callbacks.
use super::super::native_object_api::Maybe;
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
use super::api;
use super::boundary_models::Markers;
use super::canonical;
use super::err;
use super::object_binding;
use super::payload;
use super::structured_ops::checked_snapshot;
use super::structured_ops::dims_forcing;
use super::structured_ops::space;
use super::structured_ops::st;
use crate::sim::fem_simplex as fe;
#[derive(Clone)]
pub(super) struct Mesh {
    pub data: fe::SimplexMesh,
    pub frame: Frame,
}
impl NativeType for Mesh {
    const KIND: &'static str = "rssn/native/simplex-mesh/v1";
    const MAX_BYTES: usize = 4 * 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        self.frame.validate()?;
        let m = &self.data;
        if !(1..=3).contains(&m.dim)
            || m.nodes.len() > 4096
            || m.cells.len() > 8192
            || m.nodes.iter().flatten().any(|v| !v.is_finite())
            || m.cells
                .iter()
                .any(|c| c.len() != m.dim + 1 || c.iter().any(|i| *i >= m.nodes.len()))
            || m.facets
                .iter()
                .any(|(c, f, _)| *c >= m.cells.len() || *f > m.dim)
            || m.level.len() != m.cells.len()
        {
            return Err(err("simplex mesh shape/index/allocation budget invalid"));
        }
        for i in 0..m.cells.len() {
            let measure = m.measure(i);
            if !measure.is_finite() || measure <= 0.0 {
                return Err(err("simplex cells must have finite positive measure"));
            }
        }
        Ok(())
    }
}
#[derive(Clone)]
pub(super) struct Space {
    pub(super) data: fe::FemSpace,
    pub(super) frame: Frame,
}
impl NativeType for Space {
    const KIND: &'static str = "rssn/native/fem-space/v1";
    const MAX_BYTES: usize = 8 * 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        self.frame.validate()?;
        if self.data.size() > 1024 || self.data.mesh.cells.len() > 8192 {
            return Err(err("finite-element space work/allocation budget exceeded"));
        }
        Ok(())
    }
}
fn size(
    mesh: &fe::SimplexMesh,
    degree: usize,
) -> Result<(), OpError> {
    if degree != 1 && degree != 2 {
        return Err(err("FEM degree must be one or two"));
    }
    let mut edges = std::collections::BTreeSet::new();
    if degree == 2 {
        for c in &mesh.cells {
            for (i, &a) in c.iter().enumerate() {
                for &b in &c[i + 1..] {
                    edges.insert((a.min(b), a.max(b)));
                }
            }
        }
    }
    if mesh
        .nodes
        .len()
        .checked_add(edges.len())
        .is_none_or(|n| n > 1024)
    {
        return Err(err("finite-element DOF budget exceeded"));
    }
    Ok(())
}
#[derive(Clone)]
struct Constant {
    frame: Frame,
    k: f64,
    c: f64,
    source: CompiledSnapshot,
    bcs: Markers,
}
impl NativeType for Constant {
    const KIND: &'static str = "rssn/native/fem-constant-problem/v1";
    const MAX_BYTES: usize = 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        self.frame.validate()?;
        if self.k <= 0.0 || !self.k.is_finite() || !self.c.is_finite() {
            return Err(err(
                "FEM diffusivity/reaction must be finite, diffusivity positive",
            ));
        }
        Ok(())
    }
}
#[derive(Clone)]
struct Variable {
    frame: Frame,
    k: CompiledSnapshot,
    c: CompiledSnapshot,
    source: CompiledSnapshot,
    bcs: Markers,
}
impl NativeType for Variable {
    const KIND: &'static str = "rssn/native/fem-variable-problem/v1";
    const MAX_BYTES: usize = 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        self.frame.validate()
    }
}
fn boundary_finish(
    b: &Markers,
    mesh: &fe::SimplexMesh,
) -> Result<(), OpError> {
    let ids = mesh
        .facets
        .iter()
        .map(|(_, _, m)| *m)
        .collect::<std::collections::BTreeSet<_>>();
    for marker in ids {
        for datum in b.callbacks.get(&marker).unwrap_or(&b.default) {
            datum.finish()?;
        }
    }
    Ok(())
}
fn run<R>(
    p: &Constant,
    mesh: &fe::SimplexMesh,
    f: impl FnOnce(&fe::FemProblem<'_>) -> Result<R, OpError>,
) -> Result<R, OpError> {
    let source_reference = p.frame.reference([0, -1, 1])?;
    let source = |t, x| st(&p.source, &p.frame, t, x, source_reference);
    let native = fe::FemProblem {
        k: p.k,
        c: p.c,
        f: &source,
        bcs: &p.bcs.data,
    };
    let out = f(&native);
    p.source.finish()?;
    boundary_finish(&p.bcs, mesh)?;
    out
}
fn run_var<R>(
    p: &Variable,
    mesh: &fe::SimplexMesh,
    f: impl FnOnce(&fe::FemProblemVar<'_>) -> Result<R, OpError>,
) -> Result<R, OpError> {
    let diffusion_reference = p.frame.reference([2, -1, 0])?;
    let reaction_reference = p.frame.reference([0, -1, 0])?;
    let invalid_k = std::sync::atomic::AtomicBool::new(false);
    let k = |x| {
        let value = space(&p.k, &p.frame, x, diffusion_reference);
        if value <= 0.0 {
            invalid_k.store(true, std::sync::atomic::Ordering::Relaxed);
        }
        value
    };
    let c = |x| space(&p.c, &p.frame, x, reaction_reference);
    let source_reference = p.frame.reference([0, -1, 1])?;
    let source = |t, x| st(&p.source, &p.frame, t, x, source_reference);
    let native = fe::FemProblemVar {
        k: &k,
        c: &c,
        f: &source,
        bcs: &p.bcs.data,
    };
    let out = f(&native);
    p.source.finish()?;
    p.k.finish()?;
    p.c.finish()?;
    boundary_finish(&p.bcs, mesh)?;
    if invalid_k.load(std::sync::atomic::Ordering::Relaxed) {
        return Err(err(
            "variable diffusivity must be positive at every evaluated quadrature point",
        ));
    }
    out
}
fn field(
    s: &Space,
    u: &[f64],
) -> Result<(), OpError> {
    if u.len() != s.data.size() {
        return Err(err("finite-element state length differs from DOFs"));
    }
    Ok(())
}
fn indicator_result(
    values: Vec<f64>,
    s: &Space,
) -> Result<Dynamic, OpError> {
    let power = Rational::new(s.data.mesh.dim as i64 + 2, 2);
    let d = canonical(DimExpr::Mul(vec![
        s.frame.dims[2].clone(),
        DimExpr::Pow(Box::new(s.frame.dims[0].clone()), Exp::Rat(power)),
        DimExpr::Inverse(Box::new(s.frame.dims[1].clone())),
    ]))?;
    let scale = s.frame.value
        * s.frame
            .length
            .powf(f64::midpoint(s.data.mesh.dim as f64, 2.0))
        / s.frame.time;
    Ok(Dynamic(api::record(vec![
        ("indicators", values.encode()?),
        ("reference", api::quantity(scale, d)?),
    ])))
}
object_binding!(fem_simplex_interval,"num_grid_simplex_interval","crate::sim::fem_simplex::SimplexMesh::interval","Native simplex finite-element method with explicit geometry/time/field references. Diffusivity is normalized by L²/T, reaction by1/T, source by U/T; compiled callbacks use coherent coordinates and retained physical output domains. Mesh/DOF/refinement/callback work bounds are checked.",true,(frame:Handle<Frame>,n:usize,lo:f64,hi:f64)->NewObject<Mesh>,|ctx,frame,n,lo,hi|{let _=ctx;frame.with(ctx,|f|{if n==0||n>511||lo>=hi{return Err(err("interval mesh geometry/allocation budget invalid"));}Ok(NewObject::new(Mesh{data:fe::SimplexMesh::interval(n,lo,hi),frame:f.clone()},f.dims.to_vec()))})});
object_binding!(fem_simplex_rectangle,"num_grid_simplex_rectangle","crate::sim::fem_simplex::SimplexMesh::rectangle","Native simplex finite-element method with explicit geometry/time/field references. Diffusivity is normalized by L²/T, reaction by1/T, source by U/T; compiled callbacks use coherent coordinates and retained physical output domains. Mesh/DOF/refinement/callback work bounds are checked.",true,(frame:Handle<Frame>,nx:usize,ny:usize,lo:Vec<f64>,hi:Vec<f64>)->NewObject<Mesh>,|ctx,frame,nx,ny,lo,hi|{let _=ctx;frame.with(ctx,|f|{if nx==0||ny==0||nx>64||ny>64||(nx+1).checked_mul(ny+1).is_none_or(|n|n>4096){return Err(err("rectangle mesh allocation budget exceeded"));}let l:[f64;2]=lo.try_into().map_err(|_|err("rectangle needs two lower bounds"))?;let h:[f64;2]=hi.try_into().map_err(|_|err("rectangle needs two upper bounds"))?;if l.iter().zip(h).any(|(l,h)|*l>=h){return Err(err("mesh bounds must increase"));}Ok(NewObject::new(Mesh{data:fe::SimplexMesh::rectangle(nx,ny,l,h),frame:f.clone()},f.dims.to_vec()))})});
object_binding!(fem_simplex_boxed,"num_grid_simplex_boxed","crate::sim::fem_simplex::SimplexMesh::boxed","Native simplex finite-element method with explicit geometry/time/field references. Diffusivity is normalized by L²/T, reaction by1/T, source by U/T; compiled callbacks use coherent coordinates and retained physical output domains. Mesh/DOF/refinement/callback work bounds are checked.",true,(frame:Handle<Frame>,cells:Vec<usize>,lo:Vec<f64>,hi:Vec<f64>)->NewObject<Mesh>,|ctx,frame,cells,lo,hi|{let _=ctx;frame.with(ctx,|f|{let n:[usize;3]=cells.try_into().map_err(|_|err("box needs three cell counts"))?;if n.iter().any(|n|*n==0||*n>16)||n.iter().try_fold(1usize,|a,n|a.checked_mul(n+1)).is_none_or(|n|n>4096)||n.iter().product::<usize>()>1365{return Err(err("box mesh allocation budget exceeded"));}let l:[f64;3]=lo.try_into().map_err(|_|err("box needs three lower bounds"))?;let h:[f64;3]=hi.try_into().map_err(|_|err("box needs three upper bounds"))?;if l.iter().zip(h).any(|(l,h)|*l>=h){return Err(err("mesh bounds must increase"));}Ok(NewObject::new(Mesh{data:fe::SimplexMesh::boxed(n,l,h),frame:f.clone()},f.dims.to_vec()))})});
object_binding!(fem_simplex_refine_marked,"num_grid_simplex_refine_marked","crate::sim::fem_simplex::SimplexMesh::refine_marked","Native simplex finite-element method with explicit geometry/time/field references. Diffusivity is normalized by L²/T, reaction by1/T, source by U/T; compiled callbacks use coherent coordinates and retained physical output domains. Mesh/DOF/refinement/callback work bounds are checked.",true,(mesh:Handle<Mesh>,marked:Vec<bool>)->NewObject<Mesh>,|ctx,mesh,marked|{let _=ctx;mesh.with(ctx,|m|{if marked.len()!=m.data.cells.len()||marked.iter().any(|v|*v)&&(m.data.cells.len()>1024||m.data.nodes.len()+6*m.data.cells.len()>4096){return Err(err("marked simplex refinement shape/allocation budget exceeded"));}Ok(NewObject::new(Mesh{data:m.data.refine_marked(&marked),frame:m.frame.clone()},m.frame.dims.to_vec()))})});
object_binding!(fem_simplex_measure,"num_grid_simplex_measure","crate::sim::fem_simplex::SimplexMesh::measure","Native simplex finite-element method with explicit geometry/time/field references. Diffusivity is normalized by L²/T, reaction by1/T, source by U/T; compiled callbacks use coherent coordinates and retained physical output domains. Mesh/DOF/refinement/callback work bounds are checked.",true,(mesh:Handle<Mesh>,cell:usize)->Dynamic,|ctx,mesh,cell|{let _=ctx;mesh.with(ctx,|m|{if cell>=m.data.cells.len(){return Err(err("simplex cell index out of range"));}let d=canonical(DimExpr::Pow(Box::new(m.frame.dims[0].clone()),Exp::Rat(Rational::new(m.data.dim as i64,1))))?;Ok(Dynamic(api::quantity(m.data.measure(cell)*m.frame.length.powi(m.data.dim as i32),d)?))})});
object_binding!(fem_simplex_bary_gradients,"num_grid_simplex_bary_gradients","crate::sim::fem_simplex::SimplexMesh::bary_gradients","Native simplex finite-element method with explicit geometry/time/field references. Diffusivity is normalized by L²/T, reaction by1/T, source by U/T; compiled callbacks use coherent coordinates and retained physical output domains. Mesh/DOF/refinement/callback work bounds are checked.",true,(mesh:Handle<Mesh>,cell:usize)->Dynamic,|ctx,mesh,cell|{let _=ctx;mesh.with(ctx,|m|{if cell>=m.data.cells.len(){return Err(err("simplex cell index out of range"));}let d=canonical(DimExpr::Inverse(Box::new(m.frame.dims[0].clone())))?;Ok(Dynamic(api::record(vec![("gradients",m.data.bary_gradients(cell).encode()?),("reference",api::quantity(1.0/m.frame.length,d)?)])))})});
object_binding!(fem_simplex_mesh_data,"num_grid_simplex_mesh_data","","Native simplex finite-element method with explicit geometry/time/field references. Diffusivity is normalized by L²/T, reaction by1/T, source by U/T; compiled callbacks use coherent coordinates and retained physical output domains. Mesh/DOF/refinement/callback work bounds are checked.",true,(mesh:Handle<Mesh>)->Dynamic,|ctx,mesh|{let _=ctx;mesh.with(ctx,|m|Ok(Dynamic(api::record(vec![("dimension",m.data.dim.encode()?),("nodes",m.data.nodes.clone().encode()?),("cells",m.data.cells.clone().encode()?),("length_reference",api::quantity(m.frame.length,m.frame.dims[0].clone())?)]))))});
object_binding!(fem_fem_space_new,"num_grid_fem_space_new","crate::sim::fem_simplex::FemSpace::new","Native simplex finite-element method with explicit geometry/time/field references. Diffusivity is normalized by L²/T, reaction by1/T, source by U/T; compiled callbacks use coherent coordinates and retained physical output domains. Mesh/DOF/refinement/callback work bounds are checked.",true,(mesh:Handle<Mesh>,degree:usize)->NewObject<Space>,|ctx,mesh,degree|{let _=ctx;mesh.with(ctx,|m|{size(&m.data,degree)?;Ok(NewObject::new(Space{data:fe::FemSpace::new(m.data.clone(),degree).map_err(err)?,frame:m.frame.clone()},m.frame.dims.to_vec()))})});
object_binding!(fem_fem_space_size,"num_grid_fem_space_size","crate::sim::fem_simplex::FemSpace::size","Native simplex finite-element method with explicit geometry/time/field references. Diffusivity is normalized by L²/T, reaction by1/T, source by U/T; compiled callbacks use coherent coordinates and retained physical output domains. Mesh/DOF/refinement/callback work bounds are checked.",true,(space:Handle<Space>)->usize,|ctx,space|{let _=ctx;space.with(ctx,|s|Ok(s.data.size()))});
object_binding!(fem_fem_eval,"num_grid_fem_eval","crate::sim::fem_simplex::FemSpace::eval","Native simplex finite-element method with explicit geometry/time/field references. Diffusivity is normalized by L²/T, reaction by1/T, source by U/T; compiled callbacks use coherent coordinates and retained physical output domains. Mesh/DOF/refinement/callback work bounds are checked.",false,(space:Handle<Space>,values:Vec<f64>,point:Input)->Maybe<Dynamic>,|ctx,space,values,point|{let _=ctx;space.with(ctx,|s|{field(s,&values)?;api::same_dimension(&api::dimension(&point.0)?,&s.frame.dims[0])?;let p=crate::rules::numerical_abi::tensor(&point.0,1,2)?.1;if p.len()!=3{return Err(err("FEM point needs three coherent coordinates"));}Ok(Maybe(s.data.eval(&values,[p[0]/s.frame.length,p[1]/s.frame.length,p[2]/s.frame.length]).map(|v|api::quantity(v*s.frame.value,s.frame.dims[2].clone()).map(Dynamic)).transpose()?))})});
object_binding!(fem_fem_interpolate,"num_grid_fem_interpolate","crate::sim::fem_simplex::FemSpace::interpolate","Native simplex finite-element method with explicit geometry/time/field references. Diffusivity is normalized by L²/T, reaction by1/T, source by U/T; compiled callbacks use coherent coordinates and retained physical output domains. Mesh/DOF/refinement/callback work bounds are checked.",true,(space:Handle<Space>,program:Input,budget:usize)->Dynamic,|ctx,space,program,budget|{let _=ctx;space.with(ctx,|s|{let callback=checked_snapshot(ctx,&program,budget,&vec![s.frame.dims[0].clone();3],&s.frame.dims[2])?;let exact=|x|crate::rules::numerical_abi::native_grids::structured_ops::space(&callback,&s.frame,x,s.frame.value);let out=s.data.interpolate(exact);callback.finish()?;payload(&s.frame,out,vec![s.data.size()])})});
object_binding!(fem_fem_l2_error,"num_grid_fem_l2_error","crate::sim::fem_simplex::FemSpace::l2_error","Native simplex finite-element method with explicit geometry/time/field references. Diffusivity is normalized by L²/T, reaction by1/T, source by U/T; compiled callbacks use coherent coordinates and retained physical output domains. Mesh/DOF/refinement/callback work bounds are checked.",true,(space:Handle<Space>,program:Input,budget:usize,values:Vec<f64>)->Dynamic,|ctx,space,program,budget,values|{let _=ctx;space.with(ctx,|s|{let callback=checked_snapshot(ctx,&program,budget,&vec![s.frame.dims[0].clone();3],&s.frame.dims[2])?;let exact=|x|crate::rules::numerical_abi::native_grids::structured_ops::space(&callback,&s.frame,x,s.frame.value);field(s,&values)?;let out=s.data.l2_error(&values,exact);callback.finish()?;let d=canonical(DimExpr::Mul(vec![s.frame.dims[2].clone(),DimExpr::Pow(Box::new(s.frame.dims[0].clone()),Exp::Rat(Rational::new(s.data.mesh.dim as i64,2)))]))?;Ok(Dynamic(api::quantity(out*s.frame.value*s.frame.length.powf(s.data.mesh.dim as f64/2.0),d)?))})});
object_binding!(fem_fem_constant_problem,"num_grid_fem_constant_problem","","Native simplex finite-element method with explicit geometry/time/field references. Diffusivity is normalized by L²/T, reaction by1/T, source by U/T; compiled callbacks use coherent coordinates and retained physical output domains. Mesh/DOF/refinement/callback work bounds are checked.",true,(bcs:Handle<Markers>,diffusivity:f64,reaction:f64,source:Input,budget:usize)->NewObject<Constant>,|ctx,bcs,diffusivity,reaction,source,budget|{let _=ctx;bcs.with(ctx,|b|{let f=&b.frame;let callback=checked_snapshot(ctx,&source,budget,&[f.dims[1].clone(),f.dims[0].clone(),f.dims[0].clone(),f.dims[0].clone()],&dims_forcing(f,1)?)?;Ok(NewObject::dependent(Constant{frame:f.clone(),k:diffusivity,c:reaction,source:callback,bcs:b.clone()},f.dims.to_vec(),vec![bcs.value.clone(),source.0.clone()]))})});
object_binding!(fem_fem_variable_problem,"num_grid_fem_variable_problem","","Native simplex finite-element method with explicit geometry/time/field references. Diffusivity is normalized by L²/T, reaction by1/T, source by U/T; compiled callbacks use coherent coordinates and retained physical output domains. Mesh/DOF/refinement/callback work bounds are checked.",true,(bcs:Handle<Markers>,diffusivity:Input,reaction:Input,source:Input,budget:usize)->NewObject<Variable>,|ctx,bcs,diffusivity,reaction,source,budget|{let _=ctx;bcs.with(ctx,|b|{let f=&b.frame;let kdim=canonical(DimExpr::Mul(vec![DimExpr::Pow(Box::new(f.dims[0].clone()),Exp::Rat(Rational::new(2,1))),DimExpr::Inverse(Box::new(f.dims[1].clone()))]))?;let cdim=canonical(DimExpr::Inverse(Box::new(f.dims[1].clone())))?;let domain=vec![f.dims[0].clone();3];let k=checked_snapshot(ctx,&diffusivity,budget,&domain,&kdim)?;let c=checked_snapshot(ctx,&reaction,budget,&domain,&cdim)?;let src=checked_snapshot(ctx,&source,budget,&[f.dims[1].clone(),f.dims[0].clone(),f.dims[0].clone(),f.dims[0].clone()],&dims_forcing(f,1)?)?;Ok(NewObject::dependent(Variable{frame:f.clone(),k,c,source:src,bcs:b.clone()},f.dims.to_vec(),vec![bcs.value.clone(),source.0.clone(),diffusivity.0.clone(),reaction.0.clone()]))})});
object_binding!(fem_fem_solve_static,"num_grid_fem_solve_static","crate::sim::fem_simplex::solve_static","Native simplex finite-element method with explicit geometry/time/field references. Diffusivity is normalized by L²/T, reaction by1/T, source by U/T; compiled callbacks use coherent coordinates and retained physical output domains. Mesh/DOF/refinement/callback work bounds are checked.",true,(space:Handle<Space>,problem:Handle<Constant>,time:f64,solver_name:String,solver_tolerance:f64)->Dynamic,|ctx,space,problem,time,solver_name,solver_tolerance|{let _=ctx;space.with(ctx,|s|problem.with(ctx,|p|{s.frame.same(&p.frame)?;let out=run(p,&s.data.mesh,|native|fe::solve_static(&s.data,native,time,crate::rules::numerical_abi::native_grids::structured_ops::solver(&solver_name,solver_tolerance)?).map_err(err))?;payload(&s.frame,out,vec![s.data.size()])}))});
object_binding!(fem_fem_error_indicators,"num_grid_fem_error_indicators","crate::sim::fem_simplex::error_indicators","Native simplex finite-element method with explicit geometry/time/field references. Diffusivity is normalized by L²/T, reaction by1/T, source by U/T; compiled callbacks use coherent coordinates and retained physical output domains. Mesh/DOF/refinement/callback work bounds are checked.",true,(space:Handle<Space>,problem:Handle<Constant>,time:f64,values:Vec<f64>)->Dynamic,|ctx,space,problem,time,values|{let _=ctx;space.with(ctx,|s|problem.with(ctx,|p|{s.frame.same(&p.frame)?;field(s,&values)?;let out=run(p,&s.data.mesh,|native|Ok(fe::error_indicators(&s.data,&values,native,time)))?;indicator_result(out,s)}))});
object_binding!(fem_fem_solve_heat,"num_grid_fem_solve_heat","crate::sim::fem_simplex::solve_heat","Native simplex finite-element method with explicit geometry/time/field references. Diffusivity is normalized by L²/T, reaction by1/T, source by U/T; compiled callbacks use coherent coordinates and retained physical output domains. Mesh/DOF/refinement/callback work bounds are checked.",true,(space:Handle<Space>,problem:Handle<Constant>,time:f64,initial:Input,budget:usize,steps:usize,theta:f64,solver_name:String,solver_tolerance:f64)->Dynamic,|ctx,space,problem,time,initial,budget,steps,theta,solver_name,solver_tolerance|{let _=ctx;space.with(ctx,|s|problem.with(ctx,|p|{s.frame.same(&p.frame)?;if time<=0.0||steps==0||steps>8192||s.data.size().checked_mul(steps).is_none_or(|n|n>1_048_576)||!(0.0..=1.0).contains(&theta){return Err(err("FEM heat control/work budget invalid"));}let init=checked_snapshot(ctx,&initial,budget,&vec![s.frame.dims[0].clone();3],&s.frame.dims[2])?;let initial=|x|crate::rules::numerical_abi::native_grids::structured_ops::space(&init,&s.frame,x,s.frame.value);let out=run(p,&s.data.mesh,|native|fe::solve_heat(&s.data,native,initial,time,steps,theta,crate::rules::numerical_abi::native_grids::structured_ops::solver(&solver_name,solver_tolerance)?).map_err(err))?;init.finish()?;payload(&s.frame,out,vec![s.data.size()])}))});
object_binding!(fem_fem_solve_static_var,"num_grid_fem_solve_static_var","crate::sim::fem_simplex::solve_static_var","Native simplex finite-element method with explicit geometry/time/field references. Diffusivity is normalized by L²/T, reaction by1/T, source by U/T; compiled callbacks use coherent coordinates and retained physical output domains. Mesh/DOF/refinement/callback work bounds are checked.",true,(space:Handle<Space>,problem:Handle<Variable>,time:f64,solver_name:String,solver_tolerance:f64)->Dynamic,|ctx,space,problem,time,solver_name,solver_tolerance|{let _=ctx;space.with(ctx,|s|problem.with(ctx,|p|{s.frame.same(&p.frame)?;let out=run_var(p,&s.data.mesh,|native|fe::solve_static_var(&s.data,native,time,crate::rules::numerical_abi::native_grids::structured_ops::solver(&solver_name,solver_tolerance)?).map_err(err))?;payload(&s.frame,out,vec![s.data.size()])}))});
object_binding!(fem_fem_error_indicators_var,"num_grid_fem_error_indicators_var","crate::sim::fem_simplex::error_indicators_var","Native simplex finite-element method with explicit geometry/time/field references. Diffusivity is normalized by L²/T, reaction by1/T, source by U/T; compiled callbacks use coherent coordinates and retained physical output domains. Mesh/DOF/refinement/callback work bounds are checked.",true,(space:Handle<Space>,problem:Handle<Variable>,time:f64,values:Vec<f64>)->Dynamic,|ctx,space,problem,time,values|{let _=ctx;space.with(ctx,|s|problem.with(ctx,|p|{s.frame.same(&p.frame)?;field(s,&values)?;let out=run_var(p,&s.data.mesh,|native|Ok(fe::error_indicators_var(&s.data,&values,native,time)))?;indicator_result(out,s)}))});
object_binding!(fem_fem_solve_heat_var,"num_grid_fem_solve_heat_var","crate::sim::fem_simplex::solve_heat_var","Native simplex finite-element method with explicit geometry/time/field references. Diffusivity is normalized by L²/T, reaction by1/T, source by U/T; compiled callbacks use coherent coordinates and retained physical output domains. Mesh/DOF/refinement/callback work bounds are checked.",true,(space:Handle<Space>,problem:Handle<Variable>,time:f64,initial:Input,budget:usize,steps:usize,theta:f64,solver_name:String,solver_tolerance:f64)->Dynamic,|ctx,space,problem,time,initial,budget,steps,theta,solver_name,solver_tolerance|{let _=ctx;space.with(ctx,|s|problem.with(ctx,|p|{s.frame.same(&p.frame)?;if time<=0.0||steps==0||steps>8192||s.data.size().checked_mul(steps).is_none_or(|n|n>1_048_576)||!(0.0..=1.0).contains(&theta){return Err(err("FEM heat control/work budget invalid"));}let init=checked_snapshot(ctx,&initial,budget,&vec![s.frame.dims[0].clone();3],&s.frame.dims[2])?;let initial=|x|crate::rules::numerical_abi::native_grids::structured_ops::space(&init,&s.frame,x,s.frame.value);let out=run_var(p,&s.data.mesh,|native|fe::solve_heat_var(&s.data,native,initial,time,steps,theta,crate::rules::numerical_abi::native_grids::structured_ops::solver(&solver_name,solver_tolerance)?).map_err(err))?;init.finish()?;payload(&s.frame,out,vec![s.data.size()])}))});
object_binding!(fem_fem_heat_step,"num_grid_fem_heat_step","crate::sim::fem_simplex::heat_step","Native simplex finite-element method with explicit geometry/time/field references. Diffusivity is normalized by L²/T, reaction by1/T, source by U/T; compiled callbacks use coherent coordinates and retained physical output domains. Mesh/DOF/refinement/callback work bounds are checked.",true,(space:Handle<Space>,problem:Handle<Constant>,state:Vec<f64>,time:f64,dt:f64,theta:f64,solver_name:String,solver_tolerance:f64)->Dynamic,|ctx,space,problem,state,time,dt,theta,solver_name,solver_tolerance|{let _=ctx;space.with(ctx,|s|problem.with(ctx,|p|{s.frame.same(&p.frame)?;field(s,&state)?;if dt<=0.0||!(0.0..=1.0).contains(&theta){return Err(err("FEM step controls invalid"));}let out=run(p,&s.data.mesh,|native|fe::heat_step(&s.data,native,&state,time,dt,theta,crate::rules::numerical_abi::native_grids::structured_ops::solver(&solver_name,solver_tolerance)?).map_err(err))?;payload(&s.frame,out,vec![s.data.size()])}))});
object_binding!(fem_fem_heat_step_nodal_source,"num_grid_fem_heat_step_nodal_source","crate::sim::fem_simplex::heat_step_nodal_source","Native simplex finite-element method with explicit geometry/time/field references. Diffusivity is normalized by L²/T, reaction by1/T, source by U/T; compiled callbacks use coherent coordinates and retained physical output domains. Mesh/DOF/refinement/callback work bounds are checked.",true,(space:Handle<Space>,problem:Handle<Constant>,state:Vec<f64>,time:f64,dt:f64,theta:f64,before:Vec<f64>,after:Vec<f64>,solver_name:String,solver_tolerance:f64)->Dynamic,|ctx,space,problem,state,time,dt,theta,before,after,solver_name,solver_tolerance|{let _=ctx;space.with(ctx,|s|problem.with(ctx,|p|{s.frame.same(&p.frame)?;field(s,&state)?;if dt<=0.0||!(0.0..=1.0).contains(&theta){return Err(err("FEM step controls invalid"));}field(s,&before)?;field(s,&after)?;let out=run(p,&s.data.mesh,|native|fe::heat_step_nodal_source(&s.data,native,&state,time,dt,theta,crate::rules::numerical_abi::native_grids::structured_ops::solver(&solver_name,solver_tolerance)?,&before,&after).map_err(err))?;payload(&s.frame,out,vec![s.data.size()])}))});
object_binding!(fem_simplex_rule,"num_grid_simplex_rule","crate::sim::fem_simplex::simplex_rule","Native simplex finite-element method with explicit geometry/time/field references. Diffusivity is normalized by L²/T, reaction by1/T, source by U/T; compiled callbacks use coherent coordinates and retained physical output domains. Mesh/DOF/refinement/callback work bounds are checked.",true,(dimension:usize,order:usize)->(Vec<Vec<f64>>,Vec<f64>),|ctx,dimension,order|{let _=ctx;if dimension>6||order>8{return Err(err("simplex quadrature work budget exceeded"));}Ok(fe::simplex_rule(dimension,order))});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(i, BINDINGS, "numerical_vector")
}
static BINDINGS: &[&api::Binding] = &[
    &fem_solve_adaptive::BINDING,
    &fem_solve_adaptive_var::BINDING,
    &fem_simplex_interval::BINDING,
    &fem_simplex_rectangle::BINDING,
    &fem_simplex_boxed::BINDING,
    &fem_simplex_refine_marked::BINDING,
    &fem_simplex_measure::BINDING,
    &fem_simplex_bary_gradients::BINDING,
    &fem_simplex_mesh_data::BINDING,
    &fem_fem_space_new::BINDING,
    &fem_fem_space_size::BINDING,
    &fem_fem_eval::BINDING,
    &fem_fem_interpolate::BINDING,
    &fem_fem_l2_error::BINDING,
    &fem_fem_constant_problem::BINDING,
    &fem_fem_variable_problem::BINDING,
    &fem_fem_solve_static::BINDING,
    &fem_fem_error_indicators::BINDING,
    &fem_fem_solve_heat::BINDING,
    &fem_fem_solve_static_var::BINDING,
    &fem_fem_error_indicators_var::BINDING,
    &fem_fem_solve_heat_var::BINDING,
    &fem_fem_heat_step::BINDING,
    &fem_fem_heat_step_nodal_source::BINDING,
    &fem_simplex_rule::BINDING,
];
object_binding!(fem_solve_adaptive,"num_grid_fem_solve_adaptive","crate::sim::fem_simplex::solve_adaptive","Actual native spatial Dörfler refinement, independently from temporal adaptivity. Positive Length/Time/field references determine coefficient/source units and residual-indicator reference. A conservative worst-case refinement/DOF bound is checked before native allocation; actual returned space, field and estimator history are retained.",true,(mesh:Handle<Mesh>,problem:Handle<Constant>,degree:usize,time:f64,theta:f64,tolerance:f64,max_refinements:usize,solver_name:String,solver_tolerance:f64)->Dynamic,|ctx,mesh,problem,degree,time,theta,tolerance,max_refinements,solver_name,solver_tolerance|{mesh.with(ctx,|m|problem.with(ctx,|p|{m.frame.same(&p.frame)?;if max_refinements==0||max_refinements>4||theta<=0.0||theta>1.0||tolerance<=0.0{return Err(err("spatial adaptivity controls/refinement budget invalid"));}size(&m.data,degree)?;let growth=(1usize<<m.data.dim).checked_pow(max_refinements as u32).ok_or_else(||err("refinement size overflow"))?;let cells=m.data.cells.len().checked_mul(growth).ok_or_else(||err("refinement size overflow"))?;let nodes=m.data.nodes.len().checked_add(cells.checked_mul(if degree==2{12}else{6}).ok_or_else(||err("refinement size overflow"))?).ok_or_else(||err("refinement size overflow"))?;if cells>8192||nodes>1024{return Err(err("worst-case spatial refinement/DOF allocation budget exceeded"));}let out=run(p,&m.data,|native|fe::solve_adaptive(&m.data,native,degree,time,theta,tolerance,max_refinements,crate::rules::numerical_abi::native_grids::structured_ops::solver(&solver_name,solver_tolerance)?).map_err(err))?;let dofs=out.0.size();let space=Space{data:out.0,frame:m.frame.clone()};let indicator=indicator_result(out.2.iter().map(|(_,e)|*e).collect(),&space)?.0;let retained=crate::rules::numerical_abi::native_object_api::NativeOutput::encode(NewObject::new(space,m.frame.dims.to_vec()),ctx)?;Ok(Dynamic(api::record(vec![("space",retained),("solution",payload(&m.frame,out.1,vec![dofs])?.0),("history_dofs",out.2.into_iter().map(|(n,_)|n).collect::<Vec<_>>().encode()?),("history_indicator",indicator)]))) }))});
object_binding!(fem_solve_adaptive_var,"num_grid_fem_solve_adaptive_var","crate::sim::fem_simplex::solve_adaptive_var","Actual native spatial Dörfler refinement, independently from temporal adaptivity. Positive Length/Time/field references determine coefficient/source units and residual-indicator reference. A conservative worst-case refinement/DOF bound is checked before native allocation; actual returned space, field and estimator history are retained.",true,(mesh:Handle<Mesh>,problem:Handle<Variable>,degree:usize,time:f64,theta:f64,tolerance:f64,max_refinements:usize,solver_name:String,solver_tolerance:f64)->Dynamic,|ctx,mesh,problem,degree,time,theta,tolerance,max_refinements,solver_name,solver_tolerance|{mesh.with(ctx,|m|problem.with(ctx,|p|{m.frame.same(&p.frame)?;if max_refinements==0||max_refinements>4||theta<=0.0||theta>1.0||tolerance<=0.0{return Err(err("spatial adaptivity controls/refinement budget invalid"));}size(&m.data,degree)?;let growth=(1usize<<m.data.dim).checked_pow(max_refinements as u32).ok_or_else(||err("refinement size overflow"))?;let cells=m.data.cells.len().checked_mul(growth).ok_or_else(||err("refinement size overflow"))?;let nodes=m.data.nodes.len().checked_add(cells.checked_mul(if degree==2{12}else{6}).ok_or_else(||err("refinement size overflow"))?).ok_or_else(||err("refinement size overflow"))?;if cells>8192||nodes>1024{return Err(err("worst-case spatial refinement/DOF allocation budget exceeded"));}let out=run_var(p,&m.data,|native|fe::solve_adaptive_var(&m.data,native,degree,time,theta,tolerance,max_refinements,crate::rules::numerical_abi::native_grids::structured_ops::solver(&solver_name,solver_tolerance)?).map_err(err))?;let dofs=out.0.size();let space=Space{data:out.0,frame:m.frame.clone()};let indicator=indicator_result(out.2.iter().map(|(_,e)|*e).collect(),&space)?.0;let retained=crate::rules::numerical_abi::native_object_api::NativeOutput::encode(NewObject::new(space,m.frame.dims.to_vec()),ctx)?;Ok(Dynamic(api::record(vec![("space",retained),("solution",payload(&m.frame,out.1,vec![dofs])?.0),("history_dofs",out.2.into_iter().map(|(n,_)|n).collect::<Vec<_>>().encode()?),("history_indicator",indicator)]))) }))});
