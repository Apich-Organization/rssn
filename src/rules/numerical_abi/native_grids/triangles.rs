//! Native triangular P1/P2 finite elements and spatial refinement.
use super::super::native_object_api::NativeOutput;
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
use super::canonical;
use super::err;
use super::object_binding;
use super::payload;
use crate::sim::fem_triangles as tri;
#[derive(Clone)]
struct Mesh {
    data: tri::TriMesh,
    frame: Frame,
}
impl NativeType for Mesh {
    const KIND: &'static str = "rssn/native/triangle-mesh/v1";
    const MAX_BYTES: usize = 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        self.frame.validate()?;
        if self.data.nodes.len() > 1024
            || self.data.triangles.len() > 2048
            || self.data.nodes.iter().flatten().any(|v| !v.is_finite())
            || self
                .data
                .triangles
                .iter()
                .any(|t| t.iter().any(|i| *i >= self.data.nodes.len()))
        {
            return Err(err("triangle mesh shape/index/allocation invalid"));
        }
        for t in &self.data.triangles {
            let [a, b, c] = t.map(|i| self.data.nodes[i]);
            let area = (b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1]);
            if !area.is_finite() || area <= 0.0 {
                return Err(err(
                    "triangles must be counterclockwise with finite positive area",
                ));
            }
        }
        Ok(())
    }
}
#[derive(Clone)]
struct Problem {
    frame: Frame,
    callbacks: [CompiledSnapshot; 4],
}
impl NativeType for Problem {
    const KIND: &'static str = "rssn/native/triangle-problem/v1";
    const MAX_BYTES: usize = 16384;

    fn validate(&self) -> Result<(), OpError> {
        self.frame.validate()
    }
}
#[derive(Clone)]
struct Solution {
    data: tri::Solution,
    mesh: Mesh,
}
impl NativeType for Solution {
    const KIND: &'static str = "rssn/native/triangle-solution/v1";
    const MAX_BYTES: usize = 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        self.mesh.validate()?;
        if self.data.points.len() != self.data.values.len()
            || self.data.values.len() > 1024
            || self.data.values.iter().any(|v| !v.is_finite())
        {
            return Err(err("triangle solution finite storage bounds exceeded"));
        }
        Ok(())
    }
}
fn size(
    m: &tri::TriMesh,
    p: usize,
) -> Result<(), OpError> {
    if p != 1 && p != 2 {
        return Err(err("triangle degree must be 1/2"));
    }
    let mut edges = std::collections::BTreeSet::new();
    if p == 2 {
        for t in &m.triangles {
            for (a, b) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
                edges.insert((a.min(b), a.max(b)));
            }
        }
    }
    if m.nodes.len() + edges.len() > 1024 {
        return Err(err("triangle DOF budget exceeded"));
    }
    Ok(())
}
fn evaluate(
    c: &CompiledSnapshot,
    f: &Frame,
    x: f64,
    y: f64,
    out: f64,
) -> f64 {
    c.sample(&[x * f.length, y * f.length]) / out
}
fn run<T>(
    p: &Problem,
    body: impl FnOnce(&tri::Problem<'_>) -> Result<T, OpError>,
) -> Result<T, OpError> {
    let invalid = std::sync::atomic::AtomicBool::new(false);
    let f = &p.frame;
    let diffusion_reference = f.reference([2, -1, 0])?;
    let reaction_reference = f.reference([0, -1, 0])?;
    let source_reference = f.reference([0, -1, 1])?;
    let k = |x, y| {
        let v = evaluate(&p.callbacks[0], f, x, y, diffusion_reference);
        if !v.is_finite() || v <= 0.0 {
            invalid.store(true, std::sync::atomic::Ordering::Relaxed);
        }
        v
    };
    let c = |x, y| {
        let v = evaluate(&p.callbacks[1], f, x, y, reaction_reference);
        if !v.is_finite() || v < 0.0 {
            invalid.store(true, std::sync::atomic::Ordering::Relaxed);
        }
        v
    };
    let source = |x, y| evaluate(&p.callbacks[2], f, x, y, source_reference);
    let g = |x, y| evaluate(&p.callbacks[3], f, x, y, f.value);
    let result = body(&tri::Problem {
        k: &k,
        c: &c,
        f: &source,
        g: &g,
    });
    for c in &p.callbacks {
        c.finish()?;
    }
    if invalid.load(std::sync::atomic::Ordering::Relaxed) {
        return Err(err(
            "triangle diffusion/reaction must be positive/nonnegative at every evaluated point",
        ));
    }
    result
}
fn measured(
    v: f64,
    f: &Frame,
    length_power: i64,
) -> Result<Dynamic, OpError> {
    let d = canonical(DimExpr::Mul(vec![
        f.dims[2].clone(),
        DimExpr::Pow(
            Box::new(f.dims[0].clone()),
            Exp::Rat(Rational::new(length_power, 1)),
        ),
    ]))?;
    Ok(Dynamic(api::quantity(
        v * f.value * f.length.powi(length_power as i32),
        d,
    )?))
}
fn indicator(
    v: Vec<f64>,
    f: &Frame,
) -> Result<Dynamic, OpError> {
    let d = canonical(DimExpr::Mul(vec![
        f.dims[2].clone(),
        DimExpr::Pow(Box::new(f.dims[0].clone()), Exp::Rat(Rational::new(2, 1))),
        DimExpr::Inverse(Box::new(f.dims[1].clone())),
    ]))?;
    Ok(Dynamic(api::record(vec![
        ("normalized", Output::encode(v)?),
        ("reference", api::quantity(f.reference([2, -1, 1])?, d)?),
    ])))
}
object_binding!(triangle_rect,"num_triangle_rectangle","crate::sim::fem_triangles::TriMesh::rectangle","Actual triangular FEM with explicit Length/Time/field references, variable compiled diffusion/reaction/source/Dirichlet data, original CG/P1/P2/refinement algorithms and bounded native storage/work. Callback arguments and outputs are coherent physical quantities.",true,(frame:Handle<Frame>,bounds:Vec<f64>,nx:usize,ny:usize)->NewObject<Mesh>,|ctx,frame,bounds,nx,ny|{frame.with(ctx,|f|{if bounds.len()!=4||bounds[1]<=bounds[0]||bounds[3]<=bounds[2]||nx==0||ny==0||nx>31||ny>31||(nx+1)*(ny+1)>1024{return Err(err("triangle rectangle geometry/allocation invalid"));}Ok(NewObject::new(Mesh{data:tri::TriMesh::rectangle(bounds[0],bounds[1],bounds[2],bounds[3],nx,ny),frame:f.clone()},f.dims.to_vec()))})});
object_binding!(triangle_l_shape,"num_triangle_l_shape","crate::sim::fem_triangles::TriMesh::l_shape","Actual triangular FEM with explicit Length/Time/field references, variable compiled diffusion/reaction/source/Dirichlet data, original CG/P1/P2/refinement algorithms and bounded native storage/work. Callback arguments and outputs are coherent physical quantities.",true,(frame:Handle<Frame>,resolution:usize)->NewObject<Mesh>,|ctx,frame,resolution|{frame.with(ctx,|f|{if resolution==0||resolution>15{return Err(err("L-shape resolution must be 1..15"));}Ok(NewObject::new(Mesh{data:tri::TriMesh::l_shape(resolution),frame:f.clone()},f.dims.to_vec()))})});
object_binding!(triangle_compacted,"num_triangle_compacted","crate::sim::fem_triangles::TriMesh::compacted","Actual triangular FEM with explicit Length/Time/field references, variable compiled diffusion/reaction/source/Dirichlet data, original CG/P1/P2/refinement algorithms and bounded native storage/work. Callback arguments and outputs are coherent physical quantities.",true,(mesh:Handle<Mesh>)->NewObject<Mesh>,|ctx,mesh|{mesh.with(ctx,|m|Ok(NewObject::new(Mesh{data:m.data.clone().compacted(),frame:m.frame.clone()},m.frame.dims.to_vec())))});
object_binding!(triangle_boundary,"num_triangle_boundary_vertices","crate::sim::fem_triangles::TriMesh::boundary_vertices","Actual triangular FEM with explicit Length/Time/field references, variable compiled diffusion/reaction/source/Dirichlet data, original CG/P1/P2/refinement algorithms and bounded native storage/work. Callback arguments and outputs are coherent physical quantities.",true,(mesh:Handle<Mesh>)->Vec<bool>,|ctx,mesh|{mesh.with(ctx,|m|Ok(m.data.boundary_vertices()))});
object_binding!(triangle_refine,"num_triangle_refine","crate::sim::fem_triangles::TriMesh::refine","Actual triangular FEM with explicit Length/Time/field references, variable compiled diffusion/reaction/source/Dirichlet data, original CG/P1/P2/refinement algorithms and bounded native storage/work. Callback arguments and outputs are coherent physical quantities.",true,(mesh:Handle<Mesh>,marked:Vec<bool>)->NewObject<Mesh>,|ctx,mesh,marked|{mesh.with(ctx,|m|{if marked.len()!=m.data.triangles.len()||m.data.triangles.len()>256||m.data.nodes.len()+3*m.data.triangles.len()>1024{return Err(err("triangle refinement allocation/mark shape budget exceeded"));}Ok(NewObject::new(Mesh{data:m.data.refine(&marked),frame:m.frame.clone()},m.frame.dims.to_vec()))})});
object_binding!(triangle_problem,"num_triangle_problem","","Actual triangular FEM with explicit Length/Time/field references, variable compiled diffusion/reaction/source/Dirichlet data, original CG/P1/P2/refinement algorithms and bounded native storage/work. Callback arguments and outputs are coherent physical quantities.",false,(frame:Handle<Frame>,diffusion:Input,reaction:Input,source:Input,dirichlet:Input,budget:usize)->NewObject<Problem>,|ctx,frame,diffusion,reaction,source,dirichlet,budget|{frame.with(ctx,|f|{if budget==0||budget>8_388_608{return Err(err("triangle callback budget invalid"));}let domain=vec![f.dims[0].clone();2];let kd=canonical(DimExpr::Mul(vec![DimExpr::Pow(Box::new(f.dims[0].clone()),Exp::Rat(Rational::new(2,1))),DimExpr::Inverse(Box::new(f.dims[1].clone()))]))?;let cd=canonical(DimExpr::Inverse(Box::new(f.dims[1].clone())))?;let fd=crate::rules::numerical_abi::native_grids::structured_ops::dims_forcing(f,1)?;let callbacks=[crate::rules::numerical_abi::native_grids::structured_ops::checked_snapshot(ctx,&diffusion,budget,&domain,&kd)?,crate::rules::numerical_abi::native_grids::structured_ops::checked_snapshot(ctx,&reaction,budget,&domain,&cd)?,crate::rules::numerical_abi::native_grids::structured_ops::checked_snapshot(ctx,&source,budget,&domain,&fd)?,crate::rules::numerical_abi::native_grids::structured_ops::checked_snapshot(ctx,&dirichlet,budget,&domain,&f.dims[2])?];Ok(NewObject::dependent(Problem{frame:f.clone(),callbacks},f.dims.to_vec(),vec![diffusion.0,reaction.0,source.0,dirichlet.0]))})});
object_binding!(triangle_solve,"num_triangle_solve","crate::sim::fem_triangles::solve","Actual triangular FEM with explicit Length/Time/field references, variable compiled diffusion/reaction/source/Dirichlet data, original CG/P1/P2/refinement algorithms and bounded native storage/work. Callback arguments and outputs are coherent physical quantities.",true,(mesh:Handle<Mesh>,problem:Handle<Problem>,degree:usize)->NewObject<Solution>,|ctx,mesh,problem,degree|{mesh.with(ctx,|m|problem.with(ctx,|p|{m.frame.same(&p.frame)?;size(&m.data,degree)?;let data=run(p,|p|tri::solve(&m.data,if degree==1{tri::Degree::P1}else{tri::Degree::P2},p).map_err(err))?;Ok(NewObject::new(Solution{data,mesh:m.clone()},m.frame.dims.to_vec()))}))});
object_binding!(triangle_indicators,"num_triangle_error_indicators","crate::sim::fem_triangles::error_indicators","Actual triangular FEM with explicit Length/Time/field references, variable compiled diffusion/reaction/source/Dirichlet data, original CG/P1/P2/refinement algorithms and bounded native storage/work. Callback arguments and outputs are coherent physical quantities.",true,(solution:Handle<Solution>,problem:Handle<Problem>)->Dynamic,|ctx,solution,problem|{solution.with(ctx,|s|problem.with(ctx,|p|{s.mesh.frame.same(&p.frame)?;if s.data.degree!=tri::Degree::P1{return Err(err("residual indicators require a P1 solution"));}let v=run(p,|p|Ok(tri::error_indicators(&s.mesh.data,&s.data,p)))?;indicator(v,&s.mesh.frame)}))});
object_binding!(triangle_adaptive,"num_triangle_solve_adaptive","crate::sim::fem_triangles::solve_adaptive","Actual triangular FEM with explicit Length/Time/field references, variable compiled diffusion/reaction/source/Dirichlet data, original CG/P1/P2/refinement algorithms and bounded native storage/work. Callback arguments and outputs are coherent physical quantities.",true,(mesh:Handle<Mesh>,problem:Handle<Problem>,refinements:usize,theta:f64)->Dynamic,|ctx,mesh,problem,refinements,theta|{mesh.with(ctx,|m|problem.with(ctx,|p|{m.frame.same(&p.frame)?;if refinements>4||theta<=0.0||theta>1.0||m.data.triangles.len().checked_mul(4usize.pow(refinements as u32)).is_none_or(|n|n>512)||m.data.nodes.len()+3*m.data.triangles.len()*4usize.pow(refinements as u32)>1024{return Err(err("spatial refinement worst-case allocation budget exceeded"));}let(mesh,data,history)=run(p,|p|tri::solve_adaptive(&m.data,p,refinements,theta).map_err(err))?;let mesh=Mesh{data:mesh,frame:m.frame.clone()};let sh=NativeOutput::encode(NewObject::new(Solution{data,mesh:mesh.clone()},m.frame.dims.to_vec()),ctx)?;let mh=NativeOutput::encode(NewObject::new(mesh,m.frame.dims.to_vec()),ctx)?;Ok(Dynamic(api::record(vec![("mesh",mh),("solution",sh),("history_dofs",Output::encode(history.iter().map(|v|v.0).collect::<Vec<_>>())?),("history_error",indicator(history.iter().map(|v|v.1).collect(),&m.frame)?.0)])))}))});
object_binding!(triangle_max_nodal_error,"num_triangle_max_nodal_error","crate::sim::fem_triangles::Solution::max_nodal_error","Actual triangular FEM with explicit Length/Time/field references, variable compiled diffusion/reaction/source/Dirichlet data, original CG/P1/P2/refinement algorithms and bounded native storage/work. Callback arguments and outputs are coherent physical quantities.",false,(solution:Handle<Solution>,exact:Input,budget:usize)->Dynamic,|ctx,solution,exact,budget|{solution.with(ctx,|s|{let f=&s.mesh.frame;let c=crate::rules::numerical_abi::native_grids::structured_ops::checked_snapshot(ctx,&exact,budget,&vec![f.dims[0].clone();2],&f.dims[2])?;let exact=|x,y|evaluate(&c,f,x,y,f.value);let v=s.data.max_nodal_error(&exact);c.finish()?;measured(v,f,0)})});
object_binding!(triangle_l2_error,"num_triangle_l2_error","crate::sim::fem_triangles::Solution::l2_error","Actual triangular FEM with explicit Length/Time/field references, variable compiled diffusion/reaction/source/Dirichlet data, original CG/P1/P2/refinement algorithms and bounded native storage/work. Callback arguments and outputs are coherent physical quantities.",false,(solution:Handle<Solution>,exact:Input,budget:usize)->Dynamic,|ctx,solution,exact,budget|{solution.with(ctx,|s|{let f=&s.mesh.frame;let c=crate::rules::numerical_abi::native_grids::structured_ops::checked_snapshot(ctx,&exact,budget,&vec![f.dims[0].clone();2],&f.dims[2])?;let exact=|x,y|evaluate(&c,f,x,y,f.value);let v=s.data.l2_error(&s.mesh.data,&exact);c.finish()?;measured(v,f,1)})});
object_binding!(triangle_data,"num_triangle_data","","Actual triangular FEM with explicit Length/Time/field references, variable compiled diffusion/reaction/source/Dirichlet data, original CG/P1/P2/refinement algorithms and bounded native storage/work. Callback arguments and outputs are coherent physical quantities.",true,(solution:Handle<Solution>)->Dynamic,|ctx,solution|{solution.with(ctx,|s|Ok(Dynamic(api::record(vec![("points_normalized",Output::encode(s.data.points.iter().flat_map(|x|x.iter().copied()).collect::<Vec<_>>())?),("elements",Output::encode(s.data.elements.clone())?),("values",payload(&s.mesh.frame,s.data.values.clone(),vec![s.data.values.len()])?.0),("iterations",Output::encode(s.data.iterations)?)]))))});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(i, BINDINGS, "numerical_vector")
}
static BINDINGS: &[&api::Binding] = &[
    &triangle_rect::BINDING,
    &triangle_l_shape::BINDING,
    &triangle_compacted::BINDING,
    &triangle_boundary::BINDING,
    &triangle_refine::BINDING,
    &triangle_problem::BINDING,
    &triangle_solve::BINDING,
    &triangle_indicators::BINDING,
    &triangle_adaptive::BINDING,
    &triangle_max_nodal_error::BINDING,
    &triangle_l2_error::BINDING,
    &triangle_data::BINDING,
];
