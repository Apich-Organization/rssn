//! Boundary elements in an explicit length/potential reference frame.
//! Numeric arrays are normalized coordinates or field coefficients, never implicit SI.
use super::super::native_object_api::Maybe;
use super::CallCtx;
use super::DimExpr;
use super::Dynamic;
use super::Exp;
use super::Frame;
use super::Handle;
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
use super::object_binding;
use crate::sim::bem_aca as aca;
use crate::sim::bem_helmholtz as helm;
use crate::sim::physics_bem as be;
use num_complex::Complex64;
#[derive(Clone)]
struct Mesh {
    data: be::SurfaceMesh3D,
    frame: Frame,
}
impl NativeType for Mesh {
    const KIND: &'static str = "rssn/native/bem-mesh/v1";
    const MAX_BYTES: usize = 262144;

    fn validate(&self) -> Result<(), OpError> {
        self.frame.validate()?;
        if self.data.vertices.len() > 1024 || self.data.triangles.len() > 512 {
            return Err(err("BEM mesh allocation/work cap exceeded"));
        }
        self.data.validate().map_err(err)
    }
}
#[derive(Clone)]
struct Vector {
    data: Vec<f64>,
    frame: Frame,
}
impl NativeType for Vector {
    const KIND: &'static str = "rssn/native/bem-vector/v1";
    const MAX_BYTES: usize = 4096;

    fn validate(&self) -> Result<(), OpError> {
        self.frame.validate()?;
        if !(2..=3).contains(&self.data.len()) || self.data.iter().any(|x| !x.is_finite()) {
            return Err(err("BEM vector invalid"));
        }
        Ok(())
    }
}
struct Element {
    data: be::Element2D,
    frame: Frame,
}
impl NativeType for Element {
    const KIND: &'static str = "rssn/native/bem-element2/v1";
    const MAX_BYTES: usize = 4096;

    fn validate(&self) -> Result<(), OpError> {
        self.frame.validate()?;
        if !self.data.length.is_finite() || self.data.length <= 0.0 {
            return Err(err("BEM element must have positive finite length"));
        }
        Ok(())
    }
}
struct Solution {
    mesh: Mesh,
    u: Vec<f64>,
    q: Vec<f64>,
    p1: bool,
    stats: Option<aca::AcaStats>,
}
impl NativeType for Solution {
    const KIND: &'static str = "rssn/native/bem-solution/v1";
    const MAX_BYTES: usize = 524288;

    fn validate(&self) -> Result<(), OpError> {
        self.mesh.validate()?;
        let n = if self.p1 {
            self.mesh.data.vertices.len()
        } else {
            self.mesh.data.triangles.len()
        };
        if self.u.len() != n
            || self.q.len() != n
            || self.u.iter().chain(&self.q).any(|x| !x.is_finite())
        {
            return Err(err("BEM solution shape/finite invariant failed"));
        }
        Ok(())
    }
}
struct Helmholtz {
    mesh: Mesh,
    data: helm::HelmholtzSolution,
    stats: Option<helm::HelmholtzAcaStats>,
}
impl NativeType for Helmholtz {
    const KIND: &'static str = "rssn/native/helmholtz-bem/v1";
    const MAX_BYTES: usize = 524288;

    fn validate(&self) -> Result<(), OpError> {
        self.mesh.validate()?;
        if self.data.density.len() != self.mesh.data.triangles.len()
            || self
                .data
                .density
                .iter()
                .any(|z| !z.re.is_finite() || !z.im.is_finite())
        {
            return Err(err("Helmholtz solution invalid"));
        }
        Ok(())
    }
}
fn xyz(v: Vec<f64>) -> Result<[f64; 3], OpError> {
    v.try_into()
        .map_err(|_| err("three normalized coordinates required"))
}
fn bcs(
    kinds: &[String],
    values: &[f64],
    n: usize,
) -> Result<Vec<be::BoundaryCondition<f64>>, OpError> {
    if kinds.len() != n || values.len() != n {
        return Err(err(
            "one boundary kind and coefficient per element/vertex required",
        ));
    }
    kinds
        .iter()
        .zip(values)
        .map(|(k, v)| {
            match k.as_str() {
                | "Potential" => Ok(be::BoundaryCondition::Potential(*v)),
                | "Flux" => Ok(be::BoundaryCondition::Flux(*v)),
                | _ => Err(err("BEM kind must be Potential or Flux")),
            }
        })
        .collect()
}
fn options(
    leaf: usize,
    eta: f64,
    tol: f64,
    restart: usize,
    iterations: usize,
) -> Result<aca::AcaOptions, OpError> {
    if leaf == 0
        || leaf > 512
        || eta <= 0.0
        || tol <= 0.0
        || tol > 1.0
        || restart == 0
        || restart > 128
        || iterations == 0
        || iterations > 4096
    {
        return Err(err("ACA options/work cap invalid"));
    }
    Ok(aca::AcaOptions {
        leaf_size: leaf,
        eta,
        aca_tol: tol,
        gmres_tol: tol,
        restart,
        max_iter: iterations,
    })
}
fn fields(
    f: &Frame,
    u: Vec<f64>,
    q: Vec<f64>,
) -> Result<Dynamic, OpError> {
    Ok(Dynamic(api::record(vec![
        ("potential", u.encode()?),
        ("flux", q.encode()?),
        (
            "potential_reference",
            api::quantity(f.value, f.dims[2].clone())?,
        ),
        (
            "flux_reference",
            api::quantity(f.value / f.length, f.gradient()?)?,
        ),
    ])))
}
fn out(
    f: &Frame,
    value: f64,
    power: i64,
) -> Result<Dynamic, OpError> {
    let d = canonical(DimExpr::Pow(
        Box::new(f.dims[0].clone()),
        Exp::Rat(Rational::new(power, 1)),
    ))?;
    Ok(Dynamic(api::quantity(
        value * f.length.powi(i32::try_from(power).map_err(err)?),
        d,
    )?))
}
fn complex_out(
    f: &Frame,
    z: Complex64,
) -> Result<Dynamic, OpError> {
    Ok(Dynamic(api::record(vec![
        ("normalized", z.encode()?),
        ("reference", api::quantity(f.value, f.dims[2].clone())?),
    ])))
}
object_binding!(bem_vector2_new,"num_grid_bem_vector2_new","crate::sim::physics_bem::Vector2D::new","Native BEM call. Numeric coordinates are normalized by the supplied PDE frame length, potentials by its field reference, and fluxes by field/length. Outputs retain explicit physical reference witnesses.",true,(frame:Handle<Frame>,coordinates:Vec<f64>)->NewObject<Vector>,|ctx:&CallCtx<'_> ,frame:Handle<Frame>, coordinates:Vec<f64>|{frame.with(ctx,|f|{if coordinates.len()!=2{return Err(err("vector dimension mismatch"));}let native=be::Vector2D::new(coordinates[0],coordinates[1]);Ok(NewObject::new(Vector{data:vec![native.x,native.y],frame:f.clone()},f.dims.to_vec()))})});
object_binding!(bem_vector2_norm,"num_grid_bem_vector2_norm","crate::sim::physics_bem::Vector2D::norm","Native BEM call. Numeric coordinates are normalized by the supplied PDE frame length, potentials by its field reference, and fluxes by field/length. Outputs retain explicit physical reference witnesses.",true,(vector:Handle<Vector>)->Dynamic,|ctx:&CallCtx<'_> ,vector:Handle<Vector>|{vector.with(ctx,|v|{if v.data.len()!=2{return Err(err("vector dimension mismatch"));}let native=be::Vector2D::new(v.data[0],v.data[1]);out(&v.frame,native.norm(),1)})});
object_binding!(bem_vector3_new,"num_grid_bem_vector3_new","crate::sim::physics_bem::Vector3D::new","Native BEM call. Numeric coordinates are normalized by the supplied PDE frame length, potentials by its field reference, and fluxes by field/length. Outputs retain explicit physical reference witnesses.",true,(frame:Handle<Frame>,coordinates:Vec<f64>)->NewObject<Vector>,|ctx:&CallCtx<'_> ,frame:Handle<Frame>, coordinates:Vec<f64>|{frame.with(ctx,|f|{if coordinates.len()!=3{return Err(err("vector dimension mismatch"));}let native=be::Vector3D::new(coordinates[0],coordinates[1],coordinates[2]);Ok(NewObject::new(Vector{data:vec![native.x,native.y,native.z],frame:f.clone()},f.dims.to_vec()))})});
object_binding!(bem_vector3_norm,"num_grid_bem_vector3_norm","crate::sim::physics_bem::Vector3D::norm","Native BEM call. Numeric coordinates are normalized by the supplied PDE frame length, potentials by its field reference, and fluxes by field/length. Outputs retain explicit physical reference witnesses.",true,(vector:Handle<Vector>)->Dynamic,|ctx:&CallCtx<'_> ,vector:Handle<Vector>|{vector.with(ctx,|v|{if v.data.len()!=3{return Err(err("vector dimension mismatch"));}let native=be::Vector3D::new(v.data[0],v.data[1],v.data[2]);out(&v.frame,native.norm(),1)})});
object_binding!(bem_element_new,"num_grid_bem_element_new","crate::sim::physics_bem::Element2D::new","Native BEM call. Numeric coordinates are normalized by the supplied PDE frame length, potentials by its field reference, and fluxes by field/length. Outputs retain explicit physical reference witnesses.",true,(a:Handle<Vector>,b:Handle<Vector>)->NewObject<Element>,|ctx:&CallCtx<'_> ,a:Handle<Vector>, b:Handle<Vector>|{a.with(ctx,|a|b.with(ctx,|b|{a.frame.same(&b.frame)?;if a.data.len()!=2||b.data.len()!=2||a.data==b.data{return Err(err("distinct 2D endpoints required"));}Ok(NewObject::new(Element{data:be::Element2D::new(be::Vector2D::new(a.data[0],a.data[1]),be::Vector2D::new(b.data[0],b.data[1])),frame:a.frame.clone()},a.frame.dims.to_vec()))}))});
object_binding!(bem_mesh_new,"num_grid_bem_mesh_new","crate::sim::physics_bem::SurfaceMesh3D::validate","Native BEM call. Numeric coordinates are normalized by the supplied PDE frame length, potentials by its field reference, and fluxes by field/length. Outputs retain explicit physical reference witnesses.",true,(frame:Handle<Frame>,vertices:Vec<Vec<f64>>,triangles:Vec<Vec<f64>>)->NewObject<Mesh>,|ctx:&CallCtx<'_> ,frame:Handle<Frame>, vertices:Vec<Vec<f64>>, triangles:Vec<Vec<f64>>|{frame.with(ctx,|f|{if vertices.len()>1024||triangles.len()>512{return Err(err("BEM mesh cap exceeded"));}let vertices=vertices.into_iter().map(xyz).collect::<Result<Vec<_>,_>>()?;let triangles=triangles.into_iter().map(|t|{if t.len()!=3||t.iter().any(|x|*x<0.0||x.fract()!=0.0||*x>=vertices.len() as f64){return Err(err("triangle indices invalid"));}Ok([t[0] as usize,t[1] as usize,t[2] as usize])}).collect::<Result<Vec<_>,OpError>>()?;Ok(NewObject::new(Mesh{data:be::SurfaceMesh3D{vertices,triangles},frame:f.clone()},f.dims.to_vec()))})});
object_binding!(bem_mesh_icosphere,"num_grid_bem_mesh_icosphere","crate::sim::physics_bem::SurfaceMesh3D::icosphere","Native BEM call. Numeric coordinates are normalized by the supplied PDE frame length, potentials by its field reference, and fluxes by field/length. Outputs retain explicit physical reference witnesses.",true,(frame:Handle<Frame>,size:f64,center:Vec<f64>,subdivisions:usize)->NewObject<Mesh>,|ctx:&CallCtx<'_> ,frame:Handle<Frame>, size:f64, center:Vec<f64>, subdivisions:usize|{frame.with(ctx,|f|{if size<=0.0{return Err(err("mesh size must be positive"));}let data={if subdivisions>2{return Err(err("icosphere subdivision cap is two"));}be::SurfaceMesh3D::icosphere(size,xyz(center)?,subdivisions)};Ok(NewObject::new(Mesh{data,frame:f.clone()},f.dims.to_vec()))})});
object_binding!(bem_mesh_cube,"num_grid_bem_mesh_cube","crate::sim::physics_bem::SurfaceMesh3D::cube","Native BEM call. Numeric coordinates are normalized by the supplied PDE frame length, potentials by its field reference, and fluxes by field/length. Outputs retain explicit physical reference witnesses.",true,(frame:Handle<Frame>,size:f64,center:Vec<f64>,subdivisions:usize)->NewObject<Mesh>,|ctx:&CallCtx<'_> ,frame:Handle<Frame>, size:f64, center:Vec<f64>, subdivisions:usize|{frame.with(ctx,|f|{if size<=0.0{return Err(err("mesh size must be positive"));}let data={if subdivisions==0||subdivisions>6{return Err(err("cube subdivision cap is six"));}be::SurfaceMesh3D::cube(size,xyz(center)?,subdivisions)};Ok(NewObject::new(Mesh{data,frame:f.clone()},f.dims.to_vec()))})});
object_binding!(bem_mesh_corners,"num_grid_bem_mesh_corners","crate::sim::physics_bem::SurfaceMesh3D::corners","Native BEM call. Numeric coordinates are normalized by the supplied PDE frame length, potentials by its field reference, and fluxes by field/length. Outputs retain explicit physical reference witnesses.",true,(mesh:Handle<Mesh>,triangle:usize)->Maybe<Dynamic>,|ctx:&CallCtx<'_> ,mesh:Handle<Mesh>, triangle:usize|{mesh.with(ctx,|m|{let value=m.data.corners(triangle);Ok(Maybe(value.map(|v|Ok::<_,OpError>(Dynamic(api::record(vec![("normalized",v.into_iter().map(|p|p.to_vec()).collect::<Vec<_>>().encode()?),("reference",api::quantity(m.frame.length,m.frame.dims[0].clone())?)] )))).transpose()?))})});
object_binding!(bem_mesh_centroid,"num_grid_bem_mesh_centroid","crate::sim::physics_bem::SurfaceMesh3D::centroid","Native BEM call. Numeric coordinates are normalized by the supplied PDE frame length, potentials by its field reference, and fluxes by field/length. Outputs retain explicit physical reference witnesses.",true,(mesh:Handle<Mesh>,triangle:usize)->Maybe<Dynamic>,|ctx:&CallCtx<'_> ,mesh:Handle<Mesh>, triangle:usize|{mesh.with(ctx,|m|{let value=m.data.centroid(triangle);Ok(Maybe(value.map(|v|Ok::<_,OpError>(Dynamic(api::record(vec![("normalized",v.to_vec().encode()?),("reference",api::quantity(m.frame.length,m.frame.dims[0].clone())?)] )))).transpose()?))})});
object_binding!(bem_mesh_normal,"num_grid_bem_mesh_normal","crate::sim::physics_bem::SurfaceMesh3D::normal","Native BEM call. Numeric coordinates are normalized by the supplied PDE frame length, potentials by its field reference, and fluxes by field/length. Outputs retain explicit physical reference witnesses.",true,(mesh:Handle<Mesh>,triangle:usize)->Maybe<Dynamic>,|ctx:&CallCtx<'_> ,mesh:Handle<Mesh>, triangle:usize|{mesh.with(ctx,|m|{let value=m.data.normal(triangle);Ok(Maybe(value.map(|v|Ok::<_,OpError>(Dynamic(api::record(vec![("normalized",v.to_vec().encode()?),("reference",Value::Float(1.0))] )))).transpose()?))})});
object_binding!(bem_mesh_area,"num_grid_bem_mesh_area","crate::sim::physics_bem::SurfaceMesh3D::area","Native BEM call. Numeric coordinates are normalized by the supplied PDE frame length, potentials by its field reference, and fluxes by field/length. Outputs retain explicit physical reference witnesses.",true,(mesh:Handle<Mesh>,triangle:usize)->Maybe<Dynamic>,|ctx:&CallCtx<'_> ,mesh:Handle<Mesh>, triangle:usize|{mesh.with(ctx,|m|{Ok(Maybe(m.data.area(triangle).map(|v|out(&m.frame,v,2)).transpose()?))})});
object_binding!(bem_mesh_total_area,"num_grid_bem_mesh_total_area","crate::sim::physics_bem::SurfaceMesh3D::total_area","Native BEM call. Numeric coordinates are normalized by the supplied PDE frame length, potentials by its field reference, and fluxes by field/length. Outputs retain explicit physical reference witnesses.",true,(mesh:Handle<Mesh>)->Dynamic,|ctx:&CallCtx<'_> ,mesh:Handle<Mesh>|{mesh.with(ctx,|m|out(&m.frame,m.data.total_area(),2))});
object_binding!(bem_mesh_volume,"num_grid_bem_mesh_volume","crate::sim::physics_bem::SurfaceMesh3D::volume","Native BEM call. Numeric coordinates are normalized by the supplied PDE frame length, potentials by its field reference, and fluxes by field/length. Outputs retain explicit physical reference witnesses.",true,(mesh:Handle<Mesh>)->Dynamic,|ctx:&CallCtx<'_> ,mesh:Handle<Mesh>|{mesh.with(ctx,|m|out(&m.frame,m.data.volume(),3))});
object_binding!(bem_mesh_validate,"num_grid_bem_mesh_validate","crate::sim::physics_bem::SurfaceMesh3D::validate","Native BEM call. Numeric coordinates are normalized by the supplied PDE frame length, potentials by its field reference, and fluxes by field/length. Outputs retain explicit physical reference witnesses.",true,(mesh:Handle<Mesh>)->bool,|ctx:&CallCtx<'_> ,mesh:Handle<Mesh>|{mesh.with(ctx,|m|{m.data.validate().map_err(err)?;Ok(true)})});
object_binding!(bem_solve3,"num_grid_bem_solve3","crate::sim::physics_bem::solve_laplace_bem_3d_mesh","Native BEM call. Numeric coordinates are normalized by the supplied PDE frame length, potentials by its field reference, and fluxes by field/length. Outputs retain explicit physical reference witnesses.",true,(mesh:Handle<Mesh>,kinds:Vec<String>,values:Vec<f64>)->NewObject<Solution>,|ctx:&CallCtx<'_> ,mesh:Handle<Mesh>, kinds:Vec<String>, values:Vec<f64>|{mesh.with(ctx,|m|{let bc=bcs(&kinds,&values,m.data.triangles.len())?;let s=be::solve_laplace_bem_3d_mesh(&m.data,&bc).map_err(err)?;Ok(NewObject::new(Solution{mesh:m.clone(),u:s.u,q:s.q,p1:false,stats:None},m.frame.dims.to_vec()))})});
object_binding!(bem_aca,"num_grid_bem_aca","crate::sim::bem_aca::solve_laplace_bem_3d_aca","Native BEM call. Numeric coordinates are normalized by the supplied PDE frame length, potentials by its field reference, and fluxes by field/length. Outputs retain explicit physical reference witnesses.",true,(mesh:Handle<Mesh>,kinds:Vec<String>,values:Vec<f64>,leaf:usize,eta:f64,tolerance:f64,restart:usize,iterations:usize)->NewObject<Solution>,|ctx:&CallCtx<'_> ,mesh:Handle<Mesh>, kinds:Vec<String>, values:Vec<f64>, leaf:usize, eta:f64, tolerance:f64, restart:usize, iterations:usize|{mesh.with(ctx,|m|{let bc=bcs(&kinds,&values,m.data.triangles.len())?;let opts=options(leaf,eta,tolerance,restart,iterations)?;let (s,stats)=aca::solve_laplace_bem_3d_aca(&m.data,&bc,&opts).map_err(err)?;Ok(NewObject::new(Solution{mesh:m.clone(),u:s.u,q:s.q,p1:false,stats:Some(stats)},m.frame.dims.to_vec()))})});
object_binding!(bem_exterior,"num_grid_bem_exterior","crate::sim::bem_aca::solve_laplace_bem_3d_exterior","Native BEM call. Numeric coordinates are normalized by the supplied PDE frame length, potentials by its field reference, and fluxes by field/length. Outputs retain explicit physical reference witnesses.",true,(mesh:Handle<Mesh>,kinds:Vec<String>,values:Vec<f64>,leaf:usize,eta:f64,tolerance:f64,restart:usize,iterations:usize)->NewObject<Solution>,|ctx:&CallCtx<'_> ,mesh:Handle<Mesh>, kinds:Vec<String>, values:Vec<f64>, leaf:usize, eta:f64, tolerance:f64, restart:usize, iterations:usize|{mesh.with(ctx,|m|{let bc=bcs(&kinds,&values,m.data.triangles.len())?;let opts=options(leaf,eta,tolerance,restart,iterations)?;let (s,stats)=aca::solve_laplace_bem_3d_exterior(&m.data,&bc,&opts).map_err(err)?;Ok(NewObject::new(Solution{mesh:m.clone(),u:s.u,q:s.q,p1:false,stats:Some(stats)},m.frame.dims.to_vec()))})});
object_binding!(bem_p1,"num_grid_bem_p1","crate::sim::bem_aca::solve_laplace_bem_3d_p1","Native BEM call. Numeric coordinates are normalized by the supplied PDE frame length, potentials by its field reference, and fluxes by field/length. Outputs retain explicit physical reference witnesses.",true,(mesh:Handle<Mesh>,kinds:Vec<String>,values:Vec<f64>,exterior:bool)->NewObject<Solution>,|ctx:&CallCtx<'_> ,mesh:Handle<Mesh>, kinds:Vec<String>, values:Vec<f64>, exterior:bool|{mesh.with(ctx,|m|{let bc=bcs(&kinds,&values,m.data.vertices.len())?;let s=aca::solve_laplace_bem_3d_p1(&m.data,&bc,exterior).map_err(err)?;Ok(NewObject::new(Solution{mesh:m.clone(),u:s.u,q:s.q,p1:true,stats:None},m.frame.dims.to_vec()))})});
object_binding!(bem_p1_aca,"num_grid_bem_p1_aca","crate::sim::bem_aca::solve_laplace_bem_3d_p1_aca","Native BEM call. Numeric coordinates are normalized by the supplied PDE frame length, potentials by its field reference, and fluxes by field/length. Outputs retain explicit physical reference witnesses.",true,(mesh:Handle<Mesh>,kinds:Vec<String>,values:Vec<f64>,exterior:bool,leaf:usize,eta:f64,tolerance:f64,restart:usize,iterations:usize)->NewObject<Solution>,|ctx:&CallCtx<'_> ,mesh:Handle<Mesh>, kinds:Vec<String>, values:Vec<f64>, exterior:bool, leaf:usize, eta:f64, tolerance:f64, restart:usize, iterations:usize|{mesh.with(ctx,|m|{let bc=bcs(&kinds,&values,m.data.vertices.len())?;let opts=options(leaf,eta,tolerance,restart,iterations)?;let (s,stats)=aca::solve_laplace_bem_3d_p1_aca(&m.data,&bc,exterior,&opts).map_err(err)?;Ok(NewObject::new(Solution{mesh:m.clone(),u:s.u,q:s.q,p1:true,stats:Some(stats)},m.frame.dims.to_vec()))})});
object_binding!(bem_solution_fields,"num_grid_bem_solution_fields","","Native BEM call. Numeric coordinates are normalized by the supplied PDE frame length, potentials by its field reference, and fluxes by field/length. Outputs retain explicit physical reference witnesses.",true,(solution:Handle<Solution>)->Dynamic,|ctx:&CallCtx<'_> ,solution:Handle<Solution>|{solution.with(ctx,|s|fields(&s.mesh.frame,s.u.clone(),s.q.clone()))});
object_binding!(bem_eval3,"num_grid_bem_eval3","crate::sim::physics_bem::evaluate_potential_3d","Native BEM call. Numeric coordinates are normalized by the supplied PDE frame length, potentials by its field reference, and fluxes by field/length. Outputs retain explicit physical reference witnesses.",true,(solution:Handle<Solution>,point:Vec<f64>)->Dynamic,|ctx:&CallCtx<'_> ,solution:Handle<Solution>, point:Vec<f64>|{solution.with(ctx,|s|{if s.p1{return Err(err("constant-element evaluator cannot interpret P1 coefficients"));}let solution=be::Bem3dSolution{u:s.u.clone(),q:s.q.clone()};let u=be::evaluate_potential_3d(xyz(point)?,&s.mesh.data,&solution).map_err(err)?;Ok(Dynamic(api::quantity(u*s.mesh.frame.value,s.mesh.frame.dims[2].clone())?))})});
object_binding!(bem_eval3_exterior,"num_grid_bem_eval3_exterior","crate::sim::bem_aca::evaluate_potential_3d_exterior","Native BEM call. Numeric coordinates are normalized by the supplied PDE frame length, potentials by its field reference, and fluxes by field/length. Outputs retain explicit physical reference witnesses.",true,(solution:Handle<Solution>,point:Vec<f64>)->Dynamic,|ctx:&CallCtx<'_> ,solution:Handle<Solution>, point:Vec<f64>|{solution.with(ctx,|s|{if s.p1{return Err(err("constant-element evaluator cannot interpret P1 coefficients"));}let solution=be::Bem3dSolution{u:s.u.clone(),q:s.q.clone()};let u=aca::evaluate_potential_3d_exterior(xyz(point)?,&s.mesh.data,&solution).map_err(err)?;Ok(Dynamic(api::quantity(u*s.mesh.frame.value,s.mesh.frame.dims[2].clone())?))})});
object_binding!(bem_compression,"num_grid_bem_compression","crate::sim::bem_aca::AcaStats::compression","Native BEM call. Numeric coordinates are normalized by the supplied PDE frame length, potentials by its field reference, and fluxes by field/length. Outputs retain explicit physical reference witnesses.",true,(solution:Handle<Solution>)->f64,|ctx:&CallCtx<'_> ,solution:Handle<Solution>|{solution.with(ctx,|s|s.stats.as_ref().map(aca::AcaStats::compression).ok_or_else(||err("solution was not computed by ACA")))});
object_binding!(bem_cylinder,"num_grid_bem_cylinder","crate::sim::physics_bem::simulate_2d_cylinder_scenario","Native BEM call. Numeric coordinates are normalized by the supplied PDE frame length, potentials by its field reference, and fluxes by field/length. Outputs retain explicit physical reference witnesses.",true,(frame:Handle<Frame>)->Dynamic,|ctx:&CallCtx<'_> ,frame:Handle<Frame>|{frame.with(ctx,|f|{let(u,q)=be::simulate_2d_cylinder_scenario().map_err(err)?;fields(f,u,q)})});
object_binding!(bem_selfcheck3,"num_grid_bem_selfcheck3","crate::sim::physics_bem::solve_laplace_bem_3d","Native BEM call. Numeric coordinates are normalized by the supplied PDE frame length, potentials by its field reference, and fluxes by field/length. Outputs retain explicit physical reference witnesses.",true,()->bool,|ctx:&CallCtx<'_>|{let _=ctx;be::solve_laplace_bem_3d().map_err(err)?;Ok(true)});
object_binding!(bem_solve2,"num_grid_bem_solve2","crate::sim::physics_bem::solve_laplace_bem_2d","Native BEM call. Numeric coordinates are normalized by the supplied PDE frame length, potentials by its field reference, and fluxes by field/length. Outputs retain explicit physical reference witnesses.",true,(frame:Handle<Frame>,points:Vec<Vec<f64>>,kinds:Vec<String>,values:Vec<f64>)->Dynamic,|ctx:&CallCtx<'_> ,frame:Handle<Frame>, points:Vec<Vec<f64>>, kinds:Vec<String>, values:Vec<f64>|{frame.with(ctx,|f|{if points.len()<3||points.len()>512||points.iter().any(|p|p.len()!=2){return Err(err("BEM polygon needs 3..512 two-dimensional points"));}let points=points.iter().map(|p|(p[0],p[1])).collect::<Vec<_>>();for i in 0..points.len(){if points[i]==points[(i+1)%points.len()]{return Err(err("polygon has a zero-length edge"));}}let(u,q)=be::solve_laplace_bem_2d(&points,&bcs(&kinds,&values,points.len())?).map_err(err)?;fields(f,u,q)})});
object_binding!(bem_eval2,"num_grid_bem_eval2","crate::sim::physics_bem::evaluate_potential_2d","Native BEM call. Numeric coordinates are normalized by the supplied PDE frame length, potentials by its field reference, and fluxes by field/length. Outputs retain explicit physical reference witnesses.",true,(frame:Handle<Frame>,point:Vec<f64>,points:Vec<Vec<f64>>,potential:Vec<f64>,flux:Vec<f64>)->Dynamic,|ctx:&CallCtx<'_> ,frame:Handle<Frame>, point:Vec<f64>, points:Vec<Vec<f64>>, potential:Vec<f64>, flux:Vec<f64>|{frame.with(ctx,|f|{let n=points.len();if !(3..=512).contains(&n)||point.len()!=2||points.iter().any(|p|p.len()!=2)||potential.len()!=n||flux.len()!=n{return Err(err("BEM 2D data shape mismatch"));}let mut elements=Vec::with_capacity(n);for i in 0..n{let a=&points[i];let b=&points[(i+1)%n];if a==b{return Err(err("zero-length edge"));}elements.push(be::Element2D::new(be::Vector2D::new(a[0],a[1]),be::Vector2D::new(b[0],b[1])));}let u=be::evaluate_potential_2d((point[0],point[1]),&elements,&potential,&flux);Ok(Dynamic(api::quantity(u*f.value,f.dims[2].clone())?))})});
object_binding!(helmholtz_solve,"num_grid_helmholtz_solve","crate::sim::bem_helmholtz::solve_helmholtz_exterior","Native exterior Helmholtz solver. Wave number is normalized by inverse frame length, coupling by frame length, direction dimensionless. Unit incident amplitude corresponds to the frame potential reference; returned fields retain it. Mesh, dense work, ACA iterations and storage are bounded.",true,(mesh:Handle<Mesh>,wave_number:f64,kind:String,direction:Vec<f64>,coupling:Maybe<Complex64>)->NewObject<Helmholtz>,|ctx:&CallCtx<'_> ,mesh:Handle<Mesh>, wave_number:f64, kind:String, direction:Vec<f64>, coupling:Maybe<Complex64>|{mesh.with(ctx,|m|{if wave_number<=0.0{return Err(err("positive wave number times length-reference required"));}let direction=xyz(direction)?;if direction.iter().all(|x|*x==0.0){return Err(err("nonzero wave direction required"));}let kind=match kind.as_str(){"SoundHard"=>helm::Scatterer::SoundHard,"SoundSoft"=>helm::Scatterer::SoundSoft,_=>return Err(err("scatterer must be SoundHard or SoundSoft"))};let data=helm::solve_helmholtz_exterior(&m.data,wave_number,kind,helm::PlaneWave{direction},coupling.0).map_err(err)?;Ok(NewObject::new(Helmholtz{mesh:m.clone(),data,stats:None},m.frame.dims.to_vec()))})});
object_binding!(helmholtz_solve_aca,"num_grid_helmholtz_solve_aca","crate::sim::bem_helmholtz::solve_helmholtz_exterior_aca","Native exterior Helmholtz solver. Wave number is normalized by inverse frame length, coupling by frame length, direction dimensionless. Unit incident amplitude corresponds to the frame potential reference; returned fields retain it. Mesh, dense work, ACA iterations and storage are bounded.",true,(mesh:Handle<Mesh>,wave_number:f64,kind:String,direction:Vec<f64>,coupling:Maybe<Complex64>,leaf:usize,eta:f64,tolerance:f64,restart:usize,iterations:usize)->NewObject<Helmholtz>,|ctx:&CallCtx<'_> ,mesh:Handle<Mesh>, wave_number:f64, kind:String, direction:Vec<f64>, coupling:Maybe<Complex64>, leaf:usize, eta:f64, tolerance:f64, restart:usize, iterations:usize|{mesh.with(ctx,|m|{if wave_number<=0.0{return Err(err("positive wave number times length-reference required"));}let direction=xyz(direction)?;if direction.iter().all(|x|*x==0.0){return Err(err("nonzero wave direction required"));}let kind=match kind.as_str(){"SoundHard"=>helm::Scatterer::SoundHard,"SoundSoft"=>helm::Scatterer::SoundSoft,_=>return Err(err("scatterer must be SoundHard or SoundSoft"))};let o=options(leaf,eta,tolerance,restart,iterations)?;let opts=helm::HelmholtzAcaOptions{leaf_size:o.leaf_size,eta:o.eta,aca_tol:o.aca_tol,gmres_tol:o.gmres_tol,restart:o.restart,max_iter:o.max_iter};let(data,stats)=helm::solve_helmholtz_exterior_aca(&m.data,wave_number,kind,helm::PlaneWave{direction},coupling.0,&opts).map_err(err)?;Ok(NewObject::new(Helmholtz{mesh:m.clone(),data,stats:Some(stats)},m.frame.dims.to_vec()))})});
object_binding!(helmholtz_incident_field,"num_grid_helmholtz_incident_field","crate::sim::bem_helmholtz::HelmholtzSolution::incident_field","Native BEM call. Numeric coordinates are normalized by the supplied PDE frame length, potentials by its field reference, and fluxes by field/length. Outputs retain explicit physical reference witnesses.",true,(solution:Handle<Helmholtz>,point:Vec<f64>)->Dynamic,|ctx:&CallCtx<'_> ,solution:Handle<Helmholtz>, point:Vec<f64>|{solution.with(ctx,|s|complex_out(&s.mesh.frame,s.data.incident_field(xyz(point)?)))});
object_binding!(helmholtz_scattered_field,"num_grid_helmholtz_scattered_field","crate::sim::bem_helmholtz::HelmholtzSolution::scattered_field","Native BEM call. Numeric coordinates are normalized by the supplied PDE frame length, potentials by its field reference, and fluxes by field/length. Outputs retain explicit physical reference witnesses.",true,(solution:Handle<Helmholtz>,point:Vec<f64>)->Dynamic,|ctx:&CallCtx<'_> ,solution:Handle<Helmholtz>, point:Vec<f64>|{solution.with(ctx,|s|complex_out(&s.mesh.frame,s.data.scattered_field(&s.mesh.data,xyz(point)?)))});
object_binding!(helmholtz_stats,"num_grid_helmholtz_stats","","Native BEM call. Numeric coordinates are normalized by the supplied PDE frame length, potentials by its field reference, and fluxes by field/length. Outputs retain explicit physical reference witnesses.",true,(solution:Handle<Helmholtz>)->Dynamic,|ctx:&CallCtx<'_> ,solution:Handle<Helmholtz>|{solution.with(ctx,|s|{let a=s.stats.as_ref().ok_or_else(||err("solution is not compressed"))?;Ok(Dynamic(api::record(vec![("compression",a.compression.encode()?),("iterations",a.iterations.encode()?),("stored",a.stored.encode()?),("low_rank_blocks",a.low_rank_blocks.encode()?),("dense_blocks",a.dense_blocks.encode()?)])))})});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &bem_vector2_new::BINDING,
            &bem_vector2_norm::BINDING,
            &bem_vector3_new::BINDING,
            &bem_vector3_norm::BINDING,
            &bem_element_new::BINDING,
            &bem_mesh_new::BINDING,
            &bem_mesh_icosphere::BINDING,
            &bem_mesh_cube::BINDING,
            &bem_mesh_corners::BINDING,
            &bem_mesh_centroid::BINDING,
            &bem_mesh_normal::BINDING,
            &bem_mesh_area::BINDING,
            &bem_mesh_total_area::BINDING,
            &bem_mesh_volume::BINDING,
            &bem_mesh_validate::BINDING,
            &bem_solve3::BINDING,
            &bem_aca::BINDING,
            &bem_exterior::BINDING,
            &bem_p1::BINDING,
            &bem_p1_aca::BINDING,
            &bem_solution_fields::BINDING,
            &bem_eval3::BINDING,
            &bem_eval3_exterior::BINDING,
            &bem_compression::BINDING,
            &bem_cylinder::BINDING,
            &bem_selfcheck3::BINDING,
            &bem_solve2::BINDING,
            &bem_eval2::BINDING,
            &helmholtz_solve::BINDING,
            &helmholtz_solve_aca::BINDING,
            &helmholtz_incident_field::BINDING,
            &helmholtz_scattered_field::BINDING,
            &helmholtz_stats::BINDING,
        ],
        "numerical_vector",
    )
}
