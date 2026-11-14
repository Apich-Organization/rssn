//! Legacy quadrilateral plane-stress models. Thickness is one Length reference.
use super::CallCtx;
use super::Dynamic;
use super::Handle;
use super::Installer;
use super::OpError;
use super::Output;
use super::RuleError;
use super::api;
use super::err;
use super::mechanics::References;
use super::object_binding;
use crate::sim::models_linear_elasticity as quad;
fn points(v: Vec<Vec<f64>>) -> Result<Vec<(f64, f64)>, OpError> {
    v.into_iter()
        .map(|p| {
            if p.len() == 2 {
                Ok((p[0], p[1]))
            } else {
                Err(err("quadrilateral nodes need two normalized coordinates"))
            }
        })
        .collect()
}
fn material(
    e: f64,
    nu: f64,
) -> Result<(), OpError> {
    if e <= 0.0 || nu <= -1.0 || nu >= 0.5 {
        return Err(err("quadrilateral material requires E>0,-1<nu<1/2"));
    }
    Ok(())
}
fn output(
    r: &References,
    v: Vec<f64>,
) -> Result<Dynamic, OpError> {
    Ok(Dynamic(api::record(vec![
        ("normalized_displacement", v.encode()?),
        ("reference", api::quantity(r.scales[0], r.dims[0].clone())?),
    ])))
}
object_binding!(quad_element,"num_grid_quad_element","crate::sim::models_linear_elasticity::element_stiffness_matrix","Native bilinear quadrilateral plane-stress call. Coordinates/displacements are normalized by the supplied mechanical Length reference, modulus by stress Mass/(Length Time^2), nodal loads by Force. Native unit thickness means one physical Length reference; stiffness entries therefore have Force/Length units. Element index, material, quadrature geometry, solve work and allocation bounds are checked.",true,(refs:Handle<References>,coordinates:Vec<Vec<f64>>,young:f64,poisson:f64)->Dynamic,|ctx:&CallCtx<'_>,refs:Handle<References>,coordinates:Vec<Vec<f64>>,young:f64,poisson:f64|{material(young,poisson)?;let p=points(coordinates)?;if p.len()!=4{return Err(err("four quadrilateral vertices required"));}quad::validate_quadrilateral([p[0],p[1],p[2],p[3]]).map_err(err)?;refs.with(ctx,|r|{let m=quad::element_stiffness_matrix(p[0],p[1],p[2],p[3],young,poisson);let(s,d)=r.derived([0,-2,1,0])?;Ok(Dynamic(api::record(vec![("rows",8usize.encode()?),("cols",8usize.encode()?),("entries",m.iter().copied().collect::<Vec<_>>().encode()?),("entry_reference",api::quantity(s,d)?)])))})});
object_binding!(quad_solve,"num_grid_quad_solve","crate::sim::models_linear_elasticity::run_elasticity_simulation","Native bilinear quadrilateral plane-stress call. Coordinates/displacements are normalized by the supplied mechanical Length reference, modulus by stress Mass/(Length Time^2), nodal loads by Force. Native unit thickness means one physical Length reference; stiffness entries therefore have Force/Length units. Element index, material, quadrature geometry, solve work and allocation bounds are checked.",true,(refs:Handle<References>,coordinates:Vec<Vec<f64>>,elements:Vec<Vec<f64>>,young:f64,poisson:f64,fixed:Vec<usize>,loads:Vec<Vec<f64>>)->Dynamic,|ctx:&CallCtx<'_>,refs:Handle<References>,coordinates:Vec<Vec<f64>>,elements:Vec<Vec<f64>>,young:f64,poisson:f64,fixed:Vec<usize>,loads:Vec<Vec<f64>>|{material(young,poisson)?;let nodes=points(coordinates)?;if nodes.is_empty()||nodes.len()>256||elements.len()>512||loads.len()>512||fixed.iter().any(|i|*i>=nodes.len())||elements.len().checked_mul(64*5000).is_none_or(|n|n>50000000){return Err(err("quadrilateral solve work/index/allocation cap exceeded"));}let elements=elements.into_iter().map(|e|{if e.len()!=4||e.iter().any(|i|*i<0.0||i.fract()!=0.0||*i>=nodes.len() as f64){return Err(err("quadrilateral connectivity invalid"));}let e=[e[0] as usize,e[1] as usize,e[2] as usize,e[3] as usize];quad::validate_quadrilateral(e.map(|i|nodes[i])).map_err(err)?;Ok(e)}).collect::<Result<Vec<_>,OpError>>()?;let loads=loads.into_iter().map(|v|{if v.len()!=3||v[0]<0.0||v[0].fract()!=0.0||v[0]>=nodes.len() as f64{return Err(err("load row needs [node,fx,fy]"));}Ok((v[0] as usize,v[1],v[2]))}).collect::<Result<Vec<_>,OpError>>()?;refs.with(ctx,|r|{let p=quad::ElasticityParameters{nodes,elements,youngs_modulus:young,poissons_ratio:poisson,fixed_nodes:fixed,loads};output(r,quad::run_elasticity_simulation(&p).map_err(err)?)})});
object_binding!(quad_cantilever,"num_grid_quad_cantilever","crate::sim::models_linear_elasticity::simulate_cantilever_beam_scenario","Execute the actual fixed native cantilever scenario (105 nodes,80 elements, CG maximum5000 iterations), writing beam_original.csv and beam_deformed.csv to the requested directory. Native educational coefficients are explicitly normalized; no measured inputs are accepted.",true,(directory:String)->bool,|ctx:&CallCtx<'_>,directory:String|{let _=ctx;if directory.is_empty()||directory.len()>4096{return Err(err("output directory path invalid"));}quad::simulate_cantilever_beam_scenario(std::path::Path::new(&directory)).map_err(err)?;Ok(true)});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &quad_element::BINDING,
            &quad_solve::BINDING,
            &quad_cantilever::BINDING,
        ],
        "numerical_vector",
    )
}
