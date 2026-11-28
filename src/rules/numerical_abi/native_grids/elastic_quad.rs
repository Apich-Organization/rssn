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
fn output(
    r: &References,
    v: Vec<f64>,
) -> Result<Dynamic, OpError> {
    Ok(Dynamic(api::record(vec![
        ("normalized_displacement", v.encode()?),
        ("reference", api::quantity(r.scales[0], r.dims[0].clone())?),
    ])))
}
object_binding!(quad_element,"num_grid_quad_element","crate::sim::models_linear_elasticity::element_stiffness_matrix","Computes the plane-stress stiffness matrix of one bilinear quadrilateral. Supply four counterclockwise [x,y] rows scaled by the mechanical Length reference, Young's modulus scaled by the stress reference, and a dimensionless Poisson ratio.\n\nThe record contains rows, cols, the 64 normalized entries in row-major order, and an entry_reference in Force per Length. Model thickness is one Length reference.\n\nInvalid material, malformed vertices, and nonpositive or nonfinite quadrature Jacobians produce errors.",true,(refs:Handle<References>,coordinates:Vec<Vec<f64>>,young:f64,poisson:f64)->Dynamic,|ctx:&CallCtx<'_>,refs:Handle<References>,coordinates:Vec<Vec<f64>>,young:f64,poisson:f64|{let p=points(coordinates)?;if p.len()!=4{return Err(err("four quadrilateral vertices required"));}let m=quad::checked_stiffness([p[0],p[1],p[2],p[3]],young,poisson).map_err(err)?;refs.with(ctx,|r|{let(s,d)=r.derived([0,-2,1,0])?;Ok(Dynamic(api::record(vec![("rows",8usize.encode()?),("cols",8usize.encode()?),("entries",m.into_iter().flatten().collect::<Vec<_>>().encode()?),("entry_reference",api::quantity(s,d)?)])))})});
object_binding!(
    quad_solve, "num_grid_quad_solve", "crate::sim::models_linear_elasticity::run_elasticity_simulation",
    "Solves a bilinear quadrilateral plane-stress mesh. Supply [x,y] node rows and four zero-based node indices per counterclockwise element, followed by Young's modulus, the dimensionless Poisson ratio, fixed node indices and [node,fx,fy] load rows.\n\nNormalize coordinates, Young's modulus and loads by the mechanical Length, stress and Force references. Thickness is one Length reference. Both displacement components of a fixed node are zero. The return record contains interleaved normalized x/y displacements and their measured Length reference.\n\nThe adapter bounds node count, load rows and estimated solve work. It rejects invalid connectivity, geometry, material and nonconvergent solves.",
    true, (refs:Handle<References>,coordinates:Vec<Vec<f64>>,elements:Vec<Vec<f64>>,young:f64,poisson:f64,fixed:Vec<usize>,loads:Vec<Vec<f64>>) -> Dynamic,
    |ctx:&CallCtx<'_>,refs:Handle<References>,coordinates:Vec<Vec<f64>>,elements:Vec<Vec<f64>>,young:f64,poisson:f64,fixed:Vec<usize>,loads:Vec<Vec<f64>>| {
        const MAX_NODES: usize = 256;
        let nodes = points(coordinates)?;
        if nodes.is_empty() || nodes.len() > MAX_NODES || loads.len() > 512
            || elements.len().checked_mul(64 * 5000).is_none_or(|n| n > 50_000_000)
        {
            return Err(err("quadrilateral solve work/allocation cap exceeded"));
        }
        let mut constrained = [false; MAX_NODES];
        for node in fixed {
            if node >= nodes.len() {
                return Err(err("quadrilateral fixed node index out of range"));
            }
            constrained[node] = true;
        }
        // Repeated constraints fix the same degrees of freedom and must not
        // multiply the number of sparse entries reserved for assembly.
        let fixed_nodes = constrained[..nodes.len()]
            .iter()
            .enumerate()
            .filter_map(|(node, &fixed)| fixed.then_some(node))
            .collect();
        let elements = elements.into_iter().map(|e| {
            if e.len() != 4 || e.iter().any(|i| *i < 0.0 || i.fract() != 0.0 || *i >= nodes.len() as f64) {
                return Err(err("quadrilateral connectivity invalid"));
            }
            Ok([e[0] as usize, e[1] as usize, e[2] as usize, e[3] as usize])
        }).collect::<Result<_, OpError>>()?;
        let loads = loads.into_iter().map(|v| {
            if v.len() != 3 || v[0] < 0.0 || v[0].fract() != 0.0 || v[0] >= nodes.len() as f64 {
                return Err(err("load row needs [node,fx,fy]"));
            }
            Ok((v[0] as usize, v[1], v[2]))
        }).collect::<Result<_, OpError>>()?;
        refs.with(ctx, |r| {
            let params = quad::ElasticityParameters {
                nodes, elements, youngs_modulus: young, poissons_ratio: poisson, fixed_nodes, loads,
            };
            output(r, quad::run_elasticity_simulation(&params).map_err(err)?)
        })
    }
);
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
