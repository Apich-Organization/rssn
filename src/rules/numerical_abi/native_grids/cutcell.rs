//! Native level-set clipping with explicit geometric references and callback limits.
use super::Frame;
use super::Handle;
use super::Input;
use super::Installer;
use super::RuleError;
use super::api;
use super::err;
use super::object_binding;
use crate::sim::cutcell;
object_binding!(cut_simplex,"num_grid_simplex_positive_fraction","crate::sim::cutcell::simplex_positive_fraction","Native positive fraction of a segment/triangle/tetrahedron in dimensionless level-set coordinates; common positive rescaling does not change the fraction.",true,(values:Vec<f64>)->f64,|ctx,values|{let _=ctx;if !(2..=4).contains(&values.len()){return Err(err("clipping simplex needs 2..4 finite level-set samples"));}Ok(cutcell::simplex_positive_fraction(&values))});
object_binding!(cut_segment,"num_grid_open_length_fraction","crate::sim::cutcell::open_length_fraction","Actual native level-set segment clipping. Coordinates/length are normalized by the explicit Length reference; compiled signed distance takes coherent coordinates and returns Length. Root/bisection callback errors and budgets are preserved.",false,(frame:Handle<Frame>,function:Input,origin:Vec<f64>,axis:usize,length:f64,budget:usize)->f64,|ctx,frame,function,origin,axis,length,budget|{frame.with(ctx,|f|{if axis>=3||length<=0.0||budget==0||budget>8_388_608{return Err(err("segment geometry/callback budget invalid"));}let origin:[f64;3]=origin.try_into().map_err(|_|err("segment origin needs three coordinates"))?;let c=crate::rules::numerical_abi::native_grids::structured_ops::checked_snapshot(ctx,&function,budget,&vec![f.dims[0].clone();3],&f.dims[0])?;let out=cutcell::open_length_fraction(&|x|c.sample(&x.map(|v|v*f.length))/f.length,origin,axis,length);c.finish()?;Ok(out)})});
object_binding!(cut_box,"num_grid_open_fraction","crate::sim::cutcell::open_fraction","Actual native adaptive level-set box clipping. Coordinates and extents use the explicit Length reference, compiled signed distance retains Length, and axes are distinct. Subdivision and callback work are bounded independently.",false,(frame:Handle<Frame>,function:Input,origin:Vec<f64>,lengths:Vec<f64>,axes:Vec<usize>,depth:usize,budget:usize)->f64,|ctx,frame,function,origin,lengths,axes,depth,budget|{frame.with(ctx,|f|{if !(1..=3).contains(&axes.len())||axes.iter().any(|a|*a>=3)||axes.iter().enumerate().any(|(i,a)|axes[i+1..].contains(a))||depth>12||depth*axes.len()>18||budget==0||budget>8_388_608{return Err(err("clipping axis/subdivision/callback work budget invalid"));}let origin:[f64;3]=origin.try_into().map_err(|_|err("box origin needs three coordinates"))?;let lengths:[f64;3]=lengths.try_into().map_err(|_|err("box lengths need three slots"))?;if axes.iter().any(|a|lengths[*a]<=0.0){return Err(err("active clipping lengths must be positive"));}let c=crate::rules::numerical_abi::native_grids::structured_ops::checked_snapshot(ctx,&function,budget,&vec![f.dims[0].clone();3],&f.dims[0])?;let out=cutcell::open_fraction(&|x|c.sample(&x.map(|v|v*f.length))/f.length,origin,lengths,&axes,depth);c.finish()?;Ok(out)})});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &cut_simplex::BINDING,
            &cut_segment::BINDING,
            &cut_box::BINDING,
        ],
        "numerical_vector",
    )
}
