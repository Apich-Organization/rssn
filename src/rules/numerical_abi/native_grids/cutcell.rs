//! Native level-set clipping with explicit geometric references and callback limits.
use super::CallCtx;
use super::Frame;
use super::Handle;
use super::Input;
use super::Installer;
use super::OpError;
use super::RuleError;
use super::api;
use super::err;
use super::object_binding;
use super::structured_ops::checked_snapshot;
use crate::sim::cutcell;

fn clip(
    ctx: &CallCtx<'_>,
    frame: &Frame,
    function: &Input,
    budget: usize,
    fraction: impl FnOnce(cutcell::LevelSet<'_>) -> f64,
) -> Result<f64, OpError> {
    let length = &frame.dims[0];
    let coordinates = [length.clone(), length.clone(), length.clone()];
    let distance = checked_snapshot(ctx, function, budget, &coordinates, length)?;
    // Clipping uses signs and ratios of distances. Retain the coherent samples
    // so that rescaling by the geometric reference cannot erase a sign.
    let result = fraction(&|point| distance.sample(&point.map(|x| x * frame.length)));
    distance.finish()?;
    Ok(result)
}

object_binding!(cut_simplex,"num_grid_simplex_positive_fraction","crate::sim::cutcell::simplex_positive_fraction","Native positive fraction of a segment/triangle/tetrahedron in dimensionless level-set coordinates; common positive rescaling does not change the fraction.",true,(values:Vec<f64>)->f64,|ctx,values|{let _=ctx;if !(2..=4).contains(&values.len()){return Err(err("clipping simplex needs 2..4 finite level-set samples"));}Ok(cutcell::simplex_positive_fraction(&values))});
object_binding!(
    cut_segment, "num_grid_open_length_fraction", "crate::sim::cutcell::open_length_fraction",
    "Returns the dimensionless open fraction of an axis-aligned segment. Axis indices are zero based. Normalize the three origin coordinates and positive segment length by the frame's Length reference.\n\nThe compiled field takes three coherent Length coordinates and returns Length. Positive distances are open. The kernel samples the segment, then bisects intervals where its sampled sign changes.\n\nThe callback budget accepts 1 through 8,388,608 evaluations, including refinement. Invalid geometry, callback units or arity, nonfinite samples, and exhausted evaluation budgets produce errors.",
    false, (frame:Handle<Frame>,function:Input,origin:Vec<f64>,axis:usize,length:f64,budget:usize) -> f64,
    |ctx,frame,function,origin,axis,length,budget| frame.with(ctx, |f| {
        if axis >= 3 || length <= 0.0 || !(1..=8_388_608).contains(&budget) {
            return Err(err("segment geometry/callback budget invalid"));
        }
        let origin = origin.try_into().map_err(|_| err("segment origin needs three coordinates"))?;
        clip(ctx, f, &function, budget, |distance| {
            cutcell::open_length_fraction(distance, origin, axis, length)
        })
    })
);
object_binding!(
    cut_box, "num_grid_open_fraction", "crate::sim::cutcell::open_fraction",
    "Returns the dimensionless open fraction along one, two, or three distinct zero-based axes. Normalize the three origin and length slots by the frame's Length reference and supply positive lengths on the active axes.\n\nThe compiled field takes three coherent Length coordinates and returns Length. Positive signed distances are open. One active axis uses segment root refinement. Two or three active axes use adaptive subdivision and simplex clipping.\n\nThe depth cap is 12 and depth times the number of axes is at most 18. The independent callback budget accepts 1 through 8,388,608 evaluations. Invalid geometry, callback units or arity, nonfinite samples and budget exhaustion produce errors.",
    false, (frame:Handle<Frame>,function:Input,origin:Vec<f64>,lengths:Vec<f64>,axes:Vec<usize>,depth:usize,budget:usize) -> f64,
    |ctx,frame,function,origin,lengths,axes,depth,budget| frame.with(ctx, |f| {
        if !(1..=3).contains(&axes.len()) || axes.iter().any(|a| *a >= 3)
            || axes.iter().enumerate().any(|(i, a)| axes[i + 1..].contains(a))
            || depth > 12 || depth * axes.len() > 18 || !(1..=8_388_608).contains(&budget)
        {
            return Err(err("clipping axis/subdivision/callback work budget invalid"));
        }
        let origin = origin.try_into().map_err(|_| err("box origin needs three coordinates"))?;
        let lengths: [f64; 3] = lengths.try_into().map_err(|_| err("box lengths need three slots"))?;
        if axes.iter().any(|a| lengths[*a] <= 0.0) {
            return Err(err("active clipping lengths must be positive"));
        }
        clip(ctx, f, &function, budget, |distance| {
            cutcell::open_fraction(distance, origin, lengths, &axes, depth)
        })
    })
);
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
