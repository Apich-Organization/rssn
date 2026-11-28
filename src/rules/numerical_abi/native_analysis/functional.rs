//! Prepared norms for fields sampled on a measured axis.
use super::Installer;
use super::RuleError;
use super::native_object_api as api;
use super::output_scalar;
use crate::kernel::CallCtx;
use crate::kernel::OpError;
use crate::kernel::Value;
use crate::kernel::signature::DimExpr;
use crate::kernel::signature::Exp;
use crate::kernel::signature::Rational;
use crate::kernel::signature::TypeExpr;
use crate::kernel::value::TensorValue;
use crate::kernels::functional_analysis::L2Grid;
use api::Dynamic;
use api::Handle;
use api::Maybe;
use api::MeasuredVector;
use api::NativeOutput;
use api::NativeType;
use api::NewObject;
use api::object_binding;

const MAX_POINTS: usize = 65_536;

impl NativeType for L2Grid {
    const KIND: &'static str = "rssn/native/analysis/l2-grid/v1";
    const MAX_BYTES: usize = MAX_POINTS * size_of::<f64>() + 8192;

    fn validate(&self) -> Result<(), OpError> {
        Ok(())
    }

    fn retained_bytes(&self) -> usize {
        size_of::<Self>() + self.heap_bytes()
    }
}

type Samples = MeasuredVector<MAX_POINTS>;

fn root_axis_dimension(
    ctx: &CallCtx<'_>,
    grid: &Handle<L2Grid>,
) -> Result<DimExpr, OpError> {
    let dimensions = grid.dimensions(ctx)?;
    let [axis] = dimensions.as_slice() else {
        return Err(api::error("L2 grid has invalid axis dimension metadata"));
    };
    Ok(axis.clone().pow(Exp::Rat(Rational::new(1, 2))))
}

object_binding!(
    l2_grid_new,
    "num_functional_l2_grid_new",
    "crate::kernels::functional_analysis::L2Grid::new",
    "Prepares trapezoidal L2 weights for finite nondecreasing coordinates. The result can be reused for independent sampled fields on the same axis. Repeated coordinates contribute zero panel width.",
    false,
    (coordinates: Samples) -> NewObject<L2Grid>,
    |_, coordinates| {
        let grid = L2Grid::new(&coordinates.values).map_err(api::error)?;
        Ok(NewObject::new(grid, vec![coordinates.dimension]))
    }
);

object_binding!(
    l2_grid_norm,
    "num_functional_l2_grid_norm",
    "crate::kernels::functional_analysis::L2Grid::norm",
    "Returns the sampled L2 norm, with dimension equal to the ordinates' dimension times the square root of the axis dimension. Returns a domain error if the norm overflows Float.",
    false,
    (grid: Handle<L2Grid>, ordinates: Samples) -> Dynamic,
    |ctx, grid, ordinates| {
        let norm_dimension = DimExpr::Mul(vec![ordinates.dimension, root_axis_dimension(ctx, &grid)?]);
        grid.with(ctx, |grid| {
            let norm = grid.norm(&ordinates.values).map_err(api::error)?;
            Ok(Dynamic(api::quantity(norm, norm_dimension)?))
        })
    }
);

object_binding!(
    l2_grid_norm_and_normalize,
    "num_functional_l2_grid_norm_and_normalize",
    "crate::kernels::functional_analysis::L2Grid::norm_and_normalize",
    "Returns {norm, normalized} from one reduction. The measured norm is None if it overflows Float. Normalized ordinates carry the inverse square root of the axis dimension and are None for zero measure. Normalization uses the unrounded scale, including when the norm underflows or overflows.",
    false,
    (grid: Handle<L2Grid>, ordinates: Samples) -> Dynamic,
    |ctx, grid, ordinates| {
        let axis_root = root_axis_dimension(ctx, &grid)?;
        let norm_dimension = DimExpr::Mul(vec![ordinates.dimension, axis_root.clone()]);
        grid.with(ctx, |grid| {
            let (norm, normalized, has_measure) = grid.normalize_with_measure(&ordinates.values).map_err(api::error)?;
            let norm = norm.is_finite().then(|| api::quantity(norm, norm_dimension)).transpose()?;
            let normalized = if has_measure {
                let data = normalized.into_iter().map(output_scalar).collect::<Result<Vec<_>, _>>()?;
                Some(Value::Tensor(TensorValue {
                    shape: vec![data.len()],
                    data,
                    elem: TypeExpr::Float.with_dim(axis_root.pow(-1)),
                }))
            } else {
                None
            };
            Ok(Dynamic(api::record(vec![
                ("norm", Maybe(norm.map(Dynamic)).encode(ctx)?),
                ("normalized", Maybe(normalized.map(Dynamic)).encode(ctx)?),
            ])))
        })
    }
);

pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &l2_grid_new::BINDING,
            &l2_grid_norm::BINDING,
            &l2_grid_norm_and_normalize::BINDING,
        ],
        "numerical_vector",
    )
}
