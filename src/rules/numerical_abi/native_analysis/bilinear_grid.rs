//! Session-owned, measured rectilinear grids for repeated interpolation.

use crate::graph::RuleError;
use crate::graph::dim::DimRule;
use crate::graph::rule::Installer;
use crate::kernel::CallCtx;
use crate::kernel::OpError;
use crate::kernel::Value;
use crate::kernel::signature::DimExpr;
use crate::kernel::signature::Signature;
use crate::kernel::signature::TypeExpr;
use crate::kernels::interp_ext::BilinearGridData;
use crate::rules::numerical_abi::native_object_api::Binding;
use crate::rules::numerical_abi::native_object_api::Dynamic;
use crate::rules::numerical_abi::native_object_api::Handle;
use crate::rules::numerical_abi::native_object_api::Input;
use crate::rules::numerical_abi::native_object_api::NativeType;
use crate::rules::numerical_abi::native_object_api::NewObject;
use crate::rules::numerical_abi::native_object_api::object_binding;
use crate::rules::numerical_abi::native_object_api::{
    self as object,
};
use crate::rules::numerical_abi::numeric;

struct PreparedGrid(BilinearGridData<Vec<f64>>);

impl NativeType for PreparedGrid {
    const KIND: &'static str = "rssn/native/interpolation/bilinear/v1";
    const MAX_BYTES: usize = 3 * object::VALUE_LIMIT * size_of::<f64>() + size_of::<Self>();

    fn validate(&self) -> Result<(), OpError> {
        Ok(())
    }

    fn retained_bytes(&self) -> usize {
        size_of::<Self>() + self.0.capacity() * size_of::<f64>()
    }
}

fn measured_tensor(
    value: &Value,
    shape: &[usize],
    argument: usize,
) -> Result<(Vec<f64>, DimExpr), OpError> {
    let Value::Tensor(tensor) = value else {
        return Err(object::error("expected a measured grid tensor"));
    };
    if tensor.shape != shape {
        return Err(object::error(
            "grid shape must be [y axis length, x axis length]",
        ));
    }
    let units = object::dimension(value)?;
    let samples = tensor
        .data
        .iter()
        .map(|value| numeric(value, argument))
        .collect::<Result<_, _>>()?;
    Ok((samples, units))
}

fn axis(
    value: &Value,
    argument: usize,
) -> Result<(Vec<f64>, DimExpr), OpError> {
    let Value::Tensor(tensor) = value else {
        return Err(object::error("expected a measured grid axis"));
    };
    measured_tensor(value, &[tensor.data.len()], argument)
}

object_binding!(
    prepare,
    "num_interp_bilinear_grid_prepare",
    "crate::kernels::interp_ext::BilinearGrid<'grid>::new",
    "Prepares a rectilinear grid with measured axes. The value tensor is [y axis length, x axis length] in row-major order. Each axis and the values may have independent units.",
    false,
    (x_axis: Input, y_axis: Input, grid: Input) -> NewObject<PreparedGrid>,
    |_, x_axis, y_axis, grid| {
        let (xs, xd) = axis(&x_axis.0, 0)?;
        let (ys, yd) = axis(&y_axis.0, 1)?;
        let (values, vd) = measured_tensor(&grid.0, &[ys.len(), xs.len()], 2)?;
        let grid = BilinearGridData::new(xs, ys, values).map_err(|e| object::error(e.to_string()))?;
        Ok(NewObject::new(PreparedGrid(grid), vec![xd, yd, vd]))
    }
);

object_binding!(
    evaluate,
    "num_interp_bilinear_grid_evaluate",
    "crate::kernels::interp_ext::BilinearGrid<'grid>::eval",
    "Interpolates one finite point on a prepared rectilinear grid. Queries are clamped at the edges. The supplied value witness sets the return unit and must match the grid cells.",
    false,
    (grid: Handle<PreparedGrid>, x: Input, y: Input, witness: Input) -> Dynamic,
    |ctx: &CallCtx<'_>, grid: Handle<PreparedGrid>, x: Input, y: Input, witness: Input| {
        let dimensions = grid.dimensions(ctx)?;
        let [x_dim, y_dim, value_dim] = dimensions.as_slice() else {
            return Err(object::error("bilinear grid has incorrect dimension metadata"));
        };
        for (input, expected) in [(&x, x_dim), (&y, y_dim), (&witness, value_dim)] {
            object::same_dimension(&object::dimension(&input.0)?, expected)?;
        }
        let x = numeric(&x.0, 1)?;
        let y = numeric(&y.0, 2)?;
        numeric(&witness.0, 3)?;
        let sample = grid.with(ctx, |model| model.0.eval(x, y).map_err(|e| object::error(e.to_string())))?;
        object::quantity(sample, value_dim.clone()).map(Dynamic)
    }
);

fn prepare_signature() -> Signature {
    Signature::new(
        vec![
            object::tensor_type(1, TypeExpr::Float.with_dim(DimExpr::Var("X".into()))),
            object::tensor_type(1, TypeExpr::Float.with_dim(DimExpr::Var("Y".into()))),
            object::tensor_type(2, TypeExpr::Float.with_dim(DimExpr::Var("Cell".into()))),
        ],
        TypeExpr::Opaque(PreparedGrid::KIND.into()),
    )
}

fn evaluate_signature() -> Signature {
    let cell = TypeExpr::Float.with_dim(DimExpr::Var("Cell".into()));
    Signature::new(
        vec![
            TypeExpr::Opaque(PreparedGrid::KIND.into()),
            TypeExpr::Float.with_dim(DimExpr::Var("X".into())),
            TypeExpr::Float.with_dim(DimExpr::Var("Y".into())),
            cell.clone(),
        ],
        cell,
    )
}

static PREPARE: Binding = Binding {
    signature: prepare_signature,
    ..prepare::BINDING
};
static EVALUATE: Binding = Binding {
    signature: evaluate_signature,
    ..evaluate::BINDING
};

pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    object::install(i, &[&PREPARE], "numerical_vector")?;
    object::install_with_dimension(i, &[&EVALUATE], "numerical_vector", DimRule::Same(3))
}
