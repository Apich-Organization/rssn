//! Session-local multilinear interpolation with coherent physical axes and values.

use crate::graph::RuleError;
use crate::graph::rule::Installer;
use crate::kernel::CallCtx;
use crate::kernel::OpError;
use crate::kernel::Value;
use crate::kernel::signature::DimExpr;
use crate::kernel::signature::TypeExpr;
use crate::kernel::value::TensorValue;
use crate::kernels::interpolate::OwnedMultilinearGrid;
use crate::rules::numerical_abi::native_object_api::Dynamic;
use crate::rules::numerical_abi::native_object_api::Handle;
use crate::rules::numerical_abi::native_object_api::Input;
use crate::rules::numerical_abi::native_object_api::MeasuredVector;
use crate::rules::numerical_abi::native_object_api::NativeArg;
use crate::rules::numerical_abi::native_object_api::NativeType;
use crate::rules::numerical_abi::native_object_api::NewObject;
use crate::rules::numerical_abi::native_object_api::object_binding;
use crate::rules::numerical_abi::native_object_api::{
    self as object,
};
use crate::rules::numerical_abi::numeric;

const VALUE_LIMIT: usize = object::VALUE_LIMIT;
// Each axis needs at least two nodes and the value vector has a fixed wire budget.
const MAX_AXES: usize = VALUE_LIMIT.ilog2() as usize;

impl NativeType for OwnedMultilinearGrid {
    const KIND: &'static str = "rssn/native/interpolation/multilinear-grid/v1";
    const MAX_BYTES: usize = 8_388_608;

    fn validate(&self) -> Result<(), OpError> {
        // Both arrays are immutable and were validated together at construction.
        Ok(())
    }

    fn retained_bytes(&self) -> usize {
        self.allocated_bytes()
    }
}

fn rank_one(mut value: &Value) -> Result<&TensorValue, OpError> {
    while let Value::Dyn(dynamic) = value {
        value = &dynamic.value;
    }
    let Value::Tensor(tensor) = value else {
        return Err(object::error(
            "multilinear grid inputs must be rank-one tensors",
        ));
    };
    if tensor.shape != [tensor.data.len()] {
        return Err(object::error(
            "multilinear grid tensor shape must match its rank-one data",
        ));
    }
    Ok(tensor)
}

fn check_structural_envelope(mut value: &Value) -> Result<(), OpError> {
    while let Value::Dyn(dynamic) = value {
        object::same_dimension(&dynamic.dim, &DimExpr::Dimless)?;
        value = &dynamic.value;
    }
    Ok(())
}

struct DecodedAxes {
    nodes: Vec<Vec<f64>>,
    dimensions: Vec<DimExpr>,
}

impl NativeArg for DecodedAxes {
    fn ty() -> TypeExpr {
        object::dynamic_type()
    }

    fn decode(
        value: &Value,
        arg: usize,
        _: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        check_structural_envelope(value)?;
        let axes = rank_one(value)?;
        if !(1..=MAX_AXES).contains(&axes.data.len())
            || matches!(&axes.elem, TypeExpr::Dim(..) | TypeExpr::Point(_))
        {
            return Err(object::error(
                "invalid multidimensional axis collection or grid budget",
            ));
        }
        let mut nodes = Vec::with_capacity(axes.data.len());
        let mut dimensions = Vec::with_capacity(axes.data.len().saturating_add(1));
        for axis in &axes.data {
            let tensor = rank_one(axis)?;
            dimensions.push(object::dimension(axis)?);
            let mut coordinates = Vec::with_capacity(tensor.data.len());
            for coordinate in &tensor.data {
                coordinates.push(numeric(coordinate, arg)?);
            }
            nodes.push(coordinates);
        }
        Ok(Self { nodes, dimensions })
    }
}

#[allow(
    clippy::single_call_fn,
    reason = "keeps heterogeneous and homogeneous query decoding under one physical contract"
)]
fn decode_query(
    value: &Value,
    axes: &[DimExpr],
    output: &mut [f64],
) -> Result<(), OpError> {
    let point = rank_one(value)?;
    if point.data.len() != axes.len() {
        return Err(object::error("the query needs one coordinate per axis"));
    }
    let shared_dimension = if matches!(&point.elem, TypeExpr::Dyn(_)) {
        check_structural_envelope(value)?;
        None
    } else {
        Some(object::dimension(value)?)
    };
    for ((coordinate, expected), slot) in point.data.iter().zip(axes).zip(output) {
        if let Some(dimension) = &shared_dimension {
            object::same_dimension(dimension, expected)?;
        } else {
            object::same_dimension(&object::dimension(coordinate)?, expected)?;
        }
        *slot = numeric(coordinate, 1)?;
    }
    Ok(())
}

// The session owns the arrays. The borrowed Rust API uses the same grid
// validation and interpolation routines without taking ownership.
object_binding!(
    new,
    "num_native_multilinear_grid_new",
    "crate::kernels::interpolate::MultilinearGrid<'grid>::new",
    "Create an immutable session-local rectilinear grid. Axes is an ordered rank-one tensor of measured rank-one axis tensors, with independent units and possibly different lengths. Values is a measured rank-one tensor laid out with the first axis varying fastest. Axes must contain at least two finite increasing nodes each. Shape, unit, and allocation bounds are checked at construction.",
    false,
    (axes: DecodedAxes, values: MeasuredVector<{ VALUE_LIMIT }>) -> NewObject<OwnedMultilinearGrid>,
    |_ctx, mut axes: DecodedAxes, values: MeasuredVector<{ VALUE_LIMIT }>| {
        axes.dimensions.push(values.dimension);
        let grid = OwnedMultilinearGrid::new(axes.nodes, values.values).map_err(object::error)?;
        Ok(NewObject::new(grid, axes.dimensions))
    }
);

object_binding!(
    evaluate,
    "num_native_multilinear_grid_eval",
    "crate::kernels::interpolate::MultilinearGrid<'grid>::evaluate",
    "Interpolate a rank-one point in axis order using an existing session-local grid. A Dynamic-element point carries each coordinate's dimension separately. A homogeneous measured point supplies one shared dimension. Every coordinate must match its axis and be inside the grid. Returns a value in the grid's stored ordinate dimension.",
    false,
    (grid: Handle<OwnedMultilinearGrid>, point: Input) -> Dynamic,
    |ctx, grid: Handle<OwnedMultilinearGrid>, point: Input| {
        let dimensions = grid.dimensions(ctx)?;
        let (value_dim, axis_dims) = dimensions
            .split_last()
            .ok_or_else(|| object::error("multilinear grid is missing its stored dimensions"))?;
        let mut coordinates = [0.0; MAX_AXES];
        let coordinates = coordinates
            .get_mut(..axis_dims.len())
            .ok_or_else(|| object::error("multilinear grid exceeds the axis budget"))?;
        decode_query(&point.0, axis_dims, coordinates)?;
        let value = grid.with(ctx, |grid| grid.evaluate(coordinates).map_err(object::error))?;
        Ok(Dynamic(object::quantity(value, value_dim.clone())?))
    }
);

pub(super) fn install(installer: &mut Installer<'_>) -> Result<(), RuleError> {
    object::install(
        installer,
        &[&new::BINDING, &evaluate::BINDING],
        "numerical_vector",
    )
}
