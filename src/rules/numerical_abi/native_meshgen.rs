//! Measured planar meshing with portable triangle tensors and work bounds.
use crate::graph::RuleError;
use crate::graph::rule::Installer;
use crate::kernel::OpError;
use crate::kernel::Value;
use crate::kernel::signature::DimExpr;
use crate::kernel::signature::Signature;
use crate::kernel::signature::TypeExpr;
use crate::kernel::value::TensorValue;
use crate::kernels::meshgen::TriMesh2;
use crate::kernels::meshgen::advancing_front::advancing_front_bounded;
use crate::kernels::meshgen::delaunay::constrained_delaunay_bounded;
use crate::rules::numerical_abi::codecs::Arg;
use crate::rules::numerical_abi::codecs::Output;
use crate::rules::numerical_abi::native_object_api::Binding;
use crate::rules::numerical_abi::native_object_api::Dynamic;
use crate::rules::numerical_abi::native_object_api::Input;
use crate::rules::numerical_abi::native_object_api::object_binding;
use crate::rules::numerical_abi::native_object_api::{
    self as object,
};

const MAX_VERTICES: usize = 4096;
const MAX_CONSTRAINTS: usize = 512;

fn matrix_rows(
    input: &Value,
    index: usize,
    limit: usize,
) -> Result<usize, OpError> {
    let Value::Tensor(tensor) = input else {
        return Err(object::error("expected a rank-two mesh tensor"));
    };
    let [rows, 2] = tensor.shape.as_slice() else {
        return Err(object::error("mesh input tensor must have two columns"));
    };
    if *rows > limit || rows.checked_mul(2) != Some(tensor.data.len()) {
        return Err(if *rows > limit {
            OpError::BudgetExceeded
        } else {
            object::error(format!("mesh tensor shape is invalid at argument {index}"))
        });
    }
    Ok(*rows)
}

fn budget(limit: usize) -> Result<usize, OpError> {
    if !(3..=MAX_VERTICES).contains(&limit) {
        return Err(OpError::BudgetExceeded);
    }
    Ok(limit)
}

struct NormalizedMesh {
    points: Vec<[f64; 2]>,
    size: f64,
    origin: [f64; 2],
    dimension: DimExpr,
}

fn normalized_points(
    input: &Value,
    reference: &Value,
    limit: usize,
) -> Result<NormalizedMesh, OpError> {
    let rows = matrix_rows(input, 0, limit)?;
    let (size, dimension) = object::reference(reference)?;
    object::same_dimension(&object::dimension(input)?, &dimension)?;
    let mut points = Vec::<[f64; 2]>::decode(input, 0)?;
    if rows < 3 {
        return Err(object::error("a planar mesh needs at least three vertices"));
    }
    let origin = *points
        .first()
        .ok_or_else(|| object::error("missing coordinate origin"))?;
    for point in &mut points {
        for (value, offset) in point.iter_mut().zip(origin) {
            let delta = *value - offset;
            *value = if delta.is_finite() {
                delta / size
            } else {
                *value / size - offset / size
            };
            if !value.is_finite() {
                return Err(object::error(
                    "mesh coordinate exceeds normalized Float range",
                ));
            }
        }
    }
    Ok(NormalizedMesh {
        points,
        size,
        origin,
        dimension,
    })
}

fn encode_mesh(
    mesh: TriMesh2,
    size: f64,
    origin: [f64; 2],
    dimension: DimExpr,
) -> Result<Dynamic, OpError> {
    let rows = mesh.points.len();
    let points = mesh
        .points
        .into_iter()
        .flat_map(|row| {
            [
                row[0].mul_add(size, origin[0]),
                row[1].mul_add(size, origin[1]),
            ]
        })
        .map(|value| {
            if value.is_finite() {
                Ok(Value::Float(value))
            } else {
                Err(object::error("mesh coordinates exceed Float range"))
            }
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok(Dynamic(object::record(vec![
        (
            "points",
            Value::Tensor(TensorValue {
                shape: vec![rows, 2],
                data: points,
                elem: TypeExpr::Float.with_dim(dimension),
            }),
        ),
        ("triangles", mesh.triangles.encode()?),
    ])))
}

fn mesher_error(message: String) -> OpError {
    if crate::cancellation::is_cancelled() {
        return crate::cancellation::Cancelled.into();
    }
    if message.ends_with("vertex budget exceeded") {
        OpError::BudgetExceeded
    } else {
        OpError::Domain(message)
    }
}

object_binding!(
    advancing,
    "num_meshgen_advancing_front_2d",
    "crate::kernels::meshgen::advancing_front::advancing_front_2d",
    "Builds a counterclockwise triangle mesh of a simple planar boundary. Pass two coordinate columns, positive target spacing, positive coordinate reference and maximum vertices. Both length-like arguments share the boundary unit.",
    false,
    (boundary: Input, spacing: Input, reference: Input, max_vertices: usize) -> Dynamic,
    |_, boundary: Input, spacing: Input, reference: Input, max_vertices: usize| {
        let limit = budget(max_vertices)?;
        let NormalizedMesh { points, size, origin, dimension } =
            normalized_points(&boundary.0, &reference.0, limit)?;
        let (spacing,spacing_dim)=object::reference(&spacing.0)?;
        object::same_dimension(&spacing_dim,&dimension)?;
        validate_loop(&points)?;
        let mesh=advancing_front_bounded(&points,spacing/size,limit).map_err(mesher_error)?;
        encode_mesh(mesh,size,origin,dimension)
    }
);

object_binding!(
    constrained,
    "num_meshgen_constrained_delaunay_2d",
    "crate::kernels::meshgen::delaunay::constrained_delaunay_2d_mesh",
    "Triangulates measured points while preserving zero-based constrained segment endpoints as chains of edges. Crossing constraints can create additional measured vertices up to the supplied vertex budget.",
    false,
    (points: Input, segments: Input, reference: Input, max_vertices: usize) -> Dynamic,
    |_, points: Input, segments: Input, reference: Input, max_vertices: usize| {
        let limit = budget(max_vertices)?;
        let NormalizedMesh { points, size, origin, dimension } =
            normalized_points(&points.0, &reference.0, limit)?;
        matrix_rows(&segments.0,1,MAX_CONSTRAINTS)?;
        let segments=Vec::<[usize;2]>::decode(&segments.0,1)?;
        let edges:Vec<_>=segments.into_iter().map(|[a,b]|(a,b)).collect();
        let mesh=constrained_delaunay_bounded(&points,&edges,limit).map_err(mesher_error)?;
        encode_mesh(mesh,size,origin,dimension)
    }
);

#[allow(
    clippy::single_call_fn,
    reason = "keep geometric boundary validation separate from the binding"
)]
fn validate_loop(points: &[[f64; 2]]) -> Result<(), OpError> {
    use crate::kernels::geometry::segment_intersection_2d;
    let length = if points.first() == points.last() {
        points.len().saturating_sub(1)
    } else {
        points.len()
    };
    if length < 3 {
        return Err(object::error(
            "advancing front needs at least three distinct boundary vertices",
        ));
    }
    let points = points
        .get(..length)
        .ok_or_else(|| object::error("invalid closed boundary length"))?;
    let edges = points
        .iter()
        .copied()
        .zip(points.iter().copied().cycle().skip(1));
    for (i, (a, b)) in edges.clone().enumerate() {
        crate::cancellation::checkpoint()?;
        for (j, (c, d)) in edges.clone().enumerate().skip(i.saturating_add(2)) {
            if i == 0 && j.saturating_add(1) == length {
                continue;
            }
            if segment_intersection_2d(a, b, c, d).is_some() {
                return Err(object::error("advancing-front boundary edges intersect"));
            }
        }
    }
    Ok(())
}

fn mesh_signature(second: TypeExpr) -> Signature {
    let coordinate = TypeExpr::Float.with_dim(DimExpr::Var("D".into()));
    let points = object::tensor_type(2, coordinate.clone());
    let indices = object::tensor_type(2, TypeExpr::Int);
    Signature::new(
        vec![points.clone(), second, coordinate, TypeExpr::Int],
        TypeExpr::Record(vec![
            ("points".into(), points),
            ("triangles".into(), indices),
        ]),
    )
}

static ADVANCING: Binding = Binding {
    signature: || mesh_signature(TypeExpr::Float.with_dim(DimExpr::Var("D".into()))),
    ..advancing::BINDING
};
static CONSTRAINED: Binding = Binding {
    signature: || mesh_signature(object::tensor_type(2, TypeExpr::Int)),
    ..constrained::BINDING
};

#[allow(clippy::single_call_fn, reason = "numerical ABI family installer")]
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    object::install(i, &[&ADVANCING, &CONSTRAINED], "numerical_vector")
}
