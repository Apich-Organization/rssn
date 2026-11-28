//! Checked, session-owned triangle meshes and in-memory CAD interchange.
use super::codecs::Arg;
use super::native_object_api as api;
use crate::graph::RuleError;
use crate::graph::dim::DimError;
use crate::graph::dim::DimRule;
use crate::graph::rule::Installer;
use crate::kernel::CallCtx;
use crate::kernel::OpError;
use crate::kernel::Value;
use crate::kernel::signature::DimExpr;
use crate::kernel::signature::TypeExpr;
use crate::kernel::value::TensorValue;
use crate::kernels::cad::csg::BoolOp;
use crate::kernels::cad::csg::boolean_checked;
use crate::kernels::cad::repair::RepairReport;
use crate::kernels::cad::repair::Topology;
use crate::kernels::cad::repair::closed_signed_volume;
use crate::kernels::cad::repair::repair;
use crate::kernels::cad::repair::topology;
use crate::kernels::cad::repair::weld;
use crate::kernels::cad::step::read_planar_step;
use crate::kernels::cad::stl::read_stl;
use crate::kernels::cad::stl::write_stl_ascii;
use crate::kernels::cad::stl::write_stl_binary;
use crate::kernels::cad::types::Mesh;
use crate::kernels::cad::types::indices_valid;
use api::Handle;
use api::Many;
use api::NativeArg;
use api::NativeOutput;
use api::NativeType;
use api::NewObject;
use api::error;
use api::object_binding;
use num_rational::BigRational;
use std::cell::OnceCell;
use std::mem::size_of;

const MAX_MESH_ROWS: usize = 16_384;
const MAX_REPAIR_ROWS: usize = 4_096;
const MAX_STL_BYTES: usize = 1_048_576;
const MAX_CAD_COORDINATE: f64 = 1e60;

fn length() -> DimExpr {
    DimExpr::base("Length")
}
fn volume() -> DimExpr {
    length().pow(3)
}
struct CadMesh {
    mesh: Mesh,
    summary: OnceCell<Topology>,
}

impl NativeType for CadMesh {
    const KIND: &'static str = "rssn/native/cad/mesh/v1";
    const MAX_BYTES: usize =
        size_of::<Self>() + MAX_MESH_ROWS * (size_of::<[f64; 3]>() + size_of::<[usize; 3]>());

    fn validate(&self) -> Result<(), OpError> {
        Ok(())
    }

    fn retained_bytes(&self) -> usize {
        size_of::<Self>()
            .saturating_add(
                self.mesh
                    .vertices
                    .capacity()
                    .saturating_mul(size_of::<[f64; 3]>()),
            )
            .saturating_add(
                self.mesh
                    .triangles
                    .capacity()
                    .saturating_mul(size_of::<[usize; 3]>()),
            )
    }
}

impl CadMesh {
    fn prepare(mesh: Mesh) -> Result<NewObject<Self>, OpError> {
        if mesh.vertices.len().max(mesh.triangles.len()) > MAX_MESH_ROWS {
            return Err(error("CAD mesh exceeds the vertex or triangle budget"));
        }
        if !indices_valid(&mesh)
            || mesh
                .vertices
                .iter()
                .flatten()
                .any(|coordinate| !coordinate.is_finite() || coordinate.abs() > MAX_CAD_COORDINATE)
        {
            return Err(error(
                "CAD coordinates are out of range or a triangle index is invalid",
            ));
        }
        Ok(NewObject::normalized(Self {
            mesh,
            summary: OnceCell::new(),
        }))
    }

    fn topology(&self) -> &Topology {
        self.summary.get_or_init(|| topology(&self.mesh))
    }

    fn bounded_repair(
        &self,
        tolerance: Length,
    ) -> Result<f64, OpError> {
        if self.mesh.vertices.len().max(self.mesh.triangles.len()) > MAX_REPAIR_ROWS {
            return Err(error("mesh exceeds the CAD repair work budget"));
        }
        if !(0.0..=MAX_CAD_COORDINATE).contains(&tolerance.0) {
            return Err(error(
                "welding tolerance must be finite, nonnegative and in range",
            ));
        }
        Ok(tolerance.0)
    }
}

#[derive(Clone, Copy)]
struct Length(f64);
impl Length {
    fn positive(self) -> Result<f64, OpError> {
        if self.0 > 0.0 {
            Ok(self.0)
        } else {
            Err(error("STL coordinate length must be positive"))
        }
    }
}
impl NativeArg for Length {
    fn ty() -> TypeExpr {
        TypeExpr::Float.with_dim(length())
    }

    fn decode(
        value: &Value,
        index: usize,
        _: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        api::same_dimension(&api::dimension(value)?, &length())?;
        super::numeric(value, index).map(Self)
    }
}

struct Volume(f64);
impl NativeOutput for Volume {
    fn ty() -> TypeExpr {
        TypeExpr::Float.with_dim(volume())
    }

    fn encode(
        self,
        _: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        api::quantity(self.0, volume())
    }
}

fn triple_rows(value: &Value) -> Result<&[Value], OpError> {
    let Value::Tensor(tensor) = value else {
        return Err(error("expected a CAD rank-two tensor"));
    };
    let [rows, 3] = tensor.shape.as_slice() else {
        return Err(error("CAD tensors must have three columns"));
    };
    if *rows > MAX_MESH_ROWS || rows.saturating_mul(3) != tensor.data.len() {
        return Err(error(
            "CAD tensor data disagrees with its shape or exceeds its work budget",
        ));
    }
    Ok(&tensor.data)
}

struct Vertices(Vec<[f64; 3]>);
impl NativeArg for Vertices {
    fn ty() -> TypeExpr {
        api::tensor_type(2, TypeExpr::Float.with_dim(length()))
    }

    fn decode(
        value: &Value,
        index: usize,
        _: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        let data = triple_rows(value)?;
        api::same_dimension(&api::dimension(value)?, &length())?;
        data.as_chunks::<3>()
            .0
            .iter()
            .map(|[x, y, z]| {
                Ok([
                    super::numeric(x, index)?,
                    super::numeric(y, index)?,
                    super::numeric(z, index)?,
                ])
            })
            .collect::<Result<_, _>>()
            .map(Self)
    }
}

struct Faces(Vec<[usize; 3]>);
impl Arg for Faces {
    fn ty() -> TypeExpr {
        api::tensor_type(2, TypeExpr::Int)
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        triple_rows(value)?
            .as_chunks::<3>()
            .0
            .iter()
            .map(|[a, b, c]| {
                Ok([
                    <usize as Arg>::decode(a, index)?,
                    <usize as Arg>::decode(b, index)?,
                    <usize as Arg>::decode(c, index)?,
                ])
            })
            .collect::<Result<_, _>>()
            .map(Self)
    }
}

struct MeshArrays(Mesh);
impl NativeOutput for MeshArrays {
    fn ty() -> TypeExpr {
        TypeExpr::Record(vec![
            ("vertices".into(), Vertices::ty()),
            ("triangles".into(), <Faces as Arg>::ty()),
        ])
    }

    fn encode(
        self,
        _: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        let mesh = self.0;
        let len = mesh.vertices.len();
        let vertices = mesh
            .vertices
            .into_iter()
            .flatten()
            .map(|value| api::quantity(value, length()))
            .collect::<Result<_, _>>()?;
        let faces = mesh.triangles.len();
        Ok(api::record(vec![
            (
                "vertices",
                Value::Tensor(TensorValue {
                    shape: vec![len, 3],
                    elem: TypeExpr::Float.with_dim(length()),
                    data: vertices,
                }),
            ),
            (
                "triangles",
                Value::Tensor(TensorValue {
                    shape: vec![faces, 3],
                    elem: TypeExpr::Int,
                    data: mesh
                        .triangles
                        .into_iter()
                        .flatten()
                        .map(Value::int)
                        .collect(),
                }),
            ),
        ]))
    }
}

struct WeldedMesh(NewObject<CadMesh>, usize);
impl NativeOutput for WeldedMesh {
    fn ty() -> TypeExpr {
        TypeExpr::Record(vec![
            ("mesh".into(), NewObject::<CadMesh>::ty()),
            ("merged_vertices".into(), TypeExpr::Int),
        ])
    }

    fn encode(
        self,
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        Ok(api::record(vec![
            ("mesh", self.0.encode(ctx)?),
            ("merged_vertices", Value::int(self.1)),
        ]))
    }
}

struct RepairedMesh(NewObject<CadMesh>, RepairReport);
impl NativeOutput for RepairedMesh {
    fn ty() -> TypeExpr {
        TypeExpr::Record(vec![
            ("mesh".into(), NewObject::<CadMesh>::ty()),
            ("merged_vertices".into(), TypeExpr::Int),
            ("removed_degenerate".into(), TypeExpr::Int),
            ("filled_holes".into(), TypeExpr::Int),
            ("flipped".into(), TypeExpr::Int),
        ])
    }

    fn encode(
        self,
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        Ok(api::record(vec![
            ("mesh", self.0.encode(ctx)?),
            ("merged_vertices", Value::int(self.1.merged_vertices)),
            ("removed_degenerate", Value::int(self.1.removed_degenerate)),
            ("filled_holes", Value::int(self.1.filled_holes)),
            ("flipped", Value::int(self.1.flipped)),
        ]))
    }
}

object_binding!(
    mesh_new, "num_cad_mesh_new", "crate::kernels::cad::types::Mesh::new",
    "Prepares an indexed mesh from a rank-two [x, y, z] tensor of lengths and zero-based triangle indices. Both tensors have three columns.",
    false, (vertices: Vertices, faces: Faces) -> NewObject<CadMesh>,
    |_, vertices, faces| CadMesh::prepare(Mesh::new(vertices.0, faces.0))
);

object_binding!(
    mesh_arrays, "num_cad_mesh_arrays", "",
    "Copies the mesh's coherent length coordinates and zero-based triangle indices. Empty tensors retain three columns.",
    true, (mesh: Handle<CadMesh>) -> MeshArrays,
    |ctx, mesh| mesh.with(ctx, |mesh| Ok(MeshArrays(mesh.mesh.clone())))
);

macro_rules! count {
    ($name:ident, $op:literal, $source:literal, $doc:literal, $method:ident) => {
        object_binding!($name, $op, $source, $doc, true, (mesh: Handle<CadMesh>) -> usize,
            |ctx, mesh| mesh.with(ctx, |mesh| Ok(mesh.mesh.$method())));
    };
}
count!(
    vertex_count,
    "num_cad_mesh_vertex_count",
    "crate::kernels::cad::types::Mesh::vertex_count",
    "Returns the number of indexed vertices.",
    vertex_count
);
count!(
    triangle_count,
    "num_cad_mesh_triangle_count",
    "crate::kernels::cad::types::Mesh::triangle_count",
    "Returns the number of indexed triangle faces.",
    triangle_count
);

macro_rules! statistic {
    ($name:ident, $op:literal, $source:literal, $doc:literal, $field:ident, $ret:ty) => {
        object_binding!($name, $op, $source, $doc, true, (mesh: Handle<CadMesh>) -> $ret,
            |ctx, mesh| mesh.with(ctx, |mesh| Ok(mesh.topology().$field)));
    };
}
statistic!(
    is_closed,
    "num_cad_mesh_is_closed",
    "crate::kernels::cad::repair::is_closed",
    "Checks that every edge of the nonempty indexed surface has exactly two incident triangles.",
    closed,
    bool
);
statistic!(
    is_manifold,
    "num_cad_mesh_is_manifold",
    "crate::kernels::cad::repair::is_manifold",
    "Checks edge-manifoldness. This allows boundary edges and requires each triangle to use three distinct indices.",
    manifold,
    bool
);
statistic!(
    euler,
    "num_cad_mesh_euler_characteristic",
    "crate::kernels::cad::repair::euler_characteristic",
    "Returns V − E + F with isolated indexed vertices included.",
    euler,
    i64
);

object_binding!(
    signed_volume, "num_cad_mesh_signed_volume", "crate::kernels::cad::repair::signed_volume",
    "Returns the signed volume of a closed consistently wound mesh. An inward-wound surface gives negative volume.",
    true, (mesh: Handle<CadMesh>) -> Volume,
    |ctx, mesh| mesh.with(ctx, |mesh| {
        if !mesh.topology().oriented { return Err(error("signed volume requires a closed consistently wound mesh")); }
        Ok(Volume(closed_signed_volume(&mesh.mesh)))
    })
);

object_binding!(
    mesh_weld, "num_cad_mesh_weld", "crate::kernels::cad::repair::weld",
    "Merges nearby vertices using the supplied length tolerance and returns a new mesh and a vertex count. Zero tolerance merges bit-identical coordinates.",
    false, (mesh: Handle<CadMesh>, tolerance: Length) -> WeldedMesh,
    |ctx, mesh, tolerance| mesh.with(ctx, |mesh| {
        let (result, removed) = weld(&mesh.mesh, mesh.bounded_repair(tolerance)?);
        Ok(WeldedMesh(CadMesh::prepare(result)?, removed))
    })
);

object_binding!(
    mesh_repair, "num_cad_mesh_repair", "crate::kernels::cad::repair::repair",
    "Welds, reorients and fills boundary loops where the repair algorithm can resolve them. Returns a new model and counts of the changes.",
    false, (mesh: Handle<CadMesh>, tolerance: Length) -> RepairedMesh,
    |ctx, mesh, tolerance| mesh.with(ctx, |mesh| {
        let (result, report) = repair(&mesh.mesh, mesh.bounded_repair(tolerance)?);
        Ok(RepairedMesh(CadMesh::prepare(result)?, report))
    })
);

impl Arg for BoolOp {
    fn ty() -> TypeExpr {
        TypeExpr::String
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        match <String as Arg>::decode(value, index)?.as_str() {
            | "union" => Ok(Self::Union),
            | "difference" => Ok(Self::Difference),
            | "intersection" => Ok(Self::Intersection),
            | _ => {
                Err(error(
                    "CSG operation must be union, difference or intersection",
                ))
            },
        }
    }
}

object_binding!(
    mesh_boolean, "num_cad_mesh_boolean", "crate::kernels::cad::csg::boolean",
    "Combines outward-wound convex meshes by union, difference or intersection within fixed evaluation limits. Inputs and result share the same coherent length frame.",
    true, (left: Handle<CadMesh>, right: Handle<CadMesh>, operation: BoolOp) -> NewObject<CadMesh>,
    |ctx, left, right, operation| left.with(ctx, |left| right.with(ctx, |right| {
        let result = boolean_checked(&left.mesh, &right.mesh, operation).map_err(|message| {
            crate::cancellation::checkpoint().map_or_else(OpError::from, |()| error(message))
        })?;
        CadMesh::prepare(result)
    }))
);

fn decode_stl(
    bytes: &[u8],
    unit: Length,
) -> Result<NewObject<CadMesh>, OpError> {
    let scale = unit.positive()?;
    let mut mesh = read_stl(bytes).map_err(error)?;
    for value in mesh.vertices.iter_mut().flatten() {
        *value *= scale;
    }
    CadMesh::prepare(mesh)
}

object_binding!(
    stl_read, "num_cad_stl_read", "crate::kernels::cad::stl::read_stl",
    "Reads binary or ASCII STL directly from byte values 0–255. Supply the length of one unit in the file. STL triangles repeat vertex coordinates until welded.",
    false, (bytes: Many<u8>, unit: Length) -> NewObject<CadMesh>,
    |_, bytes, unit| decode_stl(&bytes.0, unit)
);

object_binding!(
    stl_read_text, "num_cad_stl_read_text", "",
    "Reads ASCII STL directly from text using the supplied coordinate unit length.",
    false, (text: String, unit: Length) -> NewObject<CadMesh>,
    |_, text, unit| decode_stl(text.as_bytes(), unit)
);

fn stl_coordinates(
    mesh: &Mesh,
    unit: Length,
    max_coordinate: f64,
) -> Result<Mesh, OpError> {
    let scale = unit.positive()?;
    let mut normalized = mesh.clone();
    for value in normalized.vertices.iter_mut().flatten() {
        *value /= scale;
        if !value.is_finite() || value.abs() > max_coordinate {
            return Err(error(
                "mesh exceeds the coordinate range of this STL encoding",
            ));
        }
    }
    Ok(normalized)
}

object_binding!(
    stl_write_binary, "num_cad_stl_write_binary", "crate::kernels::cad::stl::write_stl_binary",
    "Returns binary STL as byte values 0–255. Coordinates are divided by the supplied unit length before conversion to STL Float32.",
    false, (mesh: Handle<CadMesh>, unit: Length) -> Many<u8>,
    |ctx, mesh, unit| mesh.with(ctx, |mesh| {
        if mesh.mesh.triangle_count().saturating_mul(50).saturating_add(85) > api::VALUE_LIMIT {
            return Err(error("binary STL exceeds the native byte budget"));
        }
        let scaled = stl_coordinates(&mesh.mesh, unit, f64::from(f32::MAX))?;
        Ok(Many(write_stl_binary(&scaled)))
    })
);

object_binding!(
    stl_write_ascii, "num_cad_stl_write_ascii", "crate::kernels::cad::stl::write_stl_ascii",
    "Returns ASCII STL text with coordinates divided by the supplied unit length.",
    false, (mesh: Handle<CadMesh>, unit: Length, name: String) -> String,
    |ctx, mesh, unit, name| mesh.with(ctx, |mesh| {
        if name.len() > 128 || name.chars().any(char::is_control) || mesh.mesh.triangle_count() > MAX_REPAIR_ROWS {
            return Err(error("ASCII STL name or facet count exceeds the limits"));
        }
        let scaled = stl_coordinates(&mesh.mesh, unit, MAX_CAD_COORDINATE)?;
        let text = write_stl_ascii(&scaled, &name);
        if text.len() > MAX_STL_BYTES { return Err(error("ASCII STL exceeds the text budget")); }
        Ok(text)
    })
);

object_binding!(
    step_read_text, "num_cad_step_read_text", "crate::kernels::cad::step::read_step",
    "Reads bounded planar STEP B-reps with LINE edges directly from text. Each output is a measured mesh. Declared STEP length units are applied; otherwise one coordinate unit is one metre.",
    true, (text: String) -> Many<NewObject<CadMesh>>,
    |_, text| {
        let models = read_planar_step(text.as_bytes()).map_err(error)?;
        models.into_iter().map(|mesh| {
            if !topology(&mesh).closed {
                return Err(error("STEP triangulation did not yield a closed mesh"));
            }
            CadMesh::prepare(mesh)
        }).collect::<Result<_, _>>().map(Many)
    }
);

#[allow(
    clippy::single_call_fn,
    reason = "Keeps mesh registrations and dimensions with the adapters"
)]
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    const NAMESPACE: &str = "numerical_vector";
    api::install_with_dimension(
        i,
        &[
            &mesh_new::BINDING,
            &mesh_arrays::BINDING,
            &mesh_weld::BINDING,
            &mesh_repair::BINDING,
            &mesh_boolean::BINDING,
            &stl_read::BINDING,
            &stl_read_text::BINDING,
            &stl_write_ascii::BINDING,
            &stl_write_binary::BINDING,
            &step_read_text::BINDING,
        ],
        NAMESPACE,
        DimRule::Container,
    )?;
    api::install_with_dimension(
        i,
        &[
            &vertex_count::BINDING,
            &triangle_count::BINDING,
            &is_closed::BINDING,
            &is_manifold::BINDING,
            &euler::BINDING,
        ],
        NAMESPACE,
        DimRule::Dimensionless,
    )?;
    api::install_with_dimension(
        i,
        &[&signed_volume::BINDING],
        NAMESPACE,
        DimRule::Custom(|ctx, _| {
            ctx.base("Length", BigRational::from_integer(3.into()))
                .ok_or_else(|| {
                    DimError::UnknownUnit {
                        unit: "Length".into(),
                        context: None,
                    }
                })
        }),
    )
}
