// Portable spatial definitions and flat Voronoi topology. Numerical geometry
// remains in kernels::geometry; these adapters validate and encode its results.
portable_payload!(BvhPayload, "rssn/bvh-defining-data/v1");
struct BvhModel {
    payload: BvhPayload,
    coordinate_unit: f64,
}
record!(BvhModel {
    payload: BvhPayload,
    coordinate_unit: f64
});
impl Arg for BvhModel {
    fn ty() -> TypeExpr {
        <Self as Output>::ty()
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        let fields = ["payload", "coordinate_unit"];
        let coordinate_unit = numeric(model_field(value, "coordinate_unit", &fields)?, index)?;
        if !crate::kernels::exact_eq(coordinate_unit, 1.0) {
            return Err(model_error(
                "BVH coordinate_unit must be one canonical unit",
            ));
        }
        let payload = BvhPayload::decode(model_field(value, "payload", &fields)?, index)?;
        validate_bvh_rows(&payload.0)?;
        Ok(Self {
            payload,
            coordinate_unit,
        })
    }
}
fn validate_bvh_rows(rows: &[Vec<f64>]) -> Result<(), OpError> {
    if rows.len().checked_mul(24).is_none_or(|n| n > MODEL_LIMIT) {
        return Err(model_error(
            "BVH reconstruction exceeds portable allocation budget",
        ));
    }
    for row in rows {
        if row.len() != 6
            || row.iter().any(|v| !v.is_finite())
            || row[..3].iter().zip(&row[3..]).any(|(lo, hi)| lo > hi)
        {
            return Err(model_error(
                "BVH defining rows must contain finite ordered min/max XYZ corners",
            ));
        }
    }
    Ok(())
}
fn bvh_boxes(rows: &[Vec<f64>]) -> Result<Vec<crate::kernels::geometry::Aabb>, OpError> {
    validate_bvh_rows(rows)?;
    Ok(rows
        .iter()
        .map(|row| {
            crate::kernels::geometry::Aabb {
                min: [row[0], row[1], row[2]],
                max: [row[3], row[4], row[5]],
            }
        })
        .collect())
}
struct VoronoiResult {
    vertices: Vec<[f64; 2]>,
    edges: Vec<[usize; 2]>,
    ray_vertices: Vec<usize>,
    ray_directions: Vec<[f64; 2]>,
    cell_offsets: Vec<usize>,
    cell_indices: Vec<usize>,
}
record!(VoronoiResult{vertices:Vec<[f64;2]>,edges:Vec<[usize;2]>,ray_vertices:Vec<usize>,ray_directions:Vec<[f64;2]>,cell_offsets:Vec<usize>,cell_indices:Vec<usize>});
fn voronoi_output(points: &[[f64; 2]]) -> Result<VoronoiResult, OpError> {
    if points.len() < 3
        || points
            .len()
            .checked_mul(points.len())
            .is_none_or(|n| n > 1_048_576)
    {
        return Err(model_error(
            "Voronoi requires 3..=1024 sites within its finite pair-work budget",
        ));
    }
    let diagram = crate::kernels::geometry::voronoi_2d(points);
    if diagram.vertices.is_empty() {
        return Err(model_error(
            "Voronoi sites must span a nondegenerate planar triangulation",
        ));
    }
    let count = diagram.vertices.len();
    let mut cell_offsets = Vec::with_capacity(points.len() + 1);
    let mut cell_indices = Vec::new();
    cell_offsets.push(0);
    for cell in diagram.cells {
        if cell_indices
            .len()
            .checked_add(cell.len())
            .is_none_or(|n| n > MODEL_LIMIT)
            || cell.iter().any(|i| *i >= count)
        {
            return Err(model_error(
                "Voronoi cell topology exceeds its allocation or index bounds",
            ));
        }
        cell_indices.extend(cell);
        cell_offsets.push(cell_indices.len());
    }
    let edges = diagram.edges.into_iter().map(|(a, b)| [a, b]).collect();
    let (ray_vertices, ray_directions) = diagram.rays.into_iter().unzip();
    Ok(VoronoiResult {
        vertices: diagram.vertices,
        edges,
        ray_vertices,
        ray_directions,
        cell_offsets,
        cell_indices,
    })
}

binding!(voronoi_diagram,"num_voronoi_2d","numerical_vector","Planar Voronoi diagram dual to Delaunay, with measured vertices and dimensionless topology/unit ray directions. Input requires 3..=1024 finite nondegenerate sites. Cells are CSR-style cell_offsets/cell_indices; unbounded cells are described by ray_vertices/ray_directions rather than fabricated finite polygons.","crate::kernels::geometry::voronoi_2d",DimensionPolicy::VoronoiGeometry,
(points:Vec<[f64;2]>)->Result<crate::rules::numerical_abi::codecs::VoronoiResult,OpError>,|points:Vec<[f64;2]>|super::voronoi_output(&points));
binding!(bvh_build,"num_bvh_build","numerical_vector","Build immutable portable BVH defining data from same-unit N x 3 lower/upper corner tensors. Payload kind rssn/bvh-defining-data/v1 stores only canonical finite corner rows, never pointers or internal tree nodes. Queries validate/version-check the payload and reconstruct the native tree; portable reconstruction is limited to one million work-array coefficients.","crate::kernels::geometry::Bvh::build",DimensionPolicy::BvhBuild,
(lower:Vec<[f64;3]>,upper:Vec<[f64;3]>)->Result<crate::rules::numerical_abi::codecs::BvhModel,OpError>,|lower:Vec<[f64;3]>,upper:Vec<[f64;3]>|{if lower.len()!=upper.len()||lower.len().checked_mul(24).is_none_or(|n|n>super::MODEL_LIMIT){return Err(super::model_error("BVH corner tensors must have equal finite-budget row counts"));}let rows:Vec<_>=lower.into_iter().zip(upper).map(|(lo,hi)|lo.into_iter().chain(hi).collect()).collect();super::validate_bvh_rows(&rows)?;Ok(super::BvhModel{payload:super::BvhPayload(rows),coordinate_unit:1.0})});
binding!(bvh_nearest,"num_bvh_nearest","numerical_vector","Nearest AABB on a portable BVH. Query coordinates share the descriptor's coordinate units; index/distance tensors have zero or one entry, with measured distances. An empty definition gives empty results.","crate::kernels::geometry::Bvh::nearest_aabb",DimensionPolicy::SpatialQuery,
(model:crate::rules::numerical_abi::codecs::BvhModel,query:[f64;3])->Result<crate::rules::numerical_abi::codecs::Neighbors,OpError>,|model:super::BvhModel,query:[f64;3]|{let boxes=super::bvh_boxes(&model.payload.0)?;let bvh=crate::kernels::geometry::Bvh::build(&boxes);Ok(super::neighbors(bvh.nearest_aabb(query)))});
binding!(bvh_ray,"num_bvh_ray","numerical_vector","Indices of AABBs hit by origin+t*direction, t>=0, on an immutable portable BVH. Origin shares coordinate units; direction is dimensionless and nonzero. Payload kind/version, finite corners and allocation bounds are checked before reconstruction.","crate::kernels::geometry::Bvh::ray_intersect",DimensionPolicy::SpatialRay,
(model:crate::rules::numerical_abi::codecs::BvhModel,origin:[f64;3],direction:[f64;3])->Result<Vec<usize>,OpError>,|model:super::BvhModel,origin:[f64;3],direction:[f64;3]|{if direction.iter().all(|v|crate::kernels::exact_eq(*v,0.0)){return Err(super::model_error("BVH ray direction must be nonzero"));}let boxes=super::bvh_boxes(&model.payload.0)?;Ok(crate::kernels::geometry::Bvh::build(&boxes).ray_intersect(origin,direction))});

include!("mesh_codecs.rs");

/// Adapters for measured meshes, Voronoi topology and portable BVH data.
pub const GEOMETRY_EXPORTS: &[&CodecExport] = &[
    &voronoi_diagram::EXPORT,
    &meshgen_triangle_areas::EXPORT,
    &meshgen_tetra_volumes::EXPORT,
    &meshgen_min_angle_deg::EXPORT,
    &meshgen_boundary_faces::EXPORT,
    &meshgen_delaunay_2d_incremental::EXPORT,
    &bvh_build::EXPORT,
    &bvh_nearest::EXPORT,
    &bvh_ray::EXPORT,
];
