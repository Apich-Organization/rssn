use crate::kernels::meshgen;

fn validate_mesh_indices<const N: usize>(
    cells: &[[usize; N]],
    vertices: usize,
) -> Result<(), OpError> {
    for (row, cell) in cells.iter().enumerate() {
        if let Some(index) = cell.iter().find(|&&index| index >= vertices) {
            return Err(OpError::Domain(format!(
                "mesh cell {row} references vertex {index}, outside the {vertices} supplied vertices"
            )));
        }
    }
    Ok(())
}

binding!(
    meshgen_triangle_areas,
    "num_meshgen_triangle_areas",
    "numerical_vector",
    "`num_meshgen_triangle_areas(points, triangles)` returns an area for each triangle in input order, with squared coordinate units. Points have two columns and triangles have three zero-based vertex indices. Degenerate triangles have zero area. Out-of-range indices are errors.",
    "crate::kernels::meshgen::util::triangle_areas",
    DimensionPolicy::Square,
    (points: Vec<[f64; 2]>, triangles: Vec<[usize; 3]>) -> Result<Vec<f64>, OpError>,
    |points: Vec<[f64; 2]>, triangles: Vec<[usize; 3]>| {
        super::validate_mesh_indices(&triangles, points.len())?;
        Ok(super::meshgen::triangle_areas(&points, &triangles))
    }
);
binding!(
    meshgen_tetra_volumes,
    "num_meshgen_tetra_volumes",
    "numerical_vector",
    "`num_meshgen_tetra_volumes(points, tetrahedra)` returns nonnegative volumes in input order, with cubed coordinate units. Points have three columns and tetrahedra have four zero-based vertex indices. Degenerate tetrahedra have zero volume. Out-of-range indices are errors.",
    "crate::kernels::meshgen::tetra::tetra_volumes",
    DimensionPolicy::Cube,
    (points: Vec<[f64; 3]>, tetrahedra: Vec<[usize; 4]>) -> Result<Vec<f64>, OpError>,
    |points: Vec<[f64; 3]>, tetrahedra: Vec<[usize; 4]>| {
        super::validate_mesh_indices(&tetrahedra, points.len())?;
        Ok(super::meshgen::tetra_volumes(&points, &tetrahedra))
    }
);
binding!(
    meshgen_min_angle_deg,
    "num_meshgen_min_angle_deg",
    "numerical_vector",
    "`num_meshgen_min_angle_deg(points, triangles)` returns the smallest planar triangle angle as a scalar in degrees. Points have two columns and triangles have three zero-based vertex indices. An empty mesh or a degenerate triangle gives zero. Out-of-range indices are errors.",
    "crate::kernels::meshgen::util::min_angle_deg",
    DimensionPolicy::UnitlessResult,
    (points: Vec<[f64; 2]>, triangles: Vec<[usize; 3]>) -> Result<f64, OpError>,
    |points: Vec<[f64; 2]>, triangles: Vec<[usize; 3]>| {
        super::validate_mesh_indices(&triangles, points.len())?;
        Ok(super::meshgen::min_angle_deg(&points, &triangles))
    }
);
binding!(
    meshgen_boundary_faces,
    "num_meshgen_boundary_faces",
    "numerical_vector",
    "`num_meshgen_boundary_faces(tetrahedra)` returns the oriented triangular faces incident to exactly one tetrahedron, in first occurrence order. Tetrahedron rows have four distinct zero-based vertex indices. Positively oriented tetrahedra give outward face winding.",
    "crate::kernels::meshgen::tetra::boundary_faces",
    DimensionPolicy::Dimensionless,
    (tetrahedra: Vec<[usize; 4]>) -> Result<Vec<[usize; 3]>, OpError>,
    |tetrahedra: Vec<[usize; 4]>| {
        for (row, cell) in tetrahedra.iter().enumerate() {
            if cell.iter().enumerate().any(|(i, vertex)| cell[..i].contains(vertex)) {
                return Err(OpError::Domain(format!("tetrahedron {row} repeats a vertex")));
            }
        }
        Ok(super::meshgen::boundary_faces(&tetrahedra))
    }
);
binding!(
    meshgen_delaunay_2d_incremental,
    "num_meshgen_delaunay_2d_incremental",
    "numerical_vector",
    "`num_meshgen_delaunay_2d_incremental(points)` inserts between 3 and 4096 distinct planar points in input order and returns counter-clockwise triangles of zero-based vertex indices. Points may carry coordinate units. The triangulation contains indices only.",
    "crate::kernels::meshgen::delaunay::delaunay_2d_incremental",
    DimensionPolicy::UnitlessResult,
    (points: Vec<[f64; 2]>) -> Result<Vec<[usize; 3]>, OpError>,
    |points: Vec<[f64; 2]>| {
        if !(3..=4096).contains(&points.len()) {
            return Err(OpError::Domain(
                "incremental Delaunay triangulation requires between 3 and 4096 points".into(),
            ));
        }
        super::meshgen::delaunay_2d_incremental(&points).map_err(OpError::Domain)
    }
);
