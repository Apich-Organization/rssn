//! Mesh generation: Delaunay triangulation, constrained triangulation, quality
//! refinement and tetrahedralisation.
//!
//! | operator | value |
//! |---|---|
//! | `delaunay_2d(points)` | the counter-clockwise Delaunay triangulation of `points`, as `list(i, j, k)` index triples |
//! | `delaunay_3d(points)` | the Delaunay tetrahedralisation of `points`, as `list(i, j, k, l)` index quadruples |
//! | `cdt_2d(points, segments)` | the constrained Delaunay triangulation that keeps every segment as an edge |
//! | `refine_quality(points, triangles, min_angle, max_area)` | `list(points, triangles)` after inserting circumcentres until no triangle is below `min_angle` degrees or above `max_area` |
//!
//! `points` is a list of `list(x, y)` pairs (or `list(x, y, z)` triples for the
//! three-dimensional form), `segments` a list of `list(i, j)` index pairs, and
//! `triangles` a list of `list(i, j, k)` index triples. The indices refer to
//! the input list, so a caller keeps its own coordinates and the operators
//! return connectivity only — `refine_quality` is the exception, because it
//! inserts new points and must return them.
//!
//! The import and repair of mesh *files* (STL, STEP, welding, hole filling,
//! Boolean operations) is `crate::kernels::cad`: it reads a byte stream, which
//! a term cannot carry, so it is a library API rather than an operator.
//!
//! Degenerate input — fewer than three points, duplicate points, a segment
//! that names a missing vertex — is a [`DiagKind::Domain`] finding.

use super::bridge::Data;
use super::bridge::Entry;
use super::bridge::bad;
use super::bridge::require;
use super::bridge::scalar;
use crate::kernel::diag::DiagKind;

/// Most points a triangulation accepts.
const MAX_POINTS: usize = 50_000;

/// Every operator of this family.
pub(super) const ENTRIES: &[Entry] = &[
    Entry {
        name: "delaunay_2d",
        params: "t",
        eval: delaunay_2d,
        doc: "`delaunay_2d(points)`: the counter-clockwise Delaunay triangulation of a list of `list(x, y)` points, as `list(i, j, k)` index triples.",
    },
    Entry {
        name: "delaunay_3d",
        params: "t",
        eval: delaunay_3d,
        doc: "`delaunay_3d(points)`: the Delaunay tetrahedralisation of a list of `list(x, y, z)` points, as `list(i, j, k, l)` index quadruples.",
    },
    Entry {
        name: "cdt_2d",
        params: "tt",
        eval: cdt_2d,
        doc: "`cdt_2d(points, segments)`: the constrained Delaunay triangulation of `points` that keeps every `list(i, j)` segment as an edge.",
    },
    Entry {
        name: "refine_quality",
        params: "ttdd",
        eval: refine_quality,
        doc: "`refine_quality(points, triangles, min_angle, max_area)`: `list(points, triangles)` after Ruppert-style circumcentre insertion until no triangle is below `min_angle` degrees or above `max_area`.",
    },
];

/// The points of a two-dimensional list argument.
fn points2(
    args: &[Data],
    at: usize,
) -> Option<Vec<[f64; 2]>> {
    let Some(Data::Tree(entries)) = args.get(at) else {
        return bad(
            at,
            DiagKind::Shape,
            format_args!("expected a list of list(x, y) points"),
        );
    };
    require(
        !entries.is_empty() && entries.len() <= MAX_POINTS,
        Some(at),
        DiagKind::Shape,
        format_args!("expected 1 to {MAX_POINTS} points, got {}", entries.len()),
    )?;
    let mut points = Vec::with_capacity(entries.len());
    for (index, entry) in entries.iter().enumerate() {
        let Some(values) = numbers_of(entry) else {
            return bad(
                at,
                DiagKind::Shape,
                format_args!("point {index} must be a list of numbers"),
            );
        };
        let [x, y] = values.as_slice() else {
            return bad(
                at,
                DiagKind::Shape,
                format_args!("point {index} must have exactly two coordinates"),
            );
        };
        points.push([*x, *y]);
    }
    Some(points)
}

/// The points of a three-dimensional list argument.
fn points3(
    args: &[Data],
    at: usize,
) -> Option<Vec<[f64; 3]>> {
    let Some(Data::Tree(entries)) = args.get(at) else {
        return bad(
            at,
            DiagKind::Shape,
            format_args!("expected a list of list(x, y, z) points"),
        );
    };
    require(
        !entries.is_empty() && entries.len() <= MAX_POINTS,
        Some(at),
        DiagKind::Shape,
        format_args!("expected 1 to {MAX_POINTS} points, got {}", entries.len()),
    )?;
    let mut points = Vec::with_capacity(entries.len());
    for (index, entry) in entries.iter().enumerate() {
        let Some(values) = numbers_of(entry) else {
            return bad(
                at,
                DiagKind::Shape,
                format_args!("point {index} must be a list of numbers"),
            );
        };
        let [x, y, z] = values.as_slice() else {
            return bad(
                at,
                DiagKind::Shape,
                format_args!("point {index} must have exactly three coordinates"),
            );
        };
        points.push([*x, *y, *z]);
    }
    Some(points)
}

/// The numbers of a list literal.
fn numbers_of(data: &Data) -> Option<Vec<f64>> {
    match data {
        | Data::Array { shape, data } if shape.len() == 1 => Some(data.clone()),
        | _ => None,
    }
}

/// A list of index triples as a term.
fn triangles_out(triangles: Vec<[usize; 3]>) -> Data {
    Data::Tree(
        triangles
            .into_iter()
            .map(|corners| {
                Data::vector(
                    corners
                        .into_iter()
                        .map(|index| {
                            #[allow(clippy::cast_precision_loss)]
                            let value = index as f64;
                            value
                        })
                        .collect(),
                )
            })
            .collect(),
    )
}

fn delaunay_2d(args: &[Data]) -> Option<Data> {
    let points = points2(args, 0)?;
    match crate::kernels::meshgen::delaunay_2d(&points) {
        | Ok(triangles) => Some(triangles_out(triangles)),
        | Err(message) => bad(0, DiagKind::Domain, format_args!("{message}")),
    }
}

fn delaunay_3d(args: &[Data]) -> Option<Data> {
    let points = points3(args, 0)?;
    match crate::kernels::meshgen::delaunay_3d(&points) {
        | Ok(tetrahedra) => {
            Some(Data::Tree(
                tetrahedra
                    .into_iter()
                    .map(|corners| {
                        Data::vector(
                            corners
                                .into_iter()
                                .map(|index| {
                                    #[allow(clippy::cast_precision_loss)]
                                    let value = index as f64;
                                    value
                                })
                                .collect(),
                        )
                    })
                    .collect(),
            ))
        },
        | Err(message) => bad(0, DiagKind::Domain, format_args!("{message}")),
    }
}

fn cdt_2d(args: &[Data]) -> Option<Data> {
    let points = points2(args, 0)?;
    let Some(Data::Tree(entries)) = args.get(1) else {
        return bad(
            1,
            DiagKind::Shape,
            format_args!("expected a list of list(i, j) segments"),
        );
    };
    let mut segments = Vec::with_capacity(entries.len());
    for (index, entry) in entries.iter().enumerate() {
        let Some(values) = numbers_of(entry) else {
            return bad(
                1,
                DiagKind::Shape,
                format_args!("segment {index} must be a list of two indices"),
            );
        };
        let [from, to] = values.as_slice() else {
            return bad(
                1,
                DiagKind::Shape,
                format_args!("segment {index} must have exactly two indices"),
            );
        };
        let pair = |value: f64| {
            (value >= 0.0 && value.fract() == 0.0 && value < points.len() as f64)
                .then_some(value as usize)
        };
        let (Some(from), Some(to)) = (pair(*from), pair(*to)) else {
            return bad(
                1,
                DiagKind::Domain,
                format_args!("segment {index} names a vertex outside the point list"),
            );
        };
        segments.push((from, to));
    }
    match crate::kernels::meshgen::constrained_delaunay_2d(&points, &segments) {
        | Ok(triangles) => Some(triangles_out(triangles)),
        | Err(message) => bad(0, DiagKind::Domain, format_args!("{message}")),
    }
}

fn refine_quality(args: &[Data]) -> Option<Data> {
    let points = points2(args, 0)?;
    let Some(Data::Tree(entries)) = args.get(1) else {
        return bad(
            1,
            DiagKind::Shape,
            format_args!("expected a list of list(i, j, k) triangles"),
        );
    };
    let mut triangles = Vec::with_capacity(entries.len());
    for (index, entry) in entries.iter().enumerate() {
        let Some(values) = numbers_of(entry) else {
            return bad(
                1,
                DiagKind::Shape,
                format_args!("triangle {index} must be a list of three indices"),
            );
        };
        let [a, b, c] = values.as_slice() else {
            return bad(
                1,
                DiagKind::Shape,
                format_args!("triangle {index} must have exactly three indices"),
            );
        };
        let corner = |value: f64| {
            (value >= 0.0 && value.fract() == 0.0 && value < points.len() as f64)
                .then_some(value as usize)
        };
        let (Some(a), Some(b), Some(c)) = (corner(*a), corner(*b), corner(*c)) else {
            return bad(
                1,
                DiagKind::Domain,
                format_args!("triangle {index} names a vertex outside the point list"),
            );
        };
        triangles.push([a, b, c]);
    }
    let min_angle = scalar(args, 2)?;
    let max_area = scalar(args, 3)?;
    match crate::kernels::meshgen::refine_quality(&points, &triangles, min_angle, max_area) {
        | Ok(mesh) => {
            let refined_points = mesh
                .points
                .into_iter()
                .map(|[x, y]| Data::vector(vec![x, y]))
                .collect::<Vec<Data>>();
            Some(Data::Tree(vec![
                Data::Tree(refined_points),
                triangles_out(mesh.triangles),
            ]))
        },
        | Err(message) => bad(0, DiagKind::Domain, format_args!("{message}")),
    }
}
