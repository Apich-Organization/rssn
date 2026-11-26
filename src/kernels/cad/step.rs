//! A reader for ISO-10303-21 (STEP) boundary-representation files.
//!
//! The reader supports the AP203/AP214 entities describing solids, units,
//! curved edges, faces with holes, and trimmed analytic or B-spline surfaces:
//!
//! * Geometry & topology: `CARTESIAN_POINT`, `DIRECTION`, `AXIS2_PLACEMENT_3D`,
//!   `VERTEX_POINT`, `EDGE_CURVE`, `ORIENTED_EDGE`, `EDGE_LOOP`,
//!   `FACE_OUTER_BOUND`, `FACE_BOUND`, `ADVANCED_FACE`, `CLOSED_SHELL`, and
//!   `MANIFOLD_SOLID_BREP`.
//! * Curves: `LINE`, `CIRCLE`, `ELLIPSE`, `B_SPLINE_CURVE_WITH_KNOTS`, and
//!   rational NURBS curves (`RATIONAL_B_SPLINE_CURVE` in complex entities).
//! * Surfaces: `PLANE`, `CYLINDRICAL_SURFACE`, `SPHERICAL_SURFACE`,
//!   `CONICAL_SURFACE`, `B_SPLINE_SURFACE_WITH_KNOTS`, and rational NURBS
//!   surfaces (`RATIONAL_B_SPLINE_SURFACE`).
//! * Units: `SI_UNIT` (with SI prefixes for `.METRE.`) and
//!   `CONVERSION_BASED_UNIT` with `MEASURE_WITH_UNIT` /
//!   `LENGTH_MEASURE_WITH_UNIT`, scaling coordinates into metres when a length
//!   unit is declared (and leaving unitless files at scale `1.0`).
//!
//! Planar faces with inner `FACE_BOUND` holes are triangulated by projecting
//! onto the face plane, connecting holes to the outer loop via Eberly's bridge
//! algorithm, and ear-clipping. Curved faces are projected into surface `(u, v)`
//! parameter space, tessellated with interior grid refinement and Delaunay edge
//! flipping while respecting trimming boundaries and holes, and mapped back to
//! 3D.

use std::collections::HashMap;
use std::collections::HashSet;
use std::f64::consts::PI;
use std::path::Path;

use super::types::Mesh;
use super::types::cross;
use super::types::dot;
use super::types::edge_key;
use super::types::length;
use super::types::normalize;
use super::types::sub;

/// Number of polyline segments sampled along each curved `EDGE_CURVE`.
const CURVE_SEGMENTS: usize = 16;

/// Resolution of the interior `(u, v)` refinement grid on curved surfaces.
const SURFACE_GRID: usize = 12;

/// Reads every solid in a STEP file.
///
/// Returns one [`Mesh`] per `MANIFOLD_SOLID_BREP`, or one per `CLOSED_SHELL`
/// when the file has no solid wrapper. When the file declares a length unit
/// (`SI_UNIT` or `CONVERSION_BASED_UNIT`), coordinates are scaled into metres;
/// unitless files keep scale `1.0`.
///
/// # Errors
///
/// Returns `Err` when the bytes are not UTF-8, when there is no `DATA` section,
/// when there is no closed shell, or when an entity reference is undefined or
/// malformed.
pub fn read_step(bytes: &[u8]) -> Result<Vec<Mesh>, String> {
    let text =
        std::str::from_utf8(bytes).map_err(|_| "step file is not valid utf-8".to_string())?;

    let records = split_statements(text);
    let start = records
        .iter()
        .position(|record| record.trim().eq_ignore_ascii_case("DATA"))
        .ok_or_else(|| "step file has no data section".to_string())?;
    let end = records[start..]
        .iter()
        .position(|record| record.trim().eq_ignore_ascii_case("ENDSEC"))
        .map_or(records.len(), |offset| start + offset);

    let mut points: HashMap<usize, [f64; 3]> = HashMap::new();
    let mut directions: HashMap<usize, [f64; 3]> = HashMap::new();
    let mut placements: HashMap<usize, Axis2PlacementRaw> = HashMap::new();
    let mut vertex_points: HashMap<usize, usize> = HashMap::new();
    let mut curves: HashMap<usize, CurveGeom> = HashMap::new();
    let mut surfaces: HashMap<usize, SurfaceGeom> = HashMap::new();
    let mut edge_curves: HashMap<usize, EdgeCurve> = HashMap::new();
    let mut oriented_edges: HashMap<usize, (usize, bool)> = HashMap::new();
    let mut edge_loops: HashMap<usize, Vec<usize>> = HashMap::new();
    let mut face_bounds: HashMap<usize, FaceBound> = HashMap::new();
    let mut advanced_faces: HashMap<usize, AdvancedFace> = HashMap::new();
    let mut closed_shells: HashMap<usize, Vec<usize>> = HashMap::new();
    let mut solid_breps: Vec<usize> = Vec::new();

    // Unit-scaling tables.
    let mut si_length_units: HashMap<usize, f64> = HashMap::new();
    let mut measure_with_units: HashMap<usize, (f64, usize)> = HashMap::new();
    let mut conversion_units: Vec<usize> = Vec::new();

    for record in &records[start + 1..end] {
        let Some((id, parsed)) = parse_entity(record) else {
            continue;
        };
        match parsed {
            | ParsedEntity::Simple(kind, args) => {
                match kind.as_str() {
                    | "CARTESIAN_POINT" => {
                        if let Some(point) = parse_triple_arg(&args, 1) {
                            points.insert(id, point);
                        }
                    },
                    | "DIRECTION" => {
                        if let Some(dir) = parse_triple_arg(&args, 1) {
                            directions.insert(id, normalize(dir));
                        }
                    },
                    | "AXIS2_PLACEMENT_3D" => {
                        if let Some(location) = args.get(1).and_then(|a| parse_reference(a)) {
                            let axis = args.get(2).and_then(|a| parse_reference(a));
                            let ref_dir = args.get(3).and_then(|a| parse_reference(a));
                            placements.insert(
                                id,
                                Axis2PlacementRaw {
                                    location,
                                    axis,
                                    ref_dir,
                                },
                            );
                        }
                    },
                    | "VERTEX_POINT" => {
                        if let Some(point_id) = args.get(1).and_then(|a| parse_reference(a)) {
                            vertex_points.insert(id, point_id);
                        }
                    },
                    | "CIRCLE" => {
                        let placement = args.get(1).and_then(|a| parse_reference(a));
                        let radius = args.get(2).and_then(|a| a.trim().parse::<f64>().ok());
                        if let (Some(placement), Some(radius)) = (placement, radius) {
                            curves.insert(id, CurveGeom::Circle { placement, radius });
                        }
                    },
                    | "ELLIPSE" => {
                        let placement = args.get(1).and_then(|a| parse_reference(a));
                        let semi_axis_1 = args.get(2).and_then(|a| a.trim().parse::<f64>().ok());
                        let semi_axis_2 = args.get(3).and_then(|a| a.trim().parse::<f64>().ok());
                        if let (Some(placement), Some(semi_axis_1), Some(semi_axis_2)) =
                            (placement, semi_axis_1, semi_axis_2)
                        {
                            curves.insert(
                                id,
                                CurveGeom::Ellipse {
                                    placement,
                                    semi_axis_1,
                                    semi_axis_2,
                                },
                            );
                        }
                    },
                    | "B_SPLINE_CURVE_WITH_KNOTS" => {
                        if let Some(bspline) = parse_bspline_curve_simple(&args) {
                            curves.insert(id, CurveGeom::BSpline(bspline));
                        }
                    },
                    | "PLANE" => {
                        surfaces.insert(id, SurfaceGeom::Plane);
                    },
                    | "CYLINDRICAL_SURFACE" => {
                        let placement = args.get(1).and_then(|a| parse_reference(a));
                        let radius = args.get(2).and_then(|a| a.trim().parse::<f64>().ok());
                        if let (Some(placement), Some(radius)) = (placement, radius) {
                            surfaces.insert(id, SurfaceGeom::Cylinder { placement, radius });
                        }
                    },
                    | "SPHERICAL_SURFACE" => {
                        let placement = args.get(1).and_then(|a| parse_reference(a));
                        let radius = args.get(2).and_then(|a| a.trim().parse::<f64>().ok());
                        if let (Some(placement), Some(radius)) = (placement, radius) {
                            surfaces.insert(id, SurfaceGeom::Sphere { placement, radius });
                        }
                    },
                    | "CONICAL_SURFACE" => {
                        let placement = args.get(1).and_then(|a| parse_reference(a));
                        let radius = args.get(2).and_then(|a| a.trim().parse::<f64>().ok());
                        let semi_angle = args.get(3).and_then(|a| a.trim().parse::<f64>().ok());
                        if let (Some(placement), Some(radius), Some(semi_angle)) =
                            (placement, radius, semi_angle)
                        {
                            surfaces.insert(
                                id,
                                SurfaceGeom::Cone {
                                    placement,
                                    radius,
                                    semi_angle,
                                },
                            );
                        }
                    },
                    | "B_SPLINE_SURFACE_WITH_KNOTS" => {
                        if let Some(bspline) = parse_bspline_surface_simple(&args) {
                            surfaces.insert(id, SurfaceGeom::BSpline(bspline));
                        }
                    },
                    | "EDGE_CURVE" => {
                        let first = args.get(1).and_then(|a| parse_reference(a));
                        let second = args.get(2).and_then(|a| parse_reference(a));
                        let geom = args.get(3).and_then(|a| parse_reference(a));
                        let same_sense = args.get(4).is_none_or(|a| parse_boolean(a));
                        if let (Some(first), Some(second)) = (first, second) {
                            edge_curves.insert(
                                id,
                                EdgeCurve {
                                    start: first,
                                    end: second,
                                    geom,
                                    same_sense,
                                },
                            );
                        }
                    },
                    | "ORIENTED_EDGE" => {
                        let edge = args.get(3).and_then(|a| parse_reference(a));
                        let sense = args.get(4).is_none_or(|a| parse_boolean(a));
                        if let Some(edge) = edge {
                            oriented_edges.insert(id, (edge, sense));
                        }
                    },
                    | "EDGE_LOOP" => {
                        if let Some(list) = args.get(1).and_then(|a| parse_reference_list(a)) {
                            edge_loops.insert(id, list);
                        }
                    },
                    | "FACE_OUTER_BOUND" => {
                        if let Some(loop_id) = args.get(1).and_then(|a| parse_reference(a)) {
                            let orientation = args.get(2).is_none_or(|a| parse_boolean(a));
                            face_bounds.insert(
                                id,
                                FaceBound {
                                    loop_id,
                                    is_outer: true,
                                    orientation,
                                },
                            );
                        }
                    },
                    | "FACE_BOUND" => {
                        if let Some(loop_id) = args.get(1).and_then(|a| parse_reference(a)) {
                            let orientation = args.get(2).is_none_or(|a| parse_boolean(a));
                            face_bounds.insert(
                                id,
                                FaceBound {
                                    loop_id,
                                    is_outer: false,
                                    orientation,
                                },
                            );
                        }
                    },
                    | "ADVANCED_FACE" => {
                        let bounds = args
                            .get(1)
                            .and_then(|a| parse_reference_list(a))
                            .unwrap_or_default();
                        let surface = args.get(2).and_then(|a| parse_reference(a));
                        let sense = args.last().is_none_or(|a| parse_boolean(a));
                        advanced_faces.insert(
                            id,
                            AdvancedFace {
                                bounds,
                                surface,
                                sense,
                            },
                        );
                    },
                    | "CLOSED_SHELL" => {
                        if let Some(list) = args.get(1).and_then(|a| parse_reference_list(a)) {
                            closed_shells.insert(id, list);
                        }
                    },
                    | "MANIFOLD_SOLID_BREP" => {
                        if let Some(shell) = args.get(1).and_then(|a| parse_reference(a)) {
                            solid_breps.push(shell);
                        }
                    },
                    | "SI_UNIT" => {
                        if let Some(scale) = parse_si_length_unit(&args) {
                            si_length_units.insert(id, scale);
                        }
                    },
                    | "MEASURE_WITH_UNIT" | "LENGTH_MEASURE_WITH_UNIT" => {
                        if let (Some(val), Some(unit_id)) = (
                            args.first().and_then(|a| parse_measure_value(a)),
                            args.get(1).and_then(|a| parse_reference(a)),
                        ) {
                            measure_with_units.insert(id, (val, unit_id));
                        }
                    },
                    | "CONVERSION_BASED_UNIT" => {
                        if let Some(measure_id) = args.get(1).and_then(|a| parse_reference(a)) {
                            conversion_units.push(measure_id);
                        }
                    },
                    | _ => {},
                }
            },
            | ParsedEntity::Complex(components) => {
                // Units in complex entities:
                // e.g. `( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI., .METRE.) )`
                for (comp_kind, comp_args) in &components {
                    match comp_kind.as_str() {
                        | "SI_UNIT" => {
                            if let Some(scale) = parse_si_length_unit(comp_args) {
                                si_length_units.insert(id, scale);
                            }
                        },
                        | "CONVERSION_BASED_UNIT" => {
                            if let Some(measure_id) =
                                comp_args.get(1).and_then(|a| parse_reference(a))
                            {
                                conversion_units.push(measure_id);
                            }
                        },
                        | "MEASURE_WITH_UNIT" | "LENGTH_MEASURE_WITH_UNIT" => {
                            if let (Some(val), Some(unit_id)) = (
                                comp_args.first().and_then(|a| parse_measure_value(a)),
                                comp_args.get(1).and_then(|a| parse_reference(a)),
                            ) {
                                measure_with_units.insert(id, (val, unit_id));
                            }
                        },
                        | _ => {},
                    }
                }

                // Rational / non-rational B-spline curves in complex entities.
                if let Some(bspline) = parse_bspline_curve_complex(&components) {
                    curves.insert(id, CurveGeom::BSpline(bspline));
                }

                // Rational / non-rational B-spline surfaces in complex entities.
                if let Some(bspline) = parse_bspline_surface_complex(&components) {
                    surfaces.insert(id, SurfaceGeom::BSpline(bspline));
                }
            },
        }
    }

    // Determine length unit scale factor (metres per file unit; 1.0 when unitless).
    let unit_scale = resolve_unit_scale(&si_length_units, &measure_with_units, &conversion_units);
    if (unit_scale - 1.0).abs() > 0.0 {
        for pt in points.values_mut() {
            pt[0] *= unit_scale;
            pt[1] *= unit_scale;
            pt[2] *= unit_scale;
        }
        for curve in curves.values_mut() {
            match curve {
                | CurveGeom::Circle { radius, .. } => *radius *= unit_scale,
                | CurveGeom::Ellipse {
                    semi_axis_1,
                    semi_axis_2,
                    ..
                } => {
                    *semi_axis_1 *= unit_scale;
                    *semi_axis_2 *= unit_scale;
                },
                | CurveGeom::BSpline(_) => {},
            }
        }
        for surface in surfaces.values_mut() {
            match surface {
                | SurfaceGeom::Cylinder { radius, .. }
                | SurfaceGeom::Sphere { radius, .. }
                | SurfaceGeom::Cone { radius, .. } => *radius *= unit_scale,
                | SurfaceGeom::Plane | SurfaceGeom::BSpline(_) => {},
            }
        }
    }

    let context = Context {
        points: &points,
        directions: &directions,
        placements: &placements,
        vertex_points: &vertex_points,
        curves: &curves,
        surfaces: &surfaces,
        edge_curves: &edge_curves,
        oriented_edges: &oriented_edges,
        edge_loops: &edge_loops,
        face_bounds: &face_bounds,
        advanced_faces: &advanced_faces,
    };

    let mut meshes = Vec::new();
    if !solid_breps.is_empty() {
        for shell_id in solid_breps {
            let faces = closed_shells
                .get(&shell_id)
                .ok_or_else(|| format!("closed shell #{shell_id} is undefined"))?;
            meshes.push(build_shell(&context, faces)?);
        }
        return Ok(meshes);
    }

    if closed_shells.is_empty() {
        return Err("step file contains no closed shell".to_string());
    }

    let mut shells: Vec<&Vec<usize>> = closed_shells.values().collect();
    shells.sort_by_key(|faces| faces.first().copied().unwrap_or(usize::MAX));
    for faces in shells {
        meshes.push(build_shell(&context, faces)?);
    }

    Ok(meshes)
}

/// Reads every solid in the STEP file at `path`.
///
/// # Errors
///
/// Returns `Err` when the file cannot be read or when [`read_step`] fails.
pub fn read_step_file(path: &Path) -> Result<Vec<Mesh>, String> {
    let bytes = std::fs::read(path).map_err(|error| format!("cannot read step file: {error}"))?;
    read_step(&bytes)
}

// ---------------------------------------------------------------------------
// Internal geometry & topology representations
// ---------------------------------------------------------------------------

enum ParsedEntity {
    Simple(String, Vec<String>),
    Complex(Vec<(String, Vec<String>)>),
}

#[derive(Clone, Copy, Debug)]
struct Axis2PlacementRaw {
    location: usize,
    axis: Option<usize>,
    ref_dir: Option<usize>,
}

#[derive(Clone, Copy, Debug)]
struct Frame3D {
    origin: [f64; 3],
    x: [f64; 3],
    y: [f64; 3],
    z: [f64; 3],
}

#[derive(Clone, Debug)]
struct BSplineCurveData {
    degree: usize,
    control_points: Vec<usize>,
    multiplicities: Vec<usize>,
    knots: Vec<f64>,
    weights: Option<Vec<f64>>,
}

#[derive(Clone, Debug)]
enum CurveGeom {
    Circle {
        placement: usize,
        radius: f64,
    },
    Ellipse {
        placement: usize,
        semi_axis_1: f64,
        semi_axis_2: f64,
    },
    BSpline(BSplineCurveData),
}

#[derive(Clone, Debug)]
struct BSplineSurfaceData {
    u_degree: usize,
    v_degree: usize,
    control_points: Vec<Vec<usize>>,
    u_multiplicities: Vec<usize>,
    v_multiplicities: Vec<usize>,
    u_knots: Vec<f64>,
    v_knots: Vec<f64>,
    weights: Option<Vec<Vec<f64>>>,
}

#[derive(Clone, Debug)]
enum SurfaceGeom {
    Plane,
    Cylinder {
        placement: usize,
        radius: f64,
    },
    Sphere {
        placement: usize,
        radius: f64,
    },
    Cone {
        placement: usize,
        radius: f64,
        semi_angle: f64,
    },
    BSpline(BSplineSurfaceData),
}

#[derive(Clone, Copy, Debug)]
struct EdgeCurve {
    start: usize,
    end: usize,
    geom: Option<usize>,
    same_sense: bool,
}

#[derive(Clone, Copy, Debug)]
struct FaceBound {
    loop_id: usize,
    is_outer: bool,
    orientation: bool,
}

#[derive(Clone, Debug)]
struct AdvancedFace {
    bounds: Vec<usize>,
    surface: Option<usize>,
    sense: bool,
}

/// The entity tables a shell needs while it is being tessellated.
struct Context<'a> {
    points: &'a HashMap<usize, [f64; 3]>,
    directions: &'a HashMap<usize, [f64; 3]>,
    placements: &'a HashMap<usize, Axis2PlacementRaw>,
    vertex_points: &'a HashMap<usize, usize>,
    curves: &'a HashMap<usize, CurveGeom>,
    surfaces: &'a HashMap<usize, SurfaceGeom>,
    edge_curves: &'a HashMap<usize, EdgeCurve>,
    oriented_edges: &'a HashMap<usize, (usize, bool)>,
    edge_loops: &'a HashMap<usize, Vec<usize>>,
    face_bounds: &'a HashMap<usize, FaceBound>,
    advanced_faces: &'a HashMap<usize, AdvancedFace>,
}

impl Context<'_> {
    fn resolve_frame(
        &self,
        placement_id: usize,
    ) -> Option<Frame3D> {
        let raw = self.placements.get(&placement_id)?;
        let origin = *self.points.get(&raw.location)?;
        let z = raw
            .axis
            .and_then(|id| self.directions.get(&id).copied())
            .filter(|v| length(*v) > 1e-12)
            .unwrap_or([0.0, 0.0, 1.0]);
        let z = normalize(z);

        let ref_candidate = raw
            .ref_dir
            .and_then(|id| self.directions.get(&id).copied())
            .filter(|v| length(*v) > 1e-12)
            .unwrap_or_else(|| {
                if z[0].abs() < 0.9 {
                    [1.0, 0.0, 0.0]
                } else {
                    [0.0, 1.0, 0.0]
                }
            });
        let proj = dot(ref_candidate, z);
        let mut x = [
            ref_candidate[0] - proj * z[0],
            ref_candidate[1] - proj * z[1],
            ref_candidate[2] - proj * z[2],
        ];
        if length(x) <= 1e-12 {
            let fallback = if z[0].abs() < 0.9 {
                [1.0, 0.0, 0.0]
            } else {
                [0.0, 1.0, 0.0]
            };
            x = cross(z, fallback);
        }
        let x = normalize(x);
        let y = normalize(cross(z, x));
        Some(Frame3D { origin, x, y, z })
    }
}

// ---------------------------------------------------------------------------
// Shell & face tessellation
// ---------------------------------------------------------------------------

/// Tessellates every face of `faces` into one mesh.
fn build_shell(
    context: &Context<'_>,
    faces: &[usize],
) -> Result<Mesh, String> {
    let mut vertices: Vec<[f64; 3]> = Vec::new();
    let mut triangles: Vec<[usize; 3]> = Vec::new();
    let mut seen_points: HashMap<usize, usize> = HashMap::new();
    // Cache interior vertex indices sampled for each curved `EDGE_CURVE`
    // (stored in the edge's forward `start -> end` direction) so adjacent faces
    // share identical boundary vertices.
    let mut edge_samples: HashMap<usize, Vec<usize>> = HashMap::new();

    for face_id in faces {
        let face = context
            .advanced_faces
            .get(face_id)
            .ok_or_else(|| format!("advanced face #{face_id} is undefined"))?;
        if face.bounds.is_empty() {
            continue;
        }

        let mut outer_loop: Option<Vec<usize>> = None;
        let mut hole_loops: Vec<Vec<usize>> = Vec::new();

        for &bound_id in &face.bounds {
            let bound = context
                .face_bounds
                .get(&bound_id)
                .ok_or_else(|| format!("face bound #{bound_id} is undefined"))?;
            let mut ring = sample_edge_loop(
                context,
                bound.loop_id,
                &mut vertices,
                &mut seen_points,
                &mut edge_samples,
            )?;
            if ring.len() < 3 {
                continue;
            }
            if !bound.orientation {
                ring.reverse();
            }
            if bound.is_outer && outer_loop.is_none() {
                outer_loop = Some(ring);
            } else {
                hole_loops.push(ring);
            }
        }

        // If no bound was explicitly tagged `FACE_OUTER_BOUND`, promote the
        // first loop to outer.
        let Some(outer) = outer_loop.or_else(|| {
            if hole_loops.is_empty() {
                None
            } else {
                Some(hole_loops.remove(0))
            }
        }) else {
            continue;
        };

        let curved_surface = face
            .surface
            .and_then(|sid| context.surfaces.get(&sid))
            .and_then(|s| ResolvedSurface::from_geom(context, s));

        if let Some(surface) = curved_surface {
            tessellate_curved_face(
                &surface,
                &outer,
                &hole_loops,
                face.sense,
                &mut vertices,
                &mut triangles,
            );
        } else {
            tessellate_planar_face(&outer, &hole_loops, face.sense, &vertices, &mut triangles);
        }
    }

    if triangles.is_empty() {
        return Err("step shell has no tessellatable faces".to_string());
    }

    Ok(Mesh::new(vertices, triangles))
}

/// Samples an `EDGE_LOOP` into an ordered ring of vertex indices, tessellating
/// any curved `EDGE_CURVE`s into polylines while preserving exact corner indices.
fn sample_edge_loop(
    context: &Context<'_>,
    loop_id: usize,
    vertices: &mut Vec<[f64; 3]>,
    seen_points: &mut HashMap<usize, usize>,
    edge_samples: &mut HashMap<usize, Vec<usize>>,
) -> Result<Vec<usize>, String> {
    let oriented = context
        .edge_loops
        .get(&loop_id)
        .ok_or_else(|| format!("edge loop #{loop_id} is undefined"))?;

    let mut ring: Vec<usize> = Vec::new();

    for &oriented_id in oriented {
        let &(edge_id, forward) = context
            .oriented_edges
            .get(&oriented_id)
            .ok_or_else(|| format!("oriented edge #{oriented_id} is undefined"))?;
        let edge = context
            .edge_curves
            .get(&edge_id)
            .ok_or_else(|| format!("edge curve #{edge_id} is undefined"))?;

        let start_pt_id = *context
            .vertex_points
            .get(&edge.start)
            .ok_or_else(|| format!("vertex point #{} is undefined", edge.start))?;
        let end_pt_id = *context
            .vertex_points
            .get(&edge.end)
            .ok_or_else(|| format!("vertex point #{} is undefined", edge.end))?;

        let start_coord = *context
            .points
            .get(&start_pt_id)
            .ok_or_else(|| format!("cartesian point #{start_pt_id} is undefined"))?;
        let end_coord = *context
            .points
            .get(&end_pt_id)
            .ok_or_else(|| format!("cartesian point #{end_pt_id} is undefined"))?;

        let start_idx = intern_point(start_pt_id, start_coord, vertices, seen_points);
        let end_idx = intern_point(end_pt_id, end_coord, vertices, seen_points);

        // Compute or reuse interior curve samples for this `EDGE_CURVE`.
        let interior = if let Some(cached) = edge_samples.get(&edge_id) {
            cached.clone()
        } else {
            let sampled_pts = edge
                .geom
                .and_then(|gid| context.curves.get(&gid))
                .map(|geom| {
                    sample_curve_interior(context, geom, start_coord, end_coord, edge.same_sense)
                })
                .unwrap_or_default();
            let mut ids = Vec::with_capacity(sampled_pts.len());
            for pt in sampled_pts {
                let idx = vertices.len();
                vertices.push(pt);
                ids.push(idx);
            }
            edge_samples.insert(edge_id, ids.clone());
            ids
        };

        let from_idx = if forward {
            start_idx
        } else {
            end_idx
        };
        if ring.last() != Some(&from_idx) {
            ring.push(from_idx);
        }
        if forward {
            for &idx in &interior {
                if ring.last() != Some(&idx) {
                    ring.push(idx);
                }
            }
        } else {
            for &idx in interior.iter().rev() {
                if ring.last() != Some(&idx) {
                    ring.push(idx);
                }
            }
        }
    }

    if ring.len() >= 2 && ring[0] == ring[ring.len() - 1] {
        ring.pop();
    }
    Ok(ring)
}

fn intern_point(
    point_id: usize,
    coord: [f64; 3],
    vertices: &mut Vec<[f64; 3]>,
    seen: &mut HashMap<usize, usize>,
) -> usize {
    if let Some(&idx) = seen.get(&point_id) {
        return idx;
    }
    let idx = vertices.len();
    vertices.push(coord);
    seen.insert(point_id, idx);
    idx
}

// ---------------------------------------------------------------------------
// Curve evaluation
// ---------------------------------------------------------------------------

/// Samples `CURVE_SEGMENTS - 1` interior points along `geom` from `p0` to `p1`.
fn sample_curve_interior(
    context: &Context<'_>,
    geom: &CurveGeom,
    p0: [f64; 3],
    p1: [f64; 3],
    same_sense: bool,
) -> Vec<[f64; 3]> {
    match geom {
        | CurveGeom::Circle { placement, radius } => {
            let Some(frame) = context.resolve_frame(*placement) else {
                return Vec::new();
            };
            sample_conic_arc(frame, *radius, *radius, p0, p1, same_sense)
        },
        | CurveGeom::Ellipse {
            placement,
            semi_axis_1,
            semi_axis_2,
        } => {
            let Some(frame) = context.resolve_frame(*placement) else {
                return Vec::new();
            };
            sample_conic_arc(frame, *semi_axis_1, *semi_axis_2, p0, p1, same_sense)
        },
        | CurveGeom::BSpline(data) => sample_bspline_curve(context, data, same_sense),
    }
}

fn sample_conic_arc(
    frame: Frame3D,
    rx: f64,
    ry: f64,
    p0: [f64; 3],
    p1: [f64; 3],
    same_sense: bool,
) -> Vec<[f64; 3]> {
    if rx <= 1e-15 || ry <= 1e-15 {
        return Vec::new();
    }
    let d0 = sub(p0, frame.origin);
    let d1 = sub(p1, frame.origin);
    let theta0 = (dot(d0, frame.y) / ry).atan2(dot(d0, frame.x) / rx);
    let theta1 = (dot(d1, frame.y) / ry).atan2(dot(d1, frame.x) / rx);

    let full_circle = length(sub(p1, p0)) <= 1e-9 * rx.max(ry);
    let delta = if full_circle {
        if same_sense {
            2.0 * PI
        } else {
            -2.0 * PI
        }
    } else if same_sense {
        let mut d = (theta1 - theta0).rem_euclid(2.0 * PI);
        if d < 1e-12 {
            d = 2.0 * PI;
        }
        d
    } else {
        let mut d = (theta0 - theta1).rem_euclid(2.0 * PI);
        if d < 1e-12 {
            d = 2.0 * PI;
        }
        -d
    };

    let mut out = Vec::with_capacity(CURVE_SEGMENTS.saturating_sub(1));
    for i in 1..CURVE_SEGMENTS {
        let frac = i as f64 / CURVE_SEGMENTS as f64;
        let theta = theta0 + frac * delta;
        let cx = rx * theta.cos();
        let sy = ry * theta.sin();
        out.push([
            frame.origin[0] + cx * frame.x[0] + sy * frame.y[0],
            frame.origin[1] + cx * frame.x[1] + sy * frame.y[1],
            frame.origin[2] + cx * frame.x[2] + sy * frame.y[2],
        ]);
    }
    out
}

fn sample_bspline_curve(
    context: &Context<'_>,
    data: &BSplineCurveData,
    same_sense: bool,
) -> Vec<[f64; 3]> {
    let mut ctrl = Vec::with_capacity(data.control_points.len());
    for &pid in &data.control_points {
        let Some(&pt) = context.points.get(&pid) else {
            return Vec::new();
        };
        ctrl.push(pt);
    }
    if ctrl.len() < 2 {
        return Vec::new();
    }
    let degree = data.degree.min(ctrl.len() - 1);
    let knots = expand_knots(&data.multiplicities, &data.knots, degree, ctrl.len());
    let t_min = knots[degree];
    let t_max = knots[ctrl.len()];
    if (t_max - t_min).abs() <= 1e-15 {
        return Vec::new();
    }

    let weights = data.weights.as_deref().filter(|w| w.len() == ctrl.len());

    let mut out = Vec::with_capacity(CURVE_SEGMENTS.saturating_sub(1));
    for i in 1..CURVE_SEGMENTS {
        let frac = i as f64 / CURVE_SEGMENTS as f64;
        let u = if same_sense {
            t_min + frac * (t_max - t_min)
        } else {
            t_max - frac * (t_max - t_min)
        };
        out.push(eval_bspline_curve(&ctrl, weights, degree, &knots, u));
    }
    out
}

fn expand_knots(
    multiplicities: &[usize],
    knots: &[f64],
    degree: usize,
    n_ctrl: usize,
) -> Vec<f64> {
    let mut flat = Vec::new();
    if multiplicities.len() == knots.len() && !knots.is_empty() {
        for (&mult, &k) in multiplicities.iter().zip(knots.iter()) {
            for _ in 0..mult {
                flat.push(k);
            }
        }
    }
    if flat.len() == n_ctrl + degree + 1 {
        return flat;
    }

    let m = n_ctrl + degree + 1;
    let mut uniform = Vec::with_capacity(m);
    let internal = n_ctrl.saturating_sub(degree + 1);
    for _ in 0..=degree {
        uniform.push(0.0);
    }
    for i in 1..=internal {
        uniform.push(i as f64 / (internal + 1) as f64);
    }
    for _ in 0..=degree {
        uniform.push(1.0);
    }
    uniform
}

/// Evaluates all `n_ctrl` B-spline basis functions of `degree` at parameter `t`
/// using Cox-de Boor recursion.
fn bspline_basis(
    n_ctrl: usize,
    degree: usize,
    knots: &[f64],
    t: f64,
) -> Vec<f64> {
    let mut basis = vec![0.0; n_ctrl];
    if n_ctrl == 0 || knots.len() < n_ctrl + degree + 1 {
        return basis;
    }
    let t_min = knots[degree];
    let t_max = knots[n_ctrl];
    if t <= t_min {
        basis[0] = 1.0;
        return basis;
    }
    if t >= t_max {
        basis[n_ctrl - 1] = 1.0;
        return basis;
    }

    for i in 0..n_ctrl {
        if t >= knots[i] && t < knots[i + 1] {
            basis[i] = 1.0;
        }
    }

    for k in 1..=degree {
        let mut next = vec![0.0; n_ctrl];
        for i in 0..n_ctrl {
            let d1 = knots[i + k] - knots[i];
            let left = if d1.abs() > 1e-15 {
                ((t - knots[i]) / d1) * basis[i]
            } else {
                0.0
            };
            let right = if i + 1 < n_ctrl {
                let d2 = knots[i + k + 1] - knots[i + 1];
                if d2.abs() > 1e-15 {
                    ((knots[i + k + 1] - t) / d2) * basis[i + 1]
                } else {
                    0.0
                }
            } else {
                0.0
            };
            next[i] = left + right;
        }
        basis = next;
    }
    basis
}

fn eval_bspline_curve(
    ctrl: &[[f64; 3]],
    weights: Option<&[f64]>,
    degree: usize,
    knots: &[f64],
    t: f64,
) -> [f64; 3] {
    let basis = bspline_basis(ctrl.len(), degree, knots, t);
    let mut num = [0.0; 3];
    let mut den = 0.0;
    for (i, (&pt, &b)) in ctrl.iter().zip(basis.iter()).enumerate() {
        let w = weights.map_or(1.0, |ws| ws[i]);
        let wb = w * b;
        num[0] += wb * pt[0];
        num[1] += wb * pt[1];
        num[2] += wb * pt[2];
        den += wb;
    }
    if den.abs() <= 1e-15 {
        ctrl[0]
    } else {
        [num[0] / den, num[1] / den, num[2] / den]
    }
}

// ---------------------------------------------------------------------------
// Planar face triangulation (with FACE_BOUND holes)
// ---------------------------------------------------------------------------

fn tessellate_planar_face(
    outer: &[usize],
    holes: &[Vec<usize>],
    sense: bool,
    vertices: &[[f64; 3]],
    triangles: &mut Vec<[usize; 3]>,
) {
    let mut outer_ring = outer.to_vec();
    if !sense {
        outer_ring.reverse();
    }

    if holes.is_empty() && outer_ring.len() == 3 {
        triangles.push([outer_ring[0], outer_ring[1], outer_ring[2]]);
        return;
    }

    let mut normal = newell_normal(vertices, &outer_ring);
    if length(normal) < 0.5 {
        normal = [0.0, 0.0, 1.0];
    }
    let ref_vec = if normal[0].abs() < 0.9 {
        [1.0, 0.0, 0.0]
    } else {
        [0.0, 1.0, 0.0]
    };
    let u = normalize(cross(ref_vec, normal));
    let v = cross(normal, u);
    let origin = vertices[outer_ring[0]];

    let project = |idx: usize| -> (usize, [f64; 2]) {
        let rel = sub(vertices[idx], origin);
        (idx, [dot(rel, u), dot(rel, v)])
    };

    let mut outer_2d: Vec<(usize, [f64; 2])> = outer_ring.iter().map(|&i| project(i)).collect();
    let mut reversed_for_2d = false;
    if signed_area_2d(&outer_2d) < 0.0 {
        outer_2d.reverse();
        reversed_for_2d = true;
    }

    let mut holes_2d: Vec<Vec<(usize, [f64; 2])>> = Vec::with_capacity(holes.len());
    for hole in holes {
        if hole.len() < 3 {
            continue;
        }
        let mut h2d: Vec<(usize, [f64; 2])> = hole.iter().map(|&i| project(i)).collect();
        if signed_area_2d(&h2d) > 0.0 {
            h2d.reverse();
        }
        holes_2d.push(h2d);
    }

    let simple = bridge_holes_2d(outer_2d, holes_2d);
    let mut emitted = Vec::new();
    ear_clip_2d(&simple, &mut emitted);

    if reversed_for_2d {
        for tri in emitted {
            triangles.push([tri[0], tri[2], tri[1]]);
        }
    } else {
        triangles.extend(emitted);
    }
}

fn newell_normal(
    vertices: &[[f64; 3]],
    ring: &[usize],
) -> [f64; 3] {
    let mut n = [0.0; 3];
    let k = ring.len();
    for i in 0..k {
        let a = vertices[ring[i]];
        let b = vertices[ring[(i + 1) % k]];
        n[0] += (a[1] - b[1]) * (a[2] + b[2]);
        n[1] += (a[2] - b[2]) * (a[0] + b[0]);
        n[2] += (a[0] - b[0]) * (a[1] + b[1]);
    }
    normalize(n)
}

fn cross2(
    a: [f64; 2],
    b: [f64; 2],
    c: [f64; 2],
) -> f64 {
    (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
}

fn signed_area_2d(poly: &[(usize, [f64; 2])]) -> f64 {
    let n = poly.len();
    if n < 3 {
        return 0.0;
    }
    let mut sum = 0.0;
    for i in 0..n {
        let a = poly[i].1;
        let b = poly[(i + 1) % n].1;
        sum += a[0] * b[1] - b[0] * a[1];
    }
    0.5 * sum
}

fn point_in_tri_2d(
    p: [f64; 2],
    a: [f64; 2],
    b: [f64; 2],
    c: [f64; 2],
) -> bool {
    let area = cross2(a, b, c);
    if area.abs() <= 1e-14 {
        return false;
    }
    let sign = if area > 0.0 { 1.0 } else { -1.0 };
    cross2(a, b, p) * sign >= -1e-12
        && cross2(b, c, p) * sign >= -1e-12
        && cross2(c, a, p) * sign >= -1e-12
}

/// Connects CW inner holes to a CCW outer polygon using Eberly's horizontal-ray
/// bridge algorithm, returning a single simple polygon ready for ear-clipping.
fn bridge_holes_2d(
    mut outer: Vec<(usize, [f64; 2])>,
    mut holes: Vec<Vec<(usize, [f64; 2])>>,
) -> Vec<(usize, [f64; 2])> {
    if holes.is_empty() {
        return outer;
    }

    // Sort holes by rightmost x descending so bridges never cross unprocessed holes.
    holes.sort_by(|a, b| {
        let max_a = a.iter().map(|v| v.1[0]).fold(f64::NEG_INFINITY, f64::max);
        let max_b = b.iter().map(|v| v.1[0]).fold(f64::NEG_INFINITY, f64::max);
        max_b.total_cmp(&max_a)
    });

    for hole in holes {
        if hole.len() < 3 || outer.len() < 3 {
            continue;
        }
        let m_idx = (0..hole.len())
            .max_by(|&i, &j| {
                hole[i].1[0]
                    .total_cmp(&hole[j].1[0])
                    .then_with(|| hole[i].1[1].total_cmp(&hole[j].1[1]))
            })
            .unwrap_or(0);
        let m = hole[m_idx].1;

        // Cast horizontal ray (+x) from M against all edges of `outer`.
        let n = outer.len();
        let mut best_hit: Option<(usize, f64)> = None;
        for i in 0..n {
            let a = outer[i].1;
            let b = outer[(i + 1) % n].1;

            if (a[1] - m[1]).abs() <= 1e-12 && (b[1] - m[1]).abs() <= 1e-12 {
                let x_hi = a[0].max(b[0]);
                if x_hi >= m[0] - 1e-12 {
                    let ix = a[0].min(b[0]).max(m[0]);
                    if best_hit.is_none_or(|(_, best_x)| ix < best_x) {
                        best_hit = Some((i, ix));
                    }
                }
                continue;
            }

            if (a[1] <= m[1] && b[1] > m[1]) || (b[1] <= m[1] && a[1] > m[1]) {
                let t = (m[1] - a[1]) / (b[1] - a[1]);
                let ix = a[0] + t * (b[0] - a[0]);
                if ix >= m[0] - 1e-12 && best_hit.is_none_or(|(_, best_x)| ix < best_x) {
                    best_hit = Some((i, ix));
                }
            }
        }

        let visible_idx = if let Some((edge_i, ix)) = best_hit {
            let edge_j = (edge_i + 1) % n;
            let a = outer[edge_i].1;
            let b = outer[edge_j].1;
            let hit_pt = [ix, m[1]];

            if (a[0] - hit_pt[0]).hypot(a[1] - hit_pt[1]) <= 1e-10 {
                edge_i
            } else if (b[0] - hit_pt[0]).hypot(b[1] - hit_pt[1]) <= 1e-10 {
                edge_j
            } else {
                let p_idx = if a[0] >= b[0] {
                    edge_i
                } else {
                    edge_j
                };
                let p = outer[p_idx].1;

                // Check for reflex vertices of `outer` inside triangle (M, I, P).
                let mut best_reflex: Option<(usize, f64, f64)> = None;
                for r in 0..n {
                    if r == p_idx {
                        continue;
                    }
                    let prev = outer[(r + n - 1) % n].1;
                    let curr = outer[r].1;
                    let next = outer[(r + 1) % n].1;
                    if cross2(prev, curr, next) > 1e-12 {
                        continue;
                    }
                    if !point_in_tri_2d(curr, m, hit_pt, p) {
                        continue;
                    }
                    let dx = curr[0] - m[0];
                    let dy = curr[1] - m[1];
                    let angle = dy.abs().atan2(dx.max(0.0));
                    let dist_sq = dx * dx + dy * dy;
                    if best_reflex.is_none_or(|(_, ba, bd)| {
                        angle + 1e-12 < ba || ((angle - ba).abs() <= 1e-12 && dist_sq < bd)
                    }) {
                        best_reflex = Some((r, angle, dist_sq));
                    }
                }
                best_reflex.map_or(p_idx, |(r, _, _)| r)
            }
        } else {
            (0..n)
                .min_by(|&i, &j| {
                    let di = (outer[i].1[0] - m[0]).powi(2) + (outer[i].1[1] - m[1]).powi(2);
                    let dj = (outer[j].1[0] - m[0]).powi(2) + (outer[j].1[1] - m[1]).powi(2);
                    di.total_cmp(&dj)
                })
                .unwrap_or(0)
        };

        let mut merged = Vec::with_capacity(outer.len() + hole.len() + 2);
        merged.extend_from_slice(&outer[..=visible_idx]);
        for offset in 0..=hole.len() {
            merged.push(hole[(m_idx + offset) % hole.len()]);
        }
        merged.extend_from_slice(&outer[visible_idx..]);
        outer = merged;
    }

    outer
}

/// Ear-clipping triangulation of a CCW 2D polygon (allowing bridge duplicates).
fn ear_clip_2d(
    polygon: &[(usize, [f64; 2])],
    out: &mut Vec<[usize; 3]>,
) {
    let mut poly = polygon.to_vec();
    while poly.len() > 3 {
        let len = poly.len();
        let mut best_ear: Option<(usize, f64)> = None;

        for i in 0..len {
            let prev = poly[(i + len - 1) % len];
            let curr = poly[i];
            let next = poly[(i + 1) % len];
            if prev.0 == curr.0 || curr.0 == next.0 || next.0 == prev.0 {
                continue;
            }
            let area2 = cross2(prev.1, curr.1, next.1);
            if area2 <= 1e-12 {
                continue;
            }

            let mut empty = true;
            for j in 0..len {
                if j == i || j == (i + len - 1) % len || j == (i + 1) % len {
                    continue;
                }
                let pj = poly[j];
                if pj.0 == prev.0 || pj.0 == curr.0 || pj.0 == next.0 {
                    continue;
                }
                if (pj.1[0] - prev.1[0]).hypot(pj.1[1] - prev.1[1]) <= 1e-12
                    || (pj.1[0] - curr.1[0]).hypot(pj.1[1] - curr.1[1]) <= 1e-12
                    || (pj.1[0] - next.1[0]).hypot(pj.1[1] - next.1[1]) <= 1e-12
                {
                    continue;
                }
                if point_in_tri_2d(pj.1, prev.1, curr.1, next.1) {
                    empty = false;
                    break;
                }
            }
            if !empty {
                continue;
            }

            let perim = (curr.1[0] - prev.1[0]).hypot(curr.1[1] - prev.1[1])
                + (next.1[0] - curr.1[0]).hypot(next.1[1] - curr.1[1])
                + (prev.1[0] - next.1[0]).hypot(prev.1[1] - next.1[1]);
            let quality = area2 / (perim * perim).max(1e-18);
            if best_ear.is_none_or(|(_, q)| quality > q) {
                best_ear = Some((i, quality));
            }
        }

        let ear_idx = match best_ear {
            | Some((idx, _)) => idx,
            | None => {
                (0..len)
                    .max_by(|&a, &b| {
                        let ca = cross2(
                            poly[(a + len - 1) % len].1,
                            poly[a].1,
                            poly[(a + 1) % len].1,
                        );
                        let cb = cross2(
                            poly[(b + len - 1) % len].1,
                            poly[b].1,
                            poly[(b + 1) % len].1,
                        );
                        ca.total_cmp(&cb)
                    })
                    .unwrap_or(0)
            },
        };

        let prev_id = poly[(ear_idx + len - 1) % len].0;
        let curr_id = poly[ear_idx].0;
        let next_id = poly[(ear_idx + 1) % len].0;
        if prev_id != curr_id && curr_id != next_id && next_id != prev_id {
            out.push([prev_id, curr_id, next_id]);
        }
        poly.remove(ear_idx);
    }

    if poly.len() == 3 && poly[0].0 != poly[1].0 && poly[1].0 != poly[2].0 && poly[2].0 != poly[0].0
    {
        out.push([poly[0].0, poly[1].0, poly[2].0]);
    }
}

// ---------------------------------------------------------------------------
// Curved / trimmed surface evaluation & UV tessellation
// ---------------------------------------------------------------------------

enum ResolvedSurface {
    Cylinder {
        frame: Frame3D,
        radius: f64,
    },
    Sphere {
        frame: Frame3D,
        radius: f64,
    },
    Cone {
        frame: Frame3D,
        radius: f64,
        semi_angle: f64,
    },
    BSpline {
        u_degree: usize,
        v_degree: usize,
        ctrl: Vec<Vec<[f64; 3]>>,
        u_knots: Vec<f64>,
        v_knots: Vec<f64>,
        weights: Option<Vec<Vec<f64>>>,
    },
}

impl ResolvedSurface {
    fn from_geom(
        context: &Context<'_>,
        geom: &SurfaceGeom,
    ) -> Option<Self> {
        match geom {
            | SurfaceGeom::Plane => None,
            | SurfaceGeom::Cylinder { placement, radius } => {
                let frame = context.resolve_frame(*placement)?;
                Some(Self::Cylinder {
                    frame,
                    radius: *radius,
                })
            },
            | SurfaceGeom::Sphere { placement, radius } => {
                let frame = context.resolve_frame(*placement)?;
                Some(Self::Sphere {
                    frame,
                    radius: *radius,
                })
            },
            | SurfaceGeom::Cone {
                placement,
                radius,
                semi_angle,
            } => {
                let frame = context.resolve_frame(*placement)?;
                Some(Self::Cone {
                    frame,
                    radius: *radius,
                    semi_angle: *semi_angle,
                })
            },
            | SurfaceGeom::BSpline(data) => {
                let nu = data.control_points.len();
                if nu < 2 {
                    return None;
                }
                let nv = data.control_points[0].len();
                if nv < 2 || data.control_points.iter().any(|row| row.len() != nv) {
                    return None;
                }
                let mut ctrl = Vec::with_capacity(nu);
                for row in &data.control_points {
                    let mut pts = Vec::with_capacity(nv);
                    for &pid in row {
                        pts.push(*context.points.get(&pid)?);
                    }
                    ctrl.push(pts);
                }
                let u_degree = data.u_degree.min(nu - 1);
                let v_degree = data.v_degree.min(nv - 1);
                let u_knots = expand_knots(&data.u_multiplicities, &data.u_knots, u_degree, nu);
                let v_knots = expand_knots(&data.v_multiplicities, &data.v_knots, v_degree, nv);
                let weights = data.weights.as_ref().and_then(|w| {
                    if w.len() == nu && w.iter().all(|row| row.len() == nv) {
                        Some(w.clone())
                    } else {
                        None
                    }
                });
                Some(Self::BSpline {
                    u_degree,
                    v_degree,
                    ctrl,
                    u_knots,
                    v_knots,
                    weights,
                })
            },
        }
    }

    fn is_u_periodic(&self) -> bool {
        matches!(
            self,
            Self::Cylinder { .. } | Self::Sphere { .. } | Self::Cone { .. }
        )
    }

    fn eval(
        &self,
        u: f64,
        v: f64,
    ) -> [f64; 3] {
        match self {
            | Self::Cylinder { frame, radius } => {
                let cu = radius * u.cos();
                let su = radius * u.sin();
                [
                    frame.origin[0] + cu * frame.x[0] + su * frame.y[0] + v * frame.z[0],
                    frame.origin[1] + cu * frame.x[1] + su * frame.y[1] + v * frame.z[1],
                    frame.origin[2] + cu * frame.x[2] + su * frame.y[2] + v * frame.z[2],
                ]
            },
            | Self::Sphere { frame, radius } => {
                let cv = v.cos();
                let sv = v.sin();
                let cx = radius * cv * u.cos();
                let cy = radius * cv * u.sin();
                let cz = radius * sv;
                [
                    frame.origin[0] + cx * frame.x[0] + cy * frame.y[0] + cz * frame.z[0],
                    frame.origin[1] + cx * frame.x[1] + cy * frame.y[1] + cz * frame.z[1],
                    frame.origin[2] + cx * frame.x[2] + cy * frame.y[2] + cz * frame.z[2],
                ]
            },
            | Self::Cone {
                frame,
                radius,
                semi_angle,
            } => {
                let r = radius + v * semi_angle.sin();
                let z = v * semi_angle.cos();
                let cx = r * u.cos();
                let cy = r * u.sin();
                [
                    frame.origin[0] + cx * frame.x[0] + cy * frame.y[0] + z * frame.z[0],
                    frame.origin[1] + cx * frame.x[1] + cy * frame.y[1] + z * frame.z[1],
                    frame.origin[2] + cx * frame.x[2] + cy * frame.y[2] + z * frame.z[2],
                ]
            },
            | Self::BSpline {
                u_degree,
                v_degree,
                ctrl,
                u_knots,
                v_knots,
                weights,
            } => {
                let nu = ctrl.len();
                let nv = ctrl[0].len();
                let bu = bspline_basis(nu, *u_degree, u_knots, u);
                let bv = bspline_basis(nv, *v_degree, v_knots, v);
                let mut num = [0.0; 3];
                let mut den = 0.0;
                for i in 0..nu {
                    if bu[i].abs() <= 1e-18 {
                        continue;
                    }
                    for j in 0..nv {
                        if bv[j].abs() <= 1e-18 {
                            continue;
                        }
                        let w = weights.as_ref().map_or(1.0, |ws| ws[i][j]);
                        let wb = w * bu[i] * bv[j];
                        let pt = ctrl[i][j];
                        num[0] += wb * pt[0];
                        num[1] += wb * pt[1];
                        num[2] += wb * pt[2];
                        den += wb;
                    }
                }
                if den.abs() <= 1e-15 {
                    ctrl[0][0]
                } else {
                    [num[0] / den, num[1] / den, num[2] / den]
                }
            },
        }
    }

    fn project(
        &self,
        p: [f64; 3],
    ) -> [f64; 2] {
        match self {
            | Self::Cylinder { frame, .. } => {
                let d = sub(p, frame.origin);
                let u = dot(d, frame.y).atan2(dot(d, frame.x));
                let v = dot(d, frame.z);
                [u, v]
            },
            | Self::Sphere { frame, .. } => {
                let d = normalize(sub(p, frame.origin));
                let u = dot(d, frame.y).atan2(dot(d, frame.x));
                let v = dot(d, frame.z).clamp(-1.0, 1.0).asin();
                [u, v]
            },
            | Self::Cone {
                frame,
                radius,
                semi_angle,
            } => {
                let d = sub(p, frame.origin);
                let x = dot(d, frame.x);
                let y = dot(d, frame.y);
                let z = dot(d, frame.z);
                let u = y.atan2(x);
                let v = if semi_angle.cos().abs() > 1e-9 {
                    z / semi_angle.cos()
                } else {
                    (x.hypot(y) - radius) / semi_angle.sin().max(1e-12)
                };
                [u, v]
            },
            | Self::BSpline {
                u_degree,
                v_degree,
                ctrl,
                u_knots,
                v_knots,
                ..
            } => {
                let u_min = u_knots[*u_degree];
                let u_max = u_knots[ctrl.len()];
                let v_min = v_knots[*v_degree];
                let v_max = v_knots[ctrl[0].len()];

                // Coarse grid search for initial (u, v).
                let mut best_u = 0.5 * (u_min + u_max);
                let mut best_v = 0.5 * (v_min + v_max);
                let mut best_d2 = f64::INFINITY;
                let samples = 12_usize;
                for iu in 0..=samples {
                    let u = u_min + (iu as f64 / samples as f64) * (u_max - u_min);
                    for iv in 0..=samples {
                        let v = v_min + (iv as f64 / samples as f64) * (v_max - v_min);
                        let diff = sub(self.eval(u, v), p);
                        let d2 = dot(diff, diff);
                        if d2 < best_d2 {
                            best_d2 = d2;
                            best_u = u;
                            best_v = v;
                        }
                    }
                }

                // Damped Gauss-Newton refinement.
                let eps_u = ((u_max - u_min) * 1e-6).max(1e-9);
                let eps_v = ((v_max - v_min) * 1e-6).max(1e-9);
                for _ in 0..10 {
                    let s = self.eval(best_u, best_v);
                    let r = sub(s, p);
                    let d2 = dot(r, r);
                    if d2 <= 1e-20 {
                        break;
                    }
                    let u_hi = (best_u + eps_u).min(u_max);
                    let u_lo = (best_u - eps_u).max(u_min);
                    let v_hi = (best_v + eps_v).min(v_max);
                    let v_lo = (best_v - eps_v).max(v_min);
                    let du_span = (u_hi - u_lo).max(1e-12);
                    let dv_span = (v_hi - v_lo).max(1e-12);
                    let su = sub(self.eval(u_hi, best_v), self.eval(u_lo, best_v));
                    let sv = sub(self.eval(best_u, v_hi), self.eval(best_u, v_lo));
                    let su = [su[0] / du_span, su[1] / du_span, su[2] / du_span];
                    let sv = [sv[0] / dv_span, sv[1] / dv_span, sv[2] / dv_span];

                    let a = dot(su, su) + 1e-8;
                    let b = dot(su, sv);
                    let c = dot(sv, sv) + 1e-8;
                    let g0 = dot(su, r);
                    let g1 = dot(sv, r);
                    let det = a * c - b * b;
                    if det.abs() <= 1e-18 {
                        break;
                    }
                    let step_u = -(c * g0 - b * g1) / det;
                    let step_v = -(-b * g0 + a * g1) / det;

                    let mut accepted = false;
                    for alpha in [1.0, 0.5, 0.25, 0.1] {
                        let cand_u = (best_u + alpha * step_u).clamp(u_min, u_max);
                        let cand_v = (best_v + alpha * step_v).clamp(v_min, v_max);
                        let diff = sub(self.eval(cand_u, cand_v), p);
                        let cand_d2 = dot(diff, diff);
                        if cand_d2 < d2 {
                            best_u = cand_u;
                            best_v = cand_v;
                            accepted = true;
                            break;
                        }
                    }
                    if !accepted {
                        break;
                    }
                }

                [best_u, best_v]
            },
        }
    }
}

/// Projects boundary loops onto the curved surface's `(u, v)` domain, inserts an
/// interior grid strictly inside the trimmed region (respecting holes), runs
/// Lawson Delaunay flipping on unconstrained edges, and maps vertices back to 3D.
fn tessellate_curved_face(
    surface: &ResolvedSurface,
    outer: &[usize],
    holes: &[Vec<usize>],
    sense: bool,
    vertices: &mut Vec<[f64; 3]>,
    triangles: &mut Vec<[usize; 3]>,
) {
    let periodic = surface.is_u_periodic();

    let project_loop = |ring: &[usize]| -> Vec<(usize, [f64; 2])> {
        let mut pts: Vec<(usize, [f64; 2])> = ring
            .iter()
            .map(|&idx| (idx, surface.project(vertices[idx])))
            .collect();
        if periodic && pts.len() >= 2 {
            for i in 1..pts.len() {
                let prev_u = pts[i - 1].1[0];
                let mut u = pts[i].1[0];
                while u - prev_u > PI {
                    u -= 2.0 * PI;
                }
                while u - prev_u < -PI {
                    u += 2.0 * PI;
                }
                pts[i].1[0] = u;
            }
        }
        pts
    };

    let mut outer_uv = project_loop(outer);
    if outer_uv.len() < 3 {
        return;
    }

    let outer_u_mean = outer_uv.iter().map(|p| p.1[0]).sum::<f64>() / outer_uv.len() as f64;
    if signed_area_2d(&outer_uv) < 0.0 {
        outer_uv.reverse();
    }

    let mut holes_uv: Vec<Vec<(usize, [f64; 2])>> = Vec::with_capacity(holes.len());
    for hole in holes {
        if hole.len() < 3 {
            continue;
        }
        let mut h_uv = project_loop(hole);
        if periodic && !h_uv.is_empty() {
            let h_mean = h_uv.iter().map(|p| p.1[0]).sum::<f64>() / h_uv.len() as f64;
            let shift = ((outer_u_mean - h_mean) / (2.0 * PI)).round() * 2.0 * PI;
            for p in &mut h_uv {
                p.1[0] += shift;
            }
        }
        if signed_area_2d(&h_uv) > 0.0 {
            h_uv.reverse();
        }
        holes_uv.push(h_uv);
    }

    // Track constrained boundary edges (outer + holes) so Lawson flipping never
    // modifies the trimming boundary.
    let mut constrained: HashSet<(usize, usize)> = HashSet::new();
    for i in 0..outer_uv.len() {
        constrained.insert(edge_key(
            outer_uv[i].0,
            outer_uv[(i + 1) % outer_uv.len()].0,
        ));
    }
    for hole in &holes_uv {
        for i in 0..hole.len() {
            constrained.insert(edge_key(hole[i].0, hole[(i + 1) % hole.len()].0));
        }
    }

    // Record UV coordinates for every vertex index participating in this face.
    let mut uv_map: HashMap<usize, [f64; 2]> = HashMap::new();
    let mut u_min = f64::INFINITY;
    let mut u_max = f64::NEG_INFINITY;
    let mut v_min = f64::INFINITY;
    let mut v_max = f64::NEG_INFINITY;
    for &(id, uv) in &outer_uv {
        uv_map.insert(id, uv);
        u_min = u_min.min(uv[0]);
        u_max = u_max.max(uv[0]);
        v_min = v_min.min(uv[1]);
        v_max = v_max.max(uv[1]);
    }
    for hole in &holes_uv {
        for &(id, uv) in hole {
            uv_map.insert(id, uv);
        }
    }

    let simple_uv = bridge_holes_2d(outer_uv.clone(), holes_uv.clone());
    let mut uv_tris: Vec<[usize; 3]> = Vec::new();
    ear_clip_2d(&simple_uv, &mut uv_tris);

    let du = u_max - u_min;
    let dv = v_max - v_min;
    if du > 1e-12 && dv > 1e-12 && !uv_tris.is_empty() {
        // Generate interior UV grid points strictly inside the outer loop and
        // strictly outside all holes.
        let mut interior_pts: Vec<(usize, [f64; 2])> = Vec::new();
        for iu in 1..SURFACE_GRID {
            let u = u_min + (iu as f64 / SURFACE_GRID as f64) * du;
            for iv in 1..SURFACE_GRID {
                let v = v_min + (iv as f64 / SURFACE_GRID as f64) * dv;
                let pt = [u, v];
                if !point_in_polygon_2d(pt, &outer_uv) {
                    continue;
                }
                if holes_uv.iter().any(|h| point_in_polygon_2d(pt, h)) {
                    continue;
                }
                let vid = vertices.len();
                vertices.push(surface.eval(u, v));
                uv_map.insert(vid, pt);
                interior_pts.push((vid, pt));
            }
        }

        // Insert each interior point into the containing UV triangle.
        for (vid, pt) in interior_pts {
            let mut split_idx: Option<usize> = None;
            for (ti, &[a, b, c]) in uv_tris.iter().enumerate() {
                let (Some(&pa), Some(&pb), Some(&pc)) =
                    (uv_map.get(&a), uv_map.get(&b), uv_map.get(&c))
                else {
                    continue;
                };
                let area = cross2(pa, pb, pc);
                if area <= 1e-14 {
                    continue;
                }
                let w0 = cross2(pb, pc, pt) / area;
                let w1 = cross2(pc, pa, pt) / area;
                let w2 = cross2(pa, pb, pt) / area;
                if w0 > 1e-5 && w1 > 1e-5 && w2 > 1e-5 {
                    split_idx = Some(ti);
                    break;
                }
            }
            if let Some(ti) = split_idx {
                let [a, b, c] = uv_tris.swap_remove(ti);
                uv_tris.push([a, b, vid]);
                uv_tris.push([b, c, vid]);
                uv_tris.push([c, a, vid]);
            }
        }

        // Lawson Delaunay edge flipping in normalized UV space.
        lawson_flip_uv(
            &mut uv_tris,
            &uv_map,
            &constrained,
            [u_min, v_min],
            [du, dv],
        );
    }

    // In ISO 10303-42, `dS/du x dS/dv` points outward for cylinders, spheres,
    // cones, and standard B-spline surfaces, so CCW triangles in UV match
    // `sense == true`.
    for [a, b, c] in uv_tris {
        if sense {
            triangles.push([a, b, c]);
        } else {
            triangles.push([a, c, b]);
        }
    }
}

fn point_in_polygon_2d(
    p: [f64; 2],
    poly: &[(usize, [f64; 2])],
) -> bool {
    let n = poly.len();
    let mut inside = false;
    for i in 0..n {
        let a = poly[i].1;
        let b = poly[(i + 1) % n].1;
        if (a[1] > p[1]) != (b[1] > p[1]) {
            let x_int = a[0] + (p[1] - a[1]) * (b[0] - a[0]) / (b[1] - a[1]);
            if p[0] < x_int {
                inside = !inside;
            }
        }
    }
    inside
}

/// Performs Lawson Delaunay edge flipping on unconstrained internal UV edges.
fn lawson_flip_uv(
    tris: &mut [[usize; 3]],
    uv_map: &HashMap<usize, [f64; 2]>,
    constrained: &HashSet<(usize, usize)>,
    uv_min: [f64; 2],
    uv_span: [f64; 2],
) {
    let norm = |id: usize| -> [f64; 2] {
        let p = uv_map[&id];
        [
            (p[0] - uv_min[0]) / uv_span[0],
            (p[1] - uv_min[1]) / uv_span[1],
        ]
    };

    let max_passes = 8;
    for _ in 0..max_passes {
        let mut edge_to_tri: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
        for (ti, t) in tris.iter().enumerate() {
            for k in 0..3 {
                edge_to_tri
                    .entry(edge_key(t[k], t[(k + 1) % 3]))
                    .or_default()
                    .push(ti);
            }
        }

        let mut flipped_any = false;
        let mut modified = vec![false; tris.len()];

        for (&key, adj) in &edge_to_tri {
            if constrained.contains(&key) || adj.len() != 2 {
                continue;
            }
            let t0 = adj[0];
            let t1 = adj[1];
            if modified[t0] || modified[t1] {
                continue;
            }

            let Some((a, b, c)) = align_edge(tris[t0], key) else {
                continue;
            };
            let Some((b2, a2, d)) = align_edge(tris[t1], key) else {
                continue;
            };
            if !(a == a2 && b == b2) || c == d {
                continue;
            }

            let pa = norm(a);
            let pb = norm(b);
            let pc = norm(c);
            let pd = norm(d);

            // Quad (a, d, b, c) must be strictly convex and `d` strictly inside
            // circumcircle of (a, b, c).
            if cross2(pc, pa, pd) <= 1e-10 || cross2(pd, pb, pc) <= 1e-10 {
                continue;
            }
            if in_circle_2d(pa, pb, pc, pd) > 1e-10 {
                tris[t0] = [c, a, d];
                tris[t1] = [c, d, b];
                modified[t0] = true;
                modified[t1] = true;
                flipped_any = true;
            }
        }

        if !flipped_any {
            break;
        }
    }
}

/// Rotates `tri` so its first two vertices match the undirected edge `key`,
/// returning `(u, v, opposite)`.
fn align_edge(
    tri: [usize; 3],
    key: (usize, usize),
) -> Option<(usize, usize, usize)> {
    for k in 0..3 {
        let u = tri[k];
        let v = tri[(k + 1) % 3];
        let opp = tri[(k + 2) % 3];
        if edge_key(u, v) == key {
            return Some((u, v, opp));
        }
    }
    None
}

/// Positive when `d` lies strictly inside the circumcircle of CCW `(a, b, c)`.
fn in_circle_2d(
    a: [f64; 2],
    b: [f64; 2],
    c: [f64; 2],
    d: [f64; 2],
) -> f64 {
    let adx = a[0] - d[0];
    let ady = a[1] - d[1];
    let bdx = b[0] - d[0];
    let bdy = b[1] - d[1];
    let cdx = c[0] - d[0];
    let cdy = c[1] - d[1];

    let ad2 = adx * adx + ady * ady;
    let bd2 = bdx * bdx + bdy * bdy;
    let cd2 = cdx * cdx + cdy * cdy;

    adx * (bdy * cd2 - bd2 * cdy) - ady * (bdx * cd2 - bd2 * cdx) + ad2 * (bdx * cdy - bdy * cdx)
}

// ---------------------------------------------------------------------------
// STEP lexical & entity parsers
// ---------------------------------------------------------------------------

/// Splits a STEP file into `;`-terminated statements outside of strings.
fn split_statements(text: &str) -> Vec<&str> {
    let mut out = Vec::new();
    let bytes = text.as_bytes();
    let mut start = 0;
    let mut in_string = false;
    let mut index = 0;

    while index < bytes.len() {
        match bytes[index] {
            | b'\'' => in_string = !in_string,
            | b';' if !in_string => {
                out.push(&text[start..index]);
                start = index + 1;
            },
            | _ => {},
        }
        index += 1;
    }
    if start < text.len() {
        out.push(&text[start..]);
    }
    out
}

/// Parses `#12 = TYPE(...)` or a complex entity `#12 = ( SUB1(...) SUB2(...) )`.
fn parse_entity(record: &str) -> Option<(usize, ParsedEntity)> {
    let equals = record.find('=')?;
    let id = record[..equals]
        .trim()
        .strip_prefix('#')?
        .trim()
        .parse()
        .ok()?;
    let rest = record[equals + 1..].trim();

    if let Some(inner) = rest.strip_prefix('(').and_then(|s| s.strip_suffix(')')) {
        let components = parse_complex_components(inner);
        if components.is_empty() {
            return None;
        }
        return Some((id, ParsedEntity::Complex(components)));
    }

    let open = rest.find('(')?;
    let close = rest.rfind(')')?;
    if close <= open {
        return None;
    }
    let kind = rest[..open].trim().to_ascii_uppercase();
    if kind.is_empty() {
        return None;
    }

    Some((
        id,
        ParsedEntity::Simple(kind, split_arguments(&rest[open + 1..close])),
    ))
}

/// Parses the sub-entities inside a STEP complex entity instance.
fn parse_complex_components(raw: &str) -> Vec<(String, Vec<String>)> {
    let mut out = Vec::new();
    let bytes = raw.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        while i < bytes.len() && bytes[i].is_ascii_whitespace() {
            i += 1;
        }
        let name_start = i;
        while i < bytes.len() && (bytes[i].is_ascii_alphanumeric() || bytes[i] == b'_') {
            i += 1;
        }
        let name = raw[name_start..i].trim().to_ascii_uppercase();
        while i < bytes.len() && bytes[i].is_ascii_whitespace() {
            i += 1;
        }
        if i >= bytes.len() || bytes[i] != b'(' {
            break;
        }
        let args_start = i + 1;
        let mut depth = 1_i32;
        let mut in_string = false;
        i += 1;
        while i < bytes.len() && depth > 0 {
            match bytes[i] {
                | b'\'' => in_string = !in_string,
                | b'(' if !in_string => depth += 1,
                | b')' if !in_string => depth -= 1,
                | _ => {},
            }
            i += 1;
        }
        if depth == 0 && !name.is_empty() {
            let args_end = i - 1;
            out.push((name, split_arguments(&raw[args_start..args_end])));
        }
    }
    out
}

/// Splits the argument list of an entity at its top-level commas.
fn split_arguments(raw: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut depth = 0_i32;
    let mut in_string = false;
    let mut start = 0;

    for (index, ch) in raw.char_indices() {
        match ch {
            | '\'' => in_string = !in_string,
            | '(' if !in_string => depth += 1,
            | ')' if !in_string => depth -= 1,
            | ',' if !in_string && depth == 0 => {
                out.push(raw[start..index].trim().to_string());
                start = index + 1;
            },
            | _ => {},
        }
    }

    let tail = raw[start..].trim();
    if !tail.is_empty() {
        out.push(tail.to_string());
    }
    out
}

/// Reads a 3-element numeric tuple `(x, y, z)` at `args[index]`.
fn parse_triple_arg(
    args: &[String],
    index: usize,
) -> Option<[f64; 3]> {
    let list = args.get(index)?;
    let inner = list.trim().strip_prefix('(')?.strip_suffix(')')?;
    let parts: Vec<&str> = inner.split(',').map(str::trim).collect();
    if parts.len() < 3 {
        return None;
    }
    let x = parts[0].parse::<f64>().ok()?;
    let y = parts[1].parse::<f64>().ok()?;
    let z = parts[2].parse::<f64>().ok()?;
    Some([x, y, z])
}

/// Reads an entity reference such as `#42`.
fn parse_reference(raw: &str) -> Option<usize> {
    raw.trim().strip_prefix('#')?.trim().parse().ok()
}

/// Reads a parenthesised list of entity references such as `(#1,#2,#3)`.
fn parse_reference_list(raw: &str) -> Option<Vec<usize>> {
    let inner = raw.trim().strip_prefix('(')?.strip_suffix(')')?;
    let mut out = Vec::new();
    for part in inner.split(',') {
        let part = part.trim();
        if part.is_empty() {
            continue;
        }
        out.push(parse_reference(part)?);
    }
    Some(out)
}

/// Reads a parenthesised list of floating-point numbers.
fn parse_f64_list(raw: &str) -> Option<Vec<f64>> {
    let inner = raw.trim().strip_prefix('(')?.strip_suffix(')')?;
    let mut out = Vec::new();
    for part in inner.split(',') {
        let part = part.trim();
        if part.is_empty() {
            continue;
        }
        out.push(part.parse::<f64>().ok()?);
    }
    Some(out)
}

/// Reads a parenthesised list of unsigned integers.
fn parse_usize_list(raw: &str) -> Option<Vec<usize>> {
    let inner = raw.trim().strip_prefix('(')?.strip_suffix(')')?;
    let mut out = Vec::new();
    for part in inner.split(',') {
        let part = part.trim();
        if part.is_empty() {
            continue;
        }
        out.push(part.parse::<usize>().ok()?);
    }
    Some(out)
}

/// Reads a 2D parenthesised grid of entity references `((#1,#2),(#3,#4))`.
fn parse_reference_grid(raw: &str) -> Option<Vec<Vec<usize>>> {
    let inner = raw.trim().strip_prefix('(')?.strip_suffix(')')?;
    let rows = split_arguments(inner);
    let mut out = Vec::with_capacity(rows.len());
    for row in rows {
        out.push(parse_reference_list(&row)?);
    }
    Some(out)
}

/// Reads a 2D parenthesised grid of floats `((1.,1.),(1.,1.))`.
fn parse_f64_grid(raw: &str) -> Option<Vec<Vec<f64>>> {
    let inner = raw.trim().strip_prefix('(')?.strip_suffix(')')?;
    let rows = split_arguments(inner);
    let mut out = Vec::with_capacity(rows.len());
    for row in rows {
        out.push(parse_f64_list(&row)?);
    }
    Some(out)
}

/// Reads a STEP boolean, where `.T.` is true and anything else is false.
fn parse_boolean(raw: &str) -> bool {
    raw.trim().eq_ignore_ascii_case(".T.")
}

/// Reads `LENGTH_MEASURE(25.4)` or a bare float `25.4`.
fn parse_measure_value(raw: &str) -> Option<f64> {
    let trimmed = raw.trim();
    if let (Some(open), Some(close)) = (trimmed.find('('), trimmed.rfind(')')) {
        if close > open {
            return trimmed[open + 1..close].trim().parse::<f64>().ok();
        }
    }
    trimmed.parse::<f64>().ok()
}

/// Parses `SI_UNIT(prefix, .METRE.)` into metres-per-unit.
fn parse_si_length_unit(args: &[String]) -> Option<f64> {
    if args.len() < 2 {
        return None;
    }
    let prefix = args[args.len() - 2].trim();
    let base = args[args.len() - 1].trim();
    if !base.eq_ignore_ascii_case(".METRE.") {
        return None;
    }
    let factor = match prefix.to_ascii_uppercase().as_str() {
        | ".EXA." => 1e18,
        | ".PETA." => 1e15,
        | ".TERA." => 1e12,
        | ".GIGA." => 1e9,
        | ".MEGA." => 1e6,
        | ".KILO." => 1e3,
        | ".HECTO." => 1e2,
        | ".DECA." => 1e1,
        | "$" | "*" => 1.0,
        | ".DECI." => 1e-1,
        | ".CENTI." => 1e-2,
        | ".MILLI." => 1e-3,
        | ".MICRO." => 1e-6,
        | ".NANO." => 1e-9,
        | ".PICO." => 1e-12,
        | ".FEMTO." => 1e-15,
        | ".ATTO." => 1e-18,
        | _ => return None,
    };
    Some(factor)
}

fn resolve_unit_scale(
    si_length_units: &HashMap<usize, f64>,
    measure_with_units: &HashMap<usize, (f64, usize)>,
    conversion_units: &[usize],
) -> f64 {
    // Prefer an explicit CONVERSION_BASED_UNIT (e.g. INCH = 25.4 mm) when present.
    for &measure_id in conversion_units {
        if let Some(&(val, unit_id)) = measure_with_units.get(&measure_id) {
            if let Some(&si_factor) = si_length_units.get(&unit_id) {
                return val * si_factor;
            }
        }
    }
    if let Some((&_, &factor)) = si_length_units.iter().min_by_key(|&(&id, _)| id) {
        return factor;
    }
    1.0
}

fn parse_bspline_curve_simple(args: &[String]) -> Option<BSplineCurveData> {
    // B_SPLINE_CURVE_WITH_KNOTS(name, degree, (#cp...), form, closed, self_int, (mults...), (knots...), spec)
    let degree = args.get(1)?.trim().parse::<usize>().ok()?;
    let control_points = parse_reference_list(args.get(2)?)?;
    let multiplicities = args
        .get(6)
        .and_then(|a| parse_usize_list(a))
        .unwrap_or_default();
    let knots = args
        .get(7)
        .and_then(|a| parse_f64_list(a))
        .unwrap_or_default();
    Some(BSplineCurveData {
        degree,
        control_points,
        multiplicities,
        knots,
        weights: None,
    })
}

fn parse_bspline_curve_complex(components: &[(String, Vec<String>)]) -> Option<BSplineCurveData> {
    let mut degree: Option<usize> = None;
    let mut control_points: Option<Vec<usize>> = None;
    let mut multiplicities: Vec<usize> = Vec::new();
    let mut knots: Vec<f64> = Vec::new();
    let mut weights: Option<Vec<f64>> = None;

    for (kind, args) in components {
        match kind.as_str() {
            | "B_SPLINE_CURVE" => {
                degree = args.first().and_then(|a| a.trim().parse::<usize>().ok());
                control_points = args.get(1).and_then(|a| parse_reference_list(a));
            },
            | "B_SPLINE_CURVE_WITH_KNOTS" => {
                if let Some(m) = args.first().and_then(|a| parse_usize_list(a)) {
                    multiplicities = m;
                }
                if let Some(k) = args.get(1).and_then(|a| parse_f64_list(a)) {
                    knots = k;
                }
            },
            | "RATIONAL_B_SPLINE_CURVE" => {
                weights = args.first().and_then(|a| parse_f64_list(a));
            },
            | _ => {},
        }
    }

    Some(BSplineCurveData {
        degree: degree?,
        control_points: control_points?,
        multiplicities,
        knots,
        weights,
    })
}

fn parse_bspline_surface_simple(args: &[String]) -> Option<BSplineSurfaceData> {
    // B_SPLINE_SURFACE_WITH_KNOTS(name, u_deg, v_deg, ((#cp...)), form, u_cl, v_cl, self_int, (u_mult), (v_mult), (u_knots), (v_knots), spec)
    let u_degree = args.get(1)?.trim().parse::<usize>().ok()?;
    let v_degree = args.get(2)?.trim().parse::<usize>().ok()?;
    let control_points = parse_reference_grid(args.get(3)?)?;
    let u_multiplicities = args
        .get(8)
        .and_then(|a| parse_usize_list(a))
        .unwrap_or_default();
    let v_multiplicities = args
        .get(9)
        .and_then(|a| parse_usize_list(a))
        .unwrap_or_default();
    let u_knots = args
        .get(10)
        .and_then(|a| parse_f64_list(a))
        .unwrap_or_default();
    let v_knots = args
        .get(11)
        .and_then(|a| parse_f64_list(a))
        .unwrap_or_default();
    Some(BSplineSurfaceData {
        u_degree,
        v_degree,
        control_points,
        u_multiplicities,
        v_multiplicities,
        u_knots,
        v_knots,
        weights: None,
    })
}

fn parse_bspline_surface_complex(
    components: &[(String, Vec<String>)]
) -> Option<BSplineSurfaceData> {
    let mut u_degree: Option<usize> = None;
    let mut v_degree: Option<usize> = None;
    let mut control_points: Option<Vec<Vec<usize>>> = None;
    let mut u_multiplicities: Vec<usize> = Vec::new();
    let mut v_multiplicities: Vec<usize> = Vec::new();
    let mut u_knots: Vec<f64> = Vec::new();
    let mut v_knots: Vec<f64> = Vec::new();
    let mut weights: Option<Vec<Vec<f64>>> = None;

    for (kind, args) in components {
        match kind.as_str() {
            | "B_SPLINE_SURFACE" => {
                u_degree = args.first().and_then(|a| a.trim().parse::<usize>().ok());
                v_degree = args.get(1).and_then(|a| a.trim().parse::<usize>().ok());
                control_points = args.get(2).and_then(|a| parse_reference_grid(a));
            },
            | "B_SPLINE_SURFACE_WITH_KNOTS" => {
                if let Some(m) = args.first().and_then(|a| parse_usize_list(a)) {
                    u_multiplicities = m;
                }
                if let Some(m) = args.get(1).and_then(|a| parse_usize_list(a)) {
                    v_multiplicities = m;
                }
                if let Some(k) = args.get(2).and_then(|a| parse_f64_list(a)) {
                    u_knots = k;
                }
                if let Some(k) = args.get(3).and_then(|a| parse_f64_list(a)) {
                    v_knots = k;
                }
            },
            | "RATIONAL_B_SPLINE_SURFACE" => {
                weights = args.first().and_then(|a| parse_f64_grid(a));
            },
            | _ => {},
        }
    }

    Some(BSplineSurfaceData {
        u_degree: u_degree?,
        v_degree: v_degree?,
        control_points: control_points?,
        u_multiplicities,
        v_multiplicities,
        u_knots,
        v_knots,
        weights,
    })
}

#[cfg(test)]
mod tests {
    use super::read_step;
    use crate::kernels::cad::types::Mesh;
    use crate::kernels::cad::types::cross;
    use crate::kernels::cad::types::length;
    use crate::kernels::cad::types::sub;

    fn expect_ok<T>(result: Result<T, String>) -> T {
        match result {
            | Ok(value) => value,
            | Err(message) => panic!("unexpected error: {message}"),
        }
    }

    fn surface_area(mesh: &Mesh) -> f64 {
        mesh.triangles
            .iter()
            .map(|&[a, b, c]| {
                0.5 * length(cross(
                    sub(mesh.vertices[b], mesh.vertices[a]),
                    sub(mesh.vertices[c], mesh.vertices[a]),
                ))
            })
            .sum()
    }

    const TRIANGLE: &str = "\
ISO-10303-21;
HEADER;
FILE_DESCRIPTION((''),'2;1');
FILE_NAME('tri.step','2024-01-01',(''),(''),'','','');
FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));
ENDSEC;
DATA;
#1 = CARTESIAN_POINT('',(0.,0.,0.));
#2 = CARTESIAN_POINT('',(1.,0.,0.));
#3 = CARTESIAN_POINT('',(0.,1.,0.));
#4 = VERTEX_POINT('',#1);
#5 = VERTEX_POINT('',#2);
#6 = VERTEX_POINT('',#3);
#7 = LINE('',#1,#2);
#8 = EDGE_CURVE('',#4,#5,#7,.T.);
#9 = EDGE_CURVE('',#5,#6,#7,.T.);
#10 = EDGE_CURVE('',#6,#4,#7,.T.);
#11 = ORIENTED_EDGE('',*,*,#8,.T.);
#12 = ORIENTED_EDGE('',*,*,#9,.T.);
#13 = ORIENTED_EDGE('',*,*,#10,.T.);
#14 = EDGE_LOOP('',(#11,#12,#13));
#15 = FACE_OUTER_BOUND('',#14,.T.);
#16 = ADVANCED_FACE('',(#15),#7,.T.);
#17 = CLOSED_SHELL('',(#16));
#18 = MANIFOLD_SOLID_BREP('tri',#17);
ENDSEC;
END-ISO-10303-21;
";

    #[test]
    fn single_triangular_face_is_read() {
        let meshes = expect_ok(read_step(TRIANGLE.as_bytes()));
        assert_eq!(meshes.len(), 1);
        assert_eq!(meshes[0].triangles.len(), 1);
        assert_eq!(meshes[0].vertices.len(), 3);
        assert_eq!(meshes[0].vertices[1], [1.0, 0.0, 0.0]);
        assert_eq!(meshes[0].vertices[2], [0.0, 1.0, 0.0]);
    }

    #[test]
    fn malformed_files_are_rejected() {
        assert!(read_step(b"").is_err());
        assert!(read_step(b"ISO-10303-21;HEADER;ENDSEC;END-ISO-10303-21;").is_err());
        assert!(read_step(b"DATA;ENDSEC;").is_err());
    }

    #[test]
    fn undefined_references_are_reported() {
        let text = "\
DATA;
#1 = CARTESIAN_POINT('',(0.,0.,0.));
#2 = VERTEX_POINT('',#99);
ENDSEC;
";
        assert!(read_step(text.as_bytes()).is_err());
    }

    #[test]
    fn file_reader_reports_a_missing_path() {
        assert!(super::read_step_file(std::path::Path::new("does-not-exist.step")).is_err());
    }

    #[test]
    fn si_millimetre_unit_scales_coordinates_to_metres() {
        let text = "\
ISO-10303-21;
DATA;
#100 = ( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI., .METRE.) );
#1 = CARTESIAN_POINT('',(0.,0.,0.));
#2 = CARTESIAN_POINT('',(1000.,0.,0.));
#3 = CARTESIAN_POINT('',(0.,500.,0.));
#4 = VERTEX_POINT('',#1);
#5 = VERTEX_POINT('',#2);
#6 = VERTEX_POINT('',#3);
#7 = LINE('',#1,#2);
#8 = EDGE_CURVE('',#4,#5,#7,.T.);
#9 = EDGE_CURVE('',#5,#6,#7,.T.);
#10 = EDGE_CURVE('',#6,#4,#7,.T.);
#11 = ORIENTED_EDGE('',*,*,#8,.T.);
#12 = ORIENTED_EDGE('',*,*,#9,.T.);
#13 = ORIENTED_EDGE('',*,*,#10,.T.);
#14 = EDGE_LOOP('',(#11,#12,#13));
#15 = FACE_OUTER_BOUND('',#14,.T.);
#16 = ADVANCED_FACE('',(#15),#7,.T.);
#17 = CLOSED_SHELL('',(#16));
#18 = MANIFOLD_SOLID_BREP('mm_tri',#17);
ENDSEC;
END-ISO-10303-21;
";
        let meshes = expect_ok(read_step(text.as_bytes()));
        assert_eq!(meshes.len(), 1);
        assert!((meshes[0].vertices[1][0] - 1.0).abs() < 1e-12);
        assert!((meshes[0].vertices[2][1] - 0.5).abs() < 1e-12);
        assert!((surface_area(&meshes[0]) - 0.25).abs() < 1e-12);
    }

    #[test]
    fn b_spline_curve_edge_is_sampled_accurately() {
        // Quadratic Bezier B-spline from (0,0,0) through control (0.5, 1, 0) to (1,0,0):
        // y(t) = 2*t*(1-t), so area under the arch above y=0 is integral_0^1 2t(1-t) dt = 1/3.
        let text = "\
ISO-10303-21;
DATA;
#1 = CARTESIAN_POINT('',(0.,0.,0.));
#2 = CARTESIAN_POINT('',(0.5,1.,0.));
#3 = CARTESIAN_POINT('',(1.,0.,0.));
#4 = VERTEX_POINT('',#1);
#5 = VERTEX_POINT('',#3);
#6 = B_SPLINE_CURVE_WITH_KNOTS('',2,(#1,#2,#3),.UNSPECIFIED.,.F.,.F.,(3,3),(0.,1.),.UNSPECIFIED.);
#7 = LINE('',#3,#1);
#8 = EDGE_CURVE('',#4,#5,#6,.T.);
#9 = EDGE_CURVE('',#5,#4,#7,.T.);
#10 = ORIENTED_EDGE('',*,*,#8,.T.);
#11 = ORIENTED_EDGE('',*,*,#9,.T.);
#12 = EDGE_LOOP('',(#10,#11));
#13 = FACE_OUTER_BOUND('',#12,.T.);
#14 = ADVANCED_FACE('',(#13),#7,.T.);
#15 = CLOSED_SHELL('',(#14));
#16 = MANIFOLD_SOLID_BREP('arch',#15);
ENDSEC;
END-ISO-10303-21;
";
        let meshes = expect_ok(read_step(text.as_bytes()));
        assert_eq!(meshes.len(), 1);
        let mesh = &meshes[0];
        // 16 segments on the B-spline -> 17 boundary vertices -> 15 triangles.
        assert!(mesh.vertices.len() >= 17);
        assert!(mesh.triangles.len() >= 15);
        // Apex at t=0.5 has y = 0.5.
        let max_y = mesh
            .vertices
            .iter()
            .map(|v| v[1])
            .fold(f64::NEG_INFINITY, f64::max);
        assert!((max_y - 0.5).abs() < 1e-12, "max_y was {max_y}");
        // Area of 16-segment chordal approximation of parabola (exact = 1/3 ≈ 0.33333).
        let area = surface_area(mesh);
        assert!(
            (area - 1.0 / 3.0).abs() < 1e-2,
            "arch area was {area}, expected ~0.3333"
        );
    }

    #[test]
    fn planar_face_with_face_bound_square_hole_subtracts_hole_area() {
        // Outer 4x4 square (area 16) minus inner 2x2 square hole (area 4) = area 12.
        let text = "\
ISO-10303-21;
DATA;
#1 = CARTESIAN_POINT('',(0.,0.,0.));
#2 = CARTESIAN_POINT('',(4.,0.,0.));
#3 = CARTESIAN_POINT('',(4.,4.,0.));
#4 = CARTESIAN_POINT('',(0.,4.,0.));
#5 = CARTESIAN_POINT('',(1.,1.,0.));
#6 = CARTESIAN_POINT('',(3.,1.,0.));
#7 = CARTESIAN_POINT('',(3.,3.,0.));
#8 = CARTESIAN_POINT('',(1.,3.,0.));
#11 = VERTEX_POINT('',#1);
#12 = VERTEX_POINT('',#2);
#13 = VERTEX_POINT('',#3);
#14 = VERTEX_POINT('',#4);
#15 = VERTEX_POINT('',#5);
#16 = VERTEX_POINT('',#6);
#17 = VERTEX_POINT('',#7);
#18 = VERTEX_POINT('',#8);
#20 = LINE('',#1,#2);
#21 = EDGE_CURVE('',#11,#12,#20,.T.);
#22 = EDGE_CURVE('',#12,#13,#20,.T.);
#23 = EDGE_CURVE('',#13,#14,#20,.T.);
#24 = EDGE_CURVE('',#14,#11,#20,.T.);
#25 = EDGE_CURVE('',#15,#18,#20,.T.);
#26 = EDGE_CURVE('',#18,#17,#20,.T.);
#27 = EDGE_CURVE('',#17,#16,#20,.T.);
#28 = EDGE_CURVE('',#16,#15,#20,.T.);
#31 = ORIENTED_EDGE('',*,*,#21,.T.);
#32 = ORIENTED_EDGE('',*,*,#22,.T.);
#33 = ORIENTED_EDGE('',*,*,#23,.T.);
#34 = ORIENTED_EDGE('',*,*,#24,.T.);
#35 = ORIENTED_EDGE('',*,*,#25,.T.);
#36 = ORIENTED_EDGE('',*,*,#26,.T.);
#37 = ORIENTED_EDGE('',*,*,#27,.T.);
#38 = ORIENTED_EDGE('',*,*,#28,.T.);
#40 = EDGE_LOOP('',(#31,#32,#33,#34));
#41 = EDGE_LOOP('',(#35,#36,#37,#38));
#42 = FACE_OUTER_BOUND('',#40,.T.);
#43 = FACE_BOUND('',#41,.T.);
#44 = AXIS2_PLACEMENT_3D('',#1,#50,#51);
#50 = DIRECTION('',(0.,0.,1.));
#51 = DIRECTION('',(1.,0.,0.));
#45 = PLANE('',#44);
#46 = ADVANCED_FACE('',(#42,#43),#45,.T.);
#47 = CLOSED_SHELL('',(#46));
#48 = MANIFOLD_SOLID_BREP('plate_with_hole',#47);
ENDSEC;
END-ISO-10303-21;
";
        let meshes = expect_ok(read_step(text.as_bytes()));
        assert_eq!(meshes.len(), 1);
        let area = surface_area(&meshes[0]);
        assert!(
            (area - 12.0).abs() < 1e-9,
            "plate with hole area was {area}, expected 12.0"
        );
    }

    #[test]
    fn cylindrical_surface_trimmed_patch_has_expected_radius_and_area() {
        // Quarter-cylinder patch: radius = 2, height = 3, angle u in [0, pi/2].
        // Exact lateral surface area = (pi/2) * R * H = 3 * pi ≈ 9.42478.
        let text = "\
ISO-10303-21;
DATA;
#1 = CARTESIAN_POINT('',(0.,0.,0.));
#2 = CARTESIAN_POINT('',(0.,0.,3.));
#3 = DIRECTION('',(0.,0.,1.));
#4 = DIRECTION('',(1.,0.,0.));
#5 = AXIS2_PLACEMENT_3D('',#1,#3,#4);
#6 = AXIS2_PLACEMENT_3D('',#2,#3,#4);
#7 = CYLINDRICAL_SURFACE('',#5,2.);
#8 = CIRCLE('',#5,2.);
#9 = CIRCLE('',#6,2.);
#10 = CARTESIAN_POINT('',(2.,0.,0.));
#11 = CARTESIAN_POINT('',(0.,2.,0.));
#12 = CARTESIAN_POINT('',(0.,2.,3.));
#13 = CARTESIAN_POINT('',(2.,0.,3.));
#14 = VERTEX_POINT('',#10);
#15 = VERTEX_POINT('',#11);
#16 = VERTEX_POINT('',#12);
#17 = VERTEX_POINT('',#13);
#18 = LINE('',#11,#12);
#19 = LINE('',#13,#10);
#20 = EDGE_CURVE('',#14,#15,#8,.T.);
#21 = EDGE_CURVE('',#15,#16,#18,.T.);
#22 = EDGE_CURVE('',#16,#17,#9,.F.);
#23 = EDGE_CURVE('',#17,#14,#19,.T.);
#24 = ORIENTED_EDGE('',*,*,#20,.T.);
#25 = ORIENTED_EDGE('',*,*,#21,.T.);
#26 = ORIENTED_EDGE('',*,*,#22,.T.);
#27 = ORIENTED_EDGE('',*,*,#23,.T.);
#28 = EDGE_LOOP('',(#24,#25,#26,#27));
#29 = FACE_OUTER_BOUND('',#28,.T.);
#30 = ADVANCED_FACE('',(#29),#7,.T.);
#31 = CLOSED_SHELL('',(#30));
#32 = MANIFOLD_SOLID_BREP('cyl_patch',#31);
ENDSEC;
END-ISO-10303-21;
";
        let meshes = expect_ok(read_step(text.as_bytes()));
        assert_eq!(meshes.len(), 1);
        let mesh = &meshes[0];

        // Every vertex must lie on the cylinder x^2 + y^2 = 4.
        for v in &mesh.vertices {
            let r = v[0].hypot(v[1]);
            assert!((r - 2.0).abs() < 1e-9, "vertex {v:?} not on cylinder R=2");
            assert!(v[2] >= -1e-9 && v[2] <= 3.0 + 1e-9);
        }

        let area = surface_area(mesh);
        let expected = 3.0 * std::f64::consts::PI;
        assert!(
            (area - expected).abs() < 0.1,
            "cylinder patch area was {area}, expected ~{expected}"
        );
    }

    #[test]
    fn b_spline_surface_trimmed_patch_evaluates_interior_curvature() {
        // Bilinear/biquadratic dome: 3x3 control net with center control point
        // raised to z = 2.0, so S(0.5, 0.5) has z = 0.5 * 0.5 * 2.0 = 0.5.
        let text = "\
ISO-10303-21;
DATA;
#1 = CARTESIAN_POINT('',(0.,0.,0.));
#2 = CARTESIAN_POINT('',(0.,0.5,0.));
#3 = CARTESIAN_POINT('',(0.,1.,0.));
#4 = CARTESIAN_POINT('',(0.5,0.,0.));
#5 = CARTESIAN_POINT('',(0.5,0.5,2.));
#6 = CARTESIAN_POINT('',(0.5,1.,0.));
#7 = CARTESIAN_POINT('',(1.,0.,0.));
#8 = CARTESIAN_POINT('',(1.,0.5,0.));
#9 = CARTESIAN_POINT('',(1.,1.,0.));
#10 = B_SPLINE_SURFACE_WITH_KNOTS('',2,2,((#1,#2,#3),(#4,#5,#6),(#7,#8,#9)),.UNSPECIFIED.,.F.,.F.,.F.,(3,3),(3,3),(0.,1.),(0.,1.),.UNSPECIFIED.);
#11 = VERTEX_POINT('',#1);
#12 = VERTEX_POINT('',#7);
#13 = VERTEX_POINT('',#9);
#14 = VERTEX_POINT('',#3);
#15 = LINE('',#1,#7);
#16 = EDGE_CURVE('',#11,#12,#15,.T.);
#17 = EDGE_CURVE('',#12,#13,#15,.T.);
#18 = EDGE_CURVE('',#13,#14,#15,.T.);
#19 = EDGE_CURVE('',#14,#11,#15,.T.);
#20 = ORIENTED_EDGE('',*,*,#16,.T.);
#21 = ORIENTED_EDGE('',*,*,#17,.T.);
#22 = ORIENTED_EDGE('',*,*,#18,.T.);
#23 = ORIENTED_EDGE('',*,*,#19,.T.);
#24 = EDGE_LOOP('',(#20,#21,#22,#23));
#25 = FACE_OUTER_BOUND('',#24,.T.);
#26 = ADVANCED_FACE('',(#25),#10,.T.);
#27 = CLOSED_SHELL('',(#26));
#28 = MANIFOLD_SOLID_BREP('dome',#27);
ENDSEC;
END-ISO-10303-21;
";
        let meshes = expect_ok(read_step(text.as_bytes()));
        assert_eq!(meshes.len(), 1);
        let mesh = &meshes[0];
        assert!(
            mesh.vertices.len() > 20,
            "expected interior surface grid vertices"
        );
        let max_z = mesh
            .vertices
            .iter()
            .map(|v| v[2])
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(
            (max_z - 0.5).abs() < 0.02,
            "b-spline dome apex z was {max_z}, expected ~0.5"
        );
        // Surface area of curved dome must strictly exceed flat unit square (1.0).
        assert!(surface_area(mesh) > 1.05);
    }
}
