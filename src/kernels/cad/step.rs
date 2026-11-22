//! A subset reader for ISO-10303-21 (STEP) boundary-representation files.
//!
//! The reader understands the part of the AP203/AP214 schema that describes a
//! polygonal solid:
//!
//! `CARTESIAN_POINT`, `VERTEX_POINT`, `EDGE_CURVE`, `ORIENTED_EDGE`,
//! `EDGE_LOOP`, `FACE_OUTER_BOUND`, `FACE_BOUND`, `ADVANCED_FACE`,
//! `CLOSED_SHELL` and `MANIFOLD_SOLID_BREP`.
//!
//! Each face loop is tessellated by a triangle fan from its first vertex, which
//! is exact for the planar polygonal faces the subset targets. Curves are not
//! evaluated: an `EDGE_CURVE` contributes only its two vertices.
//!
//! # Limitations
//!
//! * Only line-like edges are supported. A file that uses
//!   `B_SPLINE_CURVE_WITH_KNOTS` is rejected.
//! * Trimmed surfaces, `FACE_BOUND` holes and non-planar faces are not
//!   honoured; a face is always the polygon of its outer loop.
//! * Units and transformations are ignored, so coordinates are returned in the
//!   file's own unit.

use std::collections::HashMap;
use std::path::Path;

use super::types::Mesh;

/// Reads every solid in a STEP file.
///
/// Returns one [`Mesh`] per `MANIFOLD_SOLID_BREP`, or one per `CLOSED_SHELL`
/// when the file has no solid wrapper.
///
/// # Errors
///
/// Returns `Err` when the bytes are not UTF-8, when there is no `DATA` section,
/// when the file uses a B-spline curve, or when an entity reference is
/// undefined or malformed.
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
    let mut vertex_points: HashMap<usize, usize> = HashMap::new();
    let mut edge_curves: HashMap<usize, (usize, usize)> = HashMap::new();
    let mut oriented_edges: HashMap<usize, (usize, bool)> = HashMap::new();
    let mut edge_loops: HashMap<usize, Vec<usize>> = HashMap::new();
    let mut face_bounds: HashMap<usize, usize> = HashMap::new();
    let mut advanced_faces: HashMap<usize, (Vec<usize>, bool)> = HashMap::new();
    let mut closed_shells: HashMap<usize, Vec<usize>> = HashMap::new();
    let mut solid_breps: Vec<usize> = Vec::new();

    for record in &records[start + 1..end] {
        let Some((id, kind, args)) = parse_record(record) else {
            continue;
        };
        match kind.as_str() {
            | "CARTESIAN_POINT" => {
                if let Some(point) = parse_point(&args) {
                    points.insert(id, point);
                }
            },
            | "VERTEX_POINT" => {
                if let Some(point_id) = args.get(1).and_then(|a| parse_reference(a)) {
                    vertex_points.insert(id, point_id);
                }
            },
            | "EDGE_CURVE" => {
                let first = args.get(1).and_then(|a| parse_reference(a));
                let second = args.get(2).and_then(|a| parse_reference(a));
                if let (Some(first), Some(second)) = (first, second) {
                    edge_curves.insert(id, (first, second));
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
            | "FACE_OUTER_BOUND" | "FACE_BOUND" => {
                if let Some(loop_id) = args.get(1).and_then(|a| parse_reference(a)) {
                    face_bounds.insert(id, loop_id);
                }
            },
            | "ADVANCED_FACE" => {
                let bounds = args
                    .get(1)
                    .and_then(|a| parse_reference_list(a))
                    .unwrap_or_default();
                let sense = args.last().is_none_or(|a| parse_boolean(a));
                advanced_faces.insert(id, (bounds, sense));
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
            | "B_SPLINE_CURVE_WITH_KNOTS" => {
                return Err("b-spline curves are not supported".to_string());
            },
            | _ => {},
        }
    }

    let context = Context {
        points: &points,
        vertex_points: &vertex_points,
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

/// The entity tables a shell needs while it is being tessellated.
struct Context<'a> {
    points: &'a HashMap<usize, [f64; 3]>,
    vertex_points: &'a HashMap<usize, usize>,
    edge_curves: &'a HashMap<usize, (usize, usize)>,
    oriented_edges: &'a HashMap<usize, (usize, bool)>,
    edge_loops: &'a HashMap<usize, Vec<usize>>,
    face_bounds: &'a HashMap<usize, usize>,
    advanced_faces: &'a HashMap<usize, (Vec<usize>, bool)>,
}

/// Tessellates every face of `faces` into one mesh.
fn build_shell(
    context: &Context<'_>,
    faces: &[usize],
) -> Result<Mesh, String> {
    let mut vertices: Vec<[f64; 3]> = Vec::new();
    let mut triangles: Vec<[usize; 3]> = Vec::new();
    let mut seen: HashMap<usize, usize> = HashMap::new();

    for face_id in faces {
        let (bounds, sense) = context
            .advanced_faces
            .get(face_id)
            .ok_or_else(|| format!("advanced face #{face_id} is undefined"))?;
        let Some(bound_id) = bounds.first() else {
            continue;
        };
        let loop_id = context
            .face_bounds
            .get(bound_id)
            .ok_or_else(|| format!("face bound #{bound_id} is undefined"))?;
        let oriented = context
            .edge_loops
            .get(loop_id)
            .ok_or_else(|| format!("edge loop #{loop_id} is undefined"))?;

        let mut ring: Vec<usize> = Vec::new();
        for oriented_id in oriented {
            let (edge_id, forward) = context
                .oriented_edges
                .get(oriented_id)
                .ok_or_else(|| format!("oriented edge #{oriented_id} is undefined"))?;
            let (first, second) = context
                .edge_curves
                .get(edge_id)
                .ok_or_else(|| format!("edge curve #{edge_id} is undefined"))?;
            let from = if *forward { *first } else { *second };

            let point_id = context
                .vertex_points
                .get(&from)
                .ok_or_else(|| format!("vertex point #{from} is undefined"))?;
            let coordinate = context
                .points
                .get(point_id)
                .ok_or_else(|| format!("cartesian point #{point_id} is undefined"))?;

            let index = match seen.get(point_id) {
                | Some(&index) => index,
                | None => {
                    let index = vertices.len();
                    vertices.push(*coordinate);
                    seen.insert(*point_id, index);
                    index
                },
            };
            if ring.last() != Some(&index) {
                ring.push(index);
            }
        }

        if ring.len() >= 2 && ring[0] == ring[ring.len() - 1] {
            ring.pop();
        }
        if ring.len() < 3 {
            continue;
        }
        if !*sense {
            ring.reverse();
        }

        for i in 1..ring.len() - 1 {
            triangles.push([ring[0], ring[i], ring[i + 1]]);
        }
    }

    if triangles.is_empty() {
        return Err("step shell has no tessellatable faces".to_string());
    }

    Ok(Mesh::new(vertices, triangles))
}

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

/// Splits `#12 = TYPE(a, b, c)` into its id, type name and raw arguments.
fn parse_record(record: &str) -> Option<(usize, String, Vec<String>)> {
    let equals = record.find('=')?;
    let id = record[..equals]
        .trim()
        .strip_prefix('#')?
        .trim()
        .parse()
        .ok()?;
    let rest = record[equals + 1..].trim();

    let open = rest.find('(')?;
    let close = rest.rfind(')')?;
    if close <= open {
        return None;
    }
    let kind = rest[..open].trim().to_ascii_uppercase();
    if kind.is_empty() {
        return None;
    }

    Some((id, kind, split_arguments(&rest[open + 1..close])))
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

/// Reads the coordinates of a `CARTESIAN_POINT` argument list.
fn parse_point(args: &[String]) -> Option<[f64; 3]> {
    let list = args.get(1)?;
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

/// Reads a STEP boolean, where `.T.` is true and anything else is false.
fn parse_boolean(raw: &str) -> bool {
    raw.trim().eq_ignore_ascii_case(".T.")
}

#[cfg(test)]
mod tests {
    use super::read_step;

    fn expect_ok<T>(result: Result<T, String>) -> T {
        match result {
            | Ok(value) => value,
            | Err(message) => panic!("unexpected error: {message}"),
        }
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
    fn b_spline_curves_are_rejected() {
        let text = "\
DATA;
#1 = B_SPLINE_CURVE_WITH_KNOTS('',3,(#2,#3),.UNSPECIFIED.,.F.,.F.,(4,4),(0.,1.),.UNSPECIFIED.);
ENDSEC;
";
        let error = match read_step(text.as_bytes()) {
            | Ok(_) => panic!("expected an error"),
            | Err(message) => message,
        };
        assert_eq!(error, "b-spline curves are not supported");
    }

    #[test]
    fn file_reader_reports_a_missing_path() {
        assert!(super::read_step_file(std::path::Path::new("does-not-exist.step")).is_err());
    }
}
