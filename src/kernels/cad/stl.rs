//! STL import and export.
//!
//! [`read_stl`] detects the binary and ASCII encodings by file size first and
//! by content second, so a binary file whose 80-byte header happens to start
//! with `solid` is still read as binary. Both encodings produce a triangle soup
//! with three fresh vertices per triangle; run [`super::repair::weld`] to build
//! shared connectivity.

use std::fmt::Write as _;

use super::types::Mesh;
use super::types::triangle_normal;

/// The size of the binary STL header in bytes.
const HEADER_LEN: usize = 80;

/// The size of one binary STL facet record in bytes.
const FACET_LEN: usize = 50;

/// Reads an STL file, detecting binary and ASCII encodings automatically.
///
/// # Errors
///
/// Returns `Err` when the data is too short to be an STL file, when a binary
/// triangle count does not match the payload length, when the ASCII body has a
/// malformed `vertex` line, or when the ASCII vertex count is not a multiple of
/// three.
pub fn read_stl(bytes: &[u8]) -> Result<Mesh, String> {
    if bytes.len() < 15 {
        return Err("stl data is too short".to_string());
    }

    if let Some(count) = binary_triangle_count(bytes) {
        return parse_binary(bytes, count);
    }

    if looks_like_ascii(bytes) {
        return parse_ascii(bytes);
    }

    Err("stl data is neither a complete binary file nor ascii".to_string())
}

/// Serialises `mesh` as a binary STL file.
///
/// Triangles whose indices are out of range are written with zero coordinates
/// rather than skipped, so the declared facet count always matches the payload.
#[must_use]
pub fn write_stl_binary(mesh: &Mesh) -> Vec<u8> {
    let count = u32::try_from(mesh.triangles.len()).unwrap_or(u32::MAX);
    let mut out = Vec::with_capacity(HEADER_LEN + 4 + FACET_LEN * mesh.triangles.len());

    let mut header = [0_u8; HEADER_LEN];
    let label = b"rssn binary stl";
    let take = label.len().min(HEADER_LEN);
    header[..take].copy_from_slice(&label[..take]);
    out.extend_from_slice(&header);
    out.extend_from_slice(&count.to_le_bytes());

    for t in &mesh.triangles {
        let normal = triangle_normal(mesh, *t).unwrap_or([0.0; 3]);
        for value in normal {
            out.extend_from_slice(&(value as f32).to_le_bytes());
        }
        for &index in t {
            let vertex = mesh.vertices.get(index).copied().unwrap_or([0.0; 3]);
            for value in vertex {
                out.extend_from_slice(&(value as f32).to_le_bytes());
            }
        }
        out.extend_from_slice(&0_u16.to_le_bytes());
    }

    out
}

/// Serialises `mesh` as an ASCII STL file with the given solid `name`.
#[must_use]
pub fn write_stl_ascii(
    mesh: &Mesh,
    name: &str,
) -> String {
    let mut out = String::with_capacity(64 + 200 * mesh.triangles.len());
    out.push_str("solid ");
    out.push_str(name);
    out.push('\n');

    for t in &mesh.triangles {
        let normal = triangle_normal(mesh, *t).unwrap_or([0.0; 3]);
        let _ = write!(
            out,
            "  facet normal {:.6} {:.6} {:.6}\n    outer loop\n",
            normal[0], normal[1], normal[2]
        );
        for &index in t {
            let v = mesh.vertices.get(index).copied().unwrap_or([0.0; 3]);
            let _ = writeln!(out, "      vertex {:.6} {:.6} {:.6}", v[0], v[1], v[2]);
        }
        out.push_str("    endloop\n  endfacet\n");
    }

    out.push_str("endsolid ");
    out.push_str(name);
    out.push('\n');
    out
}

/// The triangle count declared by a binary STL header, if the size matches.
fn binary_triangle_count(bytes: &[u8]) -> Option<usize> {
    let count = usize::try_from(read_u32_le(bytes, HEADER_LEN)?).ok()?;
    let expected = HEADER_LEN
        .checked_add(4)?
        .checked_add(count.checked_mul(FACET_LEN)?)?;
    if expected == bytes.len() {
        Some(count)
    } else {
        None
    }
}

/// Whether the payload looks like an ASCII STL file.
fn looks_like_ascii(bytes: &[u8]) -> bool {
    let head = &bytes[..bytes.len().min(2048)];
    let text = String::from_utf8_lossy(head).to_ascii_lowercase();
    let trimmed = text.trim_start();
    trimmed.starts_with("solid") || text.contains("facet") || text.contains("vertex")
}

/// Parses a binary STL payload of `count` facets.
fn parse_binary(
    bytes: &[u8],
    count: usize,
) -> Result<Mesh, String> {
    let mut vertices = Vec::with_capacity(count * 3);
    let mut triangles = Vec::with_capacity(count);

    for facet in 0..count {
        let base = HEADER_LEN + 4 + facet * FACET_LEN;
        for corner in 0..3 {
            let offset = base + 12 + corner * 12;
            let x = read_f32_le(bytes, offset)
                .ok_or_else(|| "truncated binary stl facet".to_string())?;
            let y = read_f32_le(bytes, offset + 4)
                .ok_or_else(|| "truncated binary stl facet".to_string())?;
            let z = read_f32_le(bytes, offset + 8)
                .ok_or_else(|| "truncated binary stl facet".to_string())?;
            vertices.push([f64::from(x), f64::from(y), f64::from(z)]);
        }
        let first = facet * 3;
        triangles.push([first, first + 1, first + 2]);
    }

    Ok(Mesh::new(vertices, triangles))
}

/// Parses an ASCII STL payload.
fn parse_ascii(bytes: &[u8]) -> Result<Mesh, String> {
    let text =
        std::str::from_utf8(bytes).map_err(|_| "ascii stl is not valid utf-8".to_string())?;

    let mut vertices: Vec<[f64; 3]> = Vec::new();
    for line in text.lines() {
        let mut words = line.split_whitespace();
        match words.next() {
            | Some(word) if word.eq_ignore_ascii_case("vertex") => {},
            | _ => continue,
        }

        let mut read =
            || -> Option<f64> { words.next().and_then(|token| token.parse::<f64>().ok()) };
        let Some(x) = read() else {
            return Err("malformed vertex line in ascii stl".to_string());
        };
        let Some(y) = read() else {
            return Err("malformed vertex line in ascii stl".to_string());
        };
        let Some(z) = read() else {
            return Err("malformed vertex line in ascii stl".to_string());
        };
        vertices.push([x, y, z]);
    }

    if vertices.is_empty() {
        return Err("ascii stl contains no vertices".to_string());
    }
    if !vertices.len().is_multiple_of(3) {
        return Err("ascii stl vertex count is not a multiple of three".to_string());
    }

    let triangles = (0..vertices.len() / 3)
        .map(|i| [i * 3, i * 3 + 1, i * 3 + 2])
        .collect();

    Ok(Mesh::new(vertices, triangles))
}

/// Reads a little-endian `u32` at `offset`.
fn read_u32_le(
    bytes: &[u8],
    offset: usize,
) -> Option<u32> {
    let slice = bytes.get(offset..offset + 4)?;
    let array: [u8; 4] = slice.try_into().ok()?;
    Some(u32::from_le_bytes(array))
}

/// Reads a little-endian `f32` at `offset`.
fn read_f32_le(
    bytes: &[u8],
    offset: usize,
) -> Option<f32> {
    let slice = bytes.get(offset..offset + 4)?;
    let array: [u8; 4] = slice.try_into().ok()?;
    Some(f32::from_le_bytes(array))
}

#[cfg(test)]
mod tests {
    use super::read_stl;
    use super::write_stl_ascii;
    use super::write_stl_binary;
    use crate::kernels::cad::repair::signed_volume;
    use crate::kernels::cad::repair::weld;
    use crate::kernels::cad::types::Mesh;

    fn expect_ok<T>(result: Result<T, String>) -> T {
        match result {
            | Ok(value) => value,
            | Err(message) => panic!("unexpected error: {message}"),
        }
    }

    fn tetrahedron() -> Mesh {
        Mesh::new(
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ],
            vec![[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]],
        )
    }

    #[test]
    fn binary_round_trip_preserves_the_solid() {
        let mesh = tetrahedron();
        let bytes = write_stl_binary(&mesh);
        let parsed = expect_ok(read_stl(&bytes));

        assert_eq!(parsed.triangles.len(), mesh.triangles.len());
        assert_eq!(parsed.vertices.len(), mesh.triangles.len() * 3);
        assert!((signed_volume(&parsed) - 1.0 / 6.0).abs() < 1e-6);

        let (welded, merged) = weld(&parsed, 1e-9);
        assert_eq!(welded.vertices.len(), 4);
        assert_eq!(merged, 8);
        assert!((signed_volume(&welded) - 1.0 / 6.0).abs() < 1e-6);
    }

    #[test]
    fn ascii_stl_is_parsed() {
        let text = "\
solid sample
  facet normal 0 0 1
    outer loop
      vertex 0 0 0
      vertex 1 0 0
      vertex 0 1 0
    endloop
  endfacet
  facet normal 0 0 1
    outer loop
      vertex 1 0 0
      vertex 1 1 0
      vertex 0 1 0
    endloop
  endfacet
endsolid sample
";
        let mesh = expect_ok(read_stl(text.as_bytes()));
        assert_eq!(mesh.triangles.len(), 2);
        assert_eq!(mesh.vertices.len(), 6);
        assert_eq!(mesh.vertices[1], [1.0, 0.0, 0.0]);
        assert_eq!(mesh.vertices[5], [0.0, 1.0, 0.0]);
    }

    #[test]
    fn ascii_export_round_trips() {
        let mesh = tetrahedron();
        let text = write_stl_ascii(&mesh, "tetra");
        assert!(text.starts_with("solid tetra"));
        assert!(text.trim_end().ends_with("endsolid tetra"));
        let parsed = expect_ok(read_stl(text.as_bytes()));
        assert_eq!(parsed.triangles.len(), 4);
    }

    #[test]
    fn truncated_binary_file_is_rejected() {
        let mut bytes = write_stl_binary(&tetrahedron());
        bytes.truncate(bytes.len() - 10);
        assert!(read_stl(&bytes).is_err());
    }

    #[test]
    fn short_and_malformed_inputs_are_rejected() {
        assert!(read_stl(b"").is_err());
        assert!(read_stl(b"solid empty\nendsolid empty\n").is_err());
        assert!(read_stl(b"solid bad\nvertex 1 2\nendsolid\n").is_err());
        assert!(read_stl(b"solid bad\nvertex 1 2 3\nendsolid\n").is_err());
    }
}
