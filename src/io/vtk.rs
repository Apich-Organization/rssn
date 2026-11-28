//! Complete VTK Legacy (ASCII & Binary) post-processing reader and writer.
//!
//! Supports `STRUCTURED_POINTS`, `RECTILINEAR_GRID`, `POLYDATA`, and
//! `UNSTRUCTURED_GRID` datasets with `POINT_DATA` and `CELL_DATA`
//! (`SCALARS`, `VECTORS`, `TENSORS`), in both `ASCII` and big-endian `BINARY`
//! encodings per the VTK Legacy File Format specification.

use std::fmt::Write as _;
use std::path::Path;

use ndarray::Array2;

/// VTK Legacy encoding format.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum VtkFormat {
    /// Human-readable ASCII text.
    #[default]
    Ascii,
    /// Big-endian IEEE binary payload.
    Binary,
}

/// Attribute array attached to `POINT_DATA` or `CELL_DATA`.
#[derive(Debug, Clone, PartialEq)]
pub enum VtkAttributeData {
    /// Scalar field (`1` component per point/cell).
    Scalars {
        /// Attribute name.
        name: String,
        /// Scalar values.
        values: Vec<f64>,
    },
    /// 3D vector field (`3` components per point/cell).
    Vectors {
        /// Attribute name.
        name: String,
        /// Vector values.
        values: Vec<[f64; 3]>,
    },
    /// 3x3 tensor field (`9` components per point/cell).
    Tensors {
        /// Attribute name.
        name: String,
        /// Row-major 3x3 tensor values.
        values: Vec<[f64; 9]>,
    },
}

impl VtkAttributeData {
    const fn len(&self) -> usize {
        match self {
            | Self::Scalars { values, .. } => values.len(),
            | Self::Vectors { values, .. } => values.len(),
            | Self::Tensors { values, .. } => values.len(),
        }
    }
}

/// VTK dataset geometry and topology.
#[derive(Debug, Clone, PartialEq)]
pub enum VtkDataset {
    /// Regular axis-aligned grid.
    StructuredPoints {
        /// `[nx, ny, nz]`.
        dims: [usize; 3],
        /// Grid origin `[ox, oy, oz]`.
        origin: [f64; 3],
        /// Grid spacing `[sx, sy, sz]`.
        spacing: [f64; 3],
    },
    /// Non-uniform axis-aligned grid.
    RectilinearGrid {
        /// X coordinates.
        xs: Vec<f64>,
        /// Y coordinates.
        ys: Vec<f64>,
        /// Z coordinates.
        zs: Vec<f64>,
    },
    /// Polygonal / particle data.
    PolyData {
        /// 3D vertex coordinates.
        points: Vec<[f64; 3]>,
        /// Vertex cells.
        vertices: Vec<Vec<usize>>,
        /// Polyline cells.
        lines: Vec<Vec<usize>>,
        /// Polygon cells.
        polygons: Vec<Vec<usize>>,
    },
    /// Unstructured finite-element mesh.
    UnstructuredGrid {
        /// 3D vertex coordinates.
        points: Vec<[f64; 3]>,
        /// Cell connectivity lists.
        cells: Vec<Vec<usize>>,
        /// VTK cell type IDs (`5` = tri, `9` = quad, `10` = tet, `12` = hex).
        cell_types: Vec<u8>,
    },
}

impl VtkDataset {
    /// Number of points in the dataset.
    #[must_use]
    pub const fn num_points(&self) -> usize {
        match self {
            | Self::StructuredPoints { dims, .. } => {
                dims[0].saturating_mul(dims[1]).saturating_mul(dims[2])
            },
            | Self::RectilinearGrid { xs, ys, zs } => {
                xs.len().saturating_mul(ys.len()).saturating_mul(zs.len())
            },
            | Self::PolyData { points, .. } | Self::UnstructuredGrid { points, .. } => points.len(),
        }
    }

    /// Number of cells in the dataset.
    #[must_use]
    pub fn num_cells(&self) -> usize {
        match self {
            | Self::StructuredPoints { dims, .. } => {
                dims[0].saturating_sub(1).max(1)
                    * dims[1].saturating_sub(1).max(1)
                    * dims[2].saturating_sub(1).max(1)
            },
            | Self::RectilinearGrid { xs, ys, zs } => {
                xs.len().saturating_sub(1).max(1)
                    * ys.len().saturating_sub(1).max(1)
                    * zs.len().saturating_sub(1).max(1)
            },
            | Self::PolyData {
                vertices,
                lines,
                polygons,
                ..
            } => vertices.len() + lines.len() + polygons.len(),
            | Self::UnstructuredGrid { cells, .. } => cells.len(),
        }
    }
}

/// Complete VTK Legacy file representation.
#[derive(Debug, Clone, PartialEq)]
pub struct VtkFile {
    /// Header title line (max 256 chars, single line).
    pub title: String,
    /// ASCII or Binary encoding.
    pub format: VtkFormat,
    /// Dataset geometry/topology.
    pub dataset: VtkDataset,
    /// Point-centered attributes.
    pub point_data: Vec<VtkAttributeData>,
    /// Cell-centered attributes.
    pub cell_data: Vec<VtkAttributeData>,
}

fn write_f64_slice(
    out: &mut Vec<u8>,
    format: VtkFormat,
    vals: &[f64],
) {
    match format {
        | VtkFormat::Ascii => {
            let mut s = String::new();
            for &v in vals {
                let _ = writeln!(s, "{v}");
            }
            out.extend_from_slice(s.as_bytes());
        },
        | VtkFormat::Binary => {
            for &v in vals {
                out.extend_from_slice(&v.to_be_bytes());
            }
            out.push(b'\n');
        },
    }
}

fn write_cell_list(
    out: &mut Vec<u8>,
    format: VtkFormat,
    keyword: &str,
    cells: &[Vec<usize>],
    num_points: usize,
) -> Result<(), String> {
    if cells.is_empty() {
        return Ok(());
    }
    let total_ints: usize = cells.iter().map(|c| 1 + c.len()).sum();
    for cell in cells {
        for &idx in cell {
            if idx >= num_points {
                return Err(format!(
                    "VTK cell references point {idx} >= num_points {num_points}"
                ));
            }
        }
    }
    out.extend_from_slice(format!("{keyword} {} {total_ints}\n", cells.len()).as_bytes());
    match format {
        | VtkFormat::Ascii => {
            let mut s = String::new();
            for cell in cells {
                let _ = write!(s, "{}", cell.len());
                for &idx in cell {
                    let _ = write!(s, " {idx}");
                }
                s.push('\n');
            }
            out.extend_from_slice(s.as_bytes());
        },
        | VtkFormat::Binary => {
            for cell in cells {
                out.extend_from_slice(&(cell.len() as i32).to_be_bytes());
                for &idx in cell {
                    out.extend_from_slice(&(idx as i32).to_be_bytes());
                }
            }
            out.push(b'\n');
        },
    }
    Ok(())
}

fn write_attributes(
    out: &mut Vec<u8>,
    format: VtkFormat,
    section: &str,
    expected_len: usize,
    attrs: &[VtkAttributeData],
) -> Result<(), String> {
    if attrs.is_empty() {
        return Ok(());
    }
    out.extend_from_slice(format!("{section} {expected_len}\n").as_bytes());
    for attr in attrs {
        if attr.len() != expected_len {
            return Err(format!(
                "VTK {section} attribute length {} != expected {expected_len}",
                attr.len()
            ));
        }
        match attr {
            | VtkAttributeData::Scalars { name, values } => {
                out.extend_from_slice(format!("SCALARS {name} double 1\n").as_bytes());
                out.extend_from_slice(b"LOOKUP_TABLE default\n");
                write_f64_slice(out, format, values);
            },
            | VtkAttributeData::Vectors { name, values } => {
                out.extend_from_slice(format!("VECTORS {name} double\n").as_bytes());
                let flat: Vec<f64> = values.iter().flat_map(|v| *v).collect();
                write_f64_slice(out, format, &flat);
            },
            | VtkAttributeData::Tensors { name, values } => {
                out.extend_from_slice(format!("TENSORS {name} double\n").as_bytes());
                let flat: Vec<f64> = values.iter().flat_map(|v| *v).collect();
                write_f64_slice(out, format, &flat);
            },
        }
    }
    Ok(())
}

/// Encodes a [`VtkFile`] into VTK Legacy v3.0 bytes (`ASCII` or `BINARY`).
///
/// # Errors
/// Fails when geometry dimensions are invalid or attribute counts mismatch.
pub fn encode_vtk(file: &VtkFile) -> Result<Vec<u8>, String> {
    let title = file.title.replace(['\r', '\n'], " ");
    let title = if title.is_empty() {
        "RSSN VTK Dataset"
    } else {
        title.as_str()
    };
    let mut out = Vec::new();
    out.extend_from_slice(b"# vtk DataFile Version 3.0\n");
    out.extend_from_slice(title.as_bytes());
    out.push(b'\n');
    match file.format {
        | VtkFormat::Ascii => out.extend_from_slice(b"ASCII\n"),
        | VtkFormat::Binary => out.extend_from_slice(b"BINARY\n"),
    }

    match &file.dataset {
        | VtkDataset::StructuredPoints {
            dims,
            origin,
            spacing,
        } => {
            if dims.contains(&0) {
                return Err("VTK STRUCTURED_POINTS dimensions must be positive".to_owned());
            }
            out.extend_from_slice(b"DATASET STRUCTURED_POINTS\n");
            out.extend_from_slice(
                format!("DIMENSIONS {} {} {}\n", dims[0], dims[1], dims[2]).as_bytes(),
            );
            out.extend_from_slice(
                format!("ORIGIN {} {} {}\n", origin[0], origin[1], origin[2]).as_bytes(),
            );
            out.extend_from_slice(
                format!("SPACING {} {} {}\n", spacing[0], spacing[1], spacing[2]).as_bytes(),
            );
        },
        | VtkDataset::RectilinearGrid { xs, ys, zs } => {
            if xs.is_empty() || ys.is_empty() || zs.is_empty() {
                return Err("VTK RECTILINEAR_GRID axes must be non-empty".to_owned());
            }
            out.extend_from_slice(b"DATASET RECTILINEAR_GRID\n");
            out.extend_from_slice(
                format!("DIMENSIONS {} {} {}\n", xs.len(), ys.len(), zs.len()).as_bytes(),
            );
            out.extend_from_slice(format!("X_COORDINATES {} double\n", xs.len()).as_bytes());
            write_f64_slice(&mut out, file.format, xs);
            out.extend_from_slice(format!("Y_COORDINATES {} double\n", ys.len()).as_bytes());
            write_f64_slice(&mut out, file.format, ys);
            out.extend_from_slice(format!("Z_COORDINATES {} double\n", zs.len()).as_bytes());
            write_f64_slice(&mut out, file.format, zs);
        },
        | VtkDataset::PolyData {
            points,
            vertices,
            lines,
            polygons,
        } => {
            out.extend_from_slice(b"DATASET POLYDATA\n");
            out.extend_from_slice(format!("POINTS {} double\n", points.len()).as_bytes());
            let flat: Vec<f64> = points.iter().flat_map(|p| *p).collect();
            write_f64_slice(&mut out, file.format, &flat);
            write_cell_list(&mut out, file.format, "VERTICES", vertices, points.len())?;
            write_cell_list(&mut out, file.format, "LINES", lines, points.len())?;
            write_cell_list(&mut out, file.format, "POLYGONS", polygons, points.len())?;
        },
        | VtkDataset::UnstructuredGrid {
            points,
            cells,
            cell_types,
        } => {
            if cells.len() != cell_types.len() {
                return Err("VTK cells and cell_types length mismatch".to_owned());
            }
            out.extend_from_slice(b"DATASET UNSTRUCTURED_GRID\n");
            out.extend_from_slice(format!("POINTS {} double\n", points.len()).as_bytes());
            let flat: Vec<f64> = points.iter().flat_map(|p| *p).collect();
            write_f64_slice(&mut out, file.format, &flat);
            write_cell_list(&mut out, file.format, "CELLS", cells, points.len())?;
            out.extend_from_slice(format!("CELL_TYPES {}\n", cell_types.len()).as_bytes());
            match file.format {
                | VtkFormat::Ascii => {
                    let mut s = String::new();
                    for &ct in cell_types {
                        let _ = writeln!(s, "{ct}");
                    }
                    out.extend_from_slice(s.as_bytes());
                },
                | VtkFormat::Binary => {
                    for &ct in cell_types {
                        out.extend_from_slice(&i32::from(ct).to_be_bytes());
                    }
                    out.push(b'\n');
                },
            }
        },
    }

    write_attributes(
        &mut out,
        file.format,
        "POINT_DATA",
        file.dataset.num_points(),
        &file.point_data,
    )?;
    write_attributes(
        &mut out,
        file.format,
        "CELL_DATA",
        file.dataset.num_cells(),
        &file.cell_data,
    )?;

    Ok(out)
}

/// Writes a [`VtkFile`] to `path`.
///
/// # Errors
/// As [`encode_vtk`], plus file I/O errors.
pub fn write_vtk<P: AsRef<Path>>(
    path: P,
    file: &VtkFile,
) -> Result<(), String> {
    let bytes = encode_vtk(file)?;
    std::fs::write(path.as_ref(), bytes)
        .map_err(|e| format!("cannot write {}: {e}", path.as_ref().display()))
}

// -------------------------------------------------------------------------
// Parser (ASCII & Binary)
// -------------------------------------------------------------------------

struct VtkParser<'a> {
    buf: &'a [u8],
    pos: usize,
}

impl<'a> VtkParser<'a> {
    const fn new(buf: &'a [u8]) -> Self {
        Self { buf, pos: 0 }
    }

    fn read_line(&mut self) -> Result<String, String> {
        if self.pos >= self.buf.len() {
            return Err("unexpected EOF in VTK header".to_owned());
        }
        let start = self.pos;
        while self.pos < self.buf.len() && self.buf[self.pos] != b'\n' {
            self.pos += 1;
        }
        let end = if self.pos > start && self.buf[self.pos - 1] == b'\r' {
            self.pos - 1
        } else {
            self.pos
        };
        if self.pos < self.buf.len() && self.buf[self.pos] == b'\n' {
            self.pos += 1;
        }
        String::from_utf8(self.buf[start..end].to_vec())
            .map_err(|_| "VTK header line is not valid UTF-8".to_owned())
    }

    const fn skip_ws(&mut self) {
        while self.pos < self.buf.len() && self.buf[self.pos].is_ascii_whitespace() {
            self.pos += 1;
        }
    }

    fn next_token(&mut self) -> Option<String> {
        self.skip_ws();
        if self.pos >= self.buf.len() {
            return None;
        }
        let start = self.pos;
        while self.pos < self.buf.len() && !self.buf[self.pos].is_ascii_whitespace() {
            self.pos += 1;
        }
        Some(String::from_utf8_lossy(&self.buf[start..self.pos]).into_owned())
    }

    fn read_scalars(
        &mut self,
        format: VtkFormat,
        dtype: &str,
        count: usize,
    ) -> Result<Vec<f64>, String> {
        let mut out = Vec::with_capacity(count);
        match format {
            | VtkFormat::Ascii => {
                for _ in 0..count {
                    let tok = self
                        .next_token()
                        .ok_or_else(|| "VTK ASCII numeric data truncated".to_owned())?;
                    let v = tok
                        .parse::<f64>()
                        .map_err(|_| format!("invalid VTK number `{tok}`"))?;
                    out.push(v);
                }
            },
            | VtkFormat::Binary => {
                let width = match dtype {
                    | "double" => 8,
                    | "float" | "int" | "unsigned_int" => 4,
                    | "short" | "unsigned_short" => 2,
                    | "char" | "unsigned_char" => 1,
                    | other => return Err(format!("unsupported VTK binary type `{other}`")),
                };
                let needed = count
                    .checked_mul(width)
                    .ok_or_else(|| "VTK binary payload overflows".to_owned())?;
                let slice = self
                    .buf
                    .get(self.pos..self.pos + needed)
                    .ok_or_else(|| "VTK binary payload truncated".to_owned())?;
                self.pos += needed;
                for chunk in slice.chunks_exact(width) {
                    let v = match dtype {
                        | "double" => f64::from_be_bytes(chunk.try_into().unwrap()),
                        | "float" => f64::from(f32::from_be_bytes(chunk.try_into().unwrap())),
                        | "int" => f64::from(i32::from_be_bytes(chunk.try_into().unwrap())),
                        | "unsigned_int" => {
                            f64::from(u32::from_be_bytes(chunk.try_into().unwrap()))
                        },
                        | "short" => f64::from(i16::from_be_bytes(chunk.try_into().unwrap())),
                        | "unsigned_short" => {
                            f64::from(u16::from_be_bytes(chunk.try_into().unwrap()))
                        },
                        | "char" => f64::from(chunk[0] as i8),
                        | "unsigned_char" => f64::from(chunk[0]),
                        | _ => unreachable!(),
                    };
                    out.push(v);
                }
                if self.buf.get(self.pos) == Some(&b'\n') {
                    self.pos += 1;
                }
            },
        }
        Ok(out)
    }

    fn read_ints(
        &mut self,
        format: VtkFormat,
        count: usize,
    ) -> Result<Vec<i32>, String> {
        let mut out = Vec::with_capacity(count);
        match format {
            | VtkFormat::Ascii => {
                for _ in 0..count {
                    let tok = self
                        .next_token()
                        .ok_or_else(|| "VTK ASCII integer data truncated".to_owned())?;
                    let v = tok
                        .parse::<i32>()
                        .map_err(|_| format!("invalid VTK integer `{tok}`"))?;
                    out.push(v);
                }
            },
            | VtkFormat::Binary => {
                let needed = count * 4;
                let slice = self
                    .buf
                    .get(self.pos..self.pos + needed)
                    .ok_or_else(|| "VTK binary integer payload truncated".to_owned())?;
                self.pos += needed;
                for c in slice.as_chunks::<4>().0 {
                    out.push(i32::from_be_bytes(*c));
                }
                if self.buf.get(self.pos) == Some(&b'\n') {
                    self.pos += 1;
                }
            },
        }
        Ok(out)
    }
}

fn decode_cell_list(
    raw: &[i32],
    num_cells: usize,
) -> Result<Vec<Vec<usize>>, String> {
    let mut cells = Vec::with_capacity(num_cells);
    let mut pos = 0_usize;
    for _ in 0..num_cells {
        let n = usize::try_from(*raw.get(pos).ok_or("VTK cell list truncated")?)
            .map_err(|_| "negative VTK cell vertex count".to_owned())?;
        pos += 1;
        let slice = raw
            .get(pos..pos + n)
            .ok_or_else(|| "VTK cell vertices truncated".to_owned())?;
        pos += n;
        let mut cell = Vec::with_capacity(n);
        for &idx in slice {
            cell.push(usize::try_from(idx).map_err(|_| "negative VTK vertex index".to_owned())?);
        }
        cells.push(cell);
    }
    Ok(cells)
}

/// Parses a complete VTK Legacy file (ASCII or Binary) with strict validation.
///
/// # Errors
/// Fails on invalid header, unsupported encoding/dataset, or truncated/mismatched data.
pub fn parse_vtk(bytes: &[u8]) -> Result<VtkFile, String> {
    let mut p = VtkParser::new(bytes);
    let magic = p.read_line()?;
    if !magic.starts_with("# vtk DataFile Version") {
        return Err("not a VTK Legacy file (missing # vtk DataFile Version header)".to_owned());
    }
    let title = p.read_line()?;
    let fmt_line = p.read_line()?;
    let format = match fmt_line.trim() {
        | "ASCII" => VtkFormat::Ascii,
        | "BINARY" => VtkFormat::Binary,
        | other => return Err(format!("unsupported VTK format `{other}`")),
    };

    let ds_tok = p.next_token().ok_or("missing DATASET keyword")?;
    if ds_tok != "DATASET" {
        return Err(format!("expected DATASET, got `{ds_tok}`"));
    }
    let ds_kind = p.next_token().ok_or("missing dataset type")?;

    let mut dims = [0_usize; 3];
    let mut origin = [0.0_f64; 3];
    let mut spacing = [1.0_f64; 3];
    let mut xs = Vec::new();
    let mut ys = Vec::new();
    let mut zs = Vec::new();
    let mut points = Vec::new();
    let mut vertices = Vec::new();
    let mut lines = Vec::new();
    let mut polygons = Vec::new();
    let mut cells = Vec::new();
    let mut cell_types = Vec::new();
    let mut point_data = Vec::new();
    let mut cell_data = Vec::new();
    let mut active_section: Option<(bool, usize)> = None; // (is_point_data, count)

    while let Some(tok) = p.next_token() {
        match tok.as_str() {
            | "DIMENSIONS" => {
                for d in &mut dims {
                    let s = p.next_token().ok_or("truncated DIMENSIONS")?;
                    *d = s.parse().map_err(|_| "invalid DIMENSIONS".to_owned())?;
                }
            },
            | "ORIGIN" => {
                for o in &mut origin {
                    let s = p.next_token().ok_or("truncated ORIGIN")?;
                    *o = s.parse().map_err(|_| "invalid ORIGIN".to_owned())?;
                }
            },
            | "SPACING" => {
                for sp in &mut spacing {
                    let s = p.next_token().ok_or("truncated SPACING")?;
                    *sp = s.parse().map_err(|_| "invalid SPACING".to_owned())?;
                }
            },
            | "X_COORDINATES" | "Y_COORDINATES" | "Z_COORDINATES" => {
                let n: usize = p
                    .next_token()
                    .ok_or("missing coordinate count")?
                    .parse()
                    .map_err(|_| "invalid coordinate count".to_owned())?;
                let dtype = p.next_token().ok_or("missing coordinate dtype")?;
                // Consume newline before binary block
                if format == VtkFormat::Binary && p.buf.get(p.pos) == Some(&b'\n') {
                    p.pos += 1;
                }
                let vals = p.read_scalars(format, &dtype, n)?;
                match tok.as_str() {
                    | "X_COORDINATES" => xs = vals,
                    | "Y_COORDINATES" => ys = vals,
                    | _ => zs = vals,
                }
            },
            | "POINTS" => {
                let n: usize = p
                    .next_token()
                    .ok_or("missing POINTS count")?
                    .parse()
                    .map_err(|_| "invalid POINTS count".to_owned())?;
                let dtype = p.next_token().ok_or("missing POINTS dtype")?;
                if format == VtkFormat::Binary && p.buf.get(p.pos) == Some(&b'\n') {
                    p.pos += 1;
                }
                let flat = p.read_scalars(format, &dtype, n * 3)?;
                points = flat.as_chunks::<3>().0.to_vec();
            },
            | "VERTICES" | "LINES" | "POLYGONS" | "CELLS" => {
                let n: usize = p
                    .next_token()
                    .ok_or("missing cell count")?
                    .parse()
                    .map_err(|_| "invalid cell count".to_owned())?;
                let total_ints: usize = p
                    .next_token()
                    .ok_or("missing cell total ints")?
                    .parse()
                    .map_err(|_| "invalid cell total ints".to_owned())?;
                if format == VtkFormat::Binary && p.buf.get(p.pos) == Some(&b'\n') {
                    p.pos += 1;
                }
                let raw = p.read_ints(format, total_ints)?;
                let decoded = decode_cell_list(&raw, n)?;
                match tok.as_str() {
                    | "VERTICES" => vertices = decoded,
                    | "LINES" => lines = decoded,
                    | "POLYGONS" => polygons = decoded,
                    | _ => cells = decoded,
                }
            },
            | "CELL_TYPES" => {
                let n: usize = p
                    .next_token()
                    .ok_or("missing CELL_TYPES count")?
                    .parse()
                    .map_err(|_| "invalid CELL_TYPES count".to_owned())?;
                if format == VtkFormat::Binary && p.buf.get(p.pos) == Some(&b'\n') {
                    p.pos += 1;
                }
                let raw = p.read_ints(format, n)?;
                cell_types = raw
                    .into_iter()
                    .map(|v| u8::try_from(v).map_err(|_| "invalid VTK cell type".to_owned()))
                    .collect::<Result<Vec<_>, _>>()?;
            },
            | "POINT_DATA" => {
                let n: usize = p
                    .next_token()
                    .ok_or("missing POINT_DATA count")?
                    .parse()
                    .map_err(|_| "invalid POINT_DATA count".to_owned())?;
                active_section = Some((true, n));
            },
            | "CELL_DATA" => {
                let n: usize = p
                    .next_token()
                    .ok_or("missing CELL_DATA count")?
                    .parse()
                    .map_err(|_| "invalid CELL_DATA count".to_owned())?;
                active_section = Some((false, n));
            },
            | "SCALARS" => {
                let (is_pt, count) =
                    active_section.ok_or("SCALARS outside POINT_DATA/CELL_DATA")?;
                let name = p.next_token().ok_or("missing SCALARS name")?;
                let dtype = p.next_token().ok_or("missing SCALARS dtype")?;
                // Optional numComp + LOOKUP_TABLE
                let next = p.next_token().ok_or("missing LOOKUP_TABLE")?;
                let ncomp: usize = if next == "LOOKUP_TABLE" {
                    let _table = p.next_token().ok_or("missing LOOKUP_TABLE name")?;
                    1
                } else {
                    let nc = next
                        .parse()
                        .map_err(|_| "invalid SCALARS numComp".to_owned())?;
                    let lt = p.next_token().ok_or("missing LOOKUP_TABLE")?;
                    if lt != "LOOKUP_TABLE" {
                        return Err(format!("expected LOOKUP_TABLE, got `{lt}`"));
                    }
                    let _table = p.next_token().ok_or("missing LOOKUP_TABLE name")?;
                    nc
                };
                if ncomp != 1 {
                    return Err(format!("only 1-component SCALARS supported, got {ncomp}"));
                }
                if format == VtkFormat::Binary && p.buf.get(p.pos) == Some(&b'\n') {
                    p.pos += 1;
                }
                let values = p.read_scalars(format, &dtype, count)?;
                let attr = VtkAttributeData::Scalars { name, values };
                if is_pt {
                    point_data.push(attr);
                } else {
                    cell_data.push(attr);
                }
            },
            | "VECTORS" => {
                let (is_pt, count) =
                    active_section.ok_or("VECTORS outside POINT_DATA/CELL_DATA")?;
                let name = p.next_token().ok_or("missing VECTORS name")?;
                let dtype = p.next_token().ok_or("missing VECTORS dtype")?;
                if format == VtkFormat::Binary && p.buf.get(p.pos) == Some(&b'\n') {
                    p.pos += 1;
                }
                let flat = p.read_scalars(format, &dtype, count * 3)?;
                let values = flat.as_chunks::<3>().0.to_vec();
                let attr = VtkAttributeData::Vectors { name, values };
                if is_pt {
                    point_data.push(attr);
                } else {
                    cell_data.push(attr);
                }
            },
            | "TENSORS" => {
                let (is_pt, count) =
                    active_section.ok_or("TENSORS outside POINT_DATA/CELL_DATA")?;
                let name = p.next_token().ok_or("missing TENSORS name")?;
                let dtype = p.next_token().ok_or("missing TENSORS dtype")?;
                if format == VtkFormat::Binary && p.buf.get(p.pos) == Some(&b'\n') {
                    p.pos += 1;
                }
                let flat = p.read_scalars(format, &dtype, count * 9)?;
                let values = flat.as_chunks::<9>().0.to_vec();
                let attr = VtkAttributeData::Tensors { name, values };
                if is_pt {
                    point_data.push(attr);
                } else {
                    cell_data.push(attr);
                }
            },
            | other => {
                return Err(format!("unexpected VTK keyword `{other}`"));
            },
        }
    }

    let dataset = match ds_kind.as_str() {
        | "STRUCTURED_POINTS" => {
            VtkDataset::StructuredPoints {
                dims,
                origin,
                spacing,
            }
        },
        | "RECTILINEAR_GRID" => VtkDataset::RectilinearGrid { xs, ys, zs },
        | "POLYDATA" => {
            VtkDataset::PolyData {
                points,
                vertices,
                lines,
                polygons,
            }
        },
        | "UNSTRUCTURED_GRID" => {
            VtkDataset::UnstructuredGrid {
                points,
                cells,
                cell_types,
            }
        },
        | other => return Err(format!("unsupported VTK DATASET `{other}`")),
    };

    Ok(VtkFile {
        title,
        format,
        dataset,
        point_data,
        cell_data,
    })
}

/// Reads and parses a VTK Legacy file from `path`.
///
/// # Errors
/// Fails when the file cannot be read or parsed.
pub fn read_vtk<P: AsRef<Path>>(path: P) -> Result<VtkFile, String> {
    let bytes = std::fs::read(path.as_ref())
        .map_err(|e| format!("cannot read {}: {e}", path.as_ref().display()))?;
    parse_vtk(&bytes)
}

// -------------------------------------------------------------------------
// Convenience wrappers (preserving existing public API)
// -------------------------------------------------------------------------

/// Writes a 2D scalar field as VTK Legacy ASCII `STRUCTURED_POINTS`.
///
/// # Errors
/// Fails when `data` is empty or the file cannot be written.
pub fn write_vtk_structured_points_2d<P: AsRef<Path>>(
    path: P,
    data: &Array2<f64>,
    field_name: &str,
    spacing: [f64; 2],
    origin: [f64; 2],
) -> Result<(), String> {
    let (rows, cols) = data.dim();
    if rows == 0 || cols == 0 {
        return Err("VTK 2D field cannot be empty".to_owned());
    }
    write_vtk(
        path,
        &VtkFile {
            title: "RSSN 2D Scalar Field".to_owned(),
            format: VtkFormat::Ascii,
            dataset: VtkDataset::StructuredPoints {
                dims: [cols, rows, 1],
                origin: [origin[0], origin[1], 0.0],
                spacing: [spacing[0], spacing[1], 1.0],
            },
            point_data: vec![VtkAttributeData::Scalars {
                name: field_name.to_owned(),
                values: data.iter().copied().collect(),
            }],
            cell_data: Vec::new(),
        },
    )
}

/// Writes a 3D scalar field (row-major `nx * ny * nz`) as VTK Legacy ASCII `STRUCTURED_POINTS`.
///
/// # Errors
/// Fails when dimensions are zero, `data.len() != nx * ny * nz`, or writing fails.
pub fn write_vtk_structured_points_3d<P: AsRef<Path>>(
    path: P,
    dims: [usize; 3],
    data: &[f64],
    field_name: &str,
    spacing: [f64; 3],
    origin: [f64; 3],
) -> Result<(), String> {
    let [nx, ny, nz] = dims;
    let total = nx
        .checked_mul(ny)
        .and_then(|v| v.checked_mul(nz))
        .ok_or_else(|| "VTK 3D dimensions overflow".to_owned())?;
    if total == 0 || data.len() != total {
        return Err(format!(
            "VTK 3D field expected {total} values for {nx}x{ny}x{nz}, got {}",
            data.len()
        ));
    }
    let mut ordered = Vec::with_capacity(total);
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                ordered.push(data[i * (ny * nz) + j * nz + k]);
            }
        }
    }
    write_vtk(
        path,
        &VtkFile {
            title: "RSSN 3D Scalar Field".to_owned(),
            format: VtkFormat::Ascii,
            dataset: VtkDataset::StructuredPoints {
                dims,
                origin,
                spacing,
            },
            point_data: vec![VtkAttributeData::Scalars {
                name: field_name.to_owned(),
                values: ordered,
            }],
            cell_data: Vec::new(),
        },
    )
}

/// Writes a 2D/3D rectilinear grid with non-uniform axis coordinates and scalar values.
///
/// # Errors
/// Fails when any axis is empty, `data.len() != xs.len() * ys.len() * zs.len()`, or writing fails.
pub fn write_vtk_rectilinear_grid<P: AsRef<Path>>(
    path: P,
    xs: &[f64],
    ys: &[f64],
    zs: &[f64],
    data: &[f64],
    field_name: &str,
) -> Result<(), String> {
    write_vtk(
        path,
        &VtkFile {
            title: "RSSN Rectilinear Grid".to_owned(),
            format: VtkFormat::Ascii,
            dataset: VtkDataset::RectilinearGrid {
                xs: xs.to_vec(),
                ys: ys.to_vec(),
                zs: zs.to_vec(),
            },
            point_data: vec![VtkAttributeData::Scalars {
                name: field_name.to_owned(),
                values: data.to_vec(),
            }],
            cell_data: Vec::new(),
        },
    )
}

/// Writes 3D particle positions (with optional per-particle scalar values) as `POLYDATA`.
///
/// # Errors
/// Fails when `scalars` length does not match `points.len()` or writing fails.
pub fn write_vtk_polydata_points<P: AsRef<Path>>(
    path: P,
    points: &[[f64; 3]],
    scalars: Option<(&str, &[f64])>,
) -> Result<(), String> {
    let vertices: Vec<Vec<usize>> = (0..points.len()).map(|i| vec![i]).collect();
    let point_data = match scalars {
        | Some((name, vals)) => {
            vec![VtkAttributeData::Scalars {
                name: name.to_owned(),
                values: vals.to_vec(),
            }]
        },
        | None => Vec::new(),
    };
    write_vtk(
        path,
        &VtkFile {
            title: "RSSN Particle Cloud".to_owned(),
            format: VtkFormat::Ascii,
            dataset: VtkDataset::PolyData {
                points: points.to_vec(),
                vertices,
                lines: Vec::new(),
                polygons: Vec::new(),
            },
            point_data,
            cell_data: Vec::new(),
        },
    )
}

/// Writes an unstructured mesh (`UNSTRUCTURED_GRID`) with VTK cell types
/// (e.g. `5` = triangle, `9` = quad, `10` = tetrahedron, `12` = hexahedron).
///
/// # Errors
/// Fails when `cells.len() != cell_types.len()` or writing fails.
pub fn write_vtk_unstructured_grid<P: AsRef<Path>>(
    path: P,
    points: &[[f64; 3]],
    cells: &[Vec<usize>],
    cell_types: &[u8],
    point_scalars: Option<(&str, &[f64])>,
) -> Result<(), String> {
    let point_data = match point_scalars {
        | Some((name, vals)) => {
            vec![VtkAttributeData::Scalars {
                name: name.to_owned(),
                values: vals.to_vec(),
            }]
        },
        | None => Vec::new(),
    };
    write_vtk(
        path,
        &VtkFile {
            title: "RSSN Unstructured Mesh".to_owned(),
            format: VtkFormat::Ascii,
            dataset: VtkDataset::UnstructuredGrid {
                points: points.to_vec(),
                cells: cells.to_vec(),
                cell_types: cell_types.to_vec(),
            },
            point_data,
            cell_data: Vec::new(),
        },
    )
}

/// Reads a 2D scalar field back from a VTK Legacy (ASCII or Binary) `STRUCTURED_POINTS` file.
///
/// # Errors
/// Fails when the file is not a valid 2D `STRUCTURED_POINTS` file.
pub fn read_vtk_structured_points_2d<P: AsRef<Path>>(path: P) -> Result<Array2<f64>, String> {
    let vtk = read_vtk(path)?;
    let VtkDataset::StructuredPoints { dims, .. } = vtk.dataset else {
        return Err("VTK file is not STRUCTURED_POINTS".to_owned());
    };
    let Some(VtkAttributeData::Scalars { values, .. }) = vtk.point_data.first() else {
        return Err("VTK file missing scalar POINT_DATA".to_owned());
    };
    Array2::from_shape_vec((dims[1], dims[0]), values.clone()).map_err(|e| e.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn vtk_ascii_and_binary_round_trip() {
        for format in [VtkFormat::Ascii, VtkFormat::Binary] {
            let vtk = VtkFile {
                title: "Test".to_owned(),
                format,
                dataset: VtkDataset::StructuredPoints {
                    dims: [3, 2, 1],
                    origin: [0.0, 0.0, 0.0],
                    spacing: [0.5, 0.5, 1.0],
                },
                point_data: vec![
                    VtkAttributeData::Scalars {
                        name: "u".to_owned(),
                        values: vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
                    },
                    VtkAttributeData::Vectors {
                        name: "vel".to_owned(),
                        values: vec![[1.0, 0.0, 0.0]; 6],
                    },
                ],
                cell_data: Vec::new(),
            };
            let bytes = encode_vtk(&vtk).unwrap();
            let back = parse_vtk(&bytes).unwrap();
            assert_eq!(back, vtk);
        }
    }
}
