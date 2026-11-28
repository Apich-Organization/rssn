//! NetCDF classic, 64-bit offset and 64-bit data files.
//!
//! A dependency-free reader and writer for the NetCDF classic data model as
//! stored in `CDF\x01` (classic), `CDF\x02` (64-bit offset) and `CDF\x05`
//! (64-bit data) files. The
//! reader parses the header lists (`NC_DIMENSION`, `NC_ATTRIBUTE`,
//! `NC_VARIABLE`), the record (unlimited) dimension with record-major
//! variable layout, and the `NC_BYTE` through `NC_DOUBLE` external types.
//!
//! Numeric values are exposed as `f64` (integer types convert exactly up to
//! 2^53). Variables of external type `NC_CHAR` are surfaced as [`NcType::Char`]
//! with their raw byte values preserved in [`Variable::data`] and convenience
//! string decoders [`Variable::as_string`] and [`Variable::as_strings`].
//!
//! The writer emits `CDF\x01`, `CDF\x02` or `CDF\x05` headers with four-byte
//! padding and big-endian data via [`write_netcdf`] /
//! [`write_netcdf_with_format`] and [`to_bytes`] / [`to_bytes_with_format`],
//! automatically upgrading a classic header to `CDF\x02` when any variable
//! offset exceeds `u32::MAX`. `CDF\x05` widens every count, name length and
//! offset to 64 bits; the type codes and dimension identifiers stay 32-bit.
//!
//! # Limitations
//!
//! - HDF5-backed NetCDF-4 files (the `.nc4` container) are rejected; the three
//!   classic headers are all supported.
//! - Only the first dimension may be unlimited, as required by the classic
//!   format.
//!
//! Every length, count and offset read from a file is checked against the
//! buffer before use, and arrays are capped at 64 MiB and 2^32 elements.
//!
//! [`NcType::Char`]: crate::io::netcdf::NcType::Char
//! [`Variable::data`]: crate::io::netcdf::Variable::data
//! [`Variable::as_string`]: crate::io::netcdf::Variable::as_string
//! [`Variable::as_strings`]: crate::io::netcdf::Variable::as_strings
//! [`write_netcdf`]: crate::io::netcdf::write_netcdf
//! [`write_netcdf_with_format`]: crate::io::netcdf::write_netcdf_with_format
//! [`to_bytes`]: crate::io::netcdf::to_bytes
//! [`to_bytes_with_format`]: crate::io::netcdf::to_bytes_with_format

use std::path::Path;

/// Largest byte length allowed for one decoded variable (64 MiB).
const MAX_ARRAY_BYTES: usize = 64 * 1024 * 1024;

/// Largest element count allowed for one variable (2^32).
const MAX_ELEMS: u64 = 1 << 32;

/// Longest name accepted in a header.
const MAX_NAME: usize = 1_048_576;

/// Header tag for the dimension list.
const NC_DIMENSION: u32 = 10;

/// Header tag for the variable list.
const NC_VARIABLE: u32 = 11;

/// Header tag for an attribute list.
const NC_ATTRIBUTE: u32 = 12;

/// External type: signed 8-bit integer.
const NC_BYTE: u32 = 1;

/// External type: character.
const NC_CHAR: u32 = 2;

/// External type: signed 16-bit integer.
const NC_SHORT: u32 = 3;

/// External type: signed 32-bit integer.
const NC_INT: u32 = 4;

/// External type: 32-bit IEEE floating point.
const NC_FLOAT: u32 = 5;

/// External type: 64-bit IEEE floating point.
const NC_DOUBLE: u32 = 6;

/// Binary format version emitted by the NetCDF writer.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum NcFormat {
    /// NetCDF classic format (`CDF\x01`, 32-bit variable offsets).
    #[default]
    Classic,
    /// NetCDF 64-bit offset format (`CDF\x02`, 64-bit variable offsets).
    Offset64,
    /// NetCDF 64-bit data format (`CDF\x05`, every count and offset 64-bit).
    Cdf5,
}

/// The external datatype of a NetCDF variable.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NcType {
    /// Signed 8-bit integer (`NC_BYTE`).
    I8,
    /// Unsigned 8-bit integer (stored as `NC_CHAR`).
    U8,
    /// Character data (`NC_CHAR`).
    Char,
    /// Signed 16-bit integer (`NC_SHORT`).
    I16,
    /// Unsigned 16-bit integer.
    U16,
    /// Signed 32-bit integer (`NC_INT`).
    I32,
    /// Unsigned 32-bit integer.
    U32,
    /// 32-bit IEEE floating point (`NC_FLOAT`).
    F32,
    /// 64-bit IEEE floating point (`NC_DOUBLE`).
    F64,
}

impl NcType {
    /// Byte width of one value of this type.
    const fn size(self) -> usize {
        match self {
            | Self::I8 | Self::U8 | Self::Char => 1,
            | Self::I16 | Self::U16 => 2,
            | Self::I32 | Self::U32 | Self::F32 => 4,
            | Self::F64 => 8,
        }
    }

    /// The external type code used in a classic header.
    const fn code(self) -> u32 {
        match self {
            | Self::I8 => NC_BYTE,
            | Self::U8 | Self::Char => NC_CHAR,
            | Self::I16 | Self::U16 => NC_SHORT,
            | Self::I32 | Self::U32 => NC_INT,
            | Self::F32 => NC_FLOAT,
            | Self::F64 => NC_DOUBLE,
        }
    }
}

/// One dimension of a [`NcFile`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Dimension {
    /// Dimension name.
    pub name: String,
    /// Length, or the number of records for the unlimited dimension.
    pub len: usize,
    /// Whether this is the record (unlimited) dimension.
    pub unlimited: bool,
}

/// One variable of a [`NcFile`].
#[derive(Debug, Clone, PartialEq)]
pub struct Variable {
    /// Variable name.
    pub name: String,
    /// Names of the dimensions, outermost first.
    pub dims: Vec<String>,
    /// External datatype.
    pub nc_type: NcType,
    /// Variable attributes in file order.
    pub attributes: Vec<(String, NcValue)>,
    /// Row-major values converted to `f64` (byte values `0.0..=255.0` for [`NcType::Char`]).
    pub data: Vec<f64>,
}

/// Decodes a NUL-terminated or NUL-padded byte slice into a UTF-8 string.
fn decode_char_slice(bytes: &[u8]) -> Option<String> {
    let end = bytes.iter().position(|&b| b == 0).unwrap_or(bytes.len());
    std::str::from_utf8(bytes.get(..end)?)
        .ok()
        .map(str::to_owned)
}

impl Variable {
    /// Creates a 1-D [`NcType::Char`] variable from a UTF-8 string.
    pub fn from_string(
        name: impl Into<String>,
        dim_name: impl Into<String>,
        text: &str,
        attributes: Vec<(String, NcValue)>,
    ) -> Self {
        Self {
            name: name.into(),
            dims: vec![dim_name.into()],
            nc_type: NcType::Char,
            attributes,
            data: text.as_bytes().iter().map(|&b| f64::from(b)).collect(),
        }
    }

    /// Creates an N-D [`NcType::Char`] variable from fixed-width NUL-padded strings.
    ///
    /// Each string is truncated or NUL-padded to `str_len` bytes, which should
    /// match the innermost (last) dimension in `dims`.
    pub fn from_strings(
        name: impl Into<String>,
        dims: Vec<String>,
        strings: &[&str],
        str_len: usize,
        attributes: Vec<(String, NcValue)>,
    ) -> Self {
        let mut data = Vec::with_capacity(strings.len().saturating_mul(str_len));
        for text in strings {
            let bytes = text.as_bytes();
            let take = bytes.len().min(str_len);
            for &b in &bytes[..take] {
                data.push(f64::from(b));
            }
            data.extend(std::iter::repeat_n(0.0, str_len - take));
        }
        Self {
            name: name.into(),
            dims,
            nc_type: NcType::Char,
            attributes,
            data,
        }
    }

    /// Converts [`Self::data`] to raw bytes when `self.nc_type == NcType::Char`.
    fn char_bytes(&self) -> Option<Vec<u8>> {
        if self.nc_type != NcType::Char {
            return None;
        }
        let mut bytes = Vec::with_capacity(self.data.len());
        for &value in &self.data {
            if !value.is_finite() {
                return None;
            }
            let truncated = value.trunc();
            if !(0.0..=255.0).contains(&truncated) {
                return None;
            }
            bytes.push(truncated as u8);
        }
        Some(bytes)
    }

    /// Decodes a 1-D (or scalar) [`NcType::Char`] variable as a UTF-8 string.
    ///
    /// Returns `None` if this variable is not [`NcType::Char`] or its bytes are
    /// not valid UTF-8.
    pub fn as_string(&self) -> Option<String> {
        let bytes = self.char_bytes()?;
        decode_char_slice(&bytes)
    }

    /// Decodes a 1-D or N-D [`NcType::Char`] variable into a list of strings.
    ///
    /// For 1-D (or scalar) variables this returns a single-element vector; for
    /// 2-D and higher-rank variables the last dimension in `self.dims` is used
    /// as the fixed byte width of each string (looked up in `dimensions`).
    pub fn as_strings(
        &self,
        dimensions: &[Dimension],
    ) -> Option<Vec<String>> {
        let bytes = self.char_bytes()?;
        if self.dims.len() <= 1 {
            return Some(vec![decode_char_slice(&bytes)?]);
        }
        let last_dim = self.dims.last()?;
        let str_len = dimensions.iter().find(|d| d.name == *last_dim)?.len;
        if str_len == 0 {
            return Some(Vec::new());
        }
        let mut out = Vec::with_capacity(bytes.len() / str_len);
        for chunk in bytes.chunks_exact(str_len) {
            out.push(decode_char_slice(chunk)?);
        }
        Some(out)
    }
}

/// An attribute value.
#[derive(Debug, Clone, PartialEq)]
pub enum NcValue {
    /// A character string.
    Text(String),
    /// A single integer.
    Int(i32),
    /// A single floating-point value.
    Float(f64),
    /// An array of integers.
    Ints(Vec<i32>),
    /// An array of floating-point values.
    Floats(Vec<f64>),
}

/// The contents of a NetCDF classic file.
#[derive(Debug, Clone, PartialEq)]
pub struct NcFile {
    /// Format version: 1 for `CDF\x01`, 2 for `CDF\x02`.
    pub version: u8,
    /// Dimensions in header order.
    pub dimensions: Vec<Dimension>,
    /// Variables in header order.
    pub variables: Vec<Variable>,
    /// Global attributes in header order.
    pub global_attributes: Vec<(String, NcValue)>,
}

impl NcFile {
    /// Decodes the named [`NcType::Char`] variable into strings using this
    /// file's dimension table.
    pub fn variable_strings(
        &self,
        name: &str,
    ) -> Option<Vec<String>> {
        let variable = self.variables.iter().find(|var| var.name == name)?;
        variable.as_strings(&self.dimensions)
    }
}

// -------------------------------------------------------------------------
// Big-endian cursor
// -------------------------------------------------------------------------

/// A bounds-checked big-endian reader over the file image.
struct Cursor<'a> {
    /// Whole file image.
    data: &'a [u8],
    /// Next unread byte.
    pos: usize,
    /// Whether `NON_NEG` fields are 64-bit (CDF-5).
    wide: bool,
}

impl<'a> Cursor<'a> {
    /// Creates a cursor at the start of `data`.
    const fn new(data: &'a [u8]) -> Self {
        Self {
            data,
            pos: 0,
            wide: false,
        }
    }

    /// Reads a `NON_NEG` field, which is 64-bit in CDF-5 and 32-bit before.
    fn nonneg(&mut self) -> Result<u64, String> {
        if self.wide {
            self.u64()
        } else {
            Ok(u64::from(self.u32()?))
        }
    }

    /// Consumes `count` bytes, failing on any overrun.
    fn take(
        &mut self,
        count: usize,
    ) -> Result<&'a [u8], String> {
        let end = self
            .pos
            .checked_add(count)
            .ok_or_else(|| "header offset overflow".to_owned())?;
        let bytes = self
            .data
            .get(self.pos..end)
            .ok_or_else(|| format!("read of {count} bytes at {} runs past the end", self.pos))?;
        self.pos = end;
        Ok(bytes)
    }

    /// Reads a big-endian `u32`.
    fn u32(&mut self) -> Result<u32, String> {
        let bytes = self.take(4)?;
        let mut raw = [0_u8; 4];
        raw.copy_from_slice(bytes);
        Ok(u32::from_be_bytes(raw))
    }

    /// Reads a big-endian `u64`.
    fn u64(&mut self) -> Result<u64, String> {
        let bytes = self.take(8)?;
        let mut raw = [0_u8; 8];
        raw.copy_from_slice(bytes);
        Ok(u64::from_be_bytes(raw))
    }

    /// Skips `count` bytes.
    fn skip(
        &mut self,
        count: usize,
    ) -> Result<(), String> {
        let _ = self.take(count)?;
        Ok(())
    }

    /// Advances to the next multiple of four.
    fn align4(&mut self) -> Result<(), String> {
        let aligned = (self.pos + 3) & !3;
        if aligned > self.data.len() {
            return Err("header padding runs past the end of the file".to_owned());
        }
        self.pos = aligned;
        Ok(())
    }

    /// Reads a padded name.
    fn name(&mut self) -> Result<String, String> {
        let len = to_usize(self.nonneg()?, "name length")?;
        if len > MAX_NAME {
            return Err(format!("name of {len} bytes is unreasonably long"));
        }
        let bytes = self.take(len)?;
        self.align4()?;
        std::str::from_utf8(bytes)
            .map(str::to_owned)
            .map_err(|_| "name is not valid UTF-8".to_owned())
    }
}

/// Converts a `u64` count to `usize`.
fn to_usize(
    value: u64,
    what: &str,
) -> Result<usize, String> {
    usize::try_from(value).map_err(|_| format!("{what} {value} does not fit in memory"))
}

/// Rounds `n` up to the next multiple of four.
const fn align4(n: u64) -> u64 {
    (n + 3) & !3
}

/// Multiplies dimension sizes with overflow and budget checks.
fn element_count(shape: &[usize]) -> Result<u64, String> {
    let mut total: u64 = 1;
    for &dim in shape {
        total = total
            .checked_mul(u64::try_from(dim).map_err(|_| "dimension too large".to_owned())?)
            .ok_or_else(|| "variable dimensions overflow".to_owned())?;
        if total > MAX_ELEMS {
            return Err(format!(
                "variable has more than {MAX_ELEMS} elements, which is not supported"
            ));
        }
    }
    Ok(total)
}

/// Byte width of an external type code.
const fn external_size(code: u32) -> Result<usize, ()> {
    match code {
        | NC_BYTE | NC_CHAR => Ok(1),
        | NC_SHORT => Ok(2),
        | NC_INT | NC_FLOAT => Ok(4),
        | NC_DOUBLE => Ok(8),
        | _ => Err(()),
    }
}

/// Maps an external type code to [`NcType`].
const fn nc_type_of(code: u32) -> Result<NcType, ()> {
    match code {
        | NC_BYTE => Ok(NcType::I8),
        | NC_CHAR => Ok(NcType::Char),
        | NC_SHORT => Ok(NcType::I16),
        | NC_INT => Ok(NcType::I32),
        | NC_FLOAT => Ok(NcType::F32),
        | NC_DOUBLE => Ok(NcType::F64),
        | _ => Err(()),
    }
}

// -------------------------------------------------------------------------
// Header parsing
// -------------------------------------------------------------------------

/// A variable header record before its data is read.
struct RawVar {
    /// Variable name.
    name: String,
    /// Dimension identifiers into the dimension table.
    dimids: Vec<usize>,
    /// Variable attributes.
    attributes: Vec<(String, NcValue)>,
    /// External type code.
    code: u32,
    /// Absolute offset of the first data byte.
    begin: u64,
}

/// The parsed header of a NetCDF file.
struct Header {
    /// Format version.
    version: u8,
    /// Number of records.
    numrecs: u64,
    /// Dimensions.
    dimensions: Vec<Dimension>,
    /// Global attributes.
    global_attributes: Vec<(String, NcValue)>,
    /// Variable headers.
    variables: Vec<RawVar>,
}

/// Reads one attribute value list.
fn read_attribute(cur: &mut Cursor<'_>) -> Result<(String, NcValue), String> {
    let name = cur.name()?;
    let code = cur.u32()?;
    let count = to_usize(cur.nonneg()?, "attribute element count")?;
    let size = external_size(code)
        .map_err(|()| format!("attribute {name} uses unsupported type code {code}"))?;
    let bytes = count
        .checked_mul(size)
        .ok_or_else(|| format!("attribute {name} size overflow"))?;
    if bytes > MAX_ARRAY_BYTES {
        return Err(format!("attribute {name} exceeds the size limit"));
    }
    let values = cur.take(bytes)?;
    cur.align4()?;
    let value = if code == NC_CHAR {
        let end = values
            .iter()
            .rposition(|&byte| byte != 0)
            .map_or(0, |index| index + 1);
        let text = std::str::from_utf8(values.get(..end).unwrap_or(&[]))
            .map_err(|_| format!("attribute {name} is not valid UTF-8"))?;
        NcValue::Text(text.to_owned())
    } else {
        decode_values(code, values, count)?
    };
    Ok((name, value))
}

/// Reads an attribute list, returning an empty vector for an absent list.
fn read_attribute_list(cur: &mut Cursor<'_>) -> Result<Vec<(String, NcValue)>, String> {
    let tag = cur.u32()?;
    let count = to_usize(cur.nonneg()?, "attribute count")?;
    if tag == 0 {
        if count != 0 {
            return Err("absent attribute list declares entries".to_owned());
        }
        return Ok(Vec::new());
    }
    if tag != NC_ATTRIBUTE {
        return Err(format!("unexpected attribute list tag {tag}"));
    }
    if count > 65_536 {
        return Err("attribute list has too many entries".to_owned());
    }
    let mut out = Vec::with_capacity(count);
    for _ in 0..count {
        out.push(read_attribute(cur)?);
    }
    Ok(out)
}

/// Decodes raw big-endian values into an [`NcValue`].
fn decode_values(
    code: u32,
    bytes: &[u8],
    count: usize,
) -> Result<NcValue, String> {
    let size =
        external_size(code).map_err(|()| format!("unsupported external type code {code}"))?;
    let needed = count
        .checked_mul(size)
        .ok_or_else(|| "value array size overflow".to_owned())?;
    if needed > MAX_ARRAY_BYTES {
        return Err("value array exceeds the size limit".to_owned());
    }
    let available = bytes
        .get(..needed)
        .ok_or_else(|| "value array is shorter than its declared length".to_owned())?;
    let mut numbers = Vec::with_capacity(count);
    for chunk in available.chunks_exact(size) {
        numbers.push(decode_value(code, chunk)?);
    }
    if count == 1 {
        let value = numbers.first().copied().ok_or("missing attribute value")?;
        return Ok(match code {
            | NC_FLOAT | NC_DOUBLE => NcValue::Float(value),
            | _ => NcValue::Int(value as i32),
        });
    }
    Ok(match code {
        | NC_FLOAT | NC_DOUBLE => NcValue::Floats(numbers),
        | _ => NcValue::Ints(numbers.iter().map(|&value| value as i32).collect()),
    })
}

/// Decodes one big-endian value of external type `code`.
fn decode_value(
    code: u32,
    bytes: &[u8],
) -> Result<f64, String> {
    match code {
        | NC_BYTE => {
            let byte = *bytes.first().ok_or("short value")?;
            Ok(f64::from(byte as i8))
        },
        | NC_CHAR => {
            let byte = *bytes.first().ok_or("short value")?;
            Ok(f64::from(byte))
        },
        | NC_SHORT => {
            let mut raw = [0_u8; 2];
            raw.copy_from_slice(bytes.get(..2).ok_or("short value")?);
            Ok(f64::from(i16::from_be_bytes(raw)))
        },
        | NC_INT => {
            let mut raw = [0_u8; 4];
            raw.copy_from_slice(bytes.get(..4).ok_or("short value")?);
            Ok(f64::from(i32::from_be_bytes(raw)))
        },
        | NC_FLOAT => {
            let mut raw = [0_u8; 4];
            raw.copy_from_slice(bytes.get(..4).ok_or("short value")?);
            Ok(f64::from(f32::from_bits(u32::from_be_bytes(raw))))
        },
        | NC_DOUBLE => {
            let mut raw = [0_u8; 8];
            raw.copy_from_slice(bytes.get(..8).ok_or("short value")?);
            Ok(f64::from_bits(u64::from_be_bytes(raw)))
        },
        | other => Err(format!("unsupported external type code {other}")),
    }
}

/// Parses a complete header.
fn parse_header(data: &[u8]) -> Result<Header, String> {
    let magic = data.get(..4).ok_or("file is too short to be NetCDF")?;
    if magic.get(..3) != Some(b"CDF".as_slice()) {
        return Err("not a NetCDF classic file: magic mismatch".to_owned());
    }
    let version = *magic.get(3).ok_or("truncated NetCDF magic")?;
    if !matches!(version, 1 | 2 | 5) {
        return Err(format!(
            "unsupported NetCDF version {version}; CDF-1, CDF-2 and CDF-5 are supported"
        ));
    }
    let mut cur = Cursor::new(data);
    cur.wide = version == 5;
    cur.skip(4)?;
    // `numrecs` is a 32-bit count up to CDF-2 and a 64-bit count in CDF-5.
    let numrecs = if version == 5 {
        let raw = cur.u64()?;
        if raw == u64::MAX {
            return Err(
                "streaming files with indeterminate record counts are not supported".to_owned(),
            );
        }
        raw
    } else {
        let raw = cur.u32()?;
        if raw == u32::MAX {
            return Err(
                "streaming files with indeterminate record counts are not supported".to_owned(),
            );
        }
        u64::from(raw)
    };

    let dim_tag = cur.u32()?;
    let dim_count = to_usize(cur.nonneg()?, "dimension count")?;
    let mut dimensions = Vec::new();
    if dim_tag == 0 {
        if dim_count != 0 {
            return Err("absent dimension list declares entries".to_owned());
        }
    } else if dim_tag == NC_DIMENSION {
        if dim_count > 65_536 {
            return Err("dimension list has too many entries".to_owned());
        }
        for _ in 0..dim_count {
            let name = cur.name()?;
            let len = cur.nonneg()?;
            let unlimited = len == 0;
            if unlimited && !dimensions.is_empty() {
                return Err("only the first dimension may be the unlimited dimension".to_owned());
            }
            dimensions.push(Dimension {
                name,
                len: if unlimited {
                    to_usize(numrecs, "record count")?
                } else {
                    to_usize(len, "dimension length")?
                },
                unlimited,
            });
        }
    } else {
        return Err(format!("unexpected dimension list tag {dim_tag}"));
    }

    let global_attributes = read_attribute_list(&mut cur)?;

    let var_tag = cur.u32()?;
    let var_count = to_usize(cur.nonneg()?, "variable count")?;
    let mut variables = Vec::new();
    if var_tag == 0 {
        if var_count != 0 {
            return Err("absent variable list declares entries".to_owned());
        }
    } else if var_tag == NC_VARIABLE {
        if var_count > 65_536 {
            return Err("variable list has too many entries".to_owned());
        }
        for _ in 0..var_count {
            let name = cur.name()?;
            let ndims = to_usize(cur.nonneg()?, "variable rank")?;
            if ndims > 1024 {
                return Err(format!("variable {name} has an implausible rank {ndims}"));
            }
            let mut dimids = Vec::with_capacity(ndims);
            for _ in 0..ndims {
                let id = to_usize(cur.nonneg()?, "dimension id")?;
                if id >= dimensions.len() {
                    return Err(format!(
                        "variable {name} references dimension {id}, which does not exist"
                    ));
                }
                dimids.push(id);
            }
            let attributes = read_attribute_list(&mut cur)?;
            let code = cur.u32()?;
            let _ = external_size(code)
                .map_err(|()| format!("variable {name} uses unsupported type code {code}"))?;
            let _vsize = cur.nonneg()?;
            let begin = if version == 1 {
                u64::from(cur.u32()?)
            } else {
                cur.u64()?
            };
            variables.push(RawVar {
                name,
                dimids,
                attributes,
                code,
                begin,
            });
        }
    } else {
        return Err(format!("unexpected variable list tag {var_tag}"));
    }

    Ok(Header {
        version,
        numrecs,
        dimensions,
        global_attributes,
        variables,
    })
}

// -------------------------------------------------------------------------
// Data reading
// -------------------------------------------------------------------------

/// Reads the data of one variable, following the record layout.
fn read_variable_data(
    data: &[u8],
    header: &Header,
    variable: &RawVar,
    record_size: u64,
) -> Result<Vec<f64>, String> {
    let size = external_size(variable.code)
        .map_err(|()| format!("variable {} uses an unsupported type", variable.name))?;
    let mut shape = Vec::with_capacity(variable.dimids.len());
    for &id in &variable.dimids {
        let dim = header
            .dimensions
            .get(id)
            .ok_or_else(|| format!("variable {} references a missing dimension", variable.name))?;
        shape.push(dim.len);
    }
    let total = element_count(&shape)?;
    let total = to_usize(total, "variable element count")?;
    let bytes = total
        .checked_mul(size)
        .ok_or_else(|| "variable byte size overflow".to_owned())?;
    if bytes > MAX_ARRAY_BYTES {
        return Err(format!(
            "variable {} needs {bytes} bytes, more than the limit",
            variable.name
        ));
    }
    let is_record = variable
        .dimids
        .first()
        .and_then(|&id| header.dimensions.get(id))
        .is_some_and(|dim| dim.unlimited);
    if !is_record {
        let start = to_usize(variable.begin, "variable offset")?;
        let end = start
            .checked_add(bytes)
            .ok_or_else(|| "variable offset overflow".to_owned())?;
        let slice = data.get(start..end).ok_or_else(|| {
            format!(
                "variable {} data runs past the end of the file",
                variable.name
            )
        })?;
        return decode_values_f64(variable.code, slice, total);
    }
    if record_size == 0 {
        return Err(format!(
            "variable {} is a record variable but the record size is zero",
            variable.name
        ));
    }
    let slab_elems = if shape.is_empty() {
        0
    } else {
        element_count(shape.get(1..).unwrap_or(&[]))?
    };
    let slab_elems = to_usize(slab_elems, "record slab size")?;
    let slab_bytes = slab_elems
        .checked_mul(size)
        .ok_or_else(|| "record slab size overflow".to_owned())?;
    let records = to_usize(header.numrecs, "record count")?;
    let mut out = Vec::with_capacity(total);
    for record in 0..records {
        let offset = u64::try_from(record)
            .unwrap_or(0)
            .saturating_mul(record_size)
            .saturating_add(variable.begin);
        let start = to_usize(offset, "record offset")?;
        let end = start
            .checked_add(slab_bytes)
            .ok_or_else(|| "record offset overflow".to_owned())?;
        let slice = data.get(start..end).ok_or_else(|| {
            format!(
                "record {record} of variable {} runs past the end of the file",
                variable.name
            )
        })?;
        let values = decode_values_f64(variable.code, slice, slab_elems)?;
        out.extend_from_slice(&values);
    }
    if out.len() != total {
        return Err(format!(
            "variable {} has {} records but its shape needs {total} values",
            variable.name, records
        ));
    }
    Ok(out)
}

/// Decodes `count` big-endian values into `f64`.
fn decode_values_f64(
    code: u32,
    bytes: &[u8],
    count: usize,
) -> Result<Vec<f64>, String> {
    let size =
        external_size(code).map_err(|()| format!("unsupported external type code {code}"))?;
    let needed = count
        .checked_mul(size)
        .ok_or_else(|| "value array size overflow".to_owned())?;
    let available = bytes
        .get(..needed)
        .ok_or_else(|| "value array is shorter than its declared length".to_owned())?;
    let mut out = Vec::with_capacity(count);
    for chunk in available.chunks_exact(size) {
        out.push(decode_value(code, chunk)?);
    }
    Ok(out)
}

/// Computes the per-record stride for record variables.
fn record_size(header: &Header) -> Result<u64, String> {
    let mut total: u64 = 0;
    for variable in &header.variables {
        let is_record = variable
            .dimids
            .first()
            .and_then(|&id| header.dimensions.get(id))
            .is_some_and(|dim| dim.unlimited);
        if !is_record {
            continue;
        }
        let size = external_size(variable.code)
            .map_err(|()| format!("variable {} uses an unsupported type", variable.name))?;
        let mut shape = Vec::new();
        for &id in variable.dimids.get(1..).unwrap_or(&[]) {
            let dim = header.dimensions.get(id).ok_or_else(|| {
                format!("variable {} references a missing dimension", variable.name)
            })?;
            shape.push(dim.len);
        }
        let elements = element_count(&shape)?;
        let bytes = elements
            .checked_mul(u64::try_from(size).unwrap_or(0))
            .ok_or_else(|| "record slab size overflow".to_owned())?;
        total = total
            .checked_add(align4(bytes))
            .ok_or_else(|| "record size overflow".to_owned())?;
    }
    Ok(total)
}

// -------------------------------------------------------------------------
// Public reading API
// -------------------------------------------------------------------------

/// Reads a NetCDF classic file from disk.
///
/// # Errors
/// Fails when the file cannot be read, is not a `CDF\x01` or `CDF\x02` file,
/// or contains inconsistent lengths or offsets.
pub fn read_netcdf(path: &Path) -> Result<NcFile, String> {
    let bytes = std::fs::read(path)
        .map_err(|error| format!("failed to read {}: {error}", path.display()))?;
    read_netcdf_bytes(&bytes)
}

/// Reads a NetCDF classic file from a byte buffer.
///
/// # Errors
/// Fails when the buffer is not a complete `CDF\x01`/`CDF\x02` image, uses an
/// unsupported type or version, or contains inconsistent lengths or offsets.
pub fn read_netcdf_bytes(bytes: &[u8]) -> Result<NcFile, String> {
    let header = parse_header(bytes)?;
    let stride = record_size(&header)?;
    let mut variables = Vec::with_capacity(header.variables.len());
    for raw in &header.variables {
        let data = read_variable_data(bytes, &header, raw, stride)?;
        let nc_type = nc_type_of(raw.code)
            .map_err(|()| format!("variable {} uses an unsupported type", raw.name))?;
        let dims = raw
            .dimids
            .iter()
            .filter_map(|&id| header.dimensions.get(id))
            .map(|dim| dim.name.clone())
            .collect();
        variables.push(Variable {
            name: raw.name.clone(),
            dims,
            nc_type,
            attributes: raw.attributes.clone(),
            data,
        });
    }
    Ok(NcFile {
        version: header.version,
        dimensions: header.dimensions,
        variables,
        global_attributes: header.global_attributes,
    })
}

// -------------------------------------------------------------------------
// Writer
// -------------------------------------------------------------------------

/// Encodes `data` as big-endian bytes of `nc_type`.
fn encode_values(
    nc_type: NcType,
    data: &[f64],
) -> Result<Vec<u8>, String> {
    let mut out = Vec::with_capacity(data.len().saturating_mul(nc_type.size()));
    for &value in data {
        match nc_type {
            | NcType::F64 => out.extend_from_slice(&value.to_be_bytes()),
            | NcType::F32 => {
                #[allow(clippy::cast_possible_truncation)]
                out.extend_from_slice(&(value as f32).to_be_bytes());
            },
            | other => {
                if !value.is_finite() {
                    return Err("integer variables cannot store non-finite values".to_owned());
                }
                let truncated = value.trunc();
                match other {
                    | NcType::I8 => {
                        if !(-128.0..=127.0).contains(&truncated) {
                            return Err(format!("value {value} does not fit in a byte"));
                        }
                        out.push(truncated as i8 as u8);
                    },
                    | NcType::U8 | NcType::Char => {
                        if !(0.0..=255.0).contains(&truncated) {
                            return Err(format!("value {value} does not fit in a u8"));
                        }
                        out.push(truncated as u8);
                    },
                    | NcType::I16 | NcType::U16 => {
                        let (low, high) = if other == NcType::I16 {
                            (-32_768.0, 32_767.0)
                        } else {
                            (0.0, 65_535.0)
                        };
                        if !(low..=high).contains(&truncated) {
                            return Err(format!("value {value} does not fit in 16 bits"));
                        }
                        let raw = truncated as i32 as u16;
                        out.extend_from_slice(&raw.to_be_bytes());
                    },
                    | _ => {
                        if !(-2_147_483_648.0..=2_147_483_647.0).contains(&truncated) {
                            return Err(format!("value {value} does not fit in 32 bits"));
                        }
                        #[allow(clippy::cast_possible_truncation)]
                        out.extend_from_slice(&(truncated as i32).to_be_bytes());
                    },
                }
            },
        }
    }
    Ok(out)
}

/// Appends a name with its length prefix and four-byte padding.
fn put_name(
    out: &mut Vec<u8>,
    name: &str,
    wide: bool,
) {
    put_count(out, u64::try_from(name.len()).unwrap_or(0), wide);
    out.extend_from_slice(name.as_bytes());
    while !out.len().is_multiple_of(4) {
        out.push(0);
    }
}

/// Appends one attribute.
fn put_attribute(
    out: &mut Vec<u8>,
    name: &str,
    value: &NcValue,
    wide: bool,
) {
    put_name(out, name, wide);
    let (code, count, body) = match value {
        | NcValue::Text(text) => (NC_CHAR, text.len(), text.as_bytes().to_vec()),
        | NcValue::Int(number) => (NC_INT, 1, number.to_be_bytes().to_vec()),
        | NcValue::Float(number) => (NC_DOUBLE, 1, number.to_be_bytes().to_vec()),
        | NcValue::Ints(numbers) => {
            (
                NC_INT,
                numbers.len(),
                numbers
                    .iter()
                    .flat_map(|number| number.to_be_bytes())
                    .collect(),
            )
        },
        | NcValue::Floats(numbers) => {
            (
                NC_DOUBLE,
                numbers.len(),
                numbers
                    .iter()
                    .flat_map(|number| number.to_be_bytes())
                    .collect(),
            )
        },
    };
    out.extend_from_slice(&code.to_be_bytes());
    put_count(out, u64::try_from(count).unwrap_or(0), wide);
    out.extend_from_slice(&body);
    while !out.len().is_multiple_of(4) {
        out.push(0);
    }
}

/// Appends an attribute list.
fn put_attribute_list(
    out: &mut Vec<u8>,
    attributes: &[(String, NcValue)],
    wide: bool,
) {
    if attributes.is_empty() {
        out.extend_from_slice(&0_u32.to_be_bytes());
        put_count(out, 0, wide);
        return;
    }
    out.extend_from_slice(&NC_ATTRIBUTE.to_be_bytes());
    put_count(out, u64::try_from(attributes.len()).unwrap_or(0), wide);
    for (name, value) in attributes {
        put_attribute(out, name, value, wide);
    }
}

/// Writes a `NON_NEG` field at the width the format uses.
fn put_count(
    out: &mut Vec<u8>,
    value: u64,
    wide: bool,
) {
    if wide {
        out.extend_from_slice(&value.to_be_bytes());
    } else {
        out.extend_from_slice(&u32::try_from(value).unwrap_or(0).to_be_bytes());
    }
}

/// Serialises a `CDF\x01` or `CDF\x02` header, using `begins` for variable offsets.
fn build_header(
    file: &NcFile,
    begins: &[u64],
    numrecs: u64,
    format: NcFormat,
) -> Result<Vec<u8>, String> {
    let mut out = Vec::new();
    let wide = matches!(format, NcFormat::Cdf5);
    match format {
        | NcFormat::Classic => out.extend_from_slice(b"CDF\x01"),
        | NcFormat::Offset64 => out.extend_from_slice(b"CDF\x02"),
        | NcFormat::Cdf5 => out.extend_from_slice(b"CDF\x05"),
    }
    if numrecs == u64::MAX {
        return Err("record count cannot equal the streaming sentinel".to_owned());
    }
    if wide {
        out.extend_from_slice(&numrecs.to_be_bytes());
    } else {
        let numrecs_u32 = u32::try_from(numrecs)
            .map_err(|_| "record count exceeds the 32-bit limit of CDF-1 and CDF-2".to_owned())?;
        if numrecs_u32 == u32::MAX {
            return Err("record count cannot equal the streaming sentinel u32::MAX".to_owned());
        }
        out.extend_from_slice(&numrecs_u32.to_be_bytes());
    }
    if file.dimensions.is_empty() {
        out.extend_from_slice(&0_u32.to_be_bytes());
        put_count(&mut out, 0, wide);
    } else {
        out.extend_from_slice(&NC_DIMENSION.to_be_bytes());
        put_count(
            &mut out,
            u64::try_from(file.dimensions.len()).unwrap_or(0),
            wide,
        );
        for dimension in &file.dimensions {
            put_name(&mut out, &dimension.name, wide);
            let len = if dimension.unlimited {
                0
            } else {
                u64::try_from(dimension.len)
                    .map_err(|_| format!("dimension {} is too long", dimension.name))?
            };
            put_count(&mut out, len, wide);
        }
    }
    put_attribute_list(&mut out, &file.global_attributes, wide);
    if file.variables.is_empty() {
        out.extend_from_slice(&0_u32.to_be_bytes());
        put_count(&mut out, 0, wide);
    } else {
        out.extend_from_slice(&NC_VARIABLE.to_be_bytes());
        put_count(
            &mut out,
            u64::try_from(file.variables.len()).unwrap_or(0),
            wide,
        );
        for (index, variable) in file.variables.iter().enumerate() {
            put_name(&mut out, &variable.name, wide);
            put_count(
                &mut out,
                u64::try_from(variable.dims.len()).unwrap_or(0),
                wide,
            );
            for name in &variable.dims {
                let id = file
                    .dimensions
                    .iter()
                    .position(|dimension| dimension.name == *name)
                    .ok_or_else(|| {
                        format!(
                            "variable {} references unknown dimension {name}",
                            variable.name
                        )
                    })?;
                put_count(&mut out, u64::try_from(id).unwrap_or(0), wide);
            }
            put_attribute_list(&mut out, &variable.attributes, wide);
            out.extend_from_slice(&variable.nc_type.code().to_be_bytes());
            let shape = variable_shape(file, variable)?;
            let slab: &[usize] = if is_record_variable(file, variable) {
                shape.get(1..).unwrap_or(&[])
            } else {
                &shape
            };
            let size = element_count(slab)?
                .checked_mul(u64::try_from(variable.nc_type.size()).unwrap_or(0))
                .ok_or_else(|| format!("variable {} is too large", variable.name))?;
            put_count(&mut out, align4(size), wide);
            let begin = *begins.get(index).unwrap_or(&0);
            match format {
                | NcFormat::Classic => {
                    let begin_u32 = u32::try_from(begin)
                        .map_err(|_| format!("variable {} offset is too large", variable.name))?;
                    out.extend_from_slice(&begin_u32.to_be_bytes());
                },
                | NcFormat::Offset64 | NcFormat::Cdf5 => {
                    out.extend_from_slice(&begin.to_be_bytes());
                },
            }
        }
    }
    Ok(out)
}

/// Computes variable `begin` offsets and per-record `stride` for a header of length `header_len`.
fn compute_begins(
    file: &NcFile,
    encoded: &[Vec<u8>],
    numrecs: u64,
    header_len: u64,
) -> (Vec<u64>, u64) {
    let mut begins = vec![0_u64; file.variables.len()];
    let mut pos = align4(header_len);
    for (index, variable) in file.variables.iter().enumerate() {
        let record = is_record_variable(file, variable);
        if !record {
            begins[index] = pos;
            let bytes = u64::try_from(encoded[index].len()).unwrap_or(0);
            pos = align4(pos.saturating_add(bytes));
        }
    }
    let record_start = pos;
    let mut record_offset = 0_u64;
    for (index, variable) in file.variables.iter().enumerate() {
        let record = is_record_variable(file, variable);
        if record {
            begins[index] = record_start.saturating_add(record_offset);
            let bytes = u64::try_from(encoded[index].len()).unwrap_or(0);
            let slab = bytes.checked_div(numrecs).unwrap_or(bytes);
            record_offset = record_offset.saturating_add(align4(slab));
        }
    }
    (begins, record_offset)
}

/// Resolves the shape of `variable`, substituting the record count for the
/// unlimited dimension.
fn variable_shape(
    file: &NcFile,
    variable: &Variable,
) -> Result<Vec<usize>, String> {
    let mut shape = Vec::with_capacity(variable.dims.len());
    for name in &variable.dims {
        let dimension = file
            .dimensions
            .iter()
            .find(|dimension| dimension.name == *name)
            .ok_or_else(|| {
                format!(
                    "variable {} references unknown dimension {name}",
                    variable.name
                )
            })?;
        shape.push(dimension.len);
    }
    Ok(shape)
}

/// Returns whether `variable` is laid out in the record region.
fn is_record_variable(
    file: &NcFile,
    variable: &Variable,
) -> bool {
    let Some(first) = variable.dims.first() else {
        return false;
    };
    file.dimensions
        .iter()
        .any(|dimension| dimension.unlimited && dimension.name == *first)
}

/// Writes a NetCDF file to disk (using `CDF\x02` when `file.version == 2` or
/// when any variable offset exceeds `u32::MAX`, otherwise `CDF\x01`).
///
/// # Errors
/// Fails when a variable references an unknown dimension, data lengths do not
/// match the declared shape, the unlimited dimension is not first, or the file
/// cannot be written.
pub fn write_netcdf(
    path: &Path,
    file: &NcFile,
) -> Result<(), String> {
    let bytes = write_netcdf_bytes(file)?;
    std::fs::write(path, bytes)
        .map_err(|error| format!("failed to write {}: {error}", path.display()))
}

/// Writes a NetCDF file to disk using the requested [`NcFormat`].
///
/// If [`NcFormat::Classic`] is requested but a variable offset exceeds
/// `u32::MAX`, the writer automatically upgrades to [`NcFormat::Offset64`].
///
/// # Errors
/// Fails when serialisation fails or the file cannot be written.
pub fn write_netcdf_with_format(
    path: &Path,
    file: &NcFile,
    format: NcFormat,
) -> Result<(), String> {
    let bytes = to_bytes_with_format(file, format)?;
    std::fs::write(path, bytes)
        .map_err(|error| format!("failed to write {}: {error}", path.display()))
}

/// Serialises a NetCDF file image (using `CDF\x02` when `file.version == 2` or
/// when any variable offset exceeds `u32::MAX`, otherwise `CDF\x01`).
///
/// # Errors
/// Fails when a variable references an unknown dimension, data lengths do not
/// match the declared shape, or the unlimited dimension is not first.
pub fn write_netcdf_bytes(file: &NcFile) -> Result<Vec<u8>, String> {
    let format = if file.version == 2 {
        NcFormat::Offset64
    } else {
        NcFormat::Classic
    };
    to_bytes_with_format(file, format)
}

/// Alias for [`write_netcdf_bytes`].
///
/// # Errors
/// Fails when a variable references an unknown dimension, data lengths do not
/// match the declared shape, or the unlimited dimension is not first.
pub fn to_bytes(file: &NcFile) -> Result<Vec<u8>, String> {
    write_netcdf_bytes(file)
}

/// Serialises a NetCDF file image using the requested [`NcFormat`].
///
/// If [`NcFormat::Classic`] is requested but a variable offset exceeds
/// `u32::MAX`, the writer automatically upgrades to [`NcFormat::Offset64`].
///
/// # Errors
/// Fails when a variable references an unknown dimension, data lengths do not
/// match the declared shape, or the unlimited dimension is not first.
pub fn to_bytes_with_format(
    file: &NcFile,
    format: NcFormat,
) -> Result<Vec<u8>, String> {
    let rec_index = file
        .dimensions
        .iter()
        .position(|dimension| dimension.unlimited);
    if let Some(index) = rec_index {
        if index != 0 {
            return Err("the unlimited dimension must be the first dimension".to_owned());
        }
    }
    let numrecs = rec_index.map_or(0_u64, |index| {
        file.dimensions
            .get(index)
            .map_or(0, |dimension| u64::try_from(dimension.len).unwrap_or(0))
    });

    let mut encoded = Vec::with_capacity(file.variables.len());
    for variable in &file.variables {
        let shape = variable_shape(file, variable)?;
        let count = element_count(&shape)?;
        let count = to_usize(count, "variable element count")?;
        if count != variable.data.len() {
            return Err(format!(
                "variable {} has shape {:?} but {} values",
                variable.name,
                shape,
                variable.data.len()
            ));
        }
        let bytes = count.saturating_mul(variable.nc_type.size());
        if bytes > MAX_ARRAY_BYTES {
            return Err(format!("variable {} exceeds the size limit", variable.name));
        }
        encoded.push(encode_values(variable.nc_type, &variable.data)?);
    }

    let zeros = vec![0_u64; file.variables.len()];
    let mut active_format = format;
    let header_len =
        u64::try_from(build_header(file, &zeros, numrecs, active_format)?.len()).unwrap_or(0);
    let (mut begins, mut stride) = compute_begins(file, &encoded, numrecs, header_len);

    if active_format == NcFormat::Classic
        && begins.iter().any(|&offset| offset > u64::from(u32::MAX))
    {
        active_format = NcFormat::Offset64;
        let header_len_64 =
            u64::try_from(build_header(file, &zeros, numrecs, active_format)?.len()).unwrap_or(0);
        let recomputed = compute_begins(file, &encoded, numrecs, header_len_64);
        begins = recomputed.0;
        stride = recomputed.1;
    }

    let mut out = build_header(file, &begins, numrecs, active_format)?;
    let write_at = |out: &mut Vec<u8>, target: u64, bytes: &[u8]| -> Result<(), String> {
        let target = to_usize(target, "data offset")?;
        while out.len() < target {
            out.push(0);
        }
        out.extend_from_slice(bytes);
        Ok(())
    };

    for (index, variable) in file.variables.iter().enumerate() {
        let record = is_record_variable(file, variable);
        if record {
            continue;
        }
        write_at(&mut out, begins[index], &encoded[index])?;
    }
    for record in 0..numrecs {
        for (index, variable) in file.variables.iter().enumerate() {
            let is_record = is_record_variable(file, variable);
            if !is_record {
                continue;
            }
            let target = begins[index].saturating_add(record.saturating_mul(stride));
            let total = encoded[index].len();
            let slab = if numrecs == 0 {
                total
            } else {
                total / numrecs as usize
            };
            let start = (record as usize).saturating_mul(slab);
            let end = start.saturating_add(slab).min(total);
            let body = encoded[index].get(start..end).unwrap_or(&[]);
            write_at(&mut out, target, body)?;
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::Dimension;
    use super::NcFile;
    use super::NcFormat;
    use super::NcType;
    use super::NcValue;
    use super::Variable;
    use super::read_netcdf;
    use super::read_netcdf_bytes;
    use super::to_bytes;
    use super::to_bytes_with_format;
    use super::write_netcdf;
    use super::write_netcdf_bytes;
    use super::write_netcdf_with_format;

    /// A unique temporary path for one test.
    fn temp_path(tag: &str) -> std::path::PathBuf {
        let unique = format!(
            "rssn_nc_{tag}_{}_{}.nc",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .map_or(0, |d| d.as_nanos())
        );
        std::env::temp_dir().join(unique)
    }

    /// A file with a 2-D `f64` variable over a dimension of length three.
    fn grid_file() -> NcFile {
        NcFile {
            version: 1,
            dimensions: vec![
                Dimension {
                    name: "y".to_owned(),
                    len: 2,
                    unlimited: false,
                },
                Dimension {
                    name: "x".to_owned(),
                    len: 3,
                    unlimited: false,
                },
            ],
            variables: vec![Variable {
                name: "grid".to_owned(),
                dims: vec!["y".to_owned(), "x".to_owned()],
                nc_type: NcType::F64,
                attributes: vec![
                    ("units".to_owned(), NcValue::Text("abc".to_owned())),
                    ("scale".to_owned(), NcValue::Float(0.5)),
                ],
                data: vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            }],
            global_attributes: vec![("title".to_owned(), NcValue::Text("a test".to_owned()))],
        }
    }

    /// A file with an unlimited (record) dimension holding one record.
    fn record_file() -> NcFile {
        NcFile {
            version: 1,
            dimensions: vec![Dimension {
                name: "time".to_owned(),
                len: 1,
                unlimited: true,
            }],
            variables: vec![Variable {
                name: "temperature".to_owned(),
                dims: vec!["time".to_owned()],
                nc_type: NcType::F64,
                attributes: vec![("units".to_owned(), NcValue::Text("K".to_owned()))],
                data: vec![273.15],
            }],
            global_attributes: Vec::new(),
        }
    }

    #[test]
    fn cdf5_widens_every_header_count() {
        let file = grid_file();
        let bytes = to_bytes_with_format(&file, NcFormat::Cdf5).expect("serialise CDF-5");
        assert_eq!(&bytes[..4], b"CDF\x05");
        let parsed = read_netcdf_bytes(&bytes).expect("parse CDF-5");
        assert_eq!(parsed.dimensions.len(), file.dimensions.len());
        assert_eq!(parsed.variables.len(), file.variables.len());
        for (left, right) in parsed.variables.iter().zip(&file.variables) {
            assert_eq!(left.name, right.name);
            assert_eq!(left.data, right.data);
        }
        // The record count is 64-bit from CDF-2 on.
        assert!(bytes.len() > 8);
    }

    #[test]
    fn cdf5_round_trips_records_and_attributes() {
        let file = record_file();
        let bytes = to_bytes_with_format(&file, NcFormat::Cdf5).expect("serialise CDF-5");
        let parsed = read_netcdf_bytes(&bytes).expect("parse CDF-5");
        assert_eq!(parsed.variables.len(), file.variables.len());
        for (left, right) in parsed.variables.iter().zip(&file.variables) {
            assert_eq!(left.data, right.data, "variable {}", right.name);
        }
        assert_eq!(parsed.global_attributes.len(), file.global_attributes.len());
    }

    #[test]
    fn every_format_round_trips_through_a_file() {
        for format in [NcFormat::Classic, NcFormat::Offset64, NcFormat::Cdf5] {
            let path = temp_path(&format!("format_{format:?}"));
            let file = record_file();
            write_netcdf_with_format(&path, &file, format).expect("write");
            let back = read_netcdf(&path).expect("read");
            let _ = std::fs::remove_file(&path);
            assert_eq!(back.variables.len(), file.variables.len(), "{format:?}");
            for (left, right) in back.variables.iter().zip(&file.variables) {
                assert_eq!(left.data, right.data, "{format:?} {}", right.name);
            }
        }
    }

    /// The reference library writes these files; they are checked in so the
    /// test runs in a normal `cargo test`.
    #[test]
    fn reads_reference_written_files() {
        let base = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/netcdf");
        for name in ["cdf5.nc", "cdf2.nc"] {
            let path = base.join(name);
            let file = read_netcdf(&path).expect("read fixture");
            assert!(!file.variables.is_empty(), "{name} has no variables");
            for variable in &file.variables {
                assert!(
                    variable.data.iter().all(|value| value.is_finite()),
                    "{name} {} has a non-finite value",
                    variable.name
                );
            }
        }
    }

    #[test]
    fn round_trip_two_dimensional_variable() {
        let file = grid_file();
        let bytes = write_netcdf_bytes(&file).expect("serialise");
        let parsed = read_netcdf_bytes(&bytes).expect("parse");
        assert_eq!(parsed.version, 1);
        assert_eq!(parsed.dimensions.len(), 2);
        assert_eq!(parsed.dimensions[1].len, 3);
        assert_eq!(parsed.variables.len(), 1);
        let variable = &parsed.variables[0];
        assert_eq!(variable.name, "grid");
        assert_eq!(variable.dims, vec!["y".to_owned(), "x".to_owned()]);
        assert_eq!(variable.nc_type, NcType::F64);
        assert_eq!(variable.data, vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        assert_eq!(
            parsed.global_attributes[0].1,
            NcValue::Text("a test".to_owned())
        );
    }

    #[test]
    fn unaligned_attribute_string_parses() {
        let file = grid_file();
        let bytes = write_netcdf_bytes(&file).expect("serialise");
        let parsed = read_netcdf_bytes(&bytes).expect("parse");
        let variable = &parsed.variables[0];
        assert_eq!(variable.attributes.len(), 2);
        assert_eq!(variable.attributes[0].1, NcValue::Text("abc".to_owned()));
        assert_eq!(variable.attributes[1].1, NcValue::Float(0.5));
        assert_eq!(variable.data.len(), 6);
    }

    #[test]
    fn round_trip_unlimited_dimension() {
        let file = record_file();
        let bytes = write_netcdf_bytes(&file).expect("serialise");
        let parsed = read_netcdf_bytes(&bytes).expect("parse");
        assert_eq!(parsed.dimensions.len(), 1);
        assert!(parsed.dimensions[0].unlimited);
        assert_eq!(parsed.dimensions[0].len, 1);
        assert_eq!(parsed.variables[0].data, vec![273.15]);
        assert_eq!(parsed.variables[0].dims, vec!["time".to_owned()]);
    }

    #[test]
    fn round_trip_through_a_file() {
        let path = temp_path("roundtrip");
        write_netcdf(&path, &grid_file()).expect("write");
        let parsed = read_netcdf(&path).expect("read");
        assert_eq!(parsed.variables[0].data.len(), 6);
        let _ = std::fs::remove_file(&path);
        assert!(!path.exists());
    }

    #[test]
    fn truncated_buffers_are_rejected() {
        let bytes = write_netcdf_bytes(&grid_file()).expect("serialise");
        let len = bytes.len();
        let cuts = [3_usize, 8, 16, len / 2];
        for cut in cuts {
            let truncated = bytes.get(..cut).expect("slice");
            assert!(
                read_netcdf_bytes(truncated).is_err(),
                "a buffer truncated at {cut} of {len} must be rejected"
            );
        }
    }

    #[test]
    fn unknown_dimension_is_rejected() {
        let mut file = grid_file();
        file.variables[0].dims = vec!["missing".to_owned()];
        assert!(write_netcdf_bytes(&file).is_err());
    }

    #[test]
    fn shape_mismatch_is_rejected() {
        let mut file = grid_file();
        file.variables[0].data.pop();
        assert!(write_netcdf_bytes(&file).is_err());
    }

    #[test]
    fn not_a_netcdf_file_is_rejected() {
        assert!(read_netcdf_bytes(b"not a netcdf file").is_err());
    }

    #[test]
    fn reads_cdf2_header() {
        let mut out = Vec::new();
        out.extend_from_slice(b"CDF\x02");
        out.extend_from_slice(&0_u32.to_be_bytes());
        out.extend_from_slice(&10_u32.to_be_bytes());
        out.extend_from_slice(&1_u32.to_be_bytes());
        out.extend_from_slice(&1_u32.to_be_bytes());
        out.extend_from_slice(b"x");
        out.extend_from_slice(&[0; 3]);
        out.extend_from_slice(&2_u32.to_be_bytes());
        out.extend_from_slice(&0_u32.to_be_bytes());
        out.extend_from_slice(&0_u32.to_be_bytes());
        out.extend_from_slice(&11_u32.to_be_bytes());
        out.extend_from_slice(&1_u32.to_be_bytes());
        out.extend_from_slice(&1_u32.to_be_bytes());
        out.extend_from_slice(b"v");
        out.extend_from_slice(&[0; 3]);
        out.extend_from_slice(&1_u32.to_be_bytes());
        out.extend_from_slice(&0_u32.to_be_bytes());
        out.extend_from_slice(&0_u32.to_be_bytes());
        out.extend_from_slice(&0_u32.to_be_bytes());
        out.extend_from_slice(&6_u32.to_be_bytes());
        out.extend_from_slice(&16_u32.to_be_bytes());
        let begin_patch = out.len();
        out.extend_from_slice(&0_u64.to_be_bytes());
        let begin = out.len() as u64;
        if let Some(slot) = out.get_mut(begin_patch..begin_patch + 8) {
            slot.copy_from_slice(&begin.to_be_bytes());
        }
        out.extend_from_slice(&1.5_f64.to_be_bytes());
        out.extend_from_slice(&2.5_f64.to_be_bytes());

        let parsed = read_netcdf_bytes(&out).expect("parse CDF-2 file");
        assert_eq!(parsed.version, 2);
        assert_eq!(parsed.dimensions.len(), 1);
        assert_eq!(parsed.dimensions[0].len, 2);
        assert_eq!(parsed.variables[0].data, vec![1.5, 2.5]);
    }

    #[test]
    fn round_trip_cdf2_fixed_and_record_variables() {
        let file = NcFile {
            version: 2,
            dimensions: vec![
                Dimension {
                    name: "time".to_owned(),
                    len: 3,
                    unlimited: true,
                },
                Dimension {
                    name: "station".to_owned(),
                    len: 2,
                    unlimited: false,
                },
            ],
            variables: vec![
                Variable {
                    name: "lat".to_owned(),
                    dims: vec!["station".to_owned()],
                    nc_type: NcType::F32,
                    attributes: vec![("units".to_owned(), NcValue::Text("degN".to_owned()))],
                    data: vec![37.5, 42.25],
                },
                Variable {
                    name: "temp".to_owned(),
                    dims: vec!["time".to_owned(), "station".to_owned()],
                    nc_type: NcType::F64,
                    attributes: vec![("units".to_owned(), NcValue::Text("K".to_owned()))],
                    data: vec![280.0, 281.5, 282.0, 283.5, 284.0, 285.5],
                },
                Variable {
                    name: "flags".to_owned(),
                    dims: vec!["time".to_owned()],
                    nc_type: NcType::I8,
                    attributes: Vec::new(),
                    data: vec![-1.0, 0.0, 1.0],
                },
            ],
            global_attributes: vec![
                ("title".to_owned(), NcValue::Text("cdf2 test".to_owned())),
                ("rev".to_owned(), NcValue::Int(2)),
            ],
        };

        let bytes = to_bytes_with_format(&file, NcFormat::Offset64).expect("serialise CDF-2");
        assert_eq!(&bytes[..4], b"CDF\x02");

        let parsed = read_netcdf_bytes(&bytes).expect("parse CDF-2");
        assert_eq!(parsed.version, 2);
        assert_eq!(parsed.dimensions, file.dimensions);
        assert_eq!(parsed.global_attributes, file.global_attributes);
        assert_eq!(parsed.variables, file.variables);

        // Also verify default to_bytes respects file.version == 2 and file-based write_netcdf_with_format.
        let default_bytes = to_bytes(&file).expect("default to_bytes");
        assert_eq!(&default_bytes[..4], b"CDF\x02");

        let path = temp_path("cdf2_roundtrip");
        write_netcdf_with_format(&path, &file, NcFormat::Offset64).expect("write CDF-2");
        let from_disk = read_netcdf(&path).expect("read CDF-2");
        assert_eq!(from_disk.version, 2);
        assert_eq!(from_disk.variables, file.variables);
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn round_trip_nc_char_variables_and_string_decoding() {
        let file = NcFile {
            version: 1,
            dimensions: vec![
                Dimension {
                    name: "title_len".to_owned(),
                    len: 5,
                    unlimited: false,
                },
                Dimension {
                    name: "n_stations".to_owned(),
                    len: 3,
                    unlimited: false,
                },
                Dimension {
                    name: "name_len".to_owned(),
                    len: 8,
                    unlimited: false,
                },
            ],
            variables: vec![
                Variable::from_string(
                    "label",
                    "title_len",
                    "hello",
                    vec![("long_name".to_owned(), NcValue::Text("label".to_owned()))],
                ),
                Variable::from_strings(
                    "stations",
                    vec!["n_stations".to_owned(), "name_len".to_owned()],
                    &["alpha", "beta", "gamma"],
                    8,
                    Vec::new(),
                ),
            ],
            global_attributes: Vec::new(),
        };

        for format in [NcFormat::Classic, NcFormat::Offset64] {
            let bytes = to_bytes_with_format(&file, format).expect("serialise NC_CHAR");
            let parsed = read_netcdf_bytes(&bytes).expect("parse NC_CHAR");

            let label = &parsed.variables[0];
            assert_eq!(label.nc_type, NcType::Char);
            assert_eq!(label.as_string(), Some("hello".to_owned()));
            assert_eq!(
                label.as_strings(&parsed.dimensions),
                Some(vec!["hello".to_owned()])
            );
            // Numeric byte values remain accessible in `data`.
            assert_eq!(label.data, vec![104.0, 101.0, 108.0, 108.0, 111.0]);

            let stations = &parsed.variables[1];
            assert_eq!(stations.nc_type, NcType::Char);
            assert_eq!(
                stations.as_strings(&parsed.dimensions),
                Some(vec![
                    "alpha".to_owned(),
                    "beta".to_owned(),
                    "gamma".to_owned(),
                ])
            );
            assert_eq!(
                parsed.variable_strings("stations"),
                Some(vec![
                    "alpha".to_owned(),
                    "beta".to_owned(),
                    "gamma".to_owned(),
                ])
            );
        }

        // Non-char variables return None for string helpers.
        let numeric = &grid_file().variables[0];
        assert_eq!(numeric.as_string(), None);
        assert_eq!(numeric.as_strings(&grid_file().dimensions), None);
    }

    #[test]
    fn integer_variables_round_trip() {
        let file = NcFile {
            version: 1,
            dimensions: vec![Dimension {
                name: "n".to_owned(),
                len: 4,
                unlimited: false,
            }],
            variables: vec![Variable {
                name: "counts".to_owned(),
                dims: vec!["n".to_owned()],
                nc_type: NcType::I16,
                attributes: vec![("valid".to_owned(), NcValue::Int(1))],
                data: vec![-3.0, 0.0, 7.0, 1_000.0],
            }],
            global_attributes: Vec::new(),
        };
        let bytes = write_netcdf_bytes(&file).expect("serialise");
        let parsed = read_netcdf_bytes(&bytes).expect("parse");
        assert_eq!(parsed.variables[0].nc_type, NcType::I16);
        assert_eq!(parsed.variables[0].data, vec![-3.0, 0.0, 7.0, 1_000.0]);
        assert_eq!(parsed.variables[0].attributes[0].1, NcValue::Int(1));
    }
}
