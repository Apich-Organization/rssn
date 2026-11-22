//! Apache Parquet files with a flat schema.
//!
//! The writer emits one row group with one `PLAIN`-encoded, uncompressed data
//! page per column. Optional columns carry RLE/bit-packed definition levels.
//! The reader understands that layout plus Snappy-compressed pages, and
//! validates every offset and length against the buffer before slicing.
//!
//! # Limitations
//!
//! - Flat schemas only: repeated fields and nested groups are rejected.
//! - One row group and one data page per column are written; the reader can
//!   concatenate several row groups but expects a single data page per chunk.
//! - Only `PLAIN` value encoding is read or written. Dictionary encoding,
//!   delta encoding and byte-stream split are not supported.
//! - Pages are written uncompressed. Reading supports `UNCOMPRESSED` and
//!   `SNAPPY`; other codecs return an error.
//! - Strings are read as UTF-8 and invalid sequences are rejected, so
//!   arbitrary `BYTE_ARRAY` payloads are not round-tripped byte for byte.
//! - `INT96`, `FIXED_LEN_BYTE_ARRAY` and logical types such as decimals are
//!   rejected or flattened to their physical type.

mod plain;
mod snappy;
mod thrift;

use std::path::Path;

use crate::io::frame::Column;

use thrift::ColumnChunk;
use thrift::ColumnMetaData;
use thrift::DataPageHeader;
use thrift::FileMetaData;
use thrift::PageHeader;
use thrift::RowGroup;
use thrift::SchemaElement;
use thrift::Statistics;

/// The magic bytes at both ends of a Parquet file.
const MAGIC: &[u8; 4] = b"PAR1";
/// Upper bound on declared rows, guarding against hostile metadata.
const MAX_ROWS: usize = 100_000_000;
/// Page type code for a data page.
const PAGE_DATA: i32 = 0;
/// Encoding code for `PLAIN`.
const ENCODING_PLAIN: i32 = 0;
/// Encoding code for RLE.
const ENCODING_RLE: i32 = 3;
/// Codec code for no compression.
const CODEC_UNCOMPRESSED: i32 = 0;
/// Codec code for Snappy.
const CODEC_SNAPPY: i32 = 1;
/// Repetition code for a required field.
const REPETITION_REQUIRED: i32 = 0;
/// Repetition code for an optional field.
const REPETITION_OPTIONAL: i32 = 1;
/// Repetition code for a repeated field.
const REPETITION_REPEATED: i32 = 2;
/// Converted type code for UTF-8 strings.
const CONVERTED_UTF8: i32 = 0;

/// The physical storage type of a Parquet column.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Physical {
    /// One bit per value.
    Boolean,
    /// Four-byte signed integers.
    Int32,
    /// Eight-byte signed integers.
    Int64,
    /// Four-byte floats.
    Float,
    /// Eight-byte floats.
    Double,
    /// Length-prefixed byte arrays, read back as UTF-8 strings.
    ByteArray,
}

/// One leaf column of a flat schema.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Field {
    /// Column name.
    pub name: String,
    /// Physical storage type.
    pub physical: Physical,
    /// Whether the column may contain missing values.
    pub optional: bool,
}

/// A flat Parquet schema.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Schema {
    /// Leaf fields, in file order.
    pub fields: Vec<Field>,
}

/// A decoded column, before the null mask is attached.
enum Decoded {
    /// Boolean values.
    Bool(Vec<bool>),
    /// 32-bit integers.
    I32(Vec<i32>),
    /// 64-bit integers.
    I64(Vec<i64>),
    /// 32-bit floats.
    F32(Vec<f32>),
    /// 64-bit floats.
    F64(Vec<f64>),
    /// UTF-8 strings.
    Bytes(Vec<String>),
}

/// Writes a flat frame to a Parquet file.
///
/// # Errors
///
/// Returns an error when the schema does not match the columns, a value does
/// not fit its physical type, or the file cannot be written.
pub fn write_parquet(
    path: &Path,
    schema: &Schema,
    columns: &[Column],
) -> Result<(), String> {
    let bytes = write_parquet_bytes(schema, columns)?;
    std::fs::write(path, bytes).map_err(|error| error.to_string())
}

/// Writes a flat frame to an in-memory Parquet file.
///
/// # Errors
///
/// Returns an error when the schema does not match the columns or a value
/// does not fit its physical type.
pub fn write_parquet_bytes(
    schema: &Schema,
    columns: &[Column],
) -> Result<Vec<u8>, String> {
    validate(schema, columns)?;
    let rows = columns.first().map_or(0, Column::len);
    if rows > MAX_ROWS {
        return Err("frame has too many rows for a parquet file".to_string());
    }
    for (field, column) in schema.fields.iter().zip(columns) {
        if column.len() != rows {
            return Err(format!(
                "column `{}` has {} values but the frame has {rows} rows",
                field.name,
                column.len()
            ));
        }
    }
    let mut file = Vec::new();
    file.extend_from_slice(MAGIC);
    let mut chunks = Vec::with_capacity(schema.fields.len());
    let mut total_size = 0i64;
    for (field, column) in schema.fields.iter().zip(columns) {
        let offset = i64::try_from(file.len())
            .map_err(|_| "parquet file has grown too large".to_string())?;
        let metadata = write_column_chunk(field, column, offset, &mut file)?;
        total_size = total_size
            .checked_add(metadata.total_uncompressed_size)
            .ok_or_else(|| "parquet file size overflows".to_string())?;
        chunks.push(ColumnChunk {
            file_offset: offset,
            meta_data: Some(metadata),
        });
    }
    let row_group = RowGroup {
        columns: chunks,
        total_byte_size: total_size,
        num_rows: i64::try_from(rows).map_err(|_| "row count overflows".to_string())?,
        file_offset: Some(4),
        total_compressed_size: Some(total_size),
    };
    let mut elements = Vec::with_capacity(schema.fields.len().saturating_add(1));
    elements.push(SchemaElement {
        type_: None,
        type_length: None,
        repetition_type: Some(REPETITION_REQUIRED),
        name: "schema".to_string(),
        num_children: Some(
            i32::try_from(schema.fields.len()).map_err(|_| "too many fields".to_string())?,
        ),
        converted_type: None,
        scale: None,
        precision: None,
        field_id: None,
    });
    for field in &schema.fields {
        elements.push(SchemaElement {
            type_: Some(physical_code(field.physical)),
            type_length: None,
            repetition_type: Some(if field.optional {
                REPETITION_OPTIONAL
            } else {
                REPETITION_REQUIRED
            }),
            name: field.name.clone(),
            num_children: None,
            converted_type: if field.physical == Physical::ByteArray {
                Some(CONVERTED_UTF8)
            } else {
                None
            },
            scale: None,
            precision: None,
            field_id: None,
        });
    }
    let metadata = FileMetaData {
        version: 1,
        schema: elements,
        num_rows: i64::try_from(rows).map_err(|_| "row count overflows".to_string())?,
        row_groups: vec![row_group],
    };
    let metadata_bytes = thrift::encode_file_metadata(&metadata);
    let metadata_len = u32::try_from(metadata_bytes.len())
        .map_err(|_| "parquet metadata is too large".to_string())?;
    file.extend_from_slice(&metadata_bytes);
    file.extend_from_slice(&metadata_len.to_le_bytes());
    file.extend_from_slice(MAGIC);
    Ok(file)
}

/// Reads a flat Parquet file.
///
/// # Errors
///
/// Returns an error when the file cannot be read or is malformed.
pub fn read_parquet(path: &Path) -> Result<(Schema, Vec<Column>), String> {
    let bytes = std::fs::read(path).map_err(|error| error.to_string())?;
    read_parquet_bytes(&bytes)
}

/// Reads a flat Parquet file from memory.
///
/// # Errors
///
/// Returns an error when the buffer is not a valid, flat, `PLAIN`-encoded
/// Parquet file.
pub fn read_parquet_bytes(bytes: &[u8]) -> Result<(Schema, Vec<Column>), String> {
    let len = bytes.len();
    if len < 12 {
        return Err("parquet file is too small".to_string());
    }
    if bytes.get(0..4) != Some(&MAGIC[..]) {
        return Err("parquet file does not start with the magic bytes".to_string());
    }
    if bytes.get(len - 4..len) != Some(&MAGIC[..]) {
        return Err("parquet file does not end with the magic bytes".to_string());
    }
    let length_bytes = bytes
        .get(len - 8..len - 4)
        .ok_or_else(|| "parquet footer is truncated".to_string())?;
    let mut raw = [0u8; 4];
    raw.copy_from_slice(length_bytes);
    let metadata_len = u32::from_le_bytes(raw) as usize;
    let metadata_start = len
        .checked_sub(8)
        .and_then(|value| value.checked_sub(metadata_len))
        .ok_or_else(|| "parquet metadata length exceeds the file".to_string())?;
    if metadata_start < 4 {
        return Err("parquet metadata overlaps the file header".to_string());
    }
    let metadata_bytes = bytes
        .get(metadata_start..len - 8)
        .ok_or_else(|| "parquet metadata is truncated".to_string())?;
    let metadata = thrift::decode_file_metadata(metadata_bytes)?;
    if metadata.version < 1 {
        return Err(format!("unsupported parquet version {}", metadata.version));
    }
    if metadata.num_rows < 0 || metadata.num_rows > MAX_ROWS as i64 {
        return Err("parquet file declares an impossible row count".to_string());
    }
    let fields = parse_schema(&metadata)?;
    let mut chunks: Vec<Vec<Column>> = vec![Vec::new(); fields.len()];
    let mut total_rows = 0usize;
    for row_group in &metadata.row_groups {
        if row_group.num_rows < 0 || row_group.num_rows > MAX_ROWS as i64 {
            return Err("row group declares an impossible row count".to_string());
        }
        let rows = usize::try_from(row_group.num_rows)
            .map_err(|_| "row group row count is out of range".to_string())?;
        if row_group.columns.len() != fields.len() {
            return Err(format!(
                "row group has {} columns but the schema has {}",
                row_group.columns.len(),
                fields.len()
            ));
        }
        for (index, chunk) in row_group.columns.iter().enumerate() {
            let field = fields
                .get(index)
                .ok_or_else(|| "row group column index is out of range".to_string())?;
            let chunk_metadata = chunk
                .meta_data
                .as_ref()
                .ok_or_else(|| format!("column chunk `{}` is missing metadata", field.name))?;
            if chunk_metadata.path_in_schema.len() != 1
                || chunk_metadata.path_in_schema.first().map(String::as_str)
                    != Some(field.name.as_str())
            {
                return Err(format!(
                    "column chunk path does not match field `{}`",
                    field.name
                ));
            }
            let column = read_column_chunk(bytes, chunk_metadata, field, metadata_start)?;
            if column.len() != rows {
                return Err(format!(
                    "column `{}` has {} values but its row group has {rows} rows",
                    field.name,
                    column.len()
                ));
            }
            chunks
                .get_mut(index)
                .ok_or_else(|| "column index is out of range".to_string())?
                .push(column);
        }
        total_rows = total_rows
            .checked_add(rows)
            .ok_or_else(|| "parquet row count overflows".to_string())?;
    }
    if total_rows != metadata.num_rows as usize {
        return Err(format!(
            "row groups hold {total_rows} rows but the footer declares {}",
            metadata.num_rows
        ));
    }
    let mut columns = Vec::with_capacity(fields.len());
    for (index, field) in fields.iter().enumerate() {
        let parts = chunks
            .get_mut(index)
            .ok_or_else(|| "column index is out of range".to_string())?;
        let column = if parts.is_empty() {
            empty_column(field)
        } else {
            merge_chunks(std::mem::take(parts))?
        };
        columns.push(column);
    }
    Ok((Schema { fields }, columns))
}

/// Checks that the schema and columns agree.
fn validate(
    schema: &Schema,
    columns: &[Column],
) -> Result<(), String> {
    if schema.fields.is_empty() {
        return Err("parquet schema must contain at least one field".to_string());
    }
    if schema.fields.len() != columns.len() {
        return Err(format!(
            "schema has {} fields but {} columns were given",
            schema.fields.len(),
            columns.len()
        ));
    }
    let mut names: Vec<&str> = Vec::with_capacity(schema.fields.len());
    for (field, column) in schema.fields.iter().zip(columns) {
        if field.name.is_empty() {
            return Err("parquet field names must not be empty".to_string());
        }
        if names.contains(&field.name.as_str()) {
            return Err(format!("duplicate parquet field name `{}`", field.name));
        }
        names.push(&field.name);
        if !column_matches(column, field.physical) {
            return Err(format!(
                "column `{}` does not match physical type {:?}",
                field.name, field.physical
            ));
        }
        if column.is_nullable() != field.optional {
            return Err(format!(
                "column `{}` nullability does not match the schema",
                field.name
            ));
        }
    }
    Ok(())
}

/// Checks that a column can be stored with a physical type.
fn column_matches(
    column: &Column,
    physical: Physical,
) -> bool {
    match column {
        | Column::Int(_) => matches!(physical, Physical::Int32 | Physical::Int64),
        | Column::Float(_) => matches!(physical, Physical::Float | Physical::Double),
        | Column::Bool(_) => physical == Physical::Boolean,
        | Column::Str(_) => physical == Physical::ByteArray,
        | Column::Nullable(inner, _) => column_matches(inner, physical),
    }
}

/// Maps a physical type to its Parquet type code.
fn physical_code(physical: Physical) -> i32 {
    match physical {
        | Physical::Boolean => 0,
        | Physical::Int32 => 1,
        | Physical::Int64 => 2,
        | Physical::Float => 4,
        | Physical::Double => 5,
        | Physical::ByteArray => 6,
    }
}

/// Maps a Parquet type code to a physical type.
fn physical_from_code(code: i32) -> Result<Physical, String> {
    match code {
        | 0 => Ok(Physical::Boolean),
        | 1 => Ok(Physical::Int32),
        | 2 => Ok(Physical::Int64),
        | 4 => Ok(Physical::Float),
        | 5 => Ok(Physical::Double),
        | 6 => Ok(Physical::ByteArray),
        | other => Err(format!("unsupported parquet physical type {other}")),
    }
}

/// Extracts the flat fields from the decoded footer.
fn parse_schema(metadata: &FileMetaData) -> Result<Vec<Field>, String> {
    let mut elements = metadata.schema.iter();
    let root = elements
        .next()
        .ok_or_else(|| "parquet schema is empty".to_string())?;
    if root.type_.is_some() {
        return Err("parquet schema root must be a group".to_string());
    }
    let mut fields = Vec::new();
    for element in elements {
        let code = element
            .type_
            .ok_or_else(|| format!("nested parquet field `{}` is not supported", element.name))?;
        let physical = physical_from_code(code)?;
        let repetition = element.repetition_type.unwrap_or(REPETITION_REQUIRED);
        if repetition == REPETITION_REPEATED {
            return Err(format!(
                "repeated parquet field `{}` is not supported",
                element.name
            ));
        }
        fields.push(Field {
            name: element.name.clone(),
            physical,
            optional: repetition == REPETITION_OPTIONAL,
        });
    }
    if fields.is_empty() {
        return Err("parquet schema has no columns".to_string());
    }
    if let Some(children) = root.num_children {
        if usize::try_from(children).ok() != Some(fields.len()) {
            return Err("parquet root child count does not match the fields".to_string());
        }
    }
    Ok(fields)
}

/// Returns a typed empty column for a field.
fn empty_column(field: &Field) -> Column {
    let base = match field.physical {
        | Physical::Boolean => Column::Bool(Vec::new()),
        | Physical::Int32 | Physical::Int64 => Column::Int(Vec::new()),
        | Physical::Float | Physical::Double => Column::Float(Vec::new()),
        | Physical::ByteArray => Column::Str(Vec::new()),
    };
    if field.optional {
        Column::Nullable(Box::new(base), Vec::new())
    } else {
        base
    }
}

/// Concatenates the per-row-group chunks of one column.
fn merge_chunks(mut chunks: Vec<Column>) -> Result<Column, String> {
    let first = chunks
        .first()
        .cloned()
        .ok_or_else(|| "column has no chunks".to_string())?;
    let mut merged = first;
    for chunk in chunks.drain(1..) {
        merged = concat_columns(&merged, &chunk)?;
    }
    Ok(merged)
}

/// Concatenates two columns of the same variant.
fn concat_columns(
    left: &Column,
    right: &Column,
) -> Result<Column, String> {
    match (left, right) {
        | (Column::Int(first), Column::Int(second)) => {
            let mut values = first.clone();
            values.extend_from_slice(second);
            Ok(Column::Int(values))
        },
        | (Column::Float(first), Column::Float(second)) => {
            let mut values = first.clone();
            values.extend_from_slice(second);
            Ok(Column::Float(values))
        },
        | (Column::Bool(first), Column::Bool(second)) => {
            let mut values = first.clone();
            values.extend_from_slice(second);
            Ok(Column::Bool(values))
        },
        | (Column::Str(first), Column::Str(second)) => {
            let mut values = first.clone();
            values.extend_from_slice(second);
            Ok(Column::Str(values))
        },
        | (Column::Nullable(first, first_mask), Column::Nullable(second, second_mask)) => {
            let inner = concat_columns(first, second)?;
            let mut mask = first_mask.clone();
            mask.extend_from_slice(second_mask);
            Ok(Column::Nullable(Box::new(inner), mask))
        },
        | _ => Err("parquet chunks have mismatched types".to_string()),
    }
}

/// Reads the non-missing integer at `row`, looking through null masks.
fn int_at(
    column: &Column,
    row: usize,
) -> Option<i64> {
    match column {
        | Column::Int(values) => values.get(row).copied(),
        | Column::Nullable(inner, mask) => {
            if mask.get(row).copied().unwrap_or(true) {
                None
            } else {
                int_at(inner, row)
            }
        },
        | _ => None,
    }
}

/// Reads the non-missing float at `row`, looking through null masks.
fn float_at(
    column: &Column,
    row: usize,
) -> Option<f64> {
    match column {
        | Column::Float(values) => values.get(row).copied(),
        | Column::Int(values) => values.get(row).map(|value| *value as f64),
        | Column::Nullable(inner, mask) => {
            if mask.get(row).copied().unwrap_or(true) {
                None
            } else {
                float_at(inner, row)
            }
        },
        | _ => None,
    }
}

/// Reads the non-missing boolean at `row`, looking through null masks.
fn bool_at(
    column: &Column,
    row: usize,
) -> Option<bool> {
    match column {
        | Column::Bool(values) => values.get(row).copied(),
        | Column::Nullable(inner, mask) => {
            if mask.get(row).copied().unwrap_or(true) {
                None
            } else {
                bool_at(inner, row)
            }
        },
        | _ => None,
    }
}

/// Reads the non-missing string at `row`, looking through null masks.
fn str_at(
    column: &Column,
    row: usize,
) -> Option<&str> {
    match column {
        | Column::Str(values) => values.get(row).map(String::as_str),
        | Column::Nullable(inner, mask) => {
            if mask.get(row).copied().unwrap_or(true) {
                None
            } else {
                str_at(inner, row)
            }
        },
        | _ => None,
    }
}

/// Extracts the physical values of a column, dropping missing slots.
fn extract_values(
    field: &Field,
    column: &Column,
) -> Result<Decoded, String> {
    let rows = column.len();
    let mut decoded = match field.physical {
        | Physical::Boolean => Decoded::Bool(Vec::new()),
        | Physical::Int32 => Decoded::I32(Vec::new()),
        | Physical::Int64 => Decoded::I64(Vec::new()),
        | Physical::Float => Decoded::F32(Vec::new()),
        | Physical::Double => Decoded::F64(Vec::new()),
        | Physical::ByteArray => Decoded::Bytes(Vec::new()),
    };
    for row in 0..rows {
        if column.is_null(row) {
            continue;
        }
        match &mut decoded {
            | Decoded::Bool(values) => {
                values.push(
                    bool_at(column, row)
                        .ok_or_else(|| format!("column `{}` is not boolean", field.name))?,
                )
            },
            | Decoded::I32(values) => {
                let value = int_at(column, row)
                    .ok_or_else(|| format!("column `{}` is not integral", field.name))?;
                values.push(i32::try_from(value).map_err(|_| {
                    format!(
                        "value {value} in column `{}` does not fit in int32",
                        field.name
                    )
                })?);
            },
            | Decoded::I64(values) => {
                values.push(
                    int_at(column, row)
                        .ok_or_else(|| format!("column `{}` is not integral", field.name))?,
                )
            },
            | Decoded::F32(values) => {
                values.push(
                    float_at(column, row)
                        .ok_or_else(|| format!("column `{}` is not numeric", field.name))?
                        as f32,
                )
            },
            | Decoded::F64(values) => {
                values.push(
                    float_at(column, row)
                        .ok_or_else(|| format!("column `{}` is not numeric", field.name))?,
                )
            },
            | Decoded::Bytes(values) => {
                values.push(
                    str_at(column, row)
                        .ok_or_else(|| format!("column `{}` is not a string", field.name))?
                        .to_string(),
                )
            },
        }
    }
    Ok(decoded)
}

/// Encodes already-extracted values with `PLAIN`.
fn encode_values(values: &Decoded) -> Vec<u8> {
    match values {
        | Decoded::Bool(items) => plain::encode_bool(items),
        | Decoded::I32(items) => plain::encode_i32(items),
        | Decoded::I64(items) => plain::encode_i64(items),
        | Decoded::F32(items) => plain::encode_f32(items),
        | Decoded::F64(items) => plain::encode_f64(items),
        | Decoded::Bytes(items) => {
            let bytes: Vec<Vec<u8>> = items
                .iter()
                .map(|value| value.as_bytes().to_vec())
                .collect();
            plain::encode_bytes(&bytes)
        },
    }
}

/// Builds column statistics, including the PLAIN-encoded min and max.
fn build_statistics(
    values: &Decoded,
    null_count: i64,
) -> Statistics {
    let mut statistics = Statistics {
        null_count: Some(null_count),
        ..Statistics::default()
    };
    match values {
        | Decoded::Bool(_) => {},
        | Decoded::I32(items) => {
            if let (Some(min), Some(max)) = (items.iter().min(), items.iter().max()) {
                statistics.min_value = Some(min.to_le_bytes().to_vec());
                statistics.max_value = Some(max.to_le_bytes().to_vec());
            }
        },
        | Decoded::I64(items) => {
            if let (Some(min), Some(max)) = (items.iter().min(), items.iter().max()) {
                statistics.min_value = Some(min.to_le_bytes().to_vec());
                statistics.max_value = Some(max.to_le_bytes().to_vec());
            }
        },
        | Decoded::F32(items) => {
            if let (Some(min), Some(max)) = (minimum_f32(items), maximum_f32(items)) {
                statistics.min_value = Some(min.to_le_bytes().to_vec());
                statistics.max_value = Some(max.to_le_bytes().to_vec());
            }
        },
        | Decoded::F64(items) => {
            if let (Some(min), Some(max)) = (minimum_f64(items), maximum_f64(items)) {
                statistics.min_value = Some(min.to_le_bytes().to_vec());
                statistics.max_value = Some(max.to_le_bytes().to_vec());
            }
        },
        | Decoded::Bytes(items) => {
            if let (Some(min), Some(max)) = (items.iter().min(), items.iter().max()) {
                statistics.min_value = Some(plain::encode_bytes(&[min.as_bytes().to_vec()]));
                statistics.max_value = Some(plain::encode_bytes(&[max.as_bytes().to_vec()]));
            }
        },
    }
    statistics.min = statistics.min_value.clone();
    statistics.max = statistics.max_value.clone();
    statistics
}

/// Smallest `f32`, using a total order so NaNs sort consistently.
fn minimum_f32(items: &[f32]) -> Option<f32> {
    let mut best: Option<f32> = None;
    for value in items {
        best = Some(best.map_or(*value, |current| {
            if value.total_cmp(&current) == std::cmp::Ordering::Less {
                *value
            } else {
                current
            }
        }));
    }
    best
}

/// Largest `f32`, using a total order so NaNs sort consistently.
fn maximum_f32(items: &[f32]) -> Option<f32> {
    let mut best: Option<f32> = None;
    for value in items {
        best = Some(best.map_or(*value, |current| {
            if value.total_cmp(&current) == std::cmp::Ordering::Greater {
                *value
            } else {
                current
            }
        }));
    }
    best
}

/// Smallest `f64`, using a total order so NaNs sort consistently.
fn minimum_f64(items: &[f64]) -> Option<f64> {
    let mut best: Option<f64> = None;
    for value in items {
        best = Some(best.map_or(*value, |current| {
            if value.total_cmp(&current) == std::cmp::Ordering::Less {
                *value
            } else {
                current
            }
        }));
    }
    best
}

/// Largest `f64`, using a total order so NaNs sort consistently.
fn maximum_f64(items: &[f64]) -> Option<f64> {
    let mut best: Option<f64> = None;
    for value in items {
        best = Some(best.map_or(*value, |current| {
            if value.total_cmp(&current) == std::cmp::Ordering::Greater {
                *value
            } else {
                current
            }
        }));
    }
    best
}

/// Writes one column chunk and returns its metadata.
fn write_column_chunk(
    field: &Field,
    column: &Column,
    offset: i64,
    file: &mut Vec<u8>,
) -> Result<ColumnMetaData, String> {
    let rows = column.len();
    let values = extract_values(field, column)?;
    let null_count = (0..rows).filter(|&row| column.is_null(row)).count();
    let mut payload = Vec::new();
    if field.optional {
        let levels: Vec<bool> = (0..rows).map(|row| column.is_null(row)).collect();
        let encoded = plain::encode_levels(&levels);
        let level_len = u32::try_from(encoded.len())
            .map_err(|_| "definition levels are too large".to_string())?;
        payload.extend_from_slice(&level_len.to_le_bytes());
        payload.extend_from_slice(&encoded);
    }
    payload.extend_from_slice(&encode_values(&values));
    let page_size = i32::try_from(payload.len())
        .map_err(|_| format!("column `{}` page is too large", field.name))?;
    let header = PageHeader {
        type_: PAGE_DATA,
        uncompressed_page_size: page_size,
        compressed_page_size: page_size,
        data_page: Some(DataPageHeader {
            num_values: i32::try_from(rows)
                .map_err(|_| format!("column `{}` has too many rows", field.name))?,
            encoding: ENCODING_PLAIN,
            definition_level_encoding: ENCODING_RLE,
            repetition_level_encoding: ENCODING_RLE,
        }),
    };
    let header_bytes = thrift::encode_page_header(&header);
    let total = header_bytes
        .len()
        .checked_add(payload.len())
        .ok_or_else(|| "column chunk size overflows".to_string())?;
    file.extend_from_slice(&header_bytes);
    file.extend_from_slice(&payload);
    Ok(ColumnMetaData {
        type_: physical_code(field.physical),
        encodings: vec![ENCODING_PLAIN, ENCODING_RLE],
        path_in_schema: vec![field.name.clone()],
        codec: CODEC_UNCOMPRESSED,
        num_values: i64::try_from(rows).map_err(|_| "row count overflows".to_string())?,
        total_uncompressed_size: i64::try_from(total)
            .map_err(|_| "column chunk size overflows".to_string())?,
        total_compressed_size: i64::try_from(total)
            .map_err(|_| "column chunk size overflows".to_string())?,
        data_page_offset: offset,
        dictionary_page_offset: None,
        statistics: Some(build_statistics(
            &values,
            i64::try_from(null_count).map_err(|_| "null count overflows".to_string())?,
        )),
    })
}

/// Reads one column chunk, appending its pages into a column.
fn read_column_chunk(
    bytes: &[u8],
    metadata: &ColumnMetaData,
    field: &Field,
    limit: usize,
) -> Result<Column, String> {
    if metadata.dictionary_page_offset.is_some() {
        return Err(format!(
            "column `{}` uses dictionary encoding, which is not supported",
            field.name
        ));
    }
    let total = usize::try_from(metadata.num_values)
        .map_err(|_| format!("column `{}` has a negative value count", field.name))?;
    if total > MAX_ROWS {
        return Err(format!("column `{}` declares too many values", field.name));
    }
    if total == 0 {
        return Ok(empty_column(field));
    }
    if metadata.data_page_offset < 4 {
        return Err(format!(
            "column `{}` has an invalid data page offset",
            field.name
        ));
    }
    let mut pos = usize::try_from(metadata.data_page_offset)
        .map_err(|_| format!("column `{}` has an invalid page offset", field.name))?;
    let max_definition = u32::from(field.optional);
    let mut decoded = match field.physical {
        | Physical::Boolean => Decoded::Bool(Vec::with_capacity(total)),
        | Physical::Int32 => Decoded::I32(Vec::with_capacity(total)),
        | Physical::Int64 => Decoded::I64(Vec::with_capacity(total)),
        | Physical::Float => Decoded::F32(Vec::with_capacity(total)),
        | Physical::Double => Decoded::F64(Vec::with_capacity(total)),
        | Physical::ByteArray => Decoded::Bytes(Vec::with_capacity(total)),
    };
    let mut nulls: Vec<bool> = Vec::with_capacity(total);
    let mut collected = 0usize;
    while collected < total {
        if pos >= limit || pos >= bytes.len() {
            return Err(format!(
                "column `{}` ended before all values were read",
                field.name
            ));
        }
        let slice = bytes
            .get(pos..)
            .ok_or_else(|| "page header is truncated".to_string())?;
        let (header, consumed) = thrift::decode_page_header(slice)?;
        let payload_start = pos
            .checked_add(consumed)
            .ok_or_else(|| "page header length overflows".to_string())?;
        if header.type_ != PAGE_DATA {
            return Err(format!(
                "column `{}` uses unsupported page type {}",
                field.name, header.type_
            ));
        }
        let compressed_len = usize::try_from(header.compressed_page_size)
            .map_err(|_| format!("column `{}` has a negative page size", field.name))?;
        let payload_end = payload_start
            .checked_add(compressed_len)
            .ok_or_else(|| "page length overflows".to_string())?;
        if payload_end > limit {
            return Err(format!(
                "column `{}` page extends past the data section",
                field.name
            ));
        }
        let payload = bytes
            .get(payload_start..payload_end)
            .ok_or_else(|| format!("column `{}` page is truncated", field.name))?;
        let data = match metadata.codec {
            | CODEC_UNCOMPRESSED => payload.to_vec(),
            | CODEC_SNAPPY => snappy::decompress(payload)?,
            | other => {
                return Err(format!(
                    "column `{}` uses unsupported compression codec {other}",
                    field.name
                ));
            },
        };
        let page = header
            .data_page
            .ok_or_else(|| format!("column `{}` data page has no page header", field.name))?;
        if page.encoding != ENCODING_PLAIN {
            return Err(format!(
                "column `{}` uses unsupported value encoding {}",
                field.name, page.encoding
            ));
        }
        let page_values = usize::try_from(page.num_values)
            .map_err(|_| format!("column `{}` has a negative page count", field.name))?;
        if page_values == 0 {
            return Err(format!(
                "column `{}` contains an empty data page",
                field.name
            ));
        }
        let next = collected
            .checked_add(page_values)
            .ok_or_else(|| "page value count overflows".to_string())?;
        if next > total {
            return Err(format!(
                "column `{}` page declares more values than the chunk",
                field.name
            ));
        }
        let mut offset = 0usize;
        let levels: Vec<bool> = if max_definition > 0 {
            if page.definition_level_encoding != ENCODING_RLE {
                return Err(format!(
                    "column `{}` uses unsupported definition level encoding",
                    field.name
                ));
            }
            let length_bytes = data.get(0..4).ok_or_else(|| {
                format!("column `{}` definition levels are truncated", field.name)
            })?;
            let mut raw_len = [0u8; 4];
            raw_len.copy_from_slice(length_bytes);
            let level_len = u32::from_le_bytes(raw_len) as usize;
            offset = 4usize
                .checked_add(level_len)
                .ok_or_else(|| "definition level length overflows".to_string())?;
            let level_data = data.get(4..offset).ok_or_else(|| {
                format!("column `{}` definition levels exceed the page", field.name)
            })?;
            plain::decode_levels(level_data, 1, page_values)?
        } else {
            vec![false; page_values]
        };
        let present = levels.iter().filter(|value| !**value).count();
        let rest = data
            .get(offset..)
            .ok_or_else(|| format!("column `{}` page payload is truncated", field.name))?;
        append_values(&mut decoded, field, rest, present)?;
        nulls.extend_from_slice(&levels);
        collected = next;
        pos = payload_end;
    }
    if nulls.len() != total {
        return Err(format!(
            "column `{}` has {} level entries but {total} values",
            field.name,
            nulls.len()
        ));
    }
    Ok(decoded_into_column(decoded, field, nulls))
}

/// Decodes `count` values of the field's physical type and appends them.
fn append_values(
    decoded: &mut Decoded,
    field: &Field,
    bytes: &[u8],
    count: usize,
) -> Result<(), String> {
    match decoded {
        | Decoded::Bool(values) => {
            if field.physical != Physical::Boolean {
                return Err("column type mismatch".to_string());
            }
            values.extend(plain::decode_bool(bytes, count)?);
        },
        | Decoded::I32(values) => {
            if field.physical != Physical::Int32 {
                return Err("column type mismatch".to_string());
            }
            values.extend(plain::decode_i32(bytes, count)?);
        },
        | Decoded::I64(values) => {
            if field.physical != Physical::Int64 {
                return Err("column type mismatch".to_string());
            }
            values.extend(plain::decode_i64(bytes, count)?);
        },
        | Decoded::F32(values) => {
            if field.physical != Physical::Float {
                return Err("column type mismatch".to_string());
            }
            values.extend(plain::decode_f32(bytes, count)?);
        },
        | Decoded::F64(values) => {
            if field.physical != Physical::Double {
                return Err("column type mismatch".to_string());
            }
            values.extend(plain::decode_f64(bytes, count)?);
        },
        | Decoded::Bytes(values) => {
            if field.physical != Physical::ByteArray {
                return Err("column type mismatch".to_string());
            }
            for value in plain::decode_bytes(bytes, count)? {
                let text = String::from_utf8(value.to_vec())
                    .map_err(|_| format!("column `{}` contains invalid utf-8", field.name))?;
                values.push(text);
            }
        },
    }
    Ok(())
}

/// Attaches the null mask to a decoded column.
///
/// The decoded values only cover the non-missing slots, so they are expanded
/// back to the full row count with placeholder entries under the mask.
fn decoded_into_column(
    decoded: Decoded,
    field: &Field,
    nulls: Vec<bool>,
) -> Column {
    if !field.optional {
        return decoded_column(decoded);
    }
    let mask = nulls;
    let base = match decoded {
        | Decoded::Bool(values) => Column::Bool(expand(&values, &mask)),
        | Decoded::I32(values) => {
            Column::Int(expand(&values, &mask).into_iter().map(i64::from).collect())
        },
        | Decoded::I64(values) => Column::Int(expand(&values, &mask)),
        | Decoded::F32(values) => {
            Column::Float(expand(&values, &mask).into_iter().map(f64::from).collect())
        },
        | Decoded::F64(values) => Column::Float(expand(&values, &mask)),
        | Decoded::Bytes(values) => Column::Str(expand(&values, &mask)),
    };
    Column::Nullable(Box::new(base), mask)
}

/// Expands present values into a full-length vector, filling null slots with
/// the type's default.
fn expand<T: Clone + Default>(
    values: &[T],
    mask: &[bool],
) -> Vec<T> {
    let mut out = Vec::with_capacity(mask.len());
    let mut source = 0usize;
    for missing in mask {
        if *missing {
            out.push(T::default());
        } else {
            out.push(values.get(source).cloned().unwrap_or_default());
            source += 1;
        }
    }
    out
}

/// Converts decoded values into a non-nullable column.
fn decoded_column(decoded: Decoded) -> Column {
    match decoded {
        | Decoded::Bool(values) => Column::Bool(values),
        | Decoded::I32(values) => Column::Int(values.into_iter().map(i64::from).collect()),
        | Decoded::I64(values) => Column::Int(values),
        | Decoded::F32(values) => Column::Float(values.into_iter().map(f64::from).collect()),
        | Decoded::F64(values) => Column::Float(values),
        | Decoded::Bytes(values) => Column::Str(values),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Builds a unique temporary path.
    fn temp_path(tag: &str) -> std::path::PathBuf {
        std::env::temp_dir().join(format!(
            "rssn_parquet_{tag}_{}_{}.parquet",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .map_or(0, |duration| duration.as_nanos())
        ))
    }

    /// Builds a 1000-row frame covering all five physical types.
    fn sample() -> (Schema, Vec<Column>) {
        let rows = 1000usize;
        let mask: Vec<bool> = (0..rows).map(|row| row % 7 == 0).collect();
        let booleans: Vec<bool> = (0..rows).map(|row| row % 3 == 0).collect();
        let small: Vec<i64> = (0..rows).map(|row| row as i64 - 500).collect();
        let big: Vec<i64> = (0..rows).map(|row| row as i64 * 1_000_000).collect();
        let single: Vec<f64> = (0..rows).map(|row| row as f64 * 0.5).collect();
        let double: Vec<f64> = (0..rows).map(|row| row as f64 * -1.25).collect();
        let strings: Vec<String> = (0..rows).map(|row| format!("row-{row}")).collect();
        let schema = Schema {
            fields: vec![
                Field {
                    name: "flag".to_string(),
                    physical: Physical::Boolean,
                    optional: true,
                },
                Field {
                    name: "small".to_string(),
                    physical: Physical::Int32,
                    optional: true,
                },
                Field {
                    name: "big".to_string(),
                    physical: Physical::Int64,
                    optional: false,
                },
                Field {
                    name: "single".to_string(),
                    physical: Physical::Float,
                    optional: true,
                },
                Field {
                    name: "double".to_string(),
                    physical: Physical::Double,
                    optional: false,
                },
                Field {
                    name: "text".to_string(),
                    physical: Physical::ByteArray,
                    optional: true,
                },
            ],
        };
        let columns = vec![
            Column::Nullable(Box::new(Column::Bool(booleans)), mask.clone()),
            Column::Nullable(Box::new(Column::Int(small)), mask.clone()),
            Column::Int(big),
            Column::Nullable(Box::new(Column::Float(single)), mask.clone()),
            Column::Float(double),
            Column::Nullable(Box::new(Column::Str(strings)), mask),
        ];
        (schema, columns)
    }

    #[test]
    fn round_trip_all_physical_types() {
        let (schema, columns) = sample();
        let path = temp_path("round_trip");
        write_parquet(&path, &schema, &columns).expect("write");
        let (read_schema, read_columns) = read_parquet(&path).expect("read");
        let _ = std::fs::remove_file(&path);
        assert_eq!(read_schema, schema);
        assert_eq!(read_columns.len(), columns.len());
        for (index, (expected, actual)) in columns.iter().zip(&read_columns).enumerate() {
            assert_eq!(actual.len(), expected.len(), "column {index} length");
            assert_eq!(actual.dtype(), expected.dtype(), "column {index} type");
            assert_eq!(
                actual.is_nullable(),
                expected.is_nullable(),
                "column {index} nullability"
            );
            for row in 0..expected.len() {
                assert_eq!(
                    actual.is_null(row),
                    expected.is_null(row),
                    "column {index} row {row} null"
                );
                assert_eq!(
                    actual.get_f64(row),
                    expected.get_f64(row),
                    "column {index} row {row} value"
                );
            }
        }
        assert_eq!(read_columns.first().map(Column::len), Some(1000));
        assert!(read_columns.first().is_some_and(|column| column.is_null(0)));
        assert!(
            read_columns
                .first()
                .is_some_and(|column| !column.is_null(1))
        );
    }

    #[test]
    fn utf8_column_round_trips() {
        let schema = Schema {
            fields: vec![Field {
                name: "text".to_string(),
                physical: Physical::ByteArray,
                optional: false,
            }],
        };
        let columns = vec![Column::Str(vec![
            "héllo".to_string(),
            "数据".to_string(),
            "🚀 rocket".to_string(),
        ])];
        let bytes = write_parquet_bytes(&schema, &columns).expect("write");
        let (_, read) = read_parquet_bytes(&bytes).expect("read");
        assert_eq!(read, columns);
    }

    #[test]
    fn truncation_never_panics() {
        let (schema, columns) = sample();
        let bytes = write_parquet_bytes(&schema, &columns).expect("write");
        let cuts = [
            bytes.len() / 3,
            bytes.len() / 2,
            bytes.len().saturating_sub(6),
        ];
        for cut in cuts {
            let slice = bytes.get(..cut).expect("slice");
            assert!(read_parquet_bytes(slice).is_err(), "cut at {cut}");
        }
        assert!(read_parquet_bytes(&[]).is_err());
        assert!(read_parquet_bytes(b"PAR1").is_err());
    }

    #[test]
    fn corrupt_magic_is_rejected() {
        let (schema, columns) = sample();
        let mut bytes = write_parquet_bytes(&schema, &columns).expect("write");
        if let Some(first) = bytes.first_mut() {
            *first = b'X';
        }
        assert!(read_parquet_bytes(&bytes).is_err());
    }

    #[test]
    fn schema_mismatch_is_rejected() {
        let (schema, columns) = sample();
        let wrong = Schema {
            fields: vec![Field {
                name: "flag".to_string(),
                physical: Physical::Boolean,
                optional: false,
            }],
        };
        assert!(write_parquet_bytes(&wrong, &columns).is_err());
        assert!(write_parquet_bytes(&schema, &columns[..1]).is_err());
    }

    #[test]
    fn empty_frame_round_trips() {
        let schema = Schema {
            fields: vec![Field {
                name: "value".to_string(),
                physical: Physical::Int64,
                optional: true,
            }],
        };
        let columns = vec![Column::Nullable(
            Box::new(Column::Int(Vec::new())),
            Vec::new(),
        )];
        let bytes = write_parquet_bytes(&schema, &columns).expect("write");
        let (read_schema, read_columns) = read_parquet_bytes(&bytes).expect("read");
        assert_eq!(read_schema, schema);
        assert_eq!(read_columns, columns);
    }

    #[test]
    fn file_round_trip() {
        let schema = Schema {
            fields: vec![Field {
                name: "n".to_string(),
                physical: Physical::Int32,
                optional: false,
            }],
        };
        let columns = vec![Column::Int(vec![1, 2, 3])];
        let path = temp_path("file");
        write_parquet(&path, &schema, &columns).expect("write");
        let (_, read) = read_parquet(&path).expect("read");
        let _ = std::fs::remove_file(&path);
        assert_eq!(read, columns);
    }
}
