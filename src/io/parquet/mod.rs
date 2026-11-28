//! Apache Parquet reader and writer.
//!
//! The writer emits one or more row groups with one or more data pages per
//! column, encoded as `PLAIN` or, when a column chunk repeats values often
//! enough to pay for a dictionary, as a `PLAIN` dictionary page plus
//! `RLE_DICTIONARY` data pages ([`WriteOptions::use_dictionary`]). Pages are
//! optionally compressed with Snappy, gzip, Zstandard or LZ4
//! ([`WriteOptions::compression`]).
//! Optional columns carry Parquet-spec RLE/bit-packed definition levels (`1`
//! for present, `0` for null).
//!
//! The reader validates every offset and length against the buffer before
//! slicing and supports:
//! - Multiple row groups and multiple data pages per column chunk
//! - Data Page v1 and v2 headers (including the v2 RLE boolean encoding and
//!   the split level/value compression of v2 pages)
//! - `UNCOMPRESSED`, `SNAPPY`, `GZIP`, `ZSTD` and `LZ4_RAW` page compression
//! - `PLAIN`, `PLAIN_DICTIONARY`, `RLE_DICTIONARY`, `DELTA_BINARY_PACKED`,
//!   `DELTA_LENGTH_BYTE_ARRAY`, and `BYTE_STREAM_SPLIT` encodings
//! - All Parquet physical types: `BOOLEAN`, `INT32`, `INT64`, `INT96` (mapped
//!   to Unix nanoseconds `i64`), `FLOAT`, `DOUBLE`, `BYTE_ARRAY` (UTF-8
//!   strings when annotated with `UTF8`, raw `Column::Bytes` otherwise), and
//!   `FIXED_LEN_BYTE_ARRAY`
//! - Nested optional and required struct groups (flattened to dotted column
//!   names such as `parent.child` with full multi-bit definition-level
//!   decoding)
//! - Repeated fields: `LIST`, `MAP` and bare repeated primitives become a
//!   [`Column::List`] of the element column, with null lists, empty lists and
//!   null elements kept apart. A `MAP` is exposed as its `key` and `value`
//!   lists, one column each, and lists nest to any depth.
//!
//! # Limitations
//!
//! - `DECIMAL` columns are exposed as `f64` divided by `10^scale`; the other
//!   logical types (dates, times, UUIDs) are exposed through their physical
//!   representation.
//! - The deprecated Hadoop LZ4 framing (codec 5) is rejected; the `LZ4_RAW`
//!   codec (7) that current writers emit is supported.
//!
//! [`WriteOptions::use_dictionary`]: crate::io::parquet::WriteOptions::use_dictionary
//! [`WriteOptions::compression`]: crate::io::parquet::WriteOptions::compression
//! [`Column::List`]: crate::io::frame::Column::List

pub mod bloom;
pub(crate) mod gzip;
pub(crate) mod lz4;
mod lzo;
mod plain;
pub(crate) mod snappy;
mod thrift;

pub use bloom::SplitBlockBloomFilter;
pub use thrift::ColumnIndex;
pub use thrift::OffsetIndex;
pub use thrift::PageLocation;

use std::ops::Range;
use std::path::Path;

use crate::io::frame::Column;
use crate::io::frame::DType;

use thrift::ColumnChunk;
use thrift::ColumnMetaData;
use thrift::DataPageHeader;
use thrift::DataPageHeaderV2;
use thrift::DictionaryPageHeader;
use thrift::FileMetaData;
use thrift::PageHeader;
use thrift::RowGroup;
use thrift::SchemaElement;
use thrift::Statistics;

/// The magic bytes at both ends of a Parquet file.
const MAGIC: &[u8; 4] = b"PAR1";
/// Upper bound on declared rows, guarding against hostile metadata.
const MAX_ROWS: usize = 100_000_000;
/// Page type code for a v1 data page.
const PAGE_DATA: i32 = 0;
/// Page type code for a v2 data page.
const PAGE_DATA_V2: i32 = 3;
/// Page type code for a dictionary page.
const PAGE_DICTIONARY: i32 = 2;
/// Encoding code for `PLAIN`.
const ENCODING_PLAIN: i32 = 0;
/// Encoding code for `PLAIN_DICTIONARY`.
const ENCODING_PLAIN_DICTIONARY: i32 = 2;
/// Encoding code for RLE.
const ENCODING_RLE: i32 = 3;
/// Encoding code for `DELTA_BINARY_PACKED`.
const ENCODING_DELTA_BINARY_PACKED: i32 = 5;
/// Encoding code for `DELTA_LENGTH_BYTE_ARRAY`.
const ENCODING_DELTA_LENGTH_BYTE_ARRAY: i32 = 6;
/// Encoding code for delta-encoded byte arrays.
const ENCODING_DELTA_BYTE_ARRAY: i32 = 7;
/// Encoding code for `RLE_DICTIONARY`.
const ENCODING_RLE_DICTIONARY: i32 = 8;
/// Encoding code for `BYTE_STREAM_SPLIT`.
const ENCODING_BYTE_STREAM_SPLIT: i32 = 9;
/// Codec code for no compression.
const CODEC_UNCOMPRESSED: i32 = 0;
/// Codec code for Snappy.
const CODEC_SNAPPY: i32 = 1;
/// Codec code for gzip.
const CODEC_GZIP: i32 = 2;
/// Codec code for LZO.
const CODEC_LZO: i32 = 3;
/// Codec code for Brotli.
const CODEC_BROTLI: i32 = 4;
/// Codec code for the deprecated Hadoop LZ4 framing.
const CODEC_LZ4: i32 = 5;
/// Codec code for Zstandard.
const CODEC_ZSTD: i32 = 6;
/// Codec code for a raw LZ4 block.
const CODEC_LZ4_RAW: i32 = 7;
/// Repetition code for a required field.
const REPETITION_REQUIRED: i32 = 0;
/// Repetition code for an optional field.
const REPETITION_OPTIONAL: i32 = 1;
/// Repetition code for a repeated field.
const REPETITION_REPEATED: i32 = 2;
/// Converted type code for UTF-8 strings.
const CONVERTED_UTF8: i32 = 0;
/// Converted type code for lists.
const CONVERTED_LIST: i32 = 3;
/// Converted type code for fixed-point decimals.
const CONVERTED_DECIMAL: i32 = 5;

/// Compression codec used when writing Parquet data pages.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Compression {
    /// Store page payloads uncompressed.
    #[default]
    Uncompressed,
    /// Compress each page payload with Snappy.
    Snappy,
    /// Compress each page payload with gzip (DEFLATE).
    Gzip,
    /// Compress each page payload with LZO (Hadoop-LZO framed).
    Lzo,
    /// Compress each page payload with Brotli.
    Brotli,
    /// Compress each page payload with Zstandard.
    Zstd,
    /// Compress each page payload with a raw LZ4 block.
    Lz4,
}

const fn codec_code(compression: Compression) -> i32 {
    match compression {
        | Compression::Uncompressed => CODEC_UNCOMPRESSED,
        | Compression::Snappy => CODEC_SNAPPY,
        | Compression::Gzip => CODEC_GZIP,
        | Compression::Lzo => CODEC_LZO,
        | Compression::Brotli => CODEC_BROTLI,
        | Compression::Zstd => CODEC_ZSTD,
        | Compression::Lz4 => CODEC_LZ4_RAW,
    }
}

/// Options controlling Parquet writing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct WriteOptions {
    /// Compression codec for data pages.
    pub compression: Compression,
    /// Maximum number of rows per row group (must be at least 1).
    pub max_row_group_rows: usize,
    /// Maximum number of rows per data page within a column chunk (must be at least 1).
    pub max_page_rows: usize,
    /// Emit a dictionary page and `RLE_DICTIONARY` data pages for a column
    /// chunk whose values repeat often enough to pay for the dictionary.
    pub use_dictionary: bool,
    /// Data page version to write: `1` keeps the levels inside the compressed
    /// payload, `2` stores them raw ahead of it.
    pub data_page_version: u8,
    /// Emit a Split Block Bloom Filter (`XXH64`) for each column chunk.
    pub write_bloom_filter: bool,
    /// Emit `ColumnIndex` and `OffsetIndex` (Parquet Page Index) for each column chunk.
    pub write_page_index: bool,
}

impl Default for WriteOptions {
    fn default() -> Self {
        Self {
            compression: Compression::Uncompressed,
            max_row_group_rows: 1_000_000,
            max_page_rows: 100_000,
            data_page_version: 1,
            use_dictionary: true,
            write_bloom_filter: false,
            write_page_index: false,
        }
    }
}

/// The physical storage type of a Parquet column.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Physical {
    /// One bit per value.
    Boolean,
    /// Four-byte signed integers.
    Int32,
    /// Eight-byte signed integers.
    Int64,
    /// Twelve-byte legacy Parquet timestamps, mapped to Unix nanoseconds (`i64`).
    Int96,
    /// Four-byte floats.
    Float,
    /// Eight-byte floats.
    Double,
    /// Length-prefixed byte arrays (UTF-8 strings or raw binary).
    ByteArray,
    /// Fixed-length byte arrays of the given byte width.
    FixedLenByteArray(usize),
}

/// One leaf column of a Parquet schema.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Field {
    /// Column name (dotted path for nested struct leaves).
    pub name: String,
    /// Physical storage type.
    pub physical: Physical,
    /// Whether the column may contain missing values.
    pub optional: bool,
}

/// A Parquet schema as a list of leaf fields.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Schema {
    /// Leaf fields, in file order.
    pub fields: Vec<Field>,
}

/// One repeated node between the root and a leaf.
#[derive(Debug, Clone, PartialEq, Eq)]
struct ListSpec {
    /// Definition level at which the list itself exists.
    present: u32,
    /// Whether the list may be null.
    optional: bool,
    /// Whether an element of the list may be missing.
    element_optional: bool,
}

/// Internal metadata for a leaf column parsed from the schema tree.
struct LeafInfo {
    /// Public field descriptor.
    field: Field,
    /// Path segments from root to leaf.
    path: Vec<String>,
    /// Maximum definition level along the path.
    max_def_level: u32,
    /// Maximum repetition level along the path.
    max_rep_level: u32,
    /// Repeated nodes enclosing the leaf, outermost first.
    lists: Vec<ListSpec>,
    /// Whether the element is annotated with `CONVERTED_UTF8`.
    is_utf8: bool,
    /// Scale of a `DECIMAL` annotation; values are exposed as `f64` divided by
    /// `10^scale`.
    decimal_scale: Option<u32>,
}

impl LeafInfo {
    /// Whether the leaf sits under a repeated node, making it a list.
    fn is_list(&self) -> bool {
        !self.lists.is_empty()
    }
}

/// A decoded column, before the null mask is attached.
enum Decoded {
    /// Boolean values.
    Bool(Vec<bool>),
    /// 32-bit integers.
    I32(Vec<i32>),
    /// 64-bit integers.
    I64(Vec<i64>),
    /// 96-bit timestamps converted to Unix nanoseconds.
    Int96(Vec<i64>),
    /// 32-bit floats.
    F32(Vec<f32>),
    /// 64-bit floats.
    F64(Vec<f64>),
    /// UTF-8 strings.
    Str(Vec<String>),
    /// Raw byte arrays (non-UTF-8 `BYTE_ARRAY` or `FIXED_LEN_BYTE_ARRAY`).
    RawBytes(Vec<Vec<u8>>),
}

impl Decoded {
    /// Number of values held.
    fn len(&self) -> usize {
        match self {
            | Decoded::Bool(items) => items.len(),
            | Decoded::I32(items) => items.len(),
            | Decoded::I64(items) | Decoded::Int96(items) => items.len(),
            | Decoded::F32(items) => items.len(),
            | Decoded::F64(items) => items.len(),
            | Decoded::Str(items) => items.len(),
            | Decoded::RawBytes(items) => items.len(),
        }
    }
}

/// Builds a dictionary of distinct values and one index per input value.
///
/// Returns `None` for booleans, which Parquet does not allow to be
/// dictionary-encoded.
fn dictionary_encode(values: &Decoded) -> Option<(Decoded, Vec<u32>)> {
    match values {
        | Decoded::Bool(_) => None,
        | Decoded::I32(items) => {
            let (unique, indices) = index_distinct(items);
            Some((Decoded::I32(unique), indices))
        },
        | Decoded::I64(items) => {
            let (unique, indices) = index_distinct(items);
            Some((Decoded::I64(unique), indices))
        },
        | Decoded::Int96(items) => {
            let (unique, indices) = index_distinct(items);
            Some((Decoded::Int96(unique), indices))
        },
        | Decoded::F32(items) => {
            let keys: Vec<u32> = items.iter().map(|value| value.to_bits()).collect();
            let (unique_keys, indices) = index_distinct(&keys);
            let unique = unique_keys
                .iter()
                .map(|bits| f32::from_bits(*bits))
                .collect();
            Some((Decoded::F32(unique), indices))
        },
        | Decoded::F64(items) => {
            let keys: Vec<u64> = items.iter().map(|value| value.to_bits()).collect();
            let (unique_keys, indices) = index_distinct(&keys);
            let unique = unique_keys
                .iter()
                .map(|bits| f64::from_bits(*bits))
                .collect();
            Some((Decoded::F64(unique), indices))
        },
        | Decoded::Str(items) => {
            let (unique, indices) = index_distinct(items);
            Some((Decoded::Str(unique), indices))
        },
        | Decoded::RawBytes(items) => {
            let (unique, indices) = index_distinct(items);
            Some((Decoded::RawBytes(unique), indices))
        },
    }
}

/// Returns the distinct values of `items` in first-appearance order and the
/// index of each input value in that list.
fn index_distinct<T: Clone + Eq + std::hash::Hash>(items: &[T]) -> (Vec<T>, Vec<u32>) {
    let mut position: std::collections::HashMap<&T, u32> = std::collections::HashMap::new();
    let mut unique: Vec<T> = Vec::new();
    let mut indices = Vec::with_capacity(items.len());
    for item in items {
        let next = u32::try_from(unique.len()).unwrap_or(u32::MAX);
        let slot = *position.entry(item).or_insert_with(|| {
            unique.push(item.clone());
            next
        });
        indices.push(slot);
    }
    (unique, indices)
}

/// Writes a frame to a Parquet file using default [`WriteOptions`].
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
    write_parquet_with_options(path, schema, columns, &WriteOptions::default())
}

/// Writes a frame to a Parquet file using explicit [`WriteOptions`].
///
/// # Errors
///
/// Returns an error when options are invalid, the schema does not match the
/// columns, a value does not fit its physical type, or the file cannot be
/// written.
pub fn write_parquet_with_options(
    path: &Path,
    schema: &Schema,
    columns: &[Column],
    options: &WriteOptions,
) -> Result<(), String> {
    let bytes = write_parquet_bytes_with_options(schema, columns, options)?;
    std::fs::write(path, bytes).map_err(|error| error.to_string())
}

/// Writes a frame to an in-memory Parquet file using default [`WriteOptions`].
///
/// # Errors
///
/// Returns an error when the schema does not match the columns or a value
/// does not fit its physical type.
pub fn write_parquet_bytes(
    schema: &Schema,
    columns: &[Column],
) -> Result<Vec<u8>, String> {
    write_parquet_bytes_with_options(schema, columns, &WriteOptions::default())
}

/// Writes a frame to an in-memory Parquet file using explicit [`WriteOptions`].
///
/// # Errors
///
/// Returns an error when options are invalid, the schema does not match the
/// columns, or a value does not fit its physical type.
pub fn write_parquet_bytes_with_options(
    schema: &Schema,
    columns: &[Column],
    options: &WriteOptions,
) -> Result<Vec<u8>, String> {
    if options.max_row_group_rows == 0 {
        return Err("max_row_group_rows must be at least 1".to_string());
    }
    if options.max_page_rows == 0 {
        return Err("max_page_rows must be at least 1".to_string());
    }
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

    let rg_ranges: Vec<Range<usize>> = if rows == 0 {
        std::iter::once(0..0).collect()
    } else {
        let mut ranges = Vec::new();
        let mut start = 0usize;
        while start < rows {
            let end = start.saturating_add(options.max_row_group_rows).min(rows);
            ranges.push(start..end);
            start = end;
        }
        ranges
    };

    let mut row_groups = Vec::with_capacity(rg_ranges.len());
    let mut pending_indexes: Vec<Vec<(Option<ColumnIndex>, Option<OffsetIndex>)>> =
        Vec::with_capacity(rg_ranges.len());
    for rg_range in rg_ranges {
        let rg_offset = i64::try_from(file.len())
            .map_err(|_| "parquet file has grown too large".to_string())?;
        let mut chunks = Vec::with_capacity(schema.fields.len());
        let mut rg_indexes = Vec::with_capacity(schema.fields.len());
        let mut total_uncompressed = 0i64;
        let mut total_compressed = 0i64;
        for (field, column) in schema.fields.iter().zip(columns) {
            let offset = i64::try_from(file.len())
                .map_err(|_| "parquet file has grown too large".to_string())?;
            let (metadata, col_idx, off_idx) = write_column_chunk_ext(
                field,
                column,
                rg_range.clone(),
                offset,
                options,
                &mut file,
            )?;
            total_uncompressed = total_uncompressed
                .checked_add(metadata.total_uncompressed_size)
                .ok_or_else(|| "parquet file size overflows".to_string())?;
            total_compressed = total_compressed
                .checked_add(metadata.total_compressed_size)
                .ok_or_else(|| "parquet file size overflows".to_string())?;
            chunks.push(ColumnChunk {
                file_offset: offset,
                meta_data: Some(metadata),
                offset_index_offset: None,
                offset_index_length: None,
                column_index_offset: None,
                column_index_length: None,
            });
            rg_indexes.push((col_idx, off_idx));
        }
        pending_indexes.push(rg_indexes);
        row_groups.push(RowGroup {
            columns: chunks,
            total_byte_size: total_uncompressed,
            num_rows: i64::try_from(rg_range.len())
                .map_err(|_| "row count overflows".to_string())?,
            file_offset: Some(rg_offset),
            total_compressed_size: Some(total_compressed),
        });
    }

    // Write Page Index (ColumnIndex followed by OffsetIndex) after all row groups
    if options.write_page_index {
        for (rg, rg_idx) in row_groups.iter_mut().zip(pending_indexes.iter()) {
            for (chunk, (col_idx_opt, off_idx_opt)) in rg.columns.iter_mut().zip(rg_idx.iter()) {
                if let Some(col_idx) = col_idx_opt {
                    let bytes = thrift::encode_column_index(col_idx);
                    chunk.column_index_offset = Some(file.len() as i64);
                    chunk.column_index_length = Some(bytes.len() as i32);
                    file.extend_from_slice(&bytes);
                }
                if let Some(off_idx) = off_idx_opt {
                    let bytes = thrift::encode_offset_index(off_idx);
                    chunk.offset_index_offset = Some(file.len() as i64);
                    chunk.offset_index_length = Some(bytes.len() as i32);
                    file.extend_from_slice(&bytes);
                }
            }
        }
    }

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
    for (field, column) in schema.fields.iter().zip(columns) {
        let type_length = match field.physical {
            | Physical::FixedLenByteArray(width) => {
                Some(i32::try_from(width).map_err(|_| "fixed length overflows i32".to_string())?)
            },
            | _ => None,
        };
        if column.list_values().is_some() {
            // Nested three-level LIST: one annotated group and one repeated
            // group per level, then the element leaf.
            let specs = list_specs(column, field.optional);
            let mut leaf = column;
            while let Some(values) = leaf.list_values() {
                leaf = values;
            }
            push_list_schema(
                &mut elements,
                &field.name,
                &specs,
                1,
                field,
                type_length,
                leaf.dtype(),
            );
            continue;
        }
        let converted_type =
            if field.physical == Physical::ByteArray && column.dtype() == DType::Str {
                Some(CONVERTED_UTF8)
            } else {
                None
            };
        elements.push(SchemaElement {
            type_: Some(physical_code(field.physical)),
            type_length,
            repetition_type: Some(if field.optional {
                REPETITION_OPTIONAL
            } else {
                REPETITION_REQUIRED
            }),
            name: field.name.clone(),
            num_children: None,
            converted_type,
            scale: None,
            precision: None,
            field_id: None,
        });
    }
    let metadata = FileMetaData {
        version: 1,
        schema: elements,
        num_rows: i64::try_from(rows).map_err(|_| "row count overflows".to_string())?,
        row_groups,
    };
    let metadata_bytes = thrift::encode_file_metadata(&metadata);
    let metadata_len = u32::try_from(metadata_bytes.len())
        .map_err(|_| "parquet metadata is too large".to_string())?;
    file.extend_from_slice(&metadata_bytes);
    file.extend_from_slice(&metadata_len.to_le_bytes());
    file.extend_from_slice(MAGIC);
    Ok(file)
}

/// Combined Page Index for one column chunk.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct PageIndex {
    /// Page-level min/max/null statistics, if present.
    pub column_index: Option<ColumnIndex>,
    /// Page byte offsets and row indices, if present.
    pub offset_index: Option<OffsetIndex>,
}

fn read_footer_metadata(bytes: &[u8]) -> Result<FileMetaData, String> {
    if bytes.len() < 12 {
        return Err("parquet file is too short".to_string());
    }
    if bytes.get(..4) != Some(MAGIC) || bytes.get(bytes.len() - 4..) != Some(MAGIC) {
        return Err("file does not have Parquet PAR1 magic".to_string());
    }
    let len_start = bytes.len() - 8;
    let footer_len = u32::from_le_bytes([
        bytes[len_start],
        bytes[len_start + 1],
        bytes[len_start + 2],
        bytes[len_start + 3],
    ]) as usize;
    if footer_len > len_start.saturating_sub(4) {
        return Err("parquet footer length exceeds file size".to_string());
    }
    let meta_bytes = &bytes[len_start - footer_len..len_start];
    thrift::decode_file_metadata(meta_bytes)
}

/// Reads all Bloom filters from a Parquet buffer, indexed `[row_group][column]`.
///
/// # Errors
/// Fails if the Parquet footer or a Bloom filter header/bitset is malformed.
pub fn read_parquet_bloom_filters(
    bytes: &[u8]
) -> Result<Vec<Vec<Option<SplitBlockBloomFilter>>>, String> {
    let metadata = read_footer_metadata(bytes)?;
    let mut out = Vec::with_capacity(metadata.row_groups.len());
    for rg in &metadata.row_groups {
        let mut row = Vec::with_capacity(rg.columns.len());
        for col in &rg.columns {
            let Some(meta) = &col.meta_data else {
                row.push(None);
                continue;
            };
            let Some(off) = meta.bloom_filter_offset else {
                row.push(None);
                continue;
            };
            let start =
                usize::try_from(off).map_err(|_| "negative bloom_filter_offset".to_owned())?;
            let slice = bytes
                .get(start..)
                .ok_or_else(|| "bloom_filter_offset out of range".to_owned())?;
            let (num_bytes, hdr_len) = thrift::parse_bloom_filter_header(slice)?;
            let bitset = slice
                .get(hdr_len..hdr_len + num_bytes)
                .ok_or_else(|| "bloom filter bitset truncated".to_owned())?;
            row.push(Some(SplitBlockBloomFilter::from_bytes(bitset)?));
        }
        out.push(row);
    }
    Ok(out)
}

/// Reads all Page Indexes (`ColumnIndex` and `OffsetIndex`) from a Parquet buffer,
/// indexed `[row_group][column]`.
///
/// # Errors
/// Fails if the Parquet footer or any index structure is malformed.
pub fn read_parquet_page_index(bytes: &[u8]) -> Result<Vec<Vec<PageIndex>>, String> {
    let metadata = read_footer_metadata(bytes)?;
    let mut out = Vec::with_capacity(metadata.row_groups.len());
    for rg in &metadata.row_groups {
        let mut row = Vec::with_capacity(rg.columns.len());
        for col in &rg.columns {
            let column_index = match (col.column_index_offset, col.column_index_length) {
                | (Some(off), Some(len)) => {
                    let s = usize::try_from(off).map_err(|_| "negative column_index_offset")?;
                    let l = usize::try_from(len).map_err(|_| "negative column_index_length")?;
                    let buf = bytes
                        .get(s..s + l)
                        .ok_or_else(|| "column_index out of range".to_owned())?;
                    Some(thrift::parse_column_index(buf)?)
                },
                | _ => None,
            };
            let offset_index = match (col.offset_index_offset, col.offset_index_length) {
                | (Some(off), Some(len)) => {
                    let s = usize::try_from(off).map_err(|_| "negative offset_index_offset")?;
                    let l = usize::try_from(len).map_err(|_| "negative offset_index_length")?;
                    let buf = bytes
                        .get(s..s + l)
                        .ok_or_else(|| "offset_index out of range".to_owned())?;
                    Some(thrift::parse_offset_index(buf)?)
                },
                | _ => None,
            };
            row.push(PageIndex {
                column_index,
                offset_index,
            });
        }
        out.push(row);
    }
    Ok(out)
}

/// Reads a Parquet file.
///
/// # Errors
///
/// Returns an error when the file cannot be read or is malformed.
pub fn read_parquet(path: &Path) -> Result<(Schema, Vec<Column>), String> {
    let bytes = std::fs::read(path).map_err(|error| error.to_string())?;
    read_parquet_bytes(&bytes)
}

/// Reads a Parquet file from memory.
///
/// # Errors
///
/// Returns an error when the buffer is not a valid Parquet file or uses an
/// unsupported feature (such as repeated fields or unsupported codecs).
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
    let leaves = parse_schema(&metadata)?;
    let mut chunks: Vec<Vec<Column>> = vec![Vec::new(); leaves.len()];
    let mut total_rows = 0usize;
    let mut page_decoder = PageDecoder::default();
    for row_group in &metadata.row_groups {
        if row_group.num_rows < 0 || row_group.num_rows > MAX_ROWS as i64 {
            return Err("row group declares an impossible row count".to_string());
        }
        let rows = usize::try_from(row_group.num_rows)
            .map_err(|_| "row group row count is out of range".to_string())?;
        if row_group.columns.len() != leaves.len() {
            return Err(format!(
                "row group has {} columns but the schema has {}",
                row_group.columns.len(),
                leaves.len()
            ));
        }
        for (index, chunk) in row_group.columns.iter().enumerate() {
            let leaf = leaves
                .get(index)
                .ok_or_else(|| "row group column index is out of range".to_string())?;
            let chunk_metadata = chunk
                .meta_data
                .as_ref()
                .ok_or_else(|| format!("column chunk `{}` is missing metadata", leaf.field.name))?;
            if chunk_metadata.path_in_schema != leaf.path {
                return Err(format!(
                    "column chunk path does not match field `{}`",
                    leaf.field.name
                ));
            }
            let column = read_column_chunk(
                bytes,
                chunk_metadata,
                leaf,
                metadata_start,
                &mut page_decoder,
            )?;
            if column.len() != rows {
                return Err(format!(
                    "column `{}` has {} values but its row group has {rows} rows",
                    leaf.field.name,
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
    let mut columns = Vec::with_capacity(leaves.len());
    for (index, leaf) in leaves.iter().enumerate() {
        let parts = chunks
            .get_mut(index)
            .ok_or_else(|| "column index is out of range".to_string())?;
        let column = if parts.is_empty() {
            empty_column(leaf)
        } else {
            merge_chunks(std::mem::take(parts))?
        };
        columns.push(column);
    }
    let fields = leaves.into_iter().map(|leaf| leaf.field).collect();
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
        if let Physical::FixedLenByteArray(width) = field.physical {
            if width == 0 {
                return Err(format!(
                    "column `{}` fixed-length byte array width must be positive",
                    field.name
                ));
            }
        }
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
        | Column::Int(_) => {
            matches!(
                physical,
                Physical::Int32 | Physical::Int64 | Physical::Int96
            )
        },
        | Column::Float(_) => matches!(physical, Physical::Float | Physical::Double),
        | Column::Bool(_) => physical == Physical::Boolean,
        | Column::Str(_) => physical == Physical::ByteArray,
        | Column::Bytes(_) => {
            matches!(
                physical,
                Physical::ByteArray | Physical::FixedLenByteArray(_)
            )
        },
        | Column::Nullable(inner, _) => column_matches(inner, physical),
        | Column::List { values, .. } => column_matches(values, physical),
    }
}

/// Emits the schema elements of one list level.
///
/// `name` is the group name at this level (`element` below the first), and the
/// element of a level is either the next level's group or the leaf.
fn push_list_schema(
    elements: &mut Vec<SchemaElement>,
    name: &str,
    specs: &[ListSpec],
    depth: usize,
    field: &Field,
    type_length: Option<i32>,
    leaf_dtype: DType,
) {
    let spec = &specs[depth - 1];
    elements.push(SchemaElement {
        type_: None,
        type_length: None,
        repetition_type: Some(if spec.optional {
            REPETITION_OPTIONAL
        } else {
            REPETITION_REQUIRED
        }),
        name: name.to_string(),
        num_children: Some(1),
        converted_type: Some(CONVERTED_LIST),
        scale: None,
        precision: None,
        field_id: None,
    });
    elements.push(SchemaElement {
        type_: None,
        type_length: None,
        repetition_type: Some(REPETITION_REPEATED),
        name: "list".to_string(),
        num_children: Some(1),
        converted_type: None,
        scale: None,
        precision: None,
        field_id: None,
    });
    if depth < specs.len() {
        push_list_schema(
            elements,
            "element",
            specs,
            depth + 1,
            field,
            type_length,
            leaf_dtype,
        );
    } else {
        elements.push(SchemaElement {
            type_: Some(physical_code(field.physical)),
            type_length,
            repetition_type: Some(if spec.element_optional {
                REPETITION_OPTIONAL
            } else {
                REPETITION_REQUIRED
            }),
            name: "element".to_string(),
            num_children: None,
            converted_type: if field.physical == Physical::ByteArray && leaf_dtype == DType::Str {
                Some(CONVERTED_UTF8)
            } else {
                None
            },
            scale: None,
            precision: None,
            field_id: None,
        });
    }
}

/// Maps a physical type to its Parquet type code.
fn physical_code(physical: Physical) -> i32 {
    match physical {
        | Physical::Boolean => 0,
        | Physical::Int32 => 1,
        | Physical::Int64 => 2,
        | Physical::Int96 => 3,
        | Physical::Float => 4,
        | Physical::Double => 5,
        | Physical::ByteArray => 6,
        | Physical::FixedLenByteArray(_) => 7,
    }
}

/// Maps a Parquet type code (and optional fixed length) to a physical type.
fn physical_from_code(
    code: i32,
    type_length: Option<i32>,
) -> Result<Physical, String> {
    match code {
        | 0 => Ok(Physical::Boolean),
        | 1 => Ok(Physical::Int32),
        | 2 => Ok(Physical::Int64),
        | 3 => Ok(Physical::Int96),
        | 4 => Ok(Physical::Float),
        | 5 => Ok(Physical::Double),
        | 6 => Ok(Physical::ByteArray),
        | 7 => {
            let len = type_length
                .and_then(|v| usize::try_from(v).ok())
                .filter(|&v| v > 0)
                .ok_or_else(|| "FIXED_LEN_BYTE_ARRAY requires positive type_length".to_string())?;
            Ok(Physical::FixedLenByteArray(len))
        },
        | other => Err(format!("unsupported parquet physical type {other}")),
    }
}

/// Extracts leaf fields from the schema tree, flattening nested struct groups
/// into dotted column names and tracking maximum definition levels.
fn parse_schema(metadata: &FileMetaData) -> Result<Vec<LeafInfo>, String> {
    let root = metadata
        .schema
        .first()
        .ok_or_else(|| "parquet schema is empty".to_string())?;
    if root.type_.is_some() {
        return Err("parquet schema root must be a group".to_string());
    }
    let root_children = match root.num_children {
        | Some(n) => {
            usize::try_from(n).map_err(|_| "negative parquet root child count".to_string())?
        },
        | None => metadata.schema.len().saturating_sub(1),
    };
    let mut cursor = 1usize;
    let mut leaves = Vec::new();
    let mut prefix = Vec::new();
    for _ in 0..root_children {
        visit_schema_node(
            &metadata.schema,
            &mut cursor,
            &mut prefix,
            0,
            0,
            &[],
            false,
            &mut leaves,
        )?;
    }
    if cursor != metadata.schema.len() {
        return Err("parquet root child count does not match the fields".to_string());
    }
    if leaves.is_empty() {
        return Err("parquet schema has no columns".to_string());
    }
    Ok(leaves)
}

/// Recursively traverses one node of the Parquet schema tree.
///
/// Every repeated node adds one level of nesting. A `LIST` annotation is not
/// required: a bare repeated primitive and the two-level layout some writers
/// emit both decode through the same rule, and `MAP` becomes its `key` and
/// `value` leaves.
fn visit_schema_node(
    elements: &[SchemaElement],
    cursor: &mut usize,
    prefix: &mut Vec<String>,
    def_level: u32,
    rep_level: u32,
    lists: &[ListSpec],
    parent_optional: bool,
    leaves: &mut Vec<LeafInfo>,
) -> Result<(), String> {
    let element = elements
        .get(*cursor)
        .ok_or_else(|| "parquet schema ended before declared children".to_string())?;
    *cursor += 1;

    let repetition = element.repetition_type.unwrap_or(REPETITION_REQUIRED);
    prefix.push(element.name.clone());
    let mut lists = lists.to_vec();
    let (next_def, next_rep) = if repetition == REPETITION_REPEATED {
        // The node above a repeated node is the list's group, so its
        // optionality is what decides whether the list may be null.
        if let Some(last) = lists.last_mut() {
            last.element_optional = parent_optional;
        }
        lists.push(ListSpec {
            present: def_level,
            optional: parent_optional,
            element_optional: false,
        });
        (def_level + 1, rep_level + 1)
    } else {
        (
            def_level + u32::from(repetition == REPETITION_OPTIONAL),
            rep_level,
        )
    };

    if let Some(code) = element.type_ {
        let physical = physical_from_code(code, element.type_length)?;
        let path = prefix.clone();
        let leaf_name = path.last().cloned().unwrap_or_default();
        let leaf_optional = repetition == REPETITION_OPTIONAL;
        if let Some(last) = lists.last_mut() {
            last.element_optional = leaf_optional;
        }
        let owner = prefix.first().cloned().unwrap_or_default();
        let name = if lists.is_empty() {
            path.join(".")
        } else {
            // A list's leaf is called `element` (or `item`, `array`, `list`);
            // a map's leaves are `key` and `value`.
            let bare = matches!(leaf_name.as_str(), "element" | "item" | "array" | "list")
                || leaf_name == owner;
            if bare || lists.len() > 1 {
                owner
            } else {
                format!("{owner}.{leaf_name}")
            }
        };
        let optional = lists.first().map_or(next_def > 0, |spec| spec.optional);
        leaves.push(LeafInfo {
            field: Field {
                name,
                physical,
                optional,
            },
            path,
            max_def_level: next_def,
            max_rep_level: next_rep,
            lists,
            is_utf8: element.converted_type == Some(CONVERTED_UTF8),
            decimal_scale: if element.converted_type == Some(CONVERTED_DECIMAL) {
                element.scale.and_then(|scale| u32::try_from(scale).ok())
            } else {
                None
            },
        });
    } else {
        let children = element
            .num_children
            .and_then(|n| usize::try_from(n).ok())
            .ok_or_else(|| {
                format!(
                    "parquet group `{}` is missing a valid child count",
                    prefix.join(".")
                )
            })?;
        for _ in 0..children {
            visit_schema_node(
                elements,
                cursor,
                prefix,
                next_def,
                next_rep,
                &lists,
                repetition == REPETITION_OPTIONAL,
                leaves,
            )?;
        }
    }

    prefix.pop();
    Ok(())
}

/// Returns the minimum bit width needed to represent values in `0..=max_val`.
fn bits_needed(max_val: u32) -> u32 {
    if max_val == 0 {
        0
    } else {
        32 - max_val.leading_zeros()
    }
}

/// Writes one data page and returns its uncompressed and compressed sizes.
///
/// A version-1 page stores the levels inside the compressed payload; a
/// version-2 page stores them raw before it and compresses only the values.
#[allow(clippy::too_many_arguments, reason = "one call site per column writer")]
fn write_data_page(
    file: &mut Vec<u8>,
    options: &WriteOptions,
    field_name: &str,
    repetition: &[u8],
    definition: &[u8],
    values: Vec<u8>,
    encoding: i32,
    num_values: usize,
    num_nulls: usize,
    num_rows: usize,
) -> Result<(usize, usize), String> {
    let value_len = values.len();
    let (payload, uncompressed_len, compressed_len, header) = if options.data_page_version == 2 {
        let compressed_values = compress_payload(values, options.compression);
        let levels_len = repetition.len().saturating_add(definition.len());
        let compressed_total = levels_len.saturating_add(compressed_values.len());
        let mut payload = Vec::with_capacity(compressed_total);
        payload.extend_from_slice(repetition);
        payload.extend_from_slice(definition);
        payload.extend_from_slice(&compressed_values);
        let header = PageHeader {
            type_: PAGE_DATA_V2,
            uncompressed_page_size: i32::try_from(levels_len.saturating_add(value_len))
                .map_err(|_| format!("column `{field_name}` page is too large"))?,
            compressed_page_size: i32::try_from(compressed_total)
                .map_err(|_| format!("column `{field_name}` compressed page is too large"))?,
            data_page: None,
            data_page_v2: Some(DataPageHeaderV2 {
                num_values: i32::try_from(num_values)
                    .map_err(|_| format!("column `{field_name}` has too many values"))?,
                num_nulls: i32::try_from(num_nulls)
                    .map_err(|_| format!("column `{field_name}` has too many nulls"))?,
                num_rows: i32::try_from(num_rows)
                    .map_err(|_| format!("column `{field_name}` has too many rows"))?,
                encoding,
                definition_levels_byte_length: i32::try_from(definition.len())
                    .map_err(|_| "definition levels are too large".to_string())?,
                repetition_levels_byte_length: i32::try_from(repetition.len())
                    .map_err(|_| "repetition levels are too large".to_string())?,
                is_compressed: options.compression != Compression::Uncompressed,
            }),
            dictionary_page: None,
        };
        (
            payload,
            levels_len.saturating_add(value_len),
            compressed_total,
            header,
        )
    } else {
        let mut payload = Vec::with_capacity(repetition.len() + definition.len() + value_len);
        if !repetition.is_empty() {
            payload.extend_from_slice(
                &u32::try_from(repetition.len())
                    .map_err(|_| "repetition levels are too large".to_string())?
                    .to_le_bytes(),
            );
            payload.extend_from_slice(repetition);
        }
        if !definition.is_empty() {
            payload.extend_from_slice(
                &u32::try_from(definition.len())
                    .map_err(|_| "definition levels are too large".to_string())?
                    .to_le_bytes(),
            );
            payload.extend_from_slice(definition);
        }
        payload.extend_from_slice(&values);
        let uncompressed_len = payload.len();
        let compressed_payload = compress_payload(payload, options.compression);
        let compressed_len = compressed_payload.len();
        let header = PageHeader {
            type_: PAGE_DATA,
            uncompressed_page_size: i32::try_from(uncompressed_len)
                .map_err(|_| format!("column `{field_name}` page is too large"))?,
            compressed_page_size: i32::try_from(compressed_len)
                .map_err(|_| format!("column `{field_name}` compressed page is too large"))?,
            data_page: Some(DataPageHeader {
                num_values: i32::try_from(num_values)
                    .map_err(|_| format!("column `{field_name}` has too many values"))?,
                encoding,
                definition_level_encoding: ENCODING_RLE,
                repetition_level_encoding: ENCODING_RLE,
            }),
            data_page_v2: None,
            dictionary_page: None,
        };
        (compressed_payload, uncompressed_len, compressed_len, header)
    };
    let header_bytes = thrift::encode_page_header(&header);
    file.extend_from_slice(&header_bytes);
    file.extend_from_slice(&payload);
    Ok((
        header_bytes.len().saturating_add(uncompressed_len),
        header_bytes.len().saturating_add(compressed_len),
    ))
}

/// Returns a typed empty column for a leaf field.
fn empty_column(leaf: &LeafInfo) -> Column {
    let base = match leaf.field.physical {
        | Physical::Boolean => Column::Bool(Vec::new()),
        | Physical::Int32 | Physical::Int64 | Physical::Int96 => Column::Int(Vec::new()),
        | Physical::Float | Physical::Double => Column::Float(Vec::new()),
        | Physical::ByteArray if leaf.is_utf8 => Column::Str(Vec::new()),
        | Physical::ByteArray | Physical::FixedLenByteArray(_) => Column::Bytes(Vec::new()),
    };
    let mut column = base;
    // The innermost element decides the leaf's own nullability.
    let element_optional = leaf
        .lists
        .last()
        .map_or(leaf.field.optional, |spec| spec.element_optional);
    if element_optional {
        column = Column::Nullable(Box::new(column), Vec::new());
    }
    for spec in leaf.lists.iter().rev() {
        let list = Column::List {
            values: Box::new(column),
            offsets: vec![0],
        };
        column = if spec.optional {
            Column::Nullable(Box::new(list), Vec::new())
        } else {
            list
        };
    }
    column
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
        | (Column::Bytes(first), Column::Bytes(second)) => {
            let mut values = first.clone();
            values.extend_from_slice(second);
            Ok(Column::Bytes(values))
        },
        | (Column::Nullable(first, first_mask), Column::Nullable(second, second_mask)) => {
            let inner = concat_columns(first, second)?;
            let mut mask = first_mask.clone();
            mask.extend_from_slice(second_mask);
            Ok(Column::Nullable(Box::new(inner), mask))
        },
        | (
            Column::List {
                values: first,
                offsets: first_offsets,
            },
            Column::List {
                values: second,
                offsets: second_offsets,
            },
        ) => {
            let base = first.len();
            let values = concat_columns(first, second)?;
            let mut offsets = if first_offsets.is_empty() {
                vec![0]
            } else {
                first_offsets.clone()
            };
            offsets.extend(second_offsets.iter().skip(1).map(|offset| offset + base));
            Ok(Column::List {
                values: Box::new(values),
                offsets,
            })
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

/// Reads the non-missing byte slice at `row`, looking through null masks.
fn bytes_at(
    column: &Column,
    row: usize,
) -> Option<&[u8]> {
    match column {
        | Column::Bytes(values) => values.get(row).map(Vec::as_slice),
        | Column::Nullable(inner, mask) => {
            if mask.get(row).copied().unwrap_or(true) {
                None
            } else {
                bytes_at(inner, row)
            }
        },
        | _ => None,
    }
}

/// Extracts the physical values of a column over `range`, dropping missing slots.
fn extract_values(
    field: &Field,
    column: &Column,
    range: Range<usize>,
) -> Result<Decoded, String> {
    let rows: Vec<usize> = range.collect();
    extract_selected(field, column, &rows)
}

/// Encodes the values at `rows` of `column`, skipping the missing ones.
fn extract_selected(
    field: &Field,
    column: &Column,
    rows: &[usize],
) -> Result<Decoded, String> {
    let mut decoded = match field.physical {
        | Physical::Boolean => Decoded::Bool(Vec::new()),
        | Physical::Int32 => Decoded::I32(Vec::new()),
        | Physical::Int64 => Decoded::I64(Vec::new()),
        | Physical::Int96 => Decoded::Int96(Vec::new()),
        | Physical::Float => Decoded::F32(Vec::new()),
        | Physical::Double => Decoded::F64(Vec::new()),
        | Physical::ByteArray if column.dtype() == DType::Str => Decoded::Str(Vec::new()),
        | Physical::ByteArray | Physical::FixedLenByteArray(_) => Decoded::RawBytes(Vec::new()),
    };
    for &row in rows {
        if column.is_null(row) {
            continue;
        }
        match &mut decoded {
            | Decoded::Bool(values) => {
                values.push(
                    bool_at(column, row)
                        .ok_or_else(|| format!("column `{}` is not boolean", field.name))?,
                );
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
                );
            },
            | Decoded::Int96(values) => {
                values.push(
                    int_at(column, row)
                        .ok_or_else(|| format!("column `{}` is not integral", field.name))?,
                );
            },
            | Decoded::F32(values) => {
                values.push(
                    float_at(column, row)
                        .ok_or_else(|| format!("column `{}` is not numeric", field.name))?
                        as f32,
                );
            },
            | Decoded::F64(values) => {
                values.push(
                    float_at(column, row)
                        .ok_or_else(|| format!("column `{}` is not numeric", field.name))?,
                );
            },
            | Decoded::Str(values) => {
                values.push(
                    str_at(column, row)
                        .ok_or_else(|| format!("column `{}` is not a string", field.name))?
                        .to_string(),
                );
            },
            | Decoded::RawBytes(values) => {
                let slice = bytes_at(column, row)
                    .ok_or_else(|| format!("column `{}` is not bytes", field.name))?;
                if let Physical::FixedLenByteArray(width) = field.physical {
                    if slice.len() != width {
                        return Err(format!(
                            "column `{}` expects {width}-byte values, got {} bytes",
                            field.name,
                            slice.len()
                        ));
                    }
                }
                values.push(slice.to_vec());
            },
        }
    }
    Ok(decoded)
}

/// Encodes already-extracted values with `PLAIN`.
fn encode_values(
    values: &Decoded,
    physical: Physical,
) -> Result<Vec<u8>, String> {
    match values {
        | Decoded::Bool(items) => Ok(plain::encode_bool(items)),
        | Decoded::I32(items) => Ok(plain::encode_i32(items)),
        | Decoded::I64(items) => Ok(plain::encode_i64(items)),
        | Decoded::Int96(items) => plain::encode_int96(items),
        | Decoded::F32(items) => Ok(plain::encode_f32(items)),
        | Decoded::F64(items) => Ok(plain::encode_f64(items)),
        | Decoded::Str(items) => {
            let bytes: Vec<Vec<u8>> = items
                .iter()
                .map(|value| value.as_bytes().to_vec())
                .collect();
            Ok(plain::encode_bytes(&bytes))
        },
        | Decoded::RawBytes(items) => {
            match physical {
                | Physical::FixedLenByteArray(width) => plain::encode_fixed_bytes(items, width),
                | _ => Ok(plain::encode_bytes(items)),
            }
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
        | Decoded::Bool(_) | Decoded::Int96(_) => {},
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
        | Decoded::Str(items) => {
            if let (Some(min), Some(max)) = (items.iter().min(), items.iter().max()) {
                statistics.min_value = Some(min.as_bytes().to_vec());
                statistics.max_value = Some(max.as_bytes().to_vec());
            }
        },
        | Decoded::RawBytes(items) => {
            if let (Some(min), Some(max)) = (items.iter().min(), items.iter().max()) {
                statistics.min_value = Some(min.clone());
                statistics.max_value = Some(max.clone());
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

/// Writes one column chunk (splitting across multiple data pages when needed)
/// and returns its metadata.
/// Levels of one repeated node, derived from the column it stores.
///
/// The definition levels mirror the three-level layout the writer emits: the
/// `LIST` group contributes one level when it is optional, the repeated group
/// one more, and the element one more when it can be null.
fn list_specs(
    column: &Column,
    optional: bool,
) -> Vec<ListSpec> {
    let mut specs = Vec::new();
    let mut current = column;
    let mut present = u32::from(optional);
    while let Some(values) = current.list_values() {
        let element_optional = values.is_nullable();
        specs.push(ListSpec {
            present,
            optional: current.is_nullable(),
            element_optional,
        });
        present = present + 1 + u32::from(element_optional);
        current = values;
    }
    specs
}

/// Shreds one list column into repetition and definition levels.
///
/// `depth` is the list being walked, `rep` the repetition level of the entry
/// that opened this list, and `values` collects the flat indices of the
/// present leaf values in entry order.
fn shred_list(
    column: &Column,
    specs: &[ListSpec],
    depth: usize,
    rep: u32,
    rows: Range<usize>,
    entries: &mut Vec<(u32, u32)>,
    values: &mut Vec<usize>,
) -> Result<(), String> {
    let spec = specs
        .get(depth - 1)
        .ok_or_else(|| "list depth exceeds the schema".to_string())?;
    let depth_level = u32::try_from(depth).unwrap_or(u32::MAX);
    let inner = column
        .list_values()
        .ok_or_else(|| "list column has no element column".to_string())?;
    for row in rows {
        let bounds = if column.is_null(row) {
            None
        } else {
            column.list_bounds(row)
        };
        let Some((start, end)) = bounds else {
            entries.push((rep, spec.present.saturating_sub(1)));
            continue;
        };
        if start == end {
            entries.push((rep, spec.present));
            continue;
        }
        for index in start..end {
            let element_rep = if index == start {
                rep
            } else {
                depth_level
            };
            if depth == specs.len() {
                let value_level = leaf_value_level(specs);
                if inner.is_null(index) {
                    entries.push((element_rep, value_level.saturating_sub(1)));
                } else {
                    entries.push((element_rep, value_level));
                    values.push(index);
                }
            } else {
                shred_list(
                    inner,
                    specs,
                    depth + 1,
                    element_rep,
                    index..index + 1,
                    entries,
                    values,
                )?;
            }
        }
    }
    Ok(())
}

/// Definition level of a missing element of the innermost list.
fn leaf_value_level(specs: &[ListSpec]) -> u32 {
    specs.last().map_or(0, |spec| {
        spec.present + 1 + u32::from(spec.element_optional)
    })
}

/// Writes one list column chunk, at any nesting depth.
///
/// A row contributes one entry per element, plus a single entry when the list
/// is null or empty. Pages are cut on row boundaries, and the values of a page
/// are the present leaves it contains, in order.
fn write_list_column_chunk(
    field: &Field,
    column: &Column,
    range: Range<usize>,
    offset: i64,
    options: &WriteOptions,
    file: &mut Vec<u8>,
) -> Result<ColumnMetaData, String> {
    let specs = list_specs(column, field.optional);
    if specs.is_empty() {
        return Err(format!("column `{}` is not a list", field.name));
    }
    let mut inner = column;
    while let Some(values) = inner.list_values() {
        inner = values;
    }

    // Levels for every row of the chunk, and the flat index of each present
    // leaf value in entry order.
    let mut entries: Vec<(u32, u32)> = Vec::new();
    let mut value_indices: Vec<usize> = Vec::new();
    let mut row_entries: Vec<usize> = Vec::with_capacity(range.len());
    let value_level = leaf_value_level(&specs);
    for row in range {
        let before = entries.len();
        let bounds = if column.is_null(row) {
            None
        } else {
            column.list_bounds(row)
        };
        if let Some((start, end)) = bounds {
            if start < end {
                shred_list(
                    column,
                    &specs,
                    1,
                    0,
                    row..row + 1,
                    &mut entries,
                    &mut value_indices,
                )?;
            } else {
                entries.push((0, specs[0].present));
            }
        } else {
            entries.push((0, specs[0].present.saturating_sub(1)));
        }
        row_entries.push(entries.len() - before);
    }
    let total = entries.len();
    let null_count = entries
        .iter()
        .filter(|(_, definition)| *definition < value_level)
        .count();
    let chunk_values = extract_selected(field, inner, &value_indices)?;

    let codec = codec_code(options.compression);

    let (dictionary, chunk_indices) =
        if options.use_dictionary && field.physical != Physical::Boolean {
            match dictionary_encode(&chunk_values) {
                | Some((unique, indices))
                    if unique.len() > 0 && unique.len() * 2 <= chunk_values.len() =>
                {
                    (Some(unique), Some(indices))
                },
                | _ => (None, None),
            }
        } else {
            (None, None)
        };

    let mut total_uncompressed = 0usize;
    let mut total_compressed = 0usize;
    let mut data_page_offset = offset;
    let mut dictionary_page_offset = None;

    if let Some(dict) = &dictionary {
        let payload = encode_values(dict, field.physical)?;
        let uncompressed_len = payload.len();
        let compressed_payload = compress_payload(payload, options.compression);
        let compressed_len = compressed_payload.len();
        let header = PageHeader {
            type_: PAGE_DICTIONARY,
            uncompressed_page_size: i32::try_from(uncompressed_len)
                .map_err(|_| format!("column `{}` dictionary page is too large", field.name))?,
            compressed_page_size: i32::try_from(compressed_len).map_err(|_| {
                format!(
                    "column `{}` compressed dictionary page is too large",
                    field.name
                )
            })?,
            data_page: None,
            data_page_v2: None,
            dictionary_page: Some(DictionaryPageHeader {
                num_values: i32::try_from(dict.len())
                    .map_err(|_| format!("column `{}` dictionary is too large", field.name))?,
                encoding: ENCODING_PLAIN,
            }),
        };
        let header_bytes = thrift::encode_page_header(&header);
        total_uncompressed = total_uncompressed
            .checked_add(header_bytes.len())
            .and_then(|value| value.checked_add(uncompressed_len))
            .ok_or_else(|| "column chunk size overflows".to_string())?;
        total_compressed = total_compressed
            .checked_add(header_bytes.len())
            .and_then(|value| value.checked_add(compressed_len))
            .ok_or_else(|| "column chunk size overflows".to_string())?;
        dictionary_page_offset = Some(offset);
        file.extend_from_slice(&header_bytes);
        file.extend_from_slice(&compressed_payload);
        let written = header_bytes.len().saturating_add(compressed_len);
        data_page_offset = offset.saturating_add(i64::try_from(written).unwrap_or(i64::MAX));
    }

    let mut entry_cursor = 0usize;
    let mut value_cursor = 0usize;
    let mut row_cursor = 0usize;
    while row_cursor < row_entries.len() {
        let mut page_rows = 0usize;
        let mut page_entries = 0usize;
        while row_cursor + page_rows < row_entries.len() && page_rows < options.max_page_rows {
            page_entries += row_entries[row_cursor + page_rows];
            page_rows += 1;
        }
        let entry_end = entry_cursor.saturating_add(page_entries).min(total);
        let slice = entries.get(entry_cursor..entry_end).unwrap_or(&[]);
        let repetitions: Vec<u32> = slice.iter().map(|(rep, _)| *rep).collect();
        let definitions: Vec<u32> = slice.iter().map(|(_, def)| *def).collect();
        let values_in_page = definitions
            .iter()
            .filter(|&&level| level == value_level)
            .count();

        let encoded_reps = plain::encode_level_values(&repetitions);
        let encoded_defs = plain::encode_level_values(&definitions);
        let mut value_bytes = Vec::new();
        let encoding = match (&dictionary, &chunk_indices) {
            | (Some(dict), Some(indices)) => {
                let bit_width =
                    bits_needed(u32::try_from(dict.len().saturating_sub(1)).unwrap_or(u32::MAX));
                value_bytes.push(bit_width as u8);
                let end = value_cursor
                    .saturating_add(values_in_page)
                    .min(indices.len());
                value_bytes.extend_from_slice(&plain::encode_indices(
                    indices.get(value_cursor..end).unwrap_or(&[]),
                    bit_width,
                ));
                ENCODING_RLE_DICTIONARY
            },
            | _ => {
                let end = value_cursor
                    .saturating_add(values_in_page)
                    .min(value_indices.len());
                let page_values = extract_selected(
                    field,
                    inner,
                    value_indices.get(value_cursor..end).unwrap_or(&[]),
                )?;
                value_bytes.extend_from_slice(&encode_values(&page_values, field.physical)?);
                ENCODING_PLAIN
            },
        };
        value_cursor = value_cursor.saturating_add(values_in_page);

        let (uncompressed, compressed) = write_data_page(
            file,
            options,
            &field.name,
            &encoded_reps,
            &encoded_defs,
            value_bytes,
            encoding,
            slice.len(),
            definitions
                .iter()
                .filter(|&&level| level < value_level)
                .count(),
            page_rows,
        )?;
        total_uncompressed = total_uncompressed
            .checked_add(uncompressed)
            .ok_or_else(|| "column chunk size overflows".to_string())?;
        total_compressed = total_compressed
            .checked_add(compressed)
            .ok_or_else(|| "column chunk size overflows".to_string())?;
        entry_cursor = entry_end;
        row_cursor += page_rows.max(1);
    }

    let encodings = if dictionary.is_some() {
        vec![ENCODING_RLE_DICTIONARY, ENCODING_RLE]
    } else {
        vec![ENCODING_PLAIN, ENCODING_RLE]
    };
    Ok(ColumnMetaData {
        type_: physical_code(field.physical),
        encodings,
        path_in_schema: list_path(field, &specs),
        codec,
        num_values: i64::try_from(total).map_err(|_| "value count overflows".to_string())?,
        total_uncompressed_size: i64::try_from(total_uncompressed)
            .map_err(|_| "column chunk size overflows".to_string())?,
        total_compressed_size: i64::try_from(total_compressed)
            .map_err(|_| "column chunk size overflows".to_string())?,
        data_page_offset,
        dictionary_page_offset,
        statistics: Some(build_statistics(
            &chunk_values,
            i64::try_from(null_count).map_err(|_| "null count overflows".to_string())?,
        )),
        bloom_filter_offset: None,
        bloom_filter_length: None,
    })
}

/// Schema path of a list leaf: the field, then `list` and `element` per level.
fn list_path(
    field: &Field,
    specs: &[ListSpec],
) -> Vec<String> {
    let mut path = vec![field.name.clone()];
    for _ in 0..specs.len() {
        path.push("list".to_string());
        path.push("element".to_string());
    }
    path
}

fn build_bloom_filter_for_decoded(values: &Decoded) -> SplitBlockBloomFilter {
    let mut bf = SplitBlockBloomFilter::with_ndv_fpp(values.len().max(1), 0.01);
    match values {
        | Decoded::Bool(v) => {
            for &b in v {
                bf.insert_bytes(&[u8::from(b)]);
            }
        },
        | Decoded::I32(v) => {
            for &x in v {
                bf.insert_i32(x);
            }
        },
        | Decoded::I64(v) | Decoded::Int96(v) => {
            for &x in v {
                bf.insert_i64(x);
            }
        },
        | Decoded::F32(v) => {
            for &x in v {
                bf.insert_f32(x);
            }
        },
        | Decoded::F64(v) => {
            for &x in v {
                bf.insert_f64(x);
            }
        },
        | Decoded::Str(v) => {
            for x in v {
                bf.insert_str(x);
            }
        },
        | Decoded::RawBytes(v) => {
            for x in v {
                bf.insert_bytes(x);
            }
        },
    }
    bf
}

fn write_column_chunk_ext(
    field: &Field,
    column: &Column,
    range: Range<usize>,
    offset: i64,
    options: &WriteOptions,
    file: &mut Vec<u8>,
) -> Result<(ColumnMetaData, Option<ColumnIndex>, Option<OffsetIndex>), String> {
    if column.list_values().is_some() {
        let meta = write_list_column_chunk(field, column, range, offset, options, file)?;
        return Ok((meta, None, None));
    }
    let rows = range.len();
    let chunk_values = extract_values(field, column, range.clone())?;
    let null_count = range.clone().filter(|&row| column.is_null(row)).count();
    let codec = codec_code(options.compression);

    // A dictionary pays for itself only when values repeat: build one when the
    // distinct count is at most half the value count. Booleans cannot be
    // dictionary-encoded in Parquet.
    let (dictionary, chunk_indices) =
        if options.use_dictionary && field.physical != Physical::Boolean {
            match dictionary_encode(&chunk_values) {
                | Some((unique, indices))
                    if unique.len() > 0 && unique.len() * 2 <= chunk_values.len() =>
                {
                    (Some(unique), Some(indices))
                },
                | _ => (None, None),
            }
        } else {
            (None, None)
        };

    let page_ranges: Vec<Range<usize>> = if rows == 0 {
        std::iter::once(range.start..range.start).collect()
    } else {
        let mut pages = Vec::new();
        let mut start = range.start;
        while start < range.end {
            let end = start.saturating_add(options.max_page_rows).min(range.end);
            pages.push(start..end);
            start = end;
        }
        pages
    };

    let mut total_uncompressed = 0usize;
    let mut total_compressed = 0usize;
    let mut data_page_offset = offset;
    let mut dictionary_page_offset = None;

    if let Some(dict) = &dictionary {
        let payload = encode_values(dict, field.physical)?;
        let uncompressed_len = payload.len();
        let compressed_payload = compress_payload(payload, options.compression);
        let compressed_len = compressed_payload.len();
        let header = PageHeader {
            type_: PAGE_DICTIONARY,
            uncompressed_page_size: i32::try_from(uncompressed_len)
                .map_err(|_| format!("column `{}` dictionary page is too large", field.name))?,
            compressed_page_size: i32::try_from(compressed_len).map_err(|_| {
                format!(
                    "column `{}` compressed dictionary page is too large",
                    field.name
                )
            })?,
            data_page: None,
            data_page_v2: None,
            dictionary_page: Some(DictionaryPageHeader {
                num_values: i32::try_from(dict.len())
                    .map_err(|_| format!("column `{}` dictionary is too large", field.name))?,
                encoding: ENCODING_PLAIN,
            }),
        };
        let header_bytes = thrift::encode_page_header(&header);
        total_uncompressed = total_uncompressed
            .checked_add(header_bytes.len())
            .and_then(|v| v.checked_add(uncompressed_len))
            .ok_or_else(|| "column chunk size overflows".to_string())?;
        total_compressed = total_compressed
            .checked_add(header_bytes.len())
            .and_then(|v| v.checked_add(compressed_len))
            .ok_or_else(|| "column chunk size overflows".to_string())?;
        dictionary_page_offset = Some(offset);
        file.extend_from_slice(&header_bytes);
        file.extend_from_slice(&compressed_payload);
        let written = header_bytes.len().saturating_add(compressed_len);
        data_page_offset = offset.saturating_add(i64::try_from(written).unwrap_or(i64::MAX));
    }

    let mut page_locations = Vec::with_capacity(page_ranges.len());
    let mut col_idx = ColumnIndex::default();
    let mut null_counts_vec = Vec::with_capacity(page_ranges.len());

    let mut value_cursor = 0usize;
    let mut first_row_in_rg = 0_i64;
    for page_range in page_ranges {
        let page_rows = page_range.len();
        let values = extract_values(field, column, page_range.clone())?;
        let page_nulls = page_range
            .clone()
            .filter(|&row| column.is_null(row))
            .count();
        let page_stats = build_statistics(&values, page_nulls as i64);
        let is_all_null = page_rows > 0 && page_nulls == page_rows;
        col_idx.null_pages.push(is_all_null);
        col_idx
            .min_values
            .push(page_stats.min_value.unwrap_or_default());
        col_idx
            .max_values
            .push(page_stats.max_value.unwrap_or_default());
        null_counts_vec.push(page_nulls as i64);

        let definition = if field.optional {
            // Parquet spec: definition level 1 = present, 0 = null.
            let levels: Vec<bool> = page_range.map(|row| !column.is_null(row)).collect();
            plain::encode_levels(&levels)
        } else {
            Vec::new()
        };
        let mut value_bytes = Vec::new();
        let encoding = match (&dictionary, &chunk_indices) {
            | (Some(dict), Some(indices)) => {
                let bit_width =
                    bits_needed(u32::try_from(dict.len().saturating_sub(1)).unwrap_or(u32::MAX));
                value_bytes.push(bit_width as u8);
                let end = value_cursor.saturating_add(values.len()).min(indices.len());
                value_bytes.extend_from_slice(&plain::encode_indices(
                    indices.get(value_cursor..end).unwrap_or(&[]),
                    bit_width,
                ));
                value_cursor = end;
                ENCODING_RLE_DICTIONARY
            },
            | _ => {
                value_bytes.extend_from_slice(&encode_values(&values, field.physical)?);
                ENCODING_PLAIN
            },
        };
        let page_start = file.len() as i64;
        let (uncompressed, compressed) = write_data_page(
            file,
            options,
            &field.name,
            &[],
            &definition,
            value_bytes,
            encoding,
            page_rows,
            page_nulls,
            page_rows,
        )?;
        page_locations.push(PageLocation {
            offset: page_start,
            compressed_page_size: compressed as i32,
            first_row_index: first_row_in_rg,
        });
        first_row_in_rg += page_rows as i64;
        total_uncompressed = total_uncompressed
            .checked_add(uncompressed)
            .ok_or_else(|| "column chunk size overflows".to_string())?;
        total_compressed = total_compressed
            .checked_add(compressed)
            .ok_or_else(|| "column chunk size overflows".to_string())?;
    }

    let (bloom_filter_offset, bloom_filter_length) = if options.write_bloom_filter {
        let bf = build_bloom_filter_for_decoded(&chunk_values);
        let hdr = thrift::encode_bloom_filter_header(bf.num_bytes() as i32);
        let bf_bytes = bf.to_bytes();
        let bf_off = file.len() as i64;
        let bf_len = (hdr.len() + bf_bytes.len()) as i32;
        file.extend_from_slice(&hdr);
        file.extend_from_slice(&bf_bytes);
        (Some(bf_off), Some(bf_len))
    } else {
        (None, None)
    };

    let (col_idx_opt, off_idx_opt) = if options.write_page_index {
        col_idx.null_counts = Some(null_counts_vec);
        col_idx.boundary_order = 0;
        (
            Some(col_idx),
            Some(OffsetIndex {
                page_locations,
                unencoded_byte_array_data_bytes: None,
            }),
        )
    } else {
        (None, None)
    };

    let encodings = if dictionary.is_some() {
        vec![ENCODING_RLE_DICTIONARY, ENCODING_RLE]
    } else {
        vec![ENCODING_PLAIN, ENCODING_RLE]
    };

    let meta = ColumnMetaData {
        type_: physical_code(field.physical),
        encodings,
        path_in_schema: vec![field.name.clone()],
        codec,
        num_values: i64::try_from(rows).map_err(|_| "row count overflows".to_string())?,
        total_uncompressed_size: i64::try_from(total_uncompressed)
            .map_err(|_| "column chunk size overflows".to_string())?,
        total_compressed_size: i64::try_from(total_compressed)
            .map_err(|_| "column chunk size overflows".to_string())?,
        data_page_offset,
        dictionary_page_offset,
        statistics: Some(build_statistics(
            &chunk_values,
            i64::try_from(null_count).map_err(|_| "null count overflows".to_string())?,
        )),
        bloom_filter_offset,
        bloom_filter_length,
    };
    Ok((meta, col_idx_opt, off_idx_opt))
}

/// Decodes the definition levels of one page into a null mask.
///
/// `level_data` is the RLE/bit-packed stream alone: the version-1 layout
/// prefixes it with its length, while a version-2 page stores the length in
/// the page header and passes the bytes directly.
/// Decodes a run of RLE/bit-packed definition or repetition levels.
fn decode_levels(
    level_data: &[u8],
    max_level: u32,
    count: usize,
    label: &str,
    name: &str,
) -> Result<Vec<u32>, String> {
    let bit_width = bits_needed(max_level);
    let levels = plain::decode_rle_hybrid(level_data, bit_width, count)?;
    for level in &levels {
        if *level > max_level {
            return Err(format!(
                "column `{name}` {label} level {level} exceeds max {max_level}"
            ));
        }
    }
    Ok(levels)
}

/// Reads one version-1 level section: a four-byte length and the levels.
fn read_level_section<'a>(
    data: &'a [u8],
    at: usize,
    name: &str,
) -> Result<(&'a [u8], usize), String> {
    let length_bytes = data
        .get(at..at.saturating_add(4))
        .ok_or_else(|| format!("column `{name}` levels are truncated"))?;
    let mut raw = [0u8; 4];
    raw.copy_from_slice(length_bytes);
    let len = u32::from_le_bytes(raw) as usize;
    let end = at
        .checked_add(4)
        .and_then(|start| start.checked_add(len))
        .ok_or_else(|| "level section length overflows".to_string())?;
    let levels = data
        .get(at + 4..end)
        .ok_or_else(|| format!("column `{name}` levels exceed the page"))?;
    Ok((levels, end))
}

/// Compresses a page payload with the configured codec.
fn compress_payload(
    payload: Vec<u8>,
    compression: Compression,
) -> Vec<u8> {
    match compression {
        | Compression::Uncompressed => payload,
        | Compression::Snappy => snappy::compress(&payload),
        | Compression::Gzip => gzip::compress(&payload),
        | Compression::Lzo => lzo::compress_hadoop(&payload),
        | Compression::Brotli => {
            let mut out = Vec::new();
            let mut cursor = std::io::Cursor::new(&payload);
            let params = brotli::enc::BrotliEncoderParams {
                quality: 5,
                ..Default::default()
            };
            if brotli::BrotliCompress(&mut cursor, &mut out, &params).is_ok() {
                out
            } else {
                payload
            }
        },
        | Compression::Zstd => {
            match zstd::bulk::compress(&payload, 3) {
                | Ok(compressed) => compressed,
                | Err(_) => payload,
            }
        },
        | Compression::Lz4 => lz4::compress(&payload),
    }
}

/// Reuses the Zstandard context for independent pages in a single file read.
#[derive(Default)]
struct PageDecoder {
    zstd: Option<zstd::bulk::Decompressor<'static>>,
}

impl PageDecoder {
    fn decompress_zstd(
        &mut self,
        payload: &[u8],
        expected: usize,
    ) -> std::io::Result<Vec<u8>> {
        let decoder = match &mut self.zstd {
            | Some(decoder) => decoder,
            | empty @ None => empty.insert(zstd::bulk::Decompressor::new()?),
        };
        decoder.decompress(payload, expected)
    }
}

/// Decompresses a page payload according to `codec` and validates its length.
fn decompress_page(
    payload: &[u8],
    codec: i32,
    uncompressed_size: i32,
    field_name: &str,
    page_decoder: &mut PageDecoder,
) -> Result<Vec<u8>, String> {
    let expected = usize::try_from(uncompressed_size)
        .map_err(|_| format!("column `{field_name}` has a negative uncompressed page size"))?;
    let data = match codec {
        | CODEC_UNCOMPRESSED => payload.to_vec(),
        | CODEC_SNAPPY => snappy::decompress(payload)?,
        | CODEC_GZIP => gzip::decompress(payload)?,
        | CODEC_LZO => lzo::decompress_hadoop(payload, expected)?,
        | CODEC_BROTLI => {
            let mut out = Vec::with_capacity(expected);
            let mut cursor = std::io::Cursor::new(payload);
            brotli::BrotliDecompress(&mut cursor, &mut out)
                .map_err(|error| format!("column `{field_name}` brotli page: {error}"))?;
            out
        },
        | CODEC_ZSTD => {
            page_decoder
                .decompress_zstd(payload, expected)
                .map_err(|error| format!("column `{field_name}` zstd page: {error}"))?
        },
        | CODEC_LZ4_RAW => lz4::decompress(payload, expected)?,
        | CODEC_LZ4 => lz4::decompress_hadoop(payload, expected)?,
        | other => {
            return Err(format!(
                "column `{field_name}` uses unsupported compression codec {other}"
            ));
        },
    };
    if data.len() != expected {
        return Err(format!(
            "column `{field_name}` page decompressed to {} bytes, expected {expected}",
            data.len()
        ));
    }
    Ok(data)
}

/// Creates an empty [`Decoded`] buffer for a leaf field.
fn new_decoded(
    leaf: &LeafInfo,
    capacity: usize,
) -> Decoded {
    let cap = capacity.min(1_000_000);
    // A decimal is exposed as a scaled `f64` whatever its physical storage.
    if leaf.decimal_scale.is_some() {
        return Decoded::F64(Vec::with_capacity(cap));
    }
    match leaf.field.physical {
        | Physical::Boolean => Decoded::Bool(Vec::with_capacity(cap)),
        | Physical::Int32 => Decoded::I32(Vec::with_capacity(cap)),
        | Physical::Int64 => Decoded::I64(Vec::with_capacity(cap)),
        | Physical::Int96 => Decoded::Int96(Vec::with_capacity(cap)),
        | Physical::Float => Decoded::F32(Vec::with_capacity(cap)),
        | Physical::Double => Decoded::F64(Vec::with_capacity(cap)),
        | Physical::ByteArray if leaf.is_utf8 => Decoded::Str(Vec::with_capacity(cap)),
        | Physical::ByteArray | Physical::FixedLenByteArray(_) => {
            Decoded::RawBytes(Vec::with_capacity(cap))
        },
    }
}

/// Reads and decodes a dictionary page starting at `pos`.
fn read_dictionary_page(
    bytes: &[u8],
    pos: usize,
    limit: usize,
    codec: i32,
    leaf: &LeafInfo,
    page_decoder: &mut PageDecoder,
) -> Result<(Decoded, usize), String> {
    if pos >= limit || pos >= bytes.len() {
        return Err(format!(
            "column `{}` dictionary page offset is out of bounds",
            leaf.field.name
        ));
    }
    let slice = bytes
        .get(pos..)
        .ok_or_else(|| "dictionary page header is truncated".to_string())?;
    let (header, consumed) = thrift::decode_page_header(slice)?;
    if header.type_ != PAGE_DICTIONARY {
        return Err(format!(
            "column `{}` expected dictionary page, got page type {}",
            leaf.field.name, header.type_
        ));
    }
    let payload_start = pos
        .checked_add(consumed)
        .ok_or_else(|| "dictionary page header length overflows".to_string())?;
    let compressed_len = usize::try_from(header.compressed_page_size).map_err(|_| {
        format!(
            "column `{}` has a negative dictionary size",
            leaf.field.name
        )
    })?;
    let payload_end = payload_start
        .checked_add(compressed_len)
        .ok_or_else(|| "dictionary page length overflows".to_string())?;
    if payload_end > limit {
        return Err(format!(
            "column `{}` dictionary page extends past the data section",
            leaf.field.name
        ));
    }
    let payload = bytes
        .get(payload_start..payload_end)
        .ok_or_else(|| format!("column `{}` dictionary page is truncated", leaf.field.name))?;
    let data = decompress_page(
        payload,
        codec,
        header.uncompressed_page_size,
        &leaf.field.name,
        page_decoder,
    )?;
    let dict_hdr = header.dictionary_page.ok_or_else(|| {
        format!(
            "column `{}` dictionary page is missing its header",
            leaf.field.name
        )
    })?;
    if dict_hdr.encoding != ENCODING_PLAIN && dict_hdr.encoding != ENCODING_PLAIN_DICTIONARY {
        return Err(format!(
            "column `{}` uses unsupported dictionary encoding {}",
            leaf.field.name, dict_hdr.encoding
        ));
    }
    let dict_count = usize::try_from(dict_hdr.num_values).map_err(|_| {
        format!(
            "column `{}` has a negative dictionary value count",
            leaf.field.name
        )
    })?;
    let mut dict = new_decoded(leaf, dict_count);
    append_plain_values(&mut dict, leaf, &data, dict_count)?;
    Ok((dict, payload_end))
}

/// Reads one column chunk, appending all of its pages into a column.
fn read_column_chunk(
    bytes: &[u8],
    metadata: &ColumnMetaData,
    leaf: &LeafInfo,
    limit: usize,
    page_decoder: &mut PageDecoder,
) -> Result<Column, String> {
    let total = usize::try_from(metadata.num_values)
        .map_err(|_| format!("column `{}` has a negative value count", leaf.field.name))?;
    if total > MAX_ROWS {
        return Err(format!(
            "column `{}` declares too many values",
            leaf.field.name
        ));
    }
    if total == 0 {
        return Ok(empty_column(leaf));
    }

    let mut dictionary = if let Some(dict_offset) = metadata.dictionary_page_offset {
        if dict_offset < 4 {
            return Err(format!(
                "column `{}` has an invalid dictionary page offset",
                leaf.field.name
            ));
        }
        let dpos = usize::try_from(dict_offset).map_err(|_| {
            format!(
                "column `{}` has an invalid dictionary page offset",
                leaf.field.name
            )
        })?;
        let (dict, _) =
            read_dictionary_page(bytes, dpos, limit, metadata.codec, leaf, page_decoder)?;
        Some(dict)
    } else {
        None
    };

    if metadata.data_page_offset < 4 {
        return Err(format!(
            "column `{}` has an invalid data page offset",
            leaf.field.name
        ));
    }
    let mut pos = usize::try_from(metadata.data_page_offset)
        .map_err(|_| format!("column `{}` has an invalid page offset", leaf.field.name))?;

    let mut decoded = new_decoded(leaf, total);
    let mut nulls: Vec<bool> = Vec::with_capacity(total.min(1_000_000));
    let mut reps: Vec<u32> = Vec::new();
    let mut defs: Vec<u32> = Vec::new();
    let mut collected = 0usize;

    while collected < total {
        if pos >= limit || pos >= bytes.len() {
            return Err(format!(
                "column `{}` ended before all values were read",
                leaf.field.name
            ));
        }
        let slice = bytes
            .get(pos..)
            .ok_or_else(|| "page header is truncated".to_string())?;
        let (header, consumed) = thrift::decode_page_header(slice)?;
        let payload_start = pos
            .checked_add(consumed)
            .ok_or_else(|| "page header length overflows".to_string())?;
        let compressed_len = usize::try_from(header.compressed_page_size)
            .map_err(|_| format!("column `{}` has a negative page size", leaf.field.name))?;
        let payload_end = payload_start
            .checked_add(compressed_len)
            .ok_or_else(|| "page length overflows".to_string())?;
        if payload_end > limit {
            return Err(format!(
                "column `{}` page extends past the data section",
                leaf.field.name
            ));
        }
        let payload = bytes
            .get(payload_start..payload_end)
            .ok_or_else(|| format!("column `{}` page is truncated", leaf.field.name))?;
        if header.type_ == PAGE_DICTIONARY {
            let dict_hdr = header.dictionary_page.ok_or_else(|| {
                format!(
                    "column `{}` dictionary page is missing its header",
                    leaf.field.name
                )
            })?;
            if dict_hdr.encoding != ENCODING_PLAIN && dict_hdr.encoding != ENCODING_PLAIN_DICTIONARY
            {
                return Err(format!(
                    "column `{}` uses unsupported dictionary encoding {}",
                    leaf.field.name, dict_hdr.encoding
                ));
            }
            let dict_count = usize::try_from(dict_hdr.num_values).map_err(|_| {
                format!(
                    "column `{}` has a negative dictionary value count",
                    leaf.field.name
                )
            })?;
            let data = decompress_page(
                payload,
                metadata.codec,
                header.uncompressed_page_size,
                &leaf.field.name,
                page_decoder,
            )?;
            let mut dict = new_decoded(leaf, dict_count);
            append_plain_values(&mut dict, leaf, &data, dict_count)?;
            dictionary = Some(dict);
            pos = payload_end;
            continue;
        }

        if header.type_ != PAGE_DATA && header.type_ != PAGE_DATA_V2 {
            return Err(format!(
                "column `{}` uses unsupported page type {}",
                leaf.field.name, header.type_
            ));
        }
        let page_values;
        let page_nulls: Vec<bool>;
        let page_defs: Vec<u32>;
        let page_reps: Vec<u32>;
        let data;
        let mut offset = 0usize;
        let encoding;

        if let Some(page) = header.data_page {
            data = decompress_page(
                payload,
                metadata.codec,
                header.uncompressed_page_size,
                &leaf.field.name,
                page_decoder,
            )?;
            page_values = usize::try_from(page.num_values)
                .map_err(|_| format!("column `{}` has a negative page count", leaf.field.name))?;
            if page.definition_level_encoding != ENCODING_RLE
                || page.repetition_level_encoding != ENCODING_RLE
            {
                return Err(format!(
                    "column `{}` uses unsupported level encoding",
                    leaf.field.name
                ));
            }
            let reps = if leaf.max_rep_level > 0 {
                let (level_data, next) = read_level_section(&data, 0, &leaf.field.name)?;
                offset = next;
                decode_levels(
                    level_data,
                    leaf.max_rep_level,
                    page_values,
                    "repetition",
                    &leaf.field.name,
                )?
            } else {
                vec![0u32; page_values]
            };
            let defs = if leaf.max_def_level > 0 {
                let (level_data, next) = read_level_section(&data, offset, &leaf.field.name)?;
                offset = next;
                decode_levels(
                    level_data,
                    leaf.max_def_level,
                    page_values,
                    "definition",
                    &leaf.field.name,
                )?
            } else {
                vec![0u32; page_values]
            };
            page_nulls = defs
                .iter()
                .map(|&level| level < leaf.max_def_level)
                .collect();
            page_defs = defs;
            page_reps = reps;
            encoding = page.encoding;
        } else if let Some(page) = header.data_page_v2 {
            page_values = usize::try_from(page.num_values)
                .map_err(|_| format!("column `{}` has a negative page count", leaf.field.name))?;
            let definition_len =
                usize::try_from(page.definition_levels_byte_length).map_err(|_| {
                    format!(
                        "column `{}` has a negative definition level length",
                        leaf.field.name
                    )
                })?;
            let repetition_len =
                usize::try_from(page.repetition_levels_byte_length).map_err(|_| {
                    format!(
                        "column `{}` has a negative repetition level length",
                        leaf.field.name
                    )
                })?;
            let levels_end = repetition_len
                .checked_add(definition_len)
                .ok_or_else(|| "page level lengths overflow".to_string())?;
            let repetition_data = payload.get(0..repetition_len).ok_or_else(|| {
                format!(
                    "column `{}` repetition levels exceed the page",
                    leaf.field.name
                )
            })?;
            let level_data = payload.get(repetition_len..levels_end).ok_or_else(|| {
                format!(
                    "column `{}` definition levels exceed the page",
                    leaf.field.name
                )
            })?;
            page_reps = if leaf.max_rep_level > 0 {
                decode_levels(
                    repetition_data,
                    leaf.max_rep_level,
                    page_values,
                    "repetition",
                    &leaf.field.name,
                )?
            } else {
                vec![0u32; page_values]
            };
            page_defs = if leaf.max_def_level > 0 {
                decode_levels(
                    level_data,
                    leaf.max_def_level,
                    page_values,
                    "definition",
                    &leaf.field.name,
                )?
            } else {
                vec![0u32; page_values]
            };
            page_nulls = page_defs
                .iter()
                .map(|&level| level < leaf.max_def_level)
                .collect();
            let value_bytes = payload
                .get(levels_end..)
                .ok_or_else(|| format!("column `{}` page payload is truncated", leaf.field.name))?;
            data = if page.is_compressed {
                // Only the value section is compressed; the level bytes are
                // stored raw and are not part of this length.
                let expected = usize::try_from(header.uncompressed_page_size)
                    .map_err(|_| format!("column `{}` has a negative page size", leaf.field.name))?
                    .checked_sub(levels_end)
                    .ok_or_else(|| "page level lengths exceed the page size".to_string())?;
                let expected = i32::try_from(expected)
                    .map_err(|_| "page value section is too large".to_string())?;
                decompress_page(
                    value_bytes,
                    metadata.codec,
                    expected,
                    &leaf.field.name,
                    page_decoder,
                )?
            } else {
                value_bytes.to_vec()
            };
            encoding = page.encoding;
        } else {
            return Err(format!(
                "column `{}` data page has no page header",
                leaf.field.name
            ));
        }

        if page_values == 0 {
            return Err(format!(
                "column `{}` contains an empty data page",
                leaf.field.name
            ));
        }
        let next = collected
            .checked_add(page_values)
            .ok_or_else(|| "page value count overflows".to_string())?;
        if next > total {
            return Err(format!(
                "column `{}` page declares more values than the chunk",
                leaf.field.name
            ));
        }

        let present = page_nulls.iter().filter(|&&is_null| !is_null).count();
        let rest = data
            .get(offset..)
            .ok_or_else(|| format!("column `{}` page payload is truncated", leaf.field.name))?;

        match encoding {
            | ENCODING_PLAIN => append_plain_values(&mut decoded, leaf, rest, present)?,
            | ENCODING_PLAIN_DICTIONARY | ENCODING_RLE_DICTIONARY => {
                if present > 0 {
                    let dict = dictionary.as_ref().ok_or_else(|| {
                        format!(
                            "column `{}` uses dictionary encoding without a dictionary page",
                            leaf.field.name
                        )
                    })?;
                    let bit_width = *rest.first().ok_or_else(|| {
                        format!(
                            "column `{}` dictionary index stream is truncated",
                            leaf.field.name
                        )
                    })?;
                    let indices =
                        plain::decode_rle_hybrid(&rest[1..], u32::from(bit_width), present)?;
                    append_dictionary_values(&mut decoded, dict, &indices, &leaf.field.name)?;
                }
            },
            | ENCODING_DELTA_BINARY_PACKED => {
                append_delta_binary_packed(&mut decoded, leaf, rest, present)?;
            },
            | ENCODING_DELTA_LENGTH_BYTE_ARRAY => {
                append_delta_length_byte_array(&mut decoded, leaf, rest, present)?;
            },
            | ENCODING_DELTA_BYTE_ARRAY => {
                append_delta_byte_array(&mut decoded, leaf, rest, present)?;
            },
            | ENCODING_BYTE_STREAM_SPLIT => {
                append_byte_stream_split(&mut decoded, leaf, rest, present)?;
            },
            | ENCODING_RLE => {
                append_rle_boolean(&mut decoded, rest, present)?;
            },
            | other => {
                return Err(format!(
                    "column `{}` uses unsupported value encoding {other}",
                    leaf.field.name
                ));
            },
        }

        nulls.extend_from_slice(&page_nulls);
        if leaf.is_list() {
            reps.extend_from_slice(&page_reps);
            defs.extend_from_slice(&page_defs);
        }
        collected = next;
        pos = payload_end;
    }

    if nulls.len() != total {
        return Err(format!(
            "column `{}` has {} level entries but {total} values",
            leaf.field.name,
            nulls.len()
        ));
    }
    if leaf.is_list() {
        return build_list_column(decoded, leaf, &reps, &defs);
    }
    Ok(decoded_into_column(decoded, &leaf.field, nulls))
}

/// Rebuilds a list column from the repetition and definition levels.
///
/// One entry describes one leaf slot. A repetition level of `r` starts a new
/// element of the list at depth `r`, and a definition level below a list's own
/// level means that list is null. The recursion peels one level of nesting at
/// a time, so any depth works.
fn build_list_column(
    decoded: Decoded,
    leaf: &LeafInfo,
    reps: &[u32],
    defs: &[u32],
) -> Result<Column, String> {
    if reps.len() != defs.len() {
        return Err(format!(
            "column `{}` has {} repetition levels but {} definition levels",
            leaf.field.name,
            reps.len(),
            defs.len()
        ));
    }
    if reps.first().is_some_and(|&level| level != 0) {
        return Err(format!(
            "column `{}` starts a row in the middle of a list",
            leaf.field.name
        ));
    }
    let entries: Vec<(u32, u32)> = reps.iter().copied().zip(defs.iter().copied()).collect();
    if std::env::var_os("RSSN_PARQUET_DEBUG").is_some() {
        eprintln!(
            "list {}: specs {:?} max_def {} max_rep {} entries {:?}",
            leaf.field.name,
            leaf.lists,
            leaf.max_def_level,
            leaf.max_rep_level,
            entries.iter().take(12).collect::<Vec<_>>()
        );
    }
    let mut decoded = Some(decoded);
    build_list_level(&entries, 1, leaf, &mut decoded)
}

/// Builds the column for the list at `depth` over `entries`.
fn build_list_level(
    entries: &[(u32, u32)],
    depth: usize,
    leaf: &LeafInfo,
    decoded: &mut Option<Decoded>,
) -> Result<Column, String> {
    let spec = leaf
        .lists
        .get(depth - 1)
        .ok_or_else(|| format!("column `{}` has no list at depth {depth}", leaf.field.name))?
        .clone();
    let depth_level = u32::try_from(depth).unwrap_or(u32::MAX);
    let mut offsets = vec![0usize];
    let mut mask: Vec<bool> = Vec::new();
    let mut children: Vec<(u32, u32)> = Vec::new();
    // Element slots, which is what the offsets count; `children` holds the
    // entries of every element and is longer than the number of elements.
    let mut elements = 0usize;
    let mut index = 0usize;
    while index < entries.len() {
        // One instance of this list runs until the next entry that starts an
        // outer element.
        let start = index;
        let mut end = index + 1;
        while end < entries.len() && entries[end].0 >= depth_level {
            end += 1;
        }
        let definition = entries[start].1;
        if definition < spec.present {
            mask.push(true);
        } else {
            mask.push(false);
            if definition > spec.present {
                // Elements of this list instance.
                let mut element = start;
                while element < end {
                    let mut next = element + 1;
                    while next < end && entries[next].0 > depth_level {
                        next += 1;
                    }
                    children.extend_from_slice(&entries[element..next]);
                    elements += 1;
                    element = next;
                }
            }
        }
        offsets.push(elements);
        index = end;
    }

    let values = if depth == leaf.lists.len() {
        let element_nulls: Vec<bool> = children
            .iter()
            .map(|&(_, definition)| definition < leaf.max_def_level)
            .collect();
        let decoded = decoded
            .take()
            .ok_or_else(|| format!("column `{}` values were consumed twice", leaf.field.name))?;
        if spec.element_optional {
            let element_field = Field {
                optional: true,
                ..leaf.field.clone()
            };
            decoded_into_column(decoded, &element_field, element_nulls)
        } else {
            decoded_column(decoded)
        }
    } else {
        build_list_level(&children, depth + 1, leaf, decoded)?
    };

    if std::env::var_os("RSSN_PARQUET_DEBUG").is_some() {
        eprintln!(
            "  depth {depth}: offsets {:?} mask {} children {}",
            offsets.iter().take(8).collect::<Vec<_>>(),
            mask.len(),
            children.len()
        );
    }
    let list = Column::List {
        values: Box::new(values),
        offsets,
    };
    // Nullability follows the schema, not the data, so row groups of one file
    // always concatenate.
    if spec.optional {
        Ok(Column::Nullable(Box::new(list), mask))
    } else {
        Ok(list)
    }
}

/// Decodes `count` `PLAIN`-encoded values of the field's physical type and appends them.
/// Appends boolean values stored with the RLE encoding, which version-2 data
/// pages use for `BOOLEAN` columns.
///
/// The stream is a four-byte little-endian length followed by RLE/bit-packed
/// runs at bit width one.
fn append_rle_boolean(
    decoded: &mut Decoded,
    bytes: &[u8],
    count: usize,
) -> Result<(), String> {
    if count == 0 {
        return Ok(());
    }
    let Decoded::Bool(values) = decoded else {
        return Err("rle boolean encoding used for a non-boolean column".to_string());
    };
    let header = bytes
        .get(0..4)
        .ok_or_else(|| "rle boolean page is truncated".to_string())?;
    let mut raw = [0u8; 4];
    raw.copy_from_slice(header);
    let len = u32::from_le_bytes(raw) as usize;
    let end = 4usize
        .checked_add(len)
        .ok_or_else(|| "rle boolean length overflows".to_string())?;
    let data = bytes
        .get(4..end)
        .ok_or_else(|| "rle boolean run is truncated".to_string())?;
    for bit in plain::decode_rle_hybrid(data, 1, count)? {
        values.push(bit != 0);
    }
    Ok(())
}

/// Decodes a big-endian two's-complement integer of arbitrary width.
fn decimal_from_bytes(bytes: &[u8]) -> f64 {
    if bytes.is_empty() {
        return 0.0;
    }
    let negative = bytes[0] & 0x80 != 0;
    let mut magnitude: u128 = 0;
    let mut overflow = false;
    for &byte in bytes {
        let value = if negative { !byte } else { byte };
        if magnitude > (u128::MAX >> 8) {
            overflow = true;
            break;
        }
        magnitude = (magnitude << 8) | u128::from(value);
    }
    let magnitude = if overflow {
        // Wider than 128 bits: accumulate in floating point, which is exact up
        // to 53 bits and correctly rounded beyond that.
        let mut acc = 0.0f64;
        for &byte in bytes {
            let value = if negative { !byte } else { byte };
            acc = acc * 256.0 + f64::from(value);
        }
        acc
    } else {
        magnitude as f64
    };
    if negative {
        -(magnitude + 1.0)
    } else {
        magnitude
    }
}

/// Appends `DECIMAL` values, dividing the stored integer by `10^scale`.
fn append_decimal_values(
    values: &mut Vec<f64>,
    leaf: &LeafInfo,
    bytes: &[u8],
    count: usize,
    scale: u32,
) -> Result<(), String> {
    let divisor = 10f64.powi(i32::try_from(scale).unwrap_or(i32::MAX / 2));
    match leaf.field.physical {
        | Physical::Int32 => {
            for value in plain::decode_i32(bytes, count)? {
                values.push(f64::from(value) / divisor);
            }
        },
        | Physical::Int64 => {
            for value in plain::decode_i64(bytes, count)? {
                values.push(value as f64 / divisor);
            }
        },
        | Physical::FixedLenByteArray(width) => {
            for raw in plain::decode_fixed_bytes(bytes, count, width)? {
                values.push(decimal_from_bytes(&raw) / divisor);
            }
        },
        | Physical::ByteArray => {
            for raw in plain::decode_bytes(bytes, count)? {
                values.push(decimal_from_bytes(raw) / divisor);
            }
        },
        | other => {
            return Err(format!(
                "column `{}` stores a decimal as an unsupported physical type {other:?}",
                leaf.field.name
            ));
        },
    }
    Ok(())
}

fn append_plain_values(
    decoded: &mut Decoded,
    leaf: &LeafInfo,
    bytes: &[u8],
    count: usize,
) -> Result<(), String> {
    if count == 0 {
        return Ok(());
    }
    match decoded {
        | Decoded::Bool(values) => {
            if leaf.field.physical != Physical::Boolean {
                return Err("column type mismatch".to_string());
            }
            values.extend(plain::decode_bool(bytes, count)?);
        },
        | Decoded::I32(values) => {
            if leaf.field.physical != Physical::Int32 {
                return Err("column type mismatch".to_string());
            }
            values.extend(plain::decode_i32(bytes, count)?);
        },
        | Decoded::I64(values) => {
            if leaf.field.physical != Physical::Int64 {
                return Err("column type mismatch".to_string());
            }
            values.extend(plain::decode_i64(bytes, count)?);
        },
        | Decoded::Int96(values) => {
            if leaf.field.physical != Physical::Int96 {
                return Err("column type mismatch".to_string());
            }
            values.extend(plain::decode_int96(bytes, count)?);
        },
        | Decoded::F32(values) => {
            if leaf.field.physical != Physical::Float {
                return Err("column type mismatch".to_string());
            }
            values.extend(plain::decode_f32(bytes, count)?);
        },
        | Decoded::F64(values) => {
            if let Some(scale) = leaf.decimal_scale {
                append_decimal_values(values, leaf, bytes, count, scale)?;
            } else {
                if leaf.field.physical != Physical::Double {
                    return Err("column type mismatch".to_string());
                }
                values.extend(plain::decode_f64(bytes, count)?);
            }
        },
        | Decoded::Str(values) => {
            if leaf.field.physical != Physical::ByteArray {
                return Err("column type mismatch".to_string());
            }
            for value in plain::decode_bytes(bytes, count)? {
                let text = String::from_utf8(value.to_vec())
                    .map_err(|_| format!("column `{}` contains invalid utf-8", leaf.field.name))?;
                values.push(text);
            }
        },
        | Decoded::RawBytes(values) => {
            match leaf.field.physical {
                | Physical::ByteArray => {
                    for value in plain::decode_bytes(bytes, count)? {
                        values.push(value.to_vec());
                    }
                },
                | Physical::FixedLenByteArray(width) => {
                    values.extend(plain::decode_fixed_bytes(bytes, count, width)?);
                },
                | _ => return Err("column type mismatch".to_string()),
            }
        },
    }
    Ok(())
}

/// Resolves dictionary `indices` against `dict` and appends the values to `decoded`.
fn append_dictionary_values(
    decoded: &mut Decoded,
    dict: &Decoded,
    indices: &[u32],
    field_name: &str,
) -> Result<(), String> {
    match (decoded, dict) {
        | (Decoded::Bool(dst), Decoded::Bool(src)) => {
            for &idx in indices {
                let val = *src.get(idx as usize).ok_or_else(|| {
                    format!("column `{field_name}` dictionary index {idx} is out of bounds")
                })?;
                dst.push(val);
            }
        },
        | (Decoded::I32(dst), Decoded::I32(src)) => {
            for &idx in indices {
                let val = *src.get(idx as usize).ok_or_else(|| {
                    format!("column `{field_name}` dictionary index {idx} is out of bounds")
                })?;
                dst.push(val);
            }
        },
        | (Decoded::I64(dst), Decoded::I64(src)) => {
            for &idx in indices {
                let val = *src.get(idx as usize).ok_or_else(|| {
                    format!("column `{field_name}` dictionary index {idx} is out of bounds")
                })?;
                dst.push(val);
            }
        },
        | (Decoded::Int96(dst), Decoded::Int96(src)) => {
            for &idx in indices {
                let val = *src.get(idx as usize).ok_or_else(|| {
                    format!("column `{field_name}` dictionary index {idx} is out of bounds")
                })?;
                dst.push(val);
            }
        },
        | (Decoded::F32(dst), Decoded::F32(src)) => {
            for &idx in indices {
                let val = *src.get(idx as usize).ok_or_else(|| {
                    format!("column `{field_name}` dictionary index {idx} is out of bounds")
                })?;
                dst.push(val);
            }
        },
        | (Decoded::F64(dst), Decoded::F64(src)) => {
            for &idx in indices {
                let val = *src.get(idx as usize).ok_or_else(|| {
                    format!("column `{field_name}` dictionary index {idx} is out of bounds")
                })?;
                dst.push(val);
            }
        },
        | (Decoded::Str(dst), Decoded::Str(src)) => {
            for &idx in indices {
                let val = src.get(idx as usize).ok_or_else(|| {
                    format!("column `{field_name}` dictionary index {idx} is out of bounds")
                })?;
                dst.push(val.clone());
            }
        },
        | (Decoded::RawBytes(dst), Decoded::RawBytes(src)) => {
            for &idx in indices {
                let val = src.get(idx as usize).ok_or_else(|| {
                    format!("column `{field_name}` dictionary index {idx} is out of bounds")
                })?;
                dst.push(val.clone());
            }
        },
        | _ => return Err(format!("column `{field_name}` dictionary type mismatch")),
    }
    Ok(())
}

/// Decodes a `DELTA_BINARY_PACKED` page for `Int32` or `Int64`.
fn append_delta_binary_packed(
    decoded: &mut Decoded,
    leaf: &LeafInfo,
    bytes: &[u8],
    count: usize,
) -> Result<(), String> {
    if count == 0 {
        return Ok(());
    }
    let (vals, _) = plain::decode_delta_binary_packed(bytes, count)?;
    match decoded {
        | Decoded::I32(dst) => {
            dst.extend(vals.into_iter().map(|v| v as i32));
            Ok(())
        },
        | Decoded::I64(dst) => {
            dst.extend(vals);
            Ok(())
        },
        | _ => {
            Err(format!(
                "column `{}` uses DELTA_BINARY_PACKED on non-integer physical type {:?}",
                leaf.field.name, leaf.field.physical
            ))
        },
    }
}

/// Decodes a `DELTA_LENGTH_BYTE_ARRAY` page for `ByteArray`.
fn append_delta_length_byte_array(
    decoded: &mut Decoded,
    leaf: &LeafInfo,
    bytes: &[u8],
    count: usize,
) -> Result<(), String> {
    if count == 0 {
        return Ok(());
    }
    let items = plain::decode_delta_length_byte_array(bytes, count)?;
    match decoded {
        | Decoded::Str(dst) => {
            for raw in items {
                let text = String::from_utf8(raw)
                    .map_err(|_| format!("column `{}` contains invalid utf-8", leaf.field.name))?;
                dst.push(text);
            }
            Ok(())
        },
        | Decoded::RawBytes(dst) => {
            dst.extend(items);
            Ok(())
        },
        | _ => {
            Err(format!(
                "column `{}` uses DELTA_LENGTH_BYTE_ARRAY on non-byte-array type {:?}",
                leaf.field.name, leaf.field.physical
            ))
        },
    }
}

/// Decodes a `DELTA_BYTE_ARRAY` page for `ByteArray`.
fn append_delta_byte_array(
    decoded: &mut Decoded,
    leaf: &LeafInfo,
    bytes: &[u8],
    count: usize,
) -> Result<(), String> {
    if count == 0 {
        return Ok(());
    }
    let items = plain::decode_delta_byte_array(bytes, count)?;
    match decoded {
        | Decoded::Str(dst) => {
            for raw in items {
                let text = String::from_utf8(raw)
                    .map_err(|_| format!("column `{}` contains invalid utf-8", leaf.field.name))?;
                dst.push(text);
            }
            Ok(())
        },
        | Decoded::RawBytes(dst) => {
            dst.extend(items);
            Ok(())
        },
        | _ => {
            Err(format!(
                "column `{}` uses DELTA_BYTE_ARRAY on non-byte-array type {:?}",
                leaf.field.name, leaf.field.physical
            ))
        },
    }
}

/// Decodes a `BYTE_STREAM_SPLIT` page for fixed-width numeric or byte types.
fn append_byte_stream_split(
    decoded: &mut Decoded,
    leaf: &LeafInfo,
    bytes: &[u8],
    count: usize,
) -> Result<(), String> {
    if count == 0 {
        return Ok(());
    }
    let width = match leaf.field.physical {
        | Physical::Float | Physical::Int32 => 4,
        | Physical::Double | Physical::Int64 => 8,
        | Physical::FixedLenByteArray(w) => w,
        | other => {
            return Err(format!(
                "column `{}` uses BYTE_STREAM_SPLIT on unsupported type {other:?}",
                leaf.field.name
            ));
        },
    };
    let plain_bytes = plain::decode_byte_stream_split(bytes, count, width)?;
    append_plain_values(decoded, leaf, &plain_bytes, count)
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
        | Decoded::I64(values) | Decoded::Int96(values) => Column::Int(expand(&values, &mask)),
        | Decoded::F32(values) => {
            Column::Float(expand(&values, &mask).into_iter().map(f64::from).collect())
        },
        | Decoded::F64(values) => Column::Float(expand(&values, &mask)),
        | Decoded::Str(values) => Column::Str(expand(&values, &mask)),
        | Decoded::RawBytes(values) => Column::Bytes(expand(&values, &mask)),
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
        | Decoded::I64(values) | Decoded::Int96(values) => Column::Int(values),
        | Decoded::F32(values) => Column::Float(values.into_iter().map(f64::from).collect()),
        | Decoded::F64(values) => Column::Float(values),
        | Decoded::Str(values) => Column::Str(values),
        | Decoded::RawBytes(values) => Column::Bytes(values),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use thrift::DictionaryPageHeader;

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

    /// Builds a 1000-row frame covering all basic physical types.
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

    /// Asserts two column lists match slot-by-slot (ignoring placeholders under null masks).
    fn assert_columns_match(
        expected: &[Column],
        actual: &[Column],
    ) {
        assert_eq!(actual.len(), expected.len());
        for (index, (exp, act)) in expected.iter().zip(actual).enumerate() {
            assert_eq!(act.len(), exp.len(), "column {index} length");
            assert_eq!(act.dtype(), exp.dtype(), "column {index} type");
            assert_eq!(
                act.is_nullable(),
                exp.is_nullable(),
                "column {index} nullability"
            );
            for row in 0..exp.len() {
                assert_eq!(
                    act.is_null(row),
                    exp.is_null(row),
                    "column {index} row {row} null"
                );
                if !exp.is_null(row) {
                    assert_eq!(
                        act.get_f64(row),
                        exp.get_f64(row),
                        "column {index} row {row} f64"
                    );
                    assert_eq!(
                        str_at(act, row),
                        str_at(exp, row),
                        "column {index} row {row} str"
                    );
                    assert_eq!(
                        bytes_at(act, row),
                        bytes_at(exp, row),
                        "column {index} row {row} bytes"
                    );
                }
            }
        }
    }

    /// Wraps raw page bytes and schema elements into a minimal valid Parquet file buffer.
    fn build_synthetic_parquet(
        schema: Vec<SchemaElement>,
        num_rows: i64,
        chunk_meta: ColumnMetaData,
        page_bytes: &[u8],
    ) -> Vec<u8> {
        let mut file = Vec::new();
        file.extend_from_slice(MAGIC);
        let offset = file.len() as i64;
        file.extend_from_slice(page_bytes);
        let rg = RowGroup {
            columns: vec![ColumnChunk {
                file_offset: offset,
                meta_data: Some(chunk_meta),
                ..Default::default()
            }],
            total_byte_size: page_bytes.len() as i64,
            num_rows,
            file_offset: Some(offset),
            total_compressed_size: Some(page_bytes.len() as i64),
        };
        let meta = FileMetaData {
            version: 1,
            schema,
            num_rows,
            row_groups: vec![rg],
        };
        let meta_bytes = thrift::encode_file_metadata(&meta);
        file.extend_from_slice(&meta_bytes);
        file.extend_from_slice(&(meta_bytes.len() as u32).to_le_bytes());
        file.extend_from_slice(MAGIC);
        file
    }

    #[test]
    fn round_trip_list_columns() {
        let scores = Column::new_list(
            Column::Nullable(
                Box::new(Column::Float(vec![1.5, 0.0, 2.5, 3.5])),
                vec![false, true, false, false],
            ),
            vec![0, 2, 2, 4],
        )
        .expect("scores");
        let scores = Column::Nullable(Box::new(scores), vec![false, true, false]);
        let tags = Column::new_list(
            Column::Str(vec!["a".to_string(), "b".to_string(), "c".to_string()]),
            vec![0, 2, 3, 3],
        )
        .expect("tags");
        let schema = Schema {
            fields: vec![
                Field {
                    name: "scores".to_string(),
                    physical: Physical::Double,
                    optional: true,
                },
                Field {
                    name: "tags".to_string(),
                    physical: Physical::ByteArray,
                    optional: false,
                },
            ],
        };
        let columns = vec![scores, tags];
        for options in [
            WriteOptions::default(),
            WriteOptions {
                compression: Compression::Snappy,
                max_row_group_rows: 2,
                max_page_rows: 1,
                use_dictionary: false,
                data_page_version: 2,
                ..Default::default()
            },
        ] {
            let bytes =
                write_parquet_bytes_with_options(&schema, &columns, &options).expect("write");
            let (read_schema, read_columns) = read_parquet_bytes(&bytes).expect("read");
            assert_eq!(read_schema, schema);
            assert_eq!(read_columns.len(), 2);
            assert_eq!(read_columns[0].len(), 3);
            assert_eq!(read_columns[0].dtype(), DType::List);
            assert!(read_columns[0].is_null(1));
            assert_eq!(
                read_columns[0].list_bounds(0),
                Some((0, 2)),
                "row 0 has two scores"
            );
            assert_eq!(read_columns[0].list_bounds(1), None, "row 1 is a null list");
            assert_eq!(read_columns[0].list_bounds(2), Some((2, 4)));
            let elements = read_columns[0].list_values().expect("elements");
            assert_eq!(elements.len(), 4);
            assert_eq!(elements.get_f64(0), Some(1.5));
            assert!(elements.is_null(1));
            assert_eq!(elements.get_f64(2), Some(2.5));
            assert_eq!(read_columns[1].list_bounds(0), Some((0, 2)));
            assert_eq!(read_columns[1].list_bounds(1), Some((2, 3)));
            assert_eq!(read_columns[1].list_bounds(2), Some((3, 3)));
        }
    }

    #[test]
    fn writes_version_two_data_pages() {
        let (schema, columns) = sample();
        for compression in [
            Compression::Uncompressed,
            Compression::Snappy,
            Compression::Gzip,
        ] {
            let options = WriteOptions {
                compression,
                data_page_version: 2,
                max_page_rows: 64,
                ..WriteOptions::default()
            };
            let bytes =
                write_parquet_bytes_with_options(&schema, &columns, &options).expect("write");
            let (read_schema, read_columns) = read_parquet_bytes(&bytes).expect("read");
            assert_eq!(read_schema, schema);
            assert_columns_match(&columns, &read_columns);
        }
    }

    #[test]
    fn round_trip_nested_list_columns() {
        // [[1, null], [4]], null, [], [[]]
        let inner = Column::new_list(
            Column::Nullable(
                Box::new(Column::Int(vec![1, 0, 4])),
                vec![false, true, false],
            ),
            vec![0, 2, 3, 3],
        )
        .expect("inner");
        let outer = Column::new_list(
            Column::Nullable(Box::new(inner), vec![false, false, false]),
            vec![0, 2, 2, 2, 3],
        )
        .expect("outer");
        let column = Column::Nullable(Box::new(outer), vec![false, true, false, false]);
        let schema = Schema {
            fields: vec![Field {
                name: "xs".to_string(),
                physical: Physical::Int32,
                optional: true,
            }],
        };
        for options in [
            WriteOptions::default(),
            WriteOptions {
                compression: Compression::Snappy,
                max_row_group_rows: 2,
                max_page_rows: 1,
                use_dictionary: false,
                data_page_version: 2,
                ..Default::default()
            },
        ] {
            let bytes =
                write_parquet_bytes_with_options(&schema, std::slice::from_ref(&column), &options)
                    .expect("write");
            let (read_schema, read_columns) = read_parquet_bytes(&bytes).expect("read");
            assert_eq!(read_schema, schema);
            assert_eq!(read_columns.len(), 1);
            let read = &read_columns[0];
            assert_eq!(read.len(), 4);
            assert_eq!(read.list_bounds(0), Some((0, 2)));
            assert_eq!(read.list_bounds(1), None);
            assert_eq!(read.list_bounds(2), Some((2, 2)));
            assert_eq!(read.list_bounds(3), Some((2, 3)));
            // Row 0: two inner lists.
            let inner_values = read.list_values().expect("inner column");
            assert_eq!(inner_values.dtype(), DType::List);
            assert_eq!(inner_values.list_bounds(0), Some((0, 2)));
            assert_eq!(inner_values.list_bounds(1), Some((2, 3)));
            let leaves = inner_values.list_values().expect("leaves");
            assert_eq!(leaves.get_f64(0), Some(1.0));
            assert!(leaves.is_null(1));
            assert_eq!(leaves.get_f64(2), Some(4.0));
            // Row 3 holds one empty inner list.
            assert_eq!(inner_values.list_bounds(2), Some((3, 3)));
        }
    }

    #[test]
    fn reads_arrow_written_nested_lists() {
        let base = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/parquet");
        for name in ["nested_lists.parquet", "nested_lists_v2.parquet"] {
            let (schema, columns) = read_parquet(&base.join(name)).expect("read nested list");
            assert_eq!(schema.fields[0].name, "xs");
            let read = &columns[0];
            assert_eq!(read.len(), 5, "{name}");
            assert_eq!(read.list_bounds(0), Some((0, 2)), "{name}");
            assert_eq!(read.list_bounds(1), None, "{name}");
            assert_eq!(read.list_bounds(2), Some((2, 2)), "{name}");
            assert_eq!(read.list_bounds(3), Some((2, 4)), "{name}");
            let inner = read.list_values().expect("inner");
            assert_eq!(inner.list_bounds(0), Some((0, 2)), "{name}");
            assert_eq!(inner.list_bounds(1), Some((2, 3)), "{name}");
            assert_eq!(inner.list_bounds(2), Some((3, 3)), "{name}");
            assert_eq!(inner.list_bounds(3), Some((3, 5)), "{name}");
            let leaves = inner.list_values().expect("leaves");
            assert_eq!(leaves.get_f64(0), Some(1.0), "{name}");
            assert_eq!(leaves.get_f64(1), Some(2.0), "{name}");
            assert_eq!(leaves.get_f64(2), Some(3.0), "{name}");
            assert!(leaves.is_null(4), "{name}");
            assert_eq!(leaves.get_f64(5), Some(5.0), "{name}");
        }
    }

    #[test]
    fn reads_arrow_written_value_encodings() {
        let base = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/parquet");
        // DELTA_BYTE_ARRAY: prefix lengths plus suffixes.
        let (_, columns) =
            read_parquet(&base.join("encodings_delta.parquet")).expect("read delta byte array");
        let base_column = match &columns[0] {
            | Column::Nullable(inner, _) => inner.as_ref(),
            | other => other,
        };
        let Column::Str(words) = base_column else {
            panic!("expected a string column");
        };
        assert_eq!(
            words,
            &vec![
                "alpha".to_string(),
                "alphabet".to_string(),
                "alphanumeric".to_string(),
                "beta".to_string(),
                "betamax".to_string(),
                "gamma".to_string(),
                "gammadion".to_string(),
                "delta".to_string(),
            ]
        );
        // BYTE_STREAM_SPLIT: one stream per byte of the value.
        let (_, columns) =
            read_parquet(&base.join("encodings_bss.parquet")).expect("read byte stream split");
        let base_column = match &columns[0] {
            | Column::Nullable(inner, _) => inner.as_ref(),
            | other => other,
        };
        let Column::Float(values) = base_column else {
            panic!("expected a float column");
        };
        for (index, value) in values.iter().enumerate() {
            assert!(
                (value - index as f64 * 0.5).abs() < 1e-12,
                "value {index} is {value}"
            );
        }
    }

    /// The fixture files under `tests/fixtures/parquet` are checked in, so
    /// this runs in a normal `cargo test`.
    #[test]
    fn reads_arrow_written_list_and_map_columns() {
        let base = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/parquet");
        for name in ["lists_v1.parquet", "lists_v2.parquet"] {
            let path = base.join(name);
            let (schema, columns) = read_parquet(&path).expect("read fixture");
            let names: Vec<&str> = schema.fields.iter().map(|f| f.name.as_str()).collect();
            assert_eq!(names, vec!["scores", "tags", "m.key", "m.value"], "{name}");
            assert_eq!(columns[0].len(), 4);
            assert_eq!(columns[0].list_bounds(0), Some((0, 2)));
            assert_eq!(columns[0].list_bounds(1), None);
            assert_eq!(columns[0].list_bounds(2), Some((2, 2)));
            assert_eq!(columns[0].list_bounds(3), Some((2, 4)));
            assert_eq!(columns[1].list_bounds(3), None);
            assert_eq!(columns[2].list_bounds(0), Some((0, 2)));
            assert_eq!(columns[3].list_bounds(3), Some((3, 3)));
        }
    }

    #[test]
    fn round_trip_all_physical_types() {
        let (schema, columns) = sample();
        let path = temp_path("round_trip");
        write_parquet(&path, &schema, &columns).expect("write");
        let (read_schema, read_columns) = read_parquet(&path).expect("read");
        let _ = std::fs::remove_file(&path);
        assert_eq!(read_schema, schema);
        assert_columns_match(&columns, &read_columns);
        assert_eq!(read_columns.first().map(Column::len), Some(1000));
        assert!(read_columns.first().is_some_and(|column| column.is_null(0)));
        assert!(
            read_columns
                .first()
                .is_some_and(|column| !column.is_null(1))
        );
    }

    #[test]
    fn snappy_compression_round_trip() {
        let (schema, columns) = sample();
        let uncompressed = write_parquet_bytes(&schema, &columns).expect("uncompressed");
        let opts = WriteOptions {
            compression: Compression::Snappy,
            ..WriteOptions::default()
        };
        let compressed =
            write_parquet_bytes_with_options(&schema, &columns, &opts).expect("snappy");
        assert!(compressed.len() < uncompressed.len());
        let (read_schema, read_columns) = read_parquet_bytes(&compressed).expect("read snappy");
        assert_eq!(read_schema, schema);
        assert_columns_match(&columns, &read_columns);
    }

    #[test]
    fn gzip_compression_round_trip() {
        let (schema, columns) = sample();
        let uncompressed = write_parquet_bytes(&schema, &columns).expect("uncompressed");
        let opts = WriteOptions {
            compression: Compression::Gzip,
            ..WriteOptions::default()
        };
        let compressed = write_parquet_bytes_with_options(&schema, &columns, &opts).expect("gzip");
        assert!(compressed.len() < uncompressed.len());
        let (read_schema, read_columns) = read_parquet_bytes(&compressed).expect("read gzip");
        assert_eq!(read_schema, schema);
        assert_columns_match(&columns, &read_columns);
    }

    #[test]
    fn zstd_and_lz4_compression_round_trip() {
        let (schema, columns) = sample();
        for compression in [Compression::Zstd, Compression::Lz4] {
            let opts = WriteOptions {
                compression,
                ..WriteOptions::default()
            };
            let bytes = write_parquet_bytes_with_options(&schema, &columns, &opts)
                .expect("write compressed");
            let (read_schema, read_columns) = read_parquet_bytes(&bytes).expect("read compressed");
            assert_eq!(read_schema, schema);
            assert_columns_match(&columns, &read_columns);
        }
    }

    #[test]
    fn zstd_reads_independent_pages_across_row_groups() {
        let schema = Schema {
            fields: vec![Field {
                name: "labels".into(),
                physical: Physical::ByteArray,
                optional: false,
            }],
        };
        let values = vec![Column::Str(
            (0..257).map(|i| format!("label-{}", i % 7)).collect(),
        )];
        for data_page_version in [1, 2] {
            let options = WriteOptions {
                compression: Compression::Zstd,
                max_page_rows: 17,
                max_row_group_rows: 91,
                data_page_version,
                ..Default::default()
            };
            let bytes = write_parquet_bytes_with_options(&schema, &values, &options).unwrap();
            let (actual_schema, actual_values) = read_parquet_bytes(&bytes).unwrap();
            assert_eq!(actual_schema, schema);
            assert_eq!(actual_values, values);
        }
    }

    /// Arrow writes these files; they are checked in so the test runs in a
    /// normal `cargo test`.
    #[test]
    fn reads_arrow_written_codec_fixtures() {
        let base = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/parquet");
        for name in ["codec_zstd.parquet", "codec_lz4.parquet"] {
            let path = base.join(name);
            let (_, columns) = read_parquet(&path).expect("read codec fixture");
            assert_eq!(columns.len(), 2, "{name} column count");
            let labels = match &columns[0] {
                | Column::Nullable(inner, _) => inner.as_ref(),
                | other => other,
            };
            let Column::Str(labels) = labels else {
                panic!("{name}: label column is not a string column");
            };
            assert_eq!(labels.len(), 2000);
            assert_eq!(labels[1], "kind-1");
            let values = match &columns[1] {
                | Column::Nullable(inner, _) => inner.as_ref(),
                | other => other,
            };
            let Column::Float(values) = values else {
                panic!("{name}: value column is not a float column");
            };
            assert!((values[3] - 1.5).abs() < 1e-12);
        }
    }

    #[test]
    fn reads_arrow_written_decimal_fixture() {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/parquet/decimal.parquet");
        let (_, columns) = read_parquet(&path).expect("read decimal fixture");
        let base = match &columns[0] {
            | Column::Nullable(inner, _) => inner.as_ref(),
            | other => other,
        };
        let Column::Float(values) = base else {
            panic!("decimal column is not a float column: {base:?}");
        };
        assert_eq!(values.len(), 4);
        assert!((values[0] - 123.45).abs() < 1e-9, "{}", values[0]);
        assert!((values[1] + 0.5).abs() < 1e-9, "{}", values[1]);
    }

    #[test]
    fn decimal_column_is_exposed_as_a_scaled_float() {
        // INT32 decimal with scale 2: 12345 means 123.45.
        let payload = plain::encode_i32(&[12345, -50, 0]);
        let header = thrift::encode_page_header(&PageHeader {
            type_: PAGE_DATA,
            uncompressed_page_size: payload.len() as i32,
            compressed_page_size: payload.len() as i32,
            data_page: Some(DataPageHeader {
                num_values: 3,
                encoding: ENCODING_PLAIN,
                definition_level_encoding: ENCODING_RLE,
                repetition_level_encoding: ENCODING_RLE,
            }),
            data_page_v2: None,
            dictionary_page: None,
        });
        let mut pages = header;
        pages.extend_from_slice(&payload);
        let schema = vec![
            SchemaElement {
                name: "schema".to_string(),
                num_children: Some(1),
                repetition_type: Some(REPETITION_REQUIRED),
                ..SchemaElement::default()
            },
            SchemaElement {
                type_: Some(1),
                repetition_type: Some(REPETITION_REQUIRED),
                name: "amount".to_string(),
                converted_type: Some(CONVERTED_DECIMAL),
                scale: Some(2),
                precision: Some(9),
                ..SchemaElement::default()
            },
        ];
        let chunk_meta = ColumnMetaData {
            type_: 1,
            encodings: vec![ENCODING_PLAIN],
            path_in_schema: vec!["amount".to_string()],
            codec: CODEC_UNCOMPRESSED,
            num_values: 3,
            total_uncompressed_size: pages.len() as i64,
            total_compressed_size: pages.len() as i64,
            data_page_offset: 4,
            dictionary_page_offset: None,
            statistics: None,
            ..Default::default()
        };
        let file = build_synthetic_parquet(schema, 3, chunk_meta, &pages);
        let (_, cols) = read_parquet_bytes(&file).expect("read decimal");
        let Column::Float(values) = &cols[0] else {
            panic!("decimal column is not a float column");
        };
        assert!((values[0] - 123.45).abs() < 1e-9, "{}", values[0]);
        assert!((values[1] + 0.5).abs() < 1e-9, "{}", values[1]);
        assert!(values[2].abs() < 1e-12);
    }

    #[test]
    fn dictionary_encoding_round_trip() {
        // Few distinct values, so the writer must choose a dictionary.
        let rows = 5_000usize;
        let labels: Vec<String> = (0..rows).map(|row| format!("kind-{}", row % 4)).collect();
        let codes: Vec<i64> = (0..rows).map(|row| (row % 7) as i64).collect();
        let schema = Schema {
            fields: vec![
                Field {
                    name: "label".to_string(),
                    physical: Physical::ByteArray,
                    optional: false,
                },
                Field {
                    name: "code".to_string(),
                    physical: Physical::Int64,
                    optional: true,
                },
            ],
        };
        let mask: Vec<bool> = (0..rows).map(|row| row % 11 == 0).collect();
        let columns = vec![
            Column::Str(labels),
            Column::Nullable(Box::new(Column::Int(codes)), mask),
        ];

        let plain_opts = WriteOptions {
            use_dictionary: false,
            ..WriteOptions::default()
        };
        let plain =
            write_parquet_bytes_with_options(&schema, &columns, &plain_opts).expect("plain write");
        let dict = write_parquet_bytes(&schema, &columns).expect("dictionary write");
        assert!(
            dict.len() < plain.len(),
            "dictionary {} is not smaller than plain {}",
            dict.len(),
            plain.len()
        );
        let (read_schema, read_columns) = read_parquet_bytes(&dict).expect("read dictionary");
        assert_eq!(read_schema, schema);
        assert_columns_match(&columns, &read_columns);

        // The same file with dictionary pages on top of gzip pages.
        let both = WriteOptions {
            compression: Compression::Gzip,
            ..WriteOptions::default()
        };
        let bytes = write_parquet_bytes_with_options(&schema, &columns, &both).expect("write");
        let (_, read_columns) = read_parquet_bytes(&bytes).expect("read dictionary+gzip");
        assert_columns_match(&columns, &read_columns);
    }

    #[test]
    fn single_value_dictionary_and_multi_page_round_trip() {
        // One distinct value exercises the zero-bit index stream; a small page
        // size forces the dictionary to be shared across data pages.
        let rows = 1_000usize;
        let values: Vec<String> = (0..rows).map(|_| "same".to_string()).collect();
        let schema = Schema {
            fields: vec![Field {
                name: "text".to_string(),
                physical: Physical::ByteArray,
                optional: false,
            }],
        };
        let columns = vec![Column::Str(values)];
        let opts = WriteOptions {
            max_page_rows: 37,
            ..WriteOptions::default()
        };
        let bytes = write_parquet_bytes_with_options(&schema, &columns, &opts).expect("write");
        let (_, read_columns) = read_parquet_bytes(&bytes).expect("read");
        assert_columns_match(&columns, &read_columns);
    }

    #[test]
    fn multi_row_group_and_multi_page_round_trip() {
        let (schema, columns) = sample();
        let opts = WriteOptions {
            compression: Compression::Snappy,
            max_row_group_rows: 230,
            max_page_rows: 75,
            ..WriteOptions::default()
        };
        let bytes = write_parquet_bytes_with_options(&schema, &columns, &opts).expect("write");
        let (read_schema, read_columns) = read_parquet_bytes(&bytes).expect("read");
        assert_eq!(read_schema, schema);
        assert_columns_match(&columns, &read_columns);
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
    fn raw_binary_and_fixed_len_byte_array_round_trip() {
        let schema = Schema {
            fields: vec![
                Field {
                    name: "raw".to_string(),
                    physical: Physical::ByteArray,
                    optional: true,
                },
                Field {
                    name: "fixed".to_string(),
                    physical: Physical::FixedLenByteArray(4),
                    optional: false,
                },
            ],
        };
        let columns = vec![
            Column::Nullable(
                Box::new(Column::Bytes(vec![
                    vec![0xff, 0xfe, 0x00, 0x80], // non-UTF-8
                    Vec::new(),
                    vec![0x01, 0x02],
                ])),
                vec![false, true, false],
            ),
            Column::Bytes(vec![
                vec![0xde, 0xad, 0xbe, 0xef],
                vec![0x00, 0x00, 0x00, 0x01],
                vec![0xca, 0xfe, 0xba, 0xbe],
            ]),
        ];
        let bytes = write_parquet_bytes(&schema, &columns).expect("write");
        let (read_schema, read_columns) = read_parquet_bytes(&bytes).expect("read");
        assert_eq!(read_schema, schema);
        assert_columns_match(&columns, &read_columns);
    }

    #[test]
    fn int96_round_trip() {
        let schema = Schema {
            fields: vec![Field {
                name: "ts".to_string(),
                physical: Physical::Int96,
                optional: true,
            }],
        };
        let columns = vec![Column::Nullable(
            Box::new(Column::Int(vec![
                0,
                0,
                1_700_000_000_123_456_789,
                -86_400_000_000_001,
            ])),
            vec![false, true, false, false],
        )];
        let bytes = write_parquet_bytes(&schema, &columns).expect("write");
        let (read_schema, read_columns) = read_parquet_bytes(&bytes).expect("read");
        assert_eq!(read_schema, schema);
        assert_eq!(read_columns, columns);
    }

    #[test]
    fn dictionary_encoded_page_decodes() {
        // Build a dictionary page with ["alpha", "beta", "gamma"] and a data page
        // referencing indices [2, 0, 1, 0, 2] using RLE_DICTIONARY.
        let dict_payload =
            plain::encode_bytes(&[b"alpha".to_vec(), b"beta".to_vec(), b"gamma".to_vec()]);
        let dict_header = thrift::encode_page_header(&PageHeader {
            type_: PAGE_DICTIONARY,
            uncompressed_page_size: dict_payload.len() as i32,
            compressed_page_size: dict_payload.len() as i32,
            data_page: None,
            data_page_v2: None,
            dictionary_page: Some(DictionaryPageHeader {
                num_values: 3,
                encoding: ENCODING_PLAIN,
            }),
        });
        let dict_offset = 4i64;
        let data_offset = dict_offset + (dict_header.len() + dict_payload.len()) as i64;

        // Bit width = 2, 1 bit-packed group of 8 values: [2, 0, 1, 0, 2, 0, 0, 0]
        // byte0 = 2 | (0<<2) | (1<<4) | (0<<6) = 0x12
        // byte1 = 2 | 0 = 0x02
        let data_payload = vec![2u8, (1u8 << 1) | 1, 0x12, 0x02];
        let data_header = thrift::encode_page_header(&PageHeader {
            type_: PAGE_DATA,
            uncompressed_page_size: data_payload.len() as i32,
            compressed_page_size: data_payload.len() as i32,
            data_page: Some(DataPageHeader {
                num_values: 5,
                encoding: ENCODING_RLE_DICTIONARY,
                definition_level_encoding: ENCODING_RLE,
                repetition_level_encoding: ENCODING_RLE,
            }),
            data_page_v2: None,
            dictionary_page: None,
        });

        let mut pages = Vec::new();
        pages.extend_from_slice(&dict_header);
        pages.extend_from_slice(&dict_payload);
        pages.extend_from_slice(&data_header);
        pages.extend_from_slice(&data_payload);

        let schema = vec![
            SchemaElement {
                name: "schema".to_string(),
                num_children: Some(1),
                repetition_type: Some(REPETITION_REQUIRED),
                ..SchemaElement::default()
            },
            SchemaElement {
                type_: Some(6),
                repetition_type: Some(REPETITION_REQUIRED),
                name: "word".to_string(),
                converted_type: Some(CONVERTED_UTF8),
                ..SchemaElement::default()
            },
        ];
        let chunk_meta = ColumnMetaData {
            type_: 6,
            encodings: vec![ENCODING_PLAIN_DICTIONARY, ENCODING_RLE_DICTIONARY],
            path_in_schema: vec!["word".to_string()],
            codec: CODEC_UNCOMPRESSED,
            num_values: 5,
            total_uncompressed_size: pages.len() as i64,
            total_compressed_size: pages.len() as i64,
            data_page_offset: data_offset,
            dictionary_page_offset: Some(dict_offset),
            statistics: None,
            ..Default::default()
        };
        let file = build_synthetic_parquet(schema, 5, chunk_meta, &pages);
        let (_, cols) = read_parquet_bytes(&file).expect("read dict");
        assert_eq!(
            cols,
            vec![Column::Str(vec![
                "gamma".to_string(),
                "alpha".to_string(),
                "beta".to_string(),
                "alpha".to_string(),
                "gamma".to_string(),
            ])]
        );
    }

    #[test]
    fn version_two_data_page_decodes() {
        // Required column, PLAIN values, no levels: the v2 layout puts the
        // level bytes before the (here uncompressed) value bytes.
        let payload = plain::encode_i32(&[7, -2, 42]);
        let header = thrift::encode_page_header(&PageHeader {
            type_: PAGE_DATA_V2,
            uncompressed_page_size: payload.len() as i32,
            compressed_page_size: payload.len() as i32,
            data_page: None,
            data_page_v2: Some(thrift::DataPageHeaderV2 {
                num_values: 3,
                num_nulls: 0,
                num_rows: 3,
                encoding: ENCODING_PLAIN,
                definition_levels_byte_length: 0,
                repetition_levels_byte_length: 0,
                is_compressed: false,
            }),
            dictionary_page: None,
        });
        let mut pages = header;
        pages.extend_from_slice(&payload);
        let schema = vec![
            SchemaElement {
                name: "schema".to_string(),
                num_children: Some(1),
                repetition_type: Some(REPETITION_REQUIRED),
                ..SchemaElement::default()
            },
            SchemaElement {
                type_: Some(1),
                repetition_type: Some(REPETITION_REQUIRED),
                name: "n".to_string(),
                ..SchemaElement::default()
            },
        ];
        let chunk_meta = ColumnMetaData {
            type_: 1,
            encodings: vec![ENCODING_PLAIN],
            path_in_schema: vec!["n".to_string()],
            codec: CODEC_UNCOMPRESSED,
            num_values: 3,
            total_uncompressed_size: pages.len() as i64,
            total_compressed_size: pages.len() as i64,
            data_page_offset: 4,
            dictionary_page_offset: None,
            statistics: None,
            ..Default::default()
        };
        let file = build_synthetic_parquet(schema, 3, chunk_meta, &pages);
        let (_, cols) = read_parquet_bytes(&file).expect("read v2");
        assert_eq!(cols, vec![Column::Int(vec![7, -2, 42])]);
    }

    #[test]
    fn version_two_optional_page_decodes_its_levels() {
        // Rows [1, null, 3, null]: the definition levels are stored raw (no
        // length prefix) and only the value section is compressed.
        let levels = plain::encode_levels(&[true, false, true, false]);
        let mut payload = levels.clone();
        payload.extend_from_slice(&plain::encode_i32(&[1, 3]));
        let header = thrift::encode_page_header(&PageHeader {
            type_: PAGE_DATA_V2,
            uncompressed_page_size: payload.len() as i32,
            compressed_page_size: payload.len() as i32,
            data_page: None,
            data_page_v2: Some(thrift::DataPageHeaderV2 {
                num_values: 4,
                num_nulls: 2,
                num_rows: 4,
                encoding: ENCODING_PLAIN,
                definition_levels_byte_length: levels.len() as i32,
                repetition_levels_byte_length: 0,
                is_compressed: false,
            }),
            dictionary_page: None,
        });
        let mut pages = header;
        pages.extend_from_slice(&payload);
        let schema = vec![
            SchemaElement {
                name: "schema".to_string(),
                num_children: Some(1),
                repetition_type: Some(REPETITION_REQUIRED),
                ..SchemaElement::default()
            },
            SchemaElement {
                type_: Some(1),
                repetition_type: Some(REPETITION_OPTIONAL),
                name: "n".to_string(),
                ..SchemaElement::default()
            },
        ];
        let chunk_meta = ColumnMetaData {
            type_: 1,
            encodings: vec![ENCODING_PLAIN],
            path_in_schema: vec!["n".to_string()],
            codec: CODEC_UNCOMPRESSED,
            num_values: 4,
            total_uncompressed_size: pages.len() as i64,
            total_compressed_size: pages.len() as i64,
            data_page_offset: 4,
            dictionary_page_offset: None,
            statistics: None,
            ..Default::default()
        };
        let file = build_synthetic_parquet(schema, 4, chunk_meta, &pages);
        let (_, cols) = read_parquet_bytes(&file).expect("read v2 optional");
        assert_eq!(
            cols,
            vec![Column::Nullable(
                Box::new(Column::Int(vec![1, 0, 3, 0])),
                vec![false, true, false, true]
            )]
        );
    }

    #[test]
    fn rle_boolean_values_decode() {
        // Version-2 pages store BOOLEAN columns as RLE with a four-byte length.
        let runs = plain::encode_levels(&[true, true, false, true]);
        let mut payload = Vec::new();
        payload.extend_from_slice(&(runs.len() as u32).to_le_bytes());
        payload.extend_from_slice(&runs);
        let header = thrift::encode_page_header(&PageHeader {
            type_: PAGE_DATA_V2,
            uncompressed_page_size: payload.len() as i32,
            compressed_page_size: payload.len() as i32,
            data_page: None,
            data_page_v2: Some(thrift::DataPageHeaderV2 {
                num_values: 4,
                num_nulls: 0,
                num_rows: 4,
                encoding: ENCODING_RLE,
                definition_levels_byte_length: 0,
                repetition_levels_byte_length: 0,
                is_compressed: false,
            }),
            dictionary_page: None,
        });
        let mut pages = header;
        pages.extend_from_slice(&payload);
        let schema = vec![
            SchemaElement {
                name: "schema".to_string(),
                num_children: Some(1),
                repetition_type: Some(REPETITION_REQUIRED),
                ..SchemaElement::default()
            },
            SchemaElement {
                type_: Some(0),
                repetition_type: Some(REPETITION_REQUIRED),
                name: "b".to_string(),
                ..SchemaElement::default()
            },
        ];
        let chunk_meta = ColumnMetaData {
            type_: 0,
            encodings: vec![ENCODING_RLE],
            path_in_schema: vec!["b".to_string()],
            codec: CODEC_UNCOMPRESSED,
            num_values: 4,
            total_uncompressed_size: pages.len() as i64,
            total_compressed_size: pages.len() as i64,
            data_page_offset: 4,
            dictionary_page_offset: None,
            statistics: None,
            ..Default::default()
        };
        let file = build_synthetic_parquet(schema, 4, chunk_meta, &pages);
        let (_, cols) = read_parquet_bytes(&file).expect("read rle booleans");
        assert_eq!(cols, vec![Column::Bool(vec![true, true, false, true])]);
    }

    #[test]
    fn delta_binary_packed_page_decodes() {
        // Values: [100, 105, 105, 110] -> deltas [5, 0, 5], min_delta = 0, offsets [5, 0, 5] (3 bits).
        let mut payload = Vec::new();
        thrift::write_varint(&mut payload, 128);
        thrift::write_varint(&mut payload, 4);
        thrift::write_varint(&mut payload, 4);
        thrift::write_zigzag(&mut payload, 100);
        thrift::write_zigzag(&mut payload, 0);
        payload.extend_from_slice(&[3, 0, 0, 0]);
        let mut mb = [0u8; 12];
        for (i, &v) in [5u64, 0, 5].iter().enumerate() {
            let bit_pos = i * 3;
            for b in 0..3 {
                if (v >> b) & 1 == 1 {
                    mb[(bit_pos + b) / 8] |= 1 << ((bit_pos + b) % 8);
                }
            }
        }
        payload.extend_from_slice(&mb);

        let header = thrift::encode_page_header(&PageHeader {
            type_: PAGE_DATA,
            uncompressed_page_size: payload.len() as i32,
            compressed_page_size: payload.len() as i32,
            data_page: Some(DataPageHeader {
                num_values: 4,
                encoding: ENCODING_DELTA_BINARY_PACKED,
                definition_level_encoding: ENCODING_RLE,
                repetition_level_encoding: ENCODING_RLE,
            }),
            data_page_v2: None,
            dictionary_page: None,
        });
        let mut pages = header;
        pages.extend_from_slice(&payload);

        let schema = vec![
            SchemaElement {
                name: "schema".to_string(),
                num_children: Some(1),
                repetition_type: Some(REPETITION_REQUIRED),
                ..SchemaElement::default()
            },
            SchemaElement {
                type_: Some(2),
                repetition_type: Some(REPETITION_REQUIRED),
                name: "ts".to_string(),
                ..SchemaElement::default()
            },
        ];
        let chunk_meta = ColumnMetaData {
            type_: 2,
            encodings: vec![ENCODING_DELTA_BINARY_PACKED],
            path_in_schema: vec!["ts".to_string()],
            codec: CODEC_UNCOMPRESSED,
            num_values: 4,
            total_uncompressed_size: pages.len() as i64,
            total_compressed_size: pages.len() as i64,
            data_page_offset: 4,
            dictionary_page_offset: None,
            statistics: None,
            ..Default::default()
        };
        let file = build_synthetic_parquet(schema, 4, chunk_meta, &pages);
        let (_, cols) = read_parquet_bytes(&file).expect("read delta");
        assert_eq!(cols, vec![Column::Int(vec![100, 105, 105, 110])]);
    }

    #[test]
    fn delta_length_byte_array_page_decodes() {
        // Strings ["ab", "cde", "f"] -> lengths [2, 3, 1], deltas [1, -2], min = -2, offsets [3, 0] (2 bits).
        let mut payload = Vec::new();
        thrift::write_varint(&mut payload, 128);
        thrift::write_varint(&mut payload, 4);
        thrift::write_varint(&mut payload, 3);
        thrift::write_zigzag(&mut payload, 2);
        thrift::write_zigzag(&mut payload, -2);
        payload.extend_from_slice(&[2, 0, 0, 0]);
        let mut mb = [0u8; 8];
        mb[0] = 0b00_00_00_11; // offset[0] = 3, offset[1] = 0
        payload.extend_from_slice(&mb);
        payload.extend_from_slice(b"abcdef");

        let header = thrift::encode_page_header(&PageHeader {
            type_: PAGE_DATA,
            uncompressed_page_size: payload.len() as i32,
            compressed_page_size: payload.len() as i32,
            data_page: Some(DataPageHeader {
                num_values: 3,
                encoding: ENCODING_DELTA_LENGTH_BYTE_ARRAY,
                definition_level_encoding: ENCODING_RLE,
                repetition_level_encoding: ENCODING_RLE,
            }),
            data_page_v2: None,
            dictionary_page: None,
        });
        let mut pages = header;
        pages.extend_from_slice(&payload);

        let schema = vec![
            SchemaElement {
                name: "schema".to_string(),
                num_children: Some(1),
                repetition_type: Some(REPETITION_REQUIRED),
                ..SchemaElement::default()
            },
            SchemaElement {
                type_: Some(6),
                repetition_type: Some(REPETITION_REQUIRED),
                name: "s".to_string(),
                converted_type: Some(CONVERTED_UTF8),
                ..SchemaElement::default()
            },
        ];
        let chunk_meta = ColumnMetaData {
            type_: 6,
            encodings: vec![ENCODING_DELTA_LENGTH_BYTE_ARRAY],
            path_in_schema: vec!["s".to_string()],
            codec: CODEC_UNCOMPRESSED,
            num_values: 3,
            total_uncompressed_size: pages.len() as i64,
            total_compressed_size: pages.len() as i64,
            data_page_offset: 4,
            dictionary_page_offset: None,
            statistics: None,
            ..Default::default()
        };
        let file = build_synthetic_parquet(schema, 3, chunk_meta, &pages);
        let (_, cols) = read_parquet_bytes(&file).expect("read delta len");
        assert_eq!(
            cols,
            vec![Column::Str(vec![
                "ab".to_string(),
                "cde".to_string(),
                "f".to_string(),
            ])]
        );
    }

    #[test]
    fn byte_stream_split_double_page_decodes() {
        let values = [1.25f64, -3.5, 42.0];
        let raw: Vec<[u8; 8]> = values.iter().map(|v| v.to_le_bytes()).collect();
        let mut payload = Vec::new();
        for b in 0..8 {
            for item in &raw {
                payload.push(item[b]);
            }
        }
        let header = thrift::encode_page_header(&PageHeader {
            type_: PAGE_DATA,
            uncompressed_page_size: payload.len() as i32,
            compressed_page_size: payload.len() as i32,
            data_page: Some(DataPageHeader {
                num_values: 3,
                encoding: ENCODING_BYTE_STREAM_SPLIT,
                definition_level_encoding: ENCODING_RLE,
                repetition_level_encoding: ENCODING_RLE,
            }),
            data_page_v2: None,
            dictionary_page: None,
        });
        let mut pages = header;
        pages.extend_from_slice(&payload);

        let schema = vec![
            SchemaElement {
                name: "schema".to_string(),
                num_children: Some(1),
                repetition_type: Some(REPETITION_REQUIRED),
                ..SchemaElement::default()
            },
            SchemaElement {
                type_: Some(5),
                repetition_type: Some(REPETITION_REQUIRED),
                name: "x".to_string(),
                ..SchemaElement::default()
            },
        ];
        let chunk_meta = ColumnMetaData {
            type_: 5,
            encodings: vec![ENCODING_BYTE_STREAM_SPLIT],
            path_in_schema: vec!["x".to_string()],
            codec: CODEC_UNCOMPRESSED,
            num_values: 3,
            total_uncompressed_size: pages.len() as i64,
            total_compressed_size: pages.len() as i64,
            data_page_offset: 4,
            dictionary_page_offset: None,
            statistics: None,
            ..Default::default()
        };
        let file = build_synthetic_parquet(schema, 3, chunk_meta, &pages);
        let (_, cols) = read_parquet_bytes(&file).expect("read bss");
        assert_eq!(cols, vec![Column::Float(vec![1.25, -3.5, 42.0])]);
    }

    #[test]
    fn nested_optional_struct_group_reads_with_multibit_def_levels() {
        // Schema: root -> optional group `user` -> optional int32 `age`
        // max_def_level = 2 (0 = user null, 1 = age null, 2 = present), bit_width = 2.
        // 4 rows with def levels [2, 0, 1, 2] -> values [30, 40].
        // Bit-packed 1 group of 8 2-bit values: [2, 0, 1, 2, 0, 0, 0, 0]
        // byte0 = 2 | (0<<2) | (1<<4) | (2<<6) = 0b10_01_00_10 = 0x92, byte1 = 0x00.
        let level_stream = [(1u8 << 1) | 1, 0x92, 0x00];
        let mut payload = Vec::new();
        payload.extend_from_slice(&(level_stream.len() as u32).to_le_bytes());
        payload.extend_from_slice(&level_stream);
        payload.extend_from_slice(&plain::encode_i32(&[30, 40]));

        let header = thrift::encode_page_header(&PageHeader {
            type_: PAGE_DATA,
            uncompressed_page_size: payload.len() as i32,
            compressed_page_size: payload.len() as i32,
            data_page: Some(DataPageHeader {
                num_values: 4,
                encoding: ENCODING_PLAIN,
                definition_level_encoding: ENCODING_RLE,
                repetition_level_encoding: ENCODING_RLE,
            }),
            data_page_v2: None,
            dictionary_page: None,
        });
        let mut pages = header;
        pages.extend_from_slice(&payload);

        let schema = vec![
            SchemaElement {
                name: "schema".to_string(),
                num_children: Some(1),
                repetition_type: Some(REPETITION_REQUIRED),
                ..SchemaElement::default()
            },
            SchemaElement {
                type_: None,
                repetition_type: Some(REPETITION_OPTIONAL),
                name: "user".to_string(),
                num_children: Some(1),
                ..SchemaElement::default()
            },
            SchemaElement {
                type_: Some(1),
                repetition_type: Some(REPETITION_OPTIONAL),
                name: "age".to_string(),
                ..SchemaElement::default()
            },
        ];
        let chunk_meta = ColumnMetaData {
            type_: 1,
            encodings: vec![ENCODING_PLAIN, ENCODING_RLE],
            path_in_schema: vec!["user".to_string(), "age".to_string()],
            codec: CODEC_UNCOMPRESSED,
            num_values: 4,
            total_uncompressed_size: pages.len() as i64,
            total_compressed_size: pages.len() as i64,
            data_page_offset: 4,
            dictionary_page_offset: None,
            statistics: None,
            ..Default::default()
        };
        let file = build_synthetic_parquet(schema, 4, chunk_meta, &pages);
        let (read_schema, cols) = read_parquet_bytes(&file).expect("read nested");
        assert_eq!(
            read_schema.fields,
            vec![Field {
                name: "user.age".to_string(),
                physical: Physical::Int32,
                optional: true,
            }]
        );
        assert_eq!(
            cols,
            vec![Column::Nullable(
                Box::new(Column::Int(vec![30, 0, 0, 40])),
                vec![false, true, true, false],
            )]
        );
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
