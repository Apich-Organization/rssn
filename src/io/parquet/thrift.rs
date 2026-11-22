//! Hand-written Thrift compact protocol for the Parquet metadata structures.
//!
//! Only the subset Parquet needs is implemented: structs, lists, binary
//! strings and the integer/boolean types. Every decoder validates lengths
//! against the remaining input before slicing, so a hostile file cannot make
//! it allocate or read out of bounds.

/// Compact-protocol type for `true` inside a field header.
pub(super) const T_BOOL_TRUE: u8 = 1;
/// Compact-protocol type for `false` inside a field header.
pub(super) const T_BOOL_FALSE: u8 = 2;
/// Compact-protocol type for a single byte.
pub(super) const T_BYTE: u8 = 3;
/// Compact-protocol type for a 16-bit integer.
pub(super) const T_I16: u8 = 4;
/// Compact-protocol type for a 32-bit integer.
pub(super) const T_I32: u8 = 5;
/// Compact-protocol type for a 64-bit integer.
pub(super) const T_I64: u8 = 6;
/// Compact-protocol type for a double.
pub(super) const T_DOUBLE: u8 = 7;
/// Compact-protocol type for binary data.
pub(super) const T_BINARY: u8 = 8;
/// Compact-protocol type for a list.
pub(super) const T_LIST: u8 = 9;
/// Compact-protocol type for a set.
pub(super) const T_SET: u8 = 10;
/// Compact-protocol type for a map.
pub(super) const T_MAP: u8 = 11;
/// Compact-protocol type for a nested struct.
pub(super) const T_STRUCT: u8 = 12;

/// Largest list accepted from a file, guarding against hostile lengths.
const MAX_LIST: usize = 16_000_000;

/// Writes an unsigned LEB128 varint.
pub(super) fn write_varint(
    out: &mut Vec<u8>,
    mut value: u64,
) {
    loop {
        let byte = (value & 0x7f) as u8;
        value >>= 7;
        if value == 0 {
            out.push(byte);
            break;
        }
        out.push(byte | 0x80);
    }
}

/// Writes a zigzag-encoded signed integer.
pub(super) fn write_zigzag(
    out: &mut Vec<u8>,
    value: i64,
) {
    let encoded = ((value << 1) ^ (value >> 63)) as u64;
    write_varint(out, encoded);
}

/// Writes a field header, using the short delta form when possible.
pub(super) fn write_field(
    out: &mut Vec<u8>,
    last: &mut i16,
    id: i16,
    kind: u8,
) {
    let delta = id.wrapping_sub(*last);
    if delta > 0 && delta <= 15 {
        out.push(((delta as u8) << 4) | kind);
    } else {
        out.push(kind);
        write_zigzag(out, i64::from(id));
    }
    *last = id;
}

/// Writes a length-prefixed byte string.
pub(super) fn write_binary(
    out: &mut Vec<u8>,
    bytes: &[u8],
) {
    write_varint(out, bytes.len() as u64);
    out.extend_from_slice(bytes);
}

/// Writes a list header; the caller then writes the elements.
pub(super) fn write_list(
    out: &mut Vec<u8>,
    element: u8,
    len: usize,
) {
    if len < 15 {
        out.push(((len as u8) << 4) | element);
    } else {
        out.push(0xF0 | element);
        write_varint(out, len as u64);
    }
}

/// A leaf or group entry of the Parquet schema.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub(super) struct SchemaElement {
    /// Physical type code, absent for groups.
    pub type_: Option<i32>,
    /// Fixed length for `FIXED_LEN_BYTE_ARRAY`.
    pub type_length: Option<i32>,
    /// Repetition code: 0 required, 1 optional, 2 repeated.
    pub repetition_type: Option<i32>,
    /// Field name.
    pub name: String,
    /// Number of children for a group.
    pub num_children: Option<i32>,
    /// Converted type code (0 is UTF-8).
    pub converted_type: Option<i32>,
    /// Decimal scale.
    pub scale: Option<i32>,
    /// Decimal precision.
    pub precision: Option<i32>,
    /// Field id.
    pub field_id: Option<i32>,
}

/// Column statistics.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub(super) struct Statistics {
    /// Deprecated maximum.
    pub max: Option<Vec<u8>>,
    /// Deprecated minimum.
    pub min: Option<Vec<u8>>,
    /// Number of missing values.
    pub null_count: Option<i64>,
    /// Maximum, PLAIN encoded.
    pub max_value: Option<Vec<u8>>,
    /// Minimum, PLAIN encoded.
    pub min_value: Option<Vec<u8>>,
}

/// Metadata for one column chunk.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub(super) struct ColumnMetaData {
    /// Physical type code.
    pub type_: i32,
    /// Encodings used in the chunk.
    pub encodings: Vec<i32>,
    /// Dotted path of the column.
    pub path_in_schema: Vec<String>,
    /// Compression codec code.
    pub codec: i32,
    /// Number of values, including missing ones.
    pub num_values: i64,
    /// Uncompressed size including page headers.
    pub total_uncompressed_size: i64,
    /// Compressed size including page headers.
    pub total_compressed_size: i64,
    /// Offset of the first data page.
    pub data_page_offset: i64,
    /// Offset of the dictionary page, when present.
    pub dictionary_page_offset: Option<i64>,
    /// Column statistics.
    pub statistics: Option<Statistics>,
}

/// One column chunk inside a row group.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub(super) struct ColumnChunk {
    /// Offset of the chunk.
    pub file_offset: i64,
    /// Chunk metadata.
    pub meta_data: Option<ColumnMetaData>,
}

/// A row group.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub(super) struct RowGroup {
    /// Column chunks, in schema order.
    pub columns: Vec<ColumnChunk>,
    /// Total uncompressed size.
    pub total_byte_size: i64,
    /// Number of rows.
    pub num_rows: i64,
    /// Offset of the row group.
    pub file_offset: Option<i64>,
    /// Total compressed size.
    pub total_compressed_size: Option<i64>,
}

/// The file footer.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub(super) struct FileMetaData {
    /// Format version.
    pub version: i32,
    /// Root element followed by one element per field.
    pub schema: Vec<SchemaElement>,
    /// Number of rows in the file.
    pub num_rows: i64,
    /// Row groups.
    pub row_groups: Vec<RowGroup>,
}

/// Header of a data page.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub(super) struct DataPageHeader {
    /// Number of values, including missing ones.
    pub num_values: i32,
    /// Value encoding code.
    pub encoding: i32,
    /// Definition-level encoding code.
    pub definition_level_encoding: i32,
    /// Repetition-level encoding code.
    pub repetition_level_encoding: i32,
}

/// Header preceding every page.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub(super) struct PageHeader {
    /// Page type code: 0 is a data page.
    pub type_: i32,
    /// Uncompressed payload size.
    pub uncompressed_page_size: i32,
    /// Compressed payload size.
    pub compressed_page_size: i32,
    /// Data-page details, present for data pages.
    pub data_page: Option<DataPageHeader>,
}

/// Encodes a [`FileMetaData`] value.
pub(super) fn encode_file_metadata(metadata: &FileMetaData) -> Vec<u8> {
    let mut out = Vec::new();
    let mut last = 0i16;
    write_field(&mut out, &mut last, 1, T_I32);
    write_zigzag(&mut out, i64::from(metadata.version));
    write_field(&mut out, &mut last, 2, T_LIST);
    write_list(&mut out, T_STRUCT, metadata.schema.len());
    for element in &metadata.schema {
        encode_schema_element(&mut out, element);
    }
    write_field(&mut out, &mut last, 3, T_I64);
    write_zigzag(&mut out, metadata.num_rows);
    write_field(&mut out, &mut last, 4, T_LIST);
    write_list(&mut out, T_STRUCT, metadata.row_groups.len());
    for row_group in &metadata.row_groups {
        encode_row_group(&mut out, row_group);
    }
    write_field(&mut out, &mut last, 6, T_BINARY);
    write_binary(&mut out, b"rssn");
    out.push(0);
    out
}

/// Encodes a schema element inline.
fn encode_schema_element(
    out: &mut Vec<u8>,
    element: &SchemaElement,
) {
    let mut last = 0i16;
    if let Some(value) = element.type_ {
        write_field(out, &mut last, 1, T_I32);
        write_zigzag(out, i64::from(value));
    }
    if let Some(value) = element.type_length {
        write_field(out, &mut last, 2, T_I32);
        write_zigzag(out, i64::from(value));
    }
    if let Some(value) = element.repetition_type {
        write_field(out, &mut last, 3, T_I32);
        write_zigzag(out, i64::from(value));
    }
    write_field(out, &mut last, 4, T_BINARY);
    write_binary(out, element.name.as_bytes());
    if let Some(value) = element.num_children {
        write_field(out, &mut last, 5, T_I32);
        write_zigzag(out, i64::from(value));
    }
    if let Some(value) = element.converted_type {
        write_field(out, &mut last, 6, T_I32);
        write_zigzag(out, i64::from(value));
    }
    if let Some(value) = element.scale {
        write_field(out, &mut last, 7, T_I32);
        write_zigzag(out, i64::from(value));
    }
    if let Some(value) = element.precision {
        write_field(out, &mut last, 8, T_I32);
        write_zigzag(out, i64::from(value));
    }
    if let Some(value) = element.field_id {
        write_field(out, &mut last, 9, T_I32);
        write_zigzag(out, i64::from(value));
    }
    out.push(0);
}

/// Encodes a row group inline.
fn encode_row_group(
    out: &mut Vec<u8>,
    row_group: &RowGroup,
) {
    let mut last = 0i16;
    write_field(out, &mut last, 1, T_LIST);
    write_list(out, T_STRUCT, row_group.columns.len());
    for column in &row_group.columns {
        encode_column_chunk(out, column);
    }
    write_field(out, &mut last, 2, T_I64);
    write_zigzag(out, row_group.total_byte_size);
    write_field(out, &mut last, 3, T_I64);
    write_zigzag(out, row_group.num_rows);
    if let Some(value) = row_group.file_offset {
        write_field(out, &mut last, 5, T_I64);
        write_zigzag(out, value);
    }
    if let Some(value) = row_group.total_compressed_size {
        write_field(out, &mut last, 6, T_I64);
        write_zigzag(out, value);
    }
    out.push(0);
}

/// Encodes a column chunk inline.
fn encode_column_chunk(
    out: &mut Vec<u8>,
    chunk: &ColumnChunk,
) {
    let mut last = 0i16;
    write_field(out, &mut last, 2, T_I64);
    write_zigzag(out, chunk.file_offset);
    if let Some(metadata) = &chunk.meta_data {
        write_field(out, &mut last, 3, T_STRUCT);
        encode_column_metadata(out, metadata);
    }
    out.push(0);
}

/// Encodes column metadata inline.
fn encode_column_metadata(
    out: &mut Vec<u8>,
    metadata: &ColumnMetaData,
) {
    let mut last = 0i16;
    write_field(out, &mut last, 1, T_I32);
    write_zigzag(out, i64::from(metadata.type_));
    write_field(out, &mut last, 2, T_LIST);
    write_list(out, T_I32, metadata.encodings.len());
    for encoding in &metadata.encodings {
        write_zigzag(out, i64::from(*encoding));
    }
    write_field(out, &mut last, 3, T_LIST);
    write_list(out, T_BINARY, metadata.path_in_schema.len());
    for path in &metadata.path_in_schema {
        write_binary(out, path.as_bytes());
    }
    write_field(out, &mut last, 4, T_I32);
    write_zigzag(out, i64::from(metadata.codec));
    write_field(out, &mut last, 5, T_I64);
    write_zigzag(out, metadata.num_values);
    write_field(out, &mut last, 6, T_I64);
    write_zigzag(out, metadata.total_uncompressed_size);
    write_field(out, &mut last, 7, T_I64);
    write_zigzag(out, metadata.total_compressed_size);
    write_field(out, &mut last, 9, T_I64);
    write_zigzag(out, metadata.data_page_offset);
    if let Some(value) = metadata.dictionary_page_offset {
        write_field(out, &mut last, 11, T_I64);
        write_zigzag(out, value);
    }
    if let Some(statistics) = &metadata.statistics {
        write_field(out, &mut last, 12, T_STRUCT);
        encode_statistics(out, statistics);
    }
    out.push(0);
}

/// Encodes statistics inline.
fn encode_statistics(
    out: &mut Vec<u8>,
    statistics: &Statistics,
) {
    let mut last = 0i16;
    if let Some(value) = &statistics.max {
        write_field(out, &mut last, 1, T_BINARY);
        write_binary(out, value);
    }
    if let Some(value) = &statistics.min {
        write_field(out, &mut last, 2, T_BINARY);
        write_binary(out, value);
    }
    if let Some(value) = statistics.null_count {
        write_field(out, &mut last, 3, T_I64);
        write_zigzag(out, value);
    }
    if let Some(value) = &statistics.max_value {
        write_field(out, &mut last, 5, T_BINARY);
        write_binary(out, value);
    }
    if let Some(value) = &statistics.min_value {
        write_field(out, &mut last, 6, T_BINARY);
        write_binary(out, value);
    }
    out.push(0);
}

/// Encodes a [`PageHeader`] value.
pub(super) fn encode_page_header(header: &PageHeader) -> Vec<u8> {
    let mut out = Vec::new();
    let mut last = 0i16;
    write_field(&mut out, &mut last, 1, T_I32);
    write_zigzag(&mut out, i64::from(header.type_));
    write_field(&mut out, &mut last, 2, T_I32);
    write_zigzag(&mut out, i64::from(header.uncompressed_page_size));
    write_field(&mut out, &mut last, 3, T_I32);
    write_zigzag(&mut out, i64::from(header.compressed_page_size));
    if let Some(data) = &header.data_page {
        write_field(&mut out, &mut last, 5, T_STRUCT);
        let mut inner = 0i16;
        write_field(&mut out, &mut inner, 1, T_I32);
        write_zigzag(&mut out, i64::from(data.num_values));
        write_field(&mut out, &mut inner, 2, T_I32);
        write_zigzag(&mut out, i64::from(data.encoding));
        write_field(&mut out, &mut inner, 3, T_I32);
        write_zigzag(&mut out, i64::from(data.definition_level_encoding));
        write_field(&mut out, &mut inner, 4, T_I32);
        write_zigzag(&mut out, i64::from(data.repetition_level_encoding));
        out.push(0);
    }
    out.push(0);
    out
}

/// A bounds-checked compact-protocol reader.
pub(super) struct Decoder<'a> {
    /// The buffer being decoded.
    buf: &'a [u8],
    /// Current offset.
    pos: usize,
}

impl<'a> Decoder<'a> {
    /// Starts decoding at the beginning of `buf`.
    pub(super) fn new(buf: &'a [u8]) -> Self {
        Self { buf, pos: 0 }
    }

    /// Returns the number of bytes consumed so far.
    pub(super) fn pos(&self) -> usize {
        self.pos
    }

    /// Returns the number of unread bytes.
    fn remaining(&self) -> usize {
        self.buf.len().saturating_sub(self.pos)
    }

    /// Reads one byte.
    fn byte(&mut self) -> Result<u8, String> {
        let value = *self
            .buf
            .get(self.pos)
            .ok_or_else(|| "unexpected end of thrift data".to_string())?;
        self.pos += 1;
        Ok(value)
    }

    /// Reads an unsigned varint.
    fn varint(&mut self) -> Result<u64, String> {
        let mut result = 0u64;
        let mut shift = 0u32;
        loop {
            if shift >= 64 {
                return Err("thrift varint is too long".to_string());
            }
            let byte = self.byte()?;
            result |= u64::from(byte & 0x7f) << shift;
            if byte & 0x80 == 0 {
                return Ok(result);
            }
            shift += 7;
        }
    }

    /// Reads a zigzag-encoded signed integer.
    fn zigzag(&mut self) -> Result<i64, String> {
        let value = self.varint()?;
        Ok(((value >> 1) as i64) ^ -((value & 1) as i64))
    }

    /// Reads a 32-bit integer.
    fn i32(&mut self) -> Result<i32, String> {
        i32::try_from(self.zigzag()?).map_err(|_| "thrift integer is out of range".to_string())
    }

    /// Reads a 64-bit integer.
    fn i64(&mut self) -> Result<i64, String> {
        self.zigzag()
    }

    /// Reads a length-prefixed byte string.
    fn binary(&mut self) -> Result<&'a [u8], String> {
        let len = usize::try_from(self.varint()?)
            .map_err(|_| "thrift binary length is out of range".to_string())?;
        let end = self
            .pos
            .checked_add(len)
            .ok_or_else(|| "thrift binary length overflows".to_string())?;
        let value = self
            .buf
            .get(self.pos..end)
            .ok_or_else(|| "thrift binary value is truncated".to_string())?;
        self.pos = end;
        Ok(value)
    }

    /// Reads a length-prefixed UTF-8 string.
    fn string(&mut self) -> Result<String, String> {
        String::from_utf8(self.binary()?.to_vec())
            .map_err(|_| "thrift string is not valid utf-8".to_string())
    }

    /// Reads a list header, returning the element type and length.
    fn list_header(&mut self) -> Result<(u8, usize), String> {
        let byte = self.byte()?;
        let element = byte & 0x0f;
        let short = usize::from(byte >> 4);
        let len = if short == 15 {
            usize::try_from(self.varint()?)
                .map_err(|_| "thrift list length is out of range".to_string())?
        } else {
            short
        };
        if len > MAX_LIST || len > self.remaining().saturating_add(1) {
            return Err("thrift list length exceeds the remaining data".to_string());
        }
        Ok((element, len))
    }

    /// Reads one field header, or `None` at the end of a struct.
    fn field_header(
        &mut self,
        last: &mut i16,
    ) -> Result<Option<(i16, u8)>, String> {
        let byte = self.byte()?;
        if byte == 0 {
            return Ok(None);
        }
        let kind = byte & 0x0f;
        let delta = i16::from(byte >> 4);
        let id = if delta == 0 {
            i16::try_from(self.zigzag()?)
                .map_err(|_| "thrift field id is out of range".to_string())?
        } else {
            last.checked_add(delta)
                .ok_or_else(|| "thrift field id overflows".to_string())?
        };
        *last = id;
        Ok(Some((id, kind)))
    }

    /// Skips a value of the given compact-protocol type.
    fn skip(
        &mut self,
        kind: u8,
    ) -> Result<(), String> {
        match kind {
            | T_BOOL_TRUE | T_BOOL_FALSE => Ok(()),
            | T_BYTE => {
                self.byte()?;
                Ok(())
            },
            | T_I16 | T_I32 | T_I64 => {
                self.varint()?;
                Ok(())
            },
            | T_DOUBLE => {
                let end = self
                    .pos
                    .checked_add(8)
                    .ok_or_else(|| "thrift double overflows".to_string())?;
                if end > self.buf.len() {
                    return Err("thrift double is truncated".to_string());
                }
                self.pos = end;
                Ok(())
            },
            | T_BINARY => {
                self.binary()?;
                Ok(())
            },
            | T_LIST | T_SET => {
                let (element, len) = self.list_header()?;
                for _ in 0..len {
                    self.skip(element)?;
                }
                Ok(())
            },
            | T_MAP => {
                let header = self.byte()?;
                if header == 0 {
                    return Ok(());
                }
                let key_kind = header >> 4;
                let value_kind = header & 0x0f;
                let len = usize::try_from(self.varint()?)
                    .map_err(|_| "thrift map length is out of range".to_string())?;
                if len > MAX_LIST {
                    return Err("thrift map is too large".to_string());
                }
                for _ in 0..len {
                    self.skip(key_kind)?;
                    self.skip(value_kind)?;
                }
                Ok(())
            },
            | T_STRUCT => {
                let mut last = 0i16;
                while let Some((_, field)) = self.field_header(&mut last)? {
                    self.skip(field)?;
                }
                Ok(())
            },
            | _ => Err(format!("unknown thrift type {kind}")),
        }
    }
}

/// Checks that a field has the expected compact-protocol type.
fn expect(
    kind: u8,
    wanted: u8,
) -> Result<(), String> {
    if kind == wanted {
        Ok(())
    } else {
        Err(format!("unexpected thrift type {kind}, expected {wanted}"))
    }
}

/// Decodes a [`FileMetaData`] value from the footer bytes.
pub(super) fn decode_file_metadata(bytes: &[u8]) -> Result<FileMetaData, String> {
    let mut decoder = Decoder::new(bytes);
    let metadata = parse_file_metadata(&mut decoder)?;
    Ok(metadata)
}

/// Parses the file footer.
fn parse_file_metadata(decoder: &mut Decoder<'_>) -> Result<FileMetaData, String> {
    let mut metadata = FileMetaData::default();
    let mut last = 0i16;
    while let Some((id, kind)) = decoder.field_header(&mut last)? {
        match id {
            | 1 => {
                expect(kind, T_I32)?;
                metadata.version = decoder.i32()?;
            },
            | 2 => {
                expect(kind, T_LIST)?;
                let (element, len) = decoder.list_header()?;
                if element != T_STRUCT {
                    return Err("parquet schema list must hold structs".to_string());
                }
                metadata.schema.reserve(len.min(1024));
                for _ in 0..len {
                    metadata.schema.push(parse_schema_element(decoder)?);
                }
            },
            | 3 => {
                expect(kind, T_I64)?;
                metadata.num_rows = decoder.i64()?;
            },
            | 4 => {
                expect(kind, T_LIST)?;
                let (element, len) = decoder.list_header()?;
                if element != T_STRUCT {
                    return Err("parquet row group list must hold structs".to_string());
                }
                metadata.row_groups.reserve(len.min(1024));
                for _ in 0..len {
                    metadata.row_groups.push(parse_row_group(decoder)?);
                }
            },
            | _ => decoder.skip(kind)?,
        }
    }
    Ok(metadata)
}

/// Parses one schema element.
fn parse_schema_element(decoder: &mut Decoder<'_>) -> Result<SchemaElement, String> {
    let mut element = SchemaElement::default();
    let mut last = 0i16;
    while let Some((id, kind)) = decoder.field_header(&mut last)? {
        match id {
            | 1 => {
                expect(kind, T_I32)?;
                element.type_ = Some(decoder.i32()?);
            },
            | 2 => {
                expect(kind, T_I32)?;
                element.type_length = Some(decoder.i32()?);
            },
            | 3 => {
                expect(kind, T_I32)?;
                element.repetition_type = Some(decoder.i32()?);
            },
            | 4 => {
                expect(kind, T_BINARY)?;
                element.name = decoder.string()?;
            },
            | 5 => {
                expect(kind, T_I32)?;
                element.num_children = Some(decoder.i32()?);
            },
            | 6 => {
                expect(kind, T_I32)?;
                element.converted_type = Some(decoder.i32()?);
            },
            | 7 => {
                expect(kind, T_I32)?;
                element.scale = Some(decoder.i32()?);
            },
            | 8 => {
                expect(kind, T_I32)?;
                element.precision = Some(decoder.i32()?);
            },
            | 9 => {
                expect(kind, T_I32)?;
                element.field_id = Some(decoder.i32()?);
            },
            | _ => decoder.skip(kind)?,
        }
    }
    Ok(element)
}

/// Parses one row group.
fn parse_row_group(decoder: &mut Decoder<'_>) -> Result<RowGroup, String> {
    let mut row_group = RowGroup::default();
    let mut last = 0i16;
    while let Some((id, kind)) = decoder.field_header(&mut last)? {
        match id {
            | 1 => {
                expect(kind, T_LIST)?;
                let (element, len) = decoder.list_header()?;
                if element != T_STRUCT {
                    return Err("column chunk list must hold structs".to_string());
                }
                row_group.columns.reserve(len.min(1024));
                for _ in 0..len {
                    row_group.columns.push(parse_column_chunk(decoder)?);
                }
            },
            | 2 => {
                expect(kind, T_I64)?;
                row_group.total_byte_size = decoder.i64()?;
            },
            | 3 => {
                expect(kind, T_I64)?;
                row_group.num_rows = decoder.i64()?;
            },
            | 5 => {
                expect(kind, T_I64)?;
                row_group.file_offset = Some(decoder.i64()?);
            },
            | 6 => {
                expect(kind, T_I64)?;
                row_group.total_compressed_size = Some(decoder.i64()?);
            },
            | _ => decoder.skip(kind)?,
        }
    }
    Ok(row_group)
}

/// Parses one column chunk.
fn parse_column_chunk(decoder: &mut Decoder<'_>) -> Result<ColumnChunk, String> {
    let mut chunk = ColumnChunk::default();
    let mut last = 0i16;
    while let Some((id, kind)) = decoder.field_header(&mut last)? {
        match id {
            | 2 => {
                expect(kind, T_I64)?;
                chunk.file_offset = decoder.i64()?;
            },
            | 3 => {
                expect(kind, T_STRUCT)?;
                chunk.meta_data = Some(parse_column_metadata(decoder)?);
            },
            | _ => decoder.skip(kind)?,
        }
    }
    Ok(chunk)
}

/// Parses column chunk metadata.
fn parse_column_metadata(decoder: &mut Decoder<'_>) -> Result<ColumnMetaData, String> {
    let mut metadata = ColumnMetaData::default();
    let mut last = 0i16;
    while let Some((id, kind)) = decoder.field_header(&mut last)? {
        match id {
            | 1 => {
                expect(kind, T_I32)?;
                metadata.type_ = decoder.i32()?;
            },
            | 2 => {
                expect(kind, T_LIST)?;
                let (element, len) = decoder.list_header()?;
                if element != T_I32 {
                    return Err("encoding list must hold i32 values".to_string());
                }
                metadata.encodings.reserve(len.min(1024));
                for _ in 0..len {
                    metadata.encodings.push(decoder.i32()?);
                }
            },
            | 3 => {
                expect(kind, T_LIST)?;
                let (element, len) = decoder.list_header()?;
                if element != T_BINARY {
                    return Err("column path list must hold strings".to_string());
                }
                metadata.path_in_schema.reserve(len.min(1024));
                for _ in 0..len {
                    metadata.path_in_schema.push(decoder.string()?);
                }
            },
            | 4 => {
                expect(kind, T_I32)?;
                metadata.codec = decoder.i32()?;
            },
            | 5 => {
                expect(kind, T_I64)?;
                metadata.num_values = decoder.i64()?;
            },
            | 6 => {
                expect(kind, T_I64)?;
                metadata.total_uncompressed_size = decoder.i64()?;
            },
            | 7 => {
                expect(kind, T_I64)?;
                metadata.total_compressed_size = decoder.i64()?;
            },
            | 9 => {
                expect(kind, T_I64)?;
                metadata.data_page_offset = decoder.i64()?;
            },
            | 11 => {
                expect(kind, T_I64)?;
                metadata.dictionary_page_offset = Some(decoder.i64()?);
            },
            | 12 => {
                expect(kind, T_STRUCT)?;
                metadata.statistics = Some(parse_statistics(decoder)?);
            },
            | _ => decoder.skip(kind)?,
        }
    }
    Ok(metadata)
}

/// Parses column statistics.
fn parse_statistics(decoder: &mut Decoder<'_>) -> Result<Statistics, String> {
    let mut statistics = Statistics::default();
    let mut last = 0i16;
    while let Some((id, kind)) = decoder.field_header(&mut last)? {
        match id {
            | 1 => {
                expect(kind, T_BINARY)?;
                statistics.max = Some(decoder.binary()?.to_vec());
            },
            | 2 => {
                expect(kind, T_BINARY)?;
                statistics.min = Some(decoder.binary()?.to_vec());
            },
            | 3 => {
                expect(kind, T_I64)?;
                statistics.null_count = Some(decoder.i64()?);
            },
            | 5 => {
                expect(kind, T_BINARY)?;
                statistics.max_value = Some(decoder.binary()?.to_vec());
            },
            | 6 => {
                expect(kind, T_BINARY)?;
                statistics.min_value = Some(decoder.binary()?.to_vec());
            },
            | _ => decoder.skip(kind)?,
        }
    }
    Ok(statistics)
}

/// Decodes a [`PageHeader`] and reports how many bytes it occupied.
pub(super) fn decode_page_header(bytes: &[u8]) -> Result<(PageHeader, usize), String> {
    let mut decoder = Decoder::new(bytes);
    let header = parse_page_header(&mut decoder)?;
    Ok((header, decoder.pos()))
}

/// Parses a page header.
fn parse_page_header(decoder: &mut Decoder<'_>) -> Result<PageHeader, String> {
    let mut header = PageHeader::default();
    let mut last = 0i16;
    while let Some((id, kind)) = decoder.field_header(&mut last)? {
        match id {
            | 1 => {
                expect(kind, T_I32)?;
                header.type_ = decoder.i32()?;
            },
            | 2 => {
                expect(kind, T_I32)?;
                header.uncompressed_page_size = decoder.i32()?;
            },
            | 3 => {
                expect(kind, T_I32)?;
                header.compressed_page_size = decoder.i32()?;
            },
            | 5 => {
                expect(kind, T_STRUCT)?;
                header.data_page = Some(parse_data_page_header(decoder)?);
            },
            | _ => decoder.skip(kind)?,
        }
    }
    Ok(header)
}

/// Parses a data page header.
fn parse_data_page_header(decoder: &mut Decoder<'_>) -> Result<DataPageHeader, String> {
    let mut header = DataPageHeader::default();
    let mut last = 0i16;
    while let Some((id, kind)) = decoder.field_header(&mut last)? {
        match id {
            | 1 => {
                expect(kind, T_I32)?;
                header.num_values = decoder.i32()?;
            },
            | 2 => {
                expect(kind, T_I32)?;
                header.encoding = decoder.i32()?;
            },
            | 3 => {
                expect(kind, T_I32)?;
                header.definition_level_encoding = decoder.i32()?;
            },
            | 4 => {
                expect(kind, T_I32)?;
                header.repetition_level_encoding = decoder.i32()?;
            },
            | _ => decoder.skip(kind)?,
        }
    }
    Ok(header)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_metadata() -> FileMetaData {
        FileMetaData {
            version: 1,
            schema: vec![
                SchemaElement {
                    name: "schema".to_string(),
                    num_children: Some(1),
                    repetition_type: Some(0),
                    ..SchemaElement::default()
                },
                SchemaElement {
                    type_: Some(6),
                    repetition_type: Some(1),
                    name: "text".to_string(),
                    converted_type: Some(0),
                    ..SchemaElement::default()
                },
            ],
            num_rows: 3,
            row_groups: vec![RowGroup {
                columns: vec![ColumnChunk {
                    file_offset: 4,
                    meta_data: Some(ColumnMetaData {
                        type_: 6,
                        encodings: vec![0, 3],
                        path_in_schema: vec!["text".to_string()],
                        codec: 0,
                        num_values: 3,
                        total_uncompressed_size: 42,
                        total_compressed_size: 42,
                        data_page_offset: 4,
                        dictionary_page_offset: None,
                        statistics: Some(Statistics {
                            null_count: Some(1),
                            min: Some(vec![1, 0, 0, 0, b'a']),
                            max: Some(vec![1, 0, 0, 0, b'z']),
                            min_value: Some(vec![1, 0, 0, 0, b'a']),
                            max_value: Some(vec![1, 0, 0, 0, b'z']),
                        }),
                    }),
                }],
                total_byte_size: 42,
                num_rows: 3,
                file_offset: Some(4),
                total_compressed_size: Some(42),
            }],
        }
    }

    #[test]
    fn metadata_round_trip() {
        let metadata = sample_metadata();
        let encoded = encode_file_metadata(&metadata);
        let decoded = decode_file_metadata(&encoded).expect("decode");
        assert_eq!(decoded, metadata);
    }

    #[test]
    fn page_header_round_trip() {
        let header = PageHeader {
            type_: 0,
            uncompressed_page_size: 17,
            compressed_page_size: 17,
            data_page: Some(DataPageHeader {
                num_values: 5,
                encoding: 0,
                definition_level_encoding: 3,
                repetition_level_encoding: 3,
            }),
        };
        let encoded = encode_page_header(&header);
        let (decoded, used) = decode_page_header(&encoded).expect("decode");
        assert_eq!(decoded, header);
        assert_eq!(used, encoded.len());
    }

    #[test]
    fn truncated_metadata_is_rejected() {
        let encoded = encode_file_metadata(&sample_metadata());
        for cut in [1, 5, encoded.len() / 2] {
            let slice = encoded.get(..cut).expect("slice");
            assert!(decode_file_metadata(slice).is_err());
        }
    }

    #[test]
    fn empty_struct_decodes() {
        let metadata = decode_file_metadata(&[0]).expect("decode");
        assert_eq!(metadata, FileMetaData::default());
    }

    #[test]
    fn oversized_list_is_rejected() {
        // Field 1 (list of structs) with a length that cannot fit.
        let mut bytes = vec![0x19, 0xF0 | T_STRUCT];
        write_varint(&mut bytes, 1_000_000);
        assert!(decode_file_metadata(&bytes).is_err());
    }
}
