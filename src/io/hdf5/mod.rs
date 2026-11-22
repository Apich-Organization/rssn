//! HDF5 files: groups, attributes and contiguous or chunked datasets.
//!
//! A dependency-free reader and writer for a practical subset of the HDF5 file
//! format. The reader understands superblock versions 0, 2 and 3; object
//! header versions 1 and 2; old-style groups (a symbol table with a version-1
//! B-tree and a local heap) as well as compact link messages; contiguous,
//! compact and chunked datasets; the `shuffle` and `deflate` filters; numeric
//! and fixed-length string datatypes; and attributes stored in the object
//! header.
//!
//! Numeric dataset values are exposed as `f64` (integers convert exactly up to
//! 2^53). The writer emits superblock version 0 with a single root group and
//! contiguous, uncompressed datasets.
//!
//! # Limitations
//!
//! - New-style groups that use a fractal heap and a version-2 B-tree are
//!   detected and rejected with a clear error instead of being misparsed.
//! - Chunked datasets using version-4 layout or version-2 chunk B-trees are
//!   rejected, as are the `szip`, `nbit` and `scaleoffset` filters.
//! - Compound, enumerated, variable-length, array and reference datatypes are
//!   not decoded, and string *datasets* are rejected (string *attributes* are
//!   supported).
//! - External links, user blocks and the shared message table are ignored.
//!
//! Every offset, length and count read from a file is validated against the
//! buffer before use, and arrays are capped at 64 MiB and 2^32 elements.

mod inflate;

use std::path::Path;

/// Largest byte length allowed for one decoded dataset (64 MiB).
const MAX_ARRAY_BYTES: usize = 64 * 1024 * 1024;

/// Largest element count allowed for one dataset (2^32).
const MAX_ELEMS: u64 = 1 << 32;

/// Maximum group nesting depth accepted by the reader.
const MAX_DEPTH: usize = 32;

/// Maximum number of entries accepted in one B-tree node.
const MAX_BTREE_ENTRIES: usize = 65_536;

/// The datatype of an HDF5 dataset, restricted to the numeric and string
/// classes this module understands.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum H5Type {
    /// Signed 8-bit integer.
    I8,
    /// Unsigned 8-bit integer.
    U8,
    /// Signed 16-bit integer.
    I16,
    /// Unsigned 16-bit integer.
    U16,
    /// Signed 32-bit integer.
    I32,
    /// Unsigned 32-bit integer.
    U32,
    /// Signed 64-bit integer.
    I64,
    /// Unsigned 64-bit integer.
    U64,
    /// 32-bit IEEE floating point.
    F32,
    /// 64-bit IEEE floating point.
    F64,
    /// Fixed-length string (not supported for dataset data).
    Str,
}

/// A dataset: a path, its dimensions, element type and row-major values.
#[derive(Debug, Clone, PartialEq)]
pub struct Dataset {
    /// Absolute path inside the file, for example `/temperature`.
    pub path: String,
    /// Dimension sizes, outermost first. Empty for a scalar dataset.
    pub shape: Vec<usize>,
    /// Element type of the stored data.
    pub h5_type: H5Type,
    /// Row-major values, converted to `f64`.
    pub data: Vec<f64>,
}

/// One attribute attached to a group or dataset.
#[derive(Debug, Clone, PartialEq)]
pub struct Attribute {
    /// Attribute name.
    pub name: String,
    /// Attribute payload.
    pub value: AttrValue,
}

/// The payload of an [`Attribute`].
#[derive(Debug, Clone, PartialEq)]
pub enum AttrValue {
    /// A single signed integer.
    Int(i64),
    /// A single floating-point value.
    Float(f64),
    /// A fixed-length string.
    Text(String),
    /// A one-dimensional array of signed integers.
    Ints(Vec<i64>),
    /// A one-dimensional array of floating-point values.
    Floats(Vec<f64>),
}

/// The contents of an HDF5 file.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct H5File {
    /// Paths of every group, including the root group `/`.
    pub groups: Vec<String>,
    /// Every dataset found, in traversal order.
    pub datasets: Vec<Dataset>,
    /// Attributes keyed by the path of the object that owns them.
    pub attributes: Vec<(String, Vec<Attribute>)>,
}

// -------------------------------------------------------------------------
// Byte-level reader
// -------------------------------------------------------------------------

/// A bounds-checked view of a file image with the file's address sizes.
struct Reader<'a> {
    /// Whole file image.
    buf: &'a [u8],
    /// Width of an address in bytes.
    off_size: usize,
    /// Width of a length field in bytes.
    len_size: usize,
    /// Base address subtracted from stored addresses.
    base: u64,
}

impl<'a> Reader<'a> {
    /// Reads `len` bytes at `pos`, failing on any overrun.
    fn bytes(
        &self,
        pos: u64,
        len: u64,
    ) -> Result<&'a [u8], String> {
        let start = usize::try_from(pos).map_err(|_| format!("offset {pos} is out of range"))?;
        let count = usize::try_from(len).map_err(|_| format!("length {len} is out of range"))?;
        let end = start
            .checked_add(count)
            .ok_or_else(|| "offset arithmetic overflow".to_owned())?;
        self.buf
            .get(start..end)
            .ok_or_else(|| format!("read of {len} bytes at {pos} runs past the end of the file"))
    }

    /// Reads one byte at `pos`.
    fn u8_at(
        &self,
        pos: u64,
    ) -> Result<u8, String> {
        let bytes = self.bytes(pos, 1)?;
        bytes
            .first()
            .copied()
            .ok_or_else(|| format!("missing byte at {pos}"))
    }

    /// Reads a little-endian `u16` at `pos`.
    fn u16_at(
        &self,
        pos: u64,
    ) -> Result<u16, String> {
        let b = self.bytes(pos, 2)?;
        Ok(u16::from_le_bytes([
            *b.first().ok_or("short read")?,
            *b.get(1).ok_or("short read")?,
        ]))
    }

    /// Reads a little-endian `u32` at `pos`.
    fn u32_at(
        &self,
        pos: u64,
    ) -> Result<u32, String> {
        let b = self.bytes(pos, 4)?;
        Ok(u32::from_le_bytes([
            *b.first().ok_or("short read")?,
            *b.get(1).ok_or("short read")?,
            *b.get(2).ok_or("short read")?,
            *b.get(3).ok_or("short read")?,
        ]))
    }

    /// Reads a little-endian `u64` at `pos`.
    fn u64_at(
        &self,
        pos: u64,
    ) -> Result<u64, String> {
        let b = self.bytes(pos, 8)?;
        let mut raw = [0_u8; 8];
        raw.copy_from_slice(b);
        Ok(u64::from_le_bytes(raw))
    }

    /// Reads an unsigned integer of `size` (1, 2, 4 or 8) little-endian bytes.
    fn uint_at(
        &self,
        pos: u64,
        size: usize,
    ) -> Result<u64, String> {
        match size {
            | 1 => Ok(u64::from(self.u8_at(pos)?)),
            | 2 => Ok(u64::from(self.u16_at(pos)?)),
            | 4 => Ok(u64::from(self.u32_at(pos)?)),
            | 8 => self.u64_at(pos),
            | _ => Err(format!("unsupported integer width {size}")),
        }
    }

    /// Reads a raw address without applying the base offset.
    fn off_raw(
        &self,
        pos: u64,
    ) -> Result<u64, String> {
        self.uint_at(pos, self.off_size)
    }

    /// Reads an address, mapping the all-ones sentinel to `u64::MAX` and
    /// subtracting the file base address.
    fn off_at(
        &self,
        pos: u64,
    ) -> Result<u64, String> {
        let raw = self.off_raw(pos)?;
        if self.off_size < 8 && raw == (1_u64 << (8 * self.off_size)) - 1 {
            return Ok(u64::MAX);
        }
        if raw == u64::MAX {
            return Ok(u64::MAX);
        }
        raw.checked_sub(self.base)
            .ok_or_else(|| format!("address {raw} precedes the base address"))
    }

    /// Reads a length field.
    fn len_at(
        &self,
        pos: u64,
    ) -> Result<u64, String> {
        self.uint_at(pos, self.len_size)
    }
}

/// Rounds `n` up to the next multiple of eight.
const fn align8(n: u64) -> u64 {
    (n + 7) & !7
}

/// Converts a `u64` count to `usize`, rejecting absurd values.
fn to_usize(
    value: u64,
    what: &str,
) -> Result<usize, String> {
    usize::try_from(value).map_err(|_| format!("{what} {value} does not fit in memory"))
}

/// Multiplies dimension sizes with overflow and budget checks.
fn element_count(shape: &[usize]) -> Result<u64, String> {
    let mut total: u64 = 1;
    for &dim in shape {
        total = total
            .checked_mul(u64::try_from(dim).map_err(|_| "dimension too large".to_owned())?)
            .ok_or_else(|| "dataset dimensions overflow".to_owned())?;
        if total > MAX_ELEMS {
            return Err(format!(
                "dataset has more than {MAX_ELEMS} elements, which is not supported"
            ));
        }
    }
    Ok(total)
}

/// Row-major strides for `shape`.
fn strides_of(shape: &[usize]) -> Vec<usize> {
    let mut strides = vec![0_usize; shape.len()];
    let mut acc = 1_usize;
    for (slot, &dim) in strides.iter_mut().zip(shape.iter()).rev() {
        *slot = acc;
        acc = acc.saturating_mul(dim);
    }
    strides
}

// -------------------------------------------------------------------------
// Superblock
// -------------------------------------------------------------------------

/// Parsed superblock information needed to walk the file.
struct Superblock {
    /// Address of the root group object header.
    root: u64,
    /// Width of an address in bytes.
    off_size: usize,
    /// Width of a length field in bytes.
    len_size: usize,
    /// Base address of the file.
    base: u64,
}

/// The eight-byte HDF5 signature.
const HDF5_SIGNATURE: [u8; 8] = [0x89, b'H', b'D', b'F', 0x0d, 0x0a, 0x1a, 0x0a];

/// Parses the superblock and returns the reader configured for the file.
fn parse_superblock(buf: &[u8]) -> Result<(Reader<'_>, Superblock), String> {
    let signature = buf.get(..8).ok_or("file is too short to be an HDF5 file")?;
    if signature != HDF5_SIGNATURE {
        return Err("not an HDF5 file: signature mismatch".to_owned());
    }
    let version = *buf.get(8).ok_or("truncated superblock")?;
    let (off_size, len_size, base_pos) = match version {
        | 0 | 1 => {
            (
                usize::from(*buf.get(13).ok_or("truncated superblock")?),
                usize::from(*buf.get(14).ok_or("truncated superblock")?),
                24_u64,
            )
        },
        | 2 | 3 => {
            (
                usize::from(*buf.get(9).ok_or("truncated superblock")?),
                usize::from(*buf.get(10).ok_or("truncated superblock")?),
                12_u64,
            )
        },
        | other => return Err(format!("unsupported HDF5 superblock version {other}")),
    };
    if !matches!(off_size, 2 | 4 | 8) {
        return Err(format!("unsupported HDF5 address size {off_size}"));
    }
    if !matches!(len_size, 2 | 4 | 8) {
        return Err(format!("unsupported HDF5 length size {len_size}"));
    }
    let probe = Reader {
        buf,
        off_size,
        len_size,
        base: 0,
    };
    let base = probe.off_raw(base_pos)?;
    let root_field = match version {
        | 0 | 1 => 24 + 4 * (off_size as u64) + off_size as u64,
        | _ => 12 + 3 * (off_size as u64),
    };
    let root_raw = probe.off_raw(root_field)?;
    if root_raw == u64::MAX {
        return Err("superblock has no root group".to_owned());
    }
    let root = root_raw
        .checked_sub(base)
        .ok_or_else(|| "root group address precedes the base address".to_owned())?;
    Ok((
        Reader {
            buf,
            off_size,
            len_size,
            base,
        },
        Superblock {
            root,
            off_size,
            len_size,
            base,
        },
    ))
}

// -------------------------------------------------------------------------
// Object header messages
// -------------------------------------------------------------------------

/// One object header message with its raw payload.
struct Message<'a> {
    /// Message type identifier.
    kind: u16,
    /// Message payload.
    data: &'a [u8],
}

/// Reads every message of the object header at `addr`.
fn read_messages<'a>(
    reader: &Reader<'a>,
    addr: u64,
) -> Result<Vec<Message<'a>>, String> {
    let head = reader.bytes(addr, 4)?;
    let mut messages = Vec::new();
    if head == b"OHDR" {
        read_messages_v2(reader, addr, &mut messages, 0)?;
    } else {
        read_messages_v1(reader, addr, &mut messages, 0)?;
    }
    Ok(messages)
}

/// Parses a version-1 object header.
fn read_messages_v1<'a>(
    reader: &Reader<'a>,
    addr: u64,
    out: &mut Vec<Message<'a>>,
    depth: usize,
) -> Result<(), String> {
    if depth > 8 {
        return Err("object header continuation chain is too deep".to_owned());
    }
    let version = reader.u8_at(addr)?;
    if version != 1 {
        return Err(format!("unsupported object header version {version}"));
    }
    let count = reader.u16_at(addr + 2)?;
    let size = reader.u32_at(addr + 8)?;
    let start = addr + 16;
    let end = start
        .checked_add(u64::from(size))
        .ok_or_else(|| "object header size overflow".to_owned())?;
    let mut pos = start;
    for _ in 0..count {
        if pos + 8 > end {
            return Err("object header message table overruns its block".to_owned());
        }
        let kind = reader.u16_at(pos)?;
        let msize = u64::from(reader.u16_at(pos + 2)?);
        let data_pos = pos + 8;
        let data = reader.bytes(data_pos, msize)?;
        pos = data_pos + align8(msize);
        if kind == 0x10 {
            if msize >= (reader.off_size + reader.len_size) as u64 {
                let next_addr = reader.off_at(data_pos)?;
                let next_len = reader.len_at(data_pos + reader.off_size as u64)?;
                if next_addr != u64::MAX {
                    read_continuation_v1(reader, next_addr, next_len, out, depth + 1)?;
                }
            }
        } else {
            out.push(Message { kind, data });
        }
    }
    Ok(())
}

/// Parses one version-1 continuation block of messages.
fn read_continuation_v1<'a>(
    reader: &Reader<'a>,
    addr: u64,
    len: u64,
    out: &mut Vec<Message<'a>>,
    depth: usize,
) -> Result<(), String> {
    if depth > 8 {
        return Err("object header continuation chain is too deep".to_owned());
    }
    let end = addr
        .checked_add(len)
        .ok_or_else(|| "continuation block overflow".to_owned())?;
    let mut pos = addr;
    while pos + 8 <= end {
        let kind = reader.u16_at(pos)?;
        let msize = u64::from(reader.u16_at(pos + 2)?);
        let data_pos = pos + 8;
        let data = reader.bytes(data_pos, msize)?;
        pos = data_pos + align8(msize);
        if kind == 0x10 {
            if msize >= (reader.off_size + reader.len_size) as u64 {
                let next_addr = reader.off_at(data_pos)?;
                let next_len = reader.len_at(data_pos + reader.off_size as u64)?;
                if next_addr != u64::MAX {
                    read_continuation_v1(reader, next_addr, next_len, out, depth + 1)?;
                }
            }
        } else {
            out.push(Message { kind, data });
        }
    }
    Ok(())
}

/// Parses a version-2 object header (and its `OCHK` continuation blocks).
fn read_messages_v2<'a>(
    reader: &Reader<'a>,
    addr: u64,
    out: &mut Vec<Message<'a>>,
    depth: usize,
) -> Result<(), String> {
    if depth > 8 {
        return Err("object header continuation chain is too deep".to_owned());
    }
    let version = reader.u8_at(addr + 4)?;
    if version != 2 {
        return Err(format!("unsupported object header version {version}"));
    }
    let flags = reader.u8_at(addr + 5)?;
    let mut pos = addr + 6;
    if flags & 0x20 != 0 {
        pos += 16;
    }
    if flags & 0x10 != 0 {
        pos += 4;
    }
    let size_width = 1_usize << (flags & 0x03);
    let chunk_size = reader.uint_at(pos, size_width)?;
    pos += size_width as u64;
    let end = pos
        .checked_add(chunk_size)
        .ok_or_else(|| "object header size overflow".to_owned())?;
    while pos < end {
        if pos + 4 > end {
            return Err("object header message table overruns its block".to_owned());
        }
        let kind = u16::from(reader.u8_at(pos)?);
        let msize = u64::from(reader.u16_at(pos + 1)?);
        let mflags = reader.u8_at(pos + 3)?;
        pos += 4;
        if mflags & 0x04 != 0 {
            pos += 2;
        }
        let data_pos = pos;
        let data = reader.bytes(data_pos, msize)?;
        pos += msize;
        if pos > end {
            return Err("object header message overruns its block".to_owned());
        }
        if kind == 0x10 {
            if msize >= (reader.off_size + reader.len_size) as u64 {
                let next_addr = reader.off_at(data_pos)?;
                let next_len = reader.len_at(data_pos + reader.off_size as u64)?;
                if next_addr != u64::MAX {
                    read_continuation_v2(reader, next_addr, next_len, out, depth + 1)?;
                }
            }
        } else {
            out.push(Message { kind, data });
        }
    }
    Ok(())
}

/// Parses one version-2 continuation block, which starts with `OCHK`.
fn read_continuation_v2<'a>(
    reader: &Reader<'a>,
    addr: u64,
    len: u64,
    out: &mut Vec<Message<'a>>,
    depth: usize,
) -> Result<(), String> {
    if depth > 8 {
        return Err("object header continuation chain is too deep".to_owned());
    }
    if reader.bytes(addr, 4)? != b"OCHK" {
        return Err("invalid object header continuation block signature".to_owned());
    }
    let end = addr
        .checked_add(len)
        .ok_or_else(|| "continuation block overflow".to_owned())?;
    let stop = end.saturating_sub(4);
    let mut pos = addr + 4;
    while pos + 4 <= stop {
        let kind = u16::from(reader.u8_at(pos)?);
        let msize = u64::from(reader.u16_at(pos + 1)?);
        let mflags = reader.u8_at(pos + 3)?;
        pos += 4;
        if mflags & 0x04 != 0 {
            pos += 2;
        }
        let data = reader.bytes(pos, msize)?;
        pos += msize;
        if pos > stop {
            return Err("continuation message overruns its block".to_owned());
        }
        out.push(Message { kind, data });
    }
    Ok(())
}

// -------------------------------------------------------------------------
// Dataspace, datatype, layout and filters
// -------------------------------------------------------------------------

/// Reads an unsigned integer of `size` bytes little-endian from `data`.
fn read_uint(
    data: &[u8],
    pos: usize,
    size: usize,
) -> Result<u64, String> {
    let end = pos
        .checked_add(size)
        .ok_or_else(|| "integer field offset overflow".to_owned())?;
    let bytes = data
        .get(pos..end)
        .ok_or_else(|| "truncated integer field".to_owned())?;
    let mut raw: u64 = 0;
    for (index, &byte) in bytes.iter().enumerate() {
        raw |= u64::from(byte) << (8 * index);
    }
    Ok(raw)
}

/// Parses a dataspace message into its dimension sizes.
fn parse_dataspace(
    data: &[u8],
    len_size: usize,
) -> Result<Vec<usize>, String> {
    let version = *data.first().ok_or("empty dataspace message")?;
    let rank = usize::from(*data.get(1).ok_or("truncated dataspace message")?);
    if rank > 32 {
        return Err(format!("dataspace rank {rank} is too large"));
    }
    let mut pos = match version {
        | 1 => 8_usize,
        | 2 => 4_usize,
        | other => return Err(format!("unsupported dataspace version {other}")),
    };
    let mut shape = Vec::with_capacity(rank);
    for _ in 0..rank {
        let raw = read_uint(data, pos, len_size)?;
        pos += len_size;
        shape.push(to_usize(raw, "dimension")?);
    }
    Ok(shape)
}

/// The datatype classes this module decodes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DType {
    /// Integer of `size` bytes.
    Int {
        /// Width in bytes.
        size: usize,
        /// Whether the value is signed.
        signed: bool,
        /// Whether the bytes are big-endian.
        be: bool,
    },
    /// IEEE floating point of `size` bytes.
    Float {
        /// Width in bytes (4 or 8).
        size: usize,
        /// Whether the bytes are big-endian.
        be: bool,
    },
    /// Fixed-length string of `size` bytes.
    Str {
        /// Width in bytes.
        size: usize,
    },
    /// Any other class, retained only so the caller can report it.
    Other {
        /// HDF5 datatype class number.
        class: u8,
    },
}

impl DType {
    /// Byte width of one element.
    const fn size(self) -> usize {
        match self {
            | Self::Int { size, .. } | Self::Float { size, .. } | Self::Str { size } => size,
            | Self::Other { .. } => 0,
        }
    }
}

/// Parses a datatype message.
fn parse_datatype(data: &[u8]) -> Result<DType, String> {
    let class_and_version = *data.first().ok_or("empty datatype message")?;
    let class = class_and_version & 0x0f;
    let bits = data.get(1..4).ok_or("truncated datatype message")?;
    let size = to_usize(read_uint(data, 4, 4)?, "datatype size")?;
    match class {
        | 0 => {
            let flags = *bits.first().ok_or("truncated datatype message")?;
            let be = flags & 0x01 != 0;
            let signed = flags & 0x08 != 0;
            if !matches!(size, 1 | 2 | 4 | 8) {
                return Err(format!("unsupported integer width {size}"));
            }
            Ok(DType::Int { size, signed, be })
        },
        | 1 => {
            let flags = *bits.first().ok_or("truncated datatype message")?;
            let order = flags & 0x03;
            let be = match order {
                | 0 => false,
                | 1 => true,
                | _ => return Err("VAX byte order is not supported".to_owned()),
            };
            if !matches!(size, 4 | 8) {
                return Err(format!("unsupported floating-point width {size}"));
            }
            Ok(DType::Float { size, be })
        },
        | 3 => Ok(DType::Str { size }),
        | other => Ok(DType::Other { class: other }),
    }
}

/// The storage layout of a dataset.
enum Layout {
    /// Bytes stored inside the object header.
    Compact(Vec<u8>),
    /// One contiguous block at `addr`; a zero `size` means "derive it from the
    /// dataspace".
    Contiguous {
        /// File address of the data.
        addr: u64,
        /// Byte length of the block, or zero when unrecorded.
        size: u64,
    },
    /// A chunked dataset indexed by a version-1 B-tree.
    Chunked {
        /// Address of the chunk B-tree.
        addr: u64,
        /// Chunk dimensions followed by the element size.
        dims: Vec<usize>,
    },
}

/// Parses a data layout message.
fn parse_layout(
    reader: &Reader<'_>,
    data: &[u8],
) -> Result<Layout, String> {
    let version = *data.first().ok_or("empty data layout message")?;
    match version {
        | 1 | 2 => {
            let ndims = usize::from(*data.get(1).ok_or("truncated data layout message")?);
            let class = *data.get(2).ok_or("truncated data layout message")?;
            let mut pos = 8_usize;
            let addr = if matches!(class, 1 | 2) {
                let value = read_uint(data, pos, reader.off_size)?;
                pos += reader.off_size;
                value
            } else {
                u64::MAX
            };
            let mut dims = Vec::with_capacity(ndims);
            for _ in 0..ndims {
                dims.push(to_usize(read_uint(data, pos, 4)?, "chunk dimension")?);
                pos += 4;
            }
            match class {
                | 0 => Ok(Layout::Compact(data.get(pos..).unwrap_or(&[]).to_vec())),
                | 1 => Ok(Layout::Contiguous { addr, size: 0 }),
                | 2 => Ok(Layout::Chunked { addr, dims }),
                | other => Err(format!("unsupported data layout class {other}")),
            }
        },
        | 3 => {
            let class = *data.get(1).ok_or("truncated data layout message")?;
            match class {
                | 0 => {
                    let size = to_usize(read_uint(data, 2, 2)?, "compact size")?;
                    let end = 4_usize
                        .checked_add(size)
                        .ok_or_else(|| "compact data size overflow".to_owned())?;
                    let body = data.get(4..end).ok_or("truncated compact data")?.to_vec();
                    Ok(Layout::Compact(body))
                },
                | 1 => {
                    let addr = read_uint(data, 2, reader.off_size)?;
                    let size = read_uint(data, 2 + reader.off_size, reader.len_size)?;
                    Ok(Layout::Contiguous { addr, size })
                },
                | 2 => {
                    let ndims = usize::from(*data.get(2).ok_or("truncated chunked layout")?);
                    let addr = read_uint(data, 3, reader.off_size)?;
                    let mut pos = 3 + reader.off_size;
                    let mut dims = Vec::with_capacity(ndims);
                    for _ in 0..ndims {
                        dims.push(to_usize(read_uint(data, pos, 4)?, "chunk dimension")?);
                        pos += 4;
                    }
                    if dims.is_empty() {
                        return Err("chunked layout has no dimensions".to_owned());
                    }
                    Ok(Layout::Chunked { addr, dims })
                },
                | other => Err(format!("unsupported data layout class {other}")),
            }
        },
        | 4 => {
            let class = *data.get(1).ok_or("truncated data layout message")?;
            match class {
                | 0 => {
                    let size = to_usize(read_uint(data, 2, 2)?, "compact size")?;
                    let end = 4_usize
                        .checked_add(size)
                        .ok_or_else(|| "compact data size overflow".to_owned())?;
                    let body = data.get(4..end).ok_or("truncated compact data")?.to_vec();
                    Ok(Layout::Compact(body))
                },
                | 1 => {
                    let addr = read_uint(data, 2, reader.off_size)?;
                    let size = read_uint(data, 2 + reader.off_size, reader.len_size)?;
                    Ok(Layout::Contiguous { addr, size })
                },
                | _ => Err("version-4 chunked layouts are not supported".to_owned()),
            }
        },
        | other => Err(format!("unsupported data layout version {other}")),
    }
}

/// One entry of a filter pipeline.
#[derive(Debug, Clone, PartialEq, Eq)]
struct Filter {
    /// HDF5 filter identifier.
    id: u16,
    /// Client data values for the filter.
    client: Vec<u32>,
}

/// Parses a filter pipeline message.
fn parse_filters(data: &[u8]) -> Result<Vec<Filter>, String> {
    let version = *data.first().ok_or("empty filter pipeline message")?;
    let count = usize::from(*data.get(1).ok_or("truncated filter pipeline")?);
    let mut filters = Vec::with_capacity(count);
    match version {
        | 1 => {
            let mut pos = 8_usize;
            for _ in 0..count {
                let id = u16::from_le_bytes([
                    *data.get(pos).ok_or("truncated filter")?,
                    *data.get(pos + 1).ok_or("truncated filter")?,
                ]);
                let name_len = to_usize(read_uint(data, pos + 2, 2)?, "filter name length")?;
                let nclient = to_usize(read_uint(data, pos + 6, 2)?, "filter client count")?;
                pos += 8;
                pos += to_usize(align8(u64::try_from(name_len).unwrap_or(0)), "padding")?;
                let mut client = Vec::with_capacity(nclient);
                for index in 0..nclient {
                    client.push(
                        u32::try_from(read_uint(data, pos + 4 * index, 4)?)
                            .map_err(|_| "filter client value overflow".to_owned())?,
                    );
                }
                pos += 4 * nclient;
                if nclient % 2 == 1 {
                    pos += 4;
                }
                if pos > data.len() {
                    return Err("filter pipeline overruns its message".to_owned());
                }
                filters.push(Filter { id, client });
            }
        },
        | 2 => {
            let mut pos = 2_usize;
            for _ in 0..count {
                let id = u16::from_le_bytes([
                    *data.get(pos).ok_or("truncated filter")?,
                    *data.get(pos + 1).ok_or("truncated filter")?,
                ]);
                pos += 2;
                if id >= 256 {
                    pos += 2;
                }
                let nclient = to_usize(read_uint(data, pos + 2, 2)?, "filter client count")?;
                pos += 4;
                let mut client = Vec::with_capacity(nclient);
                for index in 0..nclient {
                    client.push(
                        u32::try_from(read_uint(data, pos + 4 * index, 4)?)
                            .map_err(|_| "filter client value overflow".to_owned())?,
                    );
                }
                pos += 4 * nclient;
                if pos > data.len() {
                    return Err("filter pipeline overruns its message".to_owned());
                }
                filters.push(Filter { id, client });
            }
        },
        | other => return Err(format!("unsupported filter pipeline version {other}")),
    }
    Ok(filters)
}

// -------------------------------------------------------------------------
// Element conversion
// -------------------------------------------------------------------------

/// Converts `count` raw elements of `dtype` into `f64` values.
fn decode_elements(
    dtype: DType,
    bytes: &[u8],
    count: usize,
) -> Result<Vec<f64>, String> {
    let size = dtype.size();
    if size == 0 {
        return Err("dataset uses an unsupported datatype class".to_owned());
    }
    let needed = count
        .checked_mul(size)
        .ok_or_else(|| "dataset byte size overflow".to_owned())?;
    if needed > MAX_ARRAY_BYTES {
        return Err(format!(
            "dataset needs {needed} bytes, more than the {MAX_ARRAY_BYTES}-byte limit"
        ));
    }
    let available = bytes
        .get(..needed)
        .ok_or_else(|| "dataset data block is shorter than its dataspace".to_owned())?;
    let mut out = Vec::with_capacity(count);
    for chunk in available.chunks_exact(size) {
        out.push(decode_one(dtype, chunk)?);
    }
    Ok(out)
}

/// Converts a single raw element into an `f64`.
fn decode_one(
    dtype: DType,
    bytes: &[u8],
) -> Result<f64, String> {
    match dtype {
        | DType::Int { size, signed, be } => {
            let mut raw: u64 = 0;
            if be {
                for &byte in bytes.iter().take(size) {
                    raw = (raw << 8) | u64::from(byte);
                }
            } else {
                for (index, &byte) in bytes.iter().take(size).enumerate() {
                    raw |= u64::from(byte) << (8 * index);
                }
            }
            let value = if signed {
                match size {
                    | 1 => f64::from(raw as u8 as i8),
                    | 2 => f64::from(raw as u16 as i16),
                    | 4 => f64::from(raw as u32 as i32),
                    | _ => raw as i64 as f64,
                }
            } else {
                match size {
                    | 1 => f64::from(raw as u8),
                    | 2 => f64::from(raw as u16),
                    | 4 => f64::from(raw as u32),
                    | _ => raw as f64,
                }
            };
            Ok(value)
        },
        | DType::Float { size, be } => {
            if size == 4 {
                let mut raw = [0_u8; 4];
                raw.copy_from_slice(bytes.get(..4).ok_or("short float")?);
                let bits = if be {
                    u32::from_be_bytes(raw)
                } else {
                    u32::from_le_bytes(raw)
                };
                Ok(f64::from(f32::from_bits(bits)))
            } else {
                let mut raw = [0_u8; 8];
                raw.copy_from_slice(bytes.get(..8).ok_or("short float")?);
                let bits = if be {
                    u64::from_be_bytes(raw)
                } else {
                    u64::from_le_bytes(raw)
                };
                Ok(f64::from_bits(bits))
            }
        },
        | DType::Str { .. } => Err("string datasets are not supported".to_owned()),
        | DType::Other { class } => Err(format!("unsupported HDF5 datatype class {class}")),
    }
}

/// Maps a datatype to the public [`H5Type`] enumeration.
fn h5_type_of(dtype: DType) -> Result<H5Type, String> {
    match dtype {
        | DType::Int { size, signed, .. } => {
            match (size, signed) {
                | (1, true) => Ok(H5Type::I8),
                | (1, false) => Ok(H5Type::U8),
                | (2, true) => Ok(H5Type::I16),
                | (2, false) => Ok(H5Type::U16),
                | (4, true) => Ok(H5Type::I32),
                | (4, false) => Ok(H5Type::U32),
                | (8, true) => Ok(H5Type::I64),
                | (8, false) => Ok(H5Type::U64),
                | _ => Err(format!("unsupported integer width {size}")),
            }
        },
        | DType::Float { size, .. } => {
            match size {
                | 4 => Ok(H5Type::F32),
                | 8 => Ok(H5Type::F64),
                | _ => Err(format!("unsupported floating-point width {size}")),
            }
        },
        | DType::Str { .. } => Ok(H5Type::Str),
        | DType::Other { class } => Err(format!("unsupported HDF5 datatype class {class}")),
    }
}

// -------------------------------------------------------------------------
// Groups
// -------------------------------------------------------------------------

/// Reads a local heap and returns its data segment address and size.
fn read_local_heap(
    reader: &Reader<'_>,
    addr: u64,
) -> Result<(u64, u64), String> {
    if reader.bytes(addr, 4)? != b"HEAP" {
        return Err("invalid local heap signature".to_owned());
    }
    let size = reader.len_at(addr + 8)?;
    let data_addr = reader.off_at(addr + 8 + 2 * reader.len_size as u64)?;
    Ok((data_addr, size))
}

/// Reads a NUL-terminated name from the heap data segment.
fn heap_name(
    reader: &Reader<'_>,
    data_addr: u64,
    data_size: u64,
    offset: u64,
) -> Result<String, String> {
    if offset >= data_size {
        return Err("link name offset is outside the local heap".to_owned());
    }
    let available = data_size - offset;
    let span = available.min(4096);
    let bytes = reader.bytes(data_addr + offset, span)?;
    let end = bytes.iter().position(|&b| b == 0).unwrap_or(bytes.len());
    let name = bytes
        .get(..end)
        .ok_or_else(|| "invalid link name".to_owned())?;
    std::str::from_utf8(name)
        .map(str::to_owned)
        .map_err(|_| "link name is not valid UTF-8".to_owned())
}

/// Reads one symbol table node (`SNOD`) and appends its entries.
fn read_symbol_table_node(
    reader: &Reader<'_>,
    addr: u64,
    heap: (u64, u64),
    out: &mut Vec<(String, u64)>,
) -> Result<(), String> {
    if reader.bytes(addr, 4)? != b"SNOD" {
        return Err("invalid symbol table node signature".to_owned());
    }
    let count = usize::from(reader.u16_at(addr + 6)?);
    let entry_size = 2 * reader.off_size + 24;
    let base = addr + 8;
    for index in 0..count {
        let offset = base + (index * entry_size) as u64;
        let name_offset = reader.len_at(offset)?;
        let object = reader.off_at(offset + reader.off_size as u64)?;
        if object == u64::MAX {
            continue;
        }
        let name = heap_name(reader, heap.0, heap.1, name_offset)?;
        out.push((name, object));
    }
    Ok(())
}

/// Recursively walks a version-1 group B-tree.
fn walk_group_btree(
    reader: &Reader<'_>,
    addr: u64,
    heap: (u64, u64),
    out: &mut Vec<(String, u64)>,
    depth: usize,
) -> Result<(), String> {
    if depth > 8 {
        return Err("group B-tree is too deep".to_owned());
    }
    let head = reader.bytes(addr, 4)?;
    if head == b"BTHD" {
        return Err("fractal-heap groups are not supported".to_owned());
    }
    if head != b"TREE" {
        return Err("invalid group B-tree signature".to_owned());
    }
    let node_type = reader.u8_at(addr + 4)?;
    if node_type != 0 {
        return Err(format!("unexpected B-tree node type {node_type}"));
    }
    let level = reader.u8_at(addr + 5)?;
    let count = usize::from(reader.u16_at(addr + 6)?);
    if count > MAX_BTREE_ENTRIES {
        return Err("group B-tree node has too many entries".to_owned());
    }
    let key_size = reader.len_size;
    let mut pos = addr + 24;
    for _ in 0..count {
        let child = reader.off_at(pos + key_size as u64)?;
        if level == 0 {
            if child != u64::MAX {
                read_symbol_table_node(reader, child, heap, out)?;
            }
        } else if child != u64::MAX {
            walk_group_btree(reader, child, heap, out, depth + 1)?;
        }
        pos += key_size as u64 + reader.off_size as u64;
    }
    Ok(())
}

/// Parses a link message and returns the target when it is a hard link.
fn parse_link_message(
    reader: &Reader<'_>,
    data: &[u8],
) -> Result<Option<(String, u64)>, String> {
    let version = *data.first().ok_or("empty link message")?;
    if version != 1 {
        return Err(format!("unsupported link message version {version}"));
    }
    let flags = *data.get(1).ok_or("truncated link message")?;
    let mut pos = 2_usize;
    let link_type = if flags & 0x08 != 0 {
        let value = *data.get(pos).ok_or("truncated link message")?;
        pos += 1;
        value
    } else {
        0
    };
    if flags & 0x04 != 0 {
        pos += 8;
    }
    if flags & 0x10 != 0 {
        pos += 1;
    }
    let width = 1_usize << (flags & 0x03);
    let name_len = to_usize(read_uint(data, pos, width)?, "link name length")?;
    pos += width;
    let name_end = pos
        .checked_add(name_len)
        .ok_or_else(|| "link name length overflow".to_owned())?;
    let name_bytes = data.get(pos..name_end).ok_or("truncated link name")?;
    pos = name_end;
    let name = std::str::from_utf8(name_bytes)
        .map_err(|_| "link name is not valid UTF-8".to_owned())?
        .to_owned();
    if link_type != 0 {
        return Ok(None);
    }
    let addr = read_uint(data, pos, reader.off_size)?;
    if addr == u64::MAX {
        return Ok(None);
    }
    Ok(Some((name, addr)))
}

/// Collects the hard links declared by a group object header.
fn collect_links(
    reader: &Reader<'_>,
    messages: &[Message<'_>],
) -> Result<Vec<(String, u64)>, String> {
    let mut links = Vec::new();
    let mut symbol_table: Option<(u64, u64)> = None;
    for message in messages {
        match message.kind {
            | 0x02 => {
                let flags = *message.data.get(1).ok_or("truncated link info message")?;
                if flags & 0x01 != 0 {
                    return Err("fractal-heap groups are not supported".to_owned());
                }
            },
            | 0x11 => {
                let data = message.data;
                let btree = read_uint(data, 0, reader.off_size)?;
                let heap = read_uint(data, reader.off_size, reader.off_size)?;
                if btree != u64::MAX && heap != u64::MAX {
                    symbol_table = Some((btree, heap));
                }
            },
            | 0x06 => {
                if let Some(link) = parse_link_message(reader, message.data)? {
                    links.push(link);
                }
            },
            | _ => {},
        }
    }
    if let Some((btree, heap_addr)) = symbol_table {
        let heap = read_local_heap(reader, heap_addr)?;
        walk_group_btree(reader, btree, heap, &mut links, 0)?;
    }
    Ok(links)
}

// -------------------------------------------------------------------------
// Attributes
// -------------------------------------------------------------------------

/// Parses one attribute message.
fn parse_attribute(
    reader: &Reader<'_>,
    data: &[u8],
) -> Result<Attribute, String> {
    let version = *data.first().ok_or("empty attribute message")?;
    let (mut pos, padded) = match version {
        | 1 => (8_usize, true),
        | 2 => (8_usize, false),
        | 3 => {
            let flags = *data.get(1).ok_or("truncated attribute message")?;
            (if flags & 0x01 != 0 { 9 } else { 8 }, false)
        },
        | other => return Err(format!("unsupported attribute message version {other}")),
    };
    let name_len = to_usize(read_uint(data, 2, 2)?, "attribute name length")?;
    let dtype_len = to_usize(read_uint(data, 4, 2)?, "attribute datatype length")?;
    let dspace_len = to_usize(read_uint(data, 6, 2)?, "attribute dataspace length")?;
    let name_end = pos
        .checked_add(name_len)
        .ok_or_else(|| "attribute name length overflow".to_owned())?;
    let name_bytes = data.get(pos..name_end).ok_or("truncated attribute name")?;
    let terminator = name_bytes
        .iter()
        .position(|&b| b == 0)
        .unwrap_or(name_bytes.len());
    let name = std::str::from_utf8(name_bytes.get(..terminator).unwrap_or(&[]))
        .map_err(|_| "attribute name is not valid UTF-8".to_owned())?
        .to_owned();
    pos = advance(data, name_end, name_len, padded)?;
    let dtype_end = pos
        .checked_add(dtype_len)
        .ok_or_else(|| "attribute datatype length overflow".to_owned())?;
    let dtype_bytes = data
        .get(pos..dtype_end)
        .ok_or("truncated attribute datatype")?;
    pos = advance(data, dtype_end, dtype_len, padded)?;
    let dspace_end = pos
        .checked_add(dspace_len)
        .ok_or_else(|| "attribute dataspace length overflow".to_owned())?;
    let dspace_bytes = data
        .get(pos..dspace_end)
        .ok_or("truncated attribute dataspace")?;
    pos = advance(data, dspace_end, dspace_len, padded)?;
    let dtype = parse_datatype(dtype_bytes)?;
    let shape = parse_dataspace(dspace_bytes, reader.len_size)?;
    let count = element_count(&shape)?;
    let body = data.get(pos..).ok_or("truncated attribute data")?;
    let value = decode_attribute_value(dtype, body, count)?;
    Ok(Attribute { name, value })
}

/// Advances past a field, adding the 8-byte padding used by version-1 and
/// version-2 attribute messages.
fn advance(
    data: &[u8],
    end: usize,
    len: usize,
    padded: bool,
) -> Result<usize, String> {
    if padded {
        let pad = to_usize(align8(u64::try_from(len).unwrap_or(0)), "padding")?;
        let next = end
            .checked_add(pad - len)
            .ok_or_else(|| "attribute field overflow".to_owned())?;
        if next > data.len() {
            return Err("attribute message is truncated".to_owned());
        }
        Ok(next)
    } else {
        Ok(end)
    }
}

/// Converts raw attribute bytes into an [`AttrValue`].
fn decode_attribute_value(
    dtype: DType,
    body: &[u8],
    count: u64,
) -> Result<AttrValue, String> {
    let count = to_usize(count, "attribute element count")?;
    if let DType::Str { size } = dtype {
        let needed = count
            .checked_mul(size)
            .ok_or_else(|| "attribute size overflow".to_owned())?;
        let bytes = body
            .get(..needed)
            .ok_or_else(|| "attribute data block is too short".to_owned())?;
        let end = bytes.iter().position(|&b| b == 0).unwrap_or(bytes.len());
        let text = std::str::from_utf8(bytes.get(..end).unwrap_or(&[]))
            .map_err(|_| "attribute text is not valid UTF-8".to_owned())?;
        return Ok(AttrValue::Text(text.to_owned()));
    }
    if count == 0 {
        return Ok(AttrValue::Floats(Vec::new()));
    }
    let values = decode_elements(dtype, body, count)?;
    if count == 1 {
        let value = values.first().copied().ok_or("missing attribute value")?;
        return Ok(match dtype {
            | DType::Int { .. } => AttrValue::Int(value as i64),
            | _ => AttrValue::Float(value),
        });
    }
    Ok(match dtype {
        | DType::Int { .. } => AttrValue::Ints(values.iter().map(|&v| v as i64).collect()),
        | _ => AttrValue::Floats(values),
    })
}

// -------------------------------------------------------------------------
// Chunked storage
// -------------------------------------------------------------------------

/// One chunk record from a version-1 chunk B-tree.
struct Chunk {
    /// Element offsets of the chunk inside the dataset.
    offsets: Vec<u64>,
    /// File address of the filtered chunk bytes.
    addr: u64,
    /// Byte length of the filtered chunk.
    size: u64,
    /// Filter mask; a set bit means the filter was skipped.
    mask: u32,
}

/// Recursively walks a version-1 chunk B-tree.
fn walk_chunk_btree(
    reader: &Reader<'_>,
    addr: u64,
    ndims: usize,
    out: &mut Vec<Chunk>,
    depth: usize,
) -> Result<(), String> {
    if depth > 8 {
        return Err("chunk B-tree is too deep".to_owned());
    }
    let head = reader.bytes(addr, 4)?;
    if head == b"BTHD" {
        return Err("version-2 chunk B-trees are not supported".to_owned());
    }
    if head != b"TREE" {
        return Err("invalid chunk B-tree signature".to_owned());
    }
    let node_type = reader.u8_at(addr + 4)?;
    if node_type != 1 {
        return Err(format!("unexpected chunk B-tree node type {node_type}"));
    }
    let level = reader.u8_at(addr + 5)?;
    let count = usize::from(reader.u16_at(addr + 6)?);
    if count > MAX_BTREE_ENTRIES {
        return Err("chunk B-tree node has too many entries".to_owned());
    }
    let key_size = 8 + 4 * ndims;
    let mut pos = addr + 24;
    for _ in 0..count {
        let key = reader.bytes(pos, key_size as u64)?;
        let child = reader.off_at(pos + key_size as u64)?;
        if level == 0 {
            if child != u64::MAX {
                let size = read_uint(key, 0, 4)?;
                let mask = u32::try_from(read_uint(key, 4, 4)?)
                    .map_err(|_| "filter mask overflow".to_owned())?;
                let mut offsets = Vec::with_capacity(ndims);
                for index in 0..ndims {
                    offsets.push(read_uint(key, 8 + 4 * index, 4)?);
                }
                out.push(Chunk {
                    offsets,
                    addr: child,
                    size,
                    mask,
                });
            }
        } else if child != u64::MAX {
            walk_chunk_btree(reader, child, ndims, out, depth + 1)?;
        }
        pos += key_size as u64 + reader.off_size as u64;
    }
    Ok(())
}

/// Reverses the `shuffle` filter.
fn unshuffle(
    data: &[u8],
    element: usize,
) -> Result<Vec<u8>, String> {
    if element == 0 {
        return Err("shuffle filter has a zero element size".to_owned());
    }
    let count = data.len() / element;
    let mut out = vec![0_u8; data.len()];
    for index in 0..count {
        for byte in 0..element {
            let src = byte * count + index;
            let dst = index * element + byte;
            let value = *data.get(src).ok_or("shuffle filter input is truncated")?;
            if let Some(slot) = out.get_mut(dst) {
                *slot = value;
            }
        }
    }
    let tail = count * element;
    if tail < data.len() {
        out.get_mut(tail..)
            .ok_or("shuffle output overflow")?
            .copy_from_slice(data.get(tail..).ok_or("shuffle input truncated")?);
    }
    Ok(out)
}

/// Applies the reverse of a filter pipeline to `data`.
fn apply_filters(
    data: Vec<u8>,
    filters: &[Filter],
    mask: u32,
) -> Result<Vec<u8>, String> {
    let mut current = data;
    for (index, filter) in filters.iter().enumerate().rev() {
        if index < 32 && mask & (1_u32 << index) != 0 {
            continue;
        }
        current = match filter.id {
            | 1 => inflate::zlib_decompress(&current)?,
            | 2 => {
                let element = filter.client.first().copied().map_or(1, |v| v as usize);
                unshuffle(&current, element)?
            },
            | 3 => {
                let keep = current.len().saturating_sub(4);
                current.get(..keep).unwrap_or(&[]).to_vec()
            },
            | 4 => return Err("the szip filter is not supported".to_owned()),
            | 5 => return Err("the nbit filter is not supported".to_owned()),
            | 6 => return Err("the scaleoffset filter is not supported".to_owned()),
            | other => return Err(format!("unsupported HDF5 filter {other}")),
        };
    }
    Ok(current)
}

/// Reads and assembles a chunked dataset.
fn read_chunked(
    reader: &Reader<'_>,
    addr: u64,
    shape: &[usize],
    chunk_dims: &[usize],
    dtype: DType,
    filters: &[Filter],
) -> Result<Vec<f64>, String> {
    let total = element_count(shape)?;
    let total = to_usize(total, "dataset element count")?;
    let elem_size = dtype.size();
    if elem_size == 0 {
        return Err("chunked dataset uses an unsupported datatype".to_owned());
    }
    let rank = shape.len();
    if chunk_dims.len() < rank {
        return Err("chunked layout has fewer dimensions than the dataspace".to_owned());
    }
    let chunk_shape = chunk_dims
        .get(..rank)
        .ok_or("chunk shape is too short")?
        .to_vec();
    let ndims = chunk_dims.len();
    let mut chunks = Vec::new();
    walk_chunk_btree(reader, addr, ndims, &mut chunks, 0)?;
    let strides = strides_of(shape);
    let mut out = vec![0.0_f64; total];
    let chunk_elems = element_count(&chunk_shape)?;
    let chunk_elems = to_usize(chunk_elems, "chunk element count")?;
    for chunk in chunks {
        let raw = reader.bytes(chunk.addr, chunk.size)?;
        let filtered = apply_filters(raw.to_vec(), filters, chunk.mask)?;
        let values = decode_elements(dtype, &filtered, chunk_elems)?;
        for (flat, &value) in values.iter().enumerate() {
            let mut remainder = flat;
            let mut inside = true;
            let mut dst = 0_usize;
            for axis in (0..rank).rev() {
                let dim = chunk_shape.get(axis).copied().unwrap_or(0);
                if dim == 0 {
                    inside = false;
                    break;
                }
                let coordinate = remainder % dim;
                remainder /= dim;
                let offset = chunk.offsets.get(axis).copied().unwrap_or(0);
                let position = usize::try_from(offset)
                    .unwrap_or(usize::MAX)
                    .saturating_add(coordinate);
                if position >= shape.get(axis).copied().unwrap_or(0) {
                    inside = false;
                    break;
                }
                dst += position * strides.get(axis).copied().unwrap_or(0);
            }
            if inside {
                if let Some(slot) = out.get_mut(dst) {
                    *slot = value;
                }
            }
        }
    }
    Ok(out)
}

// -------------------------------------------------------------------------
// Datasets
// -------------------------------------------------------------------------

/// Builds a [`Dataset`] from an object header's messages.
fn read_dataset(
    reader: &Reader<'_>,
    messages: &[Message<'_>],
    path: &str,
) -> Result<Dataset, String> {
    let mut shape = Vec::new();
    let mut dtype = None;
    let mut layout = None;
    let mut filters = Vec::new();
    for message in messages {
        match message.kind {
            | 0x01 => shape = parse_dataspace(message.data, reader.len_size)?,
            | 0x03 => dtype = Some(parse_datatype(message.data)?),
            | 0x08 => layout = Some(parse_layout(reader, message.data)?),
            | 0x0B => filters = parse_filters(message.data)?,
            | _ => {},
        }
    }
    let dtype = dtype.ok_or_else(|| format!("dataset {path} has no datatype message"))?;
    let layout = layout.ok_or_else(|| format!("dataset {path} has no data layout message"))?;
    let h5_type = h5_type_of(dtype)?;
    if h5_type == H5Type::Str {
        return Err(format!(
            "dataset {path} holds strings, which are not supported"
        ));
    }
    let count = element_count(&shape)?;
    let count = to_usize(count, "dataset element count")?;
    let data = match layout {
        | Layout::Compact(body) => decode_elements(dtype, &body, count)?,
        | Layout::Contiguous { addr, size } => {
            if addr == u64::MAX {
                vec![0.0_f64; count]
            } else {
                let derived = u64::try_from(count)
                    .unwrap_or(0)
                    .saturating_mul(u64::try_from(dtype.size()).unwrap_or(0));
                let needed = size.max(derived);
                let body = reader.bytes(addr, needed)?;
                decode_elements(dtype, body, count)?
            }
        },
        | Layout::Chunked { addr, dims } => {
            if addr == u64::MAX {
                vec![0.0_f64; count]
            } else {
                read_chunked(reader, addr, &shape, &dims, dtype, &filters)?
            }
        },
    };
    Ok(Dataset {
        path: path.to_owned(),
        shape,
        h5_type,
        data,
    })
}

/// Recursively walks one object, appending groups, datasets and attributes.
fn walk_object(
    reader: &Reader<'_>,
    addr: u64,
    path: &str,
    out: &mut H5File,
    visited: &mut Vec<u64>,
    depth: usize,
) -> Result<(), String> {
    if depth > MAX_DEPTH {
        return Err(format!(
            "group nesting at {path} is deeper than {MAX_DEPTH}"
        ));
    }
    if visited.contains(&addr) {
        return Ok(());
    }
    visited.push(addr);
    let messages = read_messages(reader, addr)?;
    let mut attributes = Vec::new();
    for message in &messages {
        if message.kind == 0x0C {
            attributes.push(parse_attribute(reader, message.data)?);
        }
    }
    if !attributes.is_empty() {
        out.attributes.push((path.to_owned(), attributes));
    }
    let is_dataset = messages.iter().any(|message| message.kind == 0x01)
        && messages.iter().any(|message| message.kind == 0x03)
        && messages.iter().any(|message| message.kind == 0x08);
    if is_dataset {
        out.datasets.push(read_dataset(reader, &messages, path)?);
        return Ok(());
    }
    out.groups.push(path.to_owned());
    let mut links = collect_links(reader, &messages)?;
    links.sort_by(|left, right| left.0.cmp(&right.0));
    for (name, child) in links {
        let child_path = if path == "/" {
            format!("/{name}")
        } else {
            format!("{path}/{name}")
        };
        walk_object(reader, child, &child_path, out, visited, depth + 1)?;
    }
    Ok(())
}

// -------------------------------------------------------------------------
// Public reading API
// -------------------------------------------------------------------------

/// Reads an HDF5 file from disk.
///
/// # Errors
/// Fails when the file cannot be read, is not HDF5, uses an unsupported
/// feature, or contains inconsistent offsets.
pub fn read_hdf5(path: &Path) -> Result<H5File, String> {
    let bytes = std::fs::read(path)
        .map_err(|error| format!("failed to read {}: {error}", path.display()))?;
    read_hdf5_bytes(&bytes)
}

/// Reads an HDF5 file from a byte buffer.
///
/// # Errors
/// Fails when the buffer is not a complete HDF5 file, uses an unsupported
/// feature (fractal-heap groups, version-2 B-trees, exotic filters), or
/// contains inconsistent offsets.
pub fn read_hdf5_bytes(bytes: &[u8]) -> Result<H5File, String> {
    let (reader, superblock) = parse_superblock(bytes)?;
    let eof = superblock_end(bytes, &superblock)?;
    if eof != 0 && eof > u64::try_from(bytes.len()).unwrap_or(u64::MAX) {
        return Err(format!(
            "file is truncated: superblock declares {eof} bytes but only {} are present",
            bytes.len()
        ));
    }
    let mut file = H5File::default();
    let mut visited = Vec::new();
    walk_object(&reader, superblock.root, "/", &mut file, &mut visited, 0)?;
    Ok(file)
}

/// Returns the end-of-file address declared by the superblock, or zero when
/// the field is absent.
fn superblock_end(
    bytes: &[u8],
    superblock: &Superblock,
) -> Result<u64, String> {
    let version = *bytes.get(8).ok_or("truncated superblock")?;
    let probe = Reader {
        buf: bytes,
        off_size: superblock.off_size,
        len_size: superblock.len_size,
        base: 0,
    };
    let raw = match version {
        | 0 | 1 => probe.off_raw(24 + 2 * superblock.off_size as u64)?,
        | _ => probe.off_raw(12 + 2 * superblock.off_size as u64)?,
    };
    if raw == u64::MAX {
        return Ok(0);
    }
    Ok(raw.saturating_sub(superblock.base))
}

// -------------------------------------------------------------------------
// Writer
// -------------------------------------------------------------------------

/// Byte width of an [`H5Type`] element.
const fn type_size(h5_type: H5Type) -> usize {
    match h5_type {
        | H5Type::I8 | H5Type::U8 | H5Type::Str => 1,
        | H5Type::I16 | H5Type::U16 => 2,
        | H5Type::I32 | H5Type::U32 | H5Type::F32 => 4,
        | H5Type::I64 | H5Type::U64 | H5Type::F64 => 8,
    }
}

/// Builds a datatype message for `h5_type`.
fn datatype_message(h5_type: H5Type) -> Result<Vec<u8>, String> {
    let mut out = Vec::new();
    match h5_type {
        | H5Type::I8
        | H5Type::U8
        | H5Type::I16
        | H5Type::U16
        | H5Type::I32
        | H5Type::U32
        | H5Type::I64
        | H5Type::U64 => {
            let signed = matches!(
                h5_type,
                H5Type::I8 | H5Type::I16 | H5Type::I32 | H5Type::I64
            );
            let size = type_size(h5_type);
            out.push(0x10);
            out.push(if signed { 0x08 } else { 0x00 });
            out.extend_from_slice(&[0, 0]);
            out.extend_from_slice(&u32::try_from(size).unwrap_or(0).to_le_bytes());
            out.extend_from_slice(&0_u16.to_le_bytes());
            out.extend_from_slice(&u16::try_from(size * 8).unwrap_or(0).to_le_bytes());
        },
        | H5Type::F32 | H5Type::F64 => {
            let size = type_size(h5_type);
            out.push(0x11);
            out.push(0x20);
            out.extend_from_slice(&[0, 0]);
            out.extend_from_slice(&u32::try_from(size).unwrap_or(0).to_le_bytes());
            out.extend_from_slice(&0_u16.to_le_bytes());
            out.extend_from_slice(&u16::try_from(size * 8).unwrap_or(0).to_le_bytes());
            if size == 4 {
                out.extend_from_slice(&[23, 8, 0, 23]);
                out.extend_from_slice(&127_u32.to_le_bytes());
            } else {
                out.extend_from_slice(&[52, 11, 0, 52]);
                out.extend_from_slice(&1023_u32.to_le_bytes());
            }
        },
        | H5Type::Str => return Err("string datasets cannot be written".to_owned()),
    }
    Ok(out)
}

/// Builds a datatype message for a fixed-length string of `len` bytes.
fn string_datatype_message(len: usize) -> Vec<u8> {
    let mut out = Vec::new();
    out.push(0x13);
    out.extend_from_slice(&[0, 0, 0]);
    out.extend_from_slice(&u32::try_from(len).unwrap_or(0).to_le_bytes());
    out
}

/// Builds a dataspace message for `shape`.
fn dataspace_message(shape: &[usize]) -> Vec<u8> {
    let mut out = Vec::new();
    out.push(1);
    out.push(u8::try_from(shape.len()).unwrap_or(0));
    out.push(0);
    out.extend_from_slice(&[0; 5]);
    for &dim in shape {
        out.extend_from_slice(&u64::try_from(dim).unwrap_or(0).to_le_bytes());
    }
    out
}

/// Builds a version-3 contiguous data layout message.
fn layout_message(
    addr: u64,
    size: u64,
) -> Vec<u8> {
    let mut out = Vec::new();
    out.push(3);
    out.push(1);
    out.extend_from_slice(&addr.to_le_bytes());
    out.extend_from_slice(&size.to_le_bytes());
    out
}

/// Builds a version-2 fill value message with no stored fill value.
fn fill_value_message() -> Vec<u8> {
    vec![2, 2, 2, 0]
}

/// Builds an attribute message.
fn attribute_message(attribute: &Attribute) -> Result<Vec<u8>, String> {
    let (name_bytes, dtype, shape, body) = match &attribute.value {
        | AttrValue::Int(value) => {
            (
                attribute.name.as_bytes(),
                datatype_message(H5Type::I64)?,
                Vec::new(),
                value.to_le_bytes().to_vec(),
            )
        },
        | AttrValue::Float(value) => {
            (
                attribute.name.as_bytes(),
                datatype_message(H5Type::F64)?,
                Vec::new(),
                value.to_le_bytes().to_vec(),
            )
        },
        | AttrValue::Text(text) => {
            (
                attribute.name.as_bytes(),
                string_datatype_message(text.len()),
                Vec::new(),
                text.as_bytes().to_vec(),
            )
        },
        | AttrValue::Ints(values) => {
            (
                attribute.name.as_bytes(),
                datatype_message(H5Type::I64)?,
                vec![values.len()],
                values
                    .iter()
                    .flat_map(|value| value.to_le_bytes())
                    .collect(),
            )
        },
        | AttrValue::Floats(values) => {
            (
                attribute.name.as_bytes(),
                datatype_message(H5Type::F64)?,
                vec![values.len()],
                values
                    .iter()
                    .flat_map(|value| value.to_le_bytes())
                    .collect(),
            )
        },
    };
    let dataspace = dataspace_message(&shape);
    let name_len = name_bytes.len();
    let dtype_len = dtype.len();
    let dspace_len = dataspace.len();
    if name_len > u16::MAX as usize
        || dtype_len > u16::MAX as usize
        || dspace_len > u16::MAX as usize
    {
        return Err(format!("attribute {} is too large", attribute.name));
    }
    let mut out = Vec::new();
    out.push(3);
    out.push(0);
    out.extend_from_slice(&u16::try_from(name_len).unwrap_or(0).to_le_bytes());
    out.extend_from_slice(&u16::try_from(dtype_len).unwrap_or(0).to_le_bytes());
    out.extend_from_slice(&u16::try_from(dspace_len).unwrap_or(0).to_le_bytes());
    out.extend_from_slice(name_bytes);
    out.extend_from_slice(&dtype);
    out.extend_from_slice(&dataspace);
    out.extend_from_slice(&body);
    Ok(out)
}

/// Encodes `data` as little-endian bytes of `h5_type`.
fn encode_data(
    h5_type: H5Type,
    data: &[f64],
) -> Result<Vec<u8>, String> {
    let mut out = Vec::with_capacity(data.len().saturating_mul(type_size(h5_type)));
    match h5_type {
        | H5Type::Str => return Err("string datasets cannot be written".to_owned()),
        | H5Type::F32 => {
            for &value in data {
                #[allow(clippy::cast_possible_truncation)]
                out.extend_from_slice(&(value as f32).to_le_bytes());
            }
        },
        | H5Type::F64 => {
            for &value in data {
                out.extend_from_slice(&value.to_le_bytes());
            }
        },
        | other => {
            for &value in data {
                if !value.is_finite() {
                    return Err("integer datasets cannot store non-finite values".to_owned());
                }
                let truncated = value.trunc();
                match other {
                    | H5Type::I8 => {
                        if !(-128.0..=127.0).contains(&truncated) {
                            return Err(format!("value {value} does not fit in an i8"));
                        }
                        out.push(truncated as i8 as u8);
                    },
                    | H5Type::U8 => {
                        if !(0.0..=255.0).contains(&truncated) {
                            return Err(format!("value {value} does not fit in a u8"));
                        }
                        out.push(truncated as u8);
                    },
                    | H5Type::I16 => {
                        if !(-32_768.0..=32_767.0).contains(&truncated) {
                            return Err(format!("value {value} does not fit in an i16"));
                        }
                        out.extend_from_slice(&(truncated as i16).to_le_bytes());
                    },
                    | H5Type::U16 => {
                        if !(0.0..=65_535.0).contains(&truncated) {
                            return Err(format!("value {value} does not fit in a u16"));
                        }
                        out.extend_from_slice(&(truncated as u16).to_le_bytes());
                    },
                    | H5Type::I32 => {
                        if !(-2_147_483_648.0..=2_147_483_647.0).contains(&truncated) {
                            return Err(format!("value {value} does not fit in an i32"));
                        }
                        out.extend_from_slice(&(truncated as i32).to_le_bytes());
                    },
                    | H5Type::U32 => {
                        if !(0.0..=4_294_967_295.0).contains(&truncated) {
                            return Err(format!("value {value} does not fit in a u32"));
                        }
                        out.extend_from_slice(&(truncated as u32).to_le_bytes());
                    },
                    | H5Type::I64 => {
                        if !(-9_223_372_036_854_775_808.0..=9_223_372_036_854_775_807.0)
                            .contains(&truncated)
                        {
                            return Err(format!("value {value} does not fit in an i64"));
                        }
                        out.extend_from_slice(&(truncated as i64).to_le_bytes());
                    },
                    | _ => {
                        if !(0.0..=18_446_744_073_709_551_615.0).contains(&truncated) {
                            return Err(format!("value {value} does not fit in a u64"));
                        }
                        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
                        out.extend_from_slice(&(truncated as u64).to_le_bytes());
                    },
                }
            }
        },
    }
    Ok(out)
}

/// Serialises a version-1 object header holding `messages`.
fn object_header_v1(messages: &[(u16, Vec<u8>)]) -> Result<Vec<u8>, String> {
    let mut body = Vec::new();
    for (kind, data) in messages {
        let size = u16::try_from(data.len())
            .map_err(|_| "object header message is larger than 64 KiB".to_owned())?;
        body.extend_from_slice(&kind.to_le_bytes());
        body.extend_from_slice(&size.to_le_bytes());
        body.push(0);
        body.extend_from_slice(&[0; 3]);
        body.extend_from_slice(data);
        let padding = (8 - data.len() % 8) % 8;
        body.extend(std::iter::repeat_n(0_u8, padding));
    }
    let count =
        u16::try_from(messages.len()).map_err(|_| "too many object header messages".to_owned())?;
    let size =
        u32::try_from(body.len()).map_err(|_| "object header is larger than 4 GiB".to_owned())?;
    let mut out = Vec::with_capacity(16 + body.len());
    out.push(1);
    out.push(0);
    out.extend_from_slice(&count.to_le_bytes());
    out.extend_from_slice(&1_u32.to_le_bytes());
    out.extend_from_slice(&size.to_le_bytes());
    out.extend_from_slice(&[0; 4]);
    out.extend_from_slice(&body);
    Ok(out)
}

/// Byte length of the version-1 object header holding `messages`.
fn object_header_len(messages: &[(u16, Vec<u8>)]) -> usize {
    let body: usize = messages
        .iter()
        .map(|(_, data)| 8 + (data.len() + 7) / 8 * 8)
        .sum();
    16 + body
}

/// Normalises a dataset path to a leading slash.
fn normalise_path(path: &str) -> String {
    if path.is_empty() {
        "/".to_owned()
    } else if path.starts_with('/') {
        path.to_owned()
    } else {
        format!("/{path}")
    }
}

/// Returns the last path component, rejecting nested paths.
fn single_name(path: &str) -> Result<String, String> {
    let normalised = normalise_path(path);
    let trimmed = normalised.trim_start_matches('/');
    if trimmed.is_empty() {
        return Err("a dataset path cannot be the root group".to_owned());
    }
    if trimmed.contains('/') {
        return Err(format!(
            "nested dataset path {path} is not supported; the writer emits one root group"
        ));
    }
    Ok(trimmed.to_owned())
}

/// Collects the attributes declared for `path`.
fn attributes_for(
    attributes: &[(String, Vec<Attribute>)],
    path: &str,
) -> Vec<Attribute> {
    let wanted = normalise_path(path);
    let mut out = Vec::new();
    for (key, values) in attributes {
        if normalise_path(key) == wanted {
            out.extend(values.iter().cloned());
        }
    }
    out
}

/// Builds the object header messages for one dataset.
fn dataset_messages(
    dataset: &Dataset,
    attributes: &[Attribute],
    data_addr: u64,
    data_len: u64,
) -> Result<Vec<(u16, Vec<u8>)>, String> {
    let mut messages = vec![
        (0x01, dataspace_message(&dataset.shape)),
        (0x03, datatype_message(dataset.h5_type)?),
        (0x05, fill_value_message()),
        (0x08, layout_message(data_addr, data_len)),
    ];
    for attribute in attributes {
        messages.push((0x0C, attribute_message(attribute)?));
    }
    Ok(messages)
}

/// Serialises a local heap holding `names` and returns it with each offset.
fn local_heap(names: &[String]) -> (Vec<u8>, Vec<u64>) {
    let mut data = Vec::new();
    let mut offsets = Vec::with_capacity(names.len());
    for name in names {
        offsets.push(u64::try_from(data.len()).unwrap_or(0));
        data.extend_from_slice(name.as_bytes());
        data.push(0);
    }
    while data.len() % 8 != 0 {
        data.push(0);
    }
    let size = u64::try_from(data.len()).unwrap_or(0);
    let mut out = Vec::new();
    out.extend_from_slice(b"HEAP");
    out.push(0);
    out.extend_from_slice(&[0; 3]);
    out.extend_from_slice(&size.to_le_bytes());
    out.extend_from_slice(&1_u64.to_le_bytes());
    out.extend_from_slice(&0_u64.to_le_bytes());
    out.extend_from_slice(&data);
    (out, offsets)
}

/// Serialises a version-1 group B-tree with one leaf pointing at `snod_addr`.
fn group_btree(
    entries: &[(String, u64, u64)],
    snod_addr: u64,
) -> Vec<u8> {
    let mut out = Vec::new();
    out.extend_from_slice(b"TREE");
    out.push(0);
    out.push(0);
    let children = u16::from(!entries.is_empty());
    out.extend_from_slice(&children.to_le_bytes());
    out.extend_from_slice(&u64::MAX.to_le_bytes());
    out.extend_from_slice(&u64::MAX.to_le_bytes());
    if let Some((_, key, _)) = entries.first() {
        out.extend_from_slice(&key.to_le_bytes());
        out.extend_from_slice(&snod_addr.to_le_bytes());
    }
    out.extend_from_slice(&u64::MAX.to_le_bytes());
    out
}

/// Byte length of [`group_btree`] for `entries`.
fn group_btree_len(entries: &[(String, u64, u64)]) -> usize {
    if entries.is_empty() {
        24 + 8
    } else {
        24 + 16 + 8
    }
}

/// Serialises a symbol table node listing `entries`.
fn symbol_table_node(entries: &[(String, u64, u64)]) -> Vec<u8> {
    let mut out = Vec::new();
    out.extend_from_slice(b"SNOD");
    out.push(1);
    out.push(0);
    out.extend_from_slice(&u16::try_from(entries.len()).unwrap_or(0).to_le_bytes());
    for (_, name_offset, object) in entries {
        out.extend_from_slice(&name_offset.to_le_bytes());
        out.extend_from_slice(&object.to_le_bytes());
        out.extend_from_slice(&0_u32.to_le_bytes());
        out.extend_from_slice(&0_u32.to_le_bytes());
        out.extend_from_slice(&[0; 16]);
    }
    out
}

/// Serialises a version-0 superblock.
fn superblock_v0(
    root: u64,
    btree: u64,
    heap: u64,
    eof: u64,
) -> Vec<u8> {
    let mut out = Vec::new();
    out.extend_from_slice(&HDF5_SIGNATURE);
    out.push(0);
    out.push(0);
    out.push(0);
    out.push(0);
    out.push(0);
    out.push(8);
    out.push(8);
    out.push(0);
    out.extend_from_slice(&4_u16.to_le_bytes());
    out.extend_from_slice(&16_u16.to_le_bytes());
    out.extend_from_slice(&0_u32.to_le_bytes());
    out.extend_from_slice(&0_u64.to_le_bytes());
    out.extend_from_slice(&u64::MAX.to_le_bytes());
    out.extend_from_slice(&eof.to_le_bytes());
    out.extend_from_slice(&u64::MAX.to_le_bytes());
    out.extend_from_slice(&0_u64.to_le_bytes());
    out.extend_from_slice(&root.to_le_bytes());
    out.extend_from_slice(&1_u32.to_le_bytes());
    out.extend_from_slice(&0_u32.to_le_bytes());
    out.extend_from_slice(&btree.to_le_bytes());
    out.extend_from_slice(&heap.to_le_bytes());
    out
}

/// Writes datasets and attributes to an HDF5 file on disk.
///
/// # Errors
/// Fails when a dataset name is nested or duplicated, a string dataset is
/// requested, an integer value does not fit its declared type, or the file
/// cannot be written.
pub fn write_hdf5(
    path: &Path,
    datasets: &[Dataset],
    attributes: &[(String, Vec<Attribute>)],
) -> Result<(), String> {
    let bytes = write_hdf5_bytes(datasets, attributes)?;
    std::fs::write(path, bytes)
        .map_err(|error| format!("failed to write {}: {error}", path.display()))
}

/// Serialises datasets and attributes as an HDF5 file image.
///
/// # Errors
/// Fails when a dataset name is nested or duplicated, a string dataset is
/// requested, or an integer value does not fit its declared type.
pub fn write_hdf5_bytes(
    datasets: &[Dataset],
    attributes: &[(String, Vec<Attribute>)],
) -> Result<Vec<u8>, String> {
    let mut names = Vec::with_capacity(datasets.len());
    for dataset in datasets {
        let name = single_name(&dataset.path)?;
        if names.contains(&name) {
            return Err(format!("duplicate dataset name {name}"));
        }
        names.push(name);
    }
    let mut encoded = Vec::with_capacity(datasets.len());
    let mut lengths = Vec::with_capacity(datasets.len());
    for dataset in datasets {
        let count = element_count(&dataset.shape)?;
        let count = to_usize(count, "dataset element count")?;
        if count != dataset.data.len() {
            return Err(format!(
                "dataset {} has shape {:?} but {} values",
                dataset.path,
                dataset.shape,
                dataset.data.len()
            ));
        }
        let bytes = count.saturating_mul(type_size(dataset.h5_type));
        if bytes > MAX_ARRAY_BYTES {
            return Err(format!(
                "dataset {} exceeds the {MAX_ARRAY_BYTES}-byte limit",
                dataset.path
            ));
        }
        let body = encode_data(dataset.h5_type, &dataset.data)?;
        lengths.push(u64::try_from(body.len()).unwrap_or(0));
        encoded.push(body);
    }

    let root_attributes = attributes_for(attributes, "/");
    let mut root_messages = vec![(0x11_u16, vec![0_u8; 16])];
    for attribute in &root_attributes {
        root_messages.push((0x0C, attribute_message(attribute)?));
    }
    let root_len = object_header_len(&root_messages);

    let (mut heap, name_offsets) = local_heap(&names);
    let mut entries: Vec<(String, u64, u64)> = names
        .iter()
        .zip(name_offsets.iter())
        .map(|(name, offset)| (name.clone(), *offset, 0_u64))
        .collect();
    entries.sort_by(|left, right| left.0.cmp(&right.0));
    let btree_len = group_btree_len(&entries);
    let snod_len = 8 + entries.len() * 40;

    let mut header_lengths = Vec::with_capacity(datasets.len());
    for (index, dataset) in datasets.iter().enumerate() {
        let attrs = attributes_for(attributes, &dataset.path);
        let messages = dataset_messages(dataset, &attrs, 0, lengths[index])?;
        header_lengths.push(object_header_len(&messages));
    }

    let mut pos: u64 = 96;
    let root_addr = pos;
    pos += u64::try_from(root_len).unwrap_or(0);
    let heap_addr = pos;
    if let Some(slot) = heap.get_mut(24..32) {
        slot.copy_from_slice(&(heap_addr + 32).to_le_bytes());
    }
    pos += u64::try_from(heap.len()).unwrap_or(0);
    let btree_addr = pos;
    pos += u64::try_from(btree_len).unwrap_or(0);
    let snod_addr = pos;
    pos += u64::try_from(snod_len).unwrap_or(0);
    let mut dataset_addrs = Vec::with_capacity(datasets.len());
    for length in &lengths {
        pos = align8(pos);
        let header_addr = pos;
        let index = dataset_addrs.len();
        pos += u64::try_from(header_lengths[index]).unwrap_or(0);
        pos = align8(pos);
        let data_addr = pos;
        pos += *length;
        dataset_addrs.push((header_addr, data_addr));
    }
    let eof = align8(pos);

    for entry in &mut entries {
        if let Some(index) = names.iter().position(|name| *name == entry.0) {
            if let Some((header_addr, _)) = dataset_addrs.get(index) {
                entry.2 = *header_addr;
            }
        }
    }

    let mut out = Vec::new();
    out.extend_from_slice(&superblock_v0(root_addr, btree_addr, heap_addr, eof));
    let mut root_messages = vec![(
        0x11_u16,
        [btree_addr.to_le_bytes(), heap_addr.to_le_bytes()].concat(),
    )];
    for attribute in &root_attributes {
        root_messages.push((0x0C, attribute_message(attribute)?));
    }
    out.extend_from_slice(&object_header_v1(&root_messages)?);
    out.extend_from_slice(&heap);
    out.extend_from_slice(&group_btree(&entries, snod_addr));
    out.extend_from_slice(&symbol_table_node(&entries));
    for (index, dataset) in datasets.iter().enumerate() {
        let (header_addr, data_addr) = dataset_addrs[index];
        while u64::try_from(out.len()).unwrap_or(u64::MAX) < header_addr {
            out.push(0);
        }
        let attrs = attributes_for(attributes, &dataset.path);
        let messages = dataset_messages(dataset, &attrs, data_addr, lengths[index])?;
        out.extend_from_slice(&object_header_v1(&messages)?);
        while u64::try_from(out.len()).unwrap_or(u64::MAX) < data_addr {
            out.push(0);
        }
        out.extend_from_slice(&encoded[index]);
    }
    while u64::try_from(out.len()).unwrap_or(u64::MAX) < eof {
        out.push(0);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::AttrValue;
    use super::Attribute;
    use super::Dataset;
    use super::H5Type;
    use super::read_hdf5;
    use super::read_hdf5_bytes;
    use super::write_hdf5;
    use super::write_hdf5_bytes;

    /// A unique temporary path for one test.
    fn temp_path(tag: &str) -> std::path::PathBuf {
        let unique = format!(
            "rssn_h5_{tag}_{}_{}.h5",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .map_or(0, |d| d.as_nanos())
        );
        std::env::temp_dir().join(unique)
    }

    /// A one-dimensional `f64` dataset.
    fn line_dataset() -> Dataset {
        Dataset {
            path: "/line".to_owned(),
            shape: vec![6],
            h5_type: H5Type::F64,
            data: vec![0.0, 1.5, -2.5, 3.25, 4.0, -5.75],
        }
    }

    /// A two-dimensional `f64` dataset.
    fn grid_dataset() -> Dataset {
        Dataset {
            path: "/grid".to_owned(),
            shape: vec![2, 3],
            h5_type: H5Type::F64,
            data: vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
        }
    }

    /// An `i32` dataset.
    fn count_dataset() -> Dataset {
        Dataset {
            path: "/counts".to_owned(),
            shape: vec![4],
            h5_type: H5Type::I32,
            data: vec![-1.0, 0.0, 7.0, 1_000_000.0],
        }
    }

    /// The attributes used by the round-trip tests.
    fn sample_attributes() -> Vec<(String, Vec<Attribute>)> {
        vec![
            (
                "/".to_owned(),
                vec![
                    Attribute {
                        name: "title".to_owned(),
                        value: AttrValue::Text("round trip".to_owned()),
                    },
                    Attribute {
                        name: "version".to_owned(),
                        value: AttrValue::Int(3),
                    },
                ],
            ),
            (
                "/line".to_owned(),
                vec![
                    Attribute {
                        name: "scale".to_owned(),
                        value: AttrValue::Float(2.5),
                    },
                    Attribute {
                        name: "units".to_owned(),
                        value: AttrValue::Text("metres".to_owned()),
                    },
                ],
            ),
            (
                "/counts".to_owned(),
                vec![
                    Attribute {
                        name: "scale".to_owned(),
                        value: AttrValue::Float(2.5),
                    },
                    Attribute {
                        name: "units".to_owned(),
                        value: AttrValue::Text("metres".to_owned()),
                    },
                ],
            ),
        ]
    }

    #[test]
    fn round_trip_one_and_two_dimensional_float_datasets() {
        let datasets = vec![line_dataset(), grid_dataset()];
        let bytes = write_hdf5_bytes(&datasets, &[]).expect("serialise");
        let file = read_hdf5_bytes(&bytes).expect("parse");
        assert_eq!(file.groups, vec!["/".to_owned()]);
        assert_eq!(file.datasets.len(), 2);
        let line = file
            .datasets
            .iter()
            .find(|d| d.path == "/line")
            .expect("line dataset");
        assert_eq!(line.shape, vec![6]);
        assert_eq!(line.h5_type, H5Type::F64);
        assert_eq!(line.data, line_dataset().data);
        let grid = file
            .datasets
            .iter()
            .find(|d| d.path == "/grid")
            .expect("grid dataset");
        assert_eq!(grid.shape, vec![2, 3]);
        assert_eq!(grid.data, grid_dataset().data);
    }

    #[test]
    fn round_trip_integer_dataset_and_attributes() {
        let datasets = vec![count_dataset()];
        let attributes = sample_attributes();
        let bytes = write_hdf5_bytes(&datasets, &attributes).expect("serialise");
        let file = read_hdf5_bytes(&bytes).expect("parse");
        let counts = file
            .datasets
            .iter()
            .find(|d| d.path == "/counts")
            .expect("counts dataset");
        assert_eq!(counts.h5_type, H5Type::I32);
        assert_eq!(counts.data, vec![-1.0, 0.0, 7.0, 1_000_000.0]);
        let root_attrs = file
            .attributes
            .iter()
            .find(|(path, _)| path == "/")
            .expect("root attributes");
        assert_eq!(root_attrs.1.len(), 2);
        assert_eq!(root_attrs.1[0].name, "title");
        assert_eq!(
            root_attrs.1[0].value,
            AttrValue::Text("round trip".to_owned())
        );
        assert_eq!(root_attrs.1[1].value, AttrValue::Int(3));
        let line_attrs = file
            .attributes
            .iter()
            .find(|(path, _)| path == "/counts")
            .expect("counts attributes");
        assert_eq!(line_attrs.1[0].value, AttrValue::Float(2.5));
    }

    #[test]
    fn round_trip_through_a_file() {
        let path = temp_path("roundtrip");
        let datasets = vec![line_dataset(), count_dataset()];
        let attributes = sample_attributes();
        write_hdf5(&path, &datasets, &attributes).expect("write");
        let file = read_hdf5(&path).expect("read");
        assert_eq!(file.datasets.len(), 2);
        let _ = std::fs::remove_file(&path);
        assert!(!path.exists());
    }

    #[test]
    fn truncated_buffers_are_rejected() {
        let datasets = vec![line_dataset(), grid_dataset()];
        let attributes = sample_attributes();
        let bytes = write_hdf5_bytes(&datasets, &attributes).expect("serialise");
        let len = bytes.len();
        let cuts = [4_usize, 20, 60, 96, len / 2];
        for cut in cuts {
            let truncated = bytes.get(..cut).expect("slice");
            assert!(
                read_hdf5_bytes(truncated).is_err(),
                "a buffer truncated at {cut} of {len} must be rejected"
            );
        }
    }

    #[test]
    fn nested_and_duplicate_names_are_rejected() {
        let mut nested = line_dataset();
        nested.path = "/a/b".to_owned();
        assert!(write_hdf5_bytes(&[nested], &[]).is_err());
        let duplicate = vec![line_dataset(), line_dataset()];
        assert!(write_hdf5_bytes(&duplicate, &[]).is_err());
    }

    #[test]
    fn string_datasets_are_rejected() {
        let dataset = Dataset {
            path: "/text".to_owned(),
            shape: vec![2],
            h5_type: H5Type::Str,
            data: vec![0.0, 0.0],
        };
        assert!(write_hdf5_bytes(&[dataset], &[]).is_err());
    }

    #[test]
    fn shape_mismatch_is_rejected() {
        let dataset = Dataset {
            path: "/bad".to_owned(),
            shape: vec![3],
            h5_type: H5Type::F64,
            data: vec![1.0, 2.0],
        };
        assert!(write_hdf5_bytes(&[dataset], &[]).is_err());
    }

    #[test]
    fn reads_version_three_superblock_and_version_two_headers() {
        let mut out: Vec<u8> = Vec::new();
        out.extend_from_slice(&[0x89, b'H', b'D', b'F', 0x0d, 0x0a, 0x1a, 0x0a]);
        out.push(3);
        out.push(8);
        out.push(8);
        out.push(0);
        out.extend_from_slice(&0_u64.to_le_bytes());
        out.extend_from_slice(&u64::MAX.to_le_bytes());
        let eof_patch = out.len();
        out.extend_from_slice(&0_u64.to_le_bytes());
        let root_patch = out.len();
        out.extend_from_slice(&0_u64.to_le_bytes());
        out.extend_from_slice(&0_u32.to_le_bytes());

        let root_addr = out.len() as u64;
        patch_u64(&mut out, root_patch, root_addr);
        let mut link = vec![1_u8, 0, 1];
        link.push(b'd');
        link.extend_from_slice(&0_u64.to_le_bytes());
        let mut body = Vec::new();
        body.push(0x06_u8);
        body.extend_from_slice(&(link.len() as u16).to_le_bytes());
        body.push(0);
        body.extend_from_slice(&link);
        out.extend_from_slice(b"OHDR");
        out.push(2);
        out.push(0);
        out.push(body.len() as u8);
        let body_start = out.len();
        out.extend_from_slice(&body);
        out.extend_from_slice(&0_u32.to_le_bytes());
        let link_addr_patch = body_start + 8;

        let dataset_addr = out.len() as u64;
        patch_u64(&mut out, link_addr_patch, dataset_addr);
        let mut dataspace = vec![2_u8, 1, 0, 1];
        dataspace.extend_from_slice(&2_u64.to_le_bytes());
        let mut datatype = vec![0x11_u8, 0x20, 0, 0];
        datatype.extend_from_slice(&8_u32.to_le_bytes());
        datatype.extend_from_slice(&0_u16.to_le_bytes());
        datatype.extend_from_slice(&64_u16.to_le_bytes());
        datatype.extend_from_slice(&[52, 11, 0, 52]);
        datatype.extend_from_slice(&1023_u32.to_le_bytes());
        let mut layout = vec![3_u8, 1];
        layout.extend_from_slice(&0_u64.to_le_bytes());
        layout.extend_from_slice(&16_u64.to_le_bytes());
        let mut messages = Vec::new();
        for (kind, data) in [(0x01_u8, &dataspace), (0x03, &datatype), (0x08, &layout)] {
            messages.push(kind);
            messages.extend_from_slice(&(data.len() as u16).to_le_bytes());
            messages.push(0);
            messages.extend_from_slice(data);
        }
        out.extend_from_slice(b"OHDR");
        out.push(2);
        out.push(0);
        out.push(messages.len() as u8);
        let data_patch = out.len() + messages.len() - 18 + 2;
        out.extend_from_slice(&messages);
        out.extend_from_slice(&0_u32.to_le_bytes());
        let data_addr = out.len() as u64;
        patch_u64(&mut out, data_patch, data_addr);
        out.extend_from_slice(&1.5_f64.to_le_bytes());
        out.extend_from_slice(&(-2.25_f64).to_le_bytes());
        let eof = out.len() as u64;
        patch_u64(&mut out, eof_patch, eof);

        let file = read_hdf5_bytes(&out).expect("parse version 3 file");
        assert_eq!(file.groups, vec!["/".to_owned()]);
        assert_eq!(file.datasets.len(), 1);
        let dataset = &file.datasets[0];
        assert_eq!(dataset.path, "/d");
        assert_eq!(dataset.shape, vec![2]);
        assert_eq!(dataset.h5_type, H5Type::F64);
        assert_eq!(dataset.data, vec![1.5, -2.25]);
    }

    #[test]
    fn shuffle_filter_is_reversed() {
        let shuffled = vec![1_u8, 4, 2, 5, 3, 6];
        let restored = super::unshuffle(&shuffled, 3).expect("unshuffle");
        assert_eq!(restored, vec![1_u8, 2, 3, 4, 5, 6]);
        assert!(super::unshuffle(&shuffled, 0).is_err());
    }

    #[test]
    fn chunked_deflate_dataset_decodes() {
        let bytes = build_chunked_fixture();
        let file = read_hdf5_bytes(&bytes).expect("parse chunked fixture");
        let dataset = file
            .datasets
            .iter()
            .find(|d| d.path == "/data")
            .expect("chunked dataset");
        assert_eq!(dataset.shape, vec![4, 4]);
        assert_eq!(dataset.h5_type, H5Type::F32);
        let expected: Vec<f64> = (0..16).map(f64::from).collect();
        assert_eq!(dataset.data, expected);
    }

    // ---- synthetic chunked fixture -------------------------------------

    /// Patches a little-endian `u64` at `at`.
    fn patch_u64(
        out: &mut [u8],
        at: usize,
        value: u64,
    ) {
        if let Some(slot) = out.get_mut(at..at + 8) {
            slot.copy_from_slice(&value.to_le_bytes());
        }
    }

    /// Appends one version-1 object header message to `out`.
    fn push_message(
        out: &mut Vec<u8>,
        kind: u16,
        data: &[u8],
    ) {
        out.extend_from_slice(&kind.to_le_bytes());
        out.extend_from_slice(&(data.len() as u16).to_le_bytes());
        out.push(0);
        out.extend_from_slice(&[0; 3]);
        out.extend_from_slice(data);
        let padding = (8 - data.len() % 8) % 8;
        out.extend(std::iter::repeat_n(0_u8, padding));
    }

    /// Wraps pre-serialised `messages` in a version-1 object header.
    fn object_header(messages: &[u8]) -> Vec<u8> {
        let mut count = 0_u16;
        let mut pos = 0_usize;
        while pos + 8 <= messages.len() {
            let size = u16::from_le_bytes([messages[pos + 2], messages[pos + 3]]);
            count += 1;
            pos += 8 + (usize::from(size) + 7) / 8 * 8;
        }
        let mut out = Vec::new();
        out.push(1);
        out.push(0);
        out.extend_from_slice(&count.to_le_bytes());
        out.extend_from_slice(&1_u32.to_le_bytes());
        out.extend_from_slice(&(messages.len() as u32).to_le_bytes());
        out.extend_from_slice(&[0; 4]);
        out.extend_from_slice(messages);
        out
    }

    /// The Adler-32 checksum used by zlib.
    fn adler32(data: &[u8]) -> u32 {
        let mut s1: u32 = 1;
        let mut s2: u32 = 0;
        for &byte in data {
            s1 = (s1 + u32::from(byte)) % 65_521;
            s2 = (s2 + s1) % 65_521;
        }
        (s2 << 16) | s1
    }

    /// Wraps `raw` in a zlib stream made of one stored DEFLATE block.
    fn zlib_stored(raw: &[u8]) -> Vec<u8> {
        let mut out = vec![0x78, 0x01, 0x01];
        let len = raw.len() as u16;
        out.extend_from_slice(&len.to_le_bytes());
        out.extend_from_slice(&(!len).to_le_bytes());
        out.extend_from_slice(raw);
        out.extend_from_slice(&adler32(raw).to_be_bytes());
        out
    }

    /// Builds a synthetic 4x4 `f32` dataset stored in 2x2 chunks, each
    /// compressed with a zlib stream of one stored DEFLATE block.
    fn build_chunked_fixture() -> Vec<u8> {
        let mut out: Vec<u8> = Vec::new();
        out.extend_from_slice(&[0x89, b'H', b'D', b'F', 0x0d, 0x0a, 0x1a, 0x0a]);
        out.push(0);
        out.push(0);
        out.push(0);
        out.push(0);
        out.push(0);
        out.push(8);
        out.push(8);
        out.push(0);
        out.extend_from_slice(&4_u16.to_le_bytes());
        out.extend_from_slice(&16_u16.to_le_bytes());
        out.extend_from_slice(&0_u32.to_le_bytes());
        out.extend_from_slice(&0_u64.to_le_bytes());
        out.extend_from_slice(&u64::MAX.to_le_bytes());
        out.extend_from_slice(&0_u64.to_le_bytes());
        out.extend_from_slice(&u64::MAX.to_le_bytes());
        out.extend_from_slice(&0_u64.to_le_bytes());
        let root_patch = out.len();
        out.extend_from_slice(&0_u64.to_le_bytes());
        out.extend_from_slice(&0_u32.to_le_bytes());
        out.extend_from_slice(&0_u32.to_le_bytes());
        out.extend_from_slice(&[0; 16]);
        let root_addr = out.len() as u64;
        patch_u64(&mut out, root_patch, root_addr);

        // Root group: one compact hard link named "data".
        let mut link = Vec::new();
        link.push(1);
        link.push(0x00);
        link.push(4);
        link.extend_from_slice(b"data");
        link.extend_from_slice(&0_u64.to_le_bytes());
        let mut root_messages = Vec::new();
        push_message(&mut root_messages, 0x06, &link);
        let dataset_addr_patch = out.len() + 16 + 8 + 2 + 1 + 4;
        out.extend_from_slice(&object_header(&root_messages));
        let dataset_addr = out.len() as u64;
        patch_u64(&mut out, dataset_addr_patch, dataset_addr);

        // Dataset object header: 4x4 f32, chunked 2x2, deflate filter.
        let mut dataspace = vec![1, 2, 0, 0, 0, 0, 0, 0];
        dataspace.extend_from_slice(&4_u64.to_le_bytes());
        dataspace.extend_from_slice(&4_u64.to_le_bytes());
        let mut datatype = vec![0x11, 0x20, 0, 0];
        datatype.extend_from_slice(&4_u32.to_le_bytes());
        datatype.extend_from_slice(&0_u16.to_le_bytes());
        datatype.extend_from_slice(&32_u16.to_le_bytes());
        datatype.extend_from_slice(&[23, 8, 0, 23]);
        datatype.extend_from_slice(&127_u32.to_le_bytes());
        let mut layout = vec![3, 2, 3];
        let btree_patch = out.len() + 16 + (8 + 24) + (8 + 24) + 8 + 3;
        layout.extend_from_slice(&0_u64.to_le_bytes());
        layout.extend_from_slice(&2_u32.to_le_bytes());
        layout.extend_from_slice(&2_u32.to_le_bytes());
        layout.extend_from_slice(&4_u32.to_le_bytes());
        let filter = vec![2_u8, 1, 1, 0, 0, 0, 0, 0];
        let mut messages = Vec::new();
        push_message(&mut messages, 0x01, &dataspace);
        push_message(&mut messages, 0x03, &datatype);
        push_message(&mut messages, 0x08, &layout);
        push_message(&mut messages, 0x0B, &filter);
        out.extend_from_slice(&object_header(&messages));

        let corners = [(0_u32, 0_u32), (0, 2), (2, 0), (2, 2)];
        let mut chunk_addrs = Vec::new();
        for (row, col) in corners {
            let mut raw = Vec::new();
            for dr in 0..2_u32 {
                for dc in 0..2_u32 {
                    let value = ((row + dr) * 4 + (col + dc)) as f32;
                    raw.extend_from_slice(&value.to_le_bytes());
                }
            }
            chunk_addrs.push(out.len() as u64);
            out.extend_from_slice(&zlib_stored(&raw));
        }

        let btree_addr = out.len() as u64;
        patch_u64(&mut out, btree_patch, btree_addr);
        out.extend_from_slice(b"TREE");
        out.push(1);
        out.push(0);
        out.extend_from_slice(&4_u16.to_le_bytes());
        out.extend_from_slice(&u64::MAX.to_le_bytes());
        out.extend_from_slice(&u64::MAX.to_le_bytes());
        let chunk_len = zlib_stored(&[0_u8; 16]).len() as u32;
        for (index, (row, col)) in corners.iter().enumerate() {
            out.extend_from_slice(&chunk_len.to_le_bytes());
            out.extend_from_slice(&0_u32.to_le_bytes());
            out.extend_from_slice(&row.to_le_bytes());
            out.extend_from_slice(&col.to_le_bytes());
            out.extend_from_slice(&0_u32.to_le_bytes());
            out.extend_from_slice(&chunk_addrs[index].to_le_bytes());
        }
        out.extend_from_slice(&0_u32.to_le_bytes());
        out.extend_from_slice(&0_u32.to_le_bytes());
        out.extend_from_slice(&0_u32.to_le_bytes());
        out.extend_from_slice(&0_u32.to_le_bytes());
        out.extend_from_slice(&0_u32.to_le_bytes());
        out.extend_from_slice(&u64::MAX.to_le_bytes());
        out
    }
}
