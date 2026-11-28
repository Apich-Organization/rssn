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

pub(crate) mod checksum;
pub(crate) mod filters;
pub(crate) mod inflate;
pub mod writer;

pub use writer::DatasetWriteOptions;
pub use writer::H5Filter;
pub use writer::H5WriteOptions;
pub use writer::SoftLink;
pub use writer::write_hdf5_bytes_with_options;
pub use writer::write_hdf5_with_options;

use filters::Filter;
use filters::FilterDecoder;
use filters::apply_filters;

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
    /// Fixed-length or variable-length string.
    Str,
    /// Variable-length sequence of numeric values (stored in [`Dataset::vlen_data`]).
    VlenNumeric,
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
    /// Row-major values, converted to `f64`. Empty for string or vlen datasets.
    pub data: Vec<f64>,
    /// Row-major strings. Empty unless [`H5Type::Str`] is the element type.
    pub strings: Vec<String>,
    /// Variable-length numeric sequences. Empty unless [`H5Type::VlenNumeric`] is the element type.
    pub vlen_data: Vec<Vec<f64>>,
    /// Base numeric element type when [`H5Type::VlenNumeric`] is used.
    pub vlen_base_type: Option<H5Type>,
}

impl Dataset {
    /// A numeric dataset.
    #[must_use]
    pub fn new_numeric(
        path: impl Into<String>,
        shape: Vec<usize>,
        h5_type: H5Type,
        data: Vec<f64>,
    ) -> Self {
        Self {
            path: path.into(),
            shape,
            h5_type,
            data,
            strings: Vec::new(),
            vlen_data: Vec::new(),
            vlen_base_type: None,
        }
    }

    /// A fixed-width string dataset.
    #[must_use]
    pub fn new_strings(
        path: impl Into<String>,
        shape: Vec<usize>,
        strings: Vec<String>,
    ) -> Self {
        Self {
            path: path.into(),
            shape,
            h5_type: H5Type::Str,
            data: Vec::new(),
            strings,
            vlen_data: Vec::new(),
            vlen_base_type: None,
        }
    }

    /// A variable-length numeric dataset.
    #[must_use]
    pub fn new_vlen_numeric(
        path: impl Into<String>,
        shape: Vec<usize>,
        base_type: H5Type,
        vlen_data: Vec<Vec<f64>>,
    ) -> Self {
        Self {
            path: path.into(),
            shape,
            h5_type: H5Type::VlenNumeric,
            data: Vec::new(),
            strings: Vec::new(),
            vlen_data,
            vlen_base_type: Some(base_type),
        }
    }
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
    /// Symbolic (soft) links as `(link_path, resolved_target_path)`.
    pub soft_links: Vec<(String, String)>,
}

impl H5File {
    /// Looks up a soft link's resolved target path by its link path.
    #[must_use]
    pub fn resolve_soft_link(
        &self,
        path: &str,
    ) -> Option<&str> {
        let norm = normalise_path(path);
        self.soft_links
            .iter()
            .find(|(link, _)| *link == norm)
            .map(|(_, target)| target.as_str())
    }
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
            // The header block is padded to eight bytes, so a shorter tail is
            // slack rather than a broken message.
            break;
        }
        let kind = reader.u16_at(pos)?;
        let msize = u64::from(reader.u16_at(pos + 2)?);
        if kind == 0 {
            break;
        }
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
    // When the object tracks attribute creation order, every message header
    // carries a two-byte creation order after its four fixed bytes.
    let order_width = if flags & 0x04 != 0 {
        2_u64
    } else {
        0
    };
    let size_width = 1_usize << (flags & 0x03);
    let chunk_size = reader.uint_at(pos, size_width)?;
    pos += size_width as u64;
    let end = pos
        .checked_add(chunk_size)
        .ok_or_else(|| "object header size overflow".to_owned())?;
    while pos < end {
        if pos + 4 > end {
            // Trailing bytes below one message header are alignment padding.
            break;
        }
        let kind = u16::from(reader.u8_at(pos)?);
        let msize = u64::from(reader.u16_at(pos + 1)?);
        let mflags = reader.u8_at(pos + 3)?;
        pos += 4 + order_width;
        if kind == 0 {
            // A nil message is the filler that occupies the rest of the
            // chunk; its size field is not meaningful.
            break;
        }
        let data_pos = pos;
        let data = reader.bytes(data_pos, msize)?;
        pos += msize;
        if pos > end {
            return Err(format!(
                "object header message at {addr:#x} overruns its block"
            ));
        }
        if kind == 0x10 {
            if msize >= (reader.off_size + reader.len_size) as u64 {
                let next_addr = reader.off_at(data_pos)?;
                let next_len = reader.len_at(data_pos + reader.off_size as u64)?;
                if next_addr != u64::MAX {
                    read_continuation_v2(reader, next_addr, next_len, out, depth + 1)?;
                }
            }
        } else if mflags & 0x02 == 0 {
            // A shared message stores a reference to another object header,
            // not the message itself; skipping it keeps the rest readable.
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
    /// Enumeration: the stored value is the underlying integer, which is what
    /// this API exposes.
    Enum {
        /// Width of the underlying integer in bytes.
        size: usize,
        /// Whether the underlying integer is signed.
        signed: bool,
        /// Whether the bytes are big-endian.
        be: bool,
    },
    /// Variable-length string: each element is a 16-byte descriptor holding a
    /// length, a global-heap address and an object index.
    VlenStr,
    /// Variable-length numeric sequence: 16-byte descriptor pointing to global heap.
    VlenNumeric {
        /// Base element width in bytes.
        size: usize,
        /// Whether the base type is signed (integers).
        signed: bool,
        /// Whether the base type is big-endian.
        be: bool,
        /// Whether the base type is IEEE float.
        is_float: bool,
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
            | Self::Int { size, .. }
            | Self::Float { size, .. }
            | Self::Str { size }
            | Self::Enum { size, .. } => size,
            | Self::VlenStr | Self::VlenNumeric { .. } => 16,
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
        | 8 => {
            // Enumeration: the base integer datatype starts right after the
            // eight-byte class header; the member names are not needed.
            match parse_datatype(data.get(8..).ok_or("truncated enum datatype")?)? {
                | DType::Int { size, signed, be } => Ok(DType::Enum { size, signed, be }),
                | _ => Err("enumeration is not based on an integer type".to_owned()),
            }
        },
        | 9 => {
            // Variable-length sequence; type 1 is a string, type 0 is sequence.
            let kind = bits.first().copied().unwrap_or(0) & 0x0f;
            if kind == 1 {
                Ok(DType::VlenStr)
            } else if kind == 0 {
                match parse_datatype(data.get(8..).ok_or("truncated vlen datatype")?)? {
                    | DType::Int { size, signed, be } => {
                        Ok(DType::VlenNumeric {
                            size,
                            signed,
                            be,
                            is_float: false,
                        })
                    },
                    | DType::Float { size, be } => {
                        Ok(DType::VlenNumeric {
                            size,
                            signed: true,
                            be,
                            is_float: true,
                        })
                    },
                    | _ => Ok(DType::Other { class: 9 }),
                }
            } else {
                Ok(DType::Other { class: 9 })
            }
        },
        | other => Ok(DType::Other { class: other }),
    }
}

/// Chunk indexing method for version-4 chunked layouts.
#[derive(Debug, Clone)]
enum V4ChunkIndex {
    Single {
        filtered_size: Option<u64>,
        filter_mask: u32,
    },
    Implicit,
    FixedArray {
        page_bits: u8,
    },
    ExtensibleArray {},
    BTreeV2 {},
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
    /// A version-4 chunked dataset.
    ChunkedV4 {
        /// Address of the chunk or chunk index.
        addr: u64,
        /// Chunk dimensions (length == dataspace rank).
        dims: Vec<usize>,
        /// Indexing method.
        index: V4ChunkIndex,
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
        | 4 | 5 => {
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
                    let flags = *data.get(2).ok_or("truncated v4 chunked layout")?;
                    let ndims = usize::from(*data.get(3).ok_or("truncated v4 chunked layout")?);
                    let dim_len = usize::from(*data.get(4).ok_or("truncated v4 chunked layout")?);
                    if dim_len == 0 || dim_len > 8 {
                        return Err(format!("invalid v4 chunk dim width {dim_len}"));
                    }
                    let mut pos = 5_usize;
                    let mut dims = Vec::with_capacity(ndims);
                    for _ in 0..ndims {
                        dims.push(to_usize(read_uint(data, pos, dim_len)?, "chunk dimension")?);
                        pos += dim_len;
                    }
                    let idx_type = *data.get(pos).ok_or("truncated v4 chunk index type")?;
                    pos += 1;
                    let (index, addr) = match idx_type {
                        | 1 => {
                            // Single chunk
                            let (filtered_size, filter_mask) = if flags & 0x02 != 0 {
                                let fs = read_uint(data, pos, reader.len_size)?;
                                pos += reader.len_size;
                                let fm = u32::try_from(read_uint(data, pos, 4)?).unwrap_or(0);
                                pos += 4;
                                (Some(fs), fm)
                            } else {
                                (None, 0)
                            };
                            let addr = read_uint(data, pos, reader.off_size)?;
                            (
                                V4ChunkIndex::Single {
                                    filtered_size,
                                    filter_mask,
                                },
                                addr,
                            )
                        },
                        | 2 => {
                            // Implicit
                            let addr = read_uint(data, pos, reader.off_size)?;
                            (V4ChunkIndex::Implicit, addr)
                        },
                        | 3 => {
                            // Fixed Array
                            let page_bits = *data.get(pos).ok_or("truncated fixed array")?;
                            pos += 1;
                            let addr = read_uint(data, pos, reader.off_size)?;
                            (V4ChunkIndex::FixedArray { page_bits }, addr)
                        },
                        | 4 => {
                            // Extensible Array: 5 creation parameter bytes before addr
                            pos += 5;
                            let addr = read_uint(data, pos, reader.off_size)?;
                            (V4ChunkIndex::ExtensibleArray {}, addr)
                        },
                        | 5 => {
                            // Version-2 B-tree
                            pos += 6; // node_size(4) + split_percent(1) + merge_percent(1)
                            let addr = read_uint(data, pos, reader.off_size)?;
                            (V4ChunkIndex::BTreeV2 {}, addr)
                        },
                        | other => {
                            return Err(format!("unsupported v4 chunk index type {other}"));
                        },
                    };
                    Ok(Layout::ChunkedV4 { addr, dims, index })
                },
                | other => Err(format!("unsupported v4 data layout class {other}")),
            }
        },
        | other => Err(format!("unsupported data layout version {other}")),
    }
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
                let name_len = if id >= 256 {
                    let nl = to_usize(read_uint(data, pos, 2)?, "filter name length")?;
                    pos += 2;
                    nl
                } else {
                    0
                };
                let nclient = to_usize(read_uint(data, pos + 2, 2)?, "filter client count")?;
                pos += 4;
                if id >= 256 && name_len > 0 {
                    pos = pos
                        .checked_add(name_len)
                        .ok_or_else(|| "filter name length overflow".to_owned())?;
                }
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
/// Reads one object out of a global heap collection (`GCOL`).
fn global_heap_object(
    reader: &Reader<'_>,
    addr: u64,
    index: u32,
) -> Result<Vec<u8>, String> {
    let header = reader.bytes(addr, 16)?;
    if header.get(..4) != Some(b"GCOL") {
        return Err("vlen string does not point at a global heap collection".to_owned());
    }
    let size = reader.len_at(addr + 8)?;
    let body = reader.bytes(addr, size)?;
    let mut pos = 16usize;
    while pos + 16 <= body.len() {
        let object_index = u16::from_le_bytes([body[pos], body[pos + 1]]);
        let object_size = usize::try_from(u64::from_le_bytes(
            body.get(pos + 8..pos + 16)
                .ok_or("truncated global heap object header")?
                .try_into()
                .map_err(|_| "truncated global heap object header".to_owned())?,
        ))
        .map_err(|_| "global heap object is too large".to_owned())?;
        let start = pos + 16;
        if object_index == 0 {
            break;
        }
        if u32::from(object_index) == index {
            let end = start
                .checked_add(object_size)
                .ok_or_else(|| "global heap object size overflows".to_owned())?;
            return body
                .get(start..end)
                .map(<[u8]>::to_vec)
                .ok_or_else(|| "global heap object is truncated".to_owned());
        }
        let advance = to_usize(
            align8(u64::try_from(object_size).unwrap_or(0)),
            "heap object",
        )?;
        pos = start
            .checked_add(advance)
            .ok_or_else(|| "global heap offset overflows".to_owned())?;
    }
    Err(format!("global heap object {index} is missing"))
}

/// Decodes `count` variable-length string descriptors.
fn decode_vlen_strings(
    reader: &Reader<'_>,
    raw: &[u8],
    count: usize,
) -> Result<Vec<String>, String> {
    let mut out = Vec::with_capacity(count);
    for index in 0..count {
        let base = index.saturating_mul(16);
        let descriptor = raw
            .get(base..base + 16)
            .ok_or_else(|| "vlen string descriptor is truncated".to_owned())?;
        let length =
            u32::from_le_bytes([descriptor[0], descriptor[1], descriptor[2], descriptor[3]]);
        let mut addr_bytes = [0u8; 8];
        addr_bytes.copy_from_slice(&descriptor[4..12]);
        let addr = u64::from_le_bytes(addr_bytes);
        let mut index_bytes = [0u8; 4];
        index_bytes.copy_from_slice(&descriptor[12..16]);
        let object = u32::from_le_bytes(index_bytes);
        if length == 0 || addr == u64::MAX {
            out.push(String::new());
            continue;
        }
        let bytes = global_heap_object(reader, addr, object)?;
        let take = bytes
            .get(..usize::try_from(length).unwrap_or(usize::MAX))
            .unwrap_or(&bytes);
        out.push(String::from_utf8_lossy(take).into_owned());
    }
    Ok(out)
}

/// One member of a compound datatype.
struct CompoundMember {
    /// Member name.
    name: String,
    /// Byte offset of the member inside one compound element.
    offset: usize,
    /// Member datatype.
    dtype: DType,
}

/// Byte length of one datatype message, so a compound can walk its members.
///
/// Only the classes whose length is fixed by the format are accepted; a
/// compound member of any other class makes the whole dataset unreadable.
fn datatype_length(data: &[u8]) -> Result<usize, String> {
    let class = data.first().ok_or("empty datatype message")? & 0x0f;
    match class {
        | 0 | 1 => Ok(12),
        | 2 | 3 | 7 => Ok(8),
        | 9 => {
            let base = datatype_length(data.get(8..).ok_or("truncated vlen datatype")?)?;
            Ok(8 + base)
        },
        | other => {
            Err(format!(
                "compound member of datatype class {other} is not supported"
            ))
        },
    }
}

/// Parses the members of a compound datatype message.
///
/// The layout follows the format: a member name (NUL-terminated and, before
/// version 3, padded to eight bytes), the member byte offset, then the member
/// datatype message.
fn parse_compound_members(
    data: &[u8],
    total_size: usize,
) -> Result<Vec<CompoundMember>, String> {
    let version = *data.first().ok_or("empty compound datatype")? >> 4;
    let count = to_usize(read_uint(data, 1, 2)?, "compound member count")?;
    let mut pos = 8usize;
    let mut members = Vec::with_capacity(count);
    for _ in 0..count {
        let terminator = data
            .get(pos..)
            .and_then(|rest| rest.iter().position(|&byte| byte == 0))
            .ok_or("compound member name is not terminated")?;
        let name = std::str::from_utf8(data.get(pos..pos + terminator).unwrap_or(&[]))
            .map_err(|_| "compound member name is not valid UTF-8".to_owned())?
            .to_owned();
        let stored = terminator.saturating_add(1);
        pos += if version >= 3 {
            stored
        } else {
            usize::try_from(align8(u64::try_from(stored).unwrap_or(0))).unwrap_or(stored)
        };
        let offset = match version {
            | 1 => {
                let offset = to_usize(read_uint(data, pos, 4)?, "compound member offset")?;
                // Version 1 stores a 32-byte property descriptor.
                pos += 32;
                offset
            },
            | 2 => {
                let offset = to_usize(read_uint(data, pos, 4)?, "compound member offset")?;
                pos += 4;
                offset
            },
            | _ => {
                // Version 3 sizes the offset from the element size.
                let width = total_size.saturating_sub(1).bit_width() as usize;
                let width = width.div_ceil(8).clamp(1, 8);
                let offset = to_usize(read_uint(data, pos, width)?, "compound member offset")?;
                pos += width;
                offset
            },
        };
        let rest = data.get(pos..).ok_or("truncated compound member")?;
        let dtype = parse_datatype(rest)?;
        pos += datatype_length(rest)?;
        members.push(CompoundMember { name, offset, dtype });
    }
    Ok(members)
}

/// Decodes one compound dataset into one dataset per member.
///
/// A member is exposed when it can be represented: integers and floats become
/// values, fixed-width strings become text. Anything else is skipped, which is
/// why the returned vector can be shorter than the member list.
fn read_compound_dataset(
    dataset: &Dataset,
    dtype_bytes: &[u8],
    body: &[u8],
    count: usize,
    total_size: usize,
) -> Result<Vec<Dataset>, String> {
    let members = parse_compound_members(dtype_bytes, total_size)?;
    let mut out = Vec::with_capacity(members.len());
    for member in members {
        let width = member.dtype.size();
        if width == 0 {
            continue;
        }
        let numeric = matches!(
            member.dtype,
            DType::Int { .. } | DType::Float { .. } | DType::Enum { .. }
        );
        let string = matches!(member.dtype, DType::Str { .. });
        if !numeric && !string {
            continue;
        }
        let mut values = Vec::with_capacity(count);
        let mut texts: Vec<String> = Vec::new();
        for index in 0..count {
            let start = member
                .offset
                .saturating_add(index.saturating_mul(total_size));
            let end = start.saturating_add(width);
            let slice = body
                .get(start..end)
                .ok_or_else(|| format!("compound member {} is truncated", member.name))?;
            if numeric {
                values.push(decode_one(member.dtype, slice)?);
            } else {
                texts.extend(decode_fixed_strings(slice, 1, width));
            }
        }
        out.push(Dataset {
            path: format!("{}/{}", dataset.path, member.name),
            shape: dataset.shape.clone(),
            h5_type: if string {
                H5Type::Str
            } else {
                dataset.h5_type
            },
            data: values,
            strings: texts,
            vlen_data: Vec::new(),
            vlen_base_type: None,
        });
    }
    Ok(out)
}

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
        | DType::Enum { size, signed, be } => decode_one(DType::Int { size, signed, be }, bytes),
        | DType::VlenStr => Err("variable-length strings are decoded separately".to_owned()),
        | DType::VlenNumeric { .. } => {
            Err("variable-length numeric datasets are decoded separately".to_owned())
        },
        | DType::Other { class } => Err(format!("unsupported HDF5 datatype class {class}")),
    }
}

/// Decodes `count` variable-length numeric sequence descriptors.
fn decode_vlen_numeric(
    reader: &Reader<'_>,
    raw: &[u8],
    count: usize,
    base: DType,
) -> Result<Vec<Vec<f64>>, String> {
    let mut out = Vec::with_capacity(count);
    for index in 0..count {
        let base_off = index.saturating_mul(16);
        let descriptor = raw
            .get(base_off..base_off + 16)
            .ok_or_else(|| "vlen numeric descriptor is truncated".to_owned())?;
        let seq_len =
            u32::from_le_bytes([descriptor[0], descriptor[1], descriptor[2], descriptor[3]])
                as usize;
        let mut addr_bytes = [0_u8; 8];
        addr_bytes.copy_from_slice(&descriptor[4..12]);
        let addr = u64::from_le_bytes(addr_bytes);
        let mut index_bytes = [0_u8; 4];
        index_bytes.copy_from_slice(&descriptor[12..16]);
        let object = u32::from_le_bytes(index_bytes);
        if seq_len == 0 || addr == u64::MAX {
            out.push(Vec::new());
            continue;
        }
        let bytes = global_heap_object(reader, addr, object)?;
        let vals = decode_elements(base, &bytes, seq_len)?;
        out.push(vals);
    }
    Ok(out)
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
        | DType::Str { .. } | DType::VlenStr => Ok(H5Type::Str),
        | DType::VlenNumeric { .. } => Ok(H5Type::VlenNumeric),
        | DType::Enum { size, signed, .. } => {
            h5_type_of(DType::Int {
                size,
                signed,
                be: false,
            })
        },
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

/// Target of a group link entry.
#[derive(Debug, Clone)]
enum LinkTarget {
    Hard(u64),
    Soft(String),
}

/// Reads one symbol table node (`SNOD`) and appends its entries.
fn read_symbol_table_node(
    reader: &Reader<'_>,
    addr: u64,
    heap: (u64, u64),
    out: &mut Vec<(String, LinkTarget)>,
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
        let cache_type = reader.u32_at(offset + 2 * reader.off_size as u64)?;
        let name = heap_name(reader, heap.0, heap.1, name_offset)?;
        if cache_type == 2 {
            let target_offset = u64::from(reader.u32_at(offset + 2 * reader.off_size as u64 + 8)?);
            let target = heap_name(reader, heap.0, heap.1, target_offset)?;
            out.push((name, LinkTarget::Soft(target)));
        } else if object != u64::MAX {
            out.push((name, LinkTarget::Hard(object)));
        }
    }
    Ok(())
}

/// Recursively walks a version-1 group B-tree.
fn walk_group_btree(
    reader: &Reader<'_>,
    addr: u64,
    heap: (u64, u64),
    out: &mut Vec<(String, LinkTarget)>,
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

/// Parses a link message and returns its name and hard or soft target.
fn parse_link_message(
    reader: &Reader<'_>,
    data: &[u8],
) -> Result<Option<(String, LinkTarget)>, String> {
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
    if link_type == 1 {
        let target_len = to_usize(read_uint(data, pos, 2)?, "soft link target length")?;
        pos += 2;
        let target_bytes = data
            .get(pos..pos.saturating_add(target_len))
            .ok_or("truncated soft link target")?;
        let target = std::str::from_utf8(target_bytes)
            .map_err(|_| "soft link target is not valid UTF-8".to_owned())?
            .to_owned();
        return Ok(Some((name, LinkTarget::Soft(target))));
    }
    if link_type != 0 {
        return Ok(None);
    }
    let addr = read_uint(data, pos, reader.off_size)?;
    if addr == u64::MAX {
        return Ok(None);
    }
    Ok(Some((name, LinkTarget::Hard(addr))))
}

/// Collects the links declared by a group object header.
fn collect_links(
    reader: &Reader<'_>,
    messages: &[Message<'_>],
) -> Result<Vec<(String, LinkTarget)>, String> {
    let mut links = Vec::new();
    let mut symbol_table: Option<(u64, u64)> = None;
    let mut dense: Option<(u64, u64)> = None;
    for message in messages {
        match message.kind {
            | 0x02 => {
                let data = message.data;
                let flags = *data.get(1).ok_or("truncated link info message")?;
                let mut pos = 2usize;
                if flags & 0x01 != 0 {
                    // Maximum creation index, present when the group tracks
                    // creation order.
                    pos = pos.saturating_add(8);
                }
                let heap = read_uint(data, pos, reader.off_size)?;
                let btree = read_uint(data, pos + reader.off_size, reader.off_size)?;
                if heap != u64::MAX && btree != u64::MAX {
                    dense = Some((heap, btree));
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
    if let Some((heap_addr, btree_addr)) = dense {
        let mut records = Vec::new();
        walk_v2_btree(reader, btree_addr, &mut records)?;
        let debug = std::env::var_os("RSSN_H5_DEBUG").is_some();
        if debug {
            eprintln!(
                "dense group: {} records, heap {heap_addr:#x}",
                records.len()
            );
        }
        for id in records {
            let bytes = fractal_heap_object(reader, heap_addr, &id)?;
            match parse_link_message(reader, &bytes)? {
                | Some(link) => links.push(link),
                | None => {
                    if debug {
                        eprintln!(
                            "dropped link id {} bytes {}",
                            hex_encode(&id),
                            hex_encode(&bytes)
                        );
                    }
                },
            }
        }
    }
    Ok(links)
}

fn hex_encode(bytes: &[u8]) -> String {
    use std::fmt::Write as _;
    let mut hex = String::with_capacity(bytes.len().saturating_mul(2));
    for byte in bytes {
        let _ = write!(hex, "{byte:02x}");
    }
    hex
}

/// A fractal heap header (`FRHP`).
struct FractalHeap {
    /// Length of a heap identifier, in bytes.
    id_len: usize,
    /// Width of the object offset inside an identifier.
    off_size: usize,
    /// Width of the object length inside an identifier.
    len_size: usize,
    /// Entries per row of the doubling table.
    table_width: usize,
    /// Size of a first-row block.
    start_block_size: u64,
    /// Largest block a row may hold before it points at another indirect
    /// block.
    max_direct_block_size: u64,
    /// Width of the block offset stored in an indirect block header.
    block_offset_width: usize,
    /// File address of the root block.
    root_addr: u64,
    /// Rows in the root indirect block (zero means the root is a direct block).
    root_rows: usize,
}

/// Parses a fractal heap header.
fn parse_fractal_heap(
    reader: &Reader<'_>,
    addr: u64,
) -> Result<FractalHeap, String> {
    let head = reader.bytes(addr, 146)?;
    if head.get(..4) != Some(b"FRHP") {
        return Err("dense group storage is not a fractal heap".to_owned());
    }
    let u16_at = |at: usize| -> Result<u16, String> {
        let slice = head
            .get(at..at + 2)
            .ok_or("truncated fractal heap header")?;
        Ok(u16::from_le_bytes([slice[0], slice[1]]))
    };
    let u32_at = |at: usize| -> Result<u32, String> {
        let slice = head
            .get(at..at + 4)
            .ok_or("truncated fractal heap header")?;
        Ok(u32::from_le_bytes([slice[0], slice[1], slice[2], slice[3]]))
    };
    let u64_at = |at: usize| -> Result<u64, String> {
        let slice = head
            .get(at..at + 8)
            .ok_or("truncated fractal heap header")?;
        let mut raw = [0u8; 8];
        raw.copy_from_slice(slice);
        Ok(u64::from_le_bytes(raw))
    };
    let id_len = usize::from(u16_at(5)?);
    let max_object = u64::from(u32_at(10)?);
    let table_width = usize::from(u16_at(110)?);
    let start_block_size = u64_at(112)?;
    let max_direct_block_size = u64_at(120)?;
    let max_heap_bits = usize::from(u16_at(128)?);
    let root_addr = u64_at(132)?;
    let root_rows = usize::from(u16_at(140)?);
    if id_len < 3 || table_width == 0 {
        return Err("fractal heap header is inconsistent".to_owned());
    }
    // The identifier splits the remaining bytes between the object offset and
    // its length; both are sized from the header's limits. The block offset in
    // an indirect or direct block header is as wide as the heap's address
    // space, not as the file's offsets.
    let off_size = bytes_for_bits(max_heap_bits);
    let block_offset_width = off_size;
    let len_size = count_bytes(usize::try_from(max_object.min(max_direct_block_size)).unwrap_or(0));
    if 1 + off_size + len_size != id_len {
        return Err("fractal heap identifier length is inconsistent".to_owned());
    }
    Ok(FractalHeap {
        id_len,
        off_size,
        len_size,
        table_width,
        start_block_size,
        max_direct_block_size,
        block_offset_width,
        root_addr,
        root_rows,
    })
}

/// Resolves a heap identifier to the object bytes it names.
fn fractal_heap_object(
    reader: &Reader<'_>,
    heap_addr: u64,
    id: &[u8],
) -> Result<Vec<u8>, String> {
    let heap = parse_fractal_heap(reader, heap_addr)?;
    // Some indexes store extra fields after the identifier (the attribute name
    // index keeps the creation order and flags there), so a longer record is
    // trimmed to the identifier the heap uses.
    let id = if id.len() > heap.id_len {
        id.get(..heap.id_len).unwrap_or(id)
    } else {
        id
    };
    if id.len() != heap.id_len {
        return Err("fractal heap identifier has the wrong length".to_owned());
    }
    let kind = (id[0] >> 4) & 0x03;
    if kind != 0 {
        return Err(format!(
            "fractal heap identifier type {kind} (huge or tiny object) is not supported"
        ));
    }
    let offset = read_le(&id[1..=heap.off_size]);
    let length = read_le(&id[1 + heap.off_size..1 + heap.off_size + heap.len_size]);
    let length = usize::try_from(length).map_err(|_| "heap object is too large".to_owned())?;
    let (block_addr, block_start) = locate_heap_block(reader, &heap, offset)?;
    let within = offset
        .checked_sub(block_start)
        .ok_or_else(|| "heap object lies outside its block".to_owned())?;
    let start = block_addr
        .checked_add(within)
        .ok_or_else(|| "heap object address overflows".to_owned())?;
    reader
        .bytes(start, u64::try_from(length).unwrap_or(u64::MAX))
        .map(<[u8]>::to_vec)
}

/// Reads a little-endian unsigned integer of any width.
fn read_le(bytes: &[u8]) -> u64 {
    let mut value = 0u64;
    for (index, &byte) in bytes.iter().take(8).enumerate() {
        value |= u64::from(byte) << (8 * index);
    }
    value
}

/// Finds the direct block holding `offset` and the offset at which it starts.
///
/// The heap is a doubling table: rows zero and one hold blocks of the starting
/// size, and every later row doubles. A row past the maximum direct block size
/// points at another indirect block, which is descended into.
fn locate_heap_block(
    reader: &Reader<'_>,
    heap: &FractalHeap,
    offset: u64,
) -> Result<(u64, u64), String> {
    if heap.root_rows == 0 {
        // The root is a direct block that starts at offset zero.
        return Ok((heap.root_addr, 0));
    }
    locate_in_indirect(reader, heap, heap.root_addr, 0, heap.root_rows, offset)
}

/// Descends one indirect block of a fractal heap.
fn locate_in_indirect(
    reader: &Reader<'_>,
    heap: &FractalHeap,
    addr: u64,
    block_start: u64,
    rows: usize,
    offset: u64,
) -> Result<(u64, u64), String> {
    if rows == 0 {
        return Ok((addr, block_start));
    }
    let signature = reader.bytes(addr, 4)?;
    if signature != b"FHIB" {
        return Err("fractal heap block is not an indirect block".to_owned());
    }
    let header_len = 5usize
        .checked_add(reader.off_size)
        .and_then(|value| value.checked_add(heap.block_offset_width))
        .ok_or_else(|| "fractal heap header length overflows".to_owned())?;
    let direct_rows = direct_block_rows(heap, rows);
    let direct_entries = direct_rows.saturating_mul(heap.table_width);

    let mut row_start = block_start;
    let mut row_size = heap.start_block_size;
    for row in 0..rows {
        let row_end = row_start
            .checked_add(row_size.saturating_mul(heap.table_width as u64))
            .ok_or_else(|| "fractal heap offset overflows".to_owned())?;
        if offset < row_end {
            let column = usize::try_from((offset - row_start) / row_size.max(1)).unwrap_or(0);
            let child_start = row_start + row_size * u64::try_from(column).unwrap_or(0);
            let slot = if row < direct_rows {
                row * heap.table_width + column
            } else {
                direct_entries + (row - direct_rows) * heap.table_width + column
            };
            let entry = header_len
                .checked_add(slot.saturating_mul(8))
                .ok_or_else(|| "fractal heap entry offset overflows".to_owned())?;
            let header = reader.bytes(addr, u64::try_from(entry + 8).unwrap_or(u64::MAX))?;
            let child = read_le(
                header
                    .get(entry..entry + 8)
                    .ok_or_else(|| "fractal heap indirect block is truncated".to_owned())?,
            );
            if child == u64::MAX {
                return Err("fractal heap block is undefined".to_owned());
            }
            if row < direct_rows {
                return Ok((child, child_start));
            }
            // A nested indirect block covers its own slice of the heap and
            // repeats the doubling table from its own row count.
            let nested_rows = rows_for_block(heap, row_size);
            return locate_in_indirect(reader, heap, child, child_start, nested_rows, offset);
        }
        row_start = row_end;
        row_size = if row == 0 {
            heap.start_block_size
        } else {
            row_size.saturating_mul(2)
        };
    }
    Err("heap object offset lies outside the heap".to_owned())
}

/// Rows of a doubling table that point at direct blocks.
fn direct_block_rows(
    heap: &FractalHeap,
    rows: usize,
) -> usize {
    let start = heap.start_block_size.max(1);
    let max = heap.max_direct_block_size.max(start);
    // Rows zero and one hold the starting size, then each row doubles.
    let mut rows_for_max = 2usize;
    let mut size = start;
    while size < max && rows_for_max < 64 {
        size = size.saturating_mul(2);
        rows_for_max += 1;
    }
    rows.min(rows_for_max)
}

/// Rows in an indirect block that has to cover `block_size` heap bytes.
///
/// One row covers `start_block_size * table_width` bytes, and every extra row
/// doubles that.
fn rows_for_block(
    heap: &FractalHeap,
    block_size: u64,
) -> usize {
    let unit = heap
        .start_block_size
        .max(1)
        .saturating_mul(u64::try_from(heap.table_width.max(1)).unwrap_or(1));
    let mut rows = 1usize;
    let mut cover = unit;
    while cover < block_size && rows < 64 {
        cover = cover.saturating_mul(2);
        rows += 1;
    }
    rows
}

/// Bytes needed to hold a value of `bits` bits, or an integer of that size.
fn bytes_for_bits(bits: usize) -> usize {
    (bits / 8 + usize::from(!bits.is_multiple_of(8))).clamp(1, 32)
}

/// Bytes needed to store a count up to `value`.
fn count_bytes(value: usize) -> usize {
    let mut bytes = 1usize;
    while bytes < 8 && (1usize << (8 * bytes)) <= value {
        bytes += 1;
    }
    bytes
}

/// Field widths of the child pointers at each depth of a version-2 B-tree.
///
/// `counts[d]` is the width of the record count for a child at depth `d` and
/// `totals[d]` the width of its subtree total. The total is only stored when
/// the parent sits above the twig level.
struct BtreeWidths {
    /// Record-count widths, indexed by child depth.
    counts: Vec<usize>,
    /// Subtree-total widths, indexed by child depth.
    totals: Vec<usize>,
}

/// Computes the child-pointer field widths for a tree of `depth`.
fn btree_widths(
    depth: usize,
    node_size: usize,
    record_size: usize,
) -> BtreeWidths {
    let record_size = record_size.max(1);
    let mut counts = Vec::with_capacity(depth + 1);
    let mut totals = Vec::with_capacity(depth + 1);
    // A leaf holds records between its six-byte header and four-byte checksum.
    let mut capacity = node_size.saturating_sub(10) / record_size;
    counts.push(count_bytes(capacity));
    let mut subtree = capacity;
    totals.push(count_bytes(subtree));
    for level in 1..=depth {
        let pointer = 8
            + counts[level - 1]
            + if level > 1 {
                totals[level - 1]
            } else {
                0
            };
        capacity = node_size.saturating_sub(10 + pointer) / (record_size + pointer).max(1);
        counts.push(count_bytes(capacity));
        subtree = capacity.saturating_mul(1 + subtree);
        totals.push(count_bytes(subtree));
    }
    BtreeWidths { counts, totals }
}

/// Walks a version-2 B-tree and collects every record.
fn walk_v2_btree(
    reader: &Reader<'_>,
    addr: u64,
    out: &mut Vec<Vec<u8>>,
) -> Result<(), String> {
    let header = reader.bytes(addr, 6)?;
    match header.get(..4) {
        | Some(b"BTHD") => {
            let body = reader.bytes(addr, 26)?;
            let node_type = body[5];
            let node_size = usize::try_from(read_le(&body[6..10])).unwrap_or(512);
            let record_size = usize::from(u16::from_le_bytes([body[10], body[11]]));
            let tree_depth = usize::from(u16::from_le_bytes([body[12], body[13]]));
            let root = read_le(&body[16..24]);
            let root_count = usize::from(u16::from_le_bytes([body[24], body[25]]));
            if root == u64::MAX {
                return Ok(());
            }
            let widths = btree_widths(tree_depth, node_size, record_size);
            let first = out.len();
            walk_v2_node(
                reader,
                root,
                tree_depth,
                root_count,
                record_size,
                &widths,
                out,
            )?;
            // A link name index prefixes every record with the name hash; the
            // heap identifier follows it. An attribute name index stores the
            // identifier first and the hash last, so only type 5 is stripped.
            if node_type == 5 {
                for record in out.get_mut(first..).unwrap_or(&mut []) {
                    if record.len() > 4 {
                        record.drain(..4);
                    }
                }
            }
            Ok(())
        },
        | _ => {
            Err(format!(
                "group link index at {addr:#x} is not a version-2 B-tree"
            ))
        },
    }
}

/// Walks one version-2 B-tree node.
fn walk_v2_node(
    reader: &Reader<'_>,
    addr: u64,
    depth: usize,
    records: usize,
    record_size: usize,
    widths: &BtreeWidths,
    out: &mut Vec<Vec<u8>>,
) -> Result<(), String> {
    let signature = reader.bytes(addr, 4)?;
    match signature {
        | b"BTLF" => {
            let body = reader.bytes(
                addr,
                u64::try_from(6 + records * record_size).unwrap_or(u64::MAX),
            )?;
            for index in 0..records {
                let start = 6 + index * record_size;
                out.push(body[start..start + record_size].to_vec());
            }
            Ok(())
        },
        | b"BTIN" => {
            let count_width = *widths.counts.get(depth.saturating_sub(1)).unwrap_or(&1);
            let total_width = if depth > 1 {
                *widths.totals.get(depth.saturating_sub(1)).unwrap_or(&1)
            } else {
                0
            };
            let child_size = 8 + count_width + total_width;
            let total = 6 + records * record_size + (records + 1) * child_size;
            let body = reader.bytes(addr, u64::try_from(total).unwrap_or(u64::MAX))?;
            // A version-2 B-tree stores every record exactly once, so the
            // node's own records are part of the result.
            for index in 0..records {
                let start = 6 + index * record_size;
                out.push(body[start..start + record_size].to_vec());
            }
            // The separator records come first, then the child pointers.
            for index in 0..=records {
                let at = 6 + records * record_size + index * child_size;
                let child = read_le(&body[at..at + 8]);
                let count =
                    usize::try_from(read_le(&body[at + 8..at + 8 + count_width])).unwrap_or(0);
                if child == u64::MAX {
                    continue;
                }
                walk_v2_node(
                    reader,
                    child,
                    depth.saturating_sub(1),
                    count,
                    record_size,
                    widths,
                    out,
                )?;
            }
            Ok(())
        },
        | other => {
            Err(format!(
                "unexpected version-2 B-tree node signature {other:?}"
            ))
        },
    }
}

// -------------------------------------------------------------------------
// Attributes
// -------------------------------------------------------------------------

/// Parses an attribute info message, which points at the dense attribute
/// storage a group or dataset uses once it has many attributes.
///
/// The message is `version`, `flags`, then the maximum creation index when the
/// group tracks creation order, the fractal heap address, the name index
/// address, and the creation-order index address when it is indexed.
fn parse_attribute_info(
    reader: &Reader<'_>,
    data: &[u8],
) -> Result<Option<(u64, u64)>, String> {
    let version = *data.first().ok_or("empty attribute info message")?;
    if version != 0 {
        return Err(format!(
            "unsupported attribute info message version {version}"
        ));
    }
    let flags = *data.get(1).ok_or("truncated attribute info message")?;
    let mut pos = 2usize;
    if flags & 0x01 != 0 {
        pos += 2;
    }
    let heap = read_uint(data, pos, reader.off_size)?;
    let btree = read_uint(data, pos + reader.off_size, reader.off_size)?;
    if heap == u64::MAX || btree == u64::MAX {
        return Ok(None);
    }
    Ok(Some((heap, btree)))
}

/// Parses one attribute message.
fn parse_attribute(
    reader: &Reader<'_>,
    data: &[u8],
) -> Result<Attribute, String> {
    let version = *data.first().ok_or("empty attribute message")?;
    let (mut pos, name_padded) = match version {
        | 1 => (8_usize, true),
        | 2 => (8_usize, false),
        // Version 3 always carries the name character-set byte and never pads.
        | 3 => (9_usize, false),
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
    pos = advance(data, name_end, name_len, name_padded)?;
    let dtype_end = pos
        .checked_add(dtype_len)
        .ok_or_else(|| "attribute datatype length overflow".to_owned())?;
    let dtype_bytes = data
        .get(pos..dtype_end)
        .ok_or("truncated attribute datatype")?;
    pos = advance(data, dtype_end, dtype_len, name_padded)?;
    let dspace_end = pos
        .checked_add(dspace_len)
        .ok_or_else(|| "attribute dataspace length overflow".to_owned())?;
    let dspace_bytes = data
        .get(pos..dspace_end)
        .ok_or("truncated attribute dataspace")?;
    pos = advance(data, dspace_end, dspace_len, name_padded)?;
    let dtype = parse_datatype(dtype_bytes)?;
    let shape = parse_dataspace(dspace_bytes, reader.len_size)?;
    let count = element_count(&shape)?;
    let body = data.get(pos..).ok_or("truncated attribute data")?;
    let value = decode_attribute_value(reader, dtype, body, count)?;
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
    reader: &Reader<'_>,
    dtype: DType,
    body: &[u8],
    count: u64,
) -> Result<AttrValue, String> {
    let count = to_usize(count, "attribute element count")?;
    if dtype == DType::VlenStr {
        // A variable-length string is a descriptor into a global heap; the
        // text itself lives elsewhere in the file.
        let strings = decode_vlen_strings(reader, body, count)?;
        return Ok(AttrValue::Text(
            strings.into_iter().next().unwrap_or_default(),
        ));
    }
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
    let key_size = 8 + 8 * ndims;
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
                    offsets.push(read_uint(key, 8 + 8 * index, 8)?);
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

/// Reads and assembles a v3 chunked dataset into a contiguous raw byte buffer.
fn read_chunked_raw_v3(
    reader: &Reader<'_>,
    addr: u64,
    shape: &[usize],
    chunk_dims: &[usize],
    elem_size: usize,
    filters: &[Filter],
) -> Result<Vec<u8>, String> {
    let total = to_usize(element_count(shape)?, "dataset element count")?;
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
    assemble_chunks(
        reader,
        shape,
        &chunk_shape,
        elem_size,
        filters,
        &chunks,
        total,
    )
}

/// Assembles decoded chunks into a contiguous raw byte buffer of `total * elem_size` bytes.
fn assemble_chunks(
    reader: &Reader<'_>,
    shape: &[usize],
    chunk_shape: &[usize],
    elem_size: usize,
    filters: &[Filter],
    chunks: &[Chunk],
    total: usize,
) -> Result<Vec<u8>, String> {
    let rank = shape.len();
    let total_bytes = total
        .checked_mul(elem_size)
        .ok_or_else(|| "chunked dataset byte size overflows".to_owned())?;
    if total_bytes > MAX_ARRAY_BYTES {
        return Err("chunked dataset exceeds MAX_ARRAY_BYTES".to_owned());
    }
    let strides = strides_of(shape);
    let mut out = vec![0_u8; total_bytes];
    let chunk_elems = to_usize(element_count(chunk_shape)?, "chunk element count")?;
    let mut filter_decoder = FilterDecoder::default();
    for chunk in chunks {
        if chunk.addr == u64::MAX || chunk.size == 0 {
            continue;
        }
        let raw = reader.bytes(chunk.addr, chunk.size)?;
        let filtered = apply_filters(raw.to_vec(), filters, chunk.mask, &mut filter_decoder)?;
        for flat in 0..chunk_elems {
            let src0 = flat * elem_size;
            let Some(src_slice) = filtered.get(src0..src0 + elem_size) else {
                break;
            };
            let mut remainder = flat;
            let mut inside = true;
            let mut dst_elem = 0_usize;
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
                dst_elem += position * strides.get(axis).copied().unwrap_or(0);
            }
            if inside {
                let d0 = dst_elem * elem_size;
                if let Some(dst_slice) = out.get_mut(d0..d0 + elem_size) {
                    dst_slice.copy_from_slice(src_slice);
                }
            }
        }
    }
    Ok(out)
}

/// Computes the chunk origin coordinates for linear chunk index `idx` in row-major grid order.
fn chunk_origin_for_index(
    idx: usize,
    grid_counts: &[usize],
    chunk_shape: &[usize],
) -> Vec<u64> {
    let rank = grid_counts.len();
    let mut origins = vec![0_u64; rank];
    let mut rem = idx;
    for axis in (0..rank).rev() {
        let g = grid_counts[axis].max(1);
        origins[axis] = ((rem % g) * chunk_shape[axis]) as u64;
        rem /= g;
    }
    origins
}

/// Reads a version-4 chunked dataset (Single, Implicit, Fixed Array, Extensible Array, BTreeV2).
fn read_chunked_v4(
    reader: &Reader<'_>,
    addr: u64,
    shape: &[usize],
    chunk_dims: &[usize],
    elem_size: usize,
    filters: &[Filter],
    index: &V4ChunkIndex,
) -> Result<Vec<u8>, String> {
    let total = to_usize(element_count(shape)?, "dataset element count")?;
    let rank = shape.len();
    let chunk_dims = if chunk_dims.len() == rank + 1 {
        &chunk_dims[..rank]
    } else if chunk_dims.len() == rank {
        chunk_dims
    } else {
        return Err("v4 chunked layout dimension mismatch".to_owned());
    };
    let chunk_elems = to_usize(element_count(chunk_dims)?, "chunk element count")?;
    let raw_chunk_bytes = (chunk_elems * elem_size) as u64;
    let grid_counts: Vec<usize> = shape
        .iter()
        .zip(chunk_dims.iter())
        .map(|(&s, &c)| s.div_ceil(c.max(1)))
        .collect();
    let total_chunks: usize = grid_counts.iter().product();
    let mut chunks = Vec::new();

    match index {
        | V4ChunkIndex::Single {
            filtered_size,
            filter_mask,
        } => {
            if addr != u64::MAX {
                chunks.push(Chunk {
                    offsets: vec![0; rank],
                    addr,
                    size: filtered_size.unwrap_or(raw_chunk_bytes),
                    mask: *filter_mask,
                });
            }
        },
        | V4ChunkIndex::Implicit => {
            if addr != u64::MAX {
                for idx in 0..total_chunks {
                    chunks.push(Chunk {
                        offsets: chunk_origin_for_index(idx, &grid_counts, chunk_dims),
                        addr: addr + (idx as u64) * raw_chunk_bytes,
                        size: raw_chunk_bytes,
                        mask: 0,
                    });
                }
            }
        },
        | V4ChunkIndex::FixedArray { page_bits } => {
            if addr != u64::MAX {
                // FAHD: "FAHD"(4) + ver(1) + client_id(1) + entry_size(1) + page_bits(1) +
                // max_nelmts(len_size) + data_blk_addr(off_size) + cksum(4)
                let fahd =
                    reader.bytes(addr, 8 + reader.len_size as u64 + reader.off_size as u64)?;
                if fahd.get(..4) != Some(b"FAHD") {
                    return Err("invalid Fixed Array header signature".to_owned());
                }
                let client_id = fahd[5];
                let entry_size = usize::from(fahd[6]);
                let max_nelmts = to_usize(read_uint(fahd, 8, reader.len_size)?, "FAHD max_nelmts")?;
                let dblk_addr = read_uint(fahd, 8 + reader.len_size, reader.off_size)?;
                if dblk_addr != u64::MAX {
                    // FADB: "FADB"(4) + ver(1) + client_id(1) + header_addr(off_size)
                    let fadb_hdr_len = 6 + reader.off_size;
                    let page_size = 1_usize << (*page_bits).min(30);
                    let is_paged = max_nelmts > page_size;
                    if !is_paged {
                        let fadb = reader
                            .bytes(dblk_addr, (fadb_hdr_len + max_nelmts * entry_size) as u64)?;
                        if fadb.get(..4) != Some(b"FADB") {
                            return Err("invalid Fixed Array data block signature".to_owned());
                        }
                        for idx in 0..max_nelmts.min(total_chunks) {
                            let ep = fadb_hdr_len + idx * entry_size;
                            let caddr = read_uint(fadb, ep, reader.off_size)?;
                            if caddr == u64::MAX {
                                continue;
                            }
                            let (size, mask) = if client_id == 1 {
                                let fs = read_uint(
                                    fadb,
                                    ep + reader.off_size,
                                    entry_size - reader.off_size - 4,
                                )?;
                                let fm = u32::try_from(read_uint(fadb, ep + entry_size - 4, 4)?)
                                    .unwrap_or(0);
                                (fs, fm)
                            } else {
                                (raw_chunk_bytes, 0)
                            };
                            chunks.push(Chunk {
                                offsets: chunk_origin_for_index(idx, &grid_counts, chunk_dims),
                                addr: caddr,
                                size,
                                mask,
                            });
                        }
                    } else {
                        // Paged FADB: bitmap of ceil(npages/8) bytes after fadb_hdr, then cksum(4),
                        // then pages of `page_size * entry_size + 4` bytes.
                        let npages = max_nelmts.div_ceil(page_size);
                        let bitmap_bytes = npages.div_ceil(8);
                        let fadb_hdr =
                            reader.bytes(dblk_addr, (fadb_hdr_len + bitmap_bytes + 4) as u64)?;
                        let mut page_off = dblk_addr + (fadb_hdr_len + bitmap_bytes + 4) as u64;
                        for p in 0..npages {
                            let nelmts_in_page = if p + 1 == npages {
                                max_nelmts - p * page_size
                            } else {
                                page_size
                            };
                            let present = (fadb_hdr[fadb_hdr_len + p / 8] >> (p % 8)) & 1 != 0;
                            let pbytes = (nelmts_in_page * entry_size + 4) as u64;
                            if present {
                                let page_buf = reader.bytes(page_off, pbytes)?;
                                for i in 0..nelmts_in_page {
                                    let idx = p * page_size + i;
                                    if idx >= total_chunks {
                                        break;
                                    }
                                    let ep = i * entry_size;
                                    let caddr = read_uint(page_buf, ep, reader.off_size)?;
                                    if caddr == u64::MAX {
                                        continue;
                                    }
                                    let (size, mask) = if client_id == 1 {
                                        let fs = read_uint(
                                            page_buf,
                                            ep + reader.off_size,
                                            entry_size - reader.off_size - 4,
                                        )?;
                                        let fm = u32::try_from(read_uint(
                                            page_buf,
                                            ep + entry_size - 4,
                                            4,
                                        )?)
                                        .unwrap_or(0);
                                        (fs, fm)
                                    } else {
                                        (raw_chunk_bytes, 0)
                                    };
                                    chunks.push(Chunk {
                                        offsets: chunk_origin_for_index(
                                            idx,
                                            &grid_counts,
                                            chunk_dims,
                                        ),
                                        addr: caddr,
                                        size,
                                        mask,
                                    });
                                }
                            }
                            page_off += pbytes;
                        }
                    }
                }
            }
        },
        | V4ChunkIndex::ExtensibleArray { .. } => {
            if addr != u64::MAX {
                let eahd = reader.bytes(
                    addr,
                    16 + 6 * reader.len_size as u64 + reader.off_size as u64,
                )?;
                if eahd.get(..4) != Some(b"EAHD") {
                    return Err("invalid Extensible Array header signature".to_owned());
                }
                let client_id = eahd[5];
                let entry_size = usize::from(eahd[6]);
                let idx_blk_elmts = usize::from(eahd[8]);
                let mut pos = 12 + 4 * reader.len_size;
                let max_idx_set =
                    to_usize(read_uint(eahd, pos, reader.len_size)?, "EAHD max_idx_set")?;
                pos += reader.len_size;
                let _num_elmts = read_uint(eahd, pos, reader.len_size)?;
                pos += reader.len_size;
                let idx_blk_addr = read_uint(eahd, pos, reader.off_size)?;
                if idx_blk_addr != u64::MAX && max_idx_set > 0 {
                    // EAIB: "EAIB"(4) + ver(1) + client_id(1) + header_addr(off_size)
                    let eaib_hdr_len = 6 + reader.off_size;
                    let eaib = reader.bytes(
                        idx_blk_addr,
                        (eaib_hdr_len + idx_blk_elmts * entry_size) as u64,
                    )?;
                    if eaib.get(..4) != Some(b"EAIB") {
                        return Err("invalid Extensible Array index block signature".to_owned());
                    }
                    for idx in 0..idx_blk_elmts.min(total_chunks) {
                        let ep = eaib_hdr_len + idx * entry_size;
                        let caddr = read_uint(eaib, ep, reader.off_size)?;
                        if caddr == u64::MAX {
                            continue;
                        }
                        let (size, mask) = if client_id == 1 {
                            let fs = read_uint(
                                eaib,
                                ep + reader.off_size,
                                entry_size - reader.off_size - 4,
                            )?;
                            let fm = u32::try_from(read_uint(eaib, ep + entry_size - 4, 4)?)
                                .unwrap_or(0);
                            (fs, fm)
                        } else {
                            (raw_chunk_bytes, 0)
                        };
                        chunks.push(Chunk {
                            offsets: chunk_origin_for_index(idx, &grid_counts, chunk_dims),
                            addr: caddr,
                            size,
                            mask,
                        });
                    }
                }
            }
        },
        | V4ChunkIndex::BTreeV2 { .. } => {
            if addr != u64::MAX {
                let mut records = Vec::new();
                walk_v2_btree(reader, addr, &mut records)?;
                let bthd = reader.bytes(addr, 6)?;
                let node_type = bthd[5];
                for rec in records {
                    let caddr = read_uint(&rec, 0, reader.off_size)?;
                    if caddr == u64::MAX {
                        continue;
                    }
                    let mut rpos = reader.off_size;
                    let (size, mask) = if node_type == 11 {
                        // Filtered chunk: addr(off_size) + chunk_size(len_size? actually 8-len or var) + filter_mask(4) + coords(rank * 8)
                        let coord_bytes = rank * 8;
                        let size_len = rec.len().saturating_sub(reader.off_size + 4 + coord_bytes);
                        let fs = read_uint(&rec, rpos, size_len)?;
                        rpos += size_len;
                        let fm = u32::try_from(read_uint(&rec, rpos, 4)?).unwrap_or(0);
                        rpos += 4;
                        (fs, fm)
                    } else {
                        (raw_chunk_bytes, 0)
                    };
                    let mut offsets = Vec::with_capacity(rank);
                    for d in 0..rank {
                        let off = read_uint(&rec, rpos + d * 8, 8).unwrap_or(0);
                        offsets.push(off.saturating_mul(chunk_dims[d] as u64));
                    }
                    chunks.push(Chunk {
                        offsets,
                        addr: caddr,
                        size,
                        mask,
                    });
                }
            }
        },
    }

    assemble_chunks(
        reader, shape, chunk_dims, elem_size, filters, &chunks, total,
    )
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
    let total = to_usize(element_count(shape)?, "dataset element count")?;
    let elem_size = dtype.size();
    if elem_size == 0 {
        return Err("chunked dataset uses an unsupported datatype".to_owned());
    }
    let raw = read_chunked_raw_v3(reader, addr, shape, chunk_dims, elem_size, filters)?;
    decode_elements(dtype, &raw, total)
}

// -------------------------------------------------------------------------
// Datasets
// -------------------------------------------------------------------------

/// Builds a [`Dataset`] from an object header's messages.
fn read_dataset(
    reader: &Reader<'_>,
    messages: &[Message<'_>],
    path: &str,
) -> Result<Vec<Dataset>, String> {
    let mut shape = Vec::new();
    let mut dtype = None;
    let mut dtype_bytes: Option<&[u8]> = None;
    let mut layout = None;
    let mut filters = Vec::new();
    for message in messages {
        match message.kind {
            | 0x01 => shape = parse_dataspace(message.data, reader.len_size)?,
            | 0x03 => {
                dtype = Some(parse_datatype(message.data)?);
                dtype_bytes = Some(message.data);
            },
            | 0x08 => layout = Some(parse_layout(reader, message.data)?),
            | 0x0B => filters = parse_filters(message.data)?,
            | _ => {},
        }
    }
    let (Some(dtype), Some(layout)) = (dtype, layout) else {
        // A dataset whose datatype or layout message is shared or absent
        // cannot be represented; skipping it keeps the rest of the file.
        return Ok(Vec::new());
    };
    // A compound dataset is exposed as one dataset per representable member.
    if matches!(dtype, DType::Other { class: 6 }) {
        let Some(bytes) = dtype_bytes else {
            return Ok(Vec::new());
        };
        // A compound element is as wide as its datatype message declares.
        let total_size = to_usize(read_uint(bytes, 4, 4)?, "compound size")?;
        let count = element_count(&shape)?;
        let count = to_usize(count, "dataset element count")?;
        let body = match layout {
            | Layout::Compact(body) => body,
            | Layout::Contiguous { addr, size } => {
                if addr == u64::MAX {
                    Vec::new()
                } else {
                    let derived = u64::try_from(count)
                        .unwrap_or(0)
                        .saturating_mul(u64::try_from(total_size).unwrap_or(0));
                    reader.bytes(addr, size.max(derived))?.to_vec()
                }
            },
            | Layout::Chunked { addr, dims } => {
                if addr == u64::MAX {
                    Vec::new()
                } else {
                    read_chunked_raw_v3(reader, addr, &shape, &dims, total_size, &filters)?
                }
            },
            | Layout::ChunkedV4 {
                addr, dims, index, ..
            } => {
                if addr == u64::MAX {
                    Vec::new()
                } else {
                    read_chunked_v4(reader, addr, &shape, &dims, total_size, &filters, &index)?
                }
            },
        };
        let template = Dataset::new_numeric(path, shape, H5Type::F64, Vec::new());
        return read_compound_dataset(&template, bytes, &body, count, total_size);
    }
    let Ok(h5_type) = h5_type_of(dtype) else {
        // Array, reference and other exotic classes are skipped so a single one
        // of them cannot hide the rest of the file.
        return Ok(Vec::new());
    };
    let count = element_count(&shape)?;
    let count = to_usize(count, "dataset element count")?;
    if h5_type == H5Type::Str {
        let raw = match layout {
            | Layout::Compact(body) => body,
            | Layout::Contiguous { addr, size } => {
                if addr == u64::MAX {
                    Vec::new()
                } else {
                    let derived = u64::try_from(count)
                        .unwrap_or(0)
                        .saturating_mul(u64::try_from(dtype.size()).unwrap_or(0));
                    reader.bytes(addr, size.max(derived))?.to_vec()
                }
            },
            | Layout::Chunked { addr, dims } => {
                if addr == u64::MAX {
                    Vec::new()
                } else {
                    read_chunked_raw_v3(reader, addr, &shape, &dims, dtype.size(), &filters)?
                }
            },
            | Layout::ChunkedV4 {
                addr, dims, index, ..
            } => {
                if addr == u64::MAX {
                    Vec::new()
                } else {
                    read_chunked_v4(reader, addr, &shape, &dims, dtype.size(), &filters, &index)?
                }
            },
        };
        let strings = match dtype {
            | DType::VlenStr => decode_vlen_strings(reader, &raw, count)?,
            | DType::Str { size } => decode_fixed_strings(&raw, count, size),
            | _ => Vec::new(),
        };
        return Ok(vec![Dataset::new_strings(path, shape, strings)]);
    }
    if let DType::VlenNumeric {
        size,
        signed,
        be,
        is_float,
    } = dtype
    {
        let raw = match layout {
            | Layout::Compact(body) => body,
            | Layout::Contiguous { addr, size: bsize } => {
                if addr == u64::MAX {
                    Vec::new()
                } else {
                    let derived = u64::try_from(count).unwrap_or(0).saturating_mul(16);
                    reader.bytes(addr, bsize.max(derived))?.to_vec()
                }
            },
            | Layout::Chunked { addr, dims } => {
                if addr == u64::MAX {
                    Vec::new()
                } else {
                    read_chunked_raw_v3(reader, addr, &shape, &dims, 16, &filters)?
                }
            },
            | Layout::ChunkedV4 {
                addr, dims, index, ..
            } => {
                if addr == u64::MAX {
                    Vec::new()
                } else {
                    read_chunked_v4(reader, addr, &shape, &dims, 16, &filters, &index)?
                }
            },
        };
        let base_dtype = if is_float {
            DType::Float { size, be }
        } else {
            DType::Int { size, signed, be }
        };
        let base_h5 = h5_type_of(base_dtype)?;
        let vlen_data = if raw.is_empty() {
            vec![Vec::new(); count]
        } else {
            decode_vlen_numeric(reader, &raw, count, base_dtype)?
        };
        return Ok(vec![Dataset::new_vlen_numeric(
            path, shape, base_h5, vlen_data,
        )]);
    }
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
        | Layout::ChunkedV4 {
            addr, dims, index, ..
        } => {
            if addr == u64::MAX {
                vec![0.0_f64; count]
            } else {
                let raw =
                    read_chunked_v4(reader, addr, &shape, &dims, dtype.size(), &filters, &index)?;
                decode_elements(dtype, &raw, count)?
            }
        },
    };
    Ok(vec![Dataset::new_numeric(path, shape, h5_type, data)])
}

/// Splits fixed-width, NUL-padded strings into text.
fn decode_fixed_strings(
    raw: &[u8],
    count: usize,
    size: usize,
) -> Vec<String> {
    let size = size.max(1);
    (0..count)
        .map(|index| {
            let start = index.saturating_mul(size);
            let slice = raw.get(start..start.saturating_add(size)).unwrap_or(&[]);
            let end = slice
                .iter()
                .position(|&byte| byte == 0)
                .unwrap_or(slice.len());
            String::from_utf8_lossy(slice.get(..end).unwrap_or(&[])).into_owned()
        })
        .collect()
}

/// Recursively walks one object, appending groups, datasets and attributes.
fn walk_object(
    reader: &Reader<'_>,
    addr: u64,
    path: &str,
    out: &mut H5File,
    visited: &mut std::collections::HashSet<u64>,
    depth: usize,
) -> Result<(), String> {
    if depth > MAX_DEPTH {
        return Err(format!(
            "group nesting at {path} is deeper than {MAX_DEPTH}"
        ));
    }
    if !visited.insert(addr) {
        // A group of thousands of links to one object must not re-walk it.
        return Ok(());
    }
    let messages = read_messages(reader, addr)?;
    let mut attributes = Vec::new();
    let mut dense_attributes: Option<(u64, u64)> = None;
    for message in &messages {
        match message.kind {
            | 0x0C => attributes.push(parse_attribute(reader, message.data)?),
            | 0x15 => dense_attributes = parse_attribute_info(reader, message.data)?,
            | _ => {},
        }
    }
    if let Some((heap_addr, btree_addr)) = dense_attributes {
        let mut records = Vec::new();
        walk_v2_btree(reader, btree_addr, &mut records)?;
        let debug = std::env::var_os("RSSN_H5_DEBUG").is_some();
        for id in records {
            let bytes = fractal_heap_object(reader, heap_addr, &id)?;
            match parse_attribute(reader, &bytes) {
                | Ok(attribute) => attributes.push(attribute),
                | Err(error) => {
                    if debug {
                        eprintln!("dropped attribute at {path}: {error}");
                    }
                },
            }
        }
    }
    if !attributes.is_empty() {
        out.attributes.push((path.to_owned(), attributes));
    }
    let is_dataset = messages.iter().any(|message| message.kind == 0x01)
        && messages.iter().any(|message| message.kind == 0x03)
        && messages.iter().any(|message| message.kind == 0x08);
    if is_dataset {
        out.datasets.extend(read_dataset(reader, &messages, path)?);
        return Ok(());
    }
    out.groups.push(path.to_owned());
    let mut links = collect_links(reader, &messages)?;
    links.sort_by(|left, right| left.0.cmp(&right.0));
    for (name, target) in links {
        let child_path = if path == "/" {
            format!("/{name}")
        } else {
            format!("{path}/{name}")
        };
        match target {
            | LinkTarget::Hard(child) => {
                walk_object(reader, child, &child_path, out, visited, depth + 1)?;
            },
            | LinkTarget::Soft(raw_target) => {
                let resolved = resolve_soft_target(path, &raw_target);
                out.soft_links.push((child_path, resolved));
            },
        }
    }
    Ok(())
}

/// Resolves a soft link target path (absolute or relative to `parent_group`).
fn resolve_soft_target(
    parent_group: &str,
    raw_target: &str,
) -> String {
    let joined = if raw_target.starts_with('/') {
        raw_target.to_owned()
    } else if parent_group == "/" {
        format!("/{raw_target}")
    } else {
        format!("{parent_group}/{raw_target}")
    };
    let mut parts: Vec<&str> = Vec::new();
    for seg in joined.split('/') {
        match seg {
            | "" | "." => {},
            | ".." => {
                parts.pop();
            },
            | other => parts.push(other),
        }
    }
    if parts.is_empty() {
        "/".to_owned()
    } else {
        format!("/{}", parts.join("/"))
    }
}

/// Expands soft links in `file` so datasets/groups reachable through soft links
/// are also accessible at the soft-link path (with cycle protection).
fn resolve_soft_links_in_file(file: &mut H5File) {
    let mut extra_datasets = Vec::new();
    for (link_path, _) in &file.soft_links {
        let mut cur = link_path.as_str();
        for _ in 0..MAX_DEPTH {
            if let Some(ds) = file.datasets.iter().find(|d| d.path == cur) {
                if !file.datasets.iter().any(|d| d.path == *link_path) {
                    let mut cloned = ds.clone();
                    cloned.path.clone_from(link_path);
                    extra_datasets.push(cloned);
                }
                break;
            }
            if let Some(next) = file.resolve_soft_link(cur) {
                cur = next;
            } else {
                break;
            }
        }
    }
    file.datasets.extend(extra_datasets);
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
/// feature, or contains inconsistent offsets.
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
    let mut visited = std::collections::HashSet::new();
    walk_object(&reader, superblock.root, "/", &mut file, &mut visited, 0)?;
    resolve_soft_links_in_file(&mut file);
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
        | H5Type::VlenNumeric => 16,
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
            // Bit offset of the sign bit: the top bit of the value.
            out.push(u8::try_from(size * 8 - 1).unwrap_or(0));
            out.push(0);
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
        | H5Type::Str | H5Type::VlenNumeric => {
            return Err("string/vlen datasets use dedicated datatype builders".to_owned());
        },
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
    // `name_size` counts the NUL terminator.
    let stored_name_len = name_len.saturating_add(1);
    if stored_name_len > u16::MAX as usize
        || dtype_len > u16::MAX as usize
        || dspace_len > u16::MAX as usize
    {
        return Err(format!("attribute {} is too large", attribute.name));
    }
    let mut out = Vec::new();
    out.push(3);
    out.push(0);
    out.extend_from_slice(&u16::try_from(stored_name_len).unwrap_or(0).to_le_bytes());
    out.extend_from_slice(&u16::try_from(dtype_len).unwrap_or(0).to_le_bytes());
    out.extend_from_slice(&u16::try_from(dspace_len).unwrap_or(0).to_le_bytes());
    // Version 3 header: name character set (0 = ASCII/UTF-8), then the
    // unpadded name, datatype and dataspace.
    out.push(0);
    out.extend_from_slice(name_bytes);
    out.push(0);
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
        | H5Type::Str | H5Type::VlenNumeric => {
            return Err("string/vlen datasets cannot be encoded as flat numeric data".to_owned());
        },
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
        // A version-1 header stores the *padded* size: the reference library
        // rejects a message whose size is not a multiple of eight.
        let padding = (8 - data.len() % 8) % 8;
        let size = u16::try_from(data.len().saturating_add(padding))
            .map_err(|_| "object header message is larger than 64 KiB".to_owned())?;
        body.extend_from_slice(&kind.to_le_bytes());
        body.extend_from_slice(&size.to_le_bytes());
        body.push(0);
        body.extend_from_slice(&[0; 3]);
        body.extend_from_slice(data);
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
        .map(|(_, data)| 8 + data.len().div_ceil(8) * 8)
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

// Remove unused single_name warning

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
        (0x03, dataset_datatype_message(dataset)?),
        (0x05, fill_value_message()),
        (0x08, layout_message(data_addr, data_len)),
    ];
    for attribute in attributes {
        messages.push((0x0C, attribute_message(attribute)?));
    }
    Ok(messages)
}

/// Byte width of one element of a dataset.
///
/// A string dataset stores fixed-width `S<n>` elements, so its width is the
/// longest string in it.
fn element_size(dataset: &Dataset) -> usize {
    if dataset.h5_type == H5Type::Str {
        string_width(dataset)
    } else {
        type_size(dataset.h5_type)
    }
}

/// The fixed width used for the strings of a dataset.
fn string_width(dataset: &Dataset) -> usize {
    dataset
        .strings
        .iter()
        .map(String::len)
        .max()
        .unwrap_or(0)
        .max(1)
}

/// The datatype message of a dataset, which carries the string width.
fn dataset_datatype_message(dataset: &Dataset) -> Result<Vec<u8>, String> {
    if dataset.h5_type == H5Type::Str {
        Ok(string_datatype_message(string_width(dataset)))
    } else {
        datatype_message(dataset.h5_type)
    }
}

/// Encodes a dataset body, padding strings to their fixed width.
fn encode_dataset_data(dataset: &Dataset) -> Result<Vec<u8>, String> {
    if dataset.h5_type != H5Type::Str {
        return encode_data(dataset.h5_type, &dataset.data);
    }
    let width = string_width(dataset);
    let mut out = Vec::with_capacity(dataset.strings.len().saturating_mul(width));
    for value in &dataset.strings {
        let bytes = value.as_bytes();
        if bytes.len() > width {
            return Err(format!(
                "string {value:?} does not fit the dataset width {width}"
            ));
        }
        out.extend_from_slice(bytes);
        out.extend(std::iter::repeat_n(0_u8, width - bytes.len()));
    }
    Ok(out)
}

/// Serialises a local heap holding `names` and returns it with each offset.
fn local_heap(names: &[String]) -> (Vec<u8>, Vec<u64>) {
    // Offset zero holds the empty name: a B-tree key pointing there sorts
    // before every real name, which is how the first child's left boundary is
    // written.
    let mut data = vec![0_u8];
    let mut offsets = Vec::with_capacity(names.len());
    for name in names {
        offsets.push(u64::try_from(data.len()).unwrap_or(0));
        data.extend_from_slice(name.as_bytes());
        data.push(0);
    }
    // The data segment is at least one aligned block, even for a group with
    // no children: the free-list head lives inside it.
    while data.len() < 8 || data.len() % 8 != 0 {
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

/// Capacity of one symbol table node: two key values per child, four per node.
const SNOD_CAPACITY: usize = 8;
/// On-disk size of a symbol table node: an eight-byte header and 8 entries.
const SNOD_SIZE: usize = 8 + SNOD_CAPACITY * 40;
/// Children per version-1 B-tree node.
const BTREE_CHILDREN: usize = 8;
/// On-disk size of a version-1 group B-tree node: header, children, keys.
const BTREE_NODE_SIZE: usize = 24 + BTREE_CHILDREN * 8 + (BTREE_CHILDREN + 1) * 8;

/// What a B-tree node points at.
#[derive(Clone, Copy)]
enum GroupChild {
    /// Another B-tree node, by index.
    Node(usize),
    /// A symbol table node, by index.
    Snod(usize),
}

/// One version-1 group B-tree node.
struct GroupNode {
    /// Depth below the root; leaves are zero.
    level: u8,
    /// `children.len() + 1` boundary keys: local heap offsets of link names.
    keys: Vec<u64>,
    /// Child nodes or symbol table nodes, in name order.
    children: Vec<GroupChild>,
}

/// Builds the B-tree over `entries`, one symbol table node per eight links.
///
/// The nodes come back with children before parents, so an index always refers
/// to an already-built node.
fn build_group_nodes(entries: &[(String, u64, u64)]) -> (Vec<GroupNode>, usize) {
    let mut nodes: Vec<GroupNode> = entries
        .chunks(SNOD_CAPACITY)
        .enumerate()
        .map(|(index, chunk)| {
            GroupNode {
                level: 0,
                keys: vec![0, chunk[chunk.len() - 1].1],
                children: vec![GroupChild::Snod(index)],
            }
        })
        .collect();
    if nodes.is_empty() {
        return (
            vec![GroupNode {
                level: 0,
                keys: vec![0],
                children: Vec::new(),
            }],
            0,
        );
    }
    let mut current: Vec<usize> = (0..nodes.len()).collect();
    let mut level = 1u8;
    while current.len() > 1 {
        let mut next = Vec::new();
        for group in current.chunks(BTREE_CHILDREN) {
            // A key is the last name of the child before it; the first key is
            // the empty name so the first child is reachable.
            let mut keys: Vec<u64> = vec![0];
            for &index in group {
                keys.push(nodes[index].keys[1]);
            }
            let index = nodes.len();
            nodes.push(GroupNode {
                level,
                keys,
                children: group.iter().map(|&child| GroupChild::Node(child)).collect(),
            });
            next.push(index);
        }
        current = next;
        level = level.saturating_add(1);
    }
    let root = current.first().copied().unwrap_or(0);
    (nodes, root)
}

/// Serialises every group B-tree node, resolving child addresses.
fn group_btree(
    nodes: &[GroupNode],
    root: usize,
    node_base: u64,
    snod_base: u64,
) -> Vec<u8> {
    let stride = u64::try_from(BTREE_NODE_SIZE).unwrap_or(1);
    let snod_stride = u64::try_from(SNOD_SIZE).unwrap_or(1);
    // The address the group header points at is the root, so the root is
    // written first and every other node follows.
    let mut order: Vec<usize> = Vec::with_capacity(nodes.len());
    order.push(root);
    for index in 0..nodes.len() {
        if index != root {
            order.push(index);
        }
    }
    let address_of = |index: usize| -> u64 {
        let position = order.iter().position(|&value| value == index).unwrap_or(0);
        node_base + u64::try_from(position).unwrap_or(0) * stride
    };
    let mut out = Vec::with_capacity(nodes.len().saturating_mul(BTREE_NODE_SIZE));
    for &index in &order {
        let node = &nodes[index];
        // Siblings at the same level are linked left to right.
        let siblings: Vec<usize> = nodes
            .iter()
            .enumerate()
            .filter(|(_, other)| other.level == node.level)
            .map(|(position, _)| position)
            .collect();
        let position = siblings
            .iter()
            .position(|&value| value == index)
            .unwrap_or(0);
        let left = position
            .checked_sub(1)
            .and_then(|previous| siblings.get(previous))
            .map_or(u64::MAX, |&value| address_of(value));
        let right = siblings
            .get(position + 1)
            .map_or(u64::MAX, |&value| address_of(value));
        let mut buf = Vec::with_capacity(BTREE_NODE_SIZE);
        buf.extend_from_slice(b"TREE");
        buf.push(0);
        buf.push(node.level);
        buf.extend_from_slice(
            &u16::try_from(node.children.len())
                .unwrap_or(0)
                .to_le_bytes(),
        );
        buf.extend_from_slice(&left.to_le_bytes());
        buf.extend_from_slice(&right.to_le_bytes());
        for (slot, child) in node.children.iter().enumerate() {
            buf.extend_from_slice(&node.keys[slot].to_le_bytes());
            let address = match child {
                | GroupChild::Node(target) => address_of(*target),
                | GroupChild::Snod(target) => {
                    snod_base + u64::try_from(*target).unwrap_or(0) * snod_stride
                },
            };
            buf.extend_from_slice(&address.to_le_bytes());
        }
        let boundary = node.keys.get(node.children.len()).copied().unwrap_or(0);
        buf.extend_from_slice(&boundary.to_le_bytes());
        // The reference library reads a fixed-size node, so pad to it.
        buf.resize(BTREE_NODE_SIZE, 0);
        out.extend_from_slice(&buf);
    }
    out
}

/// Serialises one symbol table node, padded to its fixed size.
fn symbol_table_node(entries: &[(String, u64, u64)]) -> Vec<u8> {
    let mut out = Vec::with_capacity(SNOD_SIZE);
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
    out.resize(SNOD_SIZE, 0);
    out
}

/// Serialises a version-0 superblock whose root group entry points at the
/// group B-tree and local heap.
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
/// Dataset paths may be nested (`/a/b/data`); every intermediate group is
/// created. Attributes are attached to the object whose path they name.
///
/// # Errors
/// Fails when two objects share a path, a path is relative, an integer value
/// does not fit its declared type, or the file cannot be written.
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
/// Splits a dataset path into its segments, rejecting empty or relative ones.
fn path_segments(path: &str) -> Result<Vec<String>, String> {
    if !path.starts_with('/') {
        return Err(format!("dataset path `{path}` must start with `/`"));
    }
    let segments: Vec<String> = path
        .split('/')
        .filter(|segment| !segment.is_empty())
        .map(str::to_owned)
        .collect();
    if segments.is_empty() {
        return Err("dataset path must name an object".to_owned());
    }
    Ok(segments)
}

/// What one group links to.
#[derive(Clone, Copy)]
enum PlanChild {
    /// A nested group, by plan index.
    Group(usize),
    /// A dataset, by index in the caller's slice.
    Dataset(usize),
}

/// Everything needed to lay out and write one group.
struct GroupPlan {
    /// Absolute path of the group.
    path: String,
    /// Children in name order.
    children: Vec<PlanChild>,
    /// Symbol table entries: name, local-heap offset, child header address.
    entries: Vec<(String, u64, u64)>,
    /// Local heap bytes.
    heap: Vec<u8>,
    /// Group B-tree nodes.
    nodes: Vec<GroupNode>,
    /// Index of the root node in `nodes`.
    root_index: usize,
    /// Number of symbol table nodes.
    snod_count: usize,
    /// Length of the object header.
    header_len: usize,
    /// Addresses assigned during layout.
    header_addr: u64,
    /// Address of the local heap.
    heap_addr: u64,
    /// Address of the group B-tree.
    btree_addr: u64,
    /// Address of the first symbol table node.
    snod_addr: u64,
}

/// Serialises a complete HDF5 file.
///
/// Datasets may live in nested groups: every intermediate group gets its own
/// object header, local heap, group B-tree and symbol table nodes.
pub fn write_hdf5_bytes(
    datasets: &[Dataset],
    attributes: &[(String, Vec<Attribute>)],
) -> Result<Vec<u8>, String> {
    let mut encoded = Vec::with_capacity(datasets.len());
    let mut lengths = Vec::with_capacity(datasets.len());
    for dataset in datasets {
        let count = element_count(&dataset.shape)?;
        let count = to_usize(count, "dataset element count")?;
        if dataset.h5_type == H5Type::Str {
            if count != dataset.strings.len() {
                return Err(format!(
                    "dataset {} has shape {:?} but {} strings",
                    dataset.path,
                    dataset.shape,
                    dataset.strings.len()
                ));
            }
        } else if count != dataset.data.len() {
            return Err(format!(
                "dataset {} has shape {:?} but {} values",
                dataset.path,
                dataset.shape,
                dataset.data.len()
            ));
        }
        let bytes = count.saturating_mul(element_size(dataset));
        if bytes > MAX_ARRAY_BYTES {
            return Err(format!(
                "dataset {} exceeds the {MAX_ARRAY_BYTES}-byte limit",
                dataset.path
            ));
        }
        let body = encode_dataset_data(dataset)?;
        lengths.push(u64::try_from(body.len()).unwrap_or(0));
        encoded.push(body);
    }

    // --- the group tree -----------------------------------------------------
    let mut groups: Vec<GroupPlan> = vec![GroupPlan {
        path: "/".to_owned(),
        children: Vec::new(),
        entries: Vec::new(),
        heap: Vec::new(),
        nodes: Vec::new(),
        root_index: 0,
        snod_count: 1,
        header_len: 0,
        header_addr: 0,
        heap_addr: 0,
        btree_addr: 0,
        snod_addr: 0,
    }];
    let mut group_index: std::collections::BTreeMap<String, usize> =
        std::collections::BTreeMap::new();
    group_index.insert("/".to_owned(), 0);
    let mut dataset_paths: Vec<String> = Vec::with_capacity(datasets.len());

    for (position, dataset) in datasets.iter().enumerate() {
        let segments = path_segments(&dataset.path)?;
        let mut parent = 0usize;
        let mut path = String::new();
        for segment in &segments[..segments.len() - 1] {
            path.push('/');
            path.push_str(segment);
            parent = match group_index.get(&path) {
                | Some(&existing) => existing,
                | None => {
                    let index = groups.len();
                    groups.push(GroupPlan {
                        path: path.clone(),
                        children: Vec::new(),
                        entries: Vec::new(),
                        heap: Vec::new(),
                        nodes: Vec::new(),
                        root_index: 0,
                        snod_count: 1,
                        header_len: 0,
                        header_addr: 0,
                        heap_addr: 0,
                        btree_addr: 0,
                        snod_addr: 0,
                    });
                    group_index.insert(path.clone(), index);
                    groups[parent].children.push(PlanChild::Group(index));
                    index
                },
            };
        }
        let full = format!("{path}/{}", segments[segments.len() - 1]);
        if group_index.contains_key(&full) {
            return Err(format!("dataset {} collides with a group", dataset.path));
        }
        if dataset_paths.contains(&full) {
            return Err(format!("duplicate dataset name {full}"));
        }
        dataset_paths.push(full);
        groups[parent].children.push(PlanChild::Dataset(position));
    }

    // --- per-group symbol tables -------------------------------------------
    // Child names are resolved first: building a group's table needs every
    // other group's name, and the groups are mutated in place.
    let child_names: Vec<Vec<String>> = groups
        .iter()
        .map(|group| {
            group
                .children
                .iter()
                .map(|child| {
                    match child {
                        | PlanChild::Group(index) => {
                            groups[*index]
                                .path
                                .rsplit('/')
                                .next()
                                .unwrap_or("")
                                .to_owned()
                        },
                        | PlanChild::Dataset(index) => {
                            dataset_paths[*index]
                                .rsplit('/')
                                .next()
                                .unwrap_or("")
                                .to_owned()
                        },
                    }
                })
                .collect()
        })
        .collect();

    let mut duplicate = None;
    for index in 0..groups.len() {
        let mut order: Vec<usize> = (0..groups[index].children.len()).collect();
        order.sort_by(|left, right| child_names[index][*left].cmp(&child_names[index][*right]));
        let children: Vec<PlanChild> = order
            .iter()
            .map(|&position| groups[index].children[position])
            .collect();
        let names: Vec<String> = order
            .iter()
            .map(|&position| child_names[index][position].clone())
            .collect();
        if names.windows(2).any(|pair| pair[0] == pair[1]) {
            duplicate = Some(groups[index].path.clone());
        }
        let (heap, offsets) = local_heap(&names);
        let entries: Vec<(String, u64, u64)> = names
            .iter()
            .zip(offsets.iter())
            .map(|(name, offset)| (name.clone(), *offset, 0_u64))
            .collect();
        let (nodes, root_index) = build_group_nodes(&entries);
        let snod_count = entries.chunks(SNOD_CAPACITY).count().max(1);
        let group = &mut groups[index];
        group.children = children;
        group.heap = heap;
        group.entries = entries;
        group.nodes = nodes;
        group.root_index = root_index;
        group.snod_count = snod_count;
    }
    if let Some(path) = duplicate {
        return Err(format!("group {path} has two children with the same name"));
    }

    // --- header sizes -------------------------------------------------------
    for group in &mut groups {
        let attrs = attributes_for(attributes, &group.path);
        let mut messages = vec![(0x11_u16, vec![0_u8; 16])];
        for attribute in &attrs {
            messages.push((0x0C, attribute_message(attribute)?));
        }
        group.header_len = object_header_len(&messages);
    }
    let mut header_lengths = Vec::with_capacity(datasets.len());
    for (index, dataset) in datasets.iter().enumerate() {
        let attrs = attributes_for(attributes, &dataset.path);
        let messages = dataset_messages(dataset, &attrs, 0, lengths[index])?;
        header_lengths.push(object_header_len(&messages));
    }

    // --- layout -------------------------------------------------------------
    let mut pos: u64 = 96;
    for group in &mut groups {
        pos = align8(pos);
        group.header_addr = pos;
        pos += u64::try_from(group.header_len).unwrap_or(0);
        group.heap_addr = pos;
        if let Some(slot) = group.heap.get_mut(24..32) {
            slot.copy_from_slice(&(group.heap_addr + 32).to_le_bytes());
        }
        pos += u64::try_from(group.heap.len()).unwrap_or(0);
        group.btree_addr = pos;
        pos += u64::try_from(group.nodes.len().saturating_mul(BTREE_NODE_SIZE)).unwrap_or(0);
        group.snod_addr = pos;
        pos += u64::try_from(group.snod_count.saturating_mul(SNOD_SIZE)).unwrap_or(0);
    }
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

    // --- resolve the child addresses now that every object has a home -------
    for index in 0..groups.len() {
        let addresses: Vec<u64> = groups[index]
            .children
            .iter()
            .map(|child| {
                match child {
                    | PlanChild::Group(target) => groups[*target].header_addr,
                    | PlanChild::Dataset(target) => dataset_addrs[*target].0,
                }
            })
            .collect();
        for (entry, address) in groups[index].entries.iter_mut().zip(addresses) {
            entry.2 = address;
        }
    }

    // --- write --------------------------------------------------------------
    let root = &groups[0];
    let mut out = Vec::with_capacity(usize::try_from(eof).unwrap_or(0));
    out.extend_from_slice(&superblock_v0(
        root.header_addr,
        root.btree_addr,
        root.heap_addr,
        eof,
    ));
    for group in &groups {
        while u64::try_from(out.len()).unwrap_or(u64::MAX) < group.header_addr {
            out.push(0);
        }
        let attrs = attributes_for(attributes, &group.path);
        let mut messages = vec![(
            0x11_u16,
            [
                group.btree_addr.to_le_bytes(),
                group.heap_addr.to_le_bytes(),
            ]
            .concat(),
        )];
        for attribute in &attrs {
            messages.push((0x0C, attribute_message(attribute)?));
        }
        out.extend_from_slice(&object_header_v1(&messages)?);
        while u64::try_from(out.len()).unwrap_or(u64::MAX) < group.heap_addr {
            out.push(0);
        }
        out.extend_from_slice(&group.heap);
        out.extend_from_slice(&group_btree(
            &group.nodes,
            group.root_index,
            group.btree_addr,
            group.snod_addr,
        ));
        for chunk in group.entries.chunks(SNOD_CAPACITY) {
            out.extend_from_slice(&symbol_table_node(chunk));
        }
    }
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
        Dataset::new_numeric(
            "/line",
            vec![6],
            H5Type::F64,
            vec![0.0, 1.5, -2.5, 3.25, 4.0, -5.75],
        )
    }

    /// A two-dimensional `f64` dataset.
    fn grid_dataset() -> Dataset {
        Dataset::new_numeric(
            "/grid",
            vec![2, 3],
            H5Type::F64,
            vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
        )
    }

    /// An `i32` dataset.
    fn count_dataset() -> Dataset {
        Dataset::new_numeric(
            "/counts",
            vec![4],
            H5Type::I32,
            vec![-1.0, 0.0, 7.0, 1_000_000.0],
        )
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
    fn duplicate_and_relative_names_are_rejected() {
        let duplicate = vec![line_dataset(), line_dataset()];
        assert!(write_hdf5_bytes(&duplicate, &[]).is_err());
        let mut relative = line_dataset();
        relative.path = "line".to_owned();
        assert!(write_hdf5_bytes(&[relative], &[]).is_err());
    }

    #[test]
    fn string_datasets_round_trip() {
        let dataset = Dataset::new_strings(
            "/text",
            vec![3],
            vec!["alpha".to_owned(), String::new(), "gamma".to_owned()],
        );
        let bytes = write_hdf5_bytes(&[dataset], &[]).expect("write strings");
        let file = read_hdf5_bytes(&bytes).expect("read strings");
        let text = file
            .datasets
            .iter()
            .find(|dataset| dataset.path == "/text")
            .expect("string dataset");
        assert_eq!(text.h5_type, H5Type::Str);
        assert_eq!(text.strings, vec!["alpha", "", "gamma"]);
        assert_eq!(text.data, [] as [f64; 0]);

        let wrong = Dataset::new_strings("/bad", vec![2], vec!["one".to_owned()]);
        assert!(write_hdf5_bytes(&[wrong], &[]).is_err());
    }

    #[test]
    fn shape_mismatch_is_rejected() {
        let dataset = Dataset::new_numeric("/bad", vec![3], H5Type::F64, vec![1.0, 2.0]);
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
    fn enum_and_variable_length_string_datatypes_parse() {
        // Enumeration over int32: class 8, three members, base type at +8.
        let mut enum_msg = vec![0x18_u8, 0x03, 0, 0, 4, 0, 0, 0];
        enum_msg.extend_from_slice(&[0x10, 0x08, 0, 0, 4, 0, 0, 0, 0, 0, 0x20, 0]);
        let dtype = super::parse_datatype(&enum_msg).expect("enum datatype");
        assert!(matches!(
            dtype,
            super::DType::Enum {
                size: 4,
                signed: true,
                be: false
            }
        ));
        assert_eq!(super::h5_type_of(dtype).expect("enum type"), H5Type::I32);

        // Variable-length string: class 9, kind 1.
        let mut vlen_msg = vec![0x19_u8, 0x01, 0x01, 0, 16, 0, 0, 0];
        vlen_msg.extend_from_slice(&[0x10, 0, 0, 0, 1, 0, 0, 0, 0, 0, 8, 0]);
        let dtype = super::parse_datatype(&vlen_msg).expect("vlen datatype");
        assert!(matches!(dtype, super::DType::VlenStr));
        assert_eq!(super::h5_type_of(dtype).expect("vlen type"), H5Type::Str);
    }

    #[test]
    fn writes_a_group_with_several_symbol_table_nodes() {
        // Twelve datasets need two symbol table nodes and a two-level group
        // B-tree; the file has to read back exactly.
        let mut datasets = vec![
            Dataset::new_numeric(
                "/temperature",
                vec![2, 3],
                H5Type::F64,
                vec![250.0, 251.5, 252.25, 253.0, 254.5, 255.75],
            ),
            Dataset::new_strings(
                "/names",
                vec![3],
                vec!["alpha".to_owned(), "beta".to_owned(), "gamma".to_owned()],
            ),
        ];
        for index in 0..10 {
            datasets.push(Dataset::new_numeric(
                format!("/extra{index}"),
                vec![2],
                H5Type::I64,
                vec![f64::from(index), -1.0],
            ));
        }
        let attributes = vec![(
            "/temperature".to_owned(),
            vec![Attribute {
                name: "units".to_owned(),
                value: AttrValue::Text("K".to_owned()),
            }],
        )];
        let bytes = write_hdf5_bytes(&datasets, &attributes).expect("write");
        let file = read_hdf5_bytes(&bytes).expect("read");
        assert_eq!(file.datasets.len(), datasets.len());
        let temperature = file
            .datasets
            .iter()
            .find(|dataset| dataset.path == "/temperature")
            .expect("temperature");
        assert_eq!(temperature.shape, vec![2, 3]);
        assert_eq!(temperature.data[2], 252.25);
        let names = file
            .datasets
            .iter()
            .find(|dataset| dataset.path == "/names")
            .expect("names");
        assert_eq!(names.strings, vec!["alpha", "beta", "gamma"]);
        let units = file
            .attributes
            .iter()
            .find(|(path, _)| path == "/temperature")
            .map(|(_, attributes)| attributes)
            .expect("attributes");
        assert_eq!(units.len(), 1);
        assert_eq!(units[0].name, "units");
    }

    /// The fixture files under `tests/fixtures/hdf5` are checked in, so this
    /// runs in a normal `cargo test`.
    fn fixture(name: &str) -> std::path::PathBuf {
        std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/hdf5")
            .join(name)
    }

    #[test]
    fn reads_dense_groups() {
        // Small group, one heap block.
        let file = read_hdf5(&fixture("dense_20.h5")).expect("dense 20");
        assert_eq!(file.datasets.len(), 20);

        // Several heap blocks.
        let file = read_hdf5(&fixture("dense_100.h5")).expect("dense 100");
        assert_eq!(file.datasets.len(), 100);

        // Version-2 B-tree with depth two.
        let file = read_hdf5(&fixture("dense_1000.h5")).expect("dense 1000");
        assert_eq!(file.datasets.len(), 1000);

        // A group that tracks creation order: the link info message carries a
        // creation index and the link messages carry a creation order.
        let file = read_hdf5(&fixture("dense_order.h5")).expect("dense order");
        assert_eq!(file.datasets.len(), 300);
    }

    #[test]
    fn reads_dense_group_with_nested_heap_blocks() {
        // 40000 links to 100 objects: the heap needs more rows than the
        // maximum direct block size allows, so the root indirect block points
        // at further indirect blocks.
        let file = read_hdf5(&fixture("dense_nested.h5")).expect("nested heap");
        let datasets: Vec<&str> = file.datasets.iter().map(|d| d.path.as_str()).collect();
        assert_eq!(datasets.len(), 100);
        assert!(datasets.iter().all(|path| path.starts_with("/many/")));
        // Every object is reachable: the links share objects, so the reader
        // keeps the first path it sees for each one.
        assert!(
            datasets
                .iter()
                .any(|path| path.ends_with("/o0") || path.ends_with("/l0"))
        );
        let first = file.datasets.first().expect("a dataset");
        assert_eq!(first.data, vec![0.0, 1.0]);
    }

    #[test]
    fn reads_enum_and_variable_length_string_datasets() {
        let file = read_hdf5(&fixture("types.h5")).expect("typed datasets");
        let strings = file
            .datasets
            .iter()
            .find(|dataset| dataset.path == "/vlen_str")
            .expect("vlen dataset");
        assert_eq!(strings.strings, vec!["alpha", "beta", "gamma"]);
        let enumeration = file
            .datasets
            .iter()
            .find(|dataset| dataset.path == "/enum")
            .expect("enum dataset");
        assert_eq!(enumeration.data, vec![0.0, 2.0, 1.0]);
    }

    #[test]
    fn reads_dense_attributes() {
        // Fourteen attributes force the dense layout: an attribute info
        // message pointing at a fractal heap and a version-2 B-tree.
        let file = read_hdf5(&fixture("dense_attrs.h5")).expect("dense attributes");
        let attributes = file
            .attributes
            .iter()
            .find(|(path, _)| path == "/data")
            .map(|(_, attributes)| attributes)
            .expect("attributes on /data");
        assert_eq!(attributes.len(), 14);
        let mut names: Vec<&str> = attributes.iter().map(|a| a.name.as_str()).collect();
        names.sort_unstable();
        assert_eq!(names.first().copied(), Some("a00"));
        assert_eq!(names.last().copied(), Some("a13"));
        let first = attributes.iter().find(|a| a.name == "a01").expect("a01");
        match &first.value {
            | AttrValue::Float(value) => assert!((value - 1.5).abs() < 1e-12),
            | other => panic!("unexpected value {other:?}"),
        }
    }

    #[test]
    fn writes_nested_groups() {
        let datasets = vec![
            Dataset::new_numeric("/a/b/temperature", vec![2], H5Type::F64, vec![1.0, 2.5]),
            Dataset::new_numeric("/a/counts", vec![3], H5Type::I32, vec![1.0, 2.0, 3.0]),
            Dataset::new_numeric("/top", vec![1], H5Type::I64, vec![7.0]),
            Dataset::new_strings("/a/b/names", vec![2], vec!["x".to_owned(), "y".to_owned()]),
        ];
        let attributes = vec![(
            "/a".to_owned(),
            vec![Attribute {
                name: "kind".to_owned(),
                value: AttrValue::Text("group".to_owned()),
            }],
        )];
        let bytes = write_hdf5_bytes(&datasets, &attributes).expect("write");
        let file = read_hdf5_bytes(&bytes).expect("read");
        let mut paths: Vec<&str> = file.datasets.iter().map(|d| d.path.as_str()).collect();
        paths.sort_unstable();
        assert_eq!(
            paths,
            vec!["/a/b/names", "/a/b/temperature", "/a/counts", "/top"]
        );
        let nested = file
            .datasets
            .iter()
            .find(|d| d.path == "/a/b/temperature")
            .expect("temperature");
        assert_eq!(nested.data, vec![1.0, 2.5]);
        let names = file
            .datasets
            .iter()
            .find(|d| d.path == "/a/b/names")
            .expect("names");
        assert_eq!(names.strings, vec!["x", "y"]);
        assert!(
            file.groups.iter().any(|path| path == "/a/b"),
            "groups: {:?}",
            file.groups
        );
        let group_attributes = file
            .attributes
            .iter()
            .find(|(path, _)| path == "/a")
            .map(|(_, attributes)| attributes)
            .expect("attributes on /a");
        assert_eq!(group_attributes.len(), 1);
    }

    #[test]
    fn reads_compound_members_and_skips_references() {
        let file = read_hdf5(&fixture("compound.h5")).expect("compound fixture");
        // A compound dataset becomes one dataset per representable member.
        let a = file
            .datasets
            .iter()
            .find(|dataset| dataset.path == "/pair/a")
            .expect("member a");
        assert_eq!(a.data, vec![1.0, 2.0, 3.0]);
        let b = file
            .datasets
            .iter()
            .find(|dataset| dataset.path == "/pair/b")
            .expect("member b");
        assert_eq!(b.data, vec![1.5, 2.5, 3.5]);
        // A fixed-width string member becomes text.
        let labels = file
            .datasets
            .iter()
            .find(|dataset| dataset.path == "/mixed/label")
            .expect("member label");
        assert_eq!(labels.strings, vec!["aa", "bb"]);
        let n = file
            .datasets
            .iter()
            .find(|dataset| dataset.path == "/mixed/n")
            .expect("member n");
        assert_eq!(n.data, vec![7.0, 8.0]);
        // A reference dataset cannot be represented and is skipped; the plain
        // dataset beside it still reads.
        assert!(file.datasets.iter().all(|dataset| dataset.path != "/ref"));
        let plain = file
            .datasets
            .iter()
            .find(|dataset| dataset.path == "/plain")
            .expect("numeric dataset");
        assert_eq!(plain.data, vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn shuffle_filter_is_reversed() {
        let shuffled = vec![1_u8, 4, 2, 5, 3, 6];
        let restored = super::filters::unshuffle(&shuffled, 3).expect("unshuffle");
        assert_eq!(restored, vec![1_u8, 2, 3, 4, 5, 6]);
        assert!(super::filters::unshuffle(&shuffled, 0).is_err());
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
            pos += 8 + usize::from(size).div_ceil(8) * 8;
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
            out.extend_from_slice(&u64::from(*row).to_le_bytes());
            out.extend_from_slice(&u64::from(*col).to_le_bytes());
            out.extend_from_slice(&0_u64.to_le_bytes());
            out.extend_from_slice(&chunk_addrs[index].to_le_bytes());
        }
        out.extend_from_slice(&0_u32.to_le_bytes());
        out.extend_from_slice(&0_u32.to_le_bytes());
        out.extend_from_slice(&0_u64.to_le_bytes());
        out.extend_from_slice(&0_u64.to_le_bytes());
        out.extend_from_slice(&0_u64.to_le_bytes());
        out.extend_from_slice(&u64::MAX.to_le_bytes());
        out
    }
}
