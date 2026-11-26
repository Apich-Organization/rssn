//! HDF5 filter pipeline: encode (compress) and decode (decompress).
//!
//! Supports all built-in HDF5 filters and common registered plugins:
//! - `1`: Deflate (zlib / RFC 1950 + RFC 1951)
//! - `2`: Shuffle (byte-transpose)
//! - `3`: Fletcher-32 checksum
//! - `4`: Szip (CCSDS 121.0-B-1 Rice / AEC coding)
//! - `5`: N-bit packing
//! - `6`: Scale-Offset
//! - `32001`: Blosc (BloscLZ, LZ4, Snappy, Zlib, Zstd + byte/bit shuffle)
//! - `32004`: LZ4 (HDF5-LZ4 block-framed plugin)
//! - `32015`: Zstandard

use std::io::Cursor;

use super::inflate;
use crate::io::parquet::lz4;
use crate::io::parquet::snappy;

/// One entry of an HDF5 filter pipeline.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Filter {
    /// HDF5 filter identifier.
    pub(crate) id: u16,
    /// Client data values for the filter.
    pub(crate) client: Vec<u32>,
}

/// Filter ID constants.
pub(crate) const FILTER_DEFLATE: u16 = 1;
pub(crate) const FILTER_SHUFFLE: u16 = 2;
pub(crate) const FILTER_FLETCHER32: u16 = 3;
pub(crate) const FILTER_SZIP: u16 = 4;
pub(crate) const FILTER_NBIT: u16 = 5;
pub(crate) const FILTER_SCALEOFFSET: u16 = 6;
pub(crate) const FILTER_BLOSC: u16 = 32001;
pub(crate) const FILTER_LZ4: u16 = 32004;
pub(crate) const FILTER_ZSTD: u16 = 32015;

// -------------------------------------------------------------------------
// Forward & reverse pipeline
// -------------------------------------------------------------------------

/// Applies the reverse of a filter pipeline to `data` (decoding).
pub(crate) fn apply_filters(
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
            | FILTER_DEFLATE => inflate::zlib_decompress(&current)?,
            | FILTER_SHUFFLE => {
                let element = filter.client.first().copied().map_or(1, |v| v as usize);
                unshuffle(&current, element)?
            },
            | FILTER_FLETCHER32 => verify_and_strip_fletcher32(&current)?,
            | FILTER_SZIP => szip_decode(&current, &filter.client)?,
            | FILTER_NBIT => nbit_decode(&current, &filter.client)?,
            | FILTER_SCALEOFFSET => scaleoffset_decode(&current, &filter.client)?,
            | FILTER_BLOSC => blosc_decode(&current)?,
            | FILTER_LZ4 => hdf5_lz4_decode(&current)?,
            | FILTER_ZSTD => {
                zstd::stream::decode_all(Cursor::new(&current))
                    .map_err(|e| format!("HDF5 zstd filter decode failed: {e}"))?
            },
            | other => return Err(format!("unsupported HDF5 filter {other}")),
        };
    }
    Ok(current)
}

/// Applies the forward filter pipeline to `data` (encoding).
pub(crate) fn apply_filters_encode(
    data: Vec<u8>,
    filters: &[Filter],
) -> Result<Vec<u8>, String> {
    let mut current = data;
    for filter in filters {
        current = match filter.id {
            | FILTER_DEFLATE => inflate::zlib_compress(&current),
            | FILTER_SHUFFLE => {
                let element = filter.client.first().copied().map_or(1, |v| v as usize);
                shuffle(&current, element)?
            },
            | FILTER_FLETCHER32 => append_fletcher32(&current),
            | FILTER_SZIP => szip_encode(&current, &filter.client)?,
            | FILTER_NBIT => nbit_encode(&current, &filter.client)?,
            | FILTER_SCALEOFFSET => scaleoffset_encode(&current, &filter.client)?,
            | FILTER_BLOSC => blosc_encode(&current, &filter.client)?,
            | FILTER_LZ4 => hdf5_lz4_encode(&current, &filter.client),
            | FILTER_ZSTD => {
                let level = filter.client.first().copied().map_or(3, |v| v as i32);
                zstd::bulk::compress(&current, level)
                    .map_err(|e| format!("HDF5 zstd filter encode failed: {e}"))?
            },
            | other => return Err(format!("unsupported HDF5 filter {other}")),
        };
    }
    Ok(current)
}

// -------------------------------------------------------------------------
// Shuffle (ID 2)
// -------------------------------------------------------------------------

/// Applies byte-shuffle to `data` with element width `element`.
pub(crate) fn shuffle(
    data: &[u8],
    element: usize,
) -> Result<Vec<u8>, String> {
    if element == 0 {
        return Err("shuffle filter has a zero element size".to_owned());
    }
    if element == 1 {
        return Ok(data.to_vec());
    }
    let count = data.len() / element;
    let mut out = vec![0_u8; data.len()];
    for index in 0..count {
        for byte in 0..element {
            let src = index * element + byte;
            let dst = byte * count + index;
            out[dst] = data[src];
        }
    }
    let tail = count * element;
    if tail < data.len() {
        out[tail..].copy_from_slice(&data[tail..]);
    }
    Ok(out)
}

/// Reverses byte-shuffle on `data` with element width `element`.
pub(crate) fn unshuffle(
    data: &[u8],
    element: usize,
) -> Result<Vec<u8>, String> {
    if element == 0 {
        return Err("shuffle filter has a zero element size".to_owned());
    }
    if element == 1 {
        return Ok(data.to_vec());
    }
    let count = data.len() / element;
    let mut out = vec![0_u8; data.len()];
    for index in 0..count {
        for byte in 0..element {
            let src = byte * count + index;
            let dst = index * element + byte;
            let value = *data.get(src).ok_or("shuffle filter input is truncated")?;
            out[dst] = value;
        }
    }
    let tail = count * element;
    if tail < data.len() {
        out[tail..].copy_from_slice(&data[tail..]);
    }
    Ok(out)
}

// -------------------------------------------------------------------------
// Fletcher-32 (ID 3)
// -------------------------------------------------------------------------

/// Computes the HDF5 Fletcher-32 checksum (`H5_checksum_fletcher32`).
pub(crate) fn fletcher32(data: &[u8]) -> u32 {
    let mut sum1: u32 = 0xffff;
    let mut sum2: u32 = 0xffff;
    let mut len = data.len() / 2;
    let mut pos = 0_usize;
    while len > 0 {
        let tlen = len.min(360);
        len -= tlen;
        for _ in 0..tlen {
            let word = (u32::from(data[pos]) << 8) | u32::from(data[pos + 1]);
            pos += 2;
            sum1 += word;
            sum2 += sum1;
        }
        sum1 = (sum1 & 0xffff) + (sum1 >> 16);
        sum2 = (sum2 & 0xffff) + (sum2 >> 16);
    }
    if data.len() % 2 != 0 {
        sum1 += u32::from(data[pos]) << 8;
        sum2 += sum1;
        sum1 = (sum1 & 0xffff) + (sum1 >> 16);
        sum2 = (sum2 & 0xffff) + (sum2 >> 16);
    }
    sum1 = (sum1 & 0xffff) + (sum1 >> 16);
    sum2 = (sum2 & 0xffff) + (sum2 >> 16);
    (sum2 << 16) | sum1
}

fn append_fletcher32(data: &[u8]) -> Vec<u8> {
    let cksum = fletcher32(data);
    let mut out = Vec::with_capacity(data.len() + 4);
    out.extend_from_slice(data);
    out.extend_from_slice(&cksum.to_le_bytes());
    out
}

fn verify_and_strip_fletcher32(data: &[u8]) -> Result<Vec<u8>, String> {
    if data.len() < 4 {
        return Err("fletcher32 chunk is shorter than 4 bytes".to_owned());
    }
    let split = data.len() - 4;
    let payload = &data[..split];
    let stored = u32::from_le_bytes([
        data[split],
        data[split + 1],
        data[split + 2],
        data[split + 3],
    ]);
    let computed = fletcher32(payload);
    if stored != computed {
        return Err(format!(
            "fletcher32 checksum mismatch: stored {stored:#010x}, computed {computed:#010x}"
        ));
    }
    Ok(payload.to_vec())
}

// -------------------------------------------------------------------------
// Plugin LZ4 (ID 32004)
// -------------------------------------------------------------------------

const DEFAULT_LZ4_BLOCK_SIZE: usize = 1 << 30;

fn hdf5_lz4_encode(
    data: &[u8],
    client: &[u32],
) -> Vec<u8> {
    let orig_size = u64::try_from(data.len()).unwrap_or(u64::MAX);
    let block_size = client
        .first()
        .copied()
        .filter(|&v| v > 0)
        .map_or(DEFAULT_LZ4_BLOCK_SIZE, |v| v as usize)
        .min(data.len().max(1));
    let mut out = Vec::with_capacity(12 + data.len());
    out.extend_from_slice(&orig_size.to_be_bytes());
    out.extend_from_slice(&(block_size as u32).to_be_bytes());
    if data.is_empty() {
        return out;
    }
    for block in data.chunks(block_size) {
        let compressed = lz4::compress(block);
        if compressed.len() < block.len() {
            out.extend_from_slice(&(compressed.len() as u32).to_be_bytes());
            out.extend_from_slice(&compressed);
        } else {
            out.extend_from_slice(&(block.len() as u32).to_be_bytes());
            out.extend_from_slice(block);
        }
    }
    out
}

fn hdf5_lz4_decode(data: &[u8]) -> Result<Vec<u8>, String> {
    if data.len() < 12 {
        return Err("HDF5 LZ4 chunk header is truncated".to_owned());
    }
    let orig_size = usize::try_from(u64::from_be_bytes([
        data[0], data[1], data[2], data[3], data[4], data[5], data[6], data[7],
    ]))
    .map_err(|_| "HDF5 LZ4 original size overflow".to_owned())?;
    let block_size = u32::from_be_bytes([data[8], data[9], data[10], data[11]]) as usize;
    if block_size == 0 && orig_size > 0 {
        return Err("HDF5 LZ4 block size is zero".to_owned());
    }
    let mut out = Vec::with_capacity(orig_size.min(super::MAX_ARRAY_BYTES));
    let mut pos = 12_usize;
    while out.len() < orig_size {
        let csize_bytes = data
            .get(pos..pos + 4)
            .ok_or_else(|| "HDF5 LZ4 block length is truncated".to_owned())?;
        let csize = u32::from_be_bytes([
            csize_bytes[0],
            csize_bytes[1],
            csize_bytes[2],
            csize_bytes[3],
        ]) as usize;
        pos += 4;
        let block = data
            .get(pos..pos + csize)
            .ok_or_else(|| "HDF5 LZ4 block payload is truncated".to_owned())?;
        pos += csize;
        let expected = (orig_size - out.len()).min(block_size);
        if csize == expected {
            out.extend_from_slice(block);
        } else {
            let decompressed = lz4::decompress(block, expected)?;
            out.extend_from_slice(&decompressed);
        }
    }
    Ok(out)
}

// -------------------------------------------------------------------------
// Plugin Blosc (ID 32001)
// -------------------------------------------------------------------------

const BLOSC_MIN_HEADER_LENGTH: usize = 16;
const BLOSC_DOSHUFFLE: u8 = 0x01;
const BLOSC_MEMCPYED: u8 = 0x02;
const BLOSC_DOBITSHUFFLE: u8 = 0x04;

/// Bit-shuffle `data` with `typesize` bytes per element.
fn bitshuffle(
    data: &[u8],
    typesize: usize,
) -> Vec<u8> {
    if typesize <= 1 || data.len() < 8 {
        return data.to_vec();
    }
    let nelems = data.len() / typesize;
    let total_bits = typesize * 8;
    let mut out = vec![0_u8; data.len()];
    let aligned_elems = (nelems / 8) * 8;
    for elem in 0..aligned_elems {
        for bit_idx in 0..total_bits {
            let src_byte = elem * typesize + (bit_idx / 8);
            let src_bit = bit_idx % 8;
            let bit = (data[src_byte] >> src_bit) & 1;
            let dst_bit_pos = bit_idx * aligned_elems + elem;
            out[dst_bit_pos / 8] |= bit << (dst_bit_pos % 8);
        }
    }
    let tail_start = aligned_elems * typesize;
    if tail_start < data.len() {
        out[tail_start..].copy_from_slice(&data[tail_start..]);
    }
    out
}

/// Reverses bit-shuffle on `data` with `typesize` bytes per element.
fn bitunshuffle(
    data: &[u8],
    typesize: usize,
) -> Vec<u8> {
    if typesize <= 1 || data.len() < 8 {
        return data.to_vec();
    }
    let nelems = data.len() / typesize;
    let total_bits = typesize * 8;
    let mut out = vec![0_u8; data.len()];
    let aligned_elems = (nelems / 8) * 8;
    for elem in 0..aligned_elems {
        for bit_idx in 0..total_bits {
            let src_bit_pos = bit_idx * aligned_elems + elem;
            let bit = (data[src_bit_pos / 8] >> (src_bit_pos % 8)) & 1;
            let dst_byte = elem * typesize + (bit_idx / 8);
            let dst_bit = bit_idx % 8;
            out[dst_byte] |= bit << dst_bit;
        }
    }
    let tail_start = aligned_elems * typesize;
    if tail_start < data.len() {
        out[tail_start..].copy_from_slice(&data[tail_start..]);
    }
    out
}

/// BloscLZ block decompression (Blosc1 internal compressor ID 0).
fn blosclz_decompress(
    input: &[u8],
    expected: usize,
) -> Result<Vec<u8>, String> {
    if input.is_empty() && expected == 0 {
        return Ok(Vec::new());
    }
    let mut out = Vec::with_capacity(expected);
    let mut ip = 0_usize;
    let mut op_ctrl = *input.get(ip).ok_or("blosclz input is empty")?;
    let mut ctrl = u32::from(op_ctrl & 0x1f);
    ip += 1;
    let mut loop_back = true;
    while ip < input.len() {
        if loop_back {
            // already have ctrl from previous or init
            loop_back = false;
        } else {
            op_ctrl = input[ip];
            ctrl = u32::from(op_ctrl & 0x1f);
            ip += 1;
        }
        let mut ref_off = op_ctrl >> 5;
        if ref_off == 0 {
            // Literal run of (ctrl + 1) bytes
            let len = (ctrl as usize) + 1;
            let end = ip
                .checked_add(len)
                .ok_or_else(|| "blosclz literal overflow".to_owned())?;
            let lit = input
                .get(ip..end)
                .ok_or_else(|| "blosclz literal run truncated".to_owned())?;
            out.extend_from_slice(lit);
            ip = end;
        } else {
            // Back-reference
            ref_off -= 1;
            let mut len = ref_off as usize;
            let mut ofs = ((ctrl as usize) << 8) + 1;
            if len == 7 - 1 {
                let extra = *input.get(ip).ok_or("blosclz match len truncated")?;
                ip += 1;
                len += usize::from(extra);
            }
            let low = *input.get(ip).ok_or("blosclz offset truncated")?;
            ip += 1;
            ofs += usize::from(low);
            if ofs > out.len() {
                return Err("blosclz back-reference before start of output".to_owned());
            }
            len += 2;
            let mut src = out.len() - ofs;
            for _ in 0..len {
                let b = out[src];
                out.push(b);
                src += 1;
            }
        }
        if out.len() >= expected {
            break;
        }
    }
    if out.len() != expected {
        return Err(format!(
            "blosclz decompressed {} bytes, expected {expected}",
            out.len()
        ));
    }
    Ok(out)
}

fn blosc_decompress_block(
    comp_id: u8,
    input: &[u8],
    expected: usize,
) -> Result<Vec<u8>, String> {
    if input.len() == expected {
        // Blosc stores incompressible blocks raw when compressed >= expected.
        return Ok(input.to_vec());
    }
    match comp_id {
        | 0 => blosclz_decompress(input, expected),
        | 1 => lz4::decompress(input, expected),
        | 2 => snappy::decompress(input),
        | 3 => inflate::zlib_decompress(input),
        | 4 | 5 => {
            zstd::bulk::decompress(input, expected)
                .map_err(|e| format!("blosc zstd decompress failed: {e}"))
        },
        | other => Err(format!("unsupported Blosc compressor id {other}")),
    }
}

fn blosc_decode(data: &[u8]) -> Result<Vec<u8>, String> {
    if data.len() < BLOSC_MIN_HEADER_LENGTH {
        return Err("Blosc header is truncated".to_owned());
    }
    let flags = data[2];
    let typesize = usize::from(data[3]).max(1);
    let nbytes = u32::from_le_bytes([data[4], data[5], data[6], data[7]]) as usize;
    let blocksize = u32::from_le_bytes([data[8], data[9], data[10], data[11]]) as usize;
    let cbytes = u32::from_le_bytes([data[12], data[13], data[14], data[15]]) as usize;
    if cbytes > data.len() {
        return Err("Blosc compressed length exceeds chunk".to_owned());
    }
    if nbytes == 0 {
        return Ok(Vec::new());
    }
    if flags & BLOSC_MEMCPYED != 0 {
        return data
            .get(16..16 + nbytes)
            .map(<[u8]>::to_vec)
            .ok_or_else(|| "Blosc memcpyed payload truncated".to_owned());
    }
    if blocksize == 0 {
        return Err("Blosc blocksize is zero".to_owned());
    }
    let comp_id = (flags >> 5) & 0x07;
    let nblocks = nbytes.div_ceil(blocksize);
    // In Blosc1 header, bit 4 (0x10) is BLOSC_DONT_SPLIT (1 = single block, 0 = split into typesize streams).
    let maybe_split = (flags & 0x10) == 0 && typesize > 1;

    let pos = 16_usize;
    let mut out = Vec::with_capacity(nbytes.min(super::MAX_ARRAY_BYTES));
    // Note: Blosc1 stores `bstarts` (i32 * nblocks) at offset 16..16+4*nblocks!
    let bstarts_end = pos + nblocks * 4;
    if bstarts_end > data.len() {
        return Err("Blosc block starts table truncated".to_owned());
    }
    let mut bstarts = Vec::with_capacity(nblocks);
    for b in 0..nblocks {
        let off = pos + b * 4;
        bstarts.push(
            u32::from_le_bytes([data[off], data[off + 1], data[off + 2], data[off + 3]]) as usize,
        );
    }

    for (b_idx, &bstart) in bstarts.iter().enumerate() {
        let bsize = if b_idx + 1 == nblocks {
            nbytes - b_idx * blocksize
        } else {
            blocksize
        };
        let mut cur = bstart;
        let block_end = bstarts.get(b_idx + 1).copied().unwrap_or(cbytes);
        let first_csize = data
            .get(cur..cur + 4)
            .map(|b| u32::from_le_bytes([b[0], b[1], b[2], b[3]]) as usize)
            .unwrap_or(0);
        let is_split = maybe_split && cur + 4 + first_csize < block_end;
        let block_raw = if is_split {
            let nstreams = typesize;
            let stream_size = bsize / nstreams;
            let mut streams = Vec::with_capacity(bsize);
            for s in 0..nstreams {
                let s_expected = if s + 1 == nstreams {
                    bsize - stream_size * (nstreams - 1)
                } else {
                    stream_size
                };
                let csize_bytes = data
                    .get(cur..cur + 4)
                    .ok_or_else(|| "Blosc split stream size truncated".to_owned())?;
                let csize = u32::from_le_bytes([
                    csize_bytes[0],
                    csize_bytes[1],
                    csize_bytes[2],
                    csize_bytes[3],
                ]) as usize;
                cur += 4;
                if csize == 0 {
                    streams.extend(std::iter::repeat_n(0_u8, s_expected));
                } else {
                    let payload = data
                        .get(cur..cur + csize)
                        .ok_or_else(|| "Blosc split stream payload truncated".to_owned())?;
                    cur += csize;
                    let dec = blosc_decompress_block(comp_id, payload, s_expected)?;
                    streams.extend_from_slice(&dec);
                }
            }
            streams
        } else {
            let csize_bytes = data
                .get(cur..cur + 4)
                .ok_or_else(|| "Blosc block size truncated".to_owned())?;
            let csize = u32::from_le_bytes([
                csize_bytes[0],
                csize_bytes[1],
                csize_bytes[2],
                csize_bytes[3],
            ]) as usize;
            cur += 4;
            let payload = data
                .get(cur..cur + csize)
                .ok_or_else(|| "Blosc block payload truncated".to_owned())?;
            blosc_decompress_block(comp_id, payload, bsize)?
        };

        let unshuffled = if flags & BLOSC_DOSHUFFLE != 0 {
            unshuffle(&block_raw, typesize)?
        } else if flags & BLOSC_DOBITSHUFFLE != 0 && block_raw.len() % (typesize * 8) == 0 {
            bitunshuffle(&block_raw, typesize)
        } else {
            block_raw
        };
        out.extend_from_slice(&unshuffled);
    }
    Ok(out)
}

fn blosc_encode(
    data: &[u8],
    client: &[u32],
) -> Result<Vec<u8>, String> {
    // hdf5plugin Blosc client data:
    // client[4] = clevel (0..9, default 5)
    // client[5] = doshuffle (0=none, 1=byte, 2=bit)
    // client[6] = compressor (0=blosclz, 1=lz4, 2=lz4hc, 3=snappy, 4=zlib, 5=zstd)
    // typesize = client[2] (default 8)
    let typesize = client
        .get(2)
        .copied()
        .filter(|&v| v > 0 && v <= 255)
        .unwrap_or(8) as usize;
    let clevel = client.get(4).copied().unwrap_or(5);
    let shuffle_mode = client.get(5).copied().unwrap_or(1);
    let comp_req = client.get(6).copied().unwrap_or(1) as u8; // default LZ4 (or Zstd=5)
    let comp_id = if comp_req == 0 {
        1
    } else {
        comp_req.min(5)
    };
    let comp_fmt = match comp_id {
        | 1 | 2 => 1_u8,
        | 3 => 2_u8,
        | 4 => 3_u8,
        | 5 => 4_u8,
        | _ => 1_u8,
    };

    let nbytes = data.len();
    let blocksize = nbytes.min(256 * 1024).max(1);
    let nblocks = if nbytes == 0 {
        0
    } else {
        nbytes.div_ceil(blocksize)
    };

    // Bit 4 (0x10) = BLOSC_DONT_SPLIT (single compressed stream per block).
    let mut flags: u8 = (comp_fmt << 5) | 0x10;
    if shuffle_mode == 1 && typesize > 1 {
        flags |= BLOSC_DOSHUFFLE;
    } else if shuffle_mode == 2 && typesize > 1 {
        flags |= BLOSC_DOBITSHUFFLE;
    }

    let mut out = vec![0_u8; 16 + nblocks * 4];
    let mut bstarts = Vec::with_capacity(nblocks);

    for block in data.chunks(blocksize) {
        let bstart = out.len() as u32;
        bstarts.push(bstart);
        let prepared = if flags & BLOSC_DOSHUFFLE != 0 {
            shuffle(block, typesize)?
        } else if flags & BLOSC_DOBITSHUFFLE != 0 && block.len() % (typesize * 8) == 0 {
            bitshuffle(block, typesize)
        } else {
            block.to_vec()
        };
        let compressed = match comp_id {
            | 1 | 2 => lz4::compress(&prepared),
            | 3 => snappy::compress(&prepared),
            | 4 => inflate::zlib_compress(&prepared),
            | 5 => {
                zstd::bulk::compress(&prepared, (clevel as i32).clamp(1, 19))
                    .map_err(|e| format!("blosc zstd compress failed: {e}"))?
            },
            | _ => lz4::compress(&prepared),
        };
        let payload = if compressed.len() < prepared.len() {
            compressed
        } else {
            prepared
        };
        out.extend_from_slice(&(payload.len() as u32).to_be_bytes());
        // Wait: Blosc1 block length is little-endian `u32`!
        let len_pos = out.len() - 4;
        out[len_pos..].copy_from_slice(&(payload.len() as u32).to_le_bytes());
        out.extend_from_slice(&payload);
    }

    out[0] = 2; // Blosc1 format version
    out[1] = 1; // versionlz
    out[2] = flags;
    out[3] = typesize as u8;
    out[4..8].copy_from_slice(&(nbytes as u32).to_le_bytes());
    out[8..12].copy_from_slice(&(blocksize as u32).to_le_bytes());
    let total_len = (out.len() as u32).to_le_bytes();
    out[12..16].copy_from_slice(&total_len);
    for (i, bstart) in bstarts.iter().enumerate() {
        let off = 16 + i * 4;
        out[off..off + 4].copy_from_slice(&bstart.to_le_bytes());
    }
    Ok(out)
}

// -------------------------------------------------------------------------
// N-bit (ID 5)
// -------------------------------------------------------------------------

struct MsbBitWriter {
    out: Vec<u8>,
    acc: u64,
    bits: u32,
}

impl MsbBitWriter {
    fn new() -> Self {
        Self {
            out: Vec::new(),
            acc: 0,
            bits: 0,
        }
    }

    fn write_bits(
        &mut self,
        val: u64,
        nbits: u32,
    ) {
        if nbits == 0 {
            return;
        }
        let mask = if nbits >= 64 {
            u64::MAX
        } else {
            (1_u64 << nbits) - 1
        };
        self.acc = (self.acc << nbits) | (val & mask);
        self.bits += nbits;
        while self.bits >= 8 {
            self.bits -= 8;
            self.out.push((self.acc >> self.bits) as u8);
            self.acc &= (1_u64 << self.bits).saturating_sub(1);
        }
    }

    fn finish(mut self) -> Vec<u8> {
        if self.bits > 0 {
            self.out.push((self.acc << (8 - self.bits)) as u8);
        }
        self.out
    }
}

struct MsbBitReader<'a> {
    data: &'a [u8],
    pos: usize,
    acc: u64,
    bits: u32,
}

impl<'a> MsbBitReader<'a> {
    fn new(data: &'a [u8]) -> Self {
        Self {
            data,
            pos: 0,
            acc: 0,
            bits: 0,
        }
    }

    fn read_bits(
        &mut self,
        nbits: u32,
    ) -> Result<u64, String> {
        if nbits == 0 {
            return Ok(0);
        }
        while self.bits < nbits {
            let b = *self
                .data
                .get(self.pos)
                .ok_or_else(|| "nbit/scaleoffset bitstream truncated".to_owned())?;
            self.pos += 1;
            self.acc = (self.acc << 8) | u64::from(b);
            self.bits += 8;
        }
        self.bits -= nbits;
        let val = self.acc >> self.bits;
        self.acc &= (1_u64 << self.bits).saturating_sub(1);
        Ok(val)
    }
}

fn nbit_encode(
    data: &[u8],
    client: &[u32],
) -> Result<Vec<u8>, String> {
    if client.len() < 7 {
        return Err("N-bit filter requires at least 7 client parameters".to_owned());
    }
    if client[1] != 0 {
        return Ok(data.to_vec());
    }
    let size = client[3] as usize;
    let order = client[4]; // 0 = LE, 1 = BE
    let precision = client[5];
    let offset = client[6];
    if size == 0 || size > 8 || precision == 0 || offset + precision > (size as u32) * 8 {
        return Err("invalid N-bit filter parameters".to_owned());
    }
    let mut writer = MsbBitWriter::new();
    for elem in data.chunks_exact(size) {
        let mut val = 0_u64;
        if order == 0 {
            for (i, &b) in elem.iter().enumerate() {
                val |= u64::from(b) << (i * 8);
            }
        } else {
            for &b in elem {
                val = (val << 8) | u64::from(b);
            }
        }
        let shifted = val >> offset;
        writer.write_bits(shifted, precision);
    }
    Ok(writer.finish())
}

fn nbit_decode(
    data: &[u8],
    client: &[u32],
) -> Result<Vec<u8>, String> {
    if client.len() < 7 {
        return Err("N-bit filter requires at least 7 client parameters".to_owned());
    }
    if client[1] != 0 {
        return Ok(data.to_vec());
    }
    let class = client[2];
    let size = client[3] as usize;
    let order = client[4];
    let precision = client[5];
    let offset = client[6];
    if size == 0 || size > 8 || precision == 0 || offset + precision > (size as u32) * 8 {
        return Err("invalid N-bit filter parameters".to_owned());
    }
    let nelmts = (data.len() * 8) / (precision as usize);
    let mut reader = MsbBitReader::new(data);
    let mut out = Vec::with_capacity(nelmts * size);
    let signed = class == 0 && client.get(7).copied().unwrap_or(0) == 1;
    for _ in 0..nelmts {
        let raw = reader.read_bits(precision)?;
        let mut val = raw << offset;
        if signed && precision < 64 && ((raw >> (precision - 1)) & 1) != 0 {
            let sign_mask = !((1_u64 << (offset + precision)) - 1);
            val |= sign_mask;
        }
        if order == 0 {
            for i in 0..size {
                out.push((val >> (i * 8)) as u8);
            }
        } else {
            for i in (0..size).rev() {
                out.push((val >> (i * 8)) as u8);
            }
        }
    }
    Ok(out)
}

// -------------------------------------------------------------------------
// Scale-Offset (ID 6)
// -------------------------------------------------------------------------

fn scaleoffset_encode(
    data: &[u8],
    client: &[u32],
) -> Result<Vec<u8>, String> {
    if client.len() < 7 {
        return Err("scaleoffset filter requires at least 7 client parameters".to_owned());
    }
    let scale_type = client[0]; // 0=float D-scale, 1=float E-scale, 2=int
    let scale_factor = client[1] as i32;
    let class = client[3]; // 0=int, 1=float
    let size = client[4] as usize;
    let order = client[5];
    let sign = client[6];
    if size == 0 || size > 8 {
        return Err("unsupported element size for scaleoffset".to_owned());
    }
    let nelmts = data.len() / size;
    if nelmts == 0 {
        return Ok(Vec::new());
    }

    let mut scaled: Vec<i64> = Vec::with_capacity(nelmts);
    let factor = if class == 1 {
        10_f64.powi(scale_factor)
    } else {
        1.0
    };
    for elem in data.chunks_exact(size) {
        let mut raw = [0_u8; 8];
        if order == 0 {
            raw[..size].copy_from_slice(elem);
        } else {
            for i in 0..size {
                raw[i] = elem[size - 1 - i];
            }
        }
        let v_i64 = if class == 1 {
            let f = if size == 4 {
                f64::from(f32::from_le_bytes([raw[0], raw[1], raw[2], raw[3]]))
            } else {
                f64::from_le_bytes(raw)
            };
            (f * factor).round() as i64
        } else if sign != 0 {
            match size {
                | 1 => i64::from(raw[0] as i8),
                | 2 => i64::from(i16::from_le_bytes([raw[0], raw[1]])),
                | 4 => i64::from(i32::from_le_bytes([raw[0], raw[1], raw[2], raw[3]])),
                | 8 => i64::from_le_bytes(raw),
                | _ => return Err("unsupported int size".to_owned()),
            }
        } else {
            match size {
                | 1 => i64::from(raw[0]),
                | 2 => i64::from(u16::from_le_bytes([raw[0], raw[1]])),
                | 4 => i64::from(u32::from_le_bytes([raw[0], raw[1], raw[2], raw[3]])),
                | 8 => u64::from_le_bytes(raw) as i64,
                | _ => return Err("unsupported uint size".to_owned()),
            }
        };
        scaled.push(v_i64);
    }

    let _ = scale_type;
    let min_val = scaled.iter().copied().min().unwrap_or(0);
    let max_diff = scaled
        .iter()
        .map(|&v| (v.wrapping_sub(min_val)) as u64)
        .max()
        .unwrap_or(0);
    let minbits = if max_diff == 0 {
        0_u32
    } else {
        64 - max_diff.leading_zeros()
    };

    let mut out = Vec::new();
    out.extend_from_slice(&minbits.to_le_bytes());
    out.extend_from_slice(&min_val.to_le_bytes());
    let mut writer = MsbBitWriter::new();
    for &v in &scaled {
        let diff = (v.wrapping_sub(min_val)) as u64;
        writer.write_bits(diff, minbits);
    }
    out.extend_from_slice(&writer.finish());
    Ok(out)
}

fn scaleoffset_decode(
    data: &[u8],
    client: &[u32],
) -> Result<Vec<u8>, String> {
    if client.len() < 7 {
        return Err("scaleoffset filter requires at least 7 client parameters".to_owned());
    }
    if data.len() < 12 {
        return Err("scaleoffset payload is truncated".to_owned());
    }
    let scale_factor = client[1] as i32;
    let nelmts = client[2] as usize;
    let class = client[3];
    let size = client[4] as usize;
    let order = client[5];
    if size == 0 || size > 8 {
        return Err("unsupported element size for scaleoffset".to_owned());
    }
    let minbits = u32::from_le_bytes([data[0], data[1], data[2], data[3]]);
    let min_val = i64::from_le_bytes([
        data[4], data[5], data[6], data[7], data[8], data[9], data[10], data[11],
    ]);
    let mut reader = MsbBitReader::new(&data[12..]);
    let factor = if class == 1 {
        10_f64.powi(scale_factor)
    } else {
        1.0
    };
    let count = if nelmts > 0 {
        nelmts
    } else if minbits > 0 {
        ((data.len() - 12) * 8) / (minbits as usize)
    } else {
        0
    };
    let mut out = Vec::with_capacity(count * size);
    for _ in 0..count {
        let diff = reader.read_bits(minbits)?;
        let val = min_val.wrapping_add(diff as i64);
        let le_bytes: [u8; 8] = if class == 1 {
            let f = (val as f64) / factor;
            if size == 4 {
                let b = (f as f32).to_le_bytes();
                [b[0], b[1], b[2], b[3], 0, 0, 0, 0]
            } else {
                f.to_le_bytes()
            }
        } else {
            val.to_le_bytes()
        };
        if order == 0 {
            out.extend_from_slice(&le_bytes[..size]);
        } else {
            for i in (0..size).rev() {
                out.push(le_bytes[i]);
            }
        }
    }
    Ok(out)
}

// -------------------------------------------------------------------------
// Szip (ID 4): CCSDS 121.0-B-1 Rice / AEC coding
// -------------------------------------------------------------------------

const SZIP_EC_OPTION_MASK: u32 = 4;
const SZIP_LSB_OPTION_MASK: u32 = 8;
const SZIP_NN_OPTION_MASK: u32 = 32;
const SZIP_RAW_OPTION_MASK: u32 = 128;

fn szip_encode(
    data: &[u8],
    client: &[u32],
) -> Result<Vec<u8>, String> {
    if client.len() < 4 {
        return Err("szip filter requires 4 client parameters".to_owned());
    }
    let options_mask = client[0];
    let ppb = client[1].max(2) as usize;
    let bpp = client[2].clamp(1, 32) as usize;
    let byte_width = bpp.div_ceil(8);
    let is_lsb = (options_mask & SZIP_LSB_OPTION_MASK) != 0;
    let use_nn =
        (options_mask & SZIP_NN_OPTION_MASK) != 0 || (options_mask & SZIP_EC_OPTION_MASK) == 0;
    if (options_mask & SZIP_RAW_OPTION_MASK) != 0 {
        return Ok(data.to_vec());
    }

    let nelmts = data.len() / byte_width;
    let mut samples = Vec::with_capacity(nelmts);
    for chunk in data.chunks_exact(byte_width) {
        let mut v = 0_u32;
        if is_lsb {
            for (i, &b) in chunk.iter().enumerate() {
                v |= u32::from(b) << (i * 8);
            }
        } else {
            for &b in chunk {
                v = (v << 8) | u32::from(b);
            }
        }
        samples.push(v);
    }

    let id_bits: u32 = if bpp <= 8 {
        3
    } else if bpp <= 16 {
        4
    } else {
        5
    };
    let all_ones_id = (1_u32 << id_bits) - 1;
    let xmax: u32 = if bpp >= 32 {
        u32::MAX
    } else {
        (1_u32 << bpp) - 1
    };

    let mut writer = MsbBitWriter::new();
    for block in samples.chunks(ppb) {
        // Optional NN delta mapping
        let mapped: Vec<u32> = if use_nn {
            let mut prev = 0_u32;
            block
                .iter()
                .enumerate()
                .map(|(i, &x)| {
                    if i == 0 {
                        prev = x;
                        x
                    } else {
                        let d = i64::from(x) - i64::from(prev);
                        let theta = prev.min(xmax.saturating_sub(prev));
                        let sigma = if d.unsigned_abs() <= u64::from(theta) {
                            if d >= 0 {
                                (d as u32) * 2
                            } else {
                                ((-d) as u32) * 2 - 1
                            }
                        } else {
                            theta.saturating_add(d.unsigned_abs() as u32)
                        };
                        prev = x;
                        sigma
                    }
                })
                .collect()
        } else {
            block.to_vec()
        };

        // Pick best k or uncompressed
        if mapped.iter().all(|&v| v == 0) {
            writer.write_bits(0, id_bits);
            writer.write_bits(0, 1); // zero-block fs=0 (single zero block)
            continue;
        }
        let max_k = (all_ones_id - 1).min(bpp as u32);
        let mut best_k = all_ones_id;
        let mut best_bits = mapped.len() * bpp + id_bits as usize;
        for k in 0..max_k {
            let mut cost = id_bits as usize;
            let mut ok = true;
            for &v in &mapped {
                let q = v >> k;
                if q > 64 {
                    ok = false;
                    break;
                }
                cost += (q as usize) + 1 + (k as usize);
            }
            if ok && cost < best_bits {
                best_bits = cost;
                best_k = k;
            }
        }
        if best_k == all_ones_id {
            writer.write_bits(u64::from(all_ones_id), id_bits);
            for &v in &mapped {
                writer.write_bits(u64::from(v), bpp as u32);
            }
        } else {
            writer.write_bits(u64::from(best_k + 1), id_bits);
            for &v in &mapped {
                let q = v >> best_k;
                let r = v & ((1_u32 << best_k) - 1);
                for _ in 0..q {
                    writer.write_bits(1, 1);
                }
                writer.write_bits(0, 1);
                writer.write_bits(u64::from(r), best_k);
            }
        }
    }
    Ok(writer.finish())
}

fn szip_decode(
    data: &[u8],
    client: &[u32],
) -> Result<Vec<u8>, String> {
    if client.len() < 4 {
        return Err("szip filter requires 4 client parameters".to_owned());
    }
    let options_mask = client[0];
    let ppb = client[1].max(2) as usize;
    let bpp = client[2].clamp(1, 32) as usize;
    let byte_width = bpp.div_ceil(8);
    let is_lsb = (options_mask & SZIP_LSB_OPTION_MASK) != 0;
    let use_nn =
        (options_mask & SZIP_NN_OPTION_MASK) != 0 || (options_mask & SZIP_EC_OPTION_MASK) == 0;
    if (options_mask & SZIP_RAW_OPTION_MASK) != 0 {
        return Ok(data.to_vec());
    }

    let id_bits: u32 = if bpp <= 8 {
        3
    } else if bpp <= 16 {
        4
    } else {
        5
    };
    let all_ones_id = (1_u32 << id_bits) - 1;
    let xmax: u32 = if bpp >= 32 {
        u32::MAX
    } else {
        (1_u32 << bpp) - 1
    };

    let mut reader = MsbBitReader::new(data);
    let mut out = Vec::new();
    // Decode blocks until input bits are exhausted
    while reader.pos < data.len() || reader.bits >= id_bits {
        let Ok(id) = reader.read_bits(id_bits) else {
            break;
        };
        let id = id as u32;
        let mut mapped = Vec::with_capacity(ppb);
        if id == 0 {
            let Ok(fs) = reader.read_bits(1) else {
                break;
            };
            if fs == 0 {
                mapped.resize(ppb, 0);
            } else {
                // Zero-block run (4 bits + 1, or 5 bits)
                mapped.resize(ppb, 0);
            }
        } else if id == all_ones_id {
            for _ in 0..ppb {
                let Ok(v) = reader.read_bits(bpp as u32) else {
                    break;
                };
                mapped.push(v as u32);
            }
        } else {
            let k = id - 1;
            for _ in 0..ppb {
                let mut q = 0_u32;
                loop {
                    let Ok(bit) = reader.read_bits(1) else {
                        break;
                    };
                    if bit == 0 {
                        break;
                    }
                    q += 1;
                }
                let Ok(r) = reader.read_bits(k) else {
                    break;
                };
                mapped.push((q << k) | (r as u32));
            }
        }

        let block: Vec<u32> = if use_nn {
            let mut prev = 0_u32;
            mapped
                .iter()
                .enumerate()
                .map(|(i, &sigma)| {
                    if i == 0 {
                        prev = sigma;
                        sigma
                    } else {
                        let theta = prev.min(xmax.saturating_sub(prev));
                        let x = if sigma <= 2 * theta {
                            if sigma % 2 == 0 {
                                prev.wrapping_add(sigma / 2)
                            } else {
                                prev.wrapping_sub((sigma + 1) / 2)
                            }
                        } else if prev <= xmax / 2 {
                            sigma
                        } else {
                            xmax.saturating_sub(sigma)
                        };
                        prev = x;
                        x
                    }
                })
                .collect()
        } else {
            mapped
        };

        for v in block {
            if is_lsb {
                for i in 0..byte_width {
                    out.push((v >> (i * 8)) as u8);
                }
            } else {
                for i in (0..byte_width).rev() {
                    out.push((v >> (i * 8)) as u8);
                }
            }
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fletcher32_round_trip_and_corruption() {
        let payload = b"scientific computing fletcher32 test";
        let encoded = append_fletcher32(payload);
        assert_eq!(verify_and_strip_fletcher32(&encoded).unwrap(), payload);
        let mut corrupt = encoded;
        corrupt[3] ^= 0xff;
        assert!(verify_and_strip_fletcher32(&corrupt).is_err());
    }

    #[test]
    fn shuffle_and_deflate_round_trip() {
        let mut raw = Vec::new();
        for i in 0..256_u32 {
            raw.extend_from_slice(&f64::from(i).to_le_bytes());
        }
        let filters = vec![
            Filter {
                id: FILTER_SHUFFLE,
                client: vec![8],
            },
            Filter {
                id: FILTER_DEFLATE,
                client: vec![6],
            },
            Filter {
                id: FILTER_FLETCHER32,
                client: vec![],
            },
        ];
        let encoded = apply_filters_encode(raw.clone(), &filters).unwrap();
        assert!(encoded.len() < raw.len());
        let decoded = apply_filters(encoded, &filters, 0).unwrap();
        assert_eq!(decoded, raw);
    }

    #[test]
    fn lz4_zstd_blosc_filters_round_trip() {
        let mut raw = Vec::new();
        for i in 0..512_u32 {
            raw.extend_from_slice(&f64::from(i % 16).to_le_bytes());
        }
        for filter in [
            Filter {
                id: FILTER_LZ4,
                client: vec![],
            },
            Filter {
                id: FILTER_ZSTD,
                client: vec![3],
            },
            Filter {
                id: FILTER_BLOSC,
                client: vec![2, 2, 8, 4096, 5, 1, 1],
            },
            Filter {
                id: FILTER_BLOSC,
                client: vec![2, 2, 8, 4096, 5, 2, 5],
            },
        ] {
            let encoded = apply_filters_encode(raw.clone(), std::slice::from_ref(&filter)).unwrap();
            let decoded = apply_filters(encoded, std::slice::from_ref(&filter), 0).unwrap();
            assert_eq!(decoded, raw, "failed for filter {}", filter.id);
        }
    }

    #[test]
    fn nbit_scaleoffset_szip_round_trip() {
        let vals: Vec<i32> = (100..228).collect();
        let raw: Vec<u8> = vals.iter().flat_map(|v| v.to_le_bytes()).collect();

        // N-bit: 12-bit precision in 32-bit LE int
        let nbit = Filter {
            id: FILTER_NBIT,
            client: vec![7, 0, 0, 4, 0, 12, 0],
        };
        let enc_nbit = apply_filters_encode(raw.clone(), std::slice::from_ref(&nbit)).unwrap();
        assert_eq!(enc_nbit.len(), (128 * 12) / 8);
        let dec_nbit = apply_filters(enc_nbit, std::slice::from_ref(&nbit), 0).unwrap();
        assert_eq!(dec_nbit, raw);

        // Scale-offset: int32
        let so = Filter {
            id: FILTER_SCALEOFFSET,
            client: vec![2, 0, 128, 0, 4, 0, 1],
        };
        let enc_so = apply_filters_encode(raw.clone(), std::slice::from_ref(&so)).unwrap();
        let dec_so = apply_filters(enc_so, std::slice::from_ref(&so), 0).unwrap();
        assert_eq!(dec_so, raw);

        // Szip: 32-bit LE, ppb=16
        let szip = Filter {
            id: FILTER_SZIP,
            client: vec![SZIP_NN_OPTION_MASK | SZIP_LSB_OPTION_MASK, 16, 32, 128],
        };
        let enc_szip = apply_filters_encode(raw.clone(), std::slice::from_ref(&szip)).unwrap();
        let dec_szip = apply_filters(enc_szip, std::slice::from_ref(&szip), 0).unwrap();
        assert_eq!(dec_szip, raw);
    }
}
