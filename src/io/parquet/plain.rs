//! Parquet value encodings (`PLAIN`, RLE/bit-packed hybrid, `DELTA_BINARY_PACKED`,
//! `DELTA_LENGTH_BYTE_ARRAY`, `BYTE_STREAM_SPLIT`) and level codecs.

use super::thrift::write_varint;

/// Julian day number of the Unix epoch (1970-01-01).
const JULIAN_UNIX_EPOCH: i64 = 2_440_588;
/// Nanoseconds in one 24-hour day.
const NANOS_PER_DAY: i64 = 86_400_000_000_000;

/// Encodes booleans as a little-endian bit-packed byte string.
pub(super) fn encode_bool(values: &[bool]) -> Vec<u8> {
    let mut out = vec![0u8; values.len().div_ceil(8)];
    for (index, value) in values.iter().enumerate() {
        if *value {
            if let Some(byte) = out.get_mut(index / 8) {
                *byte |= 1 << (index % 8);
            }
        }
    }
    out
}

/// Decodes `count` bit-packed booleans.
pub(super) fn decode_bool(
    bytes: &[u8],
    count: usize,
) -> Result<Vec<bool>, String> {
    let needed = count.div_ceil(8);
    if bytes.len() < needed {
        return Err("boolean page is too short".to_string());
    }
    let mut out = Vec::with_capacity(count);
    for index in 0..count {
        let byte = bytes.get(index / 8).copied().unwrap_or(0);
        out.push((byte >> (index % 8)) & 1 == 1);
    }
    Ok(out)
}

/// Encodes little-endian `i32` values.
pub(super) fn encode_i32(values: &[i32]) -> Vec<u8> {
    let mut out = Vec::with_capacity(values.len().saturating_mul(4));
    for value in values {
        out.extend_from_slice(&value.to_le_bytes());
    }
    out
}

/// Encodes little-endian `i64` values.
pub(super) fn encode_i64(values: &[i64]) -> Vec<u8> {
    let mut out = Vec::with_capacity(values.len().saturating_mul(8));
    for value in values {
        out.extend_from_slice(&value.to_le_bytes());
    }
    out
}

/// Encodes little-endian `f32` values.
pub(super) fn encode_f32(values: &[f32]) -> Vec<u8> {
    let mut out = Vec::with_capacity(values.len().saturating_mul(4));
    for value in values {
        out.extend_from_slice(&value.to_le_bytes());
    }
    out
}

/// Encodes little-endian `f64` values.
pub(super) fn encode_f64(values: &[f64]) -> Vec<u8> {
    let mut out = Vec::with_capacity(values.len().saturating_mul(8));
    for value in values {
        out.extend_from_slice(&value.to_le_bytes());
    }
    out
}

/// Encodes byte arrays as a 4-byte length followed by the bytes.
pub(super) fn encode_bytes(values: &[Vec<u8>]) -> Vec<u8> {
    let mut out = Vec::new();
    for value in values {
        out.extend_from_slice(&(value.len() as u32).to_le_bytes());
        out.extend_from_slice(value);
    }
    out
}

/// Encodes Unix nanoseconds (`i64`) as 12-byte legacy Parquet `INT96` values
/// (8-byte little-endian nanoseconds-of-day followed by 4-byte little-endian Julian day).
pub(super) fn encode_int96(values: &[i64]) -> Result<Vec<u8>, String> {
    let mut out = Vec::with_capacity(values.len().saturating_mul(12));
    for &unix_nanos in values {
        let days = unix_nanos.div_euclid(NANOS_PER_DAY);
        let nanos_of_day = unix_nanos.rem_euclid(NANOS_PER_DAY) as u64;
        let julian = days
            .checked_add(JULIAN_UNIX_EPOCH)
            .and_then(|d| u32::try_from(d).ok())
            .ok_or_else(|| "int96 timestamp is out of Julian day range".to_string())?;
        out.extend_from_slice(&nanos_of_day.to_le_bytes());
        out.extend_from_slice(&julian.to_le_bytes());
    }
    Ok(out)
}

/// Decodes `count` 12-byte `INT96` values into Unix nanoseconds (`i64`).
pub(super) fn decode_int96(
    bytes: &[u8],
    count: usize,
) -> Result<Vec<i64>, String> {
    let needed = count
        .checked_mul(12)
        .ok_or_else(|| "column is too large".to_string())?;
    if bytes.len() < needed {
        return Err("int96 page is too short".to_string());
    }
    let mut out = Vec::with_capacity(count);
    for chunk in bytes.as_chunks::<12>().0.iter().take(count) {
        let nanos = u64::from_le_bytes([
            chunk[0], chunk[1], chunk[2], chunk[3], chunk[4], chunk[5], chunk[6], chunk[7],
        ]);
        let julian = u32::from_le_bytes([chunk[8], chunk[9], chunk[10], chunk[11]]);
        let total = (i128::from(julian) - i128::from(JULIAN_UNIX_EPOCH))
            * i128::from(NANOS_PER_DAY)
            + i128::from(nanos);
        let unix_nanos = i64::try_from(total)
            .map_err(|_| "int96 timestamp overflows i64 nanoseconds".to_string())?;
        out.push(unix_nanos);
    }
    Ok(out)
}

/// Encodes fixed-length byte arrays without length prefixes.
pub(super) fn encode_fixed_bytes(
    values: &[Vec<u8>],
    width: usize,
) -> Result<Vec<u8>, String> {
    if width == 0 {
        return Err("fixed-length byte array width must be positive".to_string());
    }
    let mut out = Vec::with_capacity(values.len().saturating_mul(width));
    for value in values {
        if value.len() != width {
            return Err(format!(
                "expected fixed-length byte array of {width} bytes, got {}",
                value.len()
            ));
        }
        out.extend_from_slice(value);
    }
    Ok(out)
}

/// Decodes `count` fixed-length byte arrays of `width` bytes each.
pub(super) fn decode_fixed_bytes(
    bytes: &[u8],
    count: usize,
    width: usize,
) -> Result<Vec<Vec<u8>>, String> {
    if width == 0 {
        return Err("fixed-length byte array width must be positive".to_string());
    }
    let needed = count
        .checked_mul(width)
        .ok_or_else(|| "fixed-length byte array column is too large".to_string())?;
    if bytes.len() < needed {
        return Err("fixed-length byte array page is too short".to_string());
    }
    let mut out = Vec::with_capacity(count.min(1_000_000));
    for index in 0..count {
        let start = index * width;
        out.push(bytes[start..start + width].to_vec());
    }
    Ok(out)
}

/// Reads a little-endian `u32` at `pos`.
fn read_u32(
    bytes: &[u8],
    pos: usize,
) -> Result<u32, String> {
    let end = pos
        .checked_add(4)
        .ok_or_else(|| "offset overflows".to_string())?;
    let slice = bytes
        .get(pos..end)
        .ok_or_else(|| "page is truncated".to_string())?;
    let mut raw = [0u8; 4];
    raw.copy_from_slice(slice);
    Ok(u32::from_le_bytes(raw))
}

/// Decodes `count` little-endian `i32` values.
pub(super) fn decode_i32(
    bytes: &[u8],
    count: usize,
) -> Result<Vec<i32>, String> {
    let needed = count
        .checked_mul(4)
        .ok_or_else(|| "column is too large".to_string())?;
    if bytes.len() < needed {
        return Err("int32 page is too short".to_string());
    }
    let mut out = Vec::with_capacity(count);
    for chunk in bytes.as_chunks::<4>().0.iter().take(count) {
        out.push(i32::from_le_bytes(*chunk));
    }
    Ok(out)
}

/// Decodes `count` little-endian `i64` values.
pub(super) fn decode_i64(
    bytes: &[u8],
    count: usize,
) -> Result<Vec<i64>, String> {
    let needed = count
        .checked_mul(8)
        .ok_or_else(|| "column is too large".to_string())?;
    if bytes.len() < needed {
        return Err("int64 page is too short".to_string());
    }
    let mut out = Vec::with_capacity(count);
    for chunk in bytes.as_chunks::<8>().0.iter().take(count) {
        out.push(i64::from_le_bytes(*chunk));
    }
    Ok(out)
}

/// Decodes `count` little-endian `f32` values.
pub(super) fn decode_f32(
    bytes: &[u8],
    count: usize,
) -> Result<Vec<f32>, String> {
    let needed = count
        .checked_mul(4)
        .ok_or_else(|| "column is too large".to_string())?;
    if bytes.len() < needed {
        return Err("float page is too short".to_string());
    }
    let mut out = Vec::with_capacity(count);
    for chunk in bytes.as_chunks::<4>().0.iter().take(count) {
        out.push(f32::from_le_bytes(*chunk));
    }
    Ok(out)
}

/// Decodes `count` little-endian `f64` values.
pub(super) fn decode_f64(
    bytes: &[u8],
    count: usize,
) -> Result<Vec<f64>, String> {
    let needed = count
        .checked_mul(8)
        .ok_or_else(|| "column is too large".to_string())?;
    if bytes.len() < needed {
        return Err("double page is too short".to_string());
    }
    let mut out = Vec::with_capacity(count);
    for chunk in bytes.as_chunks::<8>().0.iter().take(count) {
        out.push(f64::from_le_bytes(*chunk));
    }
    Ok(out)
}

/// Decodes `count` length-prefixed byte arrays.
pub(super) fn decode_bytes(
    bytes: &[u8],
    count: usize,
) -> Result<Vec<&[u8]>, String> {
    let mut out = Vec::with_capacity(count.min(1_000_000));
    let mut pos = 0usize;
    for _ in 0..count {
        let len = read_u32(bytes, pos)? as usize;
        pos = pos
            .checked_add(4)
            .ok_or_else(|| "byte array offset overflows".to_string())?;
        let end = pos
            .checked_add(len)
            .ok_or_else(|| "byte array length overflows".to_string())?;
        let value = bytes
            .get(pos..end)
            .ok_or_else(|| "byte array exceeds the page".to_string())?;
        out.push(value);
        pos = end;
    }
    Ok(out)
}

/// Encodes definition levels as RLE runs of a single bit.
///
/// The caller prefixes the returned bytes with their 4-byte length, as the
/// v1 data page layout requires.
pub(super) fn encode_levels(levels: &[bool]) -> Vec<u8> {
    let values: Vec<u32> = levels.iter().map(|&level| u32::from(level)).collect();
    encode_level_values(&values)
}

/// Encodes definition or repetition levels as RLE runs.
///
/// Every level fits in one byte, which is what the RLE run header expects for
/// a bit width of at most eight.
pub(super) fn encode_level_values(levels: &[u32]) -> Vec<u8> {
    let mut out = Vec::new();
    let mut index = 0usize;
    while index < levels.len() {
        let value = levels[index];
        let mut end = index + 1;
        while levels.get(end).copied() == Some(value) {
            end += 1;
        }
        let run = end.saturating_sub(index);
        write_varint(&mut out, (run as u64) << 1);
        out.push(u8::try_from(value).unwrap_or(0));
        index = end;
    }
    out
}

/// Encodes dictionary indices with the RLE/bit-packed hybrid used by
/// `RLE_DICTIONARY` data pages.
///
/// Indices are bit-packed in groups of eight; a trailing group is zero-padded,
/// which the reader discards because it stops after the declared value count.
pub(super) fn encode_indices(
    indices: &[u32],
    bit_width: u32,
) -> Vec<u8> {
    let mut out = Vec::new();
    if indices.is_empty() || bit_width == 0 {
        return out;
    }
    let width = bit_width.min(32);
    let groups = indices.len().div_ceil(8);
    write_varint(&mut out, ((groups as u64) << 1) | 1);
    let mask = if width >= 32 {
        u64::from(u32::MAX)
    } else {
        (1u64 << width) - 1
    };
    let mut accumulator = 0u64;
    let mut bits = 0u32;
    for index in 0..groups * 8 {
        let value = u64::from(indices.get(index).copied().unwrap_or(0)) & mask;
        accumulator |= value << bits;
        bits += width;
        while bits >= 8 {
            out.push((accumulator & 0xff) as u8);
            accumulator >>= 8;
            bits -= 8;
        }
    }
    if bits > 0 {
        out.push((accumulator & 0xff) as u8);
    }
    out
}

/// Decodes definition levels written with the RLE/bit-packed hybrid.
///
/// Only bit width 1 (the flat-schema case) is supported.
#[cfg(test)]
pub(super) fn decode_levels(
    bytes: &[u8],
    bit_width: u32,
    count: usize,
) -> Result<Vec<bool>, String> {
    if bit_width == 0 {
        return Ok(vec![false; count]);
    }
    if bit_width != 1 {
        return Err("only single-bit definition levels are supported".to_string());
    }
    let mut out = Vec::with_capacity(count);
    let mut pos = 0usize;
    while out.len() < count {
        let (header, next) = read_level_varint(bytes, pos)?;
        pos = next;
        if header & 1 == 0 {
            let run =
                usize::try_from(header >> 1).map_err(|_| "level run is too large".to_string())?;
            if run == 0 {
                return Err("zero-length level run".to_string());
            }
            let value = *bytes
                .get(pos)
                .ok_or_else(|| "level run value is truncated".to_string())?;
            pos += 1;
            let take = run.min(count.saturating_sub(out.len()));
            for _ in 0..take {
                out.push(value & 1 == 1);
            }
        } else {
            let groups = usize::try_from(header >> 1)
                .map_err(|_| "level group count is too large".to_string())?;
            if groups == 0 {
                return Err("zero-length bit-packed level run".to_string());
            }
            let packed = bytes
                .get(pos..pos.saturating_add(groups))
                .ok_or_else(|| "bit-packed levels are truncated".to_string())?;
            pos += groups;
            let total = groups.saturating_mul(8);
            let take = total.min(count.saturating_sub(out.len()));
            for index in 0..take {
                let byte = packed.get(index / 8).copied().unwrap_or(0);
                out.push((byte >> (index % 8)) & 1 == 1);
            }
        }
    }
    Ok(out)
}

/// Reads a varint from the level stream.
fn read_level_varint(
    bytes: &[u8],
    start: usize,
) -> Result<(u64, usize), String> {
    let mut result = 0u64;
    let mut shift = 0u32;
    let mut pos = start;
    loop {
        if shift >= 64 {
            return Err("level varint is too long".to_string());
        }
        let byte = *bytes
            .get(pos)
            .ok_or_else(|| "level data is truncated".to_string())?;
        pos += 1;
        result |= u64::from(byte & 0x7f) << shift;
        if byte & 0x80 == 0 {
            return Ok((result, pos));
        }
        shift += 7;
    }
}

/// Reads a zigzag-encoded signed varint from `bytes` starting at `start`.
fn read_zigzag_varint(
    bytes: &[u8],
    start: usize,
) -> Result<(i64, usize), String> {
    let (raw, next) = read_level_varint(bytes, start)?;
    Ok((((raw >> 1) as i64) ^ -((raw & 1) as i64), next))
}

/// Unpacks `bit_width` bits (0..=64) in LSB-first order starting at `bit_offset`.
fn unpack_bits_u64(
    bytes: &[u8],
    bit_offset: usize,
    bit_width: u32,
) -> u64 {
    let mut value = 0u64;
    let mut bits_read = 0u32;
    while bits_read < bit_width {
        let pos = bit_offset + bits_read as usize;
        let byte_idx = pos / 8;
        let bit_in_byte = (pos % 8) as u32;
        let avail = 8 - bit_in_byte;
        let take = avail.min(bit_width - bits_read);
        let mask = ((1u16 << take) - 1) as u8;
        let chunk = (bytes.get(byte_idx).copied().unwrap_or(0) >> bit_in_byte) & mask;
        value |= u64::from(chunk) << bits_read;
        bits_read += take;
    }
    value
}

/// Decodes `count` unsigned integers of `bit_width` bits (0..=32) from an
/// RLE/bit-packed hybrid stream (used for multi-bit definition levels and
/// dictionary indices).
pub(super) fn decode_rle_hybrid(
    bytes: &[u8],
    bit_width: u32,
    count: usize,
) -> Result<Vec<u32>, String> {
    if count == 0 {
        return Ok(Vec::new());
    }
    if bit_width == 0 {
        return Ok(vec![0u32; count]);
    }
    if bit_width > 32 {
        return Err(format!("rle hybrid bit width {bit_width} exceeds 32"));
    }
    let byte_width = (bit_width as usize).div_ceil(8);
    let mut out = Vec::with_capacity(count.min(1_000_000));
    let mut pos = 0usize;
    while out.len() < count {
        let (header, next) = read_level_varint(bytes, pos)?;
        pos = next;
        if header & 1 == 0 {
            let run = usize::try_from(header >> 1)
                .map_err(|_| "rle run length is too large".to_string())?;
            if run == 0 {
                return Err("zero-length rle run".to_string());
            }
            let end = pos
                .checked_add(byte_width)
                .ok_or_else(|| "rle value offset overflows".to_string())?;
            let val_bytes = bytes
                .get(pos..end)
                .ok_or_else(|| "rle run value is truncated".to_string())?;
            pos = end;
            let mut value = 0u32;
            for (index, &byte) in val_bytes.iter().enumerate() {
                value |= u32::from(byte) << (8 * index);
            }
            let take = run.min(count.saturating_sub(out.len()));
            out.extend(std::iter::repeat_n(value, take));
        } else {
            let groups = usize::try_from(header >> 1)
                .map_err(|_| "bit-packed group count is too large".to_string())?;
            if groups == 0 {
                return Err("zero-length bit-packed run".to_string());
            }
            let packed_bytes = groups
                .checked_mul(bit_width as usize)
                .ok_or_else(|| "bit-packed run length overflows".to_string())?;
            let end = pos
                .checked_add(packed_bytes)
                .ok_or_else(|| "bit-packed run offset overflows".to_string())?;
            let packed = bytes
                .get(pos..end)
                .ok_or_else(|| "bit-packed run is truncated".to_string())?;
            pos = end;
            let total = groups
                .checked_mul(8)
                .ok_or_else(|| "bit-packed value count overflows".to_string())?;
            let take = total.min(count.saturating_sub(out.len()));
            for index in 0..take {
                let val = unpack_bits_u64(packed, index * bit_width as usize, bit_width) as u32;
                out.push(val);
            }
        }
    }
    Ok(out)
}

/// Decodes a `DELTA_BINARY_PACKED` integer stream, returning the decoded
/// values and the number of bytes consumed from `bytes`.
pub(super) fn decode_delta_binary_packed(
    bytes: &[u8],
    expected_count: usize,
) -> Result<(Vec<i64>, usize), String> {
    let (block_size_u64, pos) = read_level_varint(bytes, 0)?;
    let (miniblocks_u64, pos) = read_level_varint(bytes, pos)?;
    let (total_u64, pos) = read_level_varint(bytes, pos)?;
    let (first_value, mut pos) = read_zigzag_varint(bytes, pos)?;

    let block_size = usize::try_from(block_size_u64)
        .map_err(|_| "delta block size is out of range".to_string())?;
    let miniblocks_per_block = usize::try_from(miniblocks_u64)
        .map_err(|_| "delta miniblock count is out of range".to_string())?;
    let total_count = usize::try_from(total_u64)
        .map_err(|_| "delta total value count is out of range".to_string())?;

    if total_count != expected_count {
        return Err(format!(
            "delta stream declares {total_count} values, expected {expected_count}"
        ));
    }
    if total_count == 0 {
        return Ok((Vec::new(), pos));
    }
    if block_size == 0 || miniblocks_per_block == 0 || block_size % miniblocks_per_block != 0 {
        return Err("invalid delta block geometry".to_string());
    }
    let values_per_miniblock = block_size / miniblocks_per_block;
    if values_per_miniblock == 0 || values_per_miniblock % 8 != 0 {
        return Err("invalid delta miniblock size".to_string());
    }

    let mut out = Vec::with_capacity(total_count.min(1_000_000));
    out.push(first_value);
    let mut current = first_value;

    while out.len() < total_count {
        let (min_delta, next) = read_zigzag_varint(bytes, pos)?;
        pos = next;
        let bw_end = pos
            .checked_add(miniblocks_per_block)
            .ok_or_else(|| "delta bit-width table overflows".to_string())?;
        let bit_widths = bytes
            .get(pos..bw_end)
            .ok_or_else(|| "delta bit-width table is truncated".to_string())?;
        pos = bw_end;

        for &bw in bit_widths {
            if out.len() >= total_count {
                break;
            }
            if bw > 64 {
                return Err(format!("delta miniblock bit width {bw} exceeds 64"));
            }
            let mb_bytes = values_per_miniblock
                .checked_mul(usize::from(bw))
                .ok_or_else(|| "delta miniblock size overflows".to_string())?
                / 8;
            let mb_end = pos
                .checked_add(mb_bytes)
                .ok_or_else(|| "delta miniblock offset overflows".to_string())?;
            let mb_slice = bytes
                .get(pos..mb_end)
                .ok_or_else(|| "delta miniblock data is truncated".to_string())?;
            pos = mb_end;

            let take = values_per_miniblock.min(total_count - out.len());
            for index in 0..take {
                let delta_offset =
                    unpack_bits_u64(mb_slice, index * usize::from(bw), u32::from(bw));
                let delta = min_delta.wrapping_add(delta_offset as i64);
                current = current.wrapping_add(delta);
                out.push(current);
            }
        }
    }

    Ok((out, pos))
}

/// Decodes a `DELTA_LENGTH_BYTE_ARRAY` page into raw byte arrays.
pub(super) fn decode_delta_length_byte_array(
    bytes: &[u8],
    count: usize,
) -> Result<Vec<Vec<u8>>, String> {
    let (values, _) = decode_delta_length_byte_array_at(bytes, count)?;
    Ok(values)
}

/// Decodes `DELTA_LENGTH_BYTE_ARRAY` and reports how many bytes it consumed.
pub(super) fn decode_delta_length_byte_array_at(
    bytes: &[u8],
    count: usize,
) -> Result<(Vec<Vec<u8>>, usize), String> {
    if count == 0 {
        return Ok((Vec::new(), 0));
    }
    let (lengths, mut pos) = decode_delta_binary_packed(bytes, count)?;
    let mut out = Vec::with_capacity(count.min(1_000_000));
    for len in lengths {
        let len = usize::try_from(len)
            .map_err(|_| "negative length in delta-length byte array".to_string())?;
        let end = pos
            .checked_add(len)
            .ok_or_else(|| "delta-length byte array offset overflows".to_string())?;
        let slice = bytes
            .get(pos..end)
            .ok_or_else(|| "delta-length byte array payload is truncated".to_string())?;
        out.push(slice.to_vec());
        pos = end;
    }
    Ok((out, pos))
}

/// Decodes `DELTA_BYTE_ARRAY`: prefix lengths followed by the suffixes.
///
/// Each value is the first `prefix[i]` bytes of the previous value followed by
/// `suffix[i]`, which is how sorted strings compress well.
pub(super) fn decode_delta_byte_array(
    bytes: &[u8],
    count: usize,
) -> Result<Vec<Vec<u8>>, String> {
    if count == 0 {
        return Ok(Vec::new());
    }
    let (prefixes, consumed) = decode_delta_binary_packed(bytes, count)?;
    let (suffixes, _) = decode_delta_length_byte_array_at(
        bytes
            .get(consumed..)
            .ok_or("delta byte array is truncated")?,
        count,
    )?;
    let mut out: Vec<Vec<u8>> = Vec::with_capacity(count.min(1_000_000));
    for (index, suffix) in suffixes.into_iter().enumerate() {
        let prefix = usize::try_from(prefixes.get(index).copied().unwrap_or(0)).unwrap_or(0);
        let previous: &[u8] = out.last().map_or(&[], Vec::as_slice);
        let take = prefix.min(previous.len());
        let mut value = Vec::with_capacity(take.saturating_add(suffix.len()));
        value.extend_from_slice(previous.get(..take).unwrap_or(&[]));
        value.extend_from_slice(&suffix);
        out.push(value);
    }
    Ok(out)
}

/// De-interleaves a `BYTE_STREAM_SPLIT` buffer of `count` values of `width`
/// bytes each into standard `PLAIN` little-endian layout.
pub(super) fn decode_byte_stream_split(
    bytes: &[u8],
    count: usize,
    width: usize,
) -> Result<Vec<u8>, String> {
    let total = count
        .checked_mul(width)
        .ok_or_else(|| "byte stream split size overflows".to_string())?;
    if bytes.len() < total {
        return Err("byte stream split page is too short".to_string());
    }
    let mut out = vec![0u8; total];
    for byte_idx in 0..width {
        let stream_start = byte_idx * count;
        for item_idx in 0..count {
            out[item_idx * width + byte_idx] = bytes[stream_start + item_idx];
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fixed_width_round_trips() {
        let ints = vec![1i32, -2, 3];
        assert_eq!(decode_i32(&encode_i32(&ints), 3).expect("i32"), ints);
        let longs = vec![1i64, -2, i64::MAX];
        assert_eq!(decode_i64(&encode_i64(&longs), 3).expect("i64"), longs);
        let floats = vec![1.5f32, -2.25];
        assert_eq!(decode_f32(&encode_f32(&floats), 2).expect("f32"), floats);
        let doubles = vec![1.5f64, -2.25];
        assert_eq!(decode_f64(&encode_f64(&doubles), 2).expect("f64"), doubles);
    }

    #[test]
    fn booleans_are_bit_packed() {
        let values = vec![true, false, true, true, false, false, true, false, true];
        let encoded = encode_bool(&values);
        assert_eq!(encoded.len(), 2);
        assert_eq!(decode_bool(&encoded, values.len()).expect("decode"), values);
        assert!(decode_bool(&encoded, 100).is_err());
    }

    #[test]
    fn byte_arrays_round_trip() {
        let values = vec![b"alpha".to_vec(), Vec::new(), vec![0, 255]];
        let encoded = encode_bytes(&values);
        let decoded = decode_bytes(&encoded, 3).expect("decode");
        assert_eq!(decoded, vec![&b"alpha"[..], &b""[..], &[0u8, 255][..]]);
        assert!(decode_bytes(&encoded, 4).is_err());
    }

    #[test]
    fn levels_round_trip() {
        let levels = vec![false, false, true, true, true, false, true];
        let encoded = encode_levels(&levels);
        assert_eq!(
            decode_levels(&encoded, 1, levels.len()).expect("decode"),
            levels
        );
        let all = vec![true; 13];
        assert_eq!(
            decode_levels(&encode_levels(&all), 1, all.len()).expect("decode"),
            all
        );
        assert!(decode_levels(&[], 1, 1).is_err());
        assert!(decode_levels(&encoded, 2, 1).is_err());
    }

    #[test]
    fn bit_packed_levels_are_accepted() {
        // One group of eight bit-packed levels, header = (1 << 1) | 1.
        let bytes = [3u8, 0b0000_0101];
        let levels = decode_levels(&bytes, 1, 8).expect("decode");
        assert_eq!(
            levels,
            vec![true, false, true, false, false, false, false, false]
        );
    }

    #[test]
    fn int96_and_fixed_bytes_round_trip() {
        let timestamps = vec![0i64, -1, 1_700_000_000_123_456_789, -86_400_000_000_000];
        let encoded = encode_int96(&timestamps).expect("encode");
        assert_eq!(encoded.len(), timestamps.len() * 12);
        assert_eq!(
            decode_int96(&encoded, timestamps.len()).expect("decode"),
            timestamps
        );

        let fixed = vec![vec![1, 2, 3, 4], vec![5, 6, 7, 8]];
        let enc_fixed = encode_fixed_bytes(&fixed, 4).expect("encode fixed");
        assert_eq!(
            decode_fixed_bytes(&enc_fixed, 2, 4).expect("decode fixed"),
            fixed
        );
        assert!(encode_fixed_bytes(&fixed, 3).is_err());
    }

    #[test]
    fn rle_hybrid_multibit_decodes() {
        // RLE run of 3 copies of value 2 (bit_width = 2 -> 1 byte value),
        // followed by 1 bit-packed group of 8 2-bit values: [0,1,2,3, 3,2,1,0].
        // Packed bytes: byte0 = 0 | (1<<2) | (2<<4) | (3<<6) = 0b11_10_01_00 = 0xE4
        //               byte1 = 3 | (2<<2) | (1<<4) | (0<<6) = 0b00_01_10_11 = 0x1B
        let stream = [
            (3u8 << 1),     // RLE run of 3
            2u8,            // value 2
            (1u8 << 1) | 1, // 1 bit-packed group (8 values)
            0xE4,
            0x1B,
        ];
        let decoded = decode_rle_hybrid(&stream, 2, 11).expect("decode");
        assert_eq!(decoded, vec![2, 2, 2, 0, 1, 2, 3, 3, 2, 1, 0]);
    }

    #[test]
    fn delta_binary_packed_and_delta_length_decode() {
        use super::super::thrift::write_zigzag;
        // Values: [10, 12, 15, 14, 20] -> deltas: [2, 3, -1, 6], min_delta = -1,
        // offsets: [3, 4, 0, 7] -> fits in 3 bits each.
        let mut buf = Vec::new();
        write_varint(&mut buf, 128); // block_size
        write_varint(&mut buf, 4); // miniblocks_per_block -> 32 values/miniblock
        write_varint(&mut buf, 5); // total_value_count
        write_zigzag(&mut buf, 10); // first_value
        write_zigzag(&mut buf, -1); // min_delta
        buf.extend_from_slice(&[3, 0, 0, 0]); // bit_widths for 4 miniblocks
        // Miniblock 0: 32 * 3 / 8 = 12 bytes, first 4 values are [3, 4, 0, 7].
        let mut mb = [0u8; 12];
        let vals = [3u64, 4, 0, 7];
        for (i, &v) in vals.iter().enumerate() {
            let bit_pos = i * 3;
            for b in 0..3 {
                if (v >> b) & 1 == 1 {
                    mb[(bit_pos + b) / 8] |= 1 << ((bit_pos + b) % 8);
                }
            }
        }
        buf.extend_from_slice(&mb);
        let (decoded, used) = decode_delta_binary_packed(&buf, 5).expect("delta");
        assert_eq!(decoded, vec![10, 12, 15, 14, 20]);
        assert_eq!(used, buf.len());

        // DELTA_LENGTH_BYTE_ARRAY for ["hi", "rust", "!"]: lengths [2, 4, 1], deltas [2, -3], min=-3, offsets [5, 0] (3 bits).
        let mut dl = Vec::new();
        write_varint(&mut dl, 128);
        write_varint(&mut dl, 4);
        write_varint(&mut dl, 3);
        write_zigzag(&mut dl, 2);
        write_zigzag(&mut dl, -3);
        dl.extend_from_slice(&[3, 0, 0, 0]);
        let mut mb2 = [0u8; 12];
        // offset[0] = 5 (0b101), offset[1] = 0
        mb2[0] = 0b101;
        dl.extend_from_slice(&mb2);
        dl.extend_from_slice(b"hirust!");
        let items = decode_delta_length_byte_array(&dl, 3).expect("delta len");
        assert_eq!(items, vec![b"hi".to_vec(), b"rust".to_vec(), b"!".to_vec()]);
    }

    #[test]
    fn byte_stream_split_decodes() {
        let floats = [1.5f32, -2.25, 3.75];
        let raw: Vec<[u8; 4]> = floats.iter().map(|v| v.to_le_bytes()).collect();
        let mut interleaved = Vec::new();
        for b in 0..4 {
            for item in &raw {
                interleaved.push(item[b]);
            }
        }
        let plain_bytes = decode_byte_stream_split(&interleaved, 3, 4).expect("bss");
        let decoded = decode_f32(&plain_bytes, 3).expect("f32");
        assert_eq!(decoded, floats);
    }
}
