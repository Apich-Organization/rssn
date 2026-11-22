//! Parquet `PLAIN` encoding and the RLE/bit-packed definition-level codec.

use super::thrift::write_varint;

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
    let mut out = Vec::new();
    let mut index = 0usize;
    while index < levels.len() {
        let value = levels.get(index).copied().unwrap_or(false);
        let mut end = index + 1;
        while levels.get(end).copied() == Some(value) {
            end += 1;
        }
        let run = end.saturating_sub(index);
        write_varint(&mut out, (run as u64) << 1);
        out.push(u8::from(value));
        index = end;
    }
    out
}

/// Decodes definition levels written with the RLE/bit-packed hybrid.
///
/// Only bit width 1 (the flat-schema case) is supported.
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
}
