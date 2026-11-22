//! Snappy block-format decompression.
//!
//! Parquet stores Snappy-compressed pages as one block: a varint with the
//! uncompressed length followed by literal and copy tags. Compression is not
//! implemented because the writer always emits uncompressed pages.

/// Reads a LEB128 varint from `bytes` starting at `start`.
fn read_varint(
    bytes: &[u8],
    start: usize,
) -> Result<(u64, usize), String> {
    let mut result = 0u64;
    let mut shift = 0u32;
    let mut pos = start;
    loop {
        if shift >= 64 {
            return Err("snappy varint is too long".to_string());
        }
        let byte = *bytes
            .get(pos)
            .ok_or_else(|| "snappy block is truncated".to_string())?;
        pos += 1;
        result |= u64::from(byte & 0x7f) << shift;
        if byte & 0x80 == 0 {
            return Ok((result, pos));
        }
        shift += 7;
    }
}

/// Copies `len` bytes from `offset` bytes back, byte by byte so overlapping
/// runs expand correctly.
fn copy_back(
    out: &mut Vec<u8>,
    offset: usize,
    len: usize,
    expected: usize,
) -> Result<(), String> {
    if offset == 0 || offset > out.len() {
        return Err("snappy copy offset is out of bounds".to_string());
    }
    let total = out
        .len()
        .checked_add(len)
        .ok_or_else(|| "snappy output is too large".to_string())?;
    if total > expected {
        return Err("snappy output exceeds the declared length".to_string());
    }
    let start = out.len() - offset;
    for step in 0..len {
        let byte = *out
            .get(start + step)
            .ok_or_else(|| "snappy copy reads outside the output".to_string())?;
        out.push(byte);
    }
    Ok(())
}

/// Decompresses one Snappy block.
///
/// # Errors
///
/// Returns an error when the block is truncated, declares an impossible
/// length, or contains an invalid tag.
pub(super) fn decompress(input: &[u8]) -> Result<Vec<u8>, String> {
    let (declared, mut pos) = read_varint(input, 0)?;
    let expected = usize::try_from(declared)
        .map_err(|_| "snappy length does not fit in memory".to_string())?;
    if expected > 1_000_000_000 {
        return Err("snappy block is too large".to_string());
    }
    let mut out = Vec::with_capacity(expected.min(1 << 20));
    while pos < input.len() {
        let tag = *input
            .get(pos)
            .ok_or_else(|| "snappy block is truncated".to_string())?;
        pos += 1;
        match tag & 0x03 {
            | 0 => {
                let mut len = usize::from(tag >> 2);
                if len < 60 {
                    len += 1;
                } else {
                    let extra = len - 59;
                    if extra > 4 {
                        return Err("snappy literal length tag is invalid".to_string());
                    }
                    let bytes = input
                        .get(pos..pos.saturating_add(extra))
                        .ok_or_else(|| "snappy literal length is truncated".to_string())?;
                    pos += extra;
                    let mut value = 0usize;
                    for (index, byte) in bytes.iter().enumerate() {
                        value |= usize::from(*byte) << (8 * index);
                    }
                    len = value
                        .checked_add(1)
                        .ok_or_else(|| "snappy literal length overflows".to_string())?;
                }
                let end = pos
                    .checked_add(len)
                    .ok_or_else(|| "snappy literal length overflows".to_string())?;
                let literal = input
                    .get(pos..end)
                    .ok_or_else(|| "snappy literal exceeds the block".to_string())?;
                pos = end;
                if out
                    .len()
                    .checked_add(len)
                    .ok_or_else(|| "snappy output is too large".to_string())?
                    > expected
                {
                    return Err("snappy output exceeds the declared length".to_string());
                }
                out.extend_from_slice(literal);
            },
            | 1 => {
                let len = usize::from((tag >> 2) & 0x07) + 4;
                let low = *input
                    .get(pos)
                    .ok_or_else(|| "snappy copy is truncated".to_string())?;
                pos += 1;
                let offset = (usize::from(tag >> 5) << 8) | usize::from(low);
                copy_back(&mut out, offset, len, expected)?;
            },
            | 2 => {
                let len = usize::from(tag >> 2) + 1;
                let bytes = input
                    .get(pos..pos.saturating_add(2))
                    .ok_or_else(|| "snappy copy is truncated".to_string())?;
                pos += 2;
                let offset = usize::from(bytes.first().copied().unwrap_or(0))
                    | (usize::from(bytes.get(1).copied().unwrap_or(0)) << 8);
                copy_back(&mut out, offset, len, expected)?;
            },
            | 3 => {
                let len = usize::from(tag >> 2) + 1;
                let bytes = input
                    .get(pos..pos.saturating_add(4))
                    .ok_or_else(|| "snappy copy is truncated".to_string())?;
                pos += 4;
                let offset = u32::from_le_bytes([
                    bytes.first().copied().unwrap_or(0),
                    bytes.get(1).copied().unwrap_or(0),
                    bytes.get(2).copied().unwrap_or(0),
                    bytes.get(3).copied().unwrap_or(0),
                ]) as usize;
                copy_back(&mut out, offset, len, expected)?;
            },
            | _ => return Err("snappy tag is invalid".to_string()),
        }
    }
    if out.len() != expected {
        return Err("snappy block length does not match its header".to_string());
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn literal_block() {
        let mut block = vec![5];
        block.push((5 - 1) << 2);
        block.extend_from_slice(b"hello");
        assert_eq!(decompress(&block).expect("decompress"), b"hello");
    }

    #[test]
    fn overlapping_copy_expands() {
        let mut block = vec![6];
        block.push((2 - 1) << 2);
        block.extend_from_slice(b"ab");
        block.push(((4 - 1) << 2) | 2);
        block.extend_from_slice(&[2, 0]);
        assert_eq!(decompress(&block).expect("decompress"), b"ababab");
    }

    #[test]
    fn truncated_and_corrupt_blocks_fail() {
        assert!(decompress(&[]).is_err());
        assert!(decompress(&[5, 16, b'h']).is_err());
        let mut bad_copy = vec![4];
        bad_copy.push((1 - 1) << 2);
        bad_copy.push(b'a');
        bad_copy.push(((2 - 1) << 2) | 2);
        bad_copy.extend_from_slice(&[9, 0]);
        assert!(decompress(&bad_copy).is_err());
    }

    #[test]
    fn declared_length_must_match() {
        let mut block = vec![4];
        block.push((2 - 1) << 2);
        block.extend_from_slice(b"ab");
        assert!(decompress(&block).is_err());
    }
}
