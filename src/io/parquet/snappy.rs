//! Snappy block-format compression and decompression.
//!
//! Parquet stores Snappy-compressed pages as one block: a varint with the
//! uncompressed length followed by literal and copy tags.

use super::thrift::write_varint;

/// Hash table bit width for LZ77 match finding.
const HASH_BITS: u32 = 14;
/// Hash table entry count.
const HASH_SIZE: usize = 1 << HASH_BITS;

/// Hashes four bytes at `input[pos..pos + 4]` into a table slot.
fn hash4(
    input: &[u8],
    pos: usize,
) -> usize {
    let v = u32::from_le_bytes([input[pos], input[pos + 1], input[pos + 2], input[pos + 3]]);
    (v.wrapping_mul(0x1e35_a7bd) >> (32 - HASH_BITS)) as usize
}

/// Appends a Snappy literal element to `out`.
fn emit_literal(
    out: &mut Vec<u8>,
    literal: &[u8],
) {
    if literal.is_empty() {
        return;
    }
    let n = literal.len() - 1;
    if n < 60 {
        out.push((n as u8) << 2);
    } else if n <= 0xff {
        out.push(60 << 2);
        out.push(n as u8);
    } else if n <= 0xffff {
        out.push(61 << 2);
        out.extend_from_slice(&(n as u16).to_le_bytes());
    } else if n <= 0xff_ffff {
        out.push(62 << 2);
        out.push((n & 0xff) as u8);
        out.push(((n >> 8) & 0xff) as u8);
        out.push(((n >> 16) & 0xff) as u8);
    } else {
        out.push(63 << 2);
        out.extend_from_slice(&(n as u32).to_le_bytes());
    }
    out.extend_from_slice(literal);
}

/// Appends one or more Snappy copy tags for `(offset, len)` to `out`.
fn emit_copy(
    out: &mut Vec<u8>,
    offset: usize,
    mut len: usize,
) {
    while len > 0 {
        let chunk = len.min(64);
        if (4..=11).contains(&chunk) && offset <= 2047 {
            let tag = (((offset >> 8) as u8) << 5) | (((chunk - 4) as u8) << 2) | 1;
            out.push(tag);
            out.push((offset & 0xff) as u8);
        } else if offset <= 65535 {
            let tag = (((chunk - 1) as u8) << 2) | 2;
            out.push(tag);
            out.extend_from_slice(&(offset as u16).to_le_bytes());
        } else {
            let tag = (((chunk - 1) as u8) << 2) | 3;
            out.push(tag);
            out.extend_from_slice(&(offset as u32).to_le_bytes());
        }
        len -= chunk;
    }
}

/// Compresses `input` into a single Snappy block.
pub(crate) fn compress(input: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(input.len() / 2 + 16);
    write_varint(&mut out, input.len() as u64);
    if input.is_empty() {
        return out;
    }
    let mut table = vec![usize::MAX; HASH_SIZE];
    let mut next_emit = 0usize;
    let mut ip = 0usize;
    while ip + 4 <= input.len() {
        let h = hash4(input, ip);
        let candidate = table[h];
        table[h] = ip;
        if candidate != usize::MAX
            && ip > candidate
            && ip - candidate <= 65535
            && input[candidate..candidate + 4] == input[ip..ip + 4]
        {
            emit_literal(&mut out, &input[next_emit..ip]);
            let offset = ip - candidate;
            let mut match_len = 4usize;
            while ip + match_len < input.len()
                && input[candidate + match_len] == input[ip + match_len]
            {
                match_len += 1;
            }
            emit_copy(&mut out, offset, match_len);
            ip += match_len;
            next_emit = ip;
        } else {
            ip += 1;
        }
    }
    if next_emit < input.len() {
        emit_literal(&mut out, &input[next_emit..]);
    }
    out
}

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
pub(crate) fn decompress(input: &[u8]) -> Result<Vec<u8>, String> {
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

    #[test]
    fn compress_round_trip_various_inputs() {
        // Empty, tiny, and short literal inputs.
        for sample in [
            &b""[..],
            &b"a"[..],
            &b"abc"[..],
            &b"hello world"[..],
            &b"abababababababababababab"[..],
        ] {
            let compressed = compress(sample);
            assert_eq!(decompress(&compressed).expect("decompress"), sample);
        }

        // Long repetitive input should compress significantly.
        let repeated = b"The quick brown fox jumps over the lazy dog. ".repeat(200);
        let compressed = compress(&repeated);
        assert!(compressed.len() < repeated.len() / 4);
        assert_eq!(decompress(&compressed).expect("decompress"), repeated);

        // Long non-repetitive literal (> 256 bytes to exercise 2-byte literal length).
        let mut pseudo_random = Vec::with_capacity(4096);
        let mut state = 0x1234_5678u32;
        for _ in 0..4096 {
            state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
            pseudo_random.push((state >> 24) as u8);
        }
        let compressed_rand = compress(&pseudo_random);
        assert_eq!(
            decompress(&compressed_rand).expect("decompress"),
            pseudo_random
        );
    }
}
