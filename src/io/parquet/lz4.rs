//! LZ4 block-format compression and decompression.
//!
//! Parquet's `LZ4_RAW` codec stores one raw LZ4 block per page: a sequence of
//! literal runs and back-references, without a frame header, length prefix or
//! checksum. The compressor is a greedy four-byte hash matcher, which is what
//! the reference implementation's fast mode does.

/// Minimum match length in an LZ4 sequence.
const MIN_MATCH: usize = 4;
/// The last match must start at least this many bytes before the end.
const LAST_LITERALS: usize = 5;
/// Matches may not start in the final `MATCH_SAFETY` bytes.
const MATCH_SAFETY: usize = 12;
/// Largest representable back-reference distance.
const MAX_DISTANCE: usize = 65_535;
/// Hash-table bit width for match finding.
const HASH_BITS: u32 = 14;
/// Hash-table size.
const HASH_SIZE: usize = 1 << HASH_BITS;

/// Hashes four bytes into a table slot.
fn hash4(
    input: &[u8],
    pos: usize,
) -> usize {
    let v = u32::from_le_bytes([input[pos], input[pos + 1], input[pos + 2], input[pos + 3]]);
    (v.wrapping_mul(2_654_435_761) >> (32 - HASH_BITS)) as usize
}

/// Appends a length extension (`255` chunks then the remainder).
fn write_length(
    out: &mut Vec<u8>,
    mut length: usize,
) {
    while length >= 255 {
        out.push(255);
        length -= 255;
    }
    out.push(length as u8);
}

/// Emits one sequence: `literal` bytes from `anchor`, then a match.
fn write_sequence(
    out: &mut Vec<u8>,
    literals: &[u8],
    offset: usize,
    match_len: usize,
) {
    let lit_len = literals.len();
    let encoded_match = match_len - MIN_MATCH;
    let token = ((lit_len.min(15) as u8) << 4) | encoded_match.min(15) as u8;
    out.push(token);
    if lit_len >= 15 {
        write_length(out, lit_len - 15);
    }
    out.extend_from_slice(literals);
    out.extend_from_slice(&(offset as u16).to_le_bytes());
    if encoded_match >= 15 {
        write_length(out, encoded_match - 15);
    }
}

/// Compresses `input` into one LZ4 block.
pub fn compress(input: &[u8]) -> Vec<u8> {
    let len = input.len();
    let mut out = Vec::with_capacity(len / 2 + 16);
    if len == 0 {
        return out;
    }
    // Short inputs are emitted as a single literal run: a match would violate
    // the trailing-literal rules.
    if len < MATCH_SAFETY {
        let token = (len.min(15) as u8) << 4;
        out.push(token);
        if len >= 15 {
            write_length(&mut out, len - 15);
        }
        out.extend_from_slice(input);
        return out;
    }

    let mut table = vec![usize::MAX; HASH_SIZE];
    let mut anchor = 0usize;
    let mut pos = 0usize;
    let match_limit = len - MATCH_SAFETY;

    while pos < match_limit {
        let slot = hash4(input, pos);
        let candidate = table[slot];
        table[slot] = pos;
        let matched = candidate != usize::MAX
            && pos - candidate <= MAX_DISTANCE
            && input[candidate..candidate + MIN_MATCH] == input[pos..pos + MIN_MATCH];
        if !matched {
            pos += 1;
            continue;
        }
        // Extend the match, but never into the trailing literal window.
        let max_end = len - LAST_LITERALS;
        let mut end = pos + MIN_MATCH;
        let mut probe = candidate + MIN_MATCH;
        while end < max_end && input[probe] == input[end] {
            end += 1;
            probe += 1;
        }
        write_sequence(&mut out, &input[anchor..pos], pos - candidate, end - pos);
        // Index the covered positions so a later match can start inside.
        for inner in pos + 1..end.min(match_limit) {
            table[hash4(input, inner)] = inner;
        }
        pos = end;
        anchor = end;
    }

    // Trailing literals (always at least `LAST_LITERALS` bytes).
    let tail = &input[anchor..];
    let token = (tail.len().min(15) as u8) << 4;
    out.push(token);
    if tail.len() >= 15 {
        write_length(&mut out, tail.len() - 15);
    }
    out.extend_from_slice(tail);
    out
}

/// Reads one length extension from `input` starting at `*pos`.
fn read_length(
    input: &[u8],
    pos: &mut usize,
) -> Result<usize, String> {
    let mut length = 0usize;
    loop {
        let byte = *input
            .get(*pos)
            .ok_or_else(|| "lz4 length extension is truncated".to_string())?;
        *pos += 1;
        length = length
            .checked_add(usize::from(byte))
            .ok_or_else(|| "lz4 length overflows".to_string())?;
        if byte != 255 {
            return Ok(length);
        }
    }
}

/// Decompresses a Hadoop-framed LZ4 page (`CODEC_LZ4 = 5`):
/// `[u32 BE uncompressed_len][u32 BE compressed_len][raw_lz4_block]...`
pub(super) fn decompress_hadoop(
    input: &[u8],
    expected: usize,
) -> Result<Vec<u8>, String> {
    if expected == 0 {
        return Ok(Vec::new());
    }
    let mut pos = 0_usize;
    let mut out = Vec::with_capacity(expected);
    while out.len() < expected {
        let hdr = input
            .get(pos..pos + 8)
            .ok_or_else(|| "Hadoop LZ4 block header is truncated".to_owned())?;
        let u_len = u32::from_be_bytes([hdr[0], hdr[1], hdr[2], hdr[3]]) as usize;
        let c_len = u32::from_be_bytes([hdr[4], hdr[5], hdr[6], hdr[7]]) as usize;
        pos += 8;
        let payload = input
            .get(pos..pos + c_len)
            .ok_or_else(|| "Hadoop LZ4 block payload is truncated".to_owned())?;
        pos += c_len;
        if c_len == u_len {
            out.extend_from_slice(payload);
        } else {
            let dec = decompress(payload, u_len)?;
            out.extend_from_slice(&dec);
        }
    }
    if out.len() != expected {
        return Err(format!(
            "Hadoop LZ4 decompressed {} bytes, expected {expected}",
            out.len()
        ));
    }
    Ok(out)
}

/// Decompresses one LZ4 block into exactly `expected` bytes.
pub fn decompress(
    input: &[u8],
    expected: usize,
) -> Result<Vec<u8>, String> {
    let mut out: Vec<u8> = Vec::with_capacity(expected.min(1 << 26));
    let mut pos = 0usize;
    while pos < input.len() {
        let token = input[pos];
        pos += 1;
        let mut literals = usize::from(token >> 4);
        if literals == 15 {
            literals = literals
                .checked_add(read_length(input, &mut pos)?)
                .ok_or_else(|| "lz4 literal length overflows".to_string())?;
        }
        let end = pos
            .checked_add(literals)
            .ok_or_else(|| "lz4 literal run overflows".to_string())?;
        let chunk = input
            .get(pos..end)
            .ok_or_else(|| "lz4 literal run is truncated".to_string())?;
        out.extend_from_slice(chunk);
        pos = end;
        // The last sequence has no match part.
        if pos == input.len() {
            break;
        }
        let offset_bytes = input
            .get(pos..pos + 2)
            .ok_or_else(|| "lz4 match offset is truncated".to_string())?;
        let offset = usize::from(u16::from_le_bytes([offset_bytes[0], offset_bytes[1]]));
        pos += 2;
        if offset == 0 || offset > out.len() {
            return Err("lz4 match offset is out of bounds".to_string());
        }
        let mut match_len = usize::from(token & 0x0f);
        if match_len == 15 {
            match_len = match_len
                .checked_add(read_length(input, &mut pos)?)
                .ok_or_else(|| "lz4 match length overflows".to_string())?;
        }
        let match_len = match_len
            .checked_add(MIN_MATCH)
            .ok_or_else(|| "lz4 match length overflows".to_string())?;
        let start = out.len() - offset;
        // Overlapping copies are byte-by-byte by definition.
        for index in 0..match_len {
            let byte = out[start + index];
            out.push(byte);
        }
        if out.len() > expected {
            return Err(format!(
                "lz4 block decoded to more than the declared {expected} bytes"
            ));
        }
    }
    if out.len() != expected {
        return Err(format!(
            "lz4 block decoded to {} bytes, expected {expected}",
            out.len()
        ));
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::compress;
    use super::decompress;

    #[test]
    fn round_trip_on_repetitive_and_random_data() {
        let mut data = Vec::new();
        for index in 0..20_000u32 {
            data.extend_from_slice(format!("record-{}-{};", index % 97, index).as_bytes());
        }
        let packed = compress(&data);
        assert!(packed.len() < data.len());
        assert_eq!(decompress(&packed, data.len()).expect("decode"), data);

        // A pseudo-random stream with no long matches must still round-trip.
        let mut noise = Vec::with_capacity(5000);
        let mut state = 0x1234_5678u32;
        for _ in 0..5000 {
            state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
            noise.push((state >> 24) as u8);
        }
        let packed = compress(&noise);
        assert_eq!(decompress(&packed, noise.len()).expect("decode"), noise);
    }

    #[test]
    fn short_and_empty_inputs_round_trip() {
        for data in [
            Vec::new(),
            vec![1u8],
            b"abc".to_vec(),
            b"0123456789ab".to_vec(),
        ] {
            let packed = compress(&data);
            assert_eq!(decompress(&packed, data.len()).expect("decode"), data);
        }
    }

    #[test]
    fn corrupt_blocks_are_rejected() {
        let data = vec![7u8; 100];
        let packed = compress(&data);
        assert!(decompress(&packed, 50).is_err());
        let mut bad = packed;
        bad.truncate(bad.len() / 2);
        assert!(decompress(&bad, data.len()).is_err());
    }
}
