//! Gzip (RFC 1952) page compression and decompression.
//!
//! Decompression reuses the DEFLATE decoder written for the HDF5 `deflate`
//! filter, so every DEFLATE block type is accepted. Compression emits a single
//! fixed-Huffman block over an LZ77 parse with a 32 KiB window; that keeps the
//! encoder small while still producing ordinary, universally readable `.gz`
//! streams (the container, CRC-32 and length trailer are exact).

use crate::io::hdf5::inflate::inflate;

/// Maximum match length in a DEFLATE block.
const MAX_MATCH: usize = 258;
/// Minimum match length worth encoding.
const MIN_MATCH: usize = 3;
/// Sliding-window size.
const WINDOW: usize = 32_768;
/// Hash-table bit width for match finding.
const HASH_BITS: u32 = 15;
/// Hash-table size.
const HASH_SIZE: usize = 1 << HASH_BITS;

/// Length codes 257..=285: base lengths.
const LEN_BASE: [u16; 29] = [
    3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59, 67, 83, 99, 115, 131,
    163, 195, 227, 258,
];

/// Extra bits for the length codes.
const LEN_EXTRA: [u8; 29] = [
    0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0,
];

/// Distance codes 0..=29: base distances.
const DIST_BASE: [u16; 30] = [
    1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513, 769, 1025, 1537,
    2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577,
];

/// Extra bits for the distance codes.
const DIST_EXTRA: [u8; 30] = [
    0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13,
    13,
];

/// LSB-first bit writer.
struct BitWriter {
    out: Vec<u8>,
    accumulator: u64,
    bits: u32,
}

impl BitWriter {
    fn new() -> Self {
        Self {
            out: Vec::new(),
            accumulator: 0,
            bits: 0,
        }
    }

    /// Writes the low `count` bits of `value`, least significant first.
    fn put(
        &mut self,
        value: u32,
        count: u32,
    ) {
        if count == 0 {
            return;
        }
        self.accumulator |= u64::from(value & mask(count)) << self.bits;
        self.bits += count;
        while self.bits >= 8 {
            self.out.push((self.accumulator & 0xff) as u8);
            self.accumulator >>= 8;
            self.bits -= 8;
        }
    }

    /// Writes a Huffman code, whose bits are transmitted most significant
    /// first inside the LSB-first stream.
    fn put_code(
        &mut self,
        code: u32,
        count: u32,
    ) {
        let mut reversed = 0u32;
        for bit in 0..count {
            reversed |= ((code >> bit) & 1) << (count - 1 - bit);
        }
        self.put(reversed, count);
    }

    fn finish(&mut self) {
        if self.bits > 0 {
            self.out.push((self.accumulator & 0xff) as u8);
            self.accumulator = 0;
            self.bits = 0;
        }
    }
}

/// `count` one-bits, for `count <= 32`.
fn mask(count: u32) -> u32 {
    if count >= 32 {
        u32::MAX
    } else {
        (1u32 << count) - 1
    }
}

/// Fixed-Huffman code and bit length for a literal/length symbol.
fn literal_code(symbol: u32) -> (u32, u32) {
    match symbol {
        | 0..=143 => (0x30 + symbol, 8),
        | 144..=255 => (0x190 + symbol - 144, 9),
        | 256..=279 => (symbol - 256, 7),
        | _ => (0xc0 + symbol - 280, 8),
    }
}

/// Length code (257..=285) and extra bits for a match length.
fn length_code(length: usize) -> (u32, u32, u32) {
    let mut index = LEN_BASE.len() - 1;
    for (i, &base) in LEN_BASE.iter().enumerate() {
        if length < usize::from(base) {
            continue;
        }
        if i + 1 < LEN_BASE.len() && length >= usize::from(LEN_BASE[i + 1]) {
            continue;
        }
        index = i;
        break;
    }
    let base = u32::from(LEN_BASE[index]);
    let extra_bits = u32::from(LEN_EXTRA[index]);
    let extra_value = (length as u32).saturating_sub(base) & mask(extra_bits);
    (257 + index as u32, extra_value, extra_bits)
}

/// Distance code and extra bits for a match distance.
fn distance_code(distance: usize) -> (u32, u32, u32) {
    let mut index = DIST_BASE.len() - 1;
    for (i, &base) in DIST_BASE.iter().enumerate() {
        if distance < usize::from(base) {
            continue;
        }
        if i + 1 < DIST_BASE.len() && distance >= usize::from(DIST_BASE[i + 1]) {
            continue;
        }
        index = i;
        break;
    }
    let base = u32::from(DIST_BASE[index]);
    let extra_bits = u32::from(DIST_EXTRA[index]);
    let extra_value = (distance as u32).saturating_sub(base) & mask(extra_bits);
    (index as u32, extra_value, extra_bits)
}

/// Hashes three bytes into a table slot.
fn hash3(
    data: &[u8],
    pos: usize,
) -> usize {
    let v =
        u32::from(data[pos]) | (u32::from(data[pos + 1]) << 8) | (u32::from(data[pos + 2]) << 16);
    (v.wrapping_mul(0x9e37_79b1) >> (32 - HASH_BITS)) as usize
}

/// Emits one literal or match.
fn emit_literal(
    writer: &mut BitWriter,
    byte: u8,
) {
    let (code, bits) = literal_code(u32::from(byte));
    writer.put_code(code, bits);
}

/// Emits a back-reference.
fn emit_match(
    writer: &mut BitWriter,
    length: usize,
    distance: usize,
) {
    let (symbol, extra, extra_bits) = length_code(length);
    let (code, bits) = literal_code(symbol);
    writer.put_code(code, bits);
    writer.put(extra, extra_bits);
    let (dist_code, dist_extra, dist_bits) = distance_code(distance);
    writer.put_code(dist_code, 5);
    writer.put(dist_extra, dist_bits);
}

/// LZ77 parse over `data`, emitting literals and matches.
fn deflate_block(
    writer: &mut BitWriter,
    data: &[u8],
) {
    let n = data.len();
    let mut table = vec![usize::MAX; HASH_SIZE];
    let mut pos = 0usize;
    while pos < n {
        let mut best_len = 0usize;
        let mut best_dist = 0usize;
        if pos + MIN_MATCH <= n {
            let slot = hash3(data, pos);
            let candidate = table[slot];
            table[slot] = pos;
            if candidate != usize::MAX && pos - candidate <= WINDOW {
                let limit = MAX_MATCH.min(n - pos);
                let mut len = 0usize;
                while len < limit && data[candidate + len] == data[pos + len] {
                    len += 1;
                }
                if len >= MIN_MATCH {
                    best_len = len;
                    best_dist = pos - candidate;
                }
            }
        }
        if best_len >= MIN_MATCH {
            emit_match(writer, best_len, best_dist);
            // Index the positions covered by the match so later matches can
            // start inside it.
            let end = (pos + best_len).min(n.saturating_sub(MIN_MATCH - 1));
            for p in pos + 1..end {
                table[hash3(data, p)] = p;
            }
            pos += best_len;
        } else {
            emit_literal(writer, data[pos]);
            pos += 1;
        }
    }
    // End of block: symbol 256 has the seven-bit code 0.
    writer.put_code(0, 7);
}

/// CRC-32 (IEEE 802.3), as required by the gzip trailer.
fn crc32(data: &[u8]) -> u32 {
    let mut crc = 0xffff_ffffu32;
    for &byte in data {
        crc ^= u32::from(byte);
        for _ in 0..8 {
            let odd = crc & 1;
            crc >>= 1;
            if odd != 0 {
                crc ^= 0xedb8_8320;
            }
        }
    }
    !crc
}

/// Compresses `data` into a raw DEFLATE stream (single fixed-Huffman block).
pub(crate) fn deflate_compress(data: &[u8]) -> Vec<u8> {
    let mut writer = BitWriter::new();
    writer.put(1, 1); // BFINAL
    writer.put(1, 2); // BTYPE = 01, fixed Huffman
    deflate_block(&mut writer, data);
    writer.finish();
    writer.out
}

/// Compresses `data` into a gzip member.
pub(super) fn compress(data: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(data.len() / 2 + 32);
    // Magic, deflate, no flags, mtime 0, no extra flags, OS unknown.
    out.extend_from_slice(&[0x1f, 0x8b, 0x08, 0x00, 0, 0, 0, 0, 0x00, 0xff]);
    out.extend_from_slice(&deflate_compress(data));
    out.extend_from_slice(&crc32(data).to_le_bytes());
    out.extend_from_slice(&u32::try_from(data.len()).unwrap_or(u32::MAX).to_le_bytes());
    out
}

/// Decompresses one gzip member.
pub(super) fn decompress(data: &[u8]) -> Result<Vec<u8>, String> {
    if data.len() < 18 {
        return Err("gzip page is truncated".to_string());
    }
    if data[0] != 0x1f || data[1] != 0x8b {
        return Err("gzip page has the wrong magic".to_string());
    }
    if data[2] != 8 {
        return Err(format!(
            "gzip page uses compression method {}, only deflate is supported",
            data[2]
        ));
    }
    let flags = data[3];
    let mut pos = 10usize;
    if flags & 0x04 != 0 {
        let extra = data
            .get(pos..pos + 2)
            .ok_or_else(|| "gzip extra field is truncated".to_string())?;
        let len = usize::from(u16::from_le_bytes([extra[0], extra[1]]));
        pos += 2 + len;
    }
    if flags & 0x08 != 0 {
        pos = skip_zero_terminated(data, pos, "file name")?;
    }
    if flags & 0x10 != 0 {
        pos = skip_zero_terminated(data, pos, "comment")?;
    }
    if flags & 0x02 != 0 {
        pos += 2;
    }
    let body_end = data
        .len()
        .checked_sub(8)
        .ok_or_else(|| "gzip trailer is missing".to_string())?;
    if pos > body_end {
        return Err("gzip header overruns the page".to_string());
    }
    let body = &data[pos..body_end];
    let out = inflate(body)?;
    let trailer = &data[body_end..];
    let expected_crc = u32::from_le_bytes([trailer[0], trailer[1], trailer[2], trailer[3]]);
    let expected_len = u32::from_le_bytes([trailer[4], trailer[5], trailer[6], trailer[7]]);
    if crc32(&out) != expected_crc {
        return Err("gzip page fails its CRC-32 check".to_string());
    }
    if out.len() as u32 != expected_len {
        return Err(format!(
            "gzip page holds {} bytes but declares {expected_len}",
            out.len()
        ));
    }
    Ok(out)
}

/// Advances past a zero-terminated header string.
fn skip_zero_terminated(
    data: &[u8],
    mut pos: usize,
    field: &str,
) -> Result<usize, String> {
    while pos < data.len() && data[pos] != 0 {
        pos += 1;
    }
    if pos >= data.len() {
        return Err(format!("gzip {field} is not terminated"));
    }
    Ok(pos + 1)
}

#[cfg(test)]
mod tests {
    use super::compress;
    use super::decompress;

    #[test]
    fn round_trip_on_repetitive_and_noisy_data() {
        let mut data = Vec::new();
        for i in 0..20_000u32 {
            data.extend_from_slice(format!("row-{}-{}", i % 97, i).as_bytes());
            data.push(if i % 3 == 0 { b'\n' } else { b',' });
        }
        let packed = compress(&data);
        assert_eq!(&packed[..2], &[0x1f, 0x8b]);
        assert!(
            packed.len() < data.len(),
            "compressed {} bytes into {}",
            data.len(),
            packed.len()
        );
        let back = decompress(&packed).expect("decompress");
        assert_eq!(back, data);
    }

    #[test]
    fn empty_and_single_byte_inputs_round_trip() {
        for data in [Vec::new(), vec![7u8], vec![0u8; 1000]] {
            let packed = compress(&data);
            assert_eq!(decompress(&packed).expect("decompress"), data);
        }
    }

    #[test]
    fn corrupt_streams_are_rejected() {
        let packed = compress(b"the quick brown fox jumps over the lazy dog");
        assert!(decompress(&packed[..10]).is_err());
        let mut bad_magic = packed.clone();
        bad_magic[0] = 0;
        assert!(decompress(&bad_magic).is_err());
        let mut bad_crc = packed.clone();
        let len = bad_crc.len();
        bad_crc[len - 8] ^= 0xff;
        assert!(decompress(&bad_crc).is_err());
    }
}
