//! Minimal zlib/DEFLATE decompression for the HDF5 `deflate` filter.
//!
//! Only decompression is provided: stored, fixed-Huffman and dynamic-Huffman
//! blocks are supported (RFC 1950 and RFC 1951). Preset dictionaries, gzip
//! headers and multi-member streams are rejected with an error.

/// Upper bound on a single decompressed buffer (64 MiB).
const MAX_OUTPUT: usize = 64 * 1024 * 1024;

/// Base lengths for the length codes 257..=285.
const LENGTH_BASE: [u16; 29] = [
    3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59, 67, 83, 99, 115, 131,
    163, 195, 227, 258,
];

/// Extra bits for the length codes 257..=285.
const LENGTH_EXTRA: [u8; 29] = [
    0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0,
];

/// Base distances for the distance codes 0..=29.
const DIST_BASE: [u16; 30] = [
    1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513, 769, 1025, 1537,
    2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577,
];

/// Extra bits for the distance codes 0..=29.
const DIST_EXTRA: [u8; 30] = [
    0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13,
    13,
];

/// Order in which the code-length code lengths are stored in a dynamic block.
const CODE_LENGTH_ORDER: [usize; 19] = [
    16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15,
];

/// A least-significant-bit-first bit reader over a byte slice.
struct BitReader<'a> {
    /// Source bytes.
    data: &'a [u8],
    /// Index of the next unread byte.
    pos: usize,
    /// Buffered bits, least significant bit first.
    bits: u32,
    /// Number of valid bits in `bits`.
    nbits: u32,
}

impl<'a> BitReader<'a> {
    /// Creates a reader positioned at the start of `data`.
    const fn new(data: &'a [u8]) -> Self {
        Self {
            data,
            pos: 0,
            bits: 0,
            nbits: 0,
        }
    }

    /// Reads `n` bits (`n <= 16`) and returns them right-aligned.
    fn take(
        &mut self,
        n: u32,
    ) -> Result<u32, String> {
        while self.nbits < n {
            let byte = *self
                .data
                .get(self.pos)
                .ok_or_else(|| "unexpected end of deflate stream".to_owned())?;
            self.pos += 1;
            self.bits |= u32::from(byte) << self.nbits;
            self.nbits += 8;
        }
        let mask = if n >= 32 {
            u32::MAX
        } else {
            (1_u32 << n) - 1
        };
        let value = self.bits & mask;
        self.bits >>= n;
        self.nbits -= n;
        Ok(value)
    }

    /// Discards bits up to the next byte boundary.
    fn align(&mut self) -> Result<(), String> {
        let drop = self.nbits % 8;
        let _ = self.take(drop)?;
        Ok(())
    }

    /// Reads one byte, honouring bits already buffered.
    fn byte(&mut self) -> Result<u8, String> {
        #[allow(clippy::cast_possible_truncation)]
        let value = self.take(8)? as u8;
        Ok(value)
    }
}

/// A canonical Huffman decoder built from a table of code lengths.
struct Huffman {
    /// Number of codes of each length (index 0 is unused).
    counts: [u16; 16],
    /// Symbols ordered by code.
    symbols: Vec<u16>,
}

impl Huffman {
    /// Builds a canonical decoder; an over-subscribed code is an error.
    fn new(lengths: &[u8]) -> Result<Self, String> {
        let mut counts = [0_u16; 16];
        for &len in lengths {
            if len > 15 {
                return Err("invalid Huffman code length".to_owned());
            }
            counts[usize::from(len)] += 1;
        }
        // Counts for length 0 are symbols that are never used.
        let mut left: i32 = 1;
        for len in 1..16 {
            left <<= 1;
            left -= i32::from(counts[len]);
            if left < 0 {
                return Err("over-subscribed Huffman code".to_owned());
            }
        }
        let mut offsets = [0_u16; 16];
        for len in 1..15 {
            offsets[len + 1] = offsets[len] + counts[len];
        }
        let mut symbols = vec![0_u16; lengths.len()];
        for (sym, &len) in lengths.iter().enumerate() {
            if len != 0 {
                #[allow(clippy::cast_possible_truncation)]
                let sym16 = sym as u16;
                let slot = usize::from(offsets[usize::from(len)]);
                if let Some(cell) = symbols.get_mut(slot) {
                    *cell = sym16;
                }
                offsets[usize::from(len)] += 1;
            }
        }
        Ok(Self { counts, symbols })
    }

    /// Decodes one symbol from `br`.
    fn decode(
        &self,
        br: &mut BitReader<'_>,
    ) -> Result<u16, String> {
        let mut code: i32 = 0;
        let mut first: i32 = 0;
        let mut index: i32 = 0;
        for len in 1..=15 {
            #[allow(clippy::cast_possible_wrap)]
            let bit = br.take(1)? as i32;
            code |= bit;
            let count = i32::from(self.counts[len]);
            if code - first < count {
                #[allow(clippy::cast_possible_wrap)]
                let at = (index + (code - first)) as usize;
                return self
                    .symbols
                    .get(at)
                    .copied()
                    .ok_or_else(|| "corrupt Huffman table".to_owned());
            }
            index += count;
            first = (first + count) << 1;
            code <<= 1;
        }
        Err("invalid Huffman code".to_owned())
    }
}

/// Builds the fixed literal/length and distance tables of RFC 1951.
fn fixed_tables() -> Result<(Huffman, Huffman), String> {
    let mut lit = [0_u8; 288];
    for (sym, slot) in lit.iter_mut().enumerate() {
        *slot = match sym {
            | 0..=143 => 8,
            | 144..=255 => 9,
            | 256..=279 => 7,
            | _ => 8,
        };
    }
    let dist = [5_u8; 30];
    Ok((Huffman::new(&lit)?, Huffman::new(&dist)?))
}

/// Reads a dynamic-block header and returns the two Huffman tables.
fn dynamic_tables(br: &mut BitReader<'_>) -> Result<(Huffman, Huffman), String> {
    #[allow(clippy::cast_possible_truncation)]
    let nlit = br.take(5)? as usize + 257;
    #[allow(clippy::cast_possible_truncation)]
    let ndist = br.take(5)? as usize + 1;
    #[allow(clippy::cast_possible_truncation)]
    let ncode = br.take(4)? as usize + 4;
    if nlit > 286 || ndist > 30 {
        return Err("invalid dynamic Huffman header".to_owned());
    }
    let mut code_lengths = [0_u8; 19];
    for &slot in CODE_LENGTH_ORDER.iter().take(ncode) {
        #[allow(clippy::cast_possible_truncation)]
        let value = br.take(3)? as u8;
        if let Some(cell) = code_lengths.get_mut(slot) {
            *cell = value;
        }
    }
    let code_table = Huffman::new(&code_lengths)?;
    let total = nlit + ndist;
    let mut lengths = Vec::with_capacity(total);
    while lengths.len() < total {
        let sym = code_table.decode(br)?;
        match sym {
            | 0..=15 => {
                #[allow(clippy::cast_possible_truncation)]
                lengths.push(sym as u8);
            },
            | 16 => {
                let prev = lengths
                    .last()
                    .copied()
                    .ok_or_else(|| "repeat code with no previous length".to_owned())?;
                let repeat = 3 + br.take(2)? as usize;
                for _ in 0..repeat {
                    lengths.push(prev);
                }
            },
            | 17 => {
                let repeat = 3 + br.take(3)? as usize;
                for _ in 0..repeat {
                    lengths.push(0);
                }
            },
            | 18 => {
                let repeat = 11 + br.take(7)? as usize;
                for _ in 0..repeat {
                    lengths.push(0);
                }
            },
            | _ => return Err("invalid code-length symbol".to_owned()),
        }
        if lengths.len() > total {
            return Err("code lengths overflow the dynamic header".to_owned());
        }
    }
    let (lit_lengths, dist_lengths) = lengths.split_at(nlit);
    Ok((Huffman::new(lit_lengths)?, Huffman::new(dist_lengths)?))
}

/// Copies a stored (uncompressed) block into `out`.
fn stored_block(
    br: &mut BitReader<'_>,
    out: &mut Vec<u8>,
) -> Result<(), String> {
    br.align()?;
    let lo = u16::from(br.byte()?);
    let hi = u16::from(br.byte()?);
    let len = usize::from(lo | (hi << 8));
    let nlo = u16::from(br.byte()?);
    let nhi = u16::from(br.byte()?);
    let nlen = nlo | (nhi << 8);
    if len != usize::from(!nlen) {
        return Err("stored block length check failed".to_owned());
    }
    if out.len() + len > MAX_OUTPUT {
        return Err("decompressed data exceeds the 64 MiB limit".to_owned());
    }
    for _ in 0..len {
        out.push(br.byte()?);
    }
    Ok(())
}

/// Decodes Huffman-coded symbols into `out` until the end-of-block code.
fn huffman_block(
    br: &mut BitReader<'_>,
    out: &mut Vec<u8>,
    lit: &Huffman,
    dist: &Huffman,
) -> Result<(), String> {
    loop {
        let sym = lit.decode(br)?;
        match sym {
            | 0..=255 => {
                if out.len() >= MAX_OUTPUT {
                    return Err("decompressed data exceeds the 64 MiB limit".to_owned());
                }
                #[allow(clippy::cast_possible_truncation)]
                out.push(sym as u8);
            },
            | 256 => return Ok(()),
            | 257..=285 => {
                #[allow(clippy::cast_possible_truncation)]
                let idx = usize::from(sym - 257);
                let base = usize::from(*LENGTH_BASE.get(idx).ok_or("bad length code")?);
                let extra = u32::from(*LENGTH_EXTRA.get(idx).ok_or("bad length code")?);
                let len = base + br.take(extra)? as usize;
                let dsym = dist.decode(br)?;
                #[allow(clippy::cast_possible_truncation)]
                let didx = usize::from(dsym);
                let dbase = usize::from(*DIST_BASE.get(didx).ok_or("bad distance code")?);
                let dextra = u32::from(*DIST_EXTRA.get(didx).ok_or("bad distance code")?);
                let distance = dbase + br.take(dextra)? as usize;
                if distance == 0 || distance > out.len() {
                    return Err("invalid deflate back-reference".to_owned());
                }
                if out.len() + len > MAX_OUTPUT {
                    return Err("decompressed data exceeds the 64 MiB limit".to_owned());
                }
                let start = out.len() - distance;
                for step in 0..len {
                    let byte = *out
                        .get(start + step)
                        .ok_or_else(|| "invalid deflate back-reference".to_owned())?;
                    out.push(byte);
                }
            },
            | _ => return Err("invalid literal/length symbol".to_owned()),
        }
    }
}

/// Decompresses a raw DEFLATE stream (no zlib wrapper).
pub(crate) fn inflate(data: &[u8]) -> Result<Vec<u8>, String> {
    let mut br = BitReader::new(data);
    let mut out = Vec::new();
    loop {
        let last = br.take(1)?;
        let btype = br.take(2)?;
        match btype {
            | 0 => stored_block(&mut br, &mut out)?,
            | 1 => {
                let (lit, dist) = fixed_tables()?;
                huffman_block(&mut br, &mut out, &lit, &dist)?;
            },
            | 2 => {
                let (lit, dist) = dynamic_tables(&mut br)?;
                huffman_block(&mut br, &mut out, &lit, &dist)?;
            },
            | _ => return Err("invalid deflate block type".to_owned()),
        }
        if last == 1 {
            break;
        }
    }
    Ok(out)
}

/// Computes the zlib Adler-32 checksum of `data`.
fn adler32(data: &[u8]) -> u32 {
    let mut s1: u32 = 1;
    let mut s2: u32 = 0;
    for &byte in data {
        s1 = (s1 + u32::from(byte)) % 65_521;
        s2 = (s2 + s1) % 65_521;
    }
    (s2 << 16) | s1
}

/// Decompresses a zlib-wrapped DEFLATE stream.
///
/// # Errors
/// Returns an error for a malformed header, a preset dictionary, a corrupt
/// deflate stream, a checksum mismatch, or output beyond 64 MiB.
pub(crate) fn zlib_decompress(data: &[u8]) -> Result<Vec<u8>, String> {
    if data.len() < 6 {
        return Err("zlib stream is too short".to_owned());
    }
    let cmf = *data.first().ok_or("empty zlib stream")?;
    let flg = *data.get(1).ok_or("empty zlib stream")?;
    if cmf & 0x0f != 8 {
        return Err("unsupported zlib compression method".to_owned());
    }
    if (u16::from(cmf) * 256 + u16::from(flg)) % 31 != 0 {
        return Err("zlib header checksum mismatch".to_owned());
    }
    if flg & 0x20 != 0 {
        return Err("zlib preset dictionaries are not supported".to_owned());
    }
    let body = data.get(2..).ok_or("truncated zlib stream")?;
    let out = inflate(body)?;
    let tail = data
        .get(data.len() - 4..)
        .ok_or("truncated zlib checksum")?;
    let expected = u32::from_be_bytes([
        *tail.first().ok_or("truncated zlib checksum")?,
        *tail.get(1).ok_or("truncated zlib checksum")?,
        *tail.get(2).ok_or("truncated zlib checksum")?,
        *tail.get(3).ok_or("truncated zlib checksum")?,
    ]);
    if adler32(&out) != expected {
        return Err("zlib Adler-32 checksum mismatch".to_owned());
    }
    Ok(out)
}

/// Compresses `data` into a zlib stream (RFC 1950 + RFC 1951 + Adler-32).
pub(crate) fn zlib_compress(data: &[u8]) -> Vec<u8> {
    let compressed = crate::io::parquet::gzip::deflate_compress(data);
    let mut out = Vec::with_capacity(compressed.len() + 6);
    out.push(0x78);
    out.push(0x9C);
    out.extend_from_slice(&compressed);
    out.extend_from_slice(&adler32(data).to_be_bytes());
    out
}

#[cfg(test)]
mod tests {
    use super::adler32;
    use super::inflate;
    use super::zlib_decompress;

    /// zlib stream produced by `zlib.compress(b"hello hello hello hello", 9)`.
    const HELLO: [u8; 16] = [
        120, 218, 203, 72, 205, 201, 201, 87, 200, 64, 39, 1, 104, 3, 8, 177,
    ];

    /// zlib stream produced by `zlib.compress` over the `DYNAMIC_TEXT` payload;
    /// it uses a dynamic-Huffman block.
    const DYNAMIC: [u8; 87] = [
        120, 218, 101, 143, 65, 10, 192, 32, 16, 3, 191, 226, 215, 118, 169, 212, 131, 66, 15, 253,
        63, 90, 67, 219, 129, 92, 66, 200, 198, 89, 55, 235, 29, 37, 31, 57, 99, 140, 40, 71, 237,
        203, 71, 191, 218, 171, 74, 52, 205, 175, 173, 148, 250, 207, 80, 32, 104, 199, 178, 162,
        209, 59, 207, 216, 124, 4, 162, 53, 248, 97, 66, 152, 24, 16, 247, 49, 70, 224, 235, 9,
        183, 187, 38, 30, 231, 123, 27,
    ];

    /// Expected plain text of [`DYNAMIC`].
    const DYNAMIC_TEXT: &str = "beta beta gamma delta alpha alpha delta gamma beta beta delta delta \
        delta beta beta beta delta alpha alpha beta alpha gamma alpha gamma delta delta delta delta \
        delta beta gamma alpha alpha beta delta beta gamma delta gamma delta delta gamma delta beta \
        gamma alpha gamma beta gamma alpha beta gamma gamma alpha alpha delta delta alpha gamma alpha";

    #[test]
    fn inflate_known_zlib_stream() {
        let out = zlib_decompress(&HELLO).expect("decompress");
        assert_eq!(out, b"hello hello hello hello");
    }

    #[test]
    fn inflate_dynamic_huffman_block() {
        let out = zlib_decompress(&DYNAMIC).expect("decompress");
        assert_eq!(out, DYNAMIC_TEXT.as_bytes());
    }

    #[test]
    fn inflate_stored_block() {
        // zlib header, final stored block holding "abc", then the checksum.
        let payload = b"abc";
        let mut stream = vec![0x78, 0x01, 0x01];
        #[allow(clippy::cast_possible_truncation)]
        let len = payload.len() as u16;
        stream.extend_from_slice(&len.to_le_bytes());
        stream.extend_from_slice(&(!len).to_le_bytes());
        stream.extend_from_slice(payload);
        stream.extend_from_slice(&adler32(payload).to_be_bytes());
        let out = zlib_decompress(&stream).expect("decompress");
        assert_eq!(out, payload);
    }

    #[test]
    fn inflate_rejects_bad_checksum() {
        let mut stream = HELLO;
        if let Some(last) = stream.last_mut() {
            *last ^= 0xff;
        }
        assert!(zlib_decompress(&stream).is_err());
    }

    #[test]
    fn inflate_rejects_truncated_stream() {
        for cut in [3, 8, 12] {
            let truncated = HELLO.get(..cut).expect("slice");
            assert!(zlib_decompress(truncated).is_err());
        }
    }

    #[test]
    fn inflate_rejects_invalid_block_type() {
        assert!(inflate(&[0x07]).is_err());
    }
}
