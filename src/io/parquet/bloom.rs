//! Parquet Split Block Bloom Filter (SBBF) with XXH64 hashing.
//!
//! Implements the exact Apache Parquet Bloom Filter specification:
//! - Hash: XXH64 (`seed = 0`)
//! - Algorithm: Split Block Bloom Filter (256-bit / 32-byte blocks, 8 x `u32` words)
//! - Compression: `UNCOMPRESSED`

/// Salt constants from the Parquet Bloom Filter specification.
const SALT: [u32; 8] = [
    0x47b6_137b,
    0x4497_4d91,
    0x8824_ad5b,
    0xa2b7_289d,
    0x7054_95c7,
    0x2df1_424b,
    0x9efc_4947,
    0x5c6b_fb31,
];

/// Minimum SBBF bitset size in bytes (one 32-byte block = 256 bits).
const MIN_BYTES: usize = 32;
/// Maximum SBBF bitset size in bytes (128 MiB).
const MAX_BYTES: usize = 128 * 1024 * 1024;

/// Standard XXH64 hash function (seed = `0` for Parquet).
#[must_use]
pub fn xxhash64(
    data: &[u8],
    seed: u64,
) -> u64 {
    const P1: u64 = 0x9E37_79B1_85EB_CA87;
    const P2: u64 = 0xC2B2_AE3D_27D4_EB4F;
    const P3: u64 = 0x1656_67B1_9E37_79F9;
    const P4: u64 = 0x85EB_CA77_C2B2_AE63;
    const P5: u64 = 0x27D4_EB2F_1656_67C5;

    let round = |acc: u64, lane: u64| -> u64 {
        acc.wrapping_add(lane.wrapping_mul(P2))
            .rotate_left(31)
            .wrapping_mul(P1)
    };
    let merge_round = |mut acc: u64, val: u64| -> u64 {
        acc ^= round(0, val);
        acc.wrapping_mul(P1).wrapping_add(P4)
    };

    let len = data.len();
    let mut pos = 0_usize;
    let mut h: u64;

    if len >= 32 {
        let mut v1 = seed.wrapping_add(P1).wrapping_add(P2);
        let mut v2 = seed.wrapping_add(P2);
        let mut v3 = seed;
        let mut v4 = seed.wrapping_sub(P1);

        while pos + 32 <= len {
            v1 = round(
                v1,
                u64::from_le_bytes(data[pos..pos + 8].try_into().unwrap()),
            );
            v2 = round(
                v2,
                u64::from_le_bytes(data[pos + 8..pos + 16].try_into().unwrap()),
            );
            v3 = round(
                v3,
                u64::from_le_bytes(data[pos + 16..pos + 24].try_into().unwrap()),
            );
            v4 = round(
                v4,
                u64::from_le_bytes(data[pos + 24..pos + 32].try_into().unwrap()),
            );
            pos += 32;
        }

        h = v1
            .rotate_left(1)
            .wrapping_add(v2.rotate_left(7))
            .wrapping_add(v3.rotate_left(12))
            .wrapping_add(v4.rotate_left(18));
        h = merge_round(h, v1);
        h = merge_round(h, v2);
        h = merge_round(h, v3);
        h = merge_round(h, v4);
    } else {
        h = seed.wrapping_add(P5);
    }

    h = h.wrapping_add(len as u64);

    while pos + 8 <= len {
        let k1 = round(
            0,
            u64::from_le_bytes(data[pos..pos + 8].try_into().unwrap()),
        );
        h ^= k1;
        h = h.rotate_left(27).wrapping_mul(P1).wrapping_add(P4);
        pos += 8;
    }

    if pos + 4 <= len {
        let k1 = u64::from(u32::from_le_bytes(data[pos..pos + 4].try_into().unwrap()));
        h ^= k1.wrapping_mul(P1);
        h = h.rotate_left(23).wrapping_mul(P2).wrapping_add(P3);
        pos += 4;
    }

    while pos < len {
        h ^= u64::from(data[pos]).wrapping_mul(P5);
        h = h.rotate_left(11).wrapping_mul(P1);
        pos += 1;
    }

    h ^= h >> 33;
    h = h.wrapping_mul(P2);
    h ^= h >> 29;
    h = h.wrapping_mul(P3);
    h ^= h >> 32;
    h
}

/// A Parquet Split Block Bloom Filter.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SplitBlockBloomFilter {
    /// Raw little-endian 32-byte blocks (`8 * u32` per block).
    blocks: Vec<[u32; 8]>,
}

impl SplitBlockBloomFilter {
    /// Creates an empty SBBF sized for `ndv` distinct values and false-positive
    /// probability `fpp` (clamped to `[32, 128 MiB]` and rounded up to a power of 2).
    #[must_use]
    pub fn with_ndv_fpp(
        ndv: usize,
        fpp: f64,
    ) -> Self {
        let fpp = fpp.clamp(1e-9, 0.5);
        let ndv = ndv.max(1) as f64;
        // Parquet spec formula: m = -ndv * ln(fpp) / (ln(2)^2) bits, or ~8*ndv/ln(2)
        let bits = (-ndv * fpp.ln() / (std::f64::consts::LN_2 * std::f64::consts::LN_2)).ceil();
        let bytes = ((bits / 8.0).ceil() as usize)
            .clamp(MIN_BYTES, MAX_BYTES)
            .next_power_of_two();
        Self {
            blocks: vec![[0_u32; 8]; bytes / 32],
        }
    }

    /// Reconstructs an SBBF from raw little-endian bytes.
    ///
    /// # Errors
    /// Fails if `bytes` length is not a power-of-two multiple of 32 in `[32, 128 MiB]`.
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, String> {
        if bytes.len() < MIN_BYTES || bytes.len() > MAX_BYTES || !bytes.len().is_power_of_two() {
            return Err(format!(
                "invalid SBBF byte length {}; expected power of two in {MIN_BYTES}..={MAX_BYTES}",
                bytes.len()
            ));
        }
        let mut blocks = Vec::with_capacity(bytes.len() / 32);
        for chunk in bytes.as_chunks::<32>().0 {
            let mut block = [0_u32; 8];
            for (i, slot) in block.iter_mut().enumerate() {
                let off = i * 4;
                *slot = u32::from_le_bytes([
                    chunk[off],
                    chunk[off + 1],
                    chunk[off + 2],
                    chunk[off + 3],
                ]);
            }
            blocks.push(block);
        }
        Ok(Self { blocks })
    }

    /// Serializes the SBBF bitset to raw little-endian bytes.
    #[must_use]
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(self.blocks.len() * 32);
        for block in &self.blocks {
            for &word in block {
                out.extend_from_slice(&word.to_le_bytes());
            }
        }
        out
    }

    /// Byte length of the bitset.
    #[must_use]
    pub fn num_bytes(&self) -> usize {
        self.blocks.len() * 32
    }

    fn block_index(
        &self,
        hash: u64,
    ) -> usize {
        let nblocks = self.blocks.len() as u64;
        (((hash >> 32).wrapping_mul(nblocks)) >> 32) as usize
    }

    fn mask(hash: u64) -> [u32; 8] {
        let key = hash as u32;
        let mut m = [0_u32; 8];
        for i in 0..8 {
            let bit = (key.wrapping_mul(SALT[i])) >> 27;
            m[i] = 1_u32 << bit;
        }
        m
    }

    /// Inserts a 64-bit XXH64 hash into the filter.
    pub fn insert_hash(
        &mut self,
        hash: u64,
    ) {
        let idx = self.block_index(hash);
        let mask = Self::mask(hash);
        for i in 0..8 {
            self.blocks[idx][i] |= mask[i];
        }
    }

    /// Checks whether a 64-bit XXH64 hash might be in the filter.
    #[must_use]
    pub fn check_hash(
        &self,
        hash: u64,
    ) -> bool {
        let idx = self.block_index(hash);
        let mask = Self::mask(hash);
        for i in 0..8 {
            if (self.blocks[idx][i] & mask[i]) == 0 {
                return false;
            }
        }
        true
    }

    /// Inserts a raw byte slice into the filter.
    pub fn insert_bytes(
        &mut self,
        val: &[u8],
    ) {
        self.insert_hash(xxhash64(val, 0));
    }

    /// Checks whether a raw byte slice might be in the filter.
    #[must_use]
    pub fn might_contain_bytes(
        &self,
        val: &[u8],
    ) -> bool {
        self.check_hash(xxhash64(val, 0))
    }

    /// Inserts a UTF-8 string into the filter.
    pub fn insert_str(
        &mut self,
        val: &str,
    ) {
        self.insert_bytes(val.as_bytes());
    }

    /// Checks whether a UTF-8 string might be in the filter.
    #[must_use]
    pub fn might_contain_str(
        &self,
        val: &str,
    ) -> bool {
        self.might_contain_bytes(val.as_bytes())
    }

    /// Inserts an `i32` value into the filter.
    pub fn insert_i32(
        &mut self,
        val: i32,
    ) {
        self.insert_bytes(&val.to_le_bytes());
    }

    /// Checks whether an `i32` value might be in the filter.
    #[must_use]
    pub fn might_contain_i32(
        &self,
        val: i32,
    ) -> bool {
        self.might_contain_bytes(&val.to_le_bytes())
    }

    /// Inserts an `i64` value into the filter.
    pub fn insert_i64(
        &mut self,
        val: i64,
    ) {
        self.insert_bytes(&val.to_le_bytes());
    }

    /// Checks whether an `i64` value might be in the filter.
    #[must_use]
    pub fn might_contain_i64(
        &self,
        val: i64,
    ) -> bool {
        self.might_contain_bytes(&val.to_le_bytes())
    }

    /// Inserts an `f32` value into the filter.
    pub fn insert_f32(
        &mut self,
        val: f32,
    ) {
        self.insert_bytes(&val.to_le_bytes());
    }

    /// Checks whether an `f32` value might be in the filter.
    #[must_use]
    pub fn might_contain_f32(
        &self,
        val: f32,
    ) -> bool {
        self.might_contain_bytes(&val.to_le_bytes())
    }

    /// Inserts an `f64` value into the filter.
    pub fn insert_f64(
        &mut self,
        val: f64,
    ) {
        self.insert_bytes(&val.to_le_bytes());
    }

    /// Checks whether an `f64` value might be in the filter.
    #[must_use]
    pub fn might_contain_f64(
        &self,
        val: f64,
    ) -> bool {
        self.might_contain_bytes(&val.to_le_bytes())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn xxhash64_known_vectors() {
        assert_eq!(xxhash64(b"", 0), 0xef46_db37_51d8_e999);
        assert_eq!(xxhash64(b"a", 0), 0xd24e_c4f1_a98c_6e5b);
        assert_eq!(xxhash64(b"abc", 0), 0x44bc_2cf5_ad77_0999);
    }

    #[test]
    fn sbbf_round_trip_and_membership() {
        let mut bf = SplitBlockBloomFilter::with_ndv_fpp(100, 0.01);
        for i in 0..100_i64 {
            bf.insert_i64(i);
        }
        let bytes = bf.to_bytes();
        let bf2 = SplitBlockBloomFilter::from_bytes(&bytes).unwrap();
        for i in 0..100_i64 {
            assert!(bf2.might_contain_i64(i));
        }
        let false_positives = (1000..2000_i64)
            .filter(|&i| bf2.might_contain_i64(i))
            .count();
        assert!(false_positives < 40, "too many FPs: {false_positives}");
    }
}
