//! The one place rssn serialises anything.
//!
//! Every byte rssn writes — the [kernel manifest](super::manifest), a
//! [`Value`](super::Value) crossing the C ABI, a session snapshot, a plugin
//! message — goes through this module, and this module uses **bincode-next**
//! exclusively:
//!
//! * [`encode_deterministic`] / [`decode_deterministic`] — bincode-next's
//!   stable binary format with deterministic map ordering. This is the
//!   format whose bytes are part of rssn's public contract: equal values
//!   always produce equal bytes, on every platform and in every process.
//! * [`encode_cbor`] / [`decode_cbor`] — bincode-next's CBOR implementation,
//!   also in deterministic mode, for hosts that want a self-describing
//!   format (browsers, JavaScript, anything that already speaks CBOR).
//!
//! Both are the same library, so there is exactly one encoder to trust and
//! one set of bug fixes to carry.

use serde::Serialize;
use serde::de::DeserializeOwned;

use bincode_next::config;

/// Why a value could not be encoded or decoded.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CodecError {
    /// What went wrong, in the encoder's own words.
    pub message: String,
}

impl std::fmt::Display for CodecError {
    fn fmt(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        f.write_str(&self.message)
    }
}

impl std::error::Error for CodecError {}

/// Encodes `value` with bincode-next's deterministic binary format.
///
/// # Errors
/// Fails when the value cannot be represented, which for the types in this
/// crate never happens; the result is a `Result` so that a host type with a
/// fallible `Serialize` cannot panic the process.
pub fn encode_deterministic<T: Serialize + ?Sized>(value: &T) -> Result<Vec<u8>, CodecError> {
    bincode_next::serde::encode_to_vec(value, config::standard().with_deterministic_bincode())
        .map_err(|e| CodecError {
            message: e.to_string(),
        })
}

/// Decodes a value written by [`encode_deterministic`].
///
/// # Errors
/// Fails on a truncated, corrupt or type-mismatched payload.
pub fn decode_deterministic<T: DeserializeOwned>(bytes: &[u8]) -> Result<T, CodecError> {
    bincode_next::serde::decode_from_slice(bytes, config::standard().with_deterministic_bincode())
        .map(|(value, _)| value)
        .map_err(|e| CodecError {
            message: e.to_string(),
        })
}

/// Encodes `value` as deterministic CBOR.
///
/// # Errors
/// See [`encode_deterministic`].
pub fn encode_cbor<T: Serialize + ?Sized>(value: &T) -> Result<Vec<u8>, CodecError> {
    bincode_next::serde::encode_to_vec(value, config::standard().with_deterministic_cbor())
        .map_err(|e| CodecError {
            message: e.to_string(),
        })
}

/// Decodes deterministic CBOR written by [`encode_cbor`].
///
/// # Errors
/// See [`decode_deterministic`].
pub fn decode_cbor<T: DeserializeOwned>(bytes: &[u8]) -> Result<T, CodecError> {
    bincode_next::serde::decode_from_slice(bytes, config::standard().with_deterministic_cbor())
        .map(|(value, _)| value)
        .map_err(|e| CodecError {
            message: e.to_string(),
        })
}

/// The BLAKE3 digest of a byte string.
#[must_use]
pub fn blake3_256(bytes: &[u8]) -> [u8; 32] {
    *blake3::hash(bytes).as_bytes()
}

/// The BLAKE3 digest of the deterministic encoding of `value`.
///
/// This is how every content address in rssn is formed: the manifest hash,
/// the signature digest, a cache key. It is stable across processes and
/// platforms because the encoding is.
#[must_use]
pub fn digest<T: Serialize + ?Sized>(value: &T) -> [u8; 32] {
    encode_deterministic(value).map_or([0_u8; 32], |bytes| blake3_256(&bytes))
}

/// The digest of `value`, truncated to 16 bytes, for places where 128 bits of
/// collision resistance are plenty (signature identity inside a manifest).
#[must_use]
pub fn short_digest<T: Serialize + ?Sized>(value: &T) -> [u8; 16] {
    let full = digest(value);
    let mut out = [0_u8; 16];
    out.copy_from_slice(full.get(..16).unwrap_or(&full));
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Serialize, serde::Deserialize, PartialEq, Debug)]
    struct Sample {
        a: Vec<u32>,
        b: String,
        c: Option<i64>,
    }

    fn sample() -> Sample {
        Sample {
            a: vec![1, 2, 3],
            b: "rssn".to_owned(),
            c: Some(-7),
        }
    }

    #[test]
    fn round_trip_both_formats() {
        let value = sample();
        let bin = encode_deterministic(&value).unwrap_or_else(|e| panic!("{e}"));
        let back: Sample = decode_deterministic(&bin).unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(back, value);
        let cbor = encode_cbor(&value).unwrap_or_else(|e| panic!("{e}"));
        let back: Sample = decode_cbor(&cbor).unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(back, value);
    }

    #[test]
    fn encoding_is_deterministic() {
        let first = encode_deterministic(&sample()).unwrap_or_default();
        for _ in 0..8 {
            assert_eq!(encode_deterministic(&sample()).unwrap_or_default(), first);
        }
        assert_eq!(digest(&sample()), digest(&sample()));
    }

    #[test]
    fn corrupt_input_is_an_error_not_a_panic() {
        assert!(decode_deterministic::<Sample>(&[0xff, 0x00]).is_err());
        assert!(decode_cbor::<Sample>(&[0xff]).is_err());
    }

    #[test]
    fn blake3_matches_the_reference_vector() {
        // BLAKE3 of the empty string.
        assert_eq!(
            blake3_256(b""),
            [
                0xaf, 0x13, 0x49, 0xb9, 0xf5, 0xf9, 0xa1, 0xa6, 0xa0, 0x40, 0x4d, 0xea, 0x36, 0xdc,
                0xc9, 0x49, 0x9b, 0xcb, 0x25, 0xc9, 0xad, 0xc1, 0x12, 0xb7, 0xcc, 0x9a, 0x93, 0xca,
                0xe4, 0x1f, 0x32, 0x62,
            ]
        );
    }
}
