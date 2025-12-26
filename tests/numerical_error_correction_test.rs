//! Tests for numerical error correction module.
//!
//! This module contains comprehensive unit tests and property-based tests
//! for Reed-Solomon codes, Hamming codes, CRC checksums, and other
//! error correction utilities.

use rssn::numerical::error_correction::*;

// ============================================================================
// Reed-Solomon Tests
// ============================================================================

#[test]

fn test_reed_solomon_encode_basic() {

    let message = vec![
        0x01, 0x02, 0x03, 0x04,
    ];

    let codeword = reed_solomon_encode(&message, 4).unwrap();

    assert_eq!(codeword.len(), 8); // 4 data + 4 parity
                                   // First 4 bytes should be the original message
    assert_eq!(&codeword[..4], &message);
}

#[test]

fn test_reed_solomon_encode_empty_message() {

    let message: Vec<u8> = vec![];

    let codeword = reed_solomon_encode(&message, 4).unwrap();

    assert_eq!(codeword.len(), 4); // 0 data + 4 parity
}

#[test]

fn test_reed_solomon_encode_max_length() {

    // Maximum message length is 255 - n_parity
    let message: Vec<u8> = (0..251).collect();

    let codeword = reed_solomon_encode(&message, 4).unwrap();

    assert_eq!(codeword.len(), 255);
}

#[test]

fn test_reed_solomon_encode_too_long() {

    let message: Vec<u8> = (0..252).collect();

    let result = reed_solomon_encode(&message, 4);

    assert!(result.is_err());
}

#[test]

fn test_reed_solomon_check_valid() {

    let message = vec![
        0x01, 0x02, 0x03, 0x04,
    ];

    let codeword = reed_solomon_encode(&message, 4).unwrap();

    assert!(reed_solomon_check(&codeword, 4));
}

#[test]

fn test_reed_solomon_check_corrupted() {

    let message = vec![
        0x01, 0x02, 0x03, 0x04,
    ];

    let mut codeword = reed_solomon_encode(&message, 4).unwrap();

    codeword[0] ^= 0xFF; // Corrupt first byte
    assert!(!reed_solomon_check(&codeword, 4));
}

#[test]

fn test_reed_solomon_decode_no_errors() {

    let message = vec![
        0x01, 0x02, 0x03, 0x04,
    ];

    let mut codeword = reed_solomon_encode(&message, 4).unwrap();

    reed_solomon_decode(&mut codeword, 4).unwrap();

    assert_eq!(&codeword[..4], &message);
}

// ============================================================================
// Hamming Code Tests
// ============================================================================

#[test]

fn test_hamming_encode_basic() {

    let data = vec![1, 0, 1, 1];

    let codeword = hamming_encode_numerical(&data).unwrap();

    assert_eq!(codeword.len(), 7);
}

#[test]

fn test_hamming_encode_all_zeros() {

    let data = vec![0, 0, 0, 0];

    let codeword = hamming_encode_numerical(&data).unwrap();

    assert_eq!(codeword, vec![0, 0, 0, 0, 0, 0, 0]);
}

#[test]

fn test_hamming_encode_all_ones() {

    let data = vec![1, 1, 1, 1];

    let codeword = hamming_encode_numerical(&data).unwrap();

    assert_eq!(codeword.len(), 7);
}

#[test]

fn test_hamming_encode_wrong_length() {

    let data = vec![1, 0, 1]; // Only 3 bits
    let result = hamming_encode_numerical(&data);

    assert!(result.is_none());
}

#[test]

fn test_hamming_decode_no_error() {

    let data = vec![1, 0, 1, 1];

    let codeword = hamming_encode_numerical(&data).unwrap();

    let (decoded, error_pos) = hamming_decode_numerical(&codeword).unwrap();

    assert_eq!(decoded, data);

    assert_eq!(error_pos, None);
}

#[test]

fn test_hamming_decode_single_error() {

    let data = vec![1, 0, 1, 1];

    let mut codeword = hamming_encode_numerical(&data).unwrap();

    codeword[2] ^= 1; // Introduce error at position 3 (1-indexed)
    let (decoded, error_pos) = hamming_decode_numerical(&codeword).unwrap();

    assert_eq!(decoded, data);

    assert_eq!(error_pos, Some(3));
}

#[test]

fn test_hamming_decode_parity_error() {

    let data = vec![1, 0, 1, 1];

    let mut codeword = hamming_encode_numerical(&data).unwrap();

    codeword[0] ^= 1; // Error in parity bit
    let (decoded, error_pos) = hamming_decode_numerical(&codeword).unwrap();

    assert_eq!(decoded, data);

    assert_eq!(error_pos, Some(1));
}

#[test]

fn test_hamming_decode_wrong_length() {

    let codeword = vec![1, 0, 1, 1, 0, 1]; // Only 6 bits
    let result = hamming_decode_numerical(&codeword);

    assert!(result.is_err());
}

#[test]

fn test_hamming_check_valid() {

    let data = vec![1, 0, 1, 1];

    let codeword = hamming_encode_numerical(&data).unwrap();

    assert!(hamming_check_numerical(&codeword));
}

#[test]

fn test_hamming_check_invalid() {

    let data = vec![1, 0, 1, 1];

    let mut codeword = hamming_encode_numerical(&data).unwrap();

    codeword[2] ^= 1; // Introduce error
    assert!(!hamming_check_numerical(&codeword));
}

#[test]

fn test_hamming_distance_equal() {

    let a = vec![1, 0, 1, 1];

    let b = vec![1, 0, 1, 1];

    assert_eq!(hamming_distance_numerical(&a, &b), Some(0));
}

#[test]

fn test_hamming_distance_different() {

    let a = vec![1, 0, 1, 1];

    let b = vec![0, 0, 1, 0];

    assert_eq!(hamming_distance_numerical(&a, &b), Some(2));
}

#[test]

fn test_hamming_distance_all_different() {

    let a = vec![0, 0, 0, 0];

    let b = vec![1, 1, 1, 1];

    assert_eq!(hamming_distance_numerical(&a, &b), Some(4));
}

#[test]

fn test_hamming_distance_length_mismatch() {

    let a = vec![1, 0, 1];

    let b = vec![1, 0, 1, 1];

    assert_eq!(hamming_distance_numerical(&a, &b), None);
}

#[test]

fn test_hamming_weight_all_zeros() {

    let data = vec![0, 0, 0, 0];

    assert_eq!(hamming_weight_numerical(&data), 0);
}

#[test]

fn test_hamming_weight_all_ones() {

    let data = vec![1, 1, 1, 1];

    assert_eq!(hamming_weight_numerical(&data), 4);
}

#[test]

fn test_hamming_weight_mixed() {

    let data = vec![1, 0, 1, 0, 1];

    assert_eq!(hamming_weight_numerical(&data), 3);
}

#[test]

fn test_hamming_weight_empty() {

    let data: Vec<u8> = vec![];

    assert_eq!(hamming_weight_numerical(&data), 0);
}

// ============================================================================
// BCH Code Tests
// ============================================================================

#[test]

fn test_bch_encode_basic() {

    let data = vec![1, 0, 1, 1, 0, 1];

    let codeword = bch_encode(&data, 2);

    assert!(codeword.len() > data.len());
}

#[test]

fn test_bch_roundtrip_no_errors() {

    let data = vec![1, 0, 1, 1];

    let codeword = bch_encode(&data, 2);

    let decoded = bch_decode(&codeword, 2).unwrap();

    assert_eq!(decoded, data);
}

// ============================================================================
// CRC-32 Tests
// ============================================================================

#[test]

fn test_crc32_empty() {

    let data: &[u8] = b"";

    let crc = crc32_compute_numerical(data);

    // CRC of empty string
    assert_eq!(crc, 0x00000000);
}

#[test]

fn test_crc32_hello_world() {

    let data = b"Hello, World!";

    let crc = crc32_compute_numerical(data);

    // Known CRC-32 value for "Hello, World!"
    assert_eq!(crc, 0xEC4AC3D0);
}

#[test]

fn test_crc32_verify_valid() {

    let data = b"Hello, World!";

    let crc = crc32_compute_numerical(data);

    assert!(crc32_verify_numerical(data, crc));
}

#[test]

fn test_crc32_verify_invalid() {

    let data = b"Hello, World!";

    let wrong_crc = 0x12345678;

    assert!(!crc32_verify_numerical(data, wrong_crc));
}

#[test]

fn test_crc32_streaming() {

    let data1 = b"Hello, ";

    let data2 = b"World!";

    // Incremental computation
    let crc = crc32_update_numerical(0xFFFFFFFF, data1);

    let crc = crc32_update_numerical(crc, data2);

    let crc = crc32_finalize_numerical(crc);

    // Full computation
    let full_crc = crc32_compute_numerical(b"Hello, World!");

    assert_eq!(crc, full_crc);
}

// ============================================================================
// CRC-16 Tests
// ============================================================================

#[test]

fn test_crc16_basic() {

    let data = b"123456789";

    let crc = crc16_compute(data);

    // CRC-16 (IBM/Modbus) value for "123456789"
    // Our implementation uses the reflected polynomial 0xA001
    assert_eq!(crc, 0x4B37);
}

#[test]

fn test_crc16_empty() {

    let data: &[u8] = b"";

    let crc = crc16_compute(data);

    // CRC-16 of empty string
    assert_eq!(crc, 0xFFFF);
}

// ============================================================================
// CRC-8 Tests
// ============================================================================

#[test]

fn test_crc8_basic() {

    let data = b"123456789";

    let crc = crc8_compute(data);

    // CRC-8 ITU value for "123456789"
    assert_ne!(crc, 0); // Just verify it produces a result
}

#[test]

fn test_crc8_empty() {

    let data: &[u8] = b"";

    let crc = crc8_compute(data);

    assert_eq!(crc, 0);
}

// ============================================================================
// Interleaving Tests
// ============================================================================

#[test]

fn test_interleave_basic() {

    let data = vec![1, 2, 3, 4, 5, 6];

    let depth = 3;

    let interleaved = interleave(&data, depth);

    assert_eq!(interleaved.len(), data.len());
}

#[test]

fn test_interleave_deinterleave_roundtrip() {

    let data = vec![
        1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
    ];

    let depth = 4;

    let interleaved = interleave(&data, depth);

    let deinterleaved = deinterleave(&interleaved, depth);

    assert_eq!(deinterleaved, data);
}

#[test]

fn test_interleave_depth_1() {

    let data = vec![1, 2, 3, 4];

    let interleaved = interleave(&data, 1);

    assert_eq!(interleaved, data);
}

#[test]

fn test_interleave_depth_0() {

    let data = vec![1, 2, 3, 4];

    let interleaved = interleave(&data, 0);

    assert_eq!(interleaved, data);
}

#[test]

fn test_interleave_empty() {

    let data: Vec<u8> = vec![];

    let interleaved = interleave(&data, 3);

    assert_eq!(interleaved, data);
}

// ============================================================================
// Convolutional Code Tests
// ============================================================================

#[test]

fn test_convolutional_encode_basic() {

    let data = vec![1, 0, 1, 1];

    let encoded = convolutional_encode(&data);

    // Rate 1/2, plus tail bits
    assert_eq!(encoded.len(), (data.len() + 2) * 2);
}

#[test]

fn test_convolutional_encode_all_zeros() {

    let data = vec![0, 0, 0, 0];

    let encoded = convolutional_encode(&data);

    // All zeros should produce all zeros
    assert!(encoded
        .iter()
        .all(|&x| x == 0));
}

// ============================================================================
// Code Theory Tests
// ============================================================================

#[test]

fn test_code_rate() {

    // Hamming(7,4) has rate 4/7
    assert!((code_rate(4, 7) - 4.0 / 7.0).abs() < 1e-10);
}

#[test]

fn test_code_rate_zero() {

    assert_eq!(code_rate(0, 7), 0.0);
}

#[test]

fn test_code_rate_n_zero() {

    assert_eq!(code_rate(4, 0), 0.0);
}

#[test]

fn test_error_correction_capability() {

    // Hamming (d=3) can correct 1 error
    assert_eq!(error_correction_capability(3), 1);

    // d=5 can correct 2 errors
    assert_eq!(error_correction_capability(5), 2);

    // d=7 can correct 3 errors
    assert_eq!(error_correction_capability(7), 3);
}

#[test]

fn test_error_detection_capability() {

    // Hamming (d=3) can detect 2 errors
    assert_eq!(error_detection_capability(3), 2);

    // d=5 can detect 4 errors
    assert_eq!(error_detection_capability(5), 4);
}

#[test]

fn test_minimum_distance() {

    let codewords = vec![
        vec![0, 0, 0, 0],
        vec![1, 1, 1, 1],
    ];

    assert_eq!(minimum_distance(&codewords), Some(4));
}

#[test]

fn test_minimum_distance_hamming_74() {

    // All Hamming(7,4) codewords have minimum distance 3
    let mut codewords = Vec::new();

    for d0 in 0..=1 {

        for d1 in 0..=1 {

            for d2 in 0..=1 {

                for d3 in 0..=1 {

                    let data = vec![d0, d1, d2, d3];

                    if let Some(cw) = hamming_encode_numerical(&data) {

                        codewords.push(cw);
                    }
                }
            }
        }
    }

    // Hamming(7,4) has minimum distance 3
    assert_eq!(minimum_distance(&codewords), Some(3));
}

#[test]

fn test_minimum_distance_single_codeword() {

    let codewords = vec![vec![1, 1, 1, 1]];

    assert_eq!(minimum_distance(&codewords), None);
}

// ============================================================================
// PolyGF256 Tests
// ============================================================================

#[test]

fn test_poly_gf256_new() {

    let poly = PolyGF256::new(vec![1, 2, 3]);

    assert_eq!(poly.0, vec![1, 2, 3]);
}

#[test]

fn test_poly_gf256_degree() {

    let poly = PolyGF256::new(vec![1, 2, 3]);

    assert_eq!(poly.degree(), 2);
}

#[test]

fn test_poly_gf256_degree_empty() {

    let poly = PolyGF256::new(vec![]);

    assert_eq!(poly.degree(), 0);
}

#[test]

fn test_poly_gf256_eval() {

    // p(x) = 1 + x (coefficients in descending order: [1, 1])
    let poly = PolyGF256::new(vec![1, 1]);

    // p(0) = 1
    assert_eq!(poly.eval(0), 1);
}

#[test]

fn test_poly_gf256_add() {

    let p1 = PolyGF256::new(vec![1, 2, 3]);

    let p2 = PolyGF256::new(vec![1, 1, 1]);

    let sum = p1.poly_add(&p2);

    // XOR addition
    assert_eq!(sum.0, vec![0, 3, 2]);
}

#[test]

fn test_poly_gf256_derivative() {

    // p(x) = x^3 + x^2 + x + 1
    let poly = PolyGF256::new(vec![1, 1, 1, 1]);

    let deriv = poly.derivative();

    // In GF(2), only odd powers survive: x^3 + x
    assert_eq!(deriv.degree(), 2);
}

// ============================================================================
// Property-Based Tests
// ============================================================================

mod proptests {

    use super::*;
    use proptest::prelude::*;

    proptest! {
        /// Hamming encode followed by decode (no errors) returns original data
        #[test]
        fn prop_hamming_encode_decode_roundtrip(
            d0 in 0u8..=1,
            d1 in 0u8..=1,
            d2 in 0u8..=1,
            d3 in 0u8..=1
        ) {
            let data = vec![d0, d1, d2, d3];
            let codeword = hamming_encode_numerical(&data).unwrap();
            let (decoded, error_pos) = hamming_decode_numerical(&codeword).unwrap();
            prop_assert_eq!(decoded, data);
            prop_assert_eq!(error_pos, None);
        }

        /// Hamming code can correct single bit errors
        #[test]
        fn prop_hamming_single_error_correction(
            d0 in 0u8..=1,
            d1 in 0u8..=1,
            d2 in 0u8..=1,
            d3 in 0u8..=1,
            error_bit in 0usize..7
        ) {
            let data = vec![d0, d1, d2, d3];
            let mut codeword = hamming_encode_numerical(&data).unwrap();
            codeword[error_bit] ^= 1; // Introduce single error
            let (decoded, error_pos) = hamming_decode_numerical(&codeword).unwrap();
            prop_assert_eq!(decoded, data);
            prop_assert_eq!(error_pos, Some(error_bit + 1));
        }

        /// Valid Hamming codeword passes check
        #[test]
        fn prop_hamming_valid_codeword_passes_check(
            d0 in 0u8..=1,
            d1 in 0u8..=1,
            d2 in 0u8..=1,
            d3 in 0u8..=1
        ) {
            let data = vec![d0, d1, d2, d3];
            let codeword = hamming_encode_numerical(&data).unwrap();
            prop_assert!(hamming_check_numerical(&codeword));
        }

        /// Corrupted Hamming codeword fails check
        #[test]
        fn prop_hamming_corrupted_codeword_fails_check(
            d0 in 0u8..=1,
            d1 in 0u8..=1,
            d2 in 0u8..=1,
            d3 in 0u8..=1,
            error_bit in 0usize..7
        ) {
            let data = vec![d0, d1, d2, d3];
            let mut codeword = hamming_encode_numerical(&data).unwrap();
            codeword[error_bit] ^= 1;
            prop_assert!(!hamming_check_numerical(&codeword));
        }

        /// Hamming distance is symmetric
        #[test]
        fn prop_hamming_distance_symmetric(
            a in proptest::collection::vec(0u8..=1, 10),
            b in proptest::collection::vec(0u8..=1, 10)
        ) {
            let dist_ab = hamming_distance_numerical(&a, &b);
            let dist_ba = hamming_distance_numerical(&b, &a);
            prop_assert_eq!(dist_ab, dist_ba);
        }

        /// Hamming distance to self is zero
        #[test]
        fn prop_hamming_distance_self_zero(a in proptest::collection::vec(0u8..=1, 1..20)) {
            let dist = hamming_distance_numerical(&a, &a);
            prop_assert_eq!(dist, Some(0));
        }

        /// Hamming weight equals distance from zero vector
        #[test]
        fn prop_hamming_weight_equals_distance_from_zero(
            data in proptest::collection::vec(0u8..=1, 1..20)
        ) {
            let zeros = vec![0u8; data.len()];
            let weight = hamming_weight_numerical(&data);
            let dist = hamming_distance_numerical(&data, &zeros);
            prop_assert_eq!(Some(weight), dist);
        }

        /// CRC32 is deterministic
        #[test]
        fn prop_crc32_deterministic(data in proptest::collection::vec(any::<u8>(), 0..100)) {
            let crc1 = crc32_compute_numerical(&data);
            let crc2 = crc32_compute_numerical(&data);
            prop_assert_eq!(crc1, crc2);
        }

        /// CRC32 verify accepts correct CRC
        #[test]
        fn prop_crc32_verify_accepts_correct(data in proptest::collection::vec(any::<u8>(), 0..100)) {
            let crc = crc32_compute_numerical(&data);
            prop_assert!(crc32_verify_numerical(&data, crc));
        }

        /// CRC32 verify rejects wrong CRC (with high probability)
        #[test]
        fn prop_crc32_verify_rejects_wrong(
            data in proptest::collection::vec(any::<u8>(), 1..100),
            wrong_crc in any::<u32>()
        ) {
            let correct_crc = crc32_compute_numerical(&data);
            if wrong_crc != correct_crc {
                prop_assert!(!crc32_verify_numerical(&data, wrong_crc));
            }
        }

        /// Interleave followed by deinterleave returns original data
        #[test]
        fn prop_interleave_roundtrip(
            data in proptest::collection::vec(any::<u8>(), 1..50),
            depth in 1usize..10
        ) {
            let interleaved = interleave(&data, depth);
            let deinterleaved = deinterleave(&interleaved, depth);
            prop_assert_eq!(deinterleaved, data);
        }

        /// Interleaving preserves data length
        #[test]
        fn prop_interleave_preserves_length(
            data in proptest::collection::vec(any::<u8>(), 0..50),
            depth in 1usize..10
        ) {
            let interleaved = interleave(&data, depth);
            prop_assert_eq!(interleaved.len(), data.len());
        }

        /// Code rate is in range [0, 1] for valid parameters
        #[test]
        fn prop_code_rate_range(k in 0usize..100, n in 1usize..100) {
            let k = k.min(n);
            let rate = code_rate(k, n);
            prop_assert!(rate >= 0.0 && rate <= 1.0);
        }

        /// Error correction capability is consistent with detection capability
        #[test]
        fn prop_correction_vs_detection(d in 1usize..20) {
            let t = error_correction_capability(d);
            let s = error_detection_capability(d);
            // Detection capability should be >= 2 * correction capability
            prop_assert!(s >= 2 * t);
        }

        /// Reed-Solomon encode produces codeword of correct length
        #[test]
        fn prop_rs_encode_length(
            message in proptest::collection::vec(any::<u8>(), 1..50),
            n_parity in 2usize..10
        ) {
            if message.len() + n_parity <= 255 {
                let codeword = reed_solomon_encode(&message, n_parity).unwrap();
                prop_assert_eq!(codeword.len(), message.len() + n_parity);
            }
        }

        /// Reed-Solomon valid codeword passes check
        #[test]
        fn prop_rs_valid_passes_check(
            message in proptest::collection::vec(any::<u8>(), 1..20),
            n_parity in 2usize..6
        ) {
            if message.len() + n_parity <= 255 {
                let codeword = reed_solomon_encode(&message, n_parity).unwrap();
                prop_assert!(reed_solomon_check(&codeword, n_parity));
            }
        }
    }
}
