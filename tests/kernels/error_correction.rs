//! Tests for the error-correction kernels (ported from
//! `numerical_error_correction_test.rs`).
//!
//! This module contains comprehensive unit tests and property-based tests
//! for Reed-Solomon codes, Hamming codes, CRC checksums, and other
//! error correction utilities.

use rssn::kernels::error_correction::*;

// ============================================================================
// Reed-Solomon Tests
// ============================================================================

#[test]

fn test_reed_solomon_encode_basic() {
    let message = vec![0x01, 0x02, 0x03, 0x04];

    let codeword = reed_solomon_encode(&message, 4).unwrap_or_else(|e| panic!("{e}"));

    assert_eq!(codeword.len(), 8); // 4 data + 4 parity
    // First 4 bytes should be the original message
    assert_eq!(&codeword[..4], &message);
}

#[test]

fn test_reed_solomon_encode_empty_message() {
    let message: Vec<u8> = vec![];

    let codeword = reed_solomon_encode(&message, 4).unwrap_or_else(|e| panic!("{e}"));

    assert_eq!(codeword.len(), 4); // 0 data + 4 parity
}

#[test]

fn test_reed_solomon_encode_max_length() {
    // Maximum message length is 255 - n_parity
    let message: Vec<u8> = (0..251).collect();

    let codeword = reed_solomon_encode(&message, 4).unwrap_or_else(|e| panic!("{e}"));

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
    let message = vec![0x01, 0x02, 0x03, 0x04];

    let codeword = reed_solomon_encode(&message, 4).unwrap_or_else(|e| panic!("{e}"));

    assert!(reed_solomon_check(&codeword, 4));
}

#[test]

fn test_reed_solomon_check_corrupted() {
    let message = vec![0x01, 0x02, 0x03, 0x04];

    let mut codeword = reed_solomon_encode(&message, 4).unwrap_or_else(|e| panic!("{e}"));

    codeword[0] ^= 0xFF; // Corrupt first byte
    assert!(!reed_solomon_check(&codeword, 4));
}

#[test]

fn test_reed_solomon_decode_no_errors() {
    let message = vec![0x01, 0x02, 0x03, 0x04];

    let mut codeword = reed_solomon_encode(&message, 4).unwrap_or_else(|e| panic!("{e}"));

    reed_solomon_decode(&mut codeword, 4).unwrap_or_else(|e| panic!("{e}"));

    assert_eq!(&codeword[..4], &message);
}

#[test]

fn test_reed_solomon_decode_with_errors() {
    let message = vec![0x01, 0x02, 0x03, 0x04];

    let codeword = reed_solomon_encode(&message, 4).unwrap_or_else(|e| panic!("{e}"));

    // Introduce an error
    let mut corrupted = codeword.clone();

    corrupted[0] ^= 0xFF;

    // Decode and correct
    reed_solomon_decode(&mut corrupted, 4).unwrap_or_else(|e| panic!("{e}"));

    assert_eq!(&corrupted[0..message.len()], &message[..],);
}

// ============================================================================
// Hamming Code Tests
// ============================================================================

#[test]

fn test_hamming_encode_basic() {
    let data = vec![1, 0, 1, 1];

    let codeword = hamming_encode_numerical(&data)
        .unwrap_or_else(|| panic!("hamming_encode_numerical returned None"));

    assert_eq!(codeword.len(), 7);
}

#[test]

fn test_hamming_encode_all_zeros() {
    let data = vec![0, 0, 0, 0];

    let codeword = hamming_encode_numerical(&data)
        .unwrap_or_else(|| panic!("hamming_encode_numerical returned None"));

    assert_eq!(codeword, vec![0, 0, 0, 0, 0, 0, 0]);
}

#[test]

fn test_hamming_encode_all_ones() {
    let data = vec![1, 1, 1, 1];

    let codeword = hamming_encode_numerical(&data)
        .unwrap_or_else(|| panic!("hamming_encode_numerical returned None"));

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

    let codeword = hamming_encode_numerical(&data)
        .unwrap_or_else(|| panic!("hamming_encode_numerical returned None"));

    let (decoded, error_pos) =
        hamming_decode_numerical(&codeword).unwrap_or_else(|e| panic!("{e}"));

    assert_eq!(decoded, data);

    assert_eq!(error_pos, None);
}

#[test]

fn test_hamming_decode_single_error() {
    let data = vec![1, 0, 1, 1];

    let mut codeword = hamming_encode_numerical(&data)
        .unwrap_or_else(|| panic!("hamming_encode_numerical returned None"));

    codeword[2] ^= 1; // Introduce error at position 3 (1-indexed)
    let (decoded, error_pos) =
        hamming_decode_numerical(&codeword).unwrap_or_else(|e| panic!("{e}"));

    assert_eq!(decoded, data);

    assert_eq!(error_pos, Some(3));
}

#[test]

fn test_hamming_decode_parity_error() {
    let data = vec![1, 0, 1, 1];

    let mut codeword = hamming_encode_numerical(&data)
        .unwrap_or_else(|| panic!("hamming_encode_numerical returned None"));

    codeword[0] ^= 1; // Error in parity bit
    let (decoded, error_pos) =
        hamming_decode_numerical(&codeword).unwrap_or_else(|e| panic!("{e}"));

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

    let codeword = hamming_encode_numerical(&data)
        .unwrap_or_else(|| panic!("hamming_encode_numerical returned None"));

    assert!(hamming_check_numerical(&codeword));
}

#[test]

fn test_hamming_check_invalid() {
    let data = vec![1, 0, 1, 1];

    let mut codeword = hamming_encode_numerical(&data)
        .unwrap_or_else(|| panic!("hamming_encode_numerical returned None"));

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

    let decoded = bch_decode(&codeword, 2).unwrap_or_else(|e| panic!("{e}"));

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
    assert_eq!(crc, 0xBB3D);
}

#[test]

fn test_crc16_empty() {
    let data: &[u8] = b"";

    let crc = crc16_compute(data);

    // CRC-16 of empty string
    assert_eq!(crc, 0x0000);
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
    let data = vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12];

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
    assert!(encoded.iter().all(|&x| x == 0));
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
    let codewords = vec![vec![0, 0, 0, 0], vec![1, 1, 1, 1]];

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

    use proptest::prelude::*;

    use super::*;

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
            let codeword = hamming_encode_numerical(&data).unwrap_or_else(|| panic!("hamming_encode_numerical returned None"));
            let (decoded, error_pos) = hamming_decode_numerical(&codeword).unwrap_or_else(|e| panic!("{e}"));
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
            let mut codeword = hamming_encode_numerical(&data).unwrap_or_else(|| panic!("hamming_encode_numerical returned None"));
            codeword[error_bit] ^= 1; // Introduce single error
            let (decoded, error_pos) = hamming_decode_numerical(&codeword).unwrap_or_else(|e| panic!("{e}"));
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
            let codeword = hamming_encode_numerical(&data).unwrap_or_else(|| panic!("hamming_encode_numerical returned None"));
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
            let mut codeword = hamming_encode_numerical(&data).unwrap_or_else(|| panic!("hamming_encode_numerical returned None"));
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
            prop_assert!((0.0..=1.0).contains(&rate));
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
                let codeword = reed_solomon_encode(&message, n_parity).unwrap_or_else(|e| panic!("{e}"));
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
                let codeword = reed_solomon_encode(&message, n_parity).unwrap_or_else(|e| panic!("{e}"));
                prop_assert!(reed_solomon_check(&codeword, n_parity));
            }
        }
    }
}

// ============================================================================
// Added: stronger assertions with independently known values
// ============================================================================

mod strengthened {
    use proptest::prelude::*;
    use proptest::test_runner::RngSeed;
    use rssn::kernels::error_correction::*;

    fn cfg() -> ProptestConfig {
        ProptestConfig {
            rng_seed: RngSeed::Fixed(0x5EED),
            failure_persistence: None,
            ..ProptestConfig::default()
        }
    }

    #[test]
    fn reed_solomon_corrects_a_single_error_and_restores_the_whole_codeword() {
        let message = vec![0x01, 0x02, 0x03, 0x04];
        let codeword = reed_solomon_encode(&message, 4).unwrap_or_else(|e| panic!("{e}"));
        let mut corrupted = codeword.clone();
        corrupted[0] ^= 0xFF;
        reed_solomon_decode(&mut corrupted, 4).unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(corrupted, codeword, "parity bytes must be repaired too");
        assert!(reed_solomon_check(&corrupted, 4));
    }

    #[test]
    fn reed_solomon_corrects_errors_in_parity_bytes() {
        let message = b"parity errors".to_vec();
        let codeword = reed_solomon_encode(&message, 6).unwrap_or_else(|e| panic!("{e}"));
        let mut corrupted = codeword.clone();
        let n = corrupted.len();
        corrupted[n - 1] ^= 0x5A;
        reed_solomon_decode(&mut corrupted, 6).unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(corrupted, codeword);
    }

    #[test]
    fn reed_solomon_corrects_up_to_half_the_parity_symbols() {
        let message: Vec<u8> = (10..30).collect();
        let n_parity = 4; // t = 2
        let codeword = reed_solomon_encode(&message, n_parity).unwrap_or_else(|e| panic!("{e}"));
        let mut corrupted = codeword.clone();
        corrupted[3] ^= 0x21;
        corrupted[17] ^= 0xC4;
        reed_solomon_decode(&mut corrupted, n_parity).unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(corrupted, codeword);
    }

    #[test]
    fn reed_solomon_syndromes_vanish_exactly_for_valid_codewords() {
        let codeword = reed_solomon_encode(&[9, 8, 7, 6, 5], 4).unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(calculate_syndromes(&codeword, 4), vec![0; 4]);
        let mut bad = codeword;
        bad[1] ^= 1;
        assert!(calculate_syndromes(&bad, 4).iter().any(|&s| s != 0));
    }

    #[test]
    fn hamming_encode_known_codewords() {
        // layout [p1 p2 d3 p4 d5 d6 d7] with p1 = d3^d5^d7, p2 = d3^d6^d7, p4 = d5^d6^d7
        assert_eq!(
            hamming_encode_numerical(&[1, 0, 1, 1]),
            Some(vec![0, 1, 1, 0, 0, 1, 1])
        );
        assert_eq!(
            hamming_encode_numerical(&[1, 1, 1, 1]),
            Some(vec![1, 1, 1, 1, 1, 1, 1])
        );
        assert_eq!(
            hamming_encode_numerical(&[0, 0, 0, 1]),
            Some(vec![1, 1, 0, 1, 0, 0, 1])
        );
    }

    #[test]
    fn hamming_decode_reports_one_indexed_error_position_for_every_bit() {
        let data = vec![0, 1, 1, 0];
        let clean = hamming_encode_numerical(&data).unwrap_or_default();
        for i in 0..7 {
            let mut cw = clean.clone();
            cw[i] ^= 1;
            assert_eq!(
                hamming_decode_numerical(&cw),
                Ok((data.clone(), Some(i + 1))),
                "bit {i}"
            );
        }
    }

    #[test]
    fn bch_corrects_a_single_data_bit_error() {
        let data = vec![1, 0, 1, 1, 0, 1, 0, 0];
        let mut cw = bch_encode(&data, 2);
        // 8 data bits need the (31, 21) code: 10 parity bits for t = 2
        assert_eq!(cw.len(), data.len() + 10);
        cw[3] ^= 1;
        assert_eq!(bch_decode(&cw, 2), Ok(data));
        assert!(bch_decode(&[1], 2).is_err());
    }

    #[test]
    fn crc_reference_values() {
        // CRC-32/ISO-HDLC check value.
        assert_eq!(crc32_compute_numerical(b"123456789"), 0xCBF4_3926);
        // CRC-16/ARC check value.
        assert_eq!(crc16_compute(b"123456789"), 0xBB3D);
        // CRC-8 (poly 0x07, init 0) check value.
        assert_eq!(crc8_compute(b"123456789"), 0xF4);
        assert_eq!(crc32_compute_numerical(b"a"), 0xE8B7_BE43);
    }

    #[test]
    fn crc32_detects_every_single_bit_flip() {
        let data = b"The quick brown fox".to_vec();
        let good = crc32_compute_numerical(&data);
        for byte in 0..data.len() {
            for bit in 0..8 {
                let mut d = data.clone();
                d[byte] ^= 1 << bit;
                assert_ne!(crc32_compute_numerical(&d), good, "byte {byte} bit {bit}");
            }
        }
    }

    #[test]
    fn interleave_known_permutation() {
        // Write row-wise into a depth x width grid, read column-wise.
        let interleaved = interleave(&[1, 2, 3, 4, 5, 6], 3);
        assert_eq!(interleaved, vec![1, 4, 2, 5, 3, 6]);
        assert_eq!(deinterleave(&interleaved, 3), vec![1, 2, 3, 4, 5, 6]);
    }

    #[test]
    fn convolutional_encoder_known_output_and_linearity() {
        // Single 1 followed by zeros produces the impulse response of the two generators.
        let out = convolutional_encode(&[1, 0, 0, 0]);
        assert_eq!(out, vec![1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0]);
        let a = convolutional_encode(&[1, 0, 1, 1]);
        let b = convolutional_encode(&[0, 1, 1, 0]);
        let sum = convolutional_encode(&[1, 1, 0, 1]);
        let xor: Vec<u8> = a.iter().zip(&b).map(|(x, y)| x ^ y).collect();
        assert_eq!(sum, xor, "the code is linear over GF(2)");
    }

    #[test]
    fn minimum_distance_of_repetition_code_and_capabilities() {
        assert_eq!(minimum_distance(&[vec![0; 5], vec![1; 5]]), Some(5));
        assert_eq!(error_correction_capability(5), 2);
        assert_eq!(error_detection_capability(1), 0);
        assert_eq!(code_rate(1, 3), 1.0 / 3.0);
    }

    #[test]
    fn gf256_polynomial_operations() {
        let p = PolyGF256::new(vec![1, 1]);
        assert_eq!(p.eval(2), 3); // 2 + 1 in GF(2^8) is XOR
        assert_eq!(p.eval(0), 1);
        let a = PolyGF256::new(vec![7, 7, 7]);
        assert_eq!(a.poly_add(&a).0, vec![0, 0, 0]);
        // Formal derivative in characteristic 2 keeps only odd-power terms.
        let d = PolyGF256::new(vec![1, 1, 1, 1]).derivative();
        assert_eq!(d.degree(), 2);
    }

    proptest! {
        #![proptest_config(cfg())]

        /// A single corrupted symbol is always repaired.
        #[test]
        fn prop_reed_solomon_corrects_one_random_error(
            message in proptest::collection::vec(any::<u8>(), 4..24),
            e in 1u8..=255, p in 0usize..64, n_parity in 2usize..8,
        ) {
            let codeword = reed_solomon_encode(&message, n_parity).map_err(TestCaseError::fail)?;
            let mut corrupted = codeword.clone();
            corrupted[p % codeword.len()] ^= e;
            reed_solomon_decode(&mut corrupted, n_parity).map_err(TestCaseError::fail)?;
            prop_assert_eq!(corrupted, codeword);
        }

        /// RS(n, k) with 4 parity symbols corrects any 2 symbol errors at random positions.
        #[test]
        fn prop_reed_solomon_corrects_random_errors(
            message in proptest::collection::vec(any::<u8>(), 4..24),
            e1 in 1u8..=255, e2 in 1u8..=255,
            p1 in 0usize..64, p2 in 0usize..64,
        ) {
            let n_parity = 4;
            let codeword = reed_solomon_encode(&message, n_parity).map_err(TestCaseError::fail)?;
            let n = codeword.len();
            let (i, j) = (p1 % n, p2 % n);
            prop_assume!(i != j);
            let mut corrupted = codeword.clone();
            corrupted[i] ^= e1;
            corrupted[j] ^= e2;
            reed_solomon_decode(&mut corrupted, n_parity).map_err(TestCaseError::fail)?;
            prop_assert_eq!(corrupted, codeword);
        }

        /// The encoded codeword always starts with the message (systematic code).
        #[test]
        fn prop_reed_solomon_is_systematic(message in proptest::collection::vec(any::<u8>(), 1..30), n_parity in 2usize..8) {
            let cw = reed_solomon_encode(&message, n_parity).map_err(TestCaseError::fail)?;
            prop_assert_eq!(&cw[..message.len()], &message[..]);
        }

        /// CRC-32 streaming over an arbitrary split equals one-shot CRC-32.
        #[test]
        fn prop_crc32_streaming_equals_one_shot(data in proptest::collection::vec(any::<u8>(), 0..64), split in 0usize..64) {
            let split = split.min(data.len());
            let crc = crc32_update_numerical(0xFFFF_FFFF, &data[..split]);
            let crc = crc32_update_numerical(crc, &data[split..]);
            prop_assert_eq!(crc32_finalize_numerical(crc), crc32_compute_numerical(&data));
        }

        /// Hamming(7,4) minimum distance: two different data words give codewords at distance >= 3.
        #[test]
        fn prop_hamming_distinct_messages_are_far_apart(a in 0u8..16, b in 0u8..16) {
            prop_assume!(a != b);
            let bits = |x: u8| vec![(x >> 3) & 1, (x >> 2) & 1, (x >> 1) & 1, x & 1];
            let ca = hamming_encode_numerical(&bits(a)).unwrap_or_default();
            let cb = hamming_encode_numerical(&bits(b)).unwrap_or_default();
            prop_assert!(hamming_distance_numerical(&ca, &cb).unwrap_or(0) >= 3);
        }

        /// The convolutional code is linear over GF(2).
        #[test]
        fn prop_convolutional_linear(a in proptest::collection::vec(0u8..=1, 6), b in proptest::collection::vec(0u8..=1, 6)) {
            let x: Vec<u8> = a.iter().zip(&b).map(|(p, q)| p ^ q).collect();
            let (ea, eb, ex) = (convolutional_encode(&a), convolutional_encode(&b), convolutional_encode(&x));
            let xor: Vec<u8> = ea.iter().zip(&eb).map(|(p, q)| p ^ q).collect();
            prop_assert_eq!(ex, xor);
        }
    }
}

mod ledger_fill {
    use rssn::kernels::error_correction::*;

    #[test]
    fn poly_sub_is_xor_like_poly_add() {
        let a = PolyGF256::new(vec![1, 2, 3]);
        let b = PolyGF256::new(vec![3, 2, 1]);
        assert_eq!(a.poly_sub(&b), a.poly_add(&b));
        assert_eq!(a.poly_sub(&a).0, vec![0, 0, 0]);
    }

    #[test]
    fn poly_mul_and_scale_follow_the_field_rules() {
        // (1 + x)(1 + x) = 1 + x^2 in characteristic 2 (coefficients in ascending order).
        let p = PolyGF256::new(vec![1, 1]);
        assert_eq!(p.poly_mul(&p).0, vec![1, 0, 1]);
        assert_eq!(p.scale(1), p);
        assert_eq!(p.scale(0).0, vec![0, 0]);
    }

    #[test]
    fn poly_div_by_the_zero_polynomial_is_an_error() {
        assert!(
            PolyGF256::new(vec![1, 0, 1])
                .poly_div(&PolyGF256::new(vec![]))
                .is_err()
        );
    }

    #[test]
    fn poly_div_returns_the_correct_quotient() {
        // Ascending order: x^3 = [0, 0, 0, 1], x = [0, 1]; x^3 / x = x^2 = [0, 0, 1].
        let (q, r) = PolyGF256::new(vec![0, 0, 0, 1])
            .poly_div(&PolyGF256::new(vec![0, 1]))
            .unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(q.normalize().0, vec![0, 0, 1]);
        assert!(r.normalize().0.is_empty());
    }

    #[test]
    fn chien_search_finds_the_error_locator_root_position() {
        // sigma(x) = 1 + alpha^3 x has its root at alpha^-3, i.e. error position 3.
        let sigma = PolyGF256::new(vec![1, 8]);
        assert_eq!(
            chien_search(&sigma).unwrap_or_else(|e| panic!("{e}")),
            vec![3]
        );
    }

    #[test]
    fn forney_algorithm_returns_the_error_magnitude() {
        // sigma = 1 + alpha^3 x, omega = alpha^6: e = omega(X^-1) X^-1 / sigma'(X^-1) = alpha^(6 - 3 - 3) = 1.
        let sigma = PolyGF256::new(vec![1, 8]);
        let omega = PolyGF256::new(vec![64]);
        assert_eq!(
            forney_algorithm(&omega, &sigma, &[3]).unwrap_or_else(|e| panic!("{e}")),
            vec![1]
        );
    }
}

mod ledger_rs_extra {
    use rssn::kernels::error_correction::*;

    #[test]
    fn reed_solomon_corrects_every_two_error_pair_in_a_short_code() {
        // Exhaustive over all position pairs with fixed magnitudes.
        let cw = reed_solomon_encode(&[3, 1, 4, 1, 5, 9], 4).unwrap();
        for i in 0..cw.len() {
            for j in (i + 1)..cw.len() {
                let mut bad = cw.clone();
                bad[i] ^= 0x35;
                bad[j] ^= 0xE1;
                reed_solomon_decode(&mut bad, 4).unwrap_or_else(|e| panic!("({i},{j}): {e}"));
                assert_eq!(bad, cw, "({i},{j})");
            }
        }
    }

    #[test]
    fn reed_solomon_corrects_three_errors_with_six_parity() {
        let cw = reed_solomon_encode(b"three errors here", 6).unwrap();
        let mut bad = cw.clone();
        bad[0] ^= 0x01;
        bad[7] ^= 0x80;
        bad[cw.len() - 2] ^= 0x7F;
        reed_solomon_decode(&mut bad, 6).unwrap();
        assert_eq!(bad, cw);
    }

    #[test]
    fn reed_solomon_never_returns_ok_with_a_non_codeword() {
        // t + 1 = 3 errors with 4 parity symbols: either an error is reported
        // (input untouched) or the result is at least a valid codeword.
        let cw = reed_solomon_encode(&(0u8..16).collect::<Vec<_>>(), 4).unwrap();
        let mut errs = 0;
        for k in 0..cw.len() - 2 {
            let mut bad = cw.clone();
            bad[k] ^= 0x11;
            bad[k + 1] ^= 0x22;
            bad[k + 2] ^= 0x44;
            let before = bad.clone();
            match reed_solomon_decode(&mut bad, 4) {
                | Ok(()) => assert!(reed_solomon_check(&bad, 4)),
                | Err(_) => {
                    errs += 1;
                    assert_eq!(bad, before);
                },
            }
        }
        assert!(
            errs > 0,
            "some 3-error patterns must be reported as uncorrectable"
        );
    }

    #[test]
    fn poly_div_reconstructs_the_dividend() {
        // p = q d + r, ascending order, in GF(2^8).
        let p = PolyGF256::new(vec![7, 0, 19, 200, 3, 55]);
        let d = PolyGF256::new(vec![5, 9, 1]);
        let (q, r) = p.poly_div(&d).unwrap();
        assert!(r.normalize().degree() < d.degree() || r.normalize().0.is_empty());
        assert_eq!(q.poly_mul(&d).poly_add(&r).normalize(), p.normalize());
    }

    #[test]
    fn poly_add_aligns_constant_terms() {
        // (1 + 2x) + (3) = 2 + 2x, not 1 + (2^3)x.
        let sum = PolyGF256::new(vec![1, 2]).poly_add(&PolyGF256::new(vec![3]));
        assert_eq!(sum.0, vec![2, 2]);
    }
}

mod bch_codes {
    use rssn::kernels::error_correction::{Bch, bch_decode, bch_encode};

    struct Rng(u64);

    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }

        fn below(&mut self, n: usize) -> usize {
            (self.next() % n as u64) as usize
        }

        fn bits(&mut self, n: usize) -> Vec<u8> {
            (0..n).map(|_| (self.next() & 1) as u8).collect()
        }

        /// `count` distinct positions in `0..n`.
        fn positions(&mut self, n: usize, count: usize) -> Vec<usize> {
            let mut v = Vec::new();
            while v.len() < count {
                let p = self.below(n);
                if !v.contains(&p) {
                    v.push(p);
                }
            }
            v
        }
    }

    fn gen_as_int(code: &Bch) -> u64 {
        code.generator().iter().enumerate().map(|(i, &b)| u64::from(b) << i).sum()
    }

    #[test]
    fn parameters_match_the_classical_tables() {
        for &(m, t, n, k) in &[
            (3, 1, 7, 4),
            (4, 1, 15, 11),
            (4, 2, 15, 7),
            (4, 3, 15, 5),
            (5, 1, 31, 26),
            (5, 2, 31, 21),
            (5, 3, 31, 16),
            (5, 5, 31, 11),
            (6, 3, 63, 45),
            (6, 4, 63, 39),
            (7, 5, 127, 92),
            (8, 4, 255, 223),
        ] {
            let c = Bch::new(m, t).unwrap();
            assert_eq!((c.n(), c.k(), c.t(), c.m()), (n, k, t, m), "BCH({n},{k},{t})");
            assert_eq!(c.parity_len(), n - k);
        }
    }

    #[test]
    fn generator_polynomials_match_references() {
        // octal 721, 2467 and 107657 (Lin & Costello, Table of BCH generators)
        assert_eq!(gen_as_int(&Bch::new(4, 2).unwrap()), 0o721);
        assert_eq!(gen_as_int(&Bch::new(4, 3).unwrap()), 0o2467);
        assert_eq!(gen_as_int(&Bch::new(5, 3).unwrap()), 0o107_657);
        // t = 1 gives the primitive polynomial itself (the Hamming code)
        assert_eq!(gen_as_int(&Bch::new(4, 1).unwrap()), 0b10011);
    }

    #[test]
    fn encoding_is_systematic_and_valid() {
        let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
        for &(m, t) in &[(4, 2), (5, 3), (6, 3), (7, 4)] {
            let c = Bch::new(m, t).unwrap();
            for _ in 0..20 {
                let msg = rng.bits(c.k());
                let cw = c.encode(&msg).unwrap();
                assert_eq!(cw.len(), c.n());
                assert_eq!(&cw[..c.k()], &msg[..]);
                assert!(c.is_codeword(&cw));
            }
        }
    }

    #[test]
    fn codewords_form_a_linear_code_with_designed_distance() {
        let c = Bch::new(4, 2).unwrap(); // (15, 7, d >= 5)
        let mut min_weight = usize::MAX;
        for m in 1u32..(1 << 7) {
            let msg: Vec<u8> = (0..7).map(|i| ((m >> (6 - i)) & 1) as u8).collect();
            let cw = c.encode(&msg).unwrap();
            min_weight = min_weight.min(cw.iter().map(|&b| usize::from(b)).sum());
        }
        assert_eq!(min_weight, 5);
        // linearity
        let a = c.encode(&[1, 0, 1, 0, 0, 1, 1]).unwrap();
        let b = c.encode(&[0, 1, 1, 0, 1, 0, 1]).unwrap();
        let sum = c.encode(&[1, 1, 0, 0, 1, 1, 0]).unwrap();
        assert!(a.iter().zip(&b).zip(&sum).all(|((x, y), z)| x ^ y == *z));
    }

    #[test]
    fn corrects_every_pattern_of_up_to_t_errors_exhaustively_for_15_7_2() {
        let c = Bch::new(4, 2).unwrap();
        let msg = [1, 0, 0, 1, 1, 0, 1];
        let cw = c.encode(&msg).unwrap();
        let mut patterns = 0;
        for i in 0..15 {
            let mut one = cw.clone();
            one[i] ^= 1;
            assert_eq!(c.correct(&one).unwrap(), (cw.clone(), 1));
            patterns += 1;
            for j in i + 1..15 {
                let mut two = one.clone();
                two[j] ^= 1;
                assert_eq!(c.correct(&two).unwrap(), (cw.clone(), 2), "errors at {i},{j}");
                assert_eq!(c.decode(&two).unwrap(), msg);
                patterns += 1;
            }
        }
        assert_eq!(patterns, 15 + 105);
    }

    #[test]
    fn corrects_random_errors_up_to_t() {
        let mut rng = Rng(0xDEAD_BEEF_CAFE_F00D);
        for &(m, t, k) in &[(4, 2, 7), (5, 3, 16), (6, 3, 45), (6, 5, 36), (7, 5, 92), (8, 4, 223)] {
            let c = Bch::new(m, t).unwrap();
            assert_eq!(c.k(), k);
            for trial in 0..60 {
                let msg = rng.bits(k);
                let mut cw = c.encode(&msg).unwrap();
                let errors = trial % (t + 1);
                for p in rng.positions(c.n(), errors) {
                    cw[p] ^= 1;
                }
                let (fixed, fixed_count) = c.correct(&cw).unwrap();
                assert_eq!(fixed_count, errors);
                assert_eq!(c.decode(&cw).unwrap(), msg, "BCH({},{k},{t}) with {errors} errors", c.n());
                assert!(c.is_codeword(&fixed));
            }
        }
    }

    #[test]
    fn t_plus_one_errors_are_detected_or_land_on_another_codeword() {
        let mut rng = Rng(42);
        let c = Bch::new(5, 3).unwrap();
        let mut detected = 0;
        for _ in 0..200 {
            let msg = rng.bits(c.k());
            let mut cw = c.encode(&msg).unwrap();
            for p in rng.positions(c.n(), 4) {
                cw[p] ^= 1;
            }
            match c.correct(&cw) {
                Err(_) => detected += 1,
                Ok((fixed, _)) => assert!(c.is_codeword(&fixed)),
            }
        }
        // bounded-distance decoding detects the majority of weight-(t+1) patterns
        assert!(detected > 100, "detected only {detected}");
    }

    #[test]
    fn shortened_codes_work() {
        let mut rng = Rng(7);
        let c = Bch::new(6, 3).unwrap(); // (63, 45)
        for len in [1usize, 10, 33, 45] {
            let msg = rng.bits(len);
            let mut cw = c.encode(&msg).unwrap();
            assert_eq!(cw.len(), len + 18);
            for p in rng.positions(cw.len(), 3) {
                cw[p] ^= 1;
            }
            assert_eq!(c.decode(&cw).unwrap(), msg);
        }
    }

    #[test]
    fn rejects_bad_parameters_and_words() {
        assert!(Bch::new(1, 1).is_err());
        assert!(Bch::new(17, 1).is_err());
        assert!(Bch::new(4, 0).is_err());
        assert!(Bch::new(4, 8).is_err());
        let c = Bch::new(4, 2).unwrap();
        assert!(c.encode(&[]).is_err());
        assert!(c.encode(&[0; 8]).is_err());
        assert!(c.encode(&[2, 0, 0]).is_err());
        assert!(c.decode(&[0; 16]).is_err());
        assert!(c.decode(&[0; 8]).is_err());
        assert!(c.decode(&[3; 15]).is_err());
        assert!(Bch::for_data_len(1 << 20, 3).is_err());
    }

    #[test]
    fn legacy_functions_use_the_real_code() {
        let mut rng = Rng(99);
        for &(len, t) in &[(4usize, 1usize), (8, 2), (11, 1), (30, 3), (100, 4)] {
            let data = rng.bits(len);
            let cw = bch_encode(&data, t);
            assert!(cw.len() > len);
            assert_eq!(&cw[..len], &data[..]);
            assert_eq!(bch_decode(&cw, t).unwrap(), data);
            let mut bad = cw.clone();
            for p in rng.positions(bad.len(), t) {
                bad[p] ^= 1;
            }
            assert_eq!(bch_decode(&bad, t).unwrap(), data, "len {len} t {t}");
        }
        // two errors, previously uncorrectable
        let data = vec![1, 0, 1, 1, 0, 1, 0, 0];
        let mut cw = bch_encode(&data, 2);
        cw[0] ^= 1;
        cw[9] ^= 1;
        assert_eq!(bch_decode(&cw, 2).unwrap(), data);
        assert!(bch_encode(&[], 2).is_empty());
        assert!(bch_encode(&[1, 0], 0).is_empty());
        assert!(bch_encode(&[1, 7], 2).is_empty());
        assert!(bch_decode(&[1, 0, 1], 2).is_err());
    }
}
