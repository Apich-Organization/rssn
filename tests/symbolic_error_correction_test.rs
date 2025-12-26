//! Test suite for error_correction module (Hamming, Reed-Solomon, and CRC codes).

use rssn::symbolic::error_correction::*;

#[test]

fn test_hamming_encode() {

    // Encode 4 bits: [1, 0, 1, 1]
    let data = vec![1u8, 0, 1, 1];

    let codeword = hamming_encode(&data).expect("Should encode 4 bits");

    assert_eq!(codeword.len(), 7);

    // Verify parity bits are correct
    // d3=1, d5=0, d6=1, d7=1
    // p1 = d3 ^ d5 ^ d7 = 1 ^ 0 ^ 1 = 0
    // p2 = d3 ^ d6 ^ d7 = 1 ^ 1 ^ 1 = 1
    // p4 = d5 ^ d6 ^ d7 = 0 ^ 1 ^ 1 = 0
    // Codeword: [p1, p2, d3, p4, d5, d6, d7] = [0, 1, 1, 0, 0, 1, 1]
    assert_eq!(codeword, vec![0, 1, 1, 0, 0, 1, 1]);
}

#[test]

fn test_hamming_encode_invalid_length() {

    let data = vec![1u8, 0, 1]; // Only 3 bits
    assert!(hamming_encode(&data).is_none());
}

#[test]

fn test_hamming_decode_no_error() {

    let codeword = vec![
        0u8, 1, 1, 0, 0, 1, 1,
    ];

    let (data, error_pos) = hamming_decode(&codeword).expect("Should decode");

    assert_eq!(data, vec![1, 0, 1, 1]);

    assert_eq!(error_pos, None);
}

#[test]

fn test_hamming_decode_single_error() {

    // Introduce error at position 5 (d5 flipped from 0 to 1)
    let codeword = vec![
        0u8, 1, 1, 0, 1, 1, 1,
    ]; // Error at index 4 (position 5)
    let (data, error_pos) = hamming_decode(&codeword).expect("Should decode and correct");

    assert_eq!(data, vec![1, 0, 1, 1]); // Should be corrected
    assert_eq!(error_pos, Some(5)); // 1-based position
}

#[test]

fn test_hamming_decode_parity_error() {

    // Introduce error at parity bit p1 (position 1)
    let codeword = vec![
        1u8, 1, 1, 0, 0, 1, 1,
    ]; // Error at index 0 (position 1)
    let (data, error_pos) = hamming_decode(&codeword).expect("Should decode and correct");

    assert_eq!(data, vec![1, 0, 1, 1]); // Data should be unchanged
    assert_eq!(error_pos, Some(1));
}

#[test]

fn test_rs_encode() {

    let data = vec![
        0x12, 0x34, 0x56, 0x78,
    ];

    let n_sym = 4; // 4 error correction symbols (t=2 errors)
    let codeword = rs_encode(&data, n_sym).expect("Should encode");

    // Codeword should be data + n_sym bytes
    assert_eq!(codeword.len(), data.len() + n_sym);

    // First bytes should be original data
    assert_eq!(&codeword[..4], &data[..]);
}

#[test]

fn test_rs_encode_decode_no_error() {

    let data = vec![
        0x12, 0x34, 0x56, 0x78, 0x9A,
    ];

    let n_sym = 4;

    let codeword = rs_encode(&data, n_sym).expect("Should encode");

    let decoded = rs_decode(&codeword, n_sym).expect("Should decode");

    assert_eq!(decoded, data);
}

#[test]

fn test_rs_roundtrip() {

    let data = b"Hello, World!".to_vec();

    let n_sym = 8; // Can correct up to 4 errors
    let codeword = rs_encode(&data, n_sym).expect("Should encode");

    let decoded = rs_decode(&codeword, n_sym).expect("Should decode");

    assert_eq!(decoded, data);
}

#[test]

fn test_rs_encode_max_length() {

    // Maximum length is 255 total
    let data = vec![0u8; 250];

    let n_sym = 5;

    let result = rs_encode(&data, n_sym);

    assert!(result.is_ok());

    // Exceeding length should fail
    let data2 = vec![0u8; 252];

    let result2 = rs_encode(&data2, n_sym);

    assert!(result2.is_err());
}

// ============================================================================
// Hamming Code Extensions Tests
// ============================================================================

#[test]

fn test_hamming_distance_equal() {

    let a = vec![
        0u8, 1, 0, 1, 1, 0, 1,
    ];

    let b = vec![
        0u8, 1, 0, 1, 1, 0, 1,
    ];

    assert_eq!(hamming_distance(&a, &b), Some(0));
}

#[test]

fn test_hamming_distance_different() {

    let a = vec![
        0u8, 1, 0, 1, 1, 0, 1,
    ];

    let b = vec![
        1u8, 0, 0, 1, 0, 0, 1,
    ]; // 3 positions differ
    assert_eq!(hamming_distance(&a, &b), Some(3));
}

#[test]

fn test_hamming_distance_all_different() {

    let a = vec![0u8, 0, 0, 0];

    let b = vec![1u8, 1, 1, 1];

    assert_eq!(hamming_distance(&a, &b), Some(4));
}

#[test]

fn test_hamming_distance_length_mismatch() {

    let a = vec![0u8, 1, 0];

    let b = vec![0u8, 1, 0, 1];

    assert_eq!(hamming_distance(&a, &b), None);
}

#[test]

fn test_hamming_weight_all_zeros() {

    let data = vec![0u8, 0, 0, 0, 0];

    assert_eq!(hamming_weight(&data), 0);
}

#[test]

fn test_hamming_weight_all_ones() {

    let data = vec![1u8, 1, 1, 1, 1];

    assert_eq!(hamming_weight(&data), 5);
}

#[test]

fn test_hamming_weight_mixed() {

    let data = vec![
        1u8, 0, 1, 0, 1, 0, 1,
    ];

    assert_eq!(hamming_weight(&data), 4);
}

#[test]

fn test_hamming_weight_empty() {

    let data: Vec<u8> = vec![];

    assert_eq!(hamming_weight(&data), 0);
}

#[test]

fn test_hamming_check_valid_codeword() {

    // Create a valid codeword for data [1, 0, 1, 1]
    let codeword = vec![
        0u8, 1, 1, 0, 0, 1, 1,
    ];

    assert!(hamming_check(&codeword));
}

#[test]

fn test_hamming_check_invalid_codeword() {

    // Flip one bit to create an invalid codeword
    let codeword = vec![
        1u8, 1, 1, 0, 0, 1, 1,
    ]; // p1 flipped
    assert!(!hamming_check(&codeword));
}

#[test]

fn test_hamming_check_wrong_length() {

    let codeword = vec![0u8, 1, 1, 0, 0, 1]; // Only 6 bits
    assert!(!hamming_check(&codeword));
}

#[test]

fn test_hamming_check_another_valid() {

    // Encode [0, 0, 0, 0] and verify
    let data = vec![0u8, 0, 0, 0];

    let codeword = hamming_encode(&data).unwrap();

    assert!(hamming_check(&codeword));
}

// ============================================================================
// Reed-Solomon Enhancements Tests
// ============================================================================

#[test]

fn test_rs_check_valid_codeword() {

    let data = vec![
        0x12, 0x34, 0x56, 0x78,
    ];

    let n_sym = 4;

    let codeword = rs_encode(&data, n_sym).unwrap();

    assert!(rs_check(&codeword, n_sym));
}

#[test]

fn test_rs_check_invalid_codeword() {

    let data = vec![
        0x12, 0x34, 0x56, 0x78,
    ];

    let n_sym = 4;

    let mut codeword = rs_encode(&data, n_sym).unwrap();

    codeword[0] ^= 0xFF; // Corrupt first byte
    assert!(!rs_check(&codeword, n_sym));
}

#[test]

fn test_rs_check_after_decode() {

    let data = b"Test Data".to_vec();

    let n_sym = 8;

    let codeword = rs_encode(&data, n_sym).unwrap();

    // Verify it's valid before any manipulation
    assert!(rs_check(&codeword, n_sym));
}

#[test]

fn test_rs_error_count_no_errors() {

    let data = vec![
        0x12, 0x34, 0x56, 0x78,
    ];

    let n_sym = 4;

    let codeword = rs_encode(&data, n_sym).unwrap();

    assert_eq!(rs_error_count(&codeword, n_sym), 0);
}

#[test]

fn test_rs_error_count_with_errors() {

    let data = vec![
        0x12, 0x34, 0x56, 0x78,
    ];

    let n_sym = 4;

    let mut codeword = rs_encode(&data, n_sym).unwrap();

    codeword[0] ^= 0xFF; // Introduce 1 error
    let error_count = rs_error_count(&codeword, n_sym);

    assert!(error_count > 0);
}

// ============================================================================
// CRC-32 Tests
// ============================================================================

#[test]

fn test_crc32_compute_empty() {

    let data: Vec<u8> = vec![];

    let crc = crc32_compute(&data);

    assert_eq!(crc, 0x00000000);
}

#[test]

fn test_crc32_compute_known_value() {

    // "123456789" should produce CRC-32 = 0xCBF43926 (IEEE 802.3)
    let data = b"123456789";

    let crc = crc32_compute(data);

    assert_eq!(crc, 0xCBF43926);
}

#[test]

fn test_crc32_compute_hello() {

    let data = b"Hello, World!";

    let crc = crc32_compute(data);

    // Just verify it produces a non-zero checksum
    assert!(crc != 0);
}

#[test]

fn test_crc32_verify_valid() {

    let data = b"Test data for CRC";

    let crc = crc32_compute(data);

    assert!(crc32_verify(data, crc));
}

#[test]

fn test_crc32_verify_invalid() {

    let data = b"Test data for CRC";

    let crc = crc32_compute(data);

    assert!(!crc32_verify(data, crc ^ 1)); // Wrong checksum
}

#[test]

fn test_crc32_verify_modified_data() {

    let data = b"Test data for CRC";

    let crc = crc32_compute(data);

    let mut modified = data.to_vec();

    modified[0] ^= 1; // Flip one bit
    assert!(!crc32_verify(&modified, crc));
}

#[test]

fn test_crc32_update_incremental() {

    let data1 = b"Hello, ";

    let data2 = b"World!";

    let full_data = b"Hello, World!";

    // Compute CRC incrementally
    let mut crc = 0xFFFFFFFF;

    crc = crc32_update(crc, data1);

    crc = crc32_update(crc, data2);

    let final_crc = crc32_finalize(crc);

    // Compare with direct computation
    let direct_crc = crc32_compute(full_data);

    assert_eq!(final_crc, direct_crc);
}

#[test]

fn test_crc32_update_single_byte() {

    let data = vec![0x42];

    let mut crc = 0xFFFFFFFF;

    crc = crc32_update(crc, &data);

    let final_crc = crc32_finalize(crc);

    assert_eq!(final_crc, crc32_compute(&data));
}

#[test]

fn test_crc32_finalize_identity() {

    // crc32_finalize just inverts the CRC
    let crc = 0x12345678;

    assert_eq!(crc32_finalize(crc), !crc);
}

// ============================================================================
// Integration Tests
// ============================================================================

#[test]

fn test_hamming_roundtrip_with_check() {

    let data = vec![1u8, 1, 0, 1];

    let codeword = hamming_encode(&data).unwrap();

    // Check valid
    assert!(hamming_check(&codeword));

    // Decode without errors
    let (decoded, pos) = hamming_decode(&codeword).unwrap();

    assert_eq!(decoded, data);

    assert_eq!(pos, None);
}

#[test]

fn test_hamming_error_detection_distance() {

    let data = vec![1u8, 0, 1, 1];

    let codeword = hamming_encode(&data).unwrap();

    // Create corrupted codeword (single bit flip)
    let mut corrupted = codeword.clone();

    corrupted[3] ^= 1;

    // Check Hamming distance between original and corrupted
    let dist = hamming_distance(&codeword, &corrupted);

    assert_eq!(dist, Some(1));

    // Verify check fails
    assert!(!hamming_check(&corrupted));
}

#[test]

fn test_rs_full_pipeline() {

    let original_data = b"RS Error Correction Test".to_vec();

    let n_sym = 8;

    // Encode
    let codeword = rs_encode(&original_data, n_sym).unwrap();

    // Check valid
    assert!(rs_check(&codeword, n_sym));

    assert_eq!(rs_error_count(&codeword, n_sym), 0);

    // Decode
    let decoded = rs_decode(&codeword, n_sym).unwrap();

    assert_eq!(decoded, original_data);
}

#[test]

fn test_crc32_data_integrity() {

    let data = b"Important data that must not be corrupted".to_vec();

    // Compute checksum
    let checksum = crc32_compute(&data);

    // Verify
    assert!(crc32_verify(&data, checksum));

    // Corrupt and verify fails
    let mut corrupted = data.clone();

    corrupted[10] ^= 0x01;

    assert!(!crc32_verify(&corrupted, checksum));
}
