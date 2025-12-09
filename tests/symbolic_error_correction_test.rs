//! Test suite for error_correction module (Hamming and Reed-Solomon codes).

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
    let codeword = vec![0u8, 1, 1, 0, 0, 1, 1];
    let (data, error_pos) = hamming_decode(&codeword).expect("Should decode");
    assert_eq!(data, vec![1, 0, 1, 1]);
    assert_eq!(error_pos, None);
}

#[test]
fn test_hamming_decode_single_error() {
    // Introduce error at position 5 (d5 flipped from 0 to 1)
    let codeword = vec![0u8, 1, 1, 0, 1, 1, 1]; // Error at index 4 (position 5)
    let (data, error_pos) = hamming_decode(&codeword).expect("Should decode and correct");
    assert_eq!(data, vec![1, 0, 1, 1]); // Should be corrected
    assert_eq!(error_pos, Some(5)); // 1-based position
}

#[test]
fn test_hamming_decode_parity_error() {
    // Introduce error at parity bit p1 (position 1)
    let codeword = vec![1u8, 1, 1, 0, 0, 1, 1]; // Error at index 0 (position 1)
    let (data, error_pos) = hamming_decode(&codeword).expect("Should decode and correct");
    assert_eq!(data, vec![1, 0, 1, 1]); // Data should be unchanged
    assert_eq!(error_pos, Some(1));
}

#[test]
fn test_rs_encode() {
    let data = vec![0x12, 0x34, 0x56, 0x78];
    let n_sym = 4; // 4 error correction symbols (t=2 errors)
    let codeword = rs_encode(&data, n_sym).expect("Should encode");
    
    // Codeword should be data + n_sym bytes
    assert_eq!(codeword.len(), data.len() + n_sym);
    
    // First bytes should be original data
    assert_eq!(&codeword[..4], &data[..]);
}

#[test]
fn test_rs_encode_decode_no_error() {
    let data = vec![0x12, 0x34, 0x56, 0x78, 0x9A];
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
