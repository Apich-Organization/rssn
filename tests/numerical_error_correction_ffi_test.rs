//! FFI Tests for numerical error correction module.
//!
//! This module contains tests for the Handle-based, JSON, and Bincode FFI APIs
//! for the numerical error correction functions.

use std::ffi::CString;

// ============================================================================
// JSON FFI Tests
// ============================================================================

#[test]

fn test_rs_encode_json() {

    let input = r#"{"message": [1, 2, 3, 4], "n_parity": 4}"#;

    let c_input =
        CString::new(input).unwrap();

    unsafe {

        let result = rssn::ffi_apis::numerical_error_correction_ffi::json::rssn_num_error_correction_rs_encode_json(c_input.as_ptr());

        assert!(!result.is_null());

        let result_str =
            std::ffi::CStr::from_ptr(
                result,
            )
            .to_string_lossy();

        assert!(
            result_str
                .contains("\"ok\":")
        );

        // Parse the result
        let parsed: serde_json::Value =
            serde_json::from_str(
                &result_str,
            )
            .unwrap();

        assert!(
            parsed["ok"].is_array()
        );

        assert_eq!(
            parsed["ok"]
                .as_array()
                .unwrap()
                .len(),
            8
        ); // 4 data + 4 parity
    }
}

#[test]

fn test_rs_check_json_valid() {

    // First encode
    let encode_input = r#"{"message": [1, 2, 3, 4], "n_parity": 4}"#;

    let c_encode_input =
        CString::new(encode_input)
            .unwrap();

    unsafe {

        let encode_result = rssn::ffi_apis::numerical_error_correction_ffi::json::rssn_num_error_correction_rs_encode_json(c_encode_input.as_ptr());

        let encode_str =
            std::ffi::CStr::from_ptr(
                encode_result,
            )
            .to_string_lossy();

        let encoded: serde_json::Value =
            serde_json::from_str(
                &encode_str,
            )
            .unwrap();

        let codeword = encoded["ok"]
            .as_array()
            .unwrap();

        // Now check - construct proper JSON using serde_json
        let check_obj = serde_json::json!({
            "codeword": codeword,
            "n_parity": 4
        });

        let check_input =
            serde_json::to_string(
                &check_obj,
            )
            .unwrap();

        let c_check_input =
            CString::new(check_input)
                .unwrap();

        let check_result = rssn::ffi_apis::numerical_error_correction_ffi::json::rssn_num_error_correction_rs_check_json(c_check_input.as_ptr());

        let check_str =
            std::ffi::CStr::from_ptr(
                check_result,
            )
            .to_string_lossy();

        let checked: serde_json::Value =
            serde_json::from_str(
                &check_str,
            )
            .unwrap();

        assert_eq!(checked["ok"], true);
    }
}

#[test]

fn test_hamming_encode_json() {

    let input =
        r#"{"data": [1, 0, 1, 1]}"#;

    let c_input =
        CString::new(input).unwrap();

    unsafe {

        let result = rssn::ffi_apis::numerical_error_correction_ffi::json::rssn_num_error_correction_hamming_encode_json(c_input.as_ptr());

        assert!(!result.is_null());

        let result_str =
            std::ffi::CStr::from_ptr(
                result,
            )
            .to_string_lossy();

        let parsed: serde_json::Value =
            serde_json::from_str(
                &result_str,
            )
            .unwrap();

        assert!(
            parsed["ok"].is_array()
        );

        assert_eq!(
            parsed["ok"]
                .as_array()
                .unwrap()
                .len(),
            7
        ); // Hamming(7,4)
    }
}

#[test]

fn test_hamming_decode_json() {

    // First encode
    let encode_input =
        r#"{"data": [1, 0, 1, 1]}"#;

    let c_encode_input =
        CString::new(encode_input)
            .unwrap();

    unsafe {

        let encode_result = rssn::ffi_apis::numerical_error_correction_ffi::json::rssn_num_error_correction_hamming_encode_json(c_encode_input.as_ptr());

        let encode_str =
            std::ffi::CStr::from_ptr(
                encode_result,
            )
            .to_string_lossy();

        let encoded: serde_json::Value =
            serde_json::from_str(
                &encode_str,
            )
            .unwrap();

        let codeword = encoded["ok"]
            .as_array()
            .unwrap();

        // Now decode - construct proper JSON using serde_json
        let decode_obj = serde_json::json!({
            "data": codeword
        });

        let decode_input =
            serde_json::to_string(
                &decode_obj,
            )
            .unwrap();

        let c_decode_input =
            CString::new(decode_input)
                .unwrap();

        let decode_result = rssn::ffi_apis::numerical_error_correction_ffi::json::rssn_num_error_correction_hamming_decode_json(c_decode_input.as_ptr());

        let decode_str =
            std::ffi::CStr::from_ptr(
                decode_result,
            )
            .to_string_lossy();

        let decoded: serde_json::Value =
            serde_json::from_str(
                &decode_str,
            )
            .unwrap();

        // Check decoded data
        assert!(
            decoded["ok"]["data"]
                .is_array()
        );

        assert_eq!(
            decoded["ok"]["data"],
            serde_json::json!([
                1, 0, 1, 1
            ])
        );

        assert!(
            decoded["ok"]["error_pos"]
                .is_null()
        );
    }
}

#[test]

fn test_hamming_check_json() {

    // First encode
    let encode_input =
        r#"{"data": [1, 0, 1, 1]}"#;

    let c_encode_input =
        CString::new(encode_input)
            .unwrap();

    unsafe {

        let encode_result = rssn::ffi_apis::numerical_error_correction_ffi::json::rssn_num_error_correction_hamming_encode_json(c_encode_input.as_ptr());

        let encode_str =
            std::ffi::CStr::from_ptr(
                encode_result,
            )
            .to_string_lossy();

        let encoded: serde_json::Value =
            serde_json::from_str(
                &encode_str,
            )
            .unwrap();

        let codeword = encoded["ok"]
            .as_array()
            .unwrap();

        // Now check - construct proper JSON using serde_json
        let check_obj = serde_json::json!({
            "data": codeword
        });

        let check_input =
            serde_json::to_string(
                &check_obj,
            )
            .unwrap();

        let c_check_input =
            CString::new(check_input)
                .unwrap();

        let check_result = rssn::ffi_apis::numerical_error_correction_ffi::json::rssn_num_error_correction_hamming_check_json(c_check_input.as_ptr());

        let check_str =
            std::ffi::CStr::from_ptr(
                check_result,
            )
            .to_string_lossy();

        let checked: serde_json::Value =
            serde_json::from_str(
                &check_str,
            )
            .unwrap();

        assert_eq!(checked["ok"], true);
    }
}

#[test]

fn test_hamming_distance_json() {

    let input = r#"{"a": [1, 0, 1, 0], "b": [0, 0, 1, 1]}"#;

    let c_input =
        CString::new(input).unwrap();

    unsafe {

        let result = rssn::ffi_apis::numerical_error_correction_ffi::json::rssn_num_error_correction_hamming_distance_json(c_input.as_ptr());

        assert!(!result.is_null());

        let result_str =
            std::ffi::CStr::from_ptr(
                result,
            )
            .to_string_lossy();

        let parsed: serde_json::Value =
            serde_json::from_str(
                &result_str,
            )
            .unwrap();

        assert_eq!(parsed["ok"], 2);
    }
}

#[test]

fn test_hamming_weight_json() {

    let input = r#"{"data": [1, 0, 1, 1, 0, 1]}"#;

    let c_input =
        CString::new(input).unwrap();

    unsafe {

        let result = rssn::ffi_apis::numerical_error_correction_ffi::json::rssn_num_error_correction_hamming_weight_json(c_input.as_ptr());

        assert!(!result.is_null());

        let result_str =
            std::ffi::CStr::from_ptr(
                result,
            )
            .to_string_lossy();

        let parsed: serde_json::Value =
            serde_json::from_str(
                &result_str,
            )
            .unwrap();

        assert_eq!(parsed["ok"], 4);
    }
}

#[test]

fn test_crc32_json() {

    // "Hello, World!" in bytes
    let input = r#"{"data": [72, 101, 108, 108, 111, 44, 32, 87, 111, 114, 108, 100, 33]}"#;

    let c_input =
        CString::new(input).unwrap();

    unsafe {

        let result = rssn::ffi_apis::numerical_error_correction_ffi::json::rssn_num_error_correction_crc32_json(c_input.as_ptr());

        assert!(!result.is_null());

        let result_str =
            std::ffi::CStr::from_ptr(
                result,
            )
            .to_string_lossy();

        let parsed: serde_json::Value =
            serde_json::from_str(
                &result_str,
            )
            .unwrap();

        assert_eq!(
            parsed["ok"],
            0xEC4AC3D0_u32
        );
    }
}

#[test]

fn test_crc32_verify_json() {

    // Use proper JSON construction with serde_json
    let data: Vec<u8> =
        b"Hello, World!".to_vec();

    let expected_crc: u32 = 0xEC4AC3D0;

    let input_obj = serde_json::json!({
        "data": data,
        "expected_crc": expected_crc
    });

    let input = serde_json::to_string(
        &input_obj,
    )
    .unwrap();

    let c_input =
        CString::new(input).unwrap();

    unsafe {

        let result = rssn::ffi_apis::numerical_error_correction_ffi::json::rssn_num_error_correction_crc32_verify_json(c_input.as_ptr());

        assert!(!result.is_null());

        let result_str =
            std::ffi::CStr::from_ptr(
                result,
            )
            .to_string_lossy();

        let parsed: serde_json::Value =
            serde_json::from_str(
                &result_str,
            )
            .unwrap();

        assert_eq!(parsed["ok"], true);
    }
}

#[test]

fn test_crc16_json() {

    // "123456789" in bytes
    let input = r#"{"data": [49, 50, 51, 52, 53, 54, 55, 56, 57]}"#;

    let c_input =
        CString::new(input).unwrap();

    unsafe {

        let result = rssn::ffi_apis::numerical_error_correction_ffi::json::rssn_num_error_correction_crc16_json(c_input.as_ptr());

        assert!(!result.is_null());

        let result_str =
            std::ffi::CStr::from_ptr(
                result,
            )
            .to_string_lossy();

        let parsed: serde_json::Value =
            serde_json::from_str(
                &result_str,
            )
            .unwrap();

        // CRC-16 (IBM/Modbus with reflected polynomial 0xA001)
        assert_eq!(
            parsed["ok"],
            0x4B37_u16
        );
    }
}

#[test]

fn test_crc8_json() {

    let input =
        r#"{"data": [1, 2, 3, 4]}"#;

    let c_input =
        CString::new(input).unwrap();

    unsafe {

        let result = rssn::ffi_apis::numerical_error_correction_ffi::json::rssn_num_error_correction_crc8_json(c_input.as_ptr());

        assert!(!result.is_null());

        let result_str =
            std::ffi::CStr::from_ptr(
                result,
            )
            .to_string_lossy();

        let parsed: serde_json::Value =
            serde_json::from_str(
                &result_str,
            )
            .unwrap();

        // Just verify it returns a number
        assert!(
            parsed["ok"].is_number()
        );
    }
}

#[test]

fn test_interleave_json() {

    let input = r#"{"data": [1, 2, 3, 4, 5, 6], "depth": 3}"#;

    let c_input =
        CString::new(input).unwrap();

    unsafe {

        let result = rssn::ffi_apis::numerical_error_correction_ffi::json::rssn_num_error_correction_interleave_json(c_input.as_ptr());

        assert!(!result.is_null());

        let result_str =
            std::ffi::CStr::from_ptr(
                result,
            )
            .to_string_lossy();

        let parsed: serde_json::Value =
            serde_json::from_str(
                &result_str,
            )
            .unwrap();

        assert!(
            parsed["ok"].is_array()
        );

        assert_eq!(
            parsed["ok"]
                .as_array()
                .unwrap()
                .len(),
            6
        );
    }
}

#[test]

fn test_deinterleave_json() {

    // First interleave
    let interleave_input = r#"{"data": [1, 2, 3, 4, 5, 6], "depth": 3}"#;

    let c_interleave_input =
        CString::new(interleave_input)
            .unwrap();

    unsafe {

        let interleave_result = rssn::ffi_apis::numerical_error_correction_ffi::json::rssn_num_error_correction_interleave_json(c_interleave_input.as_ptr());

        let interleave_str =
            std::ffi::CStr::from_ptr(
                interleave_result,
            )
            .to_string_lossy();

        let interleaved : serde_json::Value = serde_json::from_str(&interleave_str).unwrap();

        let interleaved_data =
            interleaved["ok"]
                .as_array()
                .unwrap();

        // Now deinterleave - construct proper JSON using serde_json
        let deinterleave_obj = serde_json::json!({
            "data": interleaved_data,
            "depth": 3
        });

        let deinterleave_input =
            serde_json::to_string(
                &deinterleave_obj,
            )
            .unwrap();

        let c_deinterleave_input =
            CString::new(
                deinterleave_input,
            )
            .unwrap();

        let deinterleave_result = rssn::ffi_apis::numerical_error_correction_ffi::json::rssn_num_error_correction_deinterleave_json(c_deinterleave_input.as_ptr());

        let deinterleave_str =
            std::ffi::CStr::from_ptr(
                deinterleave_result,
            )
            .to_string_lossy();

        let deinterleaved : serde_json::Value = serde_json::from_str(&deinterleave_str).unwrap();

        // Should get back original data
        assert_eq!(
            deinterleaved["ok"],
            serde_json::json!([
                1, 2, 3, 4, 5, 6
            ])
        );
    }
}

#[test]

fn test_code_rate_json() {

    let input = r#"{"k": 4, "n": 7}"#;

    let c_input =
        CString::new(input).unwrap();

    unsafe {

        let result = rssn::ffi_apis::numerical_error_correction_ffi::json::rssn_num_error_correction_code_rate_json(c_input.as_ptr());

        assert!(!result.is_null());

        let result_str =
            std::ffi::CStr::from_ptr(
                result,
            )
            .to_string_lossy();

        let parsed: serde_json::Value =
            serde_json::from_str(
                &result_str,
            )
            .unwrap();

        let rate = parsed["ok"]
            .as_f64()
            .unwrap();

        assert!(
            (rate - 4.0 / 7.0).abs()
                < 1e-10
        );
    }
}

#[test]

fn test_capability_json() {

    let input =
        r#"{"min_distance": 5}"#;

    let c_input =
        CString::new(input).unwrap();

    unsafe {

        let result = rssn::ffi_apis::numerical_error_correction_ffi::json::rssn_num_error_correction_capability_json(c_input.as_ptr());

        assert!(!result.is_null());

        let result_str =
            std::ffi::CStr::from_ptr(
                result,
            )
            .to_string_lossy();

        let parsed: serde_json::Value =
            serde_json::from_str(
                &result_str,
            )
            .unwrap();

        assert_eq!(parsed["ok"], 2); // (5-1)/2 = 2
    }
}

// ============================================================================
// Handle-based FFI Tests
// ============================================================================

#[test]

fn test_crc32_handle() {

    let data = b"Hello, World!";

    unsafe {

        let crc =
            rssn::ffi_apis::numerical_error_correction_ffi::handle::rssn_num_error_correction_crc32(
                data.as_ptr(),
                data.len(),
            );

        assert_eq!(crc, 0xEC4AC3D0);
    }
}

#[test]

fn test_crc32_verify_handle() {

    let data = b"Hello, World!";

    unsafe {

        let result = rssn::ffi_apis::numerical_error_correction_ffi::handle::rssn_num_error_correction_crc32_verify(
            data.as_ptr(),
            data.len(),
            0xEC4AC3D0,
        );

        assert_eq!(result, 1); // true
    }
}

#[test]

fn test_crc16_handle() {

    let data = b"123456789";

    unsafe {

        let crc =
            rssn::ffi_apis::numerical_error_correction_ffi::handle::rssn_num_error_correction_crc16(
                data.as_ptr(),
                data.len(),
            );

        // CRC-16 (IBM/Modbus with reflected polynomial 0xA001)
        assert_eq!(crc, 0x4B37);
    }
}

#[test]

fn test_hamming_encode_handle() {

    let data: [u8; 4] = [1, 0, 1, 1];

    let mut output: [u8; 7] = [0; 7];

    unsafe {

        let result = rssn::ffi_apis::numerical_error_correction_ffi::handle::rssn_num_error_correction_hamming_encode(
            data.as_ptr(),
            output.as_mut_ptr(),
        );

        assert_eq!(result, 0); // Success
    }
}

#[test]

fn test_hamming_decode_handle() {

    let data: [u8; 4] = [1, 0, 1, 1];

    let mut codeword: [u8; 7] = [0; 7];

    let mut decoded: [u8; 4] = [0; 4];

    let mut error_pos: i32 = 0;

    unsafe {

        // Encode
        rssn::ffi_apis::numerical_error_correction_ffi::handle::rssn_num_error_correction_hamming_encode(
            data.as_ptr(),
            codeword.as_mut_ptr(),
        );

        // Decode
        let result = rssn::ffi_apis::numerical_error_correction_ffi::handle::rssn_num_error_correction_hamming_decode(
            codeword.as_ptr(),
            decoded.as_mut_ptr(),
            &mut error_pos,
        );

        assert_eq!(result, 0); // Success
        assert_eq!(decoded, data);

        assert_eq!(error_pos, -1); // No error
    }
}

#[test]

fn test_hamming_distance_handle() {

    let a: [u8; 4] = [1, 0, 1, 0];

    let b: [u8; 4] = [0, 0, 1, 1];

    unsafe {

        let dist = rssn::ffi_apis::numerical_error_correction_ffi::handle::rssn_num_error_correction_hamming_distance(
            a.as_ptr(),
            b.as_ptr(),
            4,
        );

        assert_eq!(dist, 2);
    }
}

#[test]

fn test_hamming_weight_handle() {

    let data: [u8; 6] =
        [1, 0, 1, 1, 0, 1];

    unsafe {

        let weight = rssn::ffi_apis::numerical_error_correction_ffi::handle::rssn_num_error_correction_hamming_weight(
            data.as_ptr(),
            6,
        );

        assert_eq!(weight, 4);
    }
}

#[test]

fn test_code_rate_handle() {

    let rate =
        rssn::ffi_apis::numerical_error_correction_ffi::handle::rssn_num_error_correction_code_rate(
            4, 7,
        );

    assert!(
        (rate - 4.0 / 7.0).abs()
            < 1e-10
    );
}

#[test]

fn test_capability_handle() {

    let capability = rssn::ffi_apis::numerical_error_correction_ffi::handle::rssn_num_error_correction_capability(5);

    assert_eq!(capability, 2);
}

#[test]

fn test_detection_capability_handle() {

    let capability =
        rssn::ffi_apis::numerical_error_correction_ffi::handle::rssn_num_error_detection_capability(
            5,
        );

    assert_eq!(capability, 4);
}

#[test]

fn test_interleave_handle() {

    let data: [u8; 6] =
        [1, 2, 3, 4, 5, 6];

    let mut output: [u8; 6] = [0; 6];

    unsafe {

        let result = rssn::ffi_apis::numerical_error_correction_ffi::handle::rssn_num_error_correction_interleave(
            data.as_ptr(),
            6,
            3,
            output.as_mut_ptr(),
        );

        assert_eq!(result, 0);
    }
}

#[test]

fn test_rs_encode_handle() {

    let message: [u8; 4] = [1, 2, 3, 4];

    let mut output: [u8; 8] = [0; 8]; // 4 data + 4 parity
    let mut out_len: usize = 8;

    unsafe {

        let result = rssn::ffi_apis::numerical_error_correction_ffi::handle::rssn_num_error_correction_rs_encode(
            message.as_ptr(),
            4,
            4,
            output.as_mut_ptr(),
            &mut out_len,
        );

        assert_eq!(result, 0);

        assert_eq!(out_len, 8);

        // First 4 bytes should be the message
        assert_eq!(
            &output[.. 4],
            &message
        );
    }
}

#[test]

fn test_rs_check_handle() {

    let message: [u8; 4] = [1, 2, 3, 4];

    let mut codeword: [u8; 8] = [0; 8];

    let mut out_len: usize = 8;

    unsafe {

        rssn::ffi_apis::numerical_error_correction_ffi::handle::rssn_num_error_correction_rs_encode(
            message.as_ptr(),
            4,
            4,
            codeword.as_mut_ptr(),
            &mut out_len,
        );

        let result = rssn::ffi_apis::numerical_error_correction_ffi::handle::rssn_num_error_correction_rs_check(
            codeword.as_ptr(),
            8,
            4,
        );

        assert_eq!(result, 1); // Valid codeword
    }
}
