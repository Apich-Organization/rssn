//! Bincode-based FFI API for numerical error correction functions.

use crate::ffi_apis::common::{from_bincode_buffer, to_bincode_buffer, BincodeBuffer};
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::error_correction;
use serde::{Deserialize, Serialize};

#[derive(Deserialize)]

struct RsEncodeInput {
    message: Vec<u8>,
    n_parity: usize,
}

#[derive(Deserialize)]

struct RsDecodeInput {
    codeword: Vec<u8>,
    n_parity: usize,
}

#[derive(Deserialize)]

struct HammingInput {
    data: Vec<u8>,
}

#[derive(Deserialize)]

struct DistanceInput {
    a: Vec<u8>,
    b: Vec<u8>,
}

#[derive(Deserialize)]

struct CrcInput {
    data: Vec<u8>,
}

#[derive(Deserialize)]

struct Crc32VerifyInput {
    data: Vec<u8>,
    expected_crc: u32,
}

#[derive(Deserialize)]

struct InterleaveInput {
    data: Vec<u8>,
    depth: usize,
}

#[derive(Deserialize)]

struct CodeRateInput {
    k: usize,
    n: usize,
}

#[derive(Deserialize)]

struct CapabilityInput {
    min_distance: usize,
}

#[derive(Serialize)]

struct HammingDecodeResult {
    data: Vec<u8>,
    error_pos: Option<usize>,
}

// Reed-Solomon functions
#[no_mangle]

pub unsafe extern "C" fn rssn_num_error_correction_rs_encode_bincode(
    buffer: BincodeBuffer,
) -> BincodeBuffer {

    let input: RsEncodeInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<Vec<u8>, String> {
                ok: None,
                err: Some("Invalid Bincode".to_string()),
            })
        }
    };

    match error_correction::reed_solomon_encode(&input.message, input.n_parity) {
        Ok(codeword) => to_bincode_buffer(&FfiResult {
            ok: Some(codeword),
            err: None::<String>,
        }),
        Err(e) => to_bincode_buffer(&FfiResult::<Vec<u8>, String> {
            ok: None,
            err: Some(e),
        }),
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_error_correction_rs_decode_bincode(
    buffer: BincodeBuffer,
) -> BincodeBuffer {

    let input: RsDecodeInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<Vec<u8>, String> {
                ok: None,
                err: Some("Invalid Bincode".to_string()),
            })
        }
    };

    let mut codeword = input.codeword;

    match error_correction::reed_solomon_decode(&mut codeword, input.n_parity) {
        Ok(()) => to_bincode_buffer(&FfiResult {
            ok: Some(codeword),
            err: None::<String>,
        }),
        Err(e) => to_bincode_buffer(&FfiResult::<Vec<u8>, String> {
            ok: None,
            err: Some(e),
        }),
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_error_correction_rs_check_bincode(
    buffer: BincodeBuffer,
) -> BincodeBuffer {

    let input: RsDecodeInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<bool, String> {
                ok: None,
                err: Some("Invalid Bincode".to_string()),
            })
        }
    };

    let result = error_correction::reed_solomon_check(&input.codeword, input.n_parity);

    to_bincode_buffer(&FfiResult {
        ok: Some(result),
        err: None::<String>,
    })
}

// Hamming functions
#[no_mangle]

pub unsafe extern "C" fn rssn_num_error_correction_hamming_encode_bincode(
    buffer: BincodeBuffer,
) -> BincodeBuffer {

    let input: HammingInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<Vec<u8>, String> {
                ok: None,
                err: Some("Invalid Bincode".to_string()),
            })
        }
    };

    match error_correction::hamming_encode_numerical(&input.data) {
        Some(codeword) => to_bincode_buffer(&FfiResult {
            ok: Some(codeword),
            err: None::<String>,
        }),
        None => to_bincode_buffer(&FfiResult::<Vec<u8>, String> {
            ok: None,
            err: Some("Input must be exactly 4 bytes".to_string()),
        }),
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_error_correction_hamming_decode_bincode(
    buffer: BincodeBuffer,
) -> BincodeBuffer {

    let input: HammingInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<HammingDecodeResult, String> {
                ok: None,
                err: Some("Invalid Bincode".to_string()),
            })
        }
    };

    match error_correction::hamming_decode_numerical(&input.data) {
        Ok((data, error_pos)) => to_bincode_buffer(&FfiResult {
            ok: Some(HammingDecodeResult { data, error_pos }),
            err: None::<String>,
        }),
        Err(e) => to_bincode_buffer(&FfiResult::<HammingDecodeResult, String> {
            ok: None,
            err: Some(e),
        }),
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_error_correction_hamming_check_bincode(
    buffer: BincodeBuffer,
) -> BincodeBuffer {

    let input: HammingInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<bool, String> {
                ok: None,
                err: Some("Invalid Bincode".to_string()),
            })
        }
    };

    let result = error_correction::hamming_check_numerical(&input.data);

    to_bincode_buffer(&FfiResult {
        ok: Some(result),
        err: None::<String>,
    })
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_error_correction_hamming_distance_bincode(
    buffer: BincodeBuffer,
) -> BincodeBuffer {

    let input: DistanceInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<usize, String> {
                ok: None,
                err: Some("Invalid Bincode".to_string()),
            })
        }
    };

    match error_correction::hamming_distance_numerical(&input.a, &input.b) {
        Some(dist) => to_bincode_buffer(&FfiResult {
            ok: Some(dist),
            err: None::<String>,
        }),
        None => to_bincode_buffer(&FfiResult::<usize, String> {
            ok: None,
            err: Some("Vectors must have same length".to_string()),
        }),
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_error_correction_hamming_weight_bincode(
    buffer: BincodeBuffer,
) -> BincodeBuffer {

    let input: HammingInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<usize, String> {
                ok: None,
                err: Some("Invalid Bincode".to_string()),
            })
        }
    };

    let result = error_correction::hamming_weight_numerical(&input.data);

    to_bincode_buffer(&FfiResult {
        ok: Some(result),
        err: None::<String>,
    })
}

// CRC functions
#[no_mangle]

pub unsafe extern "C" fn rssn_num_error_correction_crc32_bincode(
    buffer: BincodeBuffer,
) -> BincodeBuffer {

    let input: CrcInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<u32, String> {
                ok: None,
                err: Some("Invalid Bincode".to_string()),
            })
        }
    };

    let result = error_correction::crc32_compute_numerical(&input.data);

    to_bincode_buffer(&FfiResult {
        ok: Some(result),
        err: None::<String>,
    })
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_error_correction_crc32_verify_bincode(
    buffer: BincodeBuffer,
) -> BincodeBuffer {

    let input: Crc32VerifyInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<bool, String> {
                ok: None,
                err: Some("Invalid Bincode".to_string()),
            })
        }
    };

    let result = error_correction::crc32_verify_numerical(&input.data, input.expected_crc);

    to_bincode_buffer(&FfiResult {
        ok: Some(result),
        err: None::<String>,
    })
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_error_correction_crc16_bincode(
    buffer: BincodeBuffer,
) -> BincodeBuffer {

    let input: CrcInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<u16, String> {
                ok: None,
                err: Some("Invalid Bincode".to_string()),
            })
        }
    };

    let result = error_correction::crc16_compute(&input.data);

    to_bincode_buffer(&FfiResult {
        ok: Some(result),
        err: None::<String>,
    })
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_error_correction_crc8_bincode(
    buffer: BincodeBuffer,
) -> BincodeBuffer {

    let input: CrcInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<u8, String> {
                ok: None,
                err: Some("Invalid Bincode".to_string()),
            })
        }
    };

    let result = error_correction::crc8_compute(&input.data);

    to_bincode_buffer(&FfiResult {
        ok: Some(result),
        err: None::<String>,
    })
}

// Interleaving functions
#[no_mangle]

pub unsafe extern "C" fn rssn_num_error_correction_interleave_bincode(
    buffer: BincodeBuffer,
) -> BincodeBuffer {

    let input: InterleaveInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<Vec<u8>, String> {
                ok: None,
                err: Some("Invalid Bincode".to_string()),
            })
        }
    };

    let result = error_correction::interleave(&input.data, input.depth);

    to_bincode_buffer(&FfiResult {
        ok: Some(result),
        err: None::<String>,
    })
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_error_correction_deinterleave_bincode(
    buffer: BincodeBuffer,
) -> BincodeBuffer {

    let input: InterleaveInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<Vec<u8>, String> {
                ok: None,
                err: Some("Invalid Bincode".to_string()),
            })
        }
    };

    let result = error_correction::deinterleave(&input.data, input.depth);

    to_bincode_buffer(&FfiResult {
        ok: Some(result),
        err: None::<String>,
    })
}

// Code theory functions
#[no_mangle]

pub unsafe extern "C" fn rssn_num_error_correction_code_rate_bincode(
    buffer: BincodeBuffer,
) -> BincodeBuffer {

    let input: CodeRateInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<f64, String> {
                ok: None,
                err: Some("Invalid Bincode".to_string()),
            })
        }
    };

    let result = error_correction::code_rate(input.k, input.n);

    to_bincode_buffer(&FfiResult {
        ok: Some(result),
        err: None::<String>,
    })
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_error_correction_capability_bincode(
    buffer: BincodeBuffer,
) -> BincodeBuffer {

    let input: CapabilityInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<usize, String> {
                ok: None,
                err: Some("Invalid Bincode".to_string()),
            })
        }
    };

    let result = error_correction::error_correction_capability(input.min_distance);

    to_bincode_buffer(&FfiResult {
        ok: Some(result),
        err: None::<String>,
    })
}
