//! JSON-based FFI API for numerical error correction functions.

use std::os::raw::c_char;

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_c_string;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::error_correction;

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
/// Encodes a message using Reed-Solomon codes with JSON for serialization.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_error_correction_rs_encode_json(
    input: *const c_char
) -> *mut c_char {

    let input : RsEncodeInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<Vec<u8>, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    match error_correction::reed_solomon_encode(
        &input.message,
        input.n_parity,
    ) {
        | Ok(codeword) => {
            to_c_string(
                serde_json::to_string(&FfiResult {
                    ok : Some(codeword),
                    err : None::<String>,
                })
                .unwrap(),
            )
        },
        | Err(e) => {
            to_c_string(
                serde_json::to_string(
                    &FfiResult::<Vec<u8>, String> {
                        ok : None,
                        err : Some(e),
                    },
                )
                .unwrap(),
            )
        },
    }
}

/// Decodes a Reed-Solomon codeword with JSON for serialization.

#[no_mangle]

pub unsafe extern "C" fn rssn_num_error_correction_rs_decode_json(
    input: *const c_char
) -> *mut c_char {

    let input : RsDecodeInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<Vec<u8>, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let mut codeword = input.codeword;

    match error_correction::reed_solomon_decode(
        &mut codeword,
        input.n_parity,
    ) {
        | Ok(()) => {
            to_c_string(
                serde_json::to_string(&FfiResult {
                    ok : Some(codeword),
                    err : None::<String>,
                })
                .unwrap(),
            )
        },
        | Err(e) => {
            to_c_string(
                serde_json::to_string(
                    &FfiResult::<Vec<u8>, String> {
                        ok : None,
                        err : Some(e),
                    },
                )
                .unwrap(),
            )
        },
    }
}

/// Checks if a Reed-Solomon codeword is valid with JSON for serialization.

#[no_mangle]

pub unsafe extern "C" fn rssn_num_error_correction_rs_check_json(
    input: *const c_char
) -> *mut c_char {

    let input : RsDecodeInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<bool, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result = error_correction::reed_solomon_check(
        &input.codeword,
        input.n_parity,
    );

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(result),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

// Hamming functions
/// Encodes a message using Hamming codes with JSON for serialization.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_error_correction_hamming_encode_json(
    input: *const c_char
) -> *mut c_char {

    let input : HammingInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<Vec<u8>, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    match error_correction::hamming_encode_numerical(&input.data) {
        | Some(codeword) => {
            to_c_string(
                serde_json::to_string(&FfiResult {
                    ok : Some(codeword),
                    err : None::<String>,
                })
                .unwrap(),
            )
        },
        | None => {
            to_c_string(
                serde_json::to_string(
                    &FfiResult::<Vec<u8>, String> {
                        ok : None,
                        err : Some("Input must be exactly 4 bytes".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    }
}

/// Decodes a Hamming codeword and corrects errors with JSON for serialization.

#[no_mangle]

pub unsafe extern "C" fn rssn_num_error_correction_hamming_decode_json(
    input: *const c_char
) -> *mut c_char {

    let input : HammingInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<HammingDecodeResult, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    match error_correction::hamming_decode_numerical(&input.data) {
        | Ok((data, error_pos)) => {
            to_c_string(
                serde_json::to_string(&FfiResult {
                    ok : Some(
                        HammingDecodeResult {
                            data,
                            error_pos,
                        },
                    ),
                    err : None::<String>,
                })
                .unwrap(),
            )
        },
        | Err(e) => {
            to_c_string(
                serde_json::to_string(
                    &FfiResult::<HammingDecodeResult, String> {
                        ok : None,
                        err : Some(e),
                    },
                )
                .unwrap(),
            )
        },
    }
}

/// Checks a Hamming codeword for errors with JSON for serialization.

#[no_mangle]

pub unsafe extern "C" fn rssn_num_error_correction_hamming_check_json(
    input: *const c_char
) -> *mut c_char {

    let input : HammingInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<bool, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result = error_correction::hamming_check_numerical(&input.data);

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(result),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

/// Computes the Hamming distance between two byte vectors with JSON for serialization.

#[no_mangle]

pub unsafe extern "C" fn rssn_num_error_correction_hamming_distance_json(
    input: *const c_char
) -> *mut c_char {

    let input : DistanceInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<usize, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    match error_correction::hamming_distance_numerical(&input.a, &input.b) {
        | Some(dist) => {
            to_c_string(
                serde_json::to_string(&FfiResult {
                    ok : Some(dist),
                    err : None::<String>,
                })
                .unwrap(),
            )
        },
        | None => {
            to_c_string(
                serde_json::to_string(
                    &FfiResult::<usize, String> {
                        ok : None,
                        err : Some("Vectors must have same length".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    }
}

/// Computes the Hamming weight of a byte vector with JSON for serialization.

#[no_mangle]

pub unsafe extern "C" fn rssn_num_error_correction_hamming_weight_json(
    input: *const c_char
) -> *mut c_char {

    let input : HammingInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<usize, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result = error_correction::hamming_weight_numerical(&input.data);

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(result),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

// CRC functions
/// Computes the CRC32 checksum of a byte vector with JSON for serialization.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_error_correction_crc32_json(
    input: *const c_char
) -> *mut c_char {

    let input : CrcInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<u32, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result = error_correction::crc32_compute_numerical(&input.data);

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(result),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

/// Verifies the CRC32 checksum of a byte vector with JSON for serialization.

#[no_mangle]

pub unsafe extern "C" fn rssn_num_error_correction_crc32_verify_json(
    input: *const c_char
) -> *mut c_char {

    let input : Crc32VerifyInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<bool, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result = error_correction::crc32_verify_numerical(
        &input.data,
        input.expected_crc,
    );

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(result),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

/// Computes the CRC16 checksum of a byte vector with JSON for serialization.

#[no_mangle]

pub unsafe extern "C" fn rssn_num_error_correction_crc16_json(
    input: *const c_char
) -> *mut c_char {

    let input : CrcInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<u16, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result =
        error_correction::crc16_compute(
            &input.data,
        );

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(result),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

/// Computes the CRC8 checksum of a byte vector with JSON for serialization.

#[no_mangle]

pub unsafe extern "C" fn rssn_num_error_correction_crc8_json(
    input: *const c_char
) -> *mut c_char {

    let input : CrcInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<u8, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result =
        error_correction::crc8_compute(
            &input.data,
        );

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(result),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

// Interleaving functions
/// Interleaves a byte vector with a given depth using JSON for serialization.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_error_correction_interleave_json(
    input: *const c_char
) -> *mut c_char {

    let input : InterleaveInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<Vec<u8>, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result =
        error_correction::interleave(
            &input.data,
            input.depth,
        );

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(result),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

/// Deinterleaves a byte vector with a given depth using JSON for serialization.

#[no_mangle]

pub unsafe extern "C" fn rssn_num_error_correction_deinterleave_json(
    input: *const c_char
) -> *mut c_char {

    let input : InterleaveInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<Vec<u8>, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result =
        error_correction::deinterleave(
            &input.data,
            input.depth,
        );

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(result),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

// Code theory functions
/// Computes the code rate (k/n) of a code with JSON for serialization.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_error_correction_code_rate_json(
    input: *const c_char
) -> *mut c_char {

    let input : CodeRateInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result =
        error_correction::code_rate(
            input.k,
            input.n,
        );

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(result),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

/// Computes the error correction capability of a code given its minimum distance, using JSON for serialization.

#[no_mangle]

pub unsafe extern "C" fn rssn_num_error_correction_capability_json(
    input: *const c_char
) -> *mut c_char {

    let input : CapabilityInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<usize, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result = error_correction::error_correction_capability(input.min_distance);

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(result),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}
