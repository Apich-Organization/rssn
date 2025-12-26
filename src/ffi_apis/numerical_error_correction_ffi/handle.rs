//! Handle-based FFI API for numerical error correction functions.

use crate::numerical::error_correction;

/// Reed-Solomon encode a message.
///
/// # Safety
/// `message_ptr` must be a valid pointer to `message_len` bytes.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_error_correction_rs_encode(
    message_ptr: *const u8,
    message_len: usize,
    n_parity: usize,
    out_ptr: *mut u8,
    out_len: *mut usize,
) -> i32 {

    if message_ptr.is_null()
        || out_ptr.is_null()
        || out_len.is_null()
    {

        return -1;
    }

    let message =
        std::slice::from_raw_parts(
            message_ptr,
            message_len,
        );

    match error_correction::reed_solomon_encode(message, n_parity) {
        | Ok(codeword) => {

            let copy_len = codeword
                .len()
                .min(*out_len);

            std::ptr::copy_nonoverlapping(
                codeword.as_ptr(),
                out_ptr,
                copy_len,
            );

            *out_len = codeword.len();

            0
        },
        | Err(_) => -2,
    }
}

/// Reed-Solomon decode a codeword in place.
///
/// # Safety
/// `codeword_ptr` must be a valid pointer to `codeword_len` bytes.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_error_correction_rs_decode(
    codeword_ptr: *mut u8,
    codeword_len: usize,
    n_parity: usize,
) -> i32 {

    if codeword_ptr.is_null() {

        return -1;
    }

    let codeword =
        std::slice::from_raw_parts_mut(
            codeword_ptr,
            codeword_len,
        );

    match error_correction::reed_solomon_decode(codeword, n_parity) {
        | Ok(()) => 0,
        | Err(_) => -2,
    }
}

/// Check if a Reed-Solomon codeword is valid.
///
/// # Safety
/// `codeword_ptr` must be a valid pointer to `codeword_len` bytes.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_error_correction_rs_check(
    codeword_ptr: *const u8,
    codeword_len: usize,
    n_parity: usize,
) -> i32 {

    if codeword_ptr.is_null() {

        return -1;
    }

    let codeword =
        std::slice::from_raw_parts(
            codeword_ptr,
            codeword_len,
        );

    if error_correction::reed_solomon_check(codeword, n_parity) {

        1
    } else {

        0
    }
}

/// Hamming encode a 4-bit data block.
///
/// # Safety
/// `data_ptr` must be a valid pointer to 4 bytes.
/// `out_ptr` must be a valid pointer to at least 7 bytes.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_error_correction_hamming_encode(
    data_ptr: *const u8,
    out_ptr: *mut u8,
) -> i32 {

    if data_ptr.is_null()
        || out_ptr.is_null()
    {

        return -1;
    }

    let data =
        std::slice::from_raw_parts(
            data_ptr,
            4,
        );

    match error_correction::hamming_encode_numerical(data) {
        | Some(codeword) => {

            std::ptr::copy_nonoverlapping(
                codeword.as_ptr(),
                out_ptr,
                7,
            );

            0
        },
        | None => -2,
    }
}

/// Hamming decode a 7-bit codeword.
///
/// # Safety
/// `codeword_ptr` must be a valid pointer to 7 bytes.
/// `out_ptr` must be a valid pointer to at least 4 bytes.
/// `error_pos_ptr` must be a valid pointer.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_error_correction_hamming_decode(
    codeword_ptr: *const u8,
    out_ptr: *mut u8,
    error_pos_ptr: *mut i32,
) -> i32 {

    if codeword_ptr.is_null()
        || out_ptr.is_null()
        || error_pos_ptr.is_null()
    {

        return -1;
    }

    let codeword =
        std::slice::from_raw_parts(
            codeword_ptr,
            7,
        );

    match error_correction::hamming_decode_numerical(codeword) {
        | Ok((data, error_pos)) => {

            std::ptr::copy_nonoverlapping(
                data.as_ptr(),
                out_ptr,
                4,
            );

            *error_pos_ptr = error_pos.map_or(-1, |p| p as i32);

            0
        },
        | Err(_) => -2,
    }
}

/// Check if a Hamming codeword is valid.
///
/// # Safety
/// `codeword_ptr` must be a valid pointer to 7 bytes.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_error_correction_hamming_check(
    codeword_ptr: *const u8
) -> i32 {

    if codeword_ptr.is_null() {

        return -1;
    }

    let codeword =
        std::slice::from_raw_parts(
            codeword_ptr,
            7,
        );

    if error_correction::hamming_check_numerical(codeword) {

        1
    } else {

        0
    }
}

/// Compute Hamming distance between two byte arrays.
///
/// # Safety
/// `a_ptr` and `b_ptr` must be valid pointers to `len` bytes each.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_error_correction_hamming_distance(
    a_ptr: *const u8,
    b_ptr: *const u8,
    len: usize,
) -> i32 {

    if a_ptr.is_null()
        || b_ptr.is_null()
    {

        return -1;
    }

    let a = std::slice::from_raw_parts(
        a_ptr, len,
    );

    let b = std::slice::from_raw_parts(
        b_ptr, len,
    );

    error_correction::hamming_distance_numerical(a, b).map_or(-1, |d| d as i32)
}

/// Compute Hamming weight of a byte array.
///
/// # Safety
/// `data_ptr` must be a valid pointer to `len` bytes.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_error_correction_hamming_weight(
    data_ptr: *const u8,
    len: usize,
) -> i32 {

    if data_ptr.is_null() {

        return -1;
    }

    let data =
        std::slice::from_raw_parts(
            data_ptr,
            len,
        );

    error_correction::hamming_weight_numerical(data) as i32
}

/// Compute CRC-32 checksum.
///
/// # Safety
/// `data_ptr` must be a valid pointer to `len` bytes.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_error_correction_crc32(
    data_ptr: *const u8,
    len: usize,
) -> u32 {

    if data_ptr.is_null() {

        return 0;
    }

    let data =
        std::slice::from_raw_parts(
            data_ptr,
            len,
        );

    error_correction::crc32_compute_numerical(data)
}

/// Verify CRC-32 checksum.
///
/// # Safety
/// `data_ptr` must be a valid pointer to `len` bytes.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_error_correction_crc32_verify(
    data_ptr: *const u8,
    len: usize,
    expected_crc: u32,
) -> i32 {

    if data_ptr.is_null() {

        return -1;
    }

    let data =
        std::slice::from_raw_parts(
            data_ptr,
            len,
        );

    if error_correction::crc32_verify_numerical(data, expected_crc) {

        1
    } else {

        0
    }
}

/// Compute CRC-16 checksum.
///
/// # Safety
/// `data_ptr` must be a valid pointer to `len` bytes.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_error_correction_crc16(
    data_ptr: *const u8,
    len: usize,
) -> u16 {

    if data_ptr.is_null() {

        return 0;
    }

    let data =
        std::slice::from_raw_parts(
            data_ptr,
            len,
        );

    error_correction::crc16_compute(
        data,
    )
}

/// Compute CRC-8 checksum.
///
/// # Safety
/// `data_ptr` must be a valid pointer to `len` bytes.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_error_correction_crc8(
    data_ptr: *const u8,
    len: usize,
) -> u8 {

    if data_ptr.is_null() {

        return 0;
    }

    let data =
        std::slice::from_raw_parts(
            data_ptr,
            len,
        );

    error_correction::crc8_compute(data)
}

/// Interleave data.
///
/// # Safety
/// `data_ptr` must be a valid pointer to `len` bytes.
/// `out_ptr` must be a valid pointer to at least `len` bytes.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_error_correction_interleave(
    data_ptr: *const u8,
    len: usize,
    depth: usize,
    out_ptr: *mut u8,
) -> i32 {

    if data_ptr.is_null()
        || out_ptr.is_null()
    {

        return -1;
    }

    let data =
        std::slice::from_raw_parts(
            data_ptr,
            len,
        );

    let result =
        error_correction::interleave(
            data, depth,
        );

    std::ptr::copy_nonoverlapping(
        result.as_ptr(),
        out_ptr,
        result.len(),
    );

    0
}

/// De-interleave data.
///
/// # Safety
/// `data_ptr` must be a valid pointer to `len` bytes.
/// `out_ptr` must be a valid pointer to at least `len` bytes.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_error_correction_deinterleave(
    data_ptr: *const u8,
    len: usize,
    depth: usize,
    out_ptr: *mut u8,
) -> i32 {

    if data_ptr.is_null()
        || out_ptr.is_null()
    {

        return -1;
    }

    let data =
        std::slice::from_raw_parts(
            data_ptr,
            len,
        );

    let result =
        error_correction::deinterleave(
            data, depth,
        );

    std::ptr::copy_nonoverlapping(
        result.as_ptr(),
        out_ptr,
        result.len(),
    );

    0
}

/// Compute code rate.
#[no_mangle]

pub extern "C" fn rssn_num_error_correction_code_rate(
    k: usize,
    n: usize,
) -> f64 {

    error_correction::code_rate(k, n)
}

/// Compute error correction capability from minimum distance.
#[no_mangle]

pub extern "C" fn rssn_num_error_correction_capability(
    min_distance: usize
) -> usize {

    error_correction::error_correction_capability(min_distance)
}

/// Compute error detection capability from minimum distance.
#[no_mangle]

pub extern "C" fn rssn_num_error_detection_capability(
    min_distance: usize
) -> usize {

    error_correction::error_detection_capability(min_distance)
}
