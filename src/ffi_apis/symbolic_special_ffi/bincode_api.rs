use crate::symbolic::special;
use crate::ffi_apis::common::*;

#[no_mangle]
pub extern "C" fn rssn_bincode_gamma_numerical(val_buf: BincodeBuffer) -> BincodeBuffer {
    let val: Option<f64> = from_bincode_buffer(&val_buf);
    if let Some(v) = val {
        to_bincode_buffer(&special::gamma_numerical(v))
    } else {
        BincodeBuffer::empty()
    }
}

#[no_mangle]
pub extern "C" fn rssn_bincode_ln_gamma_numerical(val_buf: BincodeBuffer) -> BincodeBuffer {
    let val: Option<f64> = from_bincode_buffer(&val_buf);
    if let Some(v) = val {
        to_bincode_buffer(&special::ln_gamma_numerical(v))
    } else {
        BincodeBuffer::empty()
    }
}

#[no_mangle]
pub extern "C" fn rssn_bincode_beta_numerical(a_buf: BincodeBuffer, b_buf: BincodeBuffer) -> BincodeBuffer {
    let a: Option<f64> = from_bincode_buffer(&a_buf);
    let b: Option<f64> = from_bincode_buffer(&b_buf);
    if let (Some(val_a), Some(val_b)) = (a, b) {
        to_bincode_buffer(&special::beta_numerical(val_a, val_b))
    } else {
        BincodeBuffer::empty()
    }
}

#[no_mangle]
pub extern "C" fn rssn_bincode_ln_beta_numerical(a_buf: BincodeBuffer, b_buf: BincodeBuffer) -> BincodeBuffer {
    let a: Option<f64> = from_bincode_buffer(&a_buf);
    let b: Option<f64> = from_bincode_buffer(&b_buf);
    if let (Some(val_a), Some(val_b)) = (a, b) {
        to_bincode_buffer(&special::ln_beta_numerical(val_a, val_b))
    } else {
        BincodeBuffer::empty()
    }
}

#[no_mangle]
pub extern "C" fn rssn_bincode_erf_numerical(val_buf: BincodeBuffer) -> BincodeBuffer {
    let val: Option<f64> = from_bincode_buffer(&val_buf);
    if let Some(v) = val {
        to_bincode_buffer(&special::erf_numerical(v))
    } else {
        BincodeBuffer::empty()
    }
}

#[no_mangle]
pub extern "C" fn rssn_bincode_erfc_numerical(val_buf: BincodeBuffer) -> BincodeBuffer {
    let val: Option<f64> = from_bincode_buffer(&val_buf);
    if let Some(v) = val {
        to_bincode_buffer(&special::erfc_numerical(v))
    } else {
        BincodeBuffer::empty()
    }
}
