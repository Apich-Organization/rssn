use crate::ffi_apis::common::*;
use crate::symbolic::vector::*;

#[no_mangle]

pub extern "C" fn rssn_bincode_vector_magnitude(v_buf: BincodeBuffer) -> BincodeBuffer {

    let v: Option<Vector> = from_bincode_buffer(&v_buf);

    if let Some(vector) = v {

        let result = vector.magnitude();

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

#[no_mangle]

pub extern "C" fn rssn_bincode_vector_dot(
    v1_buf: BincodeBuffer,
    v2_buf: BincodeBuffer,
) -> BincodeBuffer {

    let v1: Option<Vector> = from_bincode_buffer(&v1_buf);

    let v2: Option<Vector> = from_bincode_buffer(&v2_buf);

    if let (Some(vec1), Some(vec2)) = (v1, v2) {

        let result = vec1.dot(&vec2);

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

#[no_mangle]

pub extern "C" fn rssn_bincode_vector_cross(
    v1_buf: BincodeBuffer,
    v2_buf: BincodeBuffer,
) -> BincodeBuffer {

    let v1: Option<Vector> = from_bincode_buffer(&v1_buf);

    let v2: Option<Vector> = from_bincode_buffer(&v2_buf);

    if let (Some(vec1), Some(vec2)) = (v1, v2) {

        let result = vec1.cross(&vec2);

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

#[no_mangle]

pub extern "C" fn rssn_bincode_vector_normalize(v_buf: BincodeBuffer) -> BincodeBuffer {

    let v: Option<Vector> = from_bincode_buffer(&v_buf);

    if let Some(vector) = v {

        let result = vector.normalize();

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}
