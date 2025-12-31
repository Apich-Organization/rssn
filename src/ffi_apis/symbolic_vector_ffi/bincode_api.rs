use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::symbolic::vector::Vector;

/// Computes the magnitude of a vector.

///

/// Takes a bincode-serialized `Vector` as input.

/// Returns a bincode-serialized `Expr` representing its magnitude.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_vector_magnitude(
    v_buf: BincodeBuffer
) -> BincodeBuffer {

    let v: Option<Vector> =
        from_bincode_buffer(&v_buf);

    if let Some(vector) = v {

        let result = vector.magnitude();

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes the dot product of two vectors.

///

/// Takes two bincode-serialized `Vector` objects as input.

/// Returns a bincode-serialized `Expr` representing their dot product.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_vector_dot(
    v1_buf: BincodeBuffer,
    v2_buf: BincodeBuffer,
) -> BincodeBuffer {

    let v1: Option<Vector> =
        from_bincode_buffer(&v1_buf);

    let v2: Option<Vector> =
        from_bincode_buffer(&v2_buf);

    match (v1, v2) {
        | (Some(vec1), Some(vec2)) => {

            let result =
                vec1.dot(&vec2);

            to_bincode_buffer(&result)
        },
        | _ => BincodeBuffer::empty(),
    }
}

/// Computes the cross product of two vectors.

///

/// Takes two bincode-serialized `Vector` objects as input.

/// Returns a bincode-serialized `Vector` representing their cross product.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_vector_cross(
    v1_buf: BincodeBuffer,
    v2_buf: BincodeBuffer,
) -> BincodeBuffer {

    let v1: Option<Vector> =
        from_bincode_buffer(&v1_buf);

    let v2: Option<Vector> =
        from_bincode_buffer(&v2_buf);

    match (v1, v2) {
        | (Some(vec1), Some(vec2)) => {

            let result =
                vec1.cross(&vec2);

            to_bincode_buffer(&result)
        },
        | _ => BincodeBuffer::empty(),
    }
}

/// Normalizes a vector.

///

/// Takes a bincode-serialized `Vector` as input.

/// Returns a bincode-serialized `Vector` representing the normalized vector.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_vector_normalize(
    v_buf: BincodeBuffer
) -> BincodeBuffer {

    let v: Option<Vector> =
        from_bincode_buffer(&v_buf);

    if let Some(vector) = v {

        let result = vector.normalize();

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}
