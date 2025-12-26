use crate::ffi_apis::common::*;
use crate::symbolic::tensor::*;

#[no_mangle]

pub extern "C" fn rssn_bincode_tensor_add(
    t1_buf : BincodeBuffer,
    t2_buf : BincodeBuffer,
) -> BincodeBuffer {

    let t1 : Option<Tensor> = from_bincode_buffer(&t1_buf);

    let t2 : Option<Tensor> = from_bincode_buffer(&t2_buf);

    if let (Some(tensor1), Some(tensor2)) = (t1, t2) {

        match tensor1.add(&tensor2) {
            | Ok(result) => to_bincode_buffer(&result),
            | Err(_) => BincodeBuffer::empty(),
        }
    } else {

        BincodeBuffer::empty()
    }
}

#[no_mangle]

pub extern "C" fn rssn_bincode_tensor_scalar_mul(
    t_buf : BincodeBuffer,
    scalar_buf : BincodeBuffer,
) -> BincodeBuffer {

    let t : Option<Tensor> = from_bincode_buffer(&t_buf);

    let scalar : Option<crate::symbolic::core::Expr> = from_bincode_buffer(&scalar_buf);

    if let (Some(tensor), Some(s)) = (t, scalar) {

        match tensor.scalar_mul(&s) {
            | Ok(result) => to_bincode_buffer(&result),
            | Err(_) => BincodeBuffer::empty(),
        }
    } else {

        BincodeBuffer::empty()
    }
}

#[no_mangle]

pub extern "C" fn rssn_bincode_tensor_outer_product(
    t1_buf : BincodeBuffer,
    t2_buf : BincodeBuffer,
) -> BincodeBuffer {

    let t1 : Option<Tensor> = from_bincode_buffer(&t1_buf);

    let t2 : Option<Tensor> = from_bincode_buffer(&t2_buf);

    if let (Some(tensor1), Some(tensor2)) = (t1, t2) {

        match tensor1.outer_product(&tensor2) {
            | Ok(result) => to_bincode_buffer(&result),
            | Err(_) => BincodeBuffer::empty(),
        }
    } else {

        BincodeBuffer::empty()
    }
}
