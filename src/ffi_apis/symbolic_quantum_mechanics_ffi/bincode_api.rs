use crate::ffi_apis::common::{BincodeBuffer, from_bincode_buffer, to_bincode_buffer};
use crate::symbolic::quantum_mechanics::Bra;
use crate::symbolic::quantum_mechanics::Ket;
use crate::symbolic::quantum_mechanics::Operator;
use crate::symbolic::quantum_mechanics::{
    self,
};

/// Computes the expectation value using Bincode.
#[no_mangle]

pub extern "C" fn rssn_bincode_expectation_value(
    op_buf: BincodeBuffer,
    psi_buf: BincodeBuffer,
) -> BincodeBuffer {

    let op: Option<Operator> =
        from_bincode_buffer(&op_buf);

    let psi: Option<Ket> =
        from_bincode_buffer(&psi_buf);

    if let (Some(op), Some(psi)) =
        (op, psi)
    {

        to_bincode_buffer(&quantum_mechanics::expectation_value(&op, &psi))
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes the uncertainty using Bincode.
#[no_mangle]

pub extern "C" fn rssn_bincode_uncertainty(
    op_buf: BincodeBuffer,
    psi_buf: BincodeBuffer,
) -> BincodeBuffer {

    let op: Option<Operator> =
        from_bincode_buffer(&op_buf);

    let psi: Option<Ket> =
        from_bincode_buffer(&psi_buf);

    if let (Some(op), Some(psi)) =
        (op, psi)
    {

        to_bincode_buffer(&quantum_mechanics::uncertainty(&op, &psi))
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes the inner product <Bra|Ket> using Bincode.
#[no_mangle]

pub extern "C" fn rssn_bincode_bra_ket(
    bra_buf: BincodeBuffer,
    ket_buf: BincodeBuffer,
) -> BincodeBuffer {

    let bra: Option<Bra> =
        from_bincode_buffer(&bra_buf);

    let ket: Option<Ket> =
        from_bincode_buffer(&ket_buf);

    if let (Some(bra), Some(ket)) =
        (bra, ket)
    {

        to_bincode_buffer(
            &quantum_mechanics::bra_ket(
                &bra, &ket,
            ),
        )
    } else {

        BincodeBuffer::empty()
    }
}
