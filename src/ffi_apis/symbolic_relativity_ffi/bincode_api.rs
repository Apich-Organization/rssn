//! Bincode-based FFI API for relativity functions.

use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::relativity;

/// Calculates Lorentz factor using Bincode.
#[no_mangle]

pub extern "C" fn rssn_bincode_lorentz_factor(
    velocity_buf : BincodeBuffer
) -> BincodeBuffer {

    let velocity : Option<Expr> =
        from_bincode_buffer(
            &velocity_buf,
        );

    if let Some(v) = velocity {

        to_bincode_buffer(
            &relativity::lorentz_factor(
                &v,
            ),
        )
    } else {

        BincodeBuffer::empty()
    }
}

/// Calculates mass-energy equivalence using Bincode.
#[no_mangle]

pub extern "C" fn rssn_bincode_mass_energy_equivalence(
    mass_buf : BincodeBuffer
) -> BincodeBuffer {

    let mass : Option<Expr> =
        from_bincode_buffer(&mass_buf);

    if let Some(m) = mass {

        to_bincode_buffer(&relativity::mass_energy_equivalence(&m))
    } else {

        BincodeBuffer::empty()
    }
}

/// Calculates Schwarzschild radius using Bincode.
#[no_mangle]

pub extern "C" fn rssn_bincode_schwarzschild_radius(
    mass_buf : BincodeBuffer
) -> BincodeBuffer {

    let mass : Option<Expr> =
        from_bincode_buffer(&mass_buf);

    if let Some(m) = mass {

        to_bincode_buffer(&relativity::schwarzschild_radius(&m))
    } else {

        BincodeBuffer::empty()
    }
}
