//! Bincode-based FFI API for symbolic PDE functions.

use std::os::raw::c_char;

use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::pde;

/// Solves a PDE using Bincode with automatic method selection.
#[no_mangle]

pub extern "C" fn rssn_bincode_solve_pde(
    pde_buf : BincodeBuffer,
    func : *const c_char,
    vars_buf : BincodeBuffer,
) -> BincodeBuffer {

    let pde_expr : Option<Expr> =
        from_bincode_buffer(&pde_buf);

    let vars : Option<Vec<String>> =
        from_bincode_buffer(&vars_buf);

    let func_str = unsafe {

        if func.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(
                func,
            )
            .to_str()
            .ok()
        }
    };

    if let (
        Some(pde),
        Some(f),
        Some(v),
    ) = (
        pde_expr,
        func_str,
        vars,
    ) {

        let vars_refs : Vec<&str> = v
            .iter()
            .map(|s| s.as_str())
            .collect();

        let result = pde::solve_pde(
            &pde,
            &f,
            &vars_refs,
            None,
        );

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

/// Solves a PDE using the method of characteristics (Bincode).
#[no_mangle]

pub extern "C" fn rssn_bincode_solve_pde_by_characteristics(
    equation_buf : BincodeBuffer,
    func : *const c_char,
    vars_buf : BincodeBuffer,
) -> BincodeBuffer {

    let equation : Option<Expr> =
        from_bincode_buffer(
            &equation_buf,
        );

    let vars : Option<Vec<String>> =
        from_bincode_buffer(&vars_buf);

    let func_str = unsafe {

        if func.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(
                func,
            )
            .to_str()
            .ok()
        }
    };

    if let (
        Some(eq),
        Some(f),
        Some(v),
    ) = (
        equation,
        func_str,
        vars,
    ) {

        let vars_refs : Vec<&str> = v
            .iter()
            .map(|s| s.as_str())
            .collect();

        match pde::solve_pde_by_characteristics(&eq, f, &vars_refs) {
            | Some(result) => to_bincode_buffer(&result),
            | None => BincodeBuffer::empty(),
        }
    } else {

        BincodeBuffer::empty()
    }
}

/// Solves the 1D wave equation (Bincode).
#[no_mangle]

pub extern "C" fn rssn_bincode_solve_wave_equation_1d(
    equation_buf : BincodeBuffer,
    func : *const c_char,
    vars_buf : BincodeBuffer,
) -> BincodeBuffer {

    let equation : Option<Expr> =
        from_bincode_buffer(
            &equation_buf,
        );

    let vars : Option<Vec<String>> =
        from_bincode_buffer(&vars_buf);

    let func_str = unsafe {

        if func.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(
                func,
            )
            .to_str()
            .ok()
        }
    };

    if let (
        Some(eq),
        Some(f),
        Some(v),
    ) = (
        equation,
        func_str,
        vars,
    ) {

        let vars_refs : Vec<&str> = v
            .iter()
            .map(|s| s.as_str())
            .collect();

        match pde::solve_wave_equation_1d_dalembert(&eq, f, &vars_refs) {
            | Some(result) => to_bincode_buffer(&result),
            | None => BincodeBuffer::empty(),
        }
    } else {

        BincodeBuffer::empty()
    }
}

/// Solves the 1D heat equation (Bincode).
#[no_mangle]

pub extern "C" fn rssn_bincode_solve_heat_equation_1d(
    equation_buf : BincodeBuffer,
    func : *const c_char,
    vars_buf : BincodeBuffer,
) -> BincodeBuffer {

    let equation : Option<Expr> =
        from_bincode_buffer(
            &equation_buf,
        );

    let vars : Option<Vec<String>> =
        from_bincode_buffer(&vars_buf);

    let func_str = unsafe {

        if func.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(
                func,
            )
            .to_str()
            .ok()
        }
    };

    if let (
        Some(eq),
        Some(f),
        Some(v),
    ) = (
        equation,
        func_str,
        vars,
    ) {

        let vars_refs : Vec<&str> = v
            .iter()
            .map(|s| s.as_str())
            .collect();

        match pde::solve_heat_equation_1d(&eq, f, &vars_refs) {
            | Some(result) => to_bincode_buffer(&result),
            | None => BincodeBuffer::empty(),
        }
    } else {

        BincodeBuffer::empty()
    }
}

/// Solves the 2D Laplace equation (Bincode).
#[no_mangle]

pub extern "C" fn rssn_bincode_solve_laplace_equation_2d(
    equation_buf : BincodeBuffer,
    func : *const c_char,
    vars_buf : BincodeBuffer,
) -> BincodeBuffer {

    let equation : Option<Expr> =
        from_bincode_buffer(
            &equation_buf,
        );

    let vars : Option<Vec<String>> =
        from_bincode_buffer(&vars_buf);

    let func_str = unsafe {

        if func.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(
                func,
            )
            .to_str()
            .ok()
        }
    };

    if let (
        Some(eq),
        Some(f),
        Some(v),
    ) = (
        equation,
        func_str,
        vars,
    ) {

        let vars_refs : Vec<&str> = v
            .iter()
            .map(|s| s.as_str())
            .collect();

        match pde::solve_laplace_equation_2d(&eq, f, &vars_refs) {
            | Some(result) => to_bincode_buffer(&result),
            | None => BincodeBuffer::empty(),
        }
    } else {

        BincodeBuffer::empty()
    }
}

/// Classifies a PDE (Bincode).
#[no_mangle]

pub extern "C" fn rssn_bincode_classify_pde(
    equation_buf : BincodeBuffer,
    func : *const c_char,
    vars_buf : BincodeBuffer,
) -> BincodeBuffer {

    let equation : Option<Expr> =
        from_bincode_buffer(
            &equation_buf,
        );

    let vars : Option<Vec<String>> =
        from_bincode_buffer(&vars_buf);

    let func_str = unsafe {

        if func.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(
                func,
            )
            .to_str()
            .ok()
        }
    };

    if let (
        Some(eq),
        Some(f),
        Some(v),
    ) = (
        equation,
        func_str,
        vars,
    ) {

        let vars_refs : Vec<&str> = v
            .iter()
            .map(|s| s.as_str())
            .collect();

        let classification =
            pde::classify_pde_heuristic(
                &eq,
                f,
                &vars_refs,
            );

        to_bincode_buffer(
            &classification,
        )
    } else {

        BincodeBuffer::empty()
    }
}
