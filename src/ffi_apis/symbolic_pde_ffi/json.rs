//! JSON-based FFI API for symbolic PDE functions.

use std::os::raw::c_char;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_json_string;
use crate::symbolic::core::Expr;
use crate::symbolic::pde;

/// Solves a PDE using JSON with automatic method selection.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_solve_pde(
    pde_json: *const c_char,
    func: *const c_char,
    vars_json: *const c_char,
) -> *mut c_char {

    let pde_expr: Option<Expr> =
        from_json_string(pde_json);

    let vars: Option<Vec<String>> =
        from_json_string(vars_json);

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

    match (
        pde_expr,
        func_str,
        vars,
    ) { (
        Some(pde),
        Some(f),
        Some(v),
    ) => {

        let vars_refs: Vec<&str> = v
            .iter()
            .map(std::string::String::as_str)
            .collect();

        let result = pde::solve_pde(
            &pde,
            f,
            &vars_refs,
            None,
        );

        to_json_string(&result)
    } _ => {

        std::ptr::null_mut()
    }}
}

/// Solves a PDE using the method of characteristics (JSON).
#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_solve_pde_by_characteristics(
    equation_json: *const c_char,
    func: *const c_char,
    vars_json: *const c_char,
) -> *mut c_char {

    let equation: Option<Expr> =
        from_json_string(equation_json);

    let vars: Option<Vec<String>> =
        from_json_string(vars_json);

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

    match (
        equation,
        func_str,
        vars,
    ) { (
        Some(eq),
        Some(f),
        Some(v),
    ) => {

        let vars_refs: Vec<&str> = v
            .iter()
            .map(std::string::String::as_str)
            .collect();

        match pde::solve_pde_by_characteristics(&eq, f, &vars_refs) {
            | Some(result) => to_json_string(&result),
            | None => std::ptr::null_mut(),
        }
    } _ => {

        std::ptr::null_mut()
    }}
}

/// Solves the 1D wave equation using D'Alembert's formula (JSON).
#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_solve_wave_equation_1d(
    equation_json: *const c_char,
    func: *const c_char,
    vars_json: *const c_char,
) -> *mut c_char {

    let equation: Option<Expr> =
        from_json_string(equation_json);

    let vars: Option<Vec<String>> =
        from_json_string(vars_json);

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

    match (
        equation,
        func_str,
        vars,
    ) { (
        Some(eq),
        Some(f),
        Some(v),
    ) => {

        let vars_refs: Vec<&str> = v
            .iter()
            .map(std::string::String::as_str)
            .collect();

        match pde::solve_wave_equation_1d_dalembert(&eq, f, &vars_refs) {
            | Some(result) => to_json_string(&result),
            | None => std::ptr::null_mut(),
        }
    } _ => {

        std::ptr::null_mut()
    }}
}

/// Solves the 1D heat equation (JSON).
#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_solve_heat_equation_1d(
    equation_json: *const c_char,
    func: *const c_char,
    vars_json: *const c_char,
) -> *mut c_char {

    let equation: Option<Expr> =
        from_json_string(equation_json);

    let vars: Option<Vec<String>> =
        from_json_string(vars_json);

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

    match (
        equation,
        func_str,
        vars,
    ) { (
        Some(eq),
        Some(f),
        Some(v),
    ) => {

        let vars_refs: Vec<&str> = v
            .iter()
            .map(std::string::String::as_str)
            .collect();

        match pde::solve_heat_equation_1d(&eq, f, &vars_refs) {
            | Some(result) => to_json_string(&result),
            | None => std::ptr::null_mut(),
        }
    } _ => {

        std::ptr::null_mut()
    }}
}

/// Solves the 2D Laplace equation (JSON).
#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_solve_laplace_equation_2d(
    equation_json: *const c_char,
    func: *const c_char,
    vars_json: *const c_char,
) -> *mut c_char {

    let equation: Option<Expr> =
        from_json_string(equation_json);

    let vars: Option<Vec<String>> =
        from_json_string(vars_json);

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

    match (
        equation,
        func_str,
        vars,
    ) { (
        Some(eq),
        Some(f),
        Some(v),
    ) => {

        let vars_refs: Vec<&str> = v
            .iter()
            .map(std::string::String::as_str)
            .collect();

        match pde::solve_laplace_equation_2d(&eq, f, &vars_refs) {
            | Some(result) => to_json_string(&result),
            | None => std::ptr::null_mut(),
        }
    } _ => {

        std::ptr::null_mut()
    }}
}

/// Solves the 2D Poisson equation (JSON).
#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_solve_poisson_equation_2d(
    equation_json: *const c_char,
    func: *const c_char,
    vars_json: *const c_char,
) -> *mut c_char {

    let equation: Option<Expr> =
        from_json_string(equation_json);

    let vars: Option<Vec<String>> =
        from_json_string(vars_json);

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

    match (
        equation,
        func_str,
        vars,
    ) { (
        Some(eq),
        Some(f),
        Some(v),
    ) => {

        let vars_refs: Vec<&str> = v
            .iter()
            .map(std::string::String::as_str)
            .collect();

        match pde::solve_poisson_equation_2d(&eq, f, &vars_refs) {
            | Some(result) => to_json_string(&result),
            | None => std::ptr::null_mut(),
        }
    } _ => {

        std::ptr::null_mut()
    }}
}

/// Classifies a PDE and suggests solution methods (JSON).
#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_classify_pde(
    equation_json: *const c_char,
    func: *const c_char,
    vars_json: *const c_char,
) -> *mut c_char {

    let equation: Option<Expr> =
        from_json_string(equation_json);

    let vars: Option<Vec<String>> =
        from_json_string(vars_json);

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

    match (
        equation,
        func_str,
        vars,
    ) { (
        Some(eq),
        Some(f),
        Some(v),
    ) => {

        let vars_refs: Vec<&str> = v
            .iter()
            .map(std::string::String::as_str)
            .collect();

        let classification =
            pde::classify_pde_heuristic(
                &eq,
                f,
                &vars_refs,
            );

        to_json_string(&classification)
    } _ => {

        std::ptr::null_mut()
    }}
}
