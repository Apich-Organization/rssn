//! Handle-based FFI API for symbolic PDE functions.

use std::ffi::CStr;
use std::os::raw::c_char;

use crate::symbolic::core::Expr;
use crate::symbolic::pde;

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

unsafe fn c_str_to_str<'a>(
    s: *const c_char
) -> Option<&'a str> { unsafe {

    if s.is_null() {

        None
    } else {

        CStr::from_ptr(s)
            .to_str()
            .ok()
    }
}}

/// Solves a partial differential equation using automatic method selection.
///
/// # Safety
/// The caller must ensure `pde_expr` is a valid Expr pointer, `func` and `vars` are valid C strings,
/// and `vars_len` accurately represents the number of variables.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_solve_pde(
    pde_expr: *const Expr,
    func: *const c_char,
    vars: *const *const c_char,
    vars_len: usize,
) -> *mut Expr { unsafe {

    if pde_expr.is_null()
        || func.is_null()
        || vars.is_null()
    {

        return std::ptr::null_mut();
    }

    let pde_ref = &*pde_expr;

    let func_str = match c_str_to_str(
        func,
    ) {
        | Some(s) => s,
        | None => {
            return std::ptr::null_mut()
        },
    };

    let vars_slice =
        std::slice::from_raw_parts(
            vars,
            vars_len,
        );

    let mut vars_vec =
        Vec::with_capacity(vars_len);

    for &var_ptr in vars_slice {

        match c_str_to_str(var_ptr) {
            | Some(s) => vars_vec.push(s),
            | None => return std::ptr::null_mut(),
        }
    }

    let vars_refs: Vec<&str> =
        vars_vec.clone();

    let result = pde::solve_pde(
        pde_ref,
        func_str,
        &vars_refs,
        None,
    );

    Box::into_raw(Box::new(result))
}}

/// Solves a PDE using the method of characteristics.
///
/// # Safety
/// The caller must ensure all pointers are valid.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_solve_pde_by_characteristics(
    equation: *const Expr,
    func: *const c_char,
    vars: *const *const c_char,
    vars_len: usize,
) -> *mut Expr { unsafe {

    if equation.is_null()
        || func.is_null()
        || vars.is_null()
    {

        return std::ptr::null_mut();
    }

    let eq_ref = &*equation;

    let func_str = match c_str_to_str(
        func,
    ) {
        | Some(s) => s,
        | None => {
            return std::ptr::null_mut()
        },
    };

    let vars_slice =
        std::slice::from_raw_parts(
            vars,
            vars_len,
        );

    let mut vars_vec =
        Vec::with_capacity(vars_len);

    for &var_ptr in vars_slice {

        match c_str_to_str(var_ptr) {
            | Some(s) => vars_vec.push(s),
            | None => return std::ptr::null_mut(),
        }
    }

    let vars_refs: Vec<&str> =
        vars_vec.clone();

    match pde::solve_pde_by_characteristics(
        eq_ref,
        func_str,
        &vars_refs,
    ) {
        | Some(solution) => Box::into_raw(Box::new(solution)),
        | None => std::ptr::null_mut(),
    }
}}

/// Solves the 1D wave equation using D'Alembert's formula.
///
/// # Safety
/// The caller must ensure all pointers are valid.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_solve_wave_equation_1d_dalembert(
    equation: *const Expr,
    func: *const c_char,
    vars: *const *const c_char,
    vars_len: usize,
) -> *mut Expr { unsafe {

    if equation.is_null()
        || func.is_null()
        || vars.is_null()
    {

        return std::ptr::null_mut();
    }

    let eq_ref = &*equation;

    let func_str = match c_str_to_str(
        func,
    ) {
        | Some(s) => s,
        | None => {
            return std::ptr::null_mut()
        },
    };

    let vars_slice =
        std::slice::from_raw_parts(
            vars,
            vars_len,
        );

    let mut vars_vec =
        Vec::with_capacity(vars_len);

    for &var_ptr in vars_slice {

        match c_str_to_str(var_ptr) {
            | Some(s) => vars_vec.push(s),
            | None => return std::ptr::null_mut(),
        }
    }

    let vars_refs: Vec<&str> =
        vars_vec.clone();

    match pde::solve_wave_equation_1d_dalembert(
        eq_ref,
        func_str,
        &vars_refs,
    ) {
        | Some(solution) => Box::into_raw(Box::new(solution)),
        | None => std::ptr::null_mut(),
    }
}}

/// Solves the 1D heat equation.
///
/// # Safety
/// The caller must ensure all pointers are valid.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_solve_heat_equation_1d(
    equation: *const Expr,
    func: *const c_char,
    vars: *const *const c_char,
    vars_len: usize,
) -> *mut Expr { unsafe {

    if equation.is_null()
        || func.is_null()
        || vars.is_null()
    {

        return std::ptr::null_mut();
    }

    let eq_ref = &*equation;

    let func_str = match c_str_to_str(
        func,
    ) {
        | Some(s) => s,
        | None => {
            return std::ptr::null_mut()
        },
    };

    let vars_slice =
        std::slice::from_raw_parts(
            vars,
            vars_len,
        );

    let mut vars_vec =
        Vec::with_capacity(vars_len);

    for &var_ptr in vars_slice {

        match c_str_to_str(var_ptr) {
            | Some(s) => vars_vec.push(s),
            | None => return std::ptr::null_mut(),
        }
    }

    let vars_refs: Vec<&str> =
        vars_vec.clone();

    match pde::solve_heat_equation_1d(
        eq_ref,
        func_str,
        &vars_refs,
    ) {
        | Some(solution) => {
            Box::into_raw(Box::new(
                solution,
            ))
        },
        | None => std::ptr::null_mut(),
    }
}}

/// Solves the 2D Laplace equation.
///
/// # Safety
/// The caller must ensure all pointers are valid.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_solve_laplace_equation_2d(
    equation: *const Expr,
    func: *const c_char,
    vars: *const *const c_char,
    vars_len: usize,
) -> *mut Expr { unsafe {

    if equation.is_null()
        || func.is_null()
        || vars.is_null()
    {

        return std::ptr::null_mut();
    }

    let eq_ref = &*equation;

    let func_str = match c_str_to_str(
        func,
    ) {
        | Some(s) => s,
        | None => {
            return std::ptr::null_mut()
        },
    };

    let vars_slice =
        std::slice::from_raw_parts(
            vars,
            vars_len,
        );

    let mut vars_vec =
        Vec::with_capacity(vars_len);

    for &var_ptr in vars_slice {

        match c_str_to_str(var_ptr) {
            | Some(s) => vars_vec.push(s),
            | None => return std::ptr::null_mut(),
        }
    }

    let vars_refs: Vec<&str> =
        vars_vec.clone();

    match pde::solve_laplace_equation_2d(
        eq_ref,
        func_str,
        &vars_refs,
    ) {
        | Some(solution) => {
            Box::into_raw(Box::new(
                solution,
            ))
        },
        | None => std::ptr::null_mut(),
    }
}}

/// Solves the 2D Poisson equation.
///
/// # Safety
/// The caller must ensure all pointers are valid.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_solve_poisson_equation_2d(
    equation: *const Expr,
    func: *const c_char,
    vars: *const *const c_char,
    vars_len: usize,
) -> *mut Expr { unsafe {

    if equation.is_null()
        || func.is_null()
        || vars.is_null()
    {

        return std::ptr::null_mut();
    }

    let eq_ref = &*equation;

    let func_str = match c_str_to_str(
        func,
    ) {
        | Some(s) => s,
        | None => {
            return std::ptr::null_mut()
        },
    };

    let vars_slice =
        std::slice::from_raw_parts(
            vars,
            vars_len,
        );

    let mut vars_vec =
        Vec::with_capacity(vars_len);

    for &var_ptr in vars_slice {

        match c_str_to_str(var_ptr) {
            | Some(s) => vars_vec.push(s),
            | None => return std::ptr::null_mut(),
        }
    }

    let vars_refs: Vec<&str> =
        vars_vec.clone();

    match pde::solve_poisson_equation_2d(
        eq_ref,
        func_str,
        &vars_refs,
    ) {
        | Some(solution) => {
            Box::into_raw(Box::new(
                solution,
            ))
        },
        | None => std::ptr::null_mut(),
    }
}}

/// Solves the Helmholtz equation.
///
/// # Safety
/// The caller must ensure all pointers are valid.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_solve_helmholtz_equation(
    equation: *const Expr,
    func: *const c_char,
    vars: *const *const c_char,
    vars_len: usize,
) -> *mut Expr { unsafe {

    if equation.is_null()
        || func.is_null()
        || vars.is_null()
    {

        return std::ptr::null_mut();
    }

    let eq_ref = &*equation;

    let func_str = match c_str_to_str(
        func,
    ) {
        | Some(s) => s,
        | None => {
            return std::ptr::null_mut()
        },
    };

    let vars_slice =
        std::slice::from_raw_parts(
            vars,
            vars_len,
        );

    let mut vars_vec =
        Vec::with_capacity(vars_len);

    for &var_ptr in vars_slice {

        match c_str_to_str(var_ptr) {
            | Some(s) => vars_vec.push(s),
            | None => return std::ptr::null_mut(),
        }
    }

    let vars_refs: Vec<&str> =
        vars_vec.clone();

    match pde::solve_helmholtz_equation(
        eq_ref,
        func_str,
        &vars_refs,
    ) {
        | Some(solution) => {
            Box::into_raw(Box::new(
                solution,
            ))
        },
        | None => std::ptr::null_mut(),
    }
}}

/// Solves the Schrödinger equation.
///
/// # Safety
/// The caller must ensure all pointers are valid.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_solve_schrodinger_equation(
    equation: *const Expr,
    func: *const c_char,
    vars: *const *const c_char,
    vars_len: usize,
) -> *mut Expr { unsafe {

    if equation.is_null()
        || func.is_null()
        || vars.is_null()
    {

        return std::ptr::null_mut();
    }

    let eq_ref = &*equation;

    let func_str = match c_str_to_str(
        func,
    ) {
        | Some(s) => s,
        | None => {
            return std::ptr::null_mut()
        },
    };

    let vars_slice =
        std::slice::from_raw_parts(
            vars,
            vars_len,
        );

    let mut vars_vec =
        Vec::with_capacity(vars_len);

    for &var_ptr in vars_slice {

        match c_str_to_str(var_ptr) {
            | Some(s) => vars_vec.push(s),
            | None => return std::ptr::null_mut(),
        }
    }

    let vars_refs: Vec<&str> =
        vars_vec.clone();

    match pde::solve_schrodinger_equation(
        eq_ref,
        func_str,
        &vars_refs,
    ) {
        | Some(solution) => Box::into_raw(Box::new(solution)),
        | None => std::ptr::null_mut(),
    }
}}

/// Solves the Klein-Gordon equation.
///
/// # Safety
/// The caller must ensure all pointers are valid.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_solve_klein_gordon_equation(
    equation: *const Expr,
    func: *const c_char,
    vars: *const *const c_char,
    vars_len: usize,
) -> *mut Expr { unsafe {

    if equation.is_null()
        || func.is_null()
        || vars.is_null()
    {

        return std::ptr::null_mut();
    }

    let eq_ref = &*equation;

    let func_str = match c_str_to_str(
        func,
    ) {
        | Some(s) => s,
        | None => {
            return std::ptr::null_mut()
        },
    };

    let vars_slice =
        std::slice::from_raw_parts(
            vars,
            vars_len,
        );

    let mut vars_vec =
        Vec::with_capacity(vars_len);

    for &var_ptr in vars_slice {

        match c_str_to_str(var_ptr) {
            | Some(s) => vars_vec.push(s),
            | None => return std::ptr::null_mut(),
        }
    }

    let vars_refs: Vec<&str> =
        vars_vec.clone();

    match pde::solve_klein_gordon_equation(
        eq_ref,
        func_str,
        &vars_refs,
    ) {
        | Some(solution) => Box::into_raw(Box::new(solution)),
        | None => std::ptr::null_mut(),
    }
}}
