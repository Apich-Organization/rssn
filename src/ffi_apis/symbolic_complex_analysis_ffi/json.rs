use std::os::raw::c_char;

use crate::ffi_apis::common::*;
use crate::symbolic::complex_analysis::*;
use crate::symbolic::core::Expr;

#[no_mangle]

pub unsafe extern "C" fn path_continuation_new_json(
    func_json : *const c_char,
    var : *const c_char,
    start_point_json : *const c_char,
    order : usize,
) -> *mut c_char {

    let func: Expr = match from_json_string(func_json) {
        | Some(e) => e,
        | None => return std::ptr::null_mut(),
    };

    let var_str =
        std::ffi::CStr::from_ptr(var)
            .to_str()
            .unwrap();

    let start_point: Expr =
        match from_json_string(start_point_json) {
            | Some(e) => e,
            | None => return std::ptr::null_mut(),
        };

    let path_continuation =
        PathContinuation::new(
            &func,
            var_str,
            &start_point,
            order,
        );

    to_json_string(&path_continuation)
}

#[no_mangle]

pub unsafe extern "C" fn path_continuation_continue_along_path_json(
    pc_json : *const c_char,
    path_points_json : *const c_char,
) -> *mut c_char {

    let mut pc: PathContinuation =
        match from_json_string(pc_json) {
            | Some(e) => e,
            | None => return std::ptr::null_mut(),
        };

    let path_points: Vec<Expr> =
        match from_json_string(path_points_json) {
            | Some(e) => e,
            | None => return std::ptr::null_mut(),
        };

    match pc.continue_along_path(
        &path_points,
    ) {
        | Ok(_) => {
            to_c_string(
                "OK".to_string(),
            )
        },
        | Err(e) => to_c_string(e),
    }
}

#[no_mangle]

pub unsafe extern "C" fn path_continuation_get_final_expression_json(
    pc_json : *const c_char
) -> *mut c_char {

    let pc: PathContinuation =
        match from_json_string(pc_json) {
            | Some(e) => e,
            | None => return std::ptr::null_mut(),
        };

    to_json_string(
        &pc.get_final_expression(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn estimate_radius_of_convergence_json(
    series_expr_json : *const c_char,
    var : *const c_char,
    center_json : *const c_char,
    order : usize,
) -> f64 {

    let series_expr : Expr =
        match from_json_string(
            series_expr_json,
        ) {
            | Some(e) => e,
            | None => return 0.0,
        };

    let var_str =
        std::ffi::CStr::from_ptr(var)
            .to_str()
            .unwrap();

    let center : Expr =
        match from_json_string(
            center_json,
        ) {
            | Some(e) => e,
            | None => return 0.0,
        };

    crate::symbolic::complex_analysis::estimate_radius_of_convergence(
        &series_expr,
        var_str,
        &center,
        order,
    )
    .unwrap_or(0.0)
}

#[no_mangle]

pub unsafe extern "C" fn complex_distance_json(
    p1_json : *const c_char,
    p2_json : *const c_char,
) -> f64 {

    let p1 : Expr =
        match from_json_string(p1_json)
        {
            | Some(e) => e,
            | None => return 0.0,
        };

    let p2 : Expr =
        match from_json_string(p2_json)
        {
            | Some(e) => e,
            | None => return 0.0,
        };

    crate::symbolic::complex_analysis::complex_distance(
        &p1, &p2,
    )
    .unwrap_or(0.0)
}

#[no_mangle]

pub unsafe extern "C" fn classify_singularity_json(
    func_json : *const c_char,
    var : *const c_char,
    singularity_json : *const c_char,
    order : usize,
) -> *mut c_char {

    let func: Expr = match from_json_string(func_json) {
        | Some(e) => e,
        | None => return std::ptr::null_mut(),
    };

    let var_str =
        std::ffi::CStr::from_ptr(var)
            .to_str()
            .unwrap();

    let singularity: Expr =
        match from_json_string(singularity_json) {
            | Some(e) => e,
            | None => return std::ptr::null_mut(),
        };

    let singularity_type = crate::symbolic::complex_analysis::classify_singularity(
        &func,
        var_str,
        &singularity,
        order,
    );

    to_json_string(&singularity_type)
}

#[no_mangle]

pub unsafe extern "C" fn laurent_series_json(
    func_json : *const c_char,
    var : *const c_char,
    center_json : *const c_char,
    order : usize,
) -> *mut c_char {

    let func: Expr = match from_json_string(func_json) {
        | Some(e) => e,
        | None => return std::ptr::null_mut(),
    };

    let var_str =
        std::ffi::CStr::from_ptr(var)
            .to_str()
            .unwrap();

    let center: Expr = match from_json_string(center_json) {
        | Some(e) => e,
        | None => return std::ptr::null_mut(),
    };

    let series =
        crate::symbolic::complex_analysis::laurent_series(
            &func, var_str, &center, order,
        );

    to_json_string(&series)
}

#[no_mangle]

pub unsafe extern "C" fn calculate_residue_json(
    func_json : *const c_char,
    var : *const c_char,
    singularity_json : *const c_char,
) -> *mut c_char {

    let func: Expr = match from_json_string(func_json) {
        | Some(e) => e,
        | None => return std::ptr::null_mut(),
    };

    let var_str =
        std::ffi::CStr::from_ptr(var)
            .to_str()
            .unwrap();

    let singularity: Expr =
        match from_json_string(singularity_json) {
            | Some(e) => e,
            | None => return std::ptr::null_mut(),
        };

    let residue = crate::symbolic::complex_analysis::calculate_residue(
        &func,
        var_str,
        &singularity,
    );

    to_json_string(&residue)
}

#[no_mangle]

pub unsafe extern "C" fn contour_integral_residue_theorem_json(
    func_json : *const c_char,
    var : *const c_char,
    singularities_json : *const c_char,
) -> *mut c_char {

    let func: Expr = match from_json_string(func_json) {
        | Some(e) => e,
        | None => return std::ptr::null_mut(),
    };

    let var_str =
        std::ffi::CStr::from_ptr(var)
            .to_str()
            .unwrap();

    let singularities: Vec<Expr> =
        match from_json_string(singularities_json) {
            | Some(e) => e,
            | None => return std::ptr::null_mut(),
        };

    let result = crate::symbolic::complex_analysis::contour_integral_residue_theorem(
        &func,
        var_str,
        &singularities,
    );

    to_json_string(&result)
}

#[no_mangle]

pub unsafe extern "C" fn mobius_transformation_new_json(
    a_json : *const c_char,
    b_json : *const c_char,
    c_json : *const c_char,
    d_json : *const c_char,
) -> *mut c_char {

    let a: Expr = match from_json_string(
        a_json,
    ) {
        | Some(e) => e,
        | None => {
            return std::ptr::null_mut()
        },
    };

    let b: Expr = match from_json_string(
        b_json,
    ) {
        | Some(e) => e,
        | None => {
            return std::ptr::null_mut()
        },
    };

    let c: Expr = match from_json_string(
        c_json,
    ) {
        | Some(e) => e,
        | None => {
            return std::ptr::null_mut()
        },
    };

    let d: Expr = match from_json_string(
        d_json,
    ) {
        | Some(e) => e,
        | None => {
            return std::ptr::null_mut()
        },
    };

    let mobius =
        MobiusTransformation::new(
            a, b, c, d,
        );

    to_json_string(&mobius)
}

#[no_mangle]

pub extern "C" fn mobius_transformation_identity_json(
) -> *mut c_char {

    let mobius =
        MobiusTransformation::identity(
        );

    to_json_string(&mobius)
}

#[no_mangle]

pub unsafe extern "C" fn mobius_transformation_apply_json(
    mobius_json : *const c_char,
    z_json : *const c_char,
) -> *mut c_char {

    let mobius: MobiusTransformation =
        match from_json_string(mobius_json) {
            | Some(e) => e,
            | None => return std::ptr::null_mut(),
        };

    let z: Expr = match from_json_string(
        z_json,
    ) {
        | Some(e) => e,
        | None => {
            return std::ptr::null_mut()
        },
    };

    let result = mobius.apply(&z);

    to_json_string(&result)
}

#[no_mangle]

pub unsafe extern "C" fn mobius_transformation_compose_json(
    mobius1_json : *const c_char,
    mobius2_json : *const c_char,
) -> *mut c_char {

    let mobius1: MobiusTransformation =
        match from_json_string(mobius1_json) {
            | Some(e) => e,
            | None => return std::ptr::null_mut(),
        };

    let mobius2: MobiusTransformation =
        match from_json_string(mobius2_json) {
            | Some(e) => e,
            | None => return std::ptr::null_mut(),
        };

    let result =
        mobius1.compose(&mobius2);

    to_json_string(&result)
}

#[no_mangle]

pub unsafe extern "C" fn mobius_transformation_inverse_json(
    mobius_json : *const c_char
) -> *mut c_char {

    let mobius: MobiusTransformation =
        match from_json_string(mobius_json) {
            | Some(e) => e,
            | None => return std::ptr::null_mut(),
        };

    let result = mobius.inverse();

    to_json_string(&result)
}

#[no_mangle]

pub unsafe extern "C" fn cauchy_integral_formula_json(
    func_json : *const c_char,
    var : *const c_char,
    z0_json : *const c_char,
) -> *mut c_char {

    let func: Expr = match from_json_string(func_json) {
        | Some(e) => e,
        | None => return std::ptr::null_mut(),
    };

    let var_str =
        std::ffi::CStr::from_ptr(var)
            .to_str()
            .unwrap();

    let z0: Expr = match from_json_string(z0_json) {
        | Some(e) => e,
        | None => return std::ptr::null_mut(),
    };

    let result = crate::symbolic::complex_analysis::cauchy_integral_formula(&func, var_str, &z0);

    to_json_string(&result)
}

#[no_mangle]

pub unsafe extern "C" fn cauchy_derivative_formula_json(
    func_json : *const c_char,
    var : *const c_char,
    z0_json : *const c_char,
    n : usize,
) -> *mut c_char {

    let func: Expr = match from_json_string(func_json) {
        | Some(e) => e,
        | None => return std::ptr::null_mut(),
    };

    let var_str =
        std::ffi::CStr::from_ptr(var)
            .to_str()
            .unwrap();

    let z0: Expr = match from_json_string(z0_json) {
        | Some(e) => e,
        | None => return std::ptr::null_mut(),
    };

    let result = crate::symbolic::complex_analysis::cauchy_derivative_formula(
        &func, var_str, &z0, n,
    );

    to_json_string(&result)
}

#[no_mangle]

pub unsafe extern "C" fn complex_exp_json(
    z_json : *const c_char
) -> *mut c_char {

    let z: Expr = match from_json_string(
        z_json,
    ) {
        | Some(e) => e,
        | None => {
            return std::ptr::null_mut()
        },
    };

    let result =
        crate::symbolic::complex_analysis::complex_exp(&z);

    to_json_string(&result)
}

#[no_mangle]

pub unsafe extern "C" fn complex_log_json(
    z_json : *const c_char
) -> *mut c_char {

    let z: Expr = match from_json_string(
        z_json,
    ) {
        | Some(e) => e,
        | None => {
            return std::ptr::null_mut()
        },
    };

    let result =
        crate::symbolic::complex_analysis::complex_log(&z);

    to_json_string(&result)
}

#[no_mangle]

pub unsafe extern "C" fn complex_arg_json(
    z_json : *const c_char
) -> *mut c_char {

    let z: Expr = match from_json_string(
        z_json,
    ) {
        | Some(e) => e,
        | None => {
            return std::ptr::null_mut()
        },
    };

    let result =
        crate::symbolic::complex_analysis::complex_arg(&z);

    to_json_string(&result)
}

#[no_mangle]

pub unsafe extern "C" fn complex_modulus_json(
    z_json : *const c_char
) -> *mut c_char {

    let z: Expr = match from_json_string(
        z_json,
    ) {
        | Some(e) => e,
        | None => {
            return std::ptr::null_mut()
        },
    };

    let result =
        crate::symbolic::complex_analysis::complex_modulus(
            &z,
        );

    to_json_string(&result)
}
