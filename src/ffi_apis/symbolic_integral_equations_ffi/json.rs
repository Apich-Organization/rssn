//! JSON-based FFI API for symbolic integral equations.

use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::integral_equations::{
    solve_airfoil_equation,
    FredholmEquation,
    VolterraEquation,
};
use serde::{
    Deserialize,
    Serialize,
};
use std::os::raw::c_char;

#[derive(Serialize, Deserialize)]

struct FredholmNeumannInput {
    equation : FredholmEquation,
    iterations : usize,
}

#[derive(Serialize, Deserialize)]

struct FredholmSeparableInput {
    equation : FredholmEquation,
    a_funcs : Vec<Expr>,
    b_funcs : Vec<Expr>,
}

#[derive(Serialize, Deserialize)]

struct VolterraSuccessiveInput {
    equation : VolterraEquation,
    iterations : usize,
}

#[derive(Serialize, Deserialize)]

struct AirfoilInput {
    f_x : Expr,
    var_x : String,
    var_t : String,
}

/// Solves a Fredholm equation using the Neumann series method (JSON).
#[no_mangle]

pub extern "C" fn rssn_fredholm_solve_neumann_json(
    input_json : *const c_char
) -> *mut c_char {

    let input : Option<
        FredholmNeumannInput,
    > = from_json_string(input_json);

    let input = match input {
        | Some(i) => i,
        | None => {
            return std::ptr::null_mut()
        },
    };

    let result = input
        .equation
        .solve_neumann_series(
            input.iterations,
        );

    to_json_string(&result)
}

/// Solves a Fredholm equation with a separable kernel (JSON).
#[no_mangle]

pub extern "C" fn rssn_fredholm_solve_separable_json(
    input_json : *const c_char
) -> *mut c_char {

    let input : Option<
        FredholmSeparableInput,
    > = from_json_string(input_json);

    let input = match input {
        | Some(i) => i,
        | None => {
            return std::ptr::null_mut()
        },
    };

    match input
        .equation
        .solve_separable_kernel(
            input.a_funcs,
            input.b_funcs,
        ) {
        | Ok(result) => {
            to_json_string(&result)
        },
        | Err(_) => {
            std::ptr::null_mut()
        },
    }
}

/// Solves a Volterra equation using successive approximations (JSON).
#[no_mangle]

pub extern "C" fn rssn_volterra_solve_successive_json(
    input_json : *const c_char
) -> *mut c_char {

    let input : Option<
        VolterraSuccessiveInput,
    > = from_json_string(input_json);

    let input = match input {
        | Some(i) => i,
        | None => {
            return std::ptr::null_mut()
        },
    };

    let result = input
        .equation
        .solve_successive_approximations(input.iterations);

    to_json_string(&result)
}

/// Solves a Volterra equation by differentiation (JSON).
#[no_mangle]

pub extern "C" fn rssn_volterra_solve_by_differentiation_json(
    input_json : *const c_char
) -> *mut c_char {

    let equation : Option<
        VolterraEquation,
    > = from_json_string(input_json);

    let equation = match equation {
        | Some(e) => e,
        | None => {
            return std::ptr::null_mut()
        },
    };

    match equation
        .solve_by_differentiation()
    {
        | Ok(result) => {
            to_json_string(&result)
        },
        | Err(_) => {
            std::ptr::null_mut()
        },
    }
}

/// Solves the airfoil singular integral equation (JSON).
#[no_mangle]

pub extern "C" fn rssn_solve_airfoil_equation_json(
    input_json : *const c_char
) -> *mut c_char {

    let input : Option<AirfoilInput> =
        from_json_string(input_json);

    let input = match input {
        | Some(i) => i,
        | None => {
            return std::ptr::null_mut()
        },
    };

    let result = solve_airfoil_equation(
        &input.f_x,
        &input.var_x,
        &input.var_t,
    );

    to_json_string(&result)
}
