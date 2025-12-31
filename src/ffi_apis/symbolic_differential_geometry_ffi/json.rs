use crate::ffi_apis::common::{from_json_string, to_json_string};
use crate::symbolic::core::Expr;
use crate::symbolic::differential_geometry::{DifferentialForm, exterior_derivative, wedge_product, boundary, generalized_stokes_theorem, gauss_theorem, stokes_theorem, greens_theorem};
use crate::symbolic::vector::Vector;

/// Computes the exterior derivative of a differential form (JSON)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_exterior_derivative(
    form_json : *const std::os::raw::c_char,
    vars_json : *const std::os::raw::c_char,
) -> *mut std::os::raw::c_char {

    let form: Option<DifferentialForm> =
        from_json_string(form_json);

    let vars: Option<Vec<String>> =
        from_json_string(vars_json);

    match (form, vars) {
        | (Some(f), Some(v)) => {

            let vars_refs: Vec<&str> = v
            .iter()
            .map(std::string::String::as_str)
            .collect();

            let result =
                exterior_derivative(
                    &f,
                    &vars_refs,
                );

            to_json_string(&result)
        },
        | _ => std::ptr::null_mut(),
    }
}

/// Computes the wedge product of two differential forms (JSON)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_wedge_product(
    form1_json : *const std::os::raw::c_char,
    form2_json : *const std::os::raw::c_char,
) -> *mut std::os::raw::c_char {

    let form1: Option<
        DifferentialForm,
    > = from_json_string(form1_json);

    let form2: Option<
        DifferentialForm,
    > = from_json_string(form2_json);

    match (form1, form2) {
        | (Some(f1), Some(f2)) => {

            let result =
                wedge_product(&f1, &f2);

            to_json_string(&result)
        },
        | _ => std::ptr::null_mut(),
    }
}

/// Computes the boundary of a domain (JSON)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_boundary(
    domain_json : *const std::os::raw::c_char
) -> *mut std::os::raw::c_char {

    let domain: Option<Expr> =
        from_json_string(domain_json);

    if let Some(d) = domain {

        let result = boundary(&d);

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}

/// Represents the generalized Stokes' theorem (JSON)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_generalized_stokes_theorem(
    omega_json : *const std::os::raw::c_char,
    manifold_json : *const std::os::raw::c_char,
    vars_json : *const std::os::raw::c_char,
) -> *mut std::os::raw::c_char {

    let omega: Option<
        DifferentialForm,
    > = from_json_string(omega_json);

    let manifold: Option<Expr> =
        from_json_string(manifold_json);

    let vars: Option<Vec<String>> =
        from_json_string(vars_json);

    match (
        omega,
        manifold,
        vars,
    ) {
        | (
            Some(o),
            Some(m),
            Some(v),
        ) => {

            let vars_refs: Vec<&str> = v
            .iter()
            .map(std::string::String::as_str)
            .collect();

            let result =
            generalized_stokes_theorem(
                &o,
                &m,
                &vars_refs,
            );

            to_json_string(&result)
        },
        | _ => std::ptr::null_mut(),
    }
}

/// Represents Gauss's theorem (JSON)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_gauss_theorem(
    vector_field_json : *const std::os::raw::c_char,
    volume_json : *const std::os::raw::c_char,
) -> *mut std::os::raw::c_char {

    let vector_field: Option<Vector> =
        from_json_string(
            vector_field_json,
        );

    let volume: Option<Expr> =
        from_json_string(volume_json);

    match (vector_field, volume) {
        | (Some(vf), Some(vol)) => {

            let result = gauss_theorem(
                &vf, &vol,
            );

            to_json_string(&result)
        },
        | _ => std::ptr::null_mut(),
    }
}

/// Represents Stokes' theorem (JSON)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_stokes_theorem(
    vector_field_json : *const std::os::raw::c_char,
    surface_json : *const std::os::raw::c_char,
) -> *mut std::os::raw::c_char {

    let vector_field: Option<Vector> =
        from_json_string(
            vector_field_json,
        );

    let surface: Option<Expr> =
        from_json_string(surface_json);

    match (
        vector_field,
        surface,
    ) {
        | (Some(vf), Some(surf)) => {

            let result = stokes_theorem(
                &vf, &surf,
            );

            to_json_string(&result)
        },
        | _ => std::ptr::null_mut(),
    }
}

/// Represents Green's theorem (JSON)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_greens_theorem(
    p_json: *const std::os::raw::c_char,
    q_json: *const std::os::raw::c_char,
    domain_json : *const std::os::raw::c_char,
) -> *mut std::os::raw::c_char {

    let p: Option<Expr> =
        from_json_string(p_json);

    let q: Option<Expr> =
        from_json_string(q_json);

    let domain: Option<Expr> =
        from_json_string(domain_json);

    match (p, q, domain) {
        | (
            Some(p_expr),
            Some(q_expr),
            Some(d),
        ) => {

            let result = greens_theorem(
                &p_expr,
                &q_expr,
                &d,
            );

            to_json_string(&result)
        },
        | _ => std::ptr::null_mut(),
    }
}
