use crate::ffi_apis::common::{BincodeBuffer, from_bincode_buffer, to_bincode_buffer};
use crate::symbolic::core::Expr;
use crate::symbolic::differential_geometry::{DifferentialForm, exterior_derivative, wedge_product, boundary, generalized_stokes_theorem, gauss_theorem, stokes_theorem, greens_theorem};
use crate::symbolic::vector::Vector;

/// Computes the exterior derivative of a differential form (Bincode)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_exterior_derivative(
    form_buf: BincodeBuffer,
    vars_buf: BincodeBuffer,
) -> BincodeBuffer {

    let form: Option<DifferentialForm> =
        from_bincode_buffer(&form_buf);

    let vars: Option<Vec<String>> =
        from_bincode_buffer(&vars_buf);

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

            to_bincode_buffer(&result)
        },
        | _ => BincodeBuffer::empty(),
    }
}

/// Computes the wedge product of two differential forms (Bincode)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_wedge_product(
    form1_buf: BincodeBuffer,
    form2_buf: BincodeBuffer,
) -> BincodeBuffer {

    let form1: Option<
        DifferentialForm,
    > = from_bincode_buffer(&form1_buf);

    let form2: Option<
        DifferentialForm,
    > = from_bincode_buffer(&form2_buf);

    match (form1, form2) {
        | (Some(f1), Some(f2)) => {

            let result =
                wedge_product(&f1, &f2);

            to_bincode_buffer(&result)
        },
        | _ => BincodeBuffer::empty(),
    }
}

/// Computes the boundary of a domain (Bincode)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_boundary(
    domain_buf: BincodeBuffer
) -> BincodeBuffer {

    let domain: Option<Expr> =
        from_bincode_buffer(
            &domain_buf,
        );

    if let Some(d) = domain {

        let result = boundary(&d);

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

/// Represents the generalized Stokes' theorem (Bincode)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_generalized_stokes_theorem(
    omega_buf: BincodeBuffer,
    manifold_buf: BincodeBuffer,
    vars_buf: BincodeBuffer,
) -> BincodeBuffer {

    let omega: Option<
        DifferentialForm,
    > = from_bincode_buffer(&omega_buf);

    let manifold: Option<Expr> =
        from_bincode_buffer(
            &manifold_buf,
        );

    let vars: Option<Vec<String>> =
        from_bincode_buffer(&vars_buf);

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

            to_bincode_buffer(&result)
        },
        | _ => BincodeBuffer::empty(),
    }
}

/// Represents Gauss's theorem (Bincode)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_gauss_theorem(
    vector_field_buf: BincodeBuffer,
    volume_buf: BincodeBuffer,
) -> BincodeBuffer {

    let vector_field: Option<Vector> =
        from_bincode_buffer(
            &vector_field_buf,
        );

    let volume: Option<Expr> =
        from_bincode_buffer(
            &volume_buf,
        );

    match (vector_field, volume) {
        | (Some(vf), Some(vol)) => {

            let result = gauss_theorem(
                &vf, &vol,
            );

            to_bincode_buffer(&result)
        },
        | _ => BincodeBuffer::empty(),
    }
}

/// Represents Stokes' theorem (Bincode)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_stokes_theorem(
    vector_field_buf: BincodeBuffer,
    surface_buf: BincodeBuffer,
) -> BincodeBuffer {

    let vector_field: Option<Vector> =
        from_bincode_buffer(
            &vector_field_buf,
        );

    let surface: Option<Expr> =
        from_bincode_buffer(
            &surface_buf,
        );

    match (
        vector_field,
        surface,
    ) {
        | (Some(vf), Some(surf)) => {

            let result = stokes_theorem(
                &vf, &surf,
            );

            to_bincode_buffer(&result)
        },
        | _ => BincodeBuffer::empty(),
    }
}

/// Represents Green's theorem (Bincode)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_greens_theorem(
    p_buf: BincodeBuffer,
    q_buf: BincodeBuffer,
    domain_buf: BincodeBuffer,
) -> BincodeBuffer {

    let p: Option<Expr> =
        from_bincode_buffer(&p_buf);

    let q: Option<Expr> =
        from_bincode_buffer(&q_buf);

    let domain: Option<Expr> =
        from_bincode_buffer(
            &domain_buf,
        );

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

            to_bincode_buffer(&result)
        },
        | _ => BincodeBuffer::empty(),
    }
}
