use std::os::raw::c_char;

use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::quantum_field_theory;

/// Computes the Dirac adjoint using JSON.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_dirac_adjoint(
    psi_json : *const c_char
) -> *mut c_char {

    let psi : Option<Expr> =
        from_json_string(psi_json);

    if let Some(psi) = psi {

        to_json_string(&quantum_field_theory::dirac_adjoint(&psi))
    } else {

        std::ptr::null_mut()
    }
}

/// Computes the Feynman slash using JSON.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_feynman_slash(
    v_mu_json : *const c_char
) -> *mut c_char {

    let v_mu : Option<Expr> =
        from_json_string(v_mu_json);

    if let Some(v_mu) = v_mu {

        to_json_string(&quantum_field_theory::feynman_slash(&v_mu))
    } else {

        std::ptr::null_mut()
    }
}

/// Lagrangian density for a scalar field using JSON.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_scalar_field_lagrangian(
    phi_json : *const c_char,
    m_json : *const c_char,
) -> *mut c_char {

    let phi : Option<Expr> =
        from_json_string(phi_json);

    let m : Option<Expr> =
        from_json_string(m_json);

    if let (Some(phi), Some(m)) =
        (phi, m)
    {

        to_json_string(&quantum_field_theory::scalar_field_lagrangian(&phi, &m))
    } else {

        std::ptr::null_mut()
    }
}

/// Lagrangian density for QED using JSON.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_qed_lagrangian(
    psi_bar_json : *const c_char,
    psi_json : *const c_char,
    a_mu_json : *const c_char,
    m_json : *const c_char,
    e_json : *const c_char,
) -> *mut c_char {

    let psi_bar : Option<Expr> =
        from_json_string(psi_bar_json);

    let psi : Option<Expr> =
        from_json_string(psi_json);

    let a_mu : Option<Expr> =
        from_json_string(a_mu_json);

    let m : Option<Expr> =
        from_json_string(m_json);

    let e : Option<Expr> =
        from_json_string(e_json);

    if let (
        Some(psi_bar),
        Some(psi),
        Some(a_mu),
        Some(m),
        Some(e),
    ) = (
        psi_bar,
        psi,
        a_mu,
        m,
        e,
    ) {

        to_json_string(
            &quantum_field_theory::qed_lagrangian(
                &psi_bar,
                &psi,
                &a_mu,
                &m,
                &e,
            ),
        )
    } else {

        std::ptr::null_mut()
    }
}

/// Lagrangian density for QCD using JSON.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_qcd_lagrangian(
    psi_bar_json : *const c_char,
    psi_json : *const c_char,
    g_mu_json : *const c_char,
    m_json : *const c_char,
    gs_json : *const c_char,
) -> *mut c_char {

    let psi_bar : Option<Expr> =
        from_json_string(psi_bar_json);

    let psi : Option<Expr> =
        from_json_string(psi_json);

    let g_mu : Option<Expr> =
        from_json_string(g_mu_json);

    let m : Option<Expr> =
        from_json_string(m_json);

    let gs : Option<Expr> =
        from_json_string(gs_json);

    if let (
        Some(psi_bar),
        Some(psi),
        Some(g_mu),
        Some(m),
        Some(gs),
    ) = (
        psi_bar,
        psi,
        g_mu,
        m,
        gs,
    ) {

        to_json_string(
            &quantum_field_theory::qcd_lagrangian(
                &psi_bar,
                &psi,
                &g_mu,
                &m,
                &gs,
            ),
        )
    } else {

        std::ptr::null_mut()
    }
}

/// Computes a propagator using JSON.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_qft_propagator(
    p_json : *const c_char,
    m_json : *const c_char,
    is_fermion : bool,
) -> *mut c_char {

    let p : Option<Expr> =
        from_json_string(p_json);

    let m : Option<Expr> =
        from_json_string(m_json);

    if let (Some(p), Some(m)) = (p, m) {

        to_json_string(&quantum_field_theory::propagator(&p, &m, is_fermion))
    } else {

        std::ptr::null_mut()
    }
}
