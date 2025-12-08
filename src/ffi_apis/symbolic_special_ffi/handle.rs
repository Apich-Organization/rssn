use crate::symbolic::special;
use std::os::raw::c_double;

#[no_mangle]
pub extern "C" fn rssn_gamma_numerical(x: c_double) -> c_double {
    special::gamma_numerical(x)
}

#[no_mangle]
pub extern "C" fn rssn_ln_gamma_numerical(x: c_double) -> c_double {
    special::ln_gamma_numerical(x)
}

#[no_mangle]
pub extern "C" fn rssn_beta_numerical(a: c_double, b: c_double) -> c_double {
    special::beta_numerical(a, b)
}

#[no_mangle]
pub extern "C" fn rssn_ln_beta_numerical(a: c_double, b: c_double) -> c_double {
    special::ln_beta_numerical(a, b)
}

#[no_mangle]
pub extern "C" fn rssn_erf_numerical(x: c_double) -> c_double {
    special::erf_numerical(x)
}

#[no_mangle]
pub extern "C" fn rssn_erfc_numerical(x: c_double) -> c_double {
    special::erfc_numerical(x)
}
