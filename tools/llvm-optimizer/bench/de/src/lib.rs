include!(env!("RSSN_DE_SOURCE"));

#[unsafe(no_mangle)]
#[inline(never)]
pub extern "C" fn rssn_de_small(
    z: f64,
    s: f64,
    a: f64,
) -> f64 {
    kernels::special_ext::small(z, s, a)
}

#[unsafe(no_mangle)]
#[inline(never)]
pub extern "C" fn rssn_de_positive(
    z: f64,
    s: f64,
    a: f64,
) -> f64 {
    kernels::special_ext::positive(z, s, a)
}

#[unsafe(no_mangle)]
#[inline(never)]
pub extern "C" fn rssn_de_negative(
    z: f64,
    s: f64,
    a: f64,
) -> f64 {
    kernels::special_ext::negative(z, s, a)
}

#[unsafe(no_mangle)]
#[inline(never)]
pub extern "C" fn rssn_de_appell(
    a: f64,
    b1: f64,
    b2: f64,
    c: f64,
    x: f64,
    y: f64,
) -> f64 {
    kernels::special_ext::appell(a, b1, b2, c, x, y)
}
