use rssn::symbolic::core::Expr;
use rssn::symbolic::relativity::*;

#[test]

fn test_lorentz_factor() {

    let v = Expr::new_variable("v");

    let gamma = lorentz_factor(&v);

    let gamma_str =
        format!("{:?}", gamma);

    println!(
        "Lorentz Factor: {}",
        gamma_str
    );

    assert!(gamma_str.contains("v^(2)"));

    assert!(gamma_str.contains("c"));
}

#[test]

fn test_mass_energy() {

    let m = Expr::new_variable("m");

    let e = mass_energy_equivalence(&m);

    let e_str = format!("{:?}", e);

    assert!(e_str.contains("m"));

    assert!(e_str.contains("c^(2)"));
}

#[test]

fn test_schwarzschild_radius() {

    let m = Expr::new_variable("M");

    let rs = schwarzschild_radius(&m);

    let rs_str = format!("{:?}", rs);

    println!(
        "Schwarzschild Radius: {}",
        rs_str
    );

    assert!(rs_str.contains("G"));

    assert!(rs_str.contains("M"));

    assert!(rs_str.contains("c"));
}

#[test]

fn test_doppler_effect() {

    let f = Expr::new_variable("f");

    let v = Expr::new_variable("v");

    let f_obs = doppler_effect(&f, &v);

    let f_obs_str =
        format!("{:?}", f_obs);

    assert!(f_obs_str.contains("f"));

    assert!(f_obs_str.contains("v"));

    assert!(f_obs_str.contains("c"));

    assert!(f_obs_str.contains("0.5"));
}

#[test]

fn test_einstein_tensor() {

    let r_mu_nu =
        Expr::new_variable("R_mn");

    let r = Expr::new_variable("R");

    let g_mu_nu =
        Expr::new_variable("g_mn");

    let g_tensor = einstein_tensor(
        &r_mu_nu,
        &r,
        &g_mu_nu,
    );

    let g_str =
        format!("{:?}", g_tensor);

    assert!(g_str.contains("R_mn"));

    assert!(g_str.contains("0.5"));

    assert!(g_str.contains("R"));

    assert!(g_str.contains("g_mn"));
}
