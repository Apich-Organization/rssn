use rssn::symbolic::core::Expr;
use rssn::symbolic::solid_state_physics::*;
use rssn::symbolic::vector::Vector;

#[test]

fn test_crystal_lattice_volume() {

    let a1 = Vector::new(
        Expr::new_constant(1.0),
        Expr::new_constant(0.0),
        Expr::new_constant(0.0),
    );

    let a2 = Vector::new(
        Expr::new_constant(0.0),
        Expr::new_constant(1.0),
        Expr::new_constant(0.0),
    );

    let a3 = Vector::new(
        Expr::new_constant(0.0),
        Expr::new_constant(0.0),
        Expr::new_constant(1.0),
    );

    let lattice =
        CrystalLattice::new(a1, a2, a3);

    let vol = lattice.volume();

    assert_eq!(
        vol,
        Expr::new_constant(1.0)
    );
}

#[test]

fn test_reciprocal_lattice_vectors() {

    let a1 = Vector::new(
        Expr::new_constant(1.0),
        Expr::new_constant(0.0),
        Expr::new_constant(0.0),
    );

    let a2 = Vector::new(
        Expr::new_constant(0.0),
        Expr::new_constant(1.0),
        Expr::new_constant(0.0),
    );

    let a3 = Vector::new(
        Expr::new_constant(0.0),
        Expr::new_constant(0.0),
        Expr::new_constant(1.0),
    );

    let lattice =
        CrystalLattice::new(a1, a2, a3);

    let (b1, b2, b3) = lattice
        .reciprocal_lattice_vectors();

    // For unit cube, reciprocal vectors should be 2pi * unit vectors
    // 1.0 / 1.0 = 1.0, so it should be exactly 2 * pi
    assert!(
        b1.x.to_string()
            .contains("pi")
    );

    assert!(
        b1.x.to_string()
            .contains("2")
    );

    assert!(
        b2.y.to_string()
            .contains("pi")
    );

    assert!(
        b3.z.to_string()
            .contains("pi")
    );
}

#[test]

fn test_debye_frequency() {

    let v_s = Expr::new_variable("v_s");

    let n_dense =
        Expr::new_variable("n");

    let omega_d =
        debye_frequency(&v_s, &n_dense);

    assert!(
        omega_d
            .to_string()
            .contains("v_s")
    );

    assert!(
        omega_d
            .to_string()
            .contains("n")
    );
}

#[test]

fn test_plasma_frequency() {

    let n = Expr::new_variable("n");

    let e = Expr::new_variable("e");

    let epsilon_0 =
        Expr::new_variable("epsilon_0");

    let m = Expr::new_variable("m");

    let omega_p = plasma_frequency(
        &n,
        &e,
        &epsilon_0,
        &m,
    );

    assert!(
        omega_p
            .to_string()
            .contains("n")
    );

    assert!(
        omega_p
            .to_string()
            .contains("e")
    );

    assert!(
        omega_p
            .to_string()
            .contains("m")
    );

    assert!(
        omega_p
            .to_string()
            .contains("epsilon_0")
    );
}

#[test]

fn test_einstein_heat_capacity() {

    let n = Expr::new_variable("N");

    let einstein_temp =
        Expr::new_variable("Theta_E");

    let temp = Expr::new_variable("T");

    let cv = einstein_heat_capacity(
        &n,
        &einstein_temp,
        &temp,
    );

    assert!(
        cv.to_string()
            .contains("N")
    );

    assert!(
        cv.to_string()
            .contains("Theta_E")
    );

    assert!(
        cv.to_string()
            .contains("T")
    );
}

#[test]

fn test_london_penetration_depth() {

    let mass = Expr::new_variable("m");

    let mu_0 =
        Expr::new_variable("mu_0");

    let n_s = Expr::new_variable("n_s");

    let e = Expr::new_variable("e");

    let lambda =
        london_penetration_depth(
            &mass, &mu_0, &n_s, &e,
        );

    assert!(
        lambda
            .to_string()
            .contains("n_s")
    );

    assert!(
        lambda
            .to_string()
            .contains("mu_0")
    );
}

#[test]

fn test_hall_coefficient() {

    let n = Expr::new_variable("n");

    let q = Expr::new_variable("q");

    let rh = hall_coefficient(&n, &q);

    assert!(
        rh.to_string()
            .contains("n")
    );

    assert!(
        rh.to_string()
            .contains("q")
    );
}

#[test]

fn test_bloch_theorem() {

    let k = Vector::new(
        Expr::new_variable("k_x"),
        Expr::new_constant(0.0),
        Expr::new_constant(0.0),
    );

    let r = Vector::new(
        Expr::new_variable("x"),
        Expr::new_constant(0.0),
        Expr::new_constant(0.0),
    );

    let u = Expr::new_variable("u");

    let psi = bloch_theorem(&k, &r, &u);

    assert!(
        psi.to_string()
            .contains("exp")
    );

    assert!(
        psi.to_string()
            .contains("u")
    );
}

#[test]

fn test_energy_band() {

    let k_mag = Expr::new_variable("k");

    let m_star =
        Expr::new_variable("m_star");

    let e0 = Expr::new_variable("E_0");

    let energy = energy_band(
        &k_mag,
        &m_star,
        &e0,
    );

    assert!(
        energy
            .to_string()
            .contains("hbar")
    );

    assert!(
        energy
            .to_string()
            .contains("k")
    );
}

#[test]

fn test_fermi_energy_3d() {

    let n = Expr::new_variable("n");

    let m_star =
        Expr::new_variable("m_star");

    let ef =
        fermi_energy_3d(&n, &m_star);

    let ef_str = ef.to_string();

    assert!(ef_str.contains("hbar"));

    assert!(ef_str.contains("pi"));

    assert!(ef_str.contains("n"));
}

#[test]

fn test_drude_conductivity() {

    let n = Expr::new_variable("n");

    let e = Expr::new_variable("e");

    let tau = Expr::new_variable("tau");

    let m_star =
        Expr::new_variable("m_star");

    let sigma = drude_conductivity(
        &n,
        &e,
        &tau,
        &m_star,
    );

    let sigma_str = sigma.to_string();

    assert!(sigma_str.contains("n"));

    assert!(sigma_str.contains("e"));

    assert!(sigma_str.contains("tau"));

    assert!(
        sigma_str.contains("m_star")
    );
}
