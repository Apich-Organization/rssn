use rssn::symbolic::core::Expr;
use rssn::symbolic::solid_state_physics::*;
use rssn::symbolic::vector::Vector;

fn main() {

    println!(
        "--- Solid State Physics Demo \
         ---"
    );

    // 1. Crystal Lattice and Reciprocal Lattice
    let a1 = Vector::new(
        Expr::new_variable("a"),
        Expr::Constant(0.0),
        Expr::Constant(0.0),
    );

    let a2 = Vector::new(
        Expr::Constant(0.0),
        Expr::new_variable("a"),
        Expr::Constant(0.0),
    );

    let a3 = Vector::new(
        Expr::Constant(0.0),
        Expr::Constant(0.0),
        Expr::new_variable("a"),
    );

    let simple_cubic =
        CrystalLattice::new(a1, a2, a3);

    println!(
        "Unit Cell Volume: {}",
        simple_cubic.volume()
    );

    let a2_cross_a3 = simple_cubic
        .a2
        .cross(&simple_cubic.a3);

    println!(
        "a2 x a3: ({}, {}, {})",
        a2_cross_a3.x,
        a2_cross_a3.y,
        a2_cross_a3.z
    );

    let (b1, b2, b3) = simple_cubic
        .reciprocal_lattice_vectors();

    println!(
        "Reciprocal b1.x: {}",
        b1.x
    );

    println!(
        "Reciprocal b2: {}",
        b2.y
    );

    println!(
        "Reciprocal b3: {}",
        b3.z
    );

    // 2. Fermi Energy
    let n = Expr::new_variable("n");

    let m_star =
        Expr::new_variable("m_star");

    let ef =
        fermi_energy_3d(&n, &m_star);

    println!(
        "Fermi Energy (3D): {}",
        ef
    );

    // 3. Density of States
    let energy =
        Expr::new_variable("E");

    let volume =
        Expr::new_variable("V");

    let dos = density_of_states_3d(
        &energy, &m_star, &volume,
    );

    println!(
        "Density of States (3D): {}",
        dos
    );

    // 4. Drude Conductivity
    let e = Expr::new_variable("e");

    let tau = Expr::new_variable("tau");

    let sigma = drude_conductivity(
        &n, &e, &tau, &m_star,
    );

    println!(
        "Drude Conductivity: {}",
        sigma
    );
}
