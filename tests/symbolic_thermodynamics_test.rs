use rssn::symbolic::core::Expr;
use rssn::symbolic::thermodynamics::*;

#[test]

fn test_ideal_gas_law() {

    let p = Expr::new_variable("P");

    let v = Expr::new_variable("V");

    let n = Expr::new_variable("n");

    let r = Expr::new_variable("R");

    let t = Expr::new_variable("T");

    let eq = ideal_gas_law(
        &p, &v, &n, &r, &t,
    );

    let eq_str = format!("{:?}", eq);

    assert!(eq_str.contains("P"));

    assert!(eq_str.contains("V"));

    assert!(eq_str.contains("n"));

    assert!(eq_str.contains("T"));
}

#[test]

fn test_carnot_efficiency() {

    let tc = Expr::new_variable("Tc");

    let th = Expr::new_variable("Th");

    let eta =
        carnot_efficiency(&tc, &th);

    let eta_str = format!("{:?}", eta);

    // 1 - Tc/Th
    assert!(eta_str.contains("1"));

    assert!(eta_str.contains("Tc"));

    assert!(eta_str.contains("Th"));
}

#[test]

fn test_boltzmann_entropy() {

    let omega =
        Expr::new_variable("Omega");

    let s = boltzmann_entropy(&omega);

    let s_str = format!("{:?}", s);

    assert!(s_str.contains("k_B"));

    assert!(
        s_str.contains("ln(Omega)")
    );
}

#[test]

fn test_enthalpy_and_gibbs() {

    let u = Expr::new_variable("U");

    let p = Expr::new_variable("P");

    let v = Expr::new_variable("V");

    let t = Expr::new_variable("T");

    let s = Expr::new_variable("S");

    let h = enthalpy(&u, &p, &v);

    let g =
        gibbs_free_energy(&h, &t, &s);

    let g_str = format!("{:?}", g);

    assert!(g_str.contains("U"));

    assert!(g_str.contains("P"));

    assert!(g_str.contains("V"));

    assert!(g_str.contains("T"));

    assert!(g_str.contains("S"));
}

#[test]

fn test_maxwell_relation() {

    // A = U - TS
    // For an ideal gas state, we might have some expression for A.
    // Let's just test the symmetry of second derivatives.
    let t = Expr::new_variable("T");

    let v = Expr::new_variable("V");

    // A(T, V) = some function
    let a = Expr::new_mul(
        t.clone(),
        Expr::new_log(v.clone()),
    );

    let diff = verify_maxwell_relation_helmholtz(&a, "T", "V");

    assert_eq!(
        diff,
        Expr::new_constant(0.0)
    );
}
