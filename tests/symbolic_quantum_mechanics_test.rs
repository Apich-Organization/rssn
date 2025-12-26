use rssn::symbolic::core::Expr;
use rssn::symbolic::quantum_mechanics::*;

#[test]

fn test_bra_ket() {

    let psi = Ket {
        state: Expr::new_variable("psi"),
    };

    let phi = Bra {
        state: Expr::new_variable("phi"),
    };

    let inner = bra_ket(&phi, &psi);

    let inner_str = inner.to_string();

    assert!(inner_str.contains("integral"));

    assert!(inner_str.contains("phi"));

    assert!(inner_str.contains("psi"));
}

#[test]

fn test_commutator() {

    let a = Operator::new(Expr::new_variable("A"));

    let b = Operator::new(Expr::new_variable("B"));

    let psi = Ket {
        state: Expr::new_variable("psi"),
    };

    let comm = commutator(&a, &b, &psi);

    // In this symbolic representation, A and B commute as simple variables, so A*B*psi - B*A*psi = 0.
    assert_eq!(comm, Expr::Constant(0.0));
}

#[test]

fn test_pauli_matrices() {

    let (sx, sy, sz) = pauli_matrices();

    assert!(sx
        .to_string()
        .contains("[[0, 1]; [1, 0]]"));

    assert!(sy
        .to_string()
        .contains("i"));

    assert!(sz
        .to_string()
        .contains("1"));
}

#[test]

fn test_expectation_value() {

    let x = Operator::new(Expr::new_variable("x"));

    let psi = Ket {
        state: Expr::new_variable("psi"),
    };

    let exp_x = expectation_value(&x, &psi);

    assert!(exp_x
        .to_string()
        .contains("x"));

    assert!(exp_x
        .to_string()
        .contains("psi"));
}

#[test]

fn test_hamiltonian_free_particle() {

    let m = Expr::new_variable("m");

    let h = hamiltonian_free_particle(&m);

    assert!(h
        .op
        .to_string()
        .contains("hbar"));

    assert!(h
        .op
        .to_string()
        .contains("m"));

    assert!(h
        .op
        .to_string()
        .contains("d2_dx2"));
}
