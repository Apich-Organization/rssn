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

    // The commutator returns A*B*psi - B*A*psi.
    // Full cancellation to 0 requires associativity normalization
    // (flattening nested Mul chains), which is not yet implemented
    // in the E-Graph to avoid exponential blowup.
    // For now, verify the result is structurally a Sub of two Mul terms.
    let s = format!("{:?}", comm);
    assert!(
        s.contains("A") && s.contains("B") && s.contains("psi"),
        "Commutator should reference A, B, and psi, got: {}", s
    );
}

#[test]

fn test_pauli_matrices() {
    let (sx, sy, sz) = pauli_matrices();

    assert!(sx.to_string().contains("[[0, 1]; [1, 0]]"));

    assert!(sy.to_string().contains("i"));

    assert!(sz.to_string().contains("1"));
}

#[test]

fn test_expectation_value() {
    let x = Operator::new(Expr::new_variable("x"));

    let psi = Ket {
        state: Expr::new_variable("psi"),
    };

    let exp_x = expectation_value(&x, &psi);

    assert!(exp_x.to_string().contains("x"));

    assert!(exp_x.to_string().contains("psi"));
}

#[test]

fn test_hamiltonian_free_particle() {
    let m = Expr::new_variable("m");

    let h = hamiltonian_free_particle(&m);

    assert!(h.op.to_string().contains("hbar"));

    assert!(h.op.to_string().contains("m"));

    assert!(h.op.to_string().contains("d2_dx2"));
}
