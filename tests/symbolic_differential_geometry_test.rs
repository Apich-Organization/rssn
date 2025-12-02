use rssn::symbolic::core::Expr;
use rssn::symbolic::differential_geometry::*;
use rssn::symbolic::vector::Vector;
use std::collections::BTreeMap;

#[test]
fn test_exterior_derivative_0form() {
    // 0-form: f(x,y) = x^2 + y
    // df = 2x dx + dy
    let mut terms = BTreeMap::new();
    let x = Expr::Variable("x".to_string());
    let y = Expr::Variable("y".to_string());
    let f = Expr::new_add(Expr::new_pow(x.clone(), Expr::Constant(2.0)), y.clone());

    // 0-form is represented with blade 0
    terms.insert(0, f);
    let form = DifferentialForm { terms };

    let d_form = exterior_derivative(&form, &["x", "y"]);

    // Should have terms for dx (blade 1) and dy (blade 2)
    assert!(d_form.terms.contains_key(&1)); // dx
    assert!(d_form.terms.contains_key(&2)); // dy
}

#[test]
fn test_wedge_product() {
    // dx ^ dy
    let mut form1_terms = BTreeMap::new();
    form1_terms.insert(1, Expr::Constant(1.0)); // dx (blade 1 = 0b01)
    let form1 = DifferentialForm { terms: form1_terms };

    let mut form2_terms = BTreeMap::new();
    form2_terms.insert(2, Expr::Constant(1.0)); // dy (blade 2 = 0b10)
    let form2 = DifferentialForm { terms: form2_terms };

    let wedge = wedge_product(&form1, &form2);

    // dx ^ dy should give blade 3 (0b11)
    assert!(wedge.terms.contains_key(&3));
}

#[test]
fn test_wedge_product_antisymmetry() {
    // dx ^ dx = 0
    let mut form_terms = BTreeMap::new();
    form_terms.insert(1, Expr::Constant(1.0)); // dx
    let form = DifferentialForm { terms: form_terms };

    let wedge = wedge_product(&form, &form);

    // Should be empty or zero
    assert!(
        wedge.terms.is_empty()
            || wedge
                .terms
                .values()
                .all(|v| rssn::symbolic::simplify::is_zero(v))
    );
}

#[test]
fn test_boundary() {
    let domain = Expr::Variable("M".to_string());
    let boundary_expr = boundary(&domain);

    // Should create a Boundary expression
    if let Expr::Boundary(inner) = boundary_expr {
        assert_eq!(*inner, domain);
    } else {
        panic!("Expected Boundary expression");
    }
}

#[test]
fn test_generalized_stokes_theorem() {
    // Create a simple 1-form
    let mut terms = BTreeMap::new();
    terms.insert(1, Expr::Variable("f".to_string())); // f dx
    let omega = DifferentialForm { terms };

    let manifold = Expr::Variable("M".to_string());

    let theorem = generalized_stokes_theorem(&omega, &manifold, &["x", "y"]);

    // Should create an equation
    if let Expr::Eq(_, _) = theorem {
        // Success
    } else {
        panic!("Expected Eq expression");
    }
}

#[test]
fn test_gauss_theorem() {
    let vector_field = Vector {
        x: Expr::Variable("Fx".to_string()),
        y: Expr::Variable("Fy".to_string()),
        z: Expr::Variable("Fz".to_string()),
    };

    let volume = Expr::Variable("V".to_string());

    let theorem = gauss_theorem(&vector_field, &volume);

    // Should create an equation
    if let Expr::Eq(_, _) = theorem {
        // Success
    } else {
        panic!("Expected Eq expression");
    }
}

#[test]
fn test_stokes_theorem() {
    let vector_field = Vector {
        x: Expr::Variable("Fx".to_string()),
        y: Expr::Variable("Fy".to_string()),
        z: Expr::Variable("Fz".to_string()),
    };

    let surface = Expr::Variable("S".to_string());

    let theorem = stokes_theorem(&vector_field, &surface);

    // Should create an equation
    if let Expr::Eq(_, _) = theorem {
        // Success
    } else {
        panic!("Expected Eq expression");
    }
}

#[test]
fn test_greens_theorem() {
    let p = Expr::Variable("P".to_string());
    let q = Expr::Variable("Q".to_string());
    let domain = Expr::Variable("D".to_string());

    let theorem = greens_theorem(&p, &q, &domain);

    // Should create an equation
    if let Expr::Eq(_, _) = theorem {
        // Success
    } else {
        panic!("Expected Eq expression");
    }
}
