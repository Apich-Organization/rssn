use std::collections::HashMap;

use rssn::symbolic::core::Expr;
use rssn::symbolic::group_theory::*;

#[test]

fn test_cyclic_group_c3() {

    // Elements: e=0, a=1, a^2=2 (using integers for simplicity)
    let e = GroupElement(
        Expr::new_constant(0.0),
    );

    let a = GroupElement(
        Expr::new_constant(1.0),
    );

    let b = GroupElement(
        Expr::new_constant(2.0),
    ); // a^2

    let elements = vec![
        e.clone(),
        a.clone(),
        b.clone(),
    ];

    let mut table = HashMap::new();

    // e * x = x * e = x
    for x in &elements {

        table.insert(
            (e.clone(), x.clone()),
            x.clone(),
        );

        table.insert(
            (x.clone(), e.clone()),
            x.clone(),
        );
    }

    // a * a = b
    table.insert(
        (a.clone(), a.clone()),
        b.clone(),
    );

    // a * b = e
    table.insert(
        (a.clone(), b.clone()),
        e.clone(),
    );

    // b * a = e
    table.insert(
        (b.clone(), a.clone()),
        e.clone(),
    );

    // b * b = a
    table.insert(
        (b.clone(), b.clone()),
        a.clone(),
    );

    let group = Group::new(
        elements,
        table,
        e.clone(),
    );

    // Test properties
    assert!(group.is_abelian());

    // Orders
    assert_eq!(
        group.element_order(&e),
        Some(1)
    );

    assert_eq!(
        group.element_order(&a),
        Some(3)
    );

    assert_eq!(
        group.element_order(&b),
        Some(3)
    );

    // Center (should be whole group for abelian)
    let center = group.center();

    assert_eq!(center.len(), 3);

    // Conjugacy classes (should be singletons for abelian)
    let classes =
        group.conjugacy_classes();

    assert_eq!(classes.len(), 3);

    for class in classes {

        assert_eq!(class.len(), 1);
    }
}

#[test]

fn test_representation_c3() {

    // Trivial representation: rho(g) = [1]
    let e = GroupElement(
        Expr::new_constant(0.0),
    );

    let a = GroupElement(
        Expr::new_constant(1.0),
    );

    let b = GroupElement(
        Expr::new_constant(2.0),
    );

    let elements = vec![
        e.clone(),
        a.clone(),
        b.clone(),
    ];

    // Recreate group (needed for is_valid check)
    let mut table = HashMap::new();

    for x in &elements {

        table.insert(
            (e.clone(), x.clone()),
            x.clone(),
        );

        table.insert(
            (x.clone(), e.clone()),
            x.clone(),
        );
    }

    table.insert(
        (a.clone(), a.clone()),
        b.clone(),
    );

    table.insert(
        (a.clone(), b.clone()),
        e.clone(),
    );

    table.insert(
        (b.clone(), a.clone()),
        e.clone(),
    );

    table.insert(
        (b.clone(), b.clone()),
        a.clone(),
    );

    let group = Group::new(
        elements.clone(),
        table,
        e.clone(),
    );

    let mut matrices = HashMap::new();

    let one_matrix =
        Expr::Matrix(vec![vec![
            Expr::new_constant(1.0),
        ]]);

    matrices.insert(
        e.clone(),
        one_matrix.clone(),
    );

    matrices.insert(
        a.clone(),
        one_matrix.clone(),
    );

    matrices.insert(
        b.clone(),
        one_matrix.clone(),
    );

    let rep = Representation::new(
        elements,
        matrices,
    );

    assert!(rep.is_valid(&group));

    let chars = character(&rep);

    // Character of trivial rep is always dimension (1)
    for (_, val) in chars {

        assert_eq!(
            val,
            Expr::new_constant(1.0)
        );
    }
}
