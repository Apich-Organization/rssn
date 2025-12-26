use num_bigint::BigInt;
use rssn::symbolic::core::Expr;
use rssn::symbolic::topology::*;

#[test]

fn test_simplex_symbolic_boundary() {

    // 1-simplex (edge) [0, 1]
    let simplex = Simplex::new(&[0, 1]);

    let (faces, coeffs) = simplex.symbolic_boundary();

    // Boundary should be [1] - [0]
    assert_eq!(faces.len(), 2);

    assert_eq!(coeffs.len(), 2);

    // i=0: remove 0 -> face [1], coeff +1
    // i=1: remove 1 -> face [0], coeff -1

    assert_eq!(faces[0].0.iter().next(), Some(&1));

    assert_eq!(coeffs[0], Expr::BigInt(BigInt::from(1)));

    assert_eq!(faces[1].0.iter().next(), Some(&0));

    assert_eq!(coeffs[1], Expr::BigInt(BigInt::from(-1)));
}

#[test]

fn test_symbolic_boundary_matrix() {

    // Triangle [0, 1, 2]
    let mut complex = SimplicialComplex::new();

    complex.add_simplex(&[0, 1, 2]);

    // 2-simplices: [0, 1, 2] (1 simplex)
    // 1-simplices: [0, 1], [0, 2], [1, 2] (3 simplices)

    let boundary_matrix = complex.get_symbolic_boundary_matrix(2).unwrap();

    if let Expr::Matrix(rows) = boundary_matrix {

        assert_eq!(rows.len(), 3); // 3 edges
        assert_eq!(rows[0].len(), 1); // 1 triangle

        // Check entries are 1 or -1
        let mut non_zero_count = 0;

        for row in rows {

            if !rssn::symbolic::simplify::is_zero(&row[0]) {

                non_zero_count += 1;
            }
        }

        assert_eq!(non_zero_count, 3);
    } else {

        panic!("Expected matrix");
    }
}

#[test]

fn test_apply_symbolic_boundary_operator() {

    // Triangle [0, 1, 2]
    let mut complex = SimplicialComplex::new();

    complex.add_simplex(&[0, 1, 2]);

    let triangle = Simplex::new(&[0, 1, 2]);

    let mut chain = SymbolicChain::new(2);

    // Chain: a * [0, 1, 2]
    let a = Expr::Variable("a".to_string());

    chain.add_term(triangle, a.clone()).unwrap();

    let boundary_chain = complex.apply_symbolic_boundary_operator(&chain).unwrap();

    assert_eq!(boundary_chain.dimension, 1);

    assert_eq!(boundary_chain.terms.len(), 3);

    // Boundary should be a*([1, 2] - [0, 2] + [0, 1])
    // Coefficients should be a, -a, a

    for coeff in boundary_chain.terms.values() {

        // coeff should be a or -a
        // simplify(coeff^2 - a^2) should be 0
        let is_a = rssn::symbolic::simplify_dag::simplify(&Expr::new_sub(coeff.clone(), a.clone()));

        let is_neg_a =
            rssn::symbolic::simplify_dag::simplify(&Expr::new_add(coeff.clone(), a.clone()));

        let is_a_zero = rssn::symbolic::simplify::is_zero(&is_a);

        let is_neg_a_zero = rssn::symbolic::simplify::is_zero(&is_neg_a);

        assert!(is_a_zero || is_neg_a_zero);
    }
}
