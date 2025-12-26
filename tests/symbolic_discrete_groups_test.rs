use rssn::symbolic::discrete_groups::*;
use rssn::symbolic::group_theory::*;

#[test]

fn test_cyclic_group() {

    let c3 = cyclic_group(3);

    assert_eq!(c3.elements.len(), 3);

    assert!(c3.is_abelian());

    // Check orders
    // e has order 1
    assert_eq!(c3.element_order(&c3.identity), Some(1));

    // g and g^2 have order 3
    for el in &c3.elements {

        if *el != c3.identity {

            assert_eq!(c3.element_order(el), Some(3));
        }
    }
}

#[test]

fn test_dihedral_group() {

    let d3 = dihedral_group(3); // Order 2*3 = 6
    assert_eq!(d3.elements.len(), 6);

    // D3 is isomorphic to S3, non-abelian
    assert!(!d3.is_abelian());

    let d4 = dihedral_group(4); // Order 8
    assert_eq!(d4.elements.len(), 8);

    assert!(!d4.is_abelian());
}

#[test]

fn test_symmetric_group() {

    let s3 = symmetric_group(3).unwrap();

    assert_eq!(s3.elements.len(), 6); // 3! = 6
    assert!(!s3.is_abelian());

    let s2 = symmetric_group(2).unwrap();

    assert_eq!(s2.elements.len(), 2); // 2! = 2
    assert!(s2.is_abelian()); // S2 is abelian
}

#[test]

fn test_klein_four_group() {

    let v4 = klein_four_group();

    assert_eq!(v4.elements.len(), 4);

    assert!(v4.is_abelian());

    // Every element (except identity) has order 2
    for el in &v4.elements {

        if *el != v4.identity {

            assert_eq!(v4.element_order(el), Some(2));
        }
    }
}
