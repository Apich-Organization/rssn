use rssn::prelude::*;
use rssn::prelude::numerical::*;
use proptest::prelude::*;

#[test]
fn test_gf_p_basic() {
    let p = 11;
    let a = numerical_PrimeFieldElement::new(7, p);
    let b = numerical_PrimeFieldElement::new(5, p);
    
    // 7 + 5 = 12 ≡ 1 mod 11
    assert_eq!((a + b).value, 1);
    // 7 * 5 = 35 ≡ 2 mod 11
    assert_eq!((a * b).value, 2);
    // 7 - 5 = 2
    assert_eq!((a - b).value, 2);
    // 5 - 7 = -2 ≡ 9 mod 11
    assert_eq!((b - a).value, 9);
}

#[test]
fn test_gf_p_inverse() {
    let p = 11;
    let a = numerical_PrimeFieldElement::new(7, p);
    let inv = a.inverse().expect("7 is invertible mod 11");
    // 7 * 8 = 56 ≡ 1 mod 11
    assert_eq!(inv.value, 8);
    assert_eq!((a * inv).value, 1);
    
    let zero = numerical_PrimeFieldElement::new(0, p);
    assert!(zero.inverse().is_none());
}

#[test]
fn test_gf_p_pow() {
    let p = 11;
    let a = numerical_PrimeFieldElement::new(2, p);
    // 2^5 = 32 ≡ 10 mod 11
    assert_eq!(a.pow(5).value, 10);
    // Fermat's Little Theorem: a^(p-1) ≡ 1 mod p
    assert_eq!(a.pow(10).value, 1);
}

#[test]
fn test_gf256_basic() {
    // Addition is XOR
    assert_eq!(numerical_gf256_add(0x57, 0x83), 0x57 ^ 0x83);
    
    // Multiplication and inverse
    let a = 0x57;
    let inv = numerical_gf256_inv(a).expect("Non-zero element is invertible");
    assert_eq!(numerical_gf256_mul(a, inv), 1);
}

#[test]
fn test_gf256_pow() {
    let a = 2;
    // 2^1 = 2, 2^2 = 4, 2^3 = 8
    assert_eq!(numerical_gf256_pow(a, 0), 1);
    assert_eq!(numerical_gf256_pow(a, 1), 2);
    assert_eq!(numerical_gf256_pow(a, 2), 4);
    assert_eq!(numerical_gf256_pow(a, 3), 8);
    // 2^255 ≡ 1 in GF(2^8) generator order normally? Depends on generator.
    // Our generator x=1 is used to fill table, let's check a*inv
}

proptest! {
    #[test]
    fn prop_pfe_mul_inv(v in 1..100u64, p_idx in 0..10usize) {
        let primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29];
        let p = primes[p_idx];
        let val = v % p;
        if val != 0 {
            let a = numerical_PrimeFieldElement::new(val, p);
            let inv = a.inverse().unwrap();
            assert_eq!((a * inv).value, 1);
        }
    }

    #[test]
    fn prop_gf256_mul_inv(a in 1..255u8) {
        let inv = numerical_gf256_inv(a).unwrap();
        assert_eq!(numerical_gf256_mul(a, inv), 1);
    }
}
