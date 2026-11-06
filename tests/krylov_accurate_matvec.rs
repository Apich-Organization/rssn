use rssn::kernels::krylov::Csr;

#[test]
fn sparse_matvec_survives_intermediate_overflow_and_cancellation() {
    let matrix = Csr::from_triplets(
        2,
        2,
        &[(0, 0, 1e308), (0, 1, 9e307), (1, 0, 9e307), (1, 1, 1e308)],
    );
    let mut out = [0.0; 2];
    matrix.matvec(&[2.0, -2.0], &mut out);
    assert!((out[0] / 2e307 - 1.0).abs() < 1e-15);
    assert!((out[1] / -2e307 - 1.0).abs() < 1e-15);
}

#[test]
fn sparse_matvec_preserves_small_cross_terms() {
    let matrix = Csr::from_triplets(1, 2, &[(0, 0, 1e300), (0, 1, 1e-300)]);
    let mut out = [0.0];
    matrix.matvec(&[1e-300, 1e300], &mut out);
    assert_eq!(out[0], 2.0);
}
