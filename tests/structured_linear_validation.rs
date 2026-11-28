//! Every structured PDE driver uses the same checked native linear-system boundary.
use rssn::kernels::krylov::Csr;
use rssn::sim::structured::LinearSolver;
use rssn::sim::structured::LinearSystem;
#[test]
fn direct_system_rejects_bad_storage_rhs_and_unrepresentable_solution() {
    let a = Csr::from_triplets(2, 2, &[(0, 0, 2.0), (1, 1, 3.0)]);
    let system = LinearSystem::new(a.clone(), LinearSolver::Direct).unwrap();
    assert!(system.solve(&[1.0], None).is_err());
    assert!(system.solve(&[f64::NAN, 1.0], None).is_err());
    assert!(
        system
            .solve(&[1.0, 1.0], Some(&[f64::INFINITY, 0.0]))
            .is_err()
    );
    let mut bad = a.clone();
    bad.indptr = vec![0, 3, 2];
    assert!(LinearSystem::new(bad, LinearSolver::Direct).is_err());
    let mut bad = a.clone();
    bad.values[0] = f64::NAN;
    assert!(LinearSystem::new(bad, LinearSolver::Direct).is_err());
    assert!(LinearSystem::new(a, LinearSolver::Gmres(f64::NAN)).is_err());
    let tiny = Csr::from_triplets(1, 1, &[(0, 0, 1e-308)]);
    let system = LinearSystem::new(tiny, LinearSolver::Direct).unwrap();
    assert!(system.solve(&[1e308], None).is_err());
}
