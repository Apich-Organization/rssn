//! Validation and preconditioner parity for borrowed CSR solvers.

use rssn::kernels::numeric_common::NumericError;
use rssn::kernels::sparse_native::Csr;
use rssn::kernels::sparse_native::Krylov;
use rssn::kernels::sparse_native::PrecondKind;
use rssn::kernels::sparse_native::cg_precond;
use rssn::kernels::sparse_native::krylov_precond;

const METHODS: [Krylov; 2] = [Krylov::Gmres { restart: 3 }, Krylov::Bicgstab];

#[test]
fn preconditioned_solvers_report_input_errors_before_setup() {
    let identity = Csr {
        rows: 1,
        cols: 1,
        row_ptr: &[0, 1],
        col_idx: &[0],
        values: &[1.0],
    };
    for (rhs, tolerance) in [
        (f64::NAN, 1e-8),
        (f64::INFINITY, 1e-8),
        (1.0, f64::NAN),
        (1.0, f64::INFINITY),
        (1.0, -1.0),
    ] {
        assert!(matches!(
            cg_precond(&identity, &[rhs], tolerance, 0, PrecondKind::Identity, None),
            Err(NumericError::Domain(_))
        ));
        for method in METHODS {
            assert!(matches!(
                krylov_precond(
                    &identity,
                    &[rhs],
                    tolerance,
                    0,
                    method,
                    PrecondKind::Jacobi,
                    None
                ),
                Err(NumericError::Domain(_))
            ));
        }
    }
    assert!(matches!(
        cg_precond(&identity, &[], 0.0, 0, PrecondKind::Identity, None),
        Err(NumericError::ShapeMismatch(_))
    ));
    assert!(matches!(
        krylov_precond(
            &Csr { cols: 2, ..identity },
            &[1.0],
            0.0,
            0,
            Krylov::Bicgstab,
            PrecondKind::Jacobi,
            None
        ),
        Err(NumericError::NotSquare)
    ));
    assert!(matches!(
        cg_precond(
            &Csr {
                values: &[f64::NEG_INFINITY],
                ..identity
            },
            &[1.0],
            0.0,
            0,
            PrecondKind::Identity,
            None
        ),
        Err(NumericError::Domain(_))
    ));
}

#[test]
fn jacobi_uses_first_stored_diagonal_while_matvec_sums_duplicates() {
    let a = Csr {
        rows: 2,
        cols: 2,
        row_ptr: &[0, 3, 5],
        col_idx: &[1, 0, 0, 1, 0],
        values: &[0.0, 2.0, 3.0, 4.0, 0.0],
    };
    let b = a.matvec(&[2.0, -1.0]).unwrap();
    assert_eq!(b, [10.0, -4.0]);
    let first = cg_precond(&a, &b, 0.0, 1, PrecondKind::Jacobi, None).unwrap();
    let alpha = 54.0 / 129.0;
    assert_eq!(first.iterations, 1);
    assert_eq!(first.precond, "Jacobi");
    assert!((first.x[0] - 5.0 * alpha).abs() < 1e-14);
    assert!((first.x[1] + alpha).abs() < 1e-14);
    let solved = cg_precond(&a, &b, 1e-12, 10, PrecondKind::Jacobi, None).unwrap();
    assert!(solved.converged);
    assert!((solved.x[0] - 2.0).abs() < 1e-12);
    assert!((solved.x[1] + 1.0).abs() < 1e-12);
}

#[test]
fn borrowed_products_preserve_convergence_and_preconditioner_selection() {
    let a = Csr {
        rows: 2,
        cols: 2,
        row_ptr: &[0, 2, 4],
        col_idx: &[0, 1, 0, 1],
        values: &[4.0, -1.0, -1.0, 3.0],
    };
    let expected = [2.0, -1.0];
    let b = a.matvec(&expected).unwrap();
    for (kind, block, cg_name, general_name) in [
        (PrecondKind::Identity, None, "Identity", "Identity"),
        (PrecondKind::Jacobi, None, "Jacobi", "Jacobi"),
        (
            PrecondKind::BlockJacobi,
            Some((2, 2)),
            "BlockJacobi",
            "BlockJacobi",
        ),
        (PrecondKind::Ilu0, None, "Ilu0", "Ilu0"),
        (PrecondKind::Ilu0, Some((2, 2)), "BlockIlu0", "BlockIlu0"),
        (PrecondKind::Ic0, None, "Ic0", "Ilu0"),
        (PrecondKind::Auto, None, "Ic0", "Ilu0"),
    ] {
        let initial = cg_precond(&a, &b, 0.0, 0, kind, block).unwrap();
        assert_eq!(initial.iterations, 0);
        assert_eq!(initial.residual.to_bits(), 1.0_f64.to_bits());
        assert_eq!(initial.x, [0.0, 0.0]);
        assert!(!initial.converged);
        assert_eq!(initial.precond, cg_name);

        let cg = cg_precond(&a, &b, 1e-12, 12, kind, block).unwrap();
        assert_eq!(cg.precond, cg_name);
        assert!(cg.converged && cg.residual <= 1e-12);
        for (x, expected) in cg.x.iter().zip(expected) {
            assert!((x - expected).abs() < 1e-12);
        }
        for method in METHODS {
            let solved = krylov_precond(&a, &b, 1e-12, 12, method, kind, block).unwrap();
            assert_eq!(solved.precond, general_name);
            assert!(solved.converged && solved.residual <= 1e-12);
            for (x, expected) in solved.x.iter().zip(expected) {
                assert!((x - expected).abs() < 1e-12);
            }
        }
    }
}
