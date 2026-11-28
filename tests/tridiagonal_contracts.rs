//! Measured structured linear solve contracts.

use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::overload::matches_type;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::TensorValue;
use rssn::kernel::value::Value;

fn measured_vector(
    entries: &[f64],
    dimension: DimExpr,
) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![entries.len()],
        elem: TypeExpr::Float.with_dim(dimension),
        data: entries.iter().copied().map(Value::Float).collect(),
    })
}

#[test]
fn measured_tridiagonal_solves_preserve_rhs_over_matrix_units() {
    let session = Session::new();
    let ctx = CallCtx::new();
    let mass = DimExpr::base("Mass");
    let length = DimExpr::base("Length");
    let rhs_units = DimExpr::Mul(vec![mass.clone(), length.clone()]);
    let lower = measured_vector(&[-1.0, -1.0], mass.clone());
    let diag = measured_vector(&[2.0, 2.0, 2.0], mass.clone());
    let rhs = measured_vector(&[0.0, 0.0, 4.0], rhs_units);
    let cases = [
        (
            "num_dense_solve_tridiagonal",
            vec![lower.clone(), diag.clone(), lower.clone(), rhs.clone()],
        ),
        (
            "num_dense_solve_spd_tridiagonal",
            vec![diag.clone(), lower.clone(), rhs],
        ),
    ];
    for (name, inputs) in cases {
        let result = session.call_value(name, &inputs, &ctx).unwrap();
        let Value::Tensor(tensor) = result else {
            panic!("expected a tensor from {name}");
        };
        assert_eq!(tensor.shape, [3]);
        assert_eq!(tensor.data.len(), 3);
        assert!(matches_type(&TypeExpr::Float.with_dim(length.clone()), &tensor.elem).accepted());
        for (entry, expected) in tensor.data.iter().zip([1.0, 2.0, 3.0]) {
            match entry {
                | Value::Float(value) => assert!((*value - expected).abs() < 1e-14),
                | Value::Dyn(value) => {
                    let Value::Float(real) = value.value.as_ref() else {
                        panic!("expected a measured Float");
                    };
                    assert!((*real - expected).abs() < 1e-14);
                },
                | _ => panic!("expected a real value from {name}"),
            }
        }
    }
    let other_units = measured_vector(&[-1.0, -1.0], length);
    let rhs = measured_vector(&[0.0, 0.0, 4.0], mass);
    assert!(
        session
            .call_value(
                "num_dense_solve_spd_tridiagonal",
                &[diag, other_units, rhs],
                &ctx
            )
            .is_err()
    );
}
