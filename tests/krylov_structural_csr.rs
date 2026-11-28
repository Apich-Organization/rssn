//! Rectangular structural sparse algebra, measured damping and LSQR certification.

use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::Value;
use rssn::kernel::manifest::SignatureSource;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::RecordValue;
use rssn::kernel::value::TensorValue;
use rssn::kernels::krylov::Csr;
use rssn::kernels::krylov::lsqr;

fn vector(
    data: Vec<Value>,
    elem: TypeExpr,
) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![data.len()],
        data,
        elem,
    })
}

fn measured(
    value: f64,
    dim: &DimExpr,
) -> Value {
    Value::Dyn(DynValue {
        value: Box::new(Value::Float(value)),
        dim: dim.clone(),
    })
}

fn floats(
    values: &[f64],
    dim: Option<&DimExpr>,
) -> Value {
    let data = values
        .iter()
        .map(|&value| dim.map_or(Value::Float(value), |dim| measured(value, dim)))
        .collect();
    let elem = dim.map_or(TypeExpr::Float, |dim| {
        TypeExpr::Dim(Box::new(TypeExpr::Float), dim.clone())
    });
    vector(data, elem)
}

fn csr(
    matrix: &Csr,
    dim: Option<&DimExpr>,
) -> Value {
    let integers = |values: &[usize]| {
        vector(
            values.iter().copied().map(Value::int).collect(),
            TypeExpr::Int,
        )
    };
    Value::Record(RecordValue {
        fields: vec![
            ("rows".into(), Value::int(matrix.rows)),
            ("cols".into(), Value::int(matrix.cols)),
            ("indptr".into(), integers(&matrix.indptr)),
            ("indices".into(), integers(&matrix.indices)),
            ("values".into(), floats(&matrix.values, dim)),
        ],
    })
}

fn field<'a>(
    value: &'a Value,
    name: &str,
) -> &'a Value {
    let Value::Record(record) = value else {
        panic!("expected record")
    };
    record.get(name).expect("record field")
}

fn numeric(value: &Value) -> f64 {
    match value {
        | Value::Dyn(dynamic) => numeric(&dynamic.value),
        | _ => value.as_f64().expect("numeric value"),
    }
}

fn entries(value: &Value) -> Vec<f64> {
    let Value::Tensor(tensor) = value else {
        panic!("expected tensor")
    };
    tensor.data.iter().map(numeric).collect()
}

#[test]
fn rectangular_products_match_dense_reference_and_keep_dimensions() {
    let length = DimExpr::base("Length");
    let time = DimExpr::base("Time");
    let a = Csr::from_triplets(2, 3, &[(0, 0, 1.0), (0, 2, 2.0), (1, 1, 3.0)]);
    let b = Csr::from_triplets(3, 2, &[(0, 0, 4.0), (0, 1, 1.0), (1, 0, 5.0), (2, 1, 6.0)]);
    let dense = Value::Tensor(TensorValue {
        shape: vec![3, 2],
        data: [4.0, 1.0, 5.0, 0.0, 0.0, 6.0]
            .into_iter()
            .map(|x| measured(x, &time))
            .collect(),
        elem: TypeExpr::Dim(Box::new(TypeExpr::Float), time.clone()),
    });
    let session = Session::new();
    let ctx = CallCtx::new();
    let matrix = csr(&a, Some(&length));
    let transposed = session
        .call_value("num_csr_transpose", std::slice::from_ref(&matrix), &ctx)
        .unwrap();
    let Value::Tensor(t) = field(&transposed, "values") else {
        panic!("measured values")
    };
    assert_eq!(
        t.elem,
        TypeExpr::Dim(Box::new(TypeExpr::Float), length.clone())
    );
    assert_eq!(entries(field(&transposed, "values")), vec![1.0, 3.0, 2.0]);

    let sparse = session
        .call_value(
            "num_csr_matmul_sparse",
            &[matrix.clone(), csr(&b, Some(&time))],
            &ctx,
        )
        .unwrap();
    let dense = session
        .call_value("num_csr_matmul_dense", &[matrix, dense], &ctx)
        .unwrap();
    assert_eq!(entries(field(&sparse, "values")), vec![4.0, 13.0, 15.0]);
    assert_eq!(entries(&dense), vec![4.0, 13.0, 15.0, 0.0]);
    let (Value::Tensor(sparse), Value::Tensor(dense)) = (field(&sparse, "values"), dense) else {
        panic!("products")
    };
    assert_eq!(sparse.elem, dense.elem);
    assert_ne!(sparse.elem, TypeExpr::Float);
}

#[test]
fn damped_least_squares_keeps_the_normal_equation_and_reports_residual() {
    let length = DimExpr::base("Length");
    let time = DimExpr::base("Time");
    let matrix = Csr::from_triplets(2, 1, &[(0, 0, 2.0)]);
    let session = Session::new();
    let ctx = CallCtx::new();
    let name = "num_csr_lsqr";
    assert_eq!(
        session.manifest().op(name).unwrap().signature_source,
        SignatureSource::Declared
    );
    let args = [
        csr(&matrix, Some(&length)),
        floats(&[6.0, 8.0], Some(&time)),
        measured(1.0, &length),
        Value::Float(1e-12),
        Value::int(50),
    ];
    let result = session.call_value(name, &args, &ctx).unwrap();
    assert_eq!(*field(&result, "converged"), Value::Bool(true));
    let solution = entries(field(&result, "x"));
    assert!((solution[0] - 2.4).abs() < 1e-14);
    assert!((numeric(field(&result, "residual")) - 8.0_f64.hypot(1.2) / 10.0).abs() < 1e-14);
    let mut invalid = args;
    invalid[2] = measured(1.0, &time);
    assert!(session.call_value(name, &invalid, &ctx).is_err());
}

#[test]
fn nontrivial_small_rhs_and_zero_adjoint_behave_differently() {
    let result = lsqr(
        |x, y| y.copy_from_slice(x),
        |x, y| y.copy_from_slice(x),
        1,
        1,
        &[1e-300],
        0.0,
        1e-12,
        20,
    );
    assert!(result.converged);
    assert!((result.x[0] / 1e-300 - 1.0).abs() < 1e-14);
    let zero = lsqr(
        |_, y| y.fill(0.0),
        |_, y| y.fill(0.0),
        2,
        1,
        &[3.0, 4.0],
        1.0,
        1e-12,
        20,
    );
    assert!(zero.converged);
    assert_eq!(zero.x, vec![0.0]);
    assert_eq!(zero.residual, 1.0);
}

#[test]
fn spgemm_is_bitwise_deterministic_across_thread_counts() {
    let n = 512;
    let entries = (0..n)
        .flat_map(|row| {
            (0..24)
                .map(move |offset| (row, (row + offset * 7) % n, (offset as f64 - 12.5) * 0.0125))
        })
        .collect::<Vec<_>>();
    let a = Csr::from_triplets(n, n, &entries);
    let b = a.transpose();
    let run = |threads| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(|| a.spgemm(&b).unwrap())
    };
    let serial = run(1);
    let parallel = run(4);
    assert_eq!(serial.indptr, parallel.indptr);
    assert_eq!(serial.indices, parallel.indices);
    assert!(
        serial
            .values
            .iter()
            .zip(parallel.values)
            .all(|(a, b)| a.to_bits() == b.to_bits())
    );
    assert_eq!(b.transpose(), a);
}
