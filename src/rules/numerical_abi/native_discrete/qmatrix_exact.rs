//! Additional bounded exact matrix methods using normalized `QMatrix` handles.
#![allow(clippy::wildcard_imports)] // Shared exact matrix input and output contracts.
use super::*;
use num_bigint::BigInt;

struct IntegerMatrix {
    rows: Vec<Vec<BigInt>>,
    cols: usize,
}

impl IntegerMatrix {
    const fn new(
        rows: Vec<Vec<BigInt>>,
        cols: usize,
    ) -> Self {
        Self { rows, cols }
    }

    fn identity(size: usize) -> Self {
        let rows = (0..size)
            .map(|i| (0..size).map(|j| BigInt::from(u8::from(i == j))).collect())
            .collect();
        Self::new(rows, size)
    }
}

impl Output for IntegerMatrix {
    fn ty() -> TypeExpr {
        api::tensor_type(2, TypeExpr::Int)
    }

    fn encode(self) -> Result<Value, OpError> {
        limit(self.rows.len(), MATRIX_LIMIT)?;
        limit(self.cols, MATRIX_LIMIT)?;
        if self.rows.iter().any(|row| row.len() != self.cols) {
            return Err(error("exact integer matrix result is ragged"));
        }
        coefficient_storage(self.rows.iter().flatten().map(BigInt::bits))?;
        Ok(Value::Tensor(TensorValue::of(
            vec![self.rows.len(), self.cols],
            self.rows.into_iter().flatten().map(Value::int).collect(),
            TypeExpr::Int,
        )))
    }
}

struct HermiteForm {
    form: IntegerMatrix,
    transform: IntegerMatrix,
}

impl Output for HermiteForm {
    fn ty() -> TypeExpr {
        let matrix = <IntegerMatrix as Output>::ty();
        TypeExpr::Record(vec![
            ("form".into(), matrix.clone()),
            ("transform".into(), matrix),
        ])
    }

    fn encode(self) -> Result<Value, OpError> {
        Ok(api::record(vec![
            ("form", Output::encode(self.form)?),
            ("transform", Output::encode(self.transform)?),
        ]))
    }
}

struct SmithForm {
    left: IntegerMatrix,
    form: IntegerMatrix,
    right: IntegerMatrix,
}

impl Output for SmithForm {
    fn ty() -> TypeExpr {
        let matrix = <IntegerMatrix as Output>::ty();
        TypeExpr::Record(vec![
            ("left".into(), matrix.clone()),
            ("form".into(), matrix.clone()),
            ("right".into(), matrix),
        ])
    }

    fn encode(self) -> Result<Value, OpError> {
        Ok(api::record(vec![
            ("left", Output::encode(self.left)?),
            ("form", Output::encode(self.form)?),
            ("right", Output::encode(self.right)?),
        ]))
    }
}

struct RightHandSide(Vec<Q>);

impl NativeArg for RightHandSide {
    fn ty() -> TypeExpr {
        api::tensor_type(1, TypeExpr::Real)
    }

    fn decode(
        value: &Value,
        index: usize,
        _: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        let Value::Tensor(tensor) = value else {
            return Err(OpError::Type {
                expected: Box::new(Self::ty()),
                got: Box::new(value.type_expr()),
                arg: index,
            });
        };
        limit(tensor.data.len(), MATRIX_LIMIT)?;
        vector(value).map(Self)
    }
}

fn square(model: &QMatrix) -> Result<(), OpError> {
    if model.data.len() != model.cols {
        return Err(error("exact matrix operation requires square shape"));
    }
    Ok(())
}

fn integer_coefficients(model: &QMatrix) -> Result<Vec<Vec<BigInt>>, OpError> {
    matrix_growth(model)?;
    model
        .data
        .iter()
        .map(|row| {
            row.iter()
                .map(|entry| {
                    if !entry.is_integer() {
                        return Err(error("exact integer matrix requires integral coefficients"));
                    }
                    Ok(entry.to_integer())
                })
                .collect()
        })
        .collect()
}

fn checked_rationals(values: Option<Vec<Q>>) -> Result<Option<Vec<Q>>, OpError> {
    coefficient_storage(values.iter().flatten().map(rational_bits))?;
    Ok(values)
}

fn rhs_growth(
    a: &QMatrix,
    rhs_cols: usize,
    rhs: impl Iterator<Item = u64> + Clone,
) -> Result<(), OpError> {
    let lhs = a.data.iter().flatten().map(rational_bits);
    coefficient_storage(lhs.clone().chain(rhs.clone()))?;
    elimination_growth(
        a.data.len(),
        a.cols.saturating_add(rhs_cols),
        lhs.chain(rhs),
    )
}

object_binding!(
    qlinalg_charpoly, "num_native_qlinalg_charpoly", "crate::kernels::qlinalg::charpoly",
    "Return the exact characteristic polynomial coefficients in ascending degree, including the leading one. The normalized QMatrix must be square.",
    true, (a: Handle<QMatrix>) -> Option<Vec<Q>>, |ctx, a| a.with(ctx, |a| {
        square(a)?;
        matrix_growth(a)?;
        checked_rationals(qlinalg::charpoly(&a.data))
    })
);
object_binding!(
    qlinalg_det_bareiss, "num_native_qlinalg_det_bareiss", "crate::kernels::qlinalg::det_bareiss",
    "Return the exact determinant of a square normalized QMatrix with integer coefficients, using fraction-free Bareiss elimination.",
    true, (a: Handle<QMatrix>) -> BigInt, |ctx, a| a.with(ctx, |a| {
        square(a)?;
        let result = qlinalg::det_bareiss(&integer_coefficients(a)?);
        coefficient_storage(std::iter::once(result.bits()))?;
        Ok(result)
    })
);
object_binding!(
    qlinalg_det_modular_int, "num_native_qlinalg_det_modular_int", "crate::kernels::qlinalg::det_modular_int",
    "Return the exact determinant of a square normalized QMatrix with integer coefficients, using modular elimination and Chinese remaindering.",
    true, (a: Handle<QMatrix>) -> BigInt, |ctx, a| a.with(ctx, |a| {
        square(a)?;
        let result = qlinalg::det_modular_int(&integer_coefficients(a)?);
        coefficient_storage(std::iter::once(result.bits()))?;
        Ok(result)
    })
);
object_binding!(
    qlinalg_hermite_normal_form, "num_native_qlinalg_hermite_normal_form", "crate::kernels::qlinalg::hermite_normal_form",
    "Return the row-style Hermite normal form of a normalized QMatrix with integer coefficients.",
    true, (a: Handle<QMatrix>) -> IntegerMatrix, |ctx, a| a.with(ctx, |a| {
        let result = qlinalg::hermite_normal_form(&integer_coefficients(a)?);
        Ok(IntegerMatrix::new(result, a.cols))
    })
);
object_binding!(
    qlinalg_hermite_normal_form_with_transform, "num_native_qlinalg_hermite_normal_form_with_transform", "crate::kernels::qlinalg::hermite_normal_form_with_transform",
    "Return the row-style Hermite normal form and unimodular row transform satisfying transform * A = form.",
    true, (a: Handle<QMatrix>) -> HermiteForm, |ctx, a| a.with(ctx, |a| {
        let (form, transform) = qlinalg::hermite_normal_form_with_transform(&integer_coefficients(a)?);
        Ok(HermiteForm {
            form: IntegerMatrix::new(form, a.cols),
            transform: IntegerMatrix::new(transform, a.data.len()),
        })
    })
);
object_binding!(
    qlinalg_smith_normal_form, "num_native_qlinalg_smith_normal_form", "crate::kernels::qlinalg::smith_normal_form",
    "Return the Smith normal form and unimodular transforms satisfying left * A * right = form.",
    true, (a: Handle<QMatrix>) -> SmithForm, |ctx, a| a.with(ctx, |a| {
        let (left, form, right) = qlinalg::smith_normal_form(&integer_coefficients(a)?);
        let right = if a.data.is_empty() {
            IntegerMatrix::identity(a.cols)
        } else {
            IntegerMatrix::new(right, a.cols)
        };
        Ok(SmithForm {
            left: IntegerMatrix::new(left, a.data.len()),
            form: IntegerMatrix::new(form, a.cols),
            right,
        })
    })
);
object_binding!(
    qlinalg_smith_normal_form_diag, "num_native_qlinalg_smith_normal_form_diag", "crate::kernels::qlinalg::smith_normal_form_diag",
    "Return the positive nonzero Smith invariant factors in divisibility order. The normalized QMatrix must have integer coefficients.",
    true, (a: Handle<QMatrix>) -> Vec<BigInt>, |ctx, a| a.with(ctx, |a| {
        let values = qlinalg::smith_normal_form_diag(&integer_coefficients(a)?);
        coefficient_storage(values.iter().map(BigInt::bits))?;
        Ok(values)
    })
);
object_binding!(
    qlinalg_solve_dixon, "num_native_qlinalg_solve_dixon", "crate::kernels::qlinalg::solve_dixon",
    "Solve an exact square system by Dixon lifting and residual certification. A singular system returns None.",
    true, (a: Handle<QMatrix>, b: RightHandSide) -> Option<Vec<Q>>, |ctx, a, b| {
        let b = b.0;
        a.with(ctx, |a| {
            square(a)?;
            if a.data.len() != b.len() {
                return Err(error("exact solve right-hand side row count disagrees"));
            }
            rhs_growth(a, 1, b.iter().map(rational_bits))?;
            checked_rationals(qlinalg::solve_dixon(&a.data, &b))
        })
    }
);
object_binding!(
    qlinalg_solve_mat, "num_native_qlinalg_solve_mat", "crate::kernels::qlinalg::solve_mat",
    "Solve an exact system with multiple right-hand sides. Free variables are zero, and inconsistent systems return None.",
    true, (a: Handle<QMatrix>, b: Handle<QMatrix>) -> Maybe<NewObject<QMatrix>>, |ctx, a, b| {
        a.with(ctx, |a| b.with(ctx, |b| {
            if a.data.len() != b.data.len() {
                return Err(error("exact solve right-hand side row count disagrees"));
            }
            rhs_growth(a, b.cols, b.data.iter().flatten().map(rational_bits))?;
            let result = if a.data.is_empty() {
                Some(vec![qlinalg::zeros(b.cols); a.cols])
            } else {
                qlinalg::solve_mat(&a.data, a.cols, &b.data)
            };
            Ok(Maybe(result.map(|data| NewObject::normalized(QMatrix::new(data, b.cols)))))
        }))
    }
);
object_binding!(
    qlinalg_solve_unique, "num_native_qlinalg_solve_unique", "crate::kernels::qlinalg::solve_unique",
    "Solve an exact system with full column rank. Inconsistent and underdetermined systems return None.",
    true, (a: Handle<QMatrix>, b: RightHandSide) -> Option<Vec<Q>>, |ctx, a, b| {
        let b = b.0;
        a.with(ctx, |a| {
            if a.data.len() != b.len() {
                return Err(error("exact solve right-hand side row count disagrees"));
            }
            rhs_growth(a, 1, b.iter().map(rational_bits))?;
            checked_rationals(qlinalg::solve_unique(&a.data, a.cols, &b))
        })
    }
);

pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &qlinalg_charpoly::BINDING,
            &qlinalg_det_bareiss::BINDING,
            &qlinalg_det_modular_int::BINDING,
            &qlinalg_hermite_normal_form::BINDING,
            &qlinalg_hermite_normal_form_with_transform::BINDING,
            &qlinalg_smith_normal_form::BINDING,
            &qlinalg_smith_normal_form_diag::BINDING,
            &qlinalg_solve_dixon::BINDING,
            &qlinalg_solve_mat::BINDING,
            &qlinalg_solve_unique::BINDING,
        ],
        "numerical_vector",
    )
}
