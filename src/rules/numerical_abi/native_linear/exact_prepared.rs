//! Checked prime fields and sparse exact routines over retained matrices.

use super::Installer;
use super::PrimeMatrix;
use super::RuleError;
use super::api;
use super::bounded_q;
use super::object_binding;
use super::prime;
use super::se;
use crate::graph::dim::DimRule;
use crate::kernel::CallCtx;
use crate::kernel::OpError;
use crate::kernel::Value;
use crate::kernel::signature::TypeExpr;
use crate::kernel::value::TensorValue;
use crate::kernels::modular::Barrett;
use crate::rules::numerical_abi::codecs::Arg;
use crate::rules::numerical_abi::codecs::Output;
use api::Handle;
use api::Maybe;
use api::NativeArg;
use api::NativeOutput;
use api::NativeType;
use api::NewObject;
use api::error;
use num_bigint::BigInt;
use num_rational::BigRational;
use std::mem::size_of;

type RationalMatrix = se::SparseMat<BigRational>;

struct PrimeField(se::Gfp);
impl NativeType for PrimeField {
    const KIND: &'static str = "rssn/native/sparse-prime-field/v1";
    const MAX_BYTES: usize = size_of::<Self>();

    fn validate(&self) -> Result<(), OpError> {
        Ok(())
    }
}

object_binding!(
    field_new, "num_linear_sparse_exact_field_new", "crate::kernels::sparse_exact::Gfp::new",
    "Prepares a session-owned prime field after validating the full unsigned 64-bit modulus.",
    true, (modulus: u64) -> NewObject<PrimeField>,
    |_, modulus| Ok(NewObject::normalized(PrimeField(prime(modulus)?)))
);

object_binding!(
    field_reducer, "num_linear_sparse_exact_field_reducer", "crate::kernels::sparse_exact::Gfp::reducer",
    "Returns an independently owned, reusable Barrett reducer for this prime field.",
    true, (field: Handle<PrimeField>) -> NewObject<Barrett>,
    |ctx, field| field.with(ctx, |field| Ok(NewObject::normalized(field.0.reducer())))
);

macro_rules! field_scalar {
    ($id:ident, $op:literal, $source:literal, $doc:literal, $method:ident $(, $arg:ident)* => $result:ty) => {
        object_binding!(
            $id, $op, $source, $doc, true, (field: Handle<PrimeField> $(, $arg: u64)*) -> $result,
            |ctx, field $(, $arg)*| field.with(ctx, |field| Ok(field.0.$method($($arg),*)))
        );
    };
}

field_scalar!(
    field_prime, "num_linear_sparse_exact_field_prime", "crate::kernels::sparse_exact::Gfp::prime",
    "Returns the exact prime modulus.", prime => u64
);
field_scalar!(
    field_add, "num_linear_sparse_exact_field_add", "crate::kernels::sparse_exact::Gfp::add_mod",
    "Returns the canonical residue of a plus b, including for unreduced u64 inputs.", add_mod, left, right => u64
);
field_scalar!(
    field_sub, "num_linear_sparse_exact_field_sub", "crate::kernels::sparse_exact::Gfp::sub_mod",
    "Returns the canonical residue of a minus b, including for unreduced u64 inputs.", sub_mod, left, right => u64
);
field_scalar!(
    field_mul, "num_linear_sparse_exact_field_mul", "",
    "Returns the canonical residue of a times b using the prepared field.", mul_mod, left, right => u64
);
field_scalar!(
    field_pow, "num_linear_sparse_exact_field_pow", "",
    "Raises an element to an unsigned exponent using the prepared field.", pow_mod, base, exponent => u64
);
field_scalar!(
    field_inv, "num_linear_sparse_exact_field_inv", "",
    "Returns the canonical inverse of a field element, or None for a zero residue.", inv, value => Option<u64>
);

object_binding!(
    field_reduce, "num_linear_sparse_exact_field_reduce", "",
    "Reduces a signed exact integer of up to 16,384 bits to its canonical prime-field residue.",
    true, (field: Handle<PrimeField>, value: BigInt) -> u64,
    |ctx, field, value| {
        if value.bits() > 16_384 {
            return Err(OpError::BudgetExceeded);
        }
        field.with(ctx, |field| Ok(field.0.reduce(&value)))
    }
);

object_binding!(
    field_batch_inv, "num_linear_sparse_exact_field_batch_inv", "crate::kernels::sparse_exact::Gfp::batch_inv",
    "Inverts a vector of up to 16,384 residues with one field inversion. Returns None when any entry is zero modulo the prime.",
    true, (field: Handle<PrimeField>, values: Vec<u64>) -> Option<Vec<u64>>,
    |ctx, field, values| {
        if values.len() > 16_384 {
            return Err(OpError::BudgetExceeded);
        }
        field.with(ctx, |field| Ok(field.0.batch_inv(&values)))
    }
);

object_binding!(
    rational_to_prime, "num_linear_sparse_exact_rational_to_prime", "crate::kernels::sparse_exact::reduce_rational_matrix_mod_p",
    "Reduces a rational sparse matrix in the prepared field. A denominator divisible by the prime gives None.",
    true, (matrix: Handle<RationalMatrix>, field: Handle<PrimeField>) -> Maybe<NewObject<PrimeMatrix>>,
    |ctx, matrix, field| {
        let field = field.with(ctx, |field| Ok(field.0))?;
        matrix.with(ctx, |matrix| Ok(Maybe(se::reduce_rational_matrix_mod_p(matrix, field)
            .map(|data| NewObject::normalized(PrimeMatrix { data, field })))))
    }
);

struct PrimeDense(TensorValue);
impl NativeOutput for PrimeDense {
    fn ty() -> TypeExpr {
        api::tensor_type(2, TypeExpr::Int)
    }

    fn encode(
        self,
        _: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        Ok(Value::Tensor(self.0))
    }
}

object_binding!(
    prime_to_dense, "num_linear_sparse_exact_prime_to_dense", "",
    "Copies a prime-field sparse matrix into a rank-two tensor of exact residues.",
    true, (matrix: Handle<PrimeMatrix>) -> PrimeDense,
    |ctx, matrix| matrix.with(ctx, |matrix| {
        let mut field = matrix.field;
        let shape = vec![matrix.data.rows, matrix.data.cols];
        let data = matrix.data.to_dense(&mut field).into_iter().flatten().map(Value::int).collect();
        Ok(PrimeDense(TensorValue { shape, elem: TypeExpr::Int, data }))
    })
);

struct RcmPlan(std::sync::Arc<crate::kernels::reuse::SparseOrderingPlan>);
impl NativeOutput for RcmPlan {
    fn ty() -> TypeExpr {
        let indexes = api::tensor_type(1, TypeExpr::Int);
        TypeExpr::Record(vec![
            ("permutation".into(), indexes.clone()),
            ("inverse".into(), indexes),
        ])
    }

    fn encode(
        self,
        _: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        Ok(api::record(vec![
            ("permutation", Output::encode(self.0.perm.clone())?),
            ("inverse", Output::encode(self.0.inv_perm.clone())?),
        ]))
    }
}

object_binding!(
    rcm_plan, "num_linear_sparse_exact_rcm_plan", "crate::kernels::sparse_exact::SparseMat<E>::rcm_ordering_plan",
    "Returns a cached symmetric Reverse Cuthill–McKee ordering of a square sparse pattern. permutation[new] is the old index, and inverse[old] is the new index.",
    true, (matrix: Handle<RationalMatrix>) -> RcmPlan,
    |ctx, matrix| matrix.with(ctx, |matrix| {
        if matrix.rows != matrix.cols {
            return Err(error("RCM ordering requires a square matrix"));
        }
        Ok(RcmPlan(matrix.rcm_ordering_plan()))
    })
);

struct Rref<T: NativeType> {
    matrix: NewObject<T>,
    pivots: Vec<usize>,
}
impl<T: NativeType> NativeOutput for Rref<T> {
    fn ty() -> TypeExpr {
        TypeExpr::Record(vec![
            ("matrix".into(), <NewObject<T> as NativeOutput>::ty()),
            ("pivots".into(), api::tensor_type(1, TypeExpr::Int)),
        ])
    }

    fn encode(
        self,
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        Ok(api::record(vec![
            ("matrix", self.matrix.encode(ctx)?),
            ("pivots", Output::encode(self.pivots)?),
        ]))
    }
}

fn bound_rationals<'a>(mut values: impl Iterator<Item = &'a BigRational>) -> Result<(), OpError> {
    let mut budget = 8_192_u64;
    values.try_for_each(|value| {
        budget = budget
            .checked_sub(value.numer().bits().saturating_add(value.denom().bits()))
            .ok_or(OpError::BudgetExceeded)?;
        Ok(())
    })
}

object_binding!(
    rational_rref, "num_linear_sparse_exact_rational_rref", "crate::kernels::sparse_exact::rref",
    "Computes exact rational RREF and ascending zero-based pivot columns. Input coefficient bits are bounded before elimination.",
    true, (matrix: Handle<RationalMatrix>) -> Rref<RationalMatrix>,
    |ctx, matrix| matrix.with(ctx, |matrix| {
        bound_rationals(matrix.triplets().map(|(_, _, value)| value))?;
        let (matrix, pivots) = se::rref(&mut se::Rationals, matrix);
        Ok(Rref { matrix: NewObject::normalized(matrix), pivots })
    })
);

object_binding!(
    prime_rref, "num_linear_sparse_exact_prime_rref", "",
    "Computes exact RREF and ascending zero-based pivot columns in the matrix's field.",
    true, (matrix: Handle<PrimeMatrix>) -> Rref<PrimeMatrix>,
    |ctx, matrix| matrix.with(ctx, |matrix| {
        if matrix.data.rows.max(matrix.data.cols) > 128 {
            return Err(OpError::BudgetExceeded);
        }
        let mut field = matrix.field;
        let (data, pivots) = se::rref(&mut field, &matrix.data);
        Ok(Rref { matrix: NewObject::normalized(PrimeMatrix { data, field }), pivots })
    })
);

struct RationalBatch {
    rows: Vec<Vec<BigRational>>,
    columns: usize,
}

impl NativeArg for RationalBatch {
    fn ty() -> TypeExpr {
        api::tensor_type(2, <BigRational as Arg>::ty())
    }

    fn decode(
        value: &Value,
        index: usize,
        _: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        let Value::Tensor(tensor) = value else {
            return Err(error(
                "Dixon requires a rank-two tensor of right-hand sides",
            ));
        };
        let [rows, columns] = tensor.shape.as_slice() else {
            return Err(error(
                "Dixon requires a rank-two tensor of right-hand sides",
            ));
        };
        if *rows > 16 || *columns > 32 || rows.saturating_mul(*columns) != tensor.data.len() {
            return Err(OpError::BudgetExceeded);
        }
        let coefficients = tensor
            .data
            .iter()
            .map(|v| <BigRational as Arg>::decode(v, index))
            .collect::<Result<Vec<_>, _>>()?;
        bound_rationals(coefficients.iter())?;
        let mut coefficients = coefficients.into_iter();
        let rows = (0..*rows)
            .map(|_| coefficients.by_ref().take(*columns).collect())
            .collect();
        Ok(Self {
            rows,
            columns: *columns,
        })
    }
}

impl NativeOutput for RationalBatch {
    fn ty() -> TypeExpr {
        <Self as NativeArg>::ty()
    }

    fn encode(
        self,
        _: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        let shape = vec![self.rows.len(), self.columns];
        let data = self
            .rows
            .into_iter()
            .flatten()
            .map(|value| {
                bounded_q(&value)?;
                Output::encode(value)
            })
            .collect::<Result<_, _>>()?;
        Ok(Value::Tensor(TensorValue {
            shape,
            elem: <BigRational as Arg>::ty(),
            data,
        }))
    }
}

object_binding!(
    dixon_multi, "num_linear_sparse_exact_dixon_multi", "crate::kernels::sparse_exact::solve_dixon_multi",
    "Solves a square rational sparse system with one right-hand side per tensor row. Shares a modular factorization and returns None for a singular system.",
    true, (matrix: Handle<RationalMatrix>, rhs: RationalBatch) -> Maybe<RationalBatch>,
    |ctx, matrix, rhs| matrix.with(ctx, |matrix| {
        let dimension = matrix.rows;
        if matrix.cols != dimension || rhs.columns != dimension {
            return Err(error("Dixon requires a square matrix and right-hand rows with matching width"));
        }
        bound_rationals(matrix.triplets().map(|(_, _, value)| value))?;
        let result = se::solve_dixon_multi(matrix, &rhs.rows).map(|rows| RationalBatch { rows, columns: rhs.columns });
        Ok(Maybe(result))
    })
);

#[allow(
    clippy::single_call_fn,
    reason = "Keeps the exact sparse registrations with their adapters"
)]
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    const NAMESPACE: &str = "numerical_matrix";
    api::install_with_dimension(
        i,
        &[
            &field_prime::BINDING,
            &field_add::BINDING,
            &field_sub::BINDING,
            &field_mul::BINDING,
            &field_pow::BINDING,
            &field_reduce::BINDING,
        ],
        NAMESPACE,
        DimRule::Dimensionless,
    )?;
    api::install(
        i,
        &[
            &field_new::BINDING,
            &field_reducer::BINDING,
            &field_inv::BINDING,
            &field_batch_inv::BINDING,
            &rational_to_prime::BINDING,
            &prime_to_dense::BINDING,
            &rcm_plan::BINDING,
            &rational_rref::BINDING,
            &prime_rref::BINDING,
            &dixon_multi::BINDING,
        ],
        NAMESPACE,
    )
}
