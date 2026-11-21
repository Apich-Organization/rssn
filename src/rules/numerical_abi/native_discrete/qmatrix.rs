#![allow(clippy::wildcard_imports)] // Shared exact coefficient methods.
use super::*;
use crate::kernels::qlinalg::Q;
use crate::kernels::qlinalg::QMat;
use crate::kernels::qlinalg::{
    self,
};
use num_traits::One;
use num_traits::Signed;
use num_traits::Zero;
pub(super) struct QMatrix {
    pub(super) data: QMat,
    pub(super) cols: usize,
}
impl NativeType for QMatrix {
    const KIND: &'static str = "rssn/native/discrete/qmatrix/v1";
    const MAX_BYTES: usize = 1_048_576;

    fn validate(&self) -> Result<(), OpError> {
        limit(self.data.len(), 32)?;
        limit(self.cols, 32)?;
        if self.data.iter().any(|r| r.len() != self.cols) {
            return Err(error("exact matrix is ragged"));
        }
        let mut bits = 0u64;
        for q in self.data.iter().flatten() {
            bits = bits
                .checked_add(q.numer().bits() + q.denom().bits())
                .ok_or_else(|| error("exact coefficient bit count overflows"))?;
        }
        if bits > 1_048_576 {
            return Err(error("exact coefficient storage budget exhausted"));
        }
        Ok(())
    }
}
impl QMatrix {
    pub(super) const fn new(
        data: QMat,
        cols: usize,
    ) -> Self {
        Self { data, cols }
    }
}
pub(super) fn rational(v: &Value) -> Result<Q, OpError> {
    match v {
        | Value::Float(x) => Q::from_float(*x).ok_or_else(|| error("non-finite exact coefficient")),
        | Value::Dyn(d) => rational(&d.value),
        | Value::Quantity(q) => {
            if q.unit.affine
                || q.unit.offset.is_some()
                || q.unit.scale.pi != 0
                || q.unit.scale.denom.is_zero()
            {
                return Err(error(
                    "exact coefficients require rational linear unit scales",
                ));
            }
            Ok(rational(&q.value)? * q.unit.scale.rational())
        },
        | _ => <Q as Arg>::decode(v, 0),
    }
}
pub(super) fn vector(v: &Value) -> Result<Vec<Q>, OpError> {
    let Value::Tensor(t) = v else {
        return Err(error("exact vector expected"));
    };
    if t.shape.len() != 1 || t.shape[0] != t.data.len() {
        return Err(error("exact vector shape mismatch"));
    }
    limit(t.data.len(), 4096)?;
    t.data.iter().map(rational).collect()
}
pub(super) fn matrix(v: &Value) -> Result<QMatrix, OpError> {
    let Value::Tensor(t) = v else {
        return Err(error("exact matrix expected"));
    };
    if t.shape.len() != 2 || t.shape[0].checked_mul(t.shape[1]) != Some(t.data.len()) {
        return Err(error("exact matrix shape mismatch"));
    }
    limit(t.shape[0], 32)?;
    limit(t.shape[1], 32)?;
    let mut rows = Vec::new();
    for i in 0..t.shape[0] {
        rows.push(
            t.data[i * t.shape[1]..(i + 1) * t.shape[1]]
                .iter()
                .map(rational)
                .collect::<Result<Vec<_>, _>>()?,
        );
    }
    let matrix = QMatrix::new(rows, t.shape[1]);
    matrix.validate()?;
    Ok(matrix)
}
pub(super) fn packed(
    data: QMat,
    cols: usize,
) -> Result<Dynamic, OpError> {
    let model = QMatrix::new(data, cols);
    model.validate()?;
    Ok(Dynamic(Value::Tensor(TensorValue {
        shape: vec![model.data.len(), cols],
        data: model
            .data
            .into_iter()
            .flatten()
            .map(Value::rational)
            .collect(),
        elem: TypeExpr::Real,
    })))
}
fn columns(
    model: &QMatrix,
    cols: usize,
) -> Result<(), OpError> {
    limit(cols, 32)?;
    if !model.data.is_empty() && cols != model.cols {
        return Err(error(
            "explicit exact matrix column count disagrees with row storage",
        ));
    }
    Ok(())
}
fn matrix_growth(model: &QMatrix) -> Result<(), OpError> {
    let largest = model
        .data
        .iter()
        .flatten()
        .map(|q| q.numer().bits() + q.denom().bits())
        .max()
        .unwrap_or(0);
    if largest
        .saturating_mul(model.data.len().max(1) as u64)
        .saturating_mul(model.data.len().saturating_mul(model.cols) as u64)
        > 8_388_608
    {
        return Err(error(
            "exact elimination coefficient-growth work budget exceeded",
        ));
    }
    Ok(())
}
object_binding!(qmatrix_new,"num_native_qmatrix_new","","Construct exact rational coefficient algebra after explicit division by a same-unit positive reference. All subsequent Q operations are normalized coefficient-space operations; physical scales are retained by the caller rather than silently erased. Finite Float inputs are converted to their exact binary rational values.",false,(values:Input,reference:Input)->NewObject<QMatrix>,|_ctx:&CallCtx<'_>,values:Input,reference:Input|{api::same_dimension(&api::dimension(&values.0)?,&api::dimension(&reference.0)?)?;let scale=rational(&reference.0)?;if !scale.is_positive(){return Err(error("exact coefficient reference must be positive"));}let mut model=matrix(&values.0)?;for x in model.data.iter_mut().flatten(){*x/=&scale;}model.validate()?;Ok(NewObject::normalized(model))});
object_binding!(qmatrix_export,"num_native_qmatrix_export","","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<QMatrix>)->Dynamic,|ctx:&CallCtx<'_>,a:Handle<QMatrix>|{a.with(ctx,|m|packed(m.data.clone(),m.cols))});
object_binding!(qlinalg_q,"num_native_qlinalg_q","crate::kernels::qlinalg::q","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(n:i64)->Q,|_ctx:&CallCtx<'_>,n:i64|{Ok(qlinalg::q(n))});
object_binding!(qlinalg_zeros,"num_native_qlinalg_zeros","crate::kernels::qlinalg::zeros","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(n:usize)->Vec<Q>,|_ctx:&CallCtx<'_>,n:usize|{limit(n,4096)?;Ok(qlinalg::zeros(n))});
object_binding!(qlinalg_identity,"num_native_qlinalg_identity","crate::kernels::qlinalg::identity","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(n:usize)->NewObject<QMatrix>,|_ctx:&CallCtx<'_>,n:usize|{limit(n,32)?;Ok(NewObject::normalized(QMatrix::new(qlinalg::identity(n),n)))});
object_binding!(qlinalg_is_zero_vec,"num_native_qlinalg_is_zero_vec","crate::kernels::qlinalg::is_zero_vec","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(data:Input)->bool,|_ctx:&CallCtx<'_>,data:Input|{api::normalized(&data.0)?;Ok(qlinalg::is_zero_vec(&vector(&data.0)?))});
object_binding!(qlinalg_det,"num_native_qlinalg_det","crate::kernels::qlinalg::det","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<QMatrix>)->Q,|ctx:&CallCtx<'_>,a:Handle<QMatrix>|{a.with(ctx,|m|{if m.data.len()!=m.cols{return Err(error("exact determinant requires a square matrix"));}matrix_growth(m)?;Ok(qlinalg::det(&m.data))})});
object_binding!(qlinalg_inverse,"num_native_qlinalg_inverse","crate::kernels::qlinalg::inverse","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<QMatrix>)->Maybe<NewObject<QMatrix>>,|ctx:&CallCtx<'_>,a:Handle<QMatrix>|{a.with(ctx,|m|{if m.data.len()!=m.cols{return Err(error("exact inverse requires square shape"));}matrix_growth(m)?;Ok(Maybe(qlinalg::inverse(&m.data).map(|x|NewObject::normalized(QMatrix::new(x,m.cols)))))})});
object_binding!(qlinalg_matmul,"num_native_qlinalg_matmul","crate::kernels::qlinalg::matmul","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<QMatrix>,b:Handle<QMatrix>)->NewObject<QMatrix>,|ctx:&CallCtx<'_>,a:Handle<QMatrix>,b:Handle<QMatrix>|{a.with(ctx,|a|b.with(ctx,|b|{if a.cols!=b.data.len(){return Err(error("exact matrix product inner dimensions disagree"));}let data=if b.data.is_empty(){vec![qlinalg::zeros(b.cols);a.data.len()]}else{qlinalg::matmul(&a.data,&b.data)};Ok(NewObject::normalized(QMatrix::new(data,b.cols)))}))});
object_binding!(qlinalg_matvec,"num_native_qlinalg_matvec","crate::kernels::qlinalg::matvec","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<QMatrix>,v:Input)->Vec<Q>,|ctx:&CallCtx<'_>,a:Handle<QMatrix>,v:Input|{let v=vector(&v.0)?;a.with(ctx,|a|{if a.cols!=v.len(){return Err(error("exact matrix-vector dimensions disagree"));}Ok(qlinalg::matvec(&a.data,&v))})});
object_binding!(qlinalg_rank,"num_native_qlinalg_rank","crate::kernels::qlinalg::rank","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<QMatrix>,cols:usize)->usize,|ctx:&CallCtx<'_>,a:Handle<QMatrix>,cols:usize|{a.with(ctx,|m|{columns(m,cols)?;matrix_growth(m)?;Ok(qlinalg::rank(&m.data,cols))})});
object_binding!(qlinalg_nullspace,"num_native_qlinalg_nullspace","crate::kernels::qlinalg::nullspace","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<QMatrix>,cols:usize)->NewObject<QMatrix>,|ctx:&CallCtx<'_>,a:Handle<QMatrix>,cols:usize|{a.with(ctx,|m|{columns(m,cols)?;matrix_growth(m)?;Ok(NewObject::normalized(QMatrix::new(qlinalg::nullspace(&m.data,cols),cols)))})});
object_binding!(qlinalg_span_basis,"num_native_qlinalg_span_basis","crate::kernels::qlinalg::span_basis","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<QMatrix>,cols:usize)->NewObject<QMatrix>,|ctx:&CallCtx<'_>,a:Handle<QMatrix>,cols:usize|{a.with(ctx,|m|{columns(m,cols)?;matrix_growth(m)?;Ok(NewObject::normalized(QMatrix::new(qlinalg::span_basis(&m.data,cols),cols)))})});
object_binding!(qlinalg_transpose,"num_native_qlinalg_transpose","crate::kernels::qlinalg::transpose","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<QMatrix>,cols:usize)->NewObject<QMatrix>,|ctx:&CallCtx<'_>,a:Handle<QMatrix>,cols:usize|{a.with(ctx,|m|{columns(m,cols)?;matrix_growth(m)?;Ok(NewObject::normalized(QMatrix::new(qlinalg::transpose(&m.data,cols),m.data.len())))})});
object_binding!(qlinalg_rref,"num_native_qlinalg_rref","crate::kernels::qlinalg::rref","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<QMatrix>,cols:usize)->Vec<usize>,|ctx:&CallCtx<'_>,a:Handle<QMatrix>,cols:usize|{a.with_mut(ctx,|m|{columns(m,cols)?;matrix_growth(m)?;let mut data=m.data.clone();let pivots=qlinalg::rref(&mut data,cols);let candidate=QMatrix::new(data,cols);candidate.validate()?;*m=candidate;Ok(pivots)})});
object_binding!(qlinalg_solve,"num_native_qlinalg_solve","crate::kernels::qlinalg::solve","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<QMatrix>,cols:usize,b:Input)->Option<Vec<Q>>,|ctx:&CallCtx<'_>,a:Handle<QMatrix>,cols:usize,b:Input|{let b=vector(&b.0)?;a.with(ctx,|m|{columns(m,cols)?;matrix_growth(m)?;if b.len()!=m.data.len(){return Err(error("exact solve right-hand side row count disagrees"));}Ok(qlinalg::solve(&m.data,cols,&b))})});
object_binding!(qlinalg_reduce_mod,"num_native_qlinalg_reduce_mod","crate::kernels::qlinalg::reduce_mod","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(v:Input,basis:Handle<QMatrix>,pivots:Vec<usize>)->Vec<Q>,|ctx:&CallCtx<'_>,v:Input,basis:Handle<QMatrix>,pivots:Vec<usize>|{let v=vector(&v.0)?;basis.with(ctx,|m|{if m.cols!=v.len()||pivots.len()!=m.data.len()||pivots.iter().any(|p|*p>=v.len()){return Err(error("row-reduced basis/pivot dimensions disagree"));}for(i,p)in pivots.iter().enumerate(){if m.data[i][*p]!=Q::one()||m.data.iter().enumerate().any(|(j,row)|j!=i&&!row[*p].is_zero()){return Err(error("reduce_mod needs a reduced basis with unit pivot columns"));}}Ok(qlinalg::reduce_mod(&v,&m.data,&pivots))})});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &qmatrix_new::BINDING,
            &qmatrix_export::BINDING,
            &qlinalg_q::BINDING,
            &qlinalg_zeros::BINDING,
            &qlinalg_identity::BINDING,
            &qlinalg_is_zero_vec::BINDING,
            &qlinalg_det::BINDING,
            &qlinalg_inverse::BINDING,
            &qlinalg_matmul::BINDING,
            &qlinalg_matvec::BINDING,
            &qlinalg_rank::BINDING,
            &qlinalg_nullspace::BINDING,
            &qlinalg_span_basis::BINDING,
            &qlinalg_transpose::BINDING,
            &qlinalg_rref::BINDING,
            &qlinalg_solve::BINDING,
            &qlinalg_reduce_mod::BINDING,
        ],
        "numerical_vector",
    )
}
