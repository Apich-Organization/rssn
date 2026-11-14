#![allow(clippy::wildcard_imports)] // Exact algebra in explicitly normalized bases and lattice coordinates.
use super::qmatrix::QMatrix;
use super::qmatrix::{
    self,
};
use super::*;
use crate::kernels::lie_structure::Cartan;
use crate::kernels::lie_structure::Family;
use crate::kernels::lie_structure::Lie;
use crate::kernels::lie_structure::RootSystem;
use crate::kernels::lie_structure::{
    self as lie,
};
use crate::kernels::qlinalg::Q;
use crate::kernels::qlinalg::QMat;
use num_bigint::BigInt;
impl NativeType for RootSystem {
    const KIND: &'static str = "rssn/native/discrete/root-system/v1";
    const MAX_BYTES: usize = 1048576;

    fn validate(&self) -> Result<(), OpError> {
        limit(self.n, 8)?;
        limit(self.positive.len(), 256)?;
        Ok::<_, OpError>(())
    }
}
impl NativeType for Lie {
    const KIND: &'static str = "rssn/native/discrete/lie-algebra/v1";
    const MAX_BYTES: usize = 4194304;

    fn validate(&self) -> Result<(), OpError> {
        limit(self.n, 12)?;
        let mut bits = 0u64;
        for a in &self.c {
            if a.len() != self.n {
                return Err(error("invalid Lie constants shape"));
            }
            for b in a {
                if b.len() != self.n {
                    return Err(error("invalid Lie constants shape"));
                }
                for q in b {
                    bits = bits
                        .checked_add(q.numer().bits() + q.denom().bits())
                        .ok_or_else(|| error("Lie coefficient bit budget overflow"))?;
                }
            }
        }
        if bits > 1_048_576 {
            return Err(error("Lie coefficient bit budget exhausted"));
        }
        Ok::<_, OpError>(())
    }
}
struct Matrices(Vec<QMat>);
impl NativeArg for Matrices {
    fn ty() -> TypeExpr {
        api::tensor_type(1, TypeExpr::Opaque(QMatrix::KIND.into()))
    }

    fn decode(
        v: &Value,
        _: usize,
        c: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        let Value::Tensor(t) = v else {
            return Err(error("expected vector of exact matrix handles"));
        };
        if t.shape != vec![t.data.len()] || t.data.is_empty() || t.data.len() > 12 {
            return Err(error("exact matrix family size must be1..12"));
        }
        let values = t
            .data
            .iter()
            .map(|v| {
                let h = <Handle<QMatrix> as NativeArg>::decode(v, 0, c)?;
                h.with(c, |m| {
                    if m.data.len() != m.cols || m.cols == 0 || m.cols > 12 {
                        return Err(error("Lie matrix family requires square1..12matrices"));
                    }
                    Ok::<_, OpError>(m.data.clone())
                })
            })
            .collect::<Result<Vec<_>, _>>()?;
        let n = values[0].len();
        if values.iter().any(|m| m.len() != n) {
            return Err(error("Lie matrices must have common size"));
        }
        Ok::<_, OpError>(Self(values))
    }
}
fn letter(v: &str) -> Result<char, OpError> {
    let mut chars = v.chars();
    let c = chars
        .next()
        .ok_or_else(|| error("Cartan family letter missing"))?;
    if chars.next().is_some() {
        return Err(error("Cartan family must be one letter"));
    }
    Ok::<_, OpError>(c)
}
fn cartan(v: &Input) -> Result<Cartan, OpError> {
    api::normalized(&v.0)?;
    let Value::Tensor(t) = &v.0 else {
        return Err(error("Cartan integer matrix expected"));
    };
    if t.shape.len() != 2 || t.shape[0] != t.shape[1] || t.shape[0] == 0 || t.shape[0] > 8 {
        return Err(error("Cartan rank must be1..8"));
    }
    let n = t.shape[0];
    let data = t
        .data
        .iter()
        .map(|v| <i64 as Arg>::decode(v, 0))
        .collect::<Result<Vec<_>, _>>()?;
    if data.iter().any(|v| v.abs() > 8) {
        return Err(error("Cartan entry bound exceeded"));
    }
    Ok::<_, OpError>(data.chunks(n).map(<[i64]>::to_vec).collect())
}
fn ints(
    m: Vec<Vec<i64>>,
    cols: usize,
) -> Result<Dynamic, OpError> {
    if m.iter().any(|r| r.len() != cols) {
        return Err(error("integer matrix output ragged"));
    }
    Ok::<_, OpError>(Dynamic(Value::Tensor(TensorValue {
        shape: vec![m.len(), cols],
        data: m.into_iter().flatten().map(Value::int).collect(),
        elem: TypeExpr::Int,
    })))
}
fn qvec(
    v: &Input,
    n: usize,
) -> Result<Vec<Q>, OpError> {
    api::normalized(&v.0)?;
    let values = qmatrix::vector(&v.0)?;
    if values.len() != n {
        return Err(error("Lie coordinate vector has wrong dimension"));
    }
    Ok::<_, OpError>(values)
}
fn qspace(
    v: &Input,
    n: usize,
) -> Result<QMat, OpError> {
    api::normalized(&v.0)?;
    let m = qmatrix::matrix(&v.0)?;
    if m.cols != n || m.data.len() > n {
        return Err(error("Lie subspace must have <=n rows and n columns"));
    }
    Ok::<_, OpError>(m.data)
}
fn labels(
    v: &[i64],
    n: usize,
) -> Result<(), OpError> {
    if v.len() != n || v.iter().any(|x| *x < 0 || *x > 32) {
        return Err(error("Dynkin labels must have rank entries in0..32"));
    }
    Ok::<_, OpError>(())
}
fn weight(
    v: &[i64],
    n: usize,
) -> Result<(), OpError> {
    if v.len() != n || v.iter().any(|x| x.abs() > 256) {
        return Err(error(
            "weight coordinates exceed rank/bounded integer contract",
        ));
    }
    Ok::<_, OpError>(())
}
fn table(
    values: impl IntoIterator<Item = (Vec<i64>, BigInt)>,
    n: usize,
) -> Result<Dynamic, OpError> {
    let (w, m): (Vec<_>, Vec<_>) = values.into_iter().unzip();
    Ok::<_, OpError>(Dynamic(api::record(vec![
        ("weights", ints(w, n)?.0),
        ("multiplicities", Output::encode(m)?),
    ])))
}
fn weight_work(
    s: &RootSystem,
    top: &[i64],
) -> Result<(), OpError> {
    labels(top, s.n)?;
    let dim = s
        .weyl_dim(top)
        .ok_or_else(|| error("invalid highest weight"))?;
    if dim > BigInt::from(10000) {
        return Err(error(
            "representation dimension exceeds10000 weight-work bound",
        ));
    }
    Ok::<_, OpError>(())
}
object_binding!(lie_cartan_matrix,"num_native_lie_cartan_matrix","crate::kernels::lie_structure::cartan_matrix","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(family:String,rank:usize)->Maybe<Dynamic>,|_ctx:&CallCtx<'_>,family:String,rank:usize|{limit(rank,8)?;Ok::<_,OpError>(Maybe(lie::cartan_matrix(letter(&family)?,rank).map(|m|ints(m,rank)).transpose()?))});
object_binding!(lie_parse_type,"num_native_lie_parse_type","crate::kernels::lie_structure::parse_type","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(name:String)->Maybe<Dynamic>,|_ctx:&CallCtx<'_>,name:String|{if name.len()>16{return Err(error("Lie type name length bound exceeded"));}let m=lie::parse_type(&name);if m.as_ref().is_some_and(|m|m.len()>8){return Err(error("Cartan rank budget exceeds8"));}Ok::<_,OpError>(Maybe(m.map(|m|{let n=m.len();ints(m,n)}).transpose()?))});
object_binding!(lie_components,"num_native_lie_components","crate::kernels::lie_structure::components","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(matrix:Input)->Dynamic,|_ctx:&CallCtx<'_>,matrix:Input|{let c=cartan(&matrix)?;let mut offsets=vec![0usize];let mut indices=vec![];for row in lie::components(&c){indices.extend(row);offsets.push(indices.len());}Ok::<_,OpError>(Dynamic(api::record(vec![("offsets",Output::encode(offsets)?),("indices",Output::encode(indices)?)])))});
object_binding!(lie_positive_roots,"num_native_lie_positive_roots","crate::kernels::lie_structure::positive_roots","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(matrix:Input)->Dynamic,|_ctx:&CallCtx<'_>,matrix:Input|{let c=cartan(&matrix)?;if lie::recognize(&c).is_none(){return Err(error("positive roots require finite Cartan type"));}ints(lie::positive_roots(&c),c.len())});
object_binding!(lie_recognize,"num_native_lie_recognize","crate::kernels::lie_structure::recognize","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(matrix:Input)->Maybe<Many<Dynamic>>,|_ctx:&CallCtx<'_>,matrix:Input|{let c=cartan(&matrix)?;Ok::<_,OpError>(Maybe(lie::recognize(&c).map(|v|Many(v.into_iter().map(|(letter,rank)|Dynamic(api::record(vec![("family",Value::String(letter.to_string())),("rank",Value::int(rank))]))).collect()))))});
object_binding!(lie_euclidean_simple_roots,"num_native_lie_euclidean_simple_roots","crate::kernels::lie_structure::euclidean_simple_roots","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(family:String,rank:usize)->Maybe<Dynamic>,|_ctx:&CallCtx<'_>,family:String,rank:usize|{limit(rank,8)?;Ok::<_,OpError>(Maybe(lie::euclidean_simple_roots(letter(&family)?,rank).map(|m|{let n=m.first().map_or(0,Vec::len);qmatrix::packed(m,n)}).transpose()?))});
object_binding!(lie_standard_basis,"num_native_lie_standard_basis","crate::kernels::lie_structure::standard_basis","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(family:String,n:usize)->Maybe<Many<Dynamic>>,|_ctx:&CallCtx<'_>,family:String,n:usize|{if n==0||n>6{return Err(error("classical matrix basis order must be1..6"));}let f=match family.as_str(){"Gl"=>Family::Gl,"Sl"=>Family::Sl,"So"=>Family::So,"Sp"=>Family::Sp,"Su"=>Family::Su,_=>return Err(error("family must be Gl/Sl/So/Sp/Su"))};Ok::<_,OpError>(Maybe(lie::standard_basis(f,n).map(|basis|{Ok::<_,OpError>(Many(basis.into_iter().map(|m|{let n=m.len();let re=m.iter().map(|r|r.iter().map(|z|z.0).collect()).collect();let im=m.into_iter().map(|r|r.into_iter().map(|z|z.1).collect()).collect();Ok::<_,OpError>(Dynamic(api::record(vec![("real",ints(re,n)?.0),("imaginary",ints(im,n)?.0)])))}).collect::<Result<_,OpError>>()?))}).transpose()?))});
object_binding!(lie_charpoly,"num_native_lie_charpoly","crate::kernels::lie_structure::charpoly","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(matrix:Handle<QMatrix>)->Vec<Q>,|ctx:&CallCtx<'_>,matrix:Handle<QMatrix>|{matrix.with(ctx,|m|{if m.cols!=m.data.len()||m.cols>12{return Err(error("characteristic polynomial requires square<=12matrix"));}Ok::<_,OpError>(lie::charpoly(&m.data))})});
object_binding!(lie_rational_roots,"num_native_lie_rational_roots","crate::kernels::lie_structure::rational_roots","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(coefficients:Input)->Option<Vec<Q>>,|_ctx:&CallCtx<'_>,coefficients:Input|{api::normalized(&coefficients.0)?;let p=qmatrix::vector(&coefficients.0)?;limit(p.len(),33)?;if p.iter().any(|q|q.numer().bits()>24||q.denom().bits()>16){return Err(error("rational root factorization coefficient bit budget exceeded"));}Ok::<_,OpError>(lie::rational_roots(&p))});
object_binding!(root_system_new,"num_native_root_system_new","crate::kernels::lie_structure::RootSystem::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(matrix:Input)->Maybe<NewObject<RootSystem>>,|_ctx:&CallCtx<'_>,matrix:Input|{let c=cartan(&matrix)?;Ok::<_,OpError>(Maybe(RootSystem::new(c).map(NewObject::normalized)))});
object_binding!(root_system_num_positive,"num_native_root_system_num_positive","crate::kernels::lie_structure::RootSystem::num_positive","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<RootSystem>)->usize,|ctx:&CallCtx<'_>,a:Handle<RootSystem>|{a.with(ctx,|s|Ok::<_,OpError>(s.num_positive()))});
object_binding!(root_system_dimension,"num_native_root_system_dimension","crate::kernels::lie_structure::RootSystem::dimension","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<RootSystem>)->usize,|ctx:&CallCtx<'_>,a:Handle<RootSystem>|{a.with(ctx,|s|Ok::<_,OpError>(s.dimension()))});
object_binding!(root_system_is_irreducible,"num_native_root_system_is_irreducible","crate::kernels::lie_structure::RootSystem::is_irreducible","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<RootSystem>)->bool,|ctx:&CallCtx<'_>,a:Handle<RootSystem>|{a.with(ctx,|s|Ok::<_,OpError>(s.is_irreducible()))});
object_binding!(root_system_exponents,"num_native_root_system_exponents","crate::kernels::lie_structure::RootSystem::exponents","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<RootSystem>)->Vec<usize>,|ctx:&CallCtx<'_>,a:Handle<RootSystem>|{a.with(ctx,|s|Ok::<_,OpError>(s.exponents()))});
object_binding!(root_system_weyl_order,"num_native_root_system_weyl_order","crate::kernels::lie_structure::RootSystem::weyl_order","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<RootSystem>)->BigInt,|ctx:&CallCtx<'_>,a:Handle<RootSystem>|{a.with(ctx,|s|Ok::<_,OpError>(s.weyl_order()))});
object_binding!(root_system_highest_root,"num_native_root_system_highest_root","crate::kernels::lie_structure::RootSystem::highest_root","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<RootSystem>)->Option<Vec<i64>>,|ctx:&CallCtx<'_>,a:Handle<RootSystem>|{a.with(ctx,|s|Ok::<_,OpError>(s.highest_root()))});
object_binding!(root_system_coxeter_number,"num_native_root_system_coxeter_number","crate::kernels::lie_structure::RootSystem::coxeter_number","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<RootSystem>)->Option<usize>,|ctx:&CallCtx<'_>,a:Handle<RootSystem>|{a.with(ctx,|s|Ok::<_,OpError>(s.coxeter_number()))});
object_binding!(root_system_dual_coxeter_number,"num_native_root_system_dual_coxeter_number","crate::kernels::lie_structure::RootSystem::dual_coxeter_number","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<RootSystem>)->Option<Q>,|ctx:&CallCtx<'_>,a:Handle<RootSystem>|{a.with(ctx,|s|Ok::<_,OpError>(s.dual_coxeter_number()))});
object_binding!(root_system_gram,"num_native_root_system_gram","crate::kernels::lie_structure::RootSystem::gram","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<RootSystem>)->NewObject<QMatrix>,|ctx:&CallCtx<'_>,a:Handle<RootSystem>|{a.with(ctx,|s|Ok::<_,OpError>(NewObject::normalized(QMatrix::new(s.gram(),s.n))))});
object_binding!(root_system_fundamental_weights,"num_native_root_system_fundamental_weights","crate::kernels::lie_structure::RootSystem::fundamental_weights","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<RootSystem>)->NewObject<QMatrix>,|ctx:&CallCtx<'_>,a:Handle<RootSystem>|{a.with(ctx,|s|Ok::<_,OpError>(NewObject::normalized(QMatrix::new(s.fundamental_weights().clone(),s.n))))});
object_binding!(root_system_weyl_dim,"num_native_root_system_weyl_dim","crate::kernels::lie_structure::RootSystem::weyl_dim","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<RootSystem>,highest:Vec<i64>)->Option<BigInt>,|ctx:&CallCtx<'_>,a:Handle<RootSystem>,highest:Vec<i64>|{a.with(ctx,|s|{labels(&highest,s.n)?;Ok::<_,OpError>(s.weyl_dim(&highest))})});
object_binding!(root_system_casimir,"num_native_root_system_casimir","crate::kernels::lie_structure::RootSystem::casimir","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<RootSystem>,highest:Vec<i64>)->Option<Q>,|ctx:&CallCtx<'_>,a:Handle<RootSystem>,highest:Vec<i64>|{a.with(ctx,|s|{labels(&highest,s.n)?;Ok::<_,OpError>(s.casimir(&highest))})});
object_binding!(root_system_dynkin_index,"num_native_root_system_dynkin_index","crate::kernels::lie_structure::RootSystem::dynkin_index","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<RootSystem>,highest:Vec<i64>)->Option<Q>,|ctx:&CallCtx<'_>,a:Handle<RootSystem>,highest:Vec<i64>|{a.with(ctx,|s|{labels(&highest,s.n)?;Ok::<_,OpError>(s.dynkin_index(&highest))})});
object_binding!(root_system_inner,"num_native_root_system_inner","crate::kernels::lie_structure::RootSystem::inner","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<RootSystem>,x:Vec<i64>,y:Vec<i64>)->Q,|ctx:&CallCtx<'_>,a:Handle<RootSystem>,x:Vec<i64>,y:Vec<i64>|{a.with(ctx,|s|{weight(&x,s.n)?;weight(&y,s.n)?;Ok::<_,OpError>(s.inner(&x,&y))})});
object_binding!(root_system_dominant_multiplicities,"num_native_root_system_dominant_multiplicities","crate::kernels::lie_structure::RootSystem::dominant_multiplicities","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<RootSystem>,top:Vec<i64>,cap:usize)->Maybe<Dynamic>,|ctx:&CallCtx<'_>,a:Handle<RootSystem>,top:Vec<i64>,cap:usize|{limit(cap,10000)?;a.with(ctx,|s|{weight_work(s,&top)?;Ok::<_,OpError>(Maybe(s.dominant_multiplicities(&top,cap).map(|v|table(v,s.n)).transpose()?))})});
object_binding!(root_system_all_weights,"num_native_root_system_all_weights","crate::kernels::lie_structure::RootSystem::all_weights","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<RootSystem>,top:Vec<i64>,cap:usize)->Maybe<Dynamic>,|ctx:&CallCtx<'_>,a:Handle<RootSystem>,top:Vec<i64>,cap:usize|{limit(cap,10000)?;a.with(ctx,|s|{weight_work(s,&top)?;Ok::<_,OpError>(Maybe(s.all_weights(&top,cap).map(|v|table(v,s.n)).transpose()?))})});
object_binding!(root_system_weight_multiplicity,"num_native_root_system_weight_multiplicity","crate::kernels::lie_structure::RootSystem::weight_multiplicity","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<RootSystem>,top:Vec<i64>,x:Vec<i64>)->Option<BigInt>,|ctx:&CallCtx<'_>,a:Handle<RootSystem>,top:Vec<i64>,x:Vec<i64>|{a.with(ctx,|s|{weight_work(s,&top)?;weight(&x,s.n)?;Ok::<_,OpError>(s.weight_multiplicity(&top,&x))})});
object_binding!(root_system_tensor,"num_native_root_system_tensor","crate::kernels::lie_structure::RootSystem::tensor","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<RootSystem>,x:Vec<i64>,y:Vec<i64>)->Maybe<Dynamic>,|ctx:&CallCtx<'_>,a:Handle<RootSystem>,x:Vec<i64>,y:Vec<i64>|{a.with(ctx,|s|{weight_work(s,&x)?;weight_work(s,&y)?;Ok::<_,OpError>(Maybe(s.tensor(&x,&y).map(|v|table(v,s.n)).transpose()?))})});
object_binding!(lie_from_matrices,"num_native_lie_from_matrices","crate::kernels::lie_structure::Lie::from_matrices","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(basis:Matrices)->Maybe<NewObject<Lie>>,|_ctx:&CallCtx<'_>,basis:Matrices|{Ok::<_,OpError>(Maybe(Lie::from_matrices(&basis.0).map(NewObject::normalized)))});
object_binding!(lie_from_constants,"num_native_lie_from_constants","crate::kernels::lie_structure::Lie::from_constants","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(constants:Input)->Maybe<NewObject<Lie>>,|_ctx:&CallCtx<'_>,constants:Input|{api::normalized(&constants.0)?;let Value::Tensor(t)=&constants.0 else{return Err(error("structure constants require rank3 exact tensor"));};
if t.shape.len()!=3||t.shape[0]==0||t.shape[0]>12||t.shape.iter().any(|n|*n!=t.shape[0]){return Err(error("Lie structure constants shape must be[n,n,n] with1<=n<=12"));}let n=t.shape[0];let values=t.data.iter().map(qmatrix::rational).collect::<Result<Vec<_>,_>>()?;if values.iter().any(|q|q.numer().bits()>1024||q.denom().bits()>1024){return Err(error("Lie input coefficient bit bound exceeded"));}let c=values.chunks(n*n).map(|r|r.chunks(n).map(<[Q]>::to_vec).collect()).collect();Ok::<_,OpError>(Maybe(Lie::from_constants(c).map(NewObject::normalized)))});
object_binding!(lie_is_nilpotent,"num_native_lie_is_nilpotent","crate::kernels::lie_structure::Lie::is_nilpotent","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Lie>)->bool,|ctx:&CallCtx<'_>,a:Handle<Lie>|{a.with(ctx,|s|Ok::<_,OpError>(s.is_nilpotent()))});
object_binding!(lie_is_solvable,"num_native_lie_is_solvable","crate::kernels::lie_structure::Lie::is_solvable","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Lie>)->bool,|ctx:&CallCtx<'_>,a:Handle<Lie>|{a.with(ctx,|s|Ok::<_,OpError>(s.is_solvable()))});
object_binding!(lie_cartan_solvable_test,"num_native_lie_cartan_solvable_test","crate::kernels::lie_structure::Lie::cartan_solvable_test","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Lie>)->bool,|ctx:&CallCtx<'_>,a:Handle<Lie>|{a.with(ctx,|s|Ok::<_,OpError>(s.cartan_solvable_test()))});
object_binding!(lie_is_semisimple,"num_native_lie_is_semisimple","crate::kernels::lie_structure::Lie::is_semisimple","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Lie>)->bool,|ctx:&CallCtx<'_>,a:Handle<Lie>|{a.with(ctx,|s|Ok::<_,OpError>(s.is_semisimple()))});
object_binding!(lie_killing,"num_native_lie_killing","crate::kernels::lie_structure::Lie::killing","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Lie>)->NewObject<QMatrix>,|ctx:&CallCtx<'_>,a:Handle<Lie>|{a.with(ctx,|s|Ok::<_,OpError>(NewObject::normalized(QMatrix::new(s.killing(),s.n))))});
object_binding!(lie_center,"num_native_lie_center","crate::kernels::lie_structure::Lie::center","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Lie>)->NewObject<QMatrix>,|ctx:&CallCtx<'_>,a:Handle<Lie>|{a.with(ctx,|s|Ok::<_,OpError>(NewObject::normalized(QMatrix::new(s.center(),s.n))))});
object_binding!(lie_radical,"num_native_lie_radical","crate::kernels::lie_structure::Lie::radical","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Lie>)->NewObject<QMatrix>,|ctx:&CallCtx<'_>,a:Handle<Lie>|{a.with(ctx,|s|Ok::<_,OpError>(NewObject::normalized(QMatrix::new(s.radical(),s.n))))});
object_binding!(lie_cartan_subalgebra,"num_native_lie_cartan_subalgebra","crate::kernels::lie_structure::Lie::cartan_subalgebra","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Lie>)->NewObject<QMatrix>,|ctx:&CallCtx<'_>,a:Handle<Lie>|{a.with(ctx,|s|Ok::<_,OpError>(NewObject::normalized(QMatrix::new(s.cartan_subalgebra(),s.n))))});
object_binding!(lie_derived_series,"num_native_lie_derived_series","crate::kernels::lie_structure::Lie::derived_series","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Lie>)->Many<NewObject<QMatrix>>,|ctx:&CallCtx<'_>,a:Handle<Lie>|{a.with(ctx,|s|Ok::<_,OpError>(Many(s.derived_series().into_iter().map(|m|NewObject::normalized(QMatrix::new(m,s.n))).collect())))});
object_binding!(lie_lower_central_series,"num_native_lie_lower_central_series","crate::kernels::lie_structure::Lie::lower_central_series","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Lie>)->Many<NewObject<QMatrix>>,|ctx:&CallCtx<'_>,a:Handle<Lie>|{a.with(ctx,|s|Ok::<_,OpError>(Many(s.lower_central_series().into_iter().map(|m|NewObject::normalized(QMatrix::new(m,s.n))).collect())))});
object_binding!(lie_upper_central_series,"num_native_lie_upper_central_series","crate::kernels::lie_structure::Lie::upper_central_series","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Lie>)->Many<NewObject<QMatrix>>,|ctx:&CallCtx<'_>,a:Handle<Lie>|{a.with(ctx,|s|Ok::<_,OpError>(Many(s.upper_central_series().into_iter().map(|m|NewObject::normalized(QMatrix::new(m,s.n))).collect())))});
object_binding!(lie_bracket,"num_native_lie_bracket","crate::kernels::lie_structure::Lie::bracket","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Lie>,x:Input,y:Input)->Vec<Q>,|ctx:&CallCtx<'_>,a:Handle<Lie>,x:Input,y:Input|{a.with(ctx,|s|Ok::<_,OpError>(s.bracket(&qvec(&x,s.n)?,&qvec(&y,s.n)?)))});
object_binding!(lie_ad,"num_native_lie_ad","crate::kernels::lie_structure::Lie::ad","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Lie>,x:Input)->NewObject<QMatrix>,|ctx:&CallCtx<'_>,a:Handle<Lie>,x:Input|{a.with(ctx,|s|Ok::<_,OpError>(NewObject::normalized(QMatrix::new(s.ad(&qvec(&x,s.n)?),s.n))))});
object_binding!(lie_bracket_spaces,"num_native_lie_bracket_spaces","crate::kernels::lie_structure::Lie::bracket_spaces","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Lie>,x:Input,y:Input)->NewObject<QMatrix>,|ctx:&CallCtx<'_>,a:Handle<Lie>,x:Input,y:Input|{a.with(ctx,|s|Ok::<_,OpError>(NewObject::normalized(QMatrix::new(s.bracket_spaces(&qspace(&x,s.n)?,&qspace(&y,s.n)?),s.n))))});
object_binding!(lie_is_subalgebra,"num_native_lie_is_subalgebra","crate::kernels::lie_structure::Lie::is_subalgebra","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Lie>,space:Input)->bool,|ctx:&CallCtx<'_>,a:Handle<Lie>,space:Input|{a.with(ctx,|s|Ok::<_,OpError>(s.is_subalgebra(&qspace(&space,s.n)?)))});
object_binding!(lie_levi,"num_native_lie_levi","crate::kernels::lie_structure::Lie::levi","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Lie>)->Maybe<Dynamic>,|ctx:&CallCtx<'_>,a:Handle<Lie>|{a.with(ctx,|s|Ok::<_,OpError>(Maybe(s.levi().map(|(a,b)|Ok::<_,OpError>(Dynamic(api::record(vec![("semisimple",qmatrix::packed(a,s.n)?.0),("radical",qmatrix::packed(b,s.n)?.0)])))).transpose()?)))});
object_binding!(lie_casimir,"num_native_lie_casimir","crate::kernels::lie_structure::Lie::casimir","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Lie>,representation:Matrices)->Maybe<NewObject<QMatrix>>,|ctx:&CallCtx<'_>,a:Handle<Lie>,representation:Matrices|{a.with(ctx,|s|{if representation.0.len()!=s.n{return Err(error("representation needs one matrix per Lie basis vector"));}let n=representation.0[0].len();Ok::<_,OpError>(Maybe(s.casimir(&representation.0).map(|m|NewObject::normalized(QMatrix::new(m,n)))))})});
object_binding!(lie_root_decomposition,"num_native_lie_root_decomposition","crate::kernels::lie_structure::Lie::root_decomposition","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Lie>,cartan_basis:Input)->Maybe<Dynamic>,|ctx:&CallCtx<'_>,a:Handle<Lie>,cartan_basis:Input|{a.with(ctx,|s|{let h=qspace(&cartan_basis,s.n)?;let out=s.root_decomposition(&h).map(|(spaces,other)|{let spaces=Many(spaces.into_iter().map(|r|Ok::<_,OpError>(Dynamic(api::record(vec![("root",Output::encode(r.root)?),("space",qmatrix::packed(r.space,s.n)?.0)])))).collect::<Result<Vec<_>,OpError>>()?).encode(ctx)?;Ok::<_,OpError>(Dynamic(api::record(vec![("root_spaces",spaces),("other_dimensions",Value::int(other))])))}).transpose()?;Ok::<_,OpError>(Maybe(out))})});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &lie_cartan_matrix::BINDING,
            &lie_parse_type::BINDING,
            &lie_components::BINDING,
            &lie_positive_roots::BINDING,
            &lie_recognize::BINDING,
            &lie_euclidean_simple_roots::BINDING,
            &lie_standard_basis::BINDING,
            &lie_charpoly::BINDING,
            &lie_rational_roots::BINDING,
            &root_system_new::BINDING,
            &root_system_num_positive::BINDING,
            &root_system_dimension::BINDING,
            &root_system_is_irreducible::BINDING,
            &root_system_exponents::BINDING,
            &root_system_weyl_order::BINDING,
            &root_system_highest_root::BINDING,
            &root_system_coxeter_number::BINDING,
            &root_system_dual_coxeter_number::BINDING,
            &root_system_gram::BINDING,
            &root_system_fundamental_weights::BINDING,
            &root_system_weyl_dim::BINDING,
            &root_system_casimir::BINDING,
            &root_system_dynkin_index::BINDING,
            &root_system_inner::BINDING,
            &root_system_dominant_multiplicities::BINDING,
            &root_system_all_weights::BINDING,
            &root_system_weight_multiplicity::BINDING,
            &root_system_tensor::BINDING,
            &lie_from_matrices::BINDING,
            &lie_from_constants::BINDING,
            &lie_is_nilpotent::BINDING,
            &lie_is_solvable::BINDING,
            &lie_cartan_solvable_test::BINDING,
            &lie_is_semisimple::BINDING,
            &lie_killing::BINDING,
            &lie_center::BINDING,
            &lie_radical::BINDING,
            &lie_cartan_subalgebra::BINDING,
            &lie_derived_series::BINDING,
            &lie_lower_central_series::BINDING,
            &lie_upper_central_series::BINDING,
            &lie_bracket::BINDING,
            &lie_ad::BINDING,
            &lie_bracket_spaces::BINDING,
            &lie_is_subalgebra::BINDING,
            &lie_levi::BINDING,
            &lie_casimir::BINDING,
            &lie_root_decomposition::BINDING,
        ],
        "numerical_vector",
    )
}
