//! Exact native holonomic objects over explicitly normalized rational coordinates.
use crate::graph::RuleError;
use crate::graph::rule::Installer;
use crate::kernel::OpError;
use crate::kernels::holonomic as h;
use crate::rules::numerical_abi::codecs::Output;
use crate::rules::numerical_abi::native_object_api::Dynamic;
use crate::rules::numerical_abi::native_object_api::Handle;
use crate::rules::numerical_abi::native_object_api::Many;
use crate::rules::numerical_abi::native_object_api::Maybe;
use crate::rules::numerical_abi::native_object_api::NativeType;
use crate::rules::numerical_abi::native_object_api::NewObject;
use crate::rules::numerical_abi::native_object_api::object_binding;
use crate::rules::numerical_abi::native_object_api::{
    self as api,
};
use num_rational::BigRational as Q;
#[derive(Clone)]
struct Rational(h::RatFn);
#[derive(Clone)]
struct Atom(h::Atom);
#[derive(Clone)]
struct Expr(h::Expr);
#[derive(Clone)]
struct Operator(h::Operator);
fn poly(p: &[Q]) -> Result<(), OpError> {
    if p.len() > 65
        || p.iter()
            .any(|q| q.numer().bits() > 2048 || q.denom().bits() > 2048)
    {
        Err(api::error(
            "holonomic polynomial degree64/coefficient2048bit budget exceeded",
        ))
    } else {
        Ok(())
    }
}
fn rat(r: &h::RatFn) -> Result<(), OpError> {
    poly(&r.num)?;
    poly(&r.den)
}
fn atom(a: &h::Atom) -> Result<(), OpError> {
    match a {
        | h::Atom::Power(r, q) => {
            poly(std::slice::from_ref(q))?;
            rat(r)
        },
        | h::Atom::Rational(r)
        | h::Atom::Exp(r)
        | h::Atom::Log(r)
        | h::Atom::Atan(r)
        | h::Atom::Atanh(r)
        | h::Atom::Asin(r) => rat(r),
        | h::Atom::Sin(p) | h::Atom::Cos(p) | h::Atom::Sinh(p) | h::Atom::Cosh(p) => poly(p),
    }
}
fn expression(e: &h::Expr) -> Result<(), OpError> {
    let mut stack = vec![(e, 0)];
    let mut count = 0;
    while let Some((e, d)) = stack.pop() {
        count += 1;
        if count > 128 || d > 24 {
            return Err(api::error(
                "holonomic expression node/depth budget exceeded",
            ));
        }
        match e {
            | h::Expr::Atom(a) => atom(a)?,
            | h::Expr::Add(a, b) | h::Expr::Mul(a, b) => {
                stack.push((a, d + 1));
                stack.push((b, d + 1));
            },
            | h::Expr::Scale(q, a) => {
                poly(std::slice::from_ref(q))?;
                stack.push((a, d + 1));
            },
            | h::Expr::Pow(a, n) => {
                if *n > 16 {
                    return Err(api::error("holonomic power budget exceeds16"));
                }
                stack.push((a, d + 1));
            },
            | h::Expr::Diff(a) => stack.push((a, d + 1)),
        }
    }
    Ok(())
}
impl NativeType for Rational {
    const KIND: &'static str = "rssn/native/analysis/holonomic/rational/v1";
    const MAX_BYTES: usize = 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        rat(&self.0)
    }
}
impl NativeType for Atom {
    const KIND: &'static str = "rssn/native/analysis/holonomic/atom/v1";
    const MAX_BYTES: usize = 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        atom(&self.0)
    }
}
impl NativeType for Expr {
    const KIND: &'static str = "rssn/native/analysis/holonomic/expression/v1";
    const MAX_BYTES: usize = 8 * 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        expression(&self.0)
    }
}
impl NativeType for Operator {
    const KIND: &'static str = "rssn/native/analysis/holonomic/operator/v1";
    const MAX_BYTES: usize = 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        if self.0.len() > 9 {
            return Err(api::error("holonomic operator order budget exceeds8"));
        }
        for p in &self.0 {
            poly(p)?;
        }
        Ok(())
    }
}
#[allow(
    clippy::needless_pass_by_value,
    reason = "The typed native ABI transfers decoded arguments by value; this shared dispatcher mirrors its explicit source signature."
)]
fn limits(v: Vec<usize>) -> Result<h::Limits, OpError> {
    if v.len() != 3 || v[0] == 0 || v[0] > 8 || v[1] > 64 || v[2] == 0 || v[2] > 512 {
        return Err(api::error(
            "holonomic limits require order1..8,degree0..64,terms1..512",
        ));
    }
    Ok(h::Limits {
        max_order: v[0],
        max_degree: v[1],
        max_terms: v[2],
    })
}
fn terms(n: usize) -> Result<(), OpError> {
    if n > 512 {
        Err(api::error("holonomic series term budget exceeds512"))
    } else {
        Ok(())
    }
}
object_binding!(rational_new,"num_analysis_holonomic_rational_new","crate::kernels::holonomic::RatFn::new","Native exact holonomic operation in explicitly normalized rational coordinates. Session/kind and bounded exact coefficient/expression work are checked.",true,(num:Vec<Q>,den:Vec<Q>)->Maybe<NewObject<Rational>>,|_ctx,num,den|{poly(&num)?;poly(&den)?;Ok(Maybe(h::RatFn::new(num,den).map(|r|NewObject::normalized(Rational(r))))) });
object_binding!(rational_constant,"num_analysis_holonomic_rational_constant","crate::kernels::holonomic::RatFn::constant","Native exact holonomic operation in explicitly normalized rational coordinates. Session/kind and bounded exact coefficient/expression work are checked.",true,(q:Q)->NewObject<Rational>,|_ctx,q|{poly(std::slice::from_ref(&q))?;Ok(NewObject::normalized(Rational(h::RatFn::constant(q))))});
object_binding!(rational_poly,"num_analysis_holonomic_rational_poly","crate::kernels::holonomic::RatFn::poly","Native exact holonomic operation in explicitly normalized rational coordinates. Session/kind and bounded exact coefficient/expression work are checked.",true,(p:Vec<Q>)->NewObject<Rational>,|_ctx,p|{poly(&p)?;Ok(NewObject::normalized(Rational(h::RatFn::poly(p))))});
object_binding!(rational_zero,"num_analysis_holonomic_rational_zero","crate::kernels::holonomic::RatFn::zero","Native exact holonomic operation in explicitly normalized rational coordinates. Session/kind and bounded exact coefficient/expression work are checked.",true,()->NewObject<Rational>,|_ctx|Ok(NewObject::normalized(Rational(h::RatFn::zero()))));
object_binding!(rational_add,"num_analysis_holonomic_rational_add","crate::kernels::holonomic::RatFn::add","Native exact holonomic operation in explicitly normalized rational coordinates. Session/kind and bounded exact coefficient/expression work are checked.",true,(a:Handle<Rational>,b:Handle<Rational>)->NewObject<Rational>,|ctx,a,b|{let b=b.with(ctx,|b|Ok(b.0.clone()))?;a.with(ctx,|a|Ok(NewObject::normalized(Rational(a.0.add(&b)))))});
object_binding!(rational_sub,"num_analysis_holonomic_rational_sub","crate::kernels::holonomic::RatFn::sub","Native exact holonomic operation in explicitly normalized rational coordinates. Session/kind and bounded exact coefficient/expression work are checked.",true,(a:Handle<Rational>,b:Handle<Rational>)->NewObject<Rational>,|ctx,a,b|{let b=b.with(ctx,|b|Ok(b.0.clone()))?;a.with(ctx,|a|Ok(NewObject::normalized(Rational(a.0.sub(&b)))))});
object_binding!(rational_mul,"num_analysis_holonomic_rational_mul","crate::kernels::holonomic::RatFn::mul","Native exact holonomic operation in explicitly normalized rational coordinates. Session/kind and bounded exact coefficient/expression work are checked.",true,(a:Handle<Rational>,b:Handle<Rational>)->NewObject<Rational>,|ctx,a,b|{let b=b.with(ctx,|b|Ok(b.0.clone()))?;a.with(ctx,|a|Ok(NewObject::normalized(Rational(a.0.mul(&b)))))});
object_binding!(rational_div,"num_analysis_holonomic_rational_div","crate::kernels::holonomic::RatFn::div","Native exact holonomic operation in explicitly normalized rational coordinates. Session/kind and bounded exact coefficient/expression work are checked.",true,(a:Handle<Rational>,b:Handle<Rational>)->Maybe<NewObject<Rational>>,|ctx,a,b|{let b=b.with(ctx,|b|Ok(b.0.clone()))?;a.with(ctx,|a|Ok(Maybe(a.0.div(&b).map(|v|NewObject::normalized(Rational(v))))))});
object_binding!(rational_neg,"num_analysis_holonomic_rational_neg","crate::kernels::holonomic::RatFn::neg","Native exact holonomic operation in explicitly normalized rational coordinates. Session/kind and bounded exact coefficient/expression work are checked.",true,(a:Handle<Rational>)->NewObject<Rational>,|ctx,a|a.with(ctx,|a|Ok(NewObject::normalized(Rational(a.0.neg())))));
object_binding!(rational_deriv,"num_analysis_holonomic_rational_deriv","crate::kernels::holonomic::RatFn::deriv","Native exact holonomic operation in explicitly normalized rational coordinates. Session/kind and bounded exact coefficient/expression work are checked.",true,(a:Handle<Rational>)->NewObject<Rational>,|ctx,a|a.with(ctx,|a|Ok(NewObject::normalized(Rational(a.0.deriv())))));
object_binding!(rational_is_zero,"num_analysis_holonomic_rational_is_zero","crate::kernels::holonomic::RatFn::is_zero","Native exact holonomic operation in explicitly normalized rational coordinates. Session/kind and bounded exact coefficient/expression work are checked.",true,(a:Handle<Rational>)->bool,|ctx,a|a.with(ctx,|a|Ok(a.0.is_zero())));
object_binding!(rational_series,"num_analysis_holonomic_rational_series","crate::kernels::holonomic::RatFn::series","Native exact holonomic operation in explicitly normalized rational coordinates. Session/kind and bounded exact coefficient/expression work are checked.",true,(a:Handle<Rational>,n:usize)->Maybe<Vec<Q>>,|ctx,a,n|{terms(n)?;a.with(ctx,|a|Ok(Maybe(a.0.series(n))))});
object_binding!(atom_new,"num_analysis_holonomic_atom_new","","Native exact holonomic operation in explicitly normalized rational coordinates. Session/kind and bounded exact coefficient/expression work are checked.",true,(kind:String,rational:Handle<Rational>,polynomial:Vec<Q>,power:Q)->NewObject<Atom>,|ctx,kind,rational,polynomial,power|{poly(&polynomial)?;poly(std::slice::from_ref(&power))?;let r=rational.with(ctx,|r|Ok(r.0.clone()))?;let a=match kind.as_str(){"Rational"=>h::Atom::Rational(r),"Exp"=>h::Atom::Exp(r),"Power"=>h::Atom::Power(r,power),"Log"=>h::Atom::Log(r),"Sin"=>h::Atom::Sin(polynomial),"Cos"=>h::Atom::Cos(polynomial),"Sinh"=>h::Atom::Sinh(polynomial),"Cosh"=>h::Atom::Cosh(polynomial),"Atan"=>h::Atom::Atan(r),"Atanh"=>h::Atom::Atanh(r),"Asin"=>h::Atom::Asin(r),_=>return Err(api::error("unknown holonomic Atom variant"))};Ok(NewObject::normalized(Atom(a)))});
object_binding!(atom_is_valid,"num_analysis_holonomic_atom_is_valid","crate::kernels::holonomic::Atom::is_valid","Native exact holonomic operation in explicitly normalized rational coordinates. Session/kind and bounded exact coefficient/expression work are checked.",true,(a:Handle<Atom>)->bool,|ctx,a|a.with(ctx,|a|Ok(a.0.is_valid())));
object_binding!(atom_series,"num_analysis_holonomic_atom_series","crate::kernels::holonomic::Atom::series","Native exact holonomic operation in explicitly normalized rational coordinates. Session/kind and bounded exact coefficient/expression work are checked.",true,(a:Handle<Atom>,n:usize)->Maybe<Vec<Q>>,|ctx,a,n|{terms(n)?;a.with(ctx,|a|Ok(Maybe(a.0.series(n))))});
object_binding!(atom_annihilator,"num_analysis_holonomic_atom_annihilator","crate::kernels::holonomic::Atom::annihilator","Native exact holonomic operation in explicitly normalized rational coordinates. Session/kind and bounded exact coefficient/expression work are checked.",true,(a:Handle<Atom>)->Maybe<NewObject<Operator>>,|ctx,a|a.with(ctx,|a|Ok(Maybe(a.0.annihilator().map(|p|NewObject::normalized(Operator(p)))))));
object_binding!(expression_atom,"num_analysis_holonomic_expression_atom","","Native exact holonomic operation in explicitly normalized rational coordinates. Session/kind and bounded exact coefficient/expression work are checked.",true,(a:Handle<Atom>)->NewObject<Expr>,|ctx,a|a.with(ctx,|a|Ok(NewObject::normalized(Expr(h::Expr::Atom(a.0.clone()))))));
object_binding!(expression_binary,"num_analysis_holonomic_expression_binary","","Native exact holonomic operation in explicitly normalized rational coordinates. Session/kind and bounded exact coefficient/expression work are checked.",true,(a:Handle<Expr>,b:Handle<Expr>,product:bool)->NewObject<Expr>,|ctx,a,b,product|{let b=b.with(ctx,|b|Ok(b.0.clone()))?;a.with(ctx,|a|Ok(NewObject::normalized(Expr(if product{h::Expr::Mul(Box::new(a.0.clone()),Box::new(b))}else{h::Expr::Add(Box::new(a.0.clone()),Box::new(b))}))))});
object_binding!(expression_unary,"num_analysis_holonomic_expression_unary","","Native exact holonomic operation in explicitly normalized rational coordinates. Session/kind and bounded exact coefficient/expression work are checked.",true,(a:Handle<Expr>,kind:String,q:Q,n:u32)->NewObject<Expr>,|ctx,a,kind,q,n|{poly(std::slice::from_ref(&q))?;if n>16{return Err(api::error("holonomic expression exponent exceeds16"));}a.with(ctx,|a|{let a=Box::new(a.0.clone());Ok(NewObject::normalized(Expr(match kind.as_str(){"Scale"=>h::Expr::Scale(q,a),"Pow"=>h::Expr::Pow(a,n),"Diff"=>h::Expr::Diff(a),_=>return Err(api::error("unknown holonomic Expr unary variant"))})))})});
object_binding!(expression_series,"num_analysis_holonomic_expression_series","crate::kernels::holonomic::Expr::series","Native exact holonomic operation in explicitly normalized rational coordinates. Session/kind and bounded exact coefficient/expression work are checked.",true,(a:Handle<Expr>,n:usize)->Maybe<Vec<Q>>,|ctx,a,n|{terms(n)?;a.with(ctx,|a|Ok(Maybe(a.0.series(n))))});
object_binding!(annihilator,"num_analysis_holonomic_annihilator","crate::kernels::holonomic::annihilator","Native exact holonomic operation in explicitly normalized rational coordinates. Session/kind and bounded exact coefficient/expression work are checked.",true,(a:Handle<Expr>,bounds:Vec<usize>)->Maybe<NewObject<Operator>>,|ctx,a,bounds|{let limits=limits(bounds)?;a.with(ctx,|a|Ok(Maybe(h::annihilator(&a.0,&limits).map(|p|NewObject::normalized(Operator(p))))))});
object_binding!(is_zero,"num_analysis_holonomic_is_zero","crate::kernels::holonomic::is_zero","Native exact holonomic operation in explicitly normalized rational coordinates. Session/kind and bounded exact coefficient/expression work are checked.",true,(a:Handle<Expr>,bounds:Vec<usize>)->Maybe<bool>,|ctx,a,bounds|{let limits=limits(bounds)?;a.with(ctx,|a|Ok(Maybe(h::is_zero(&a.0,&limits))))});
object_binding!(deciding_terms,"num_analysis_holonomic_deciding_terms","crate::kernels::holonomic::deciding_terms","Native exact holonomic operation in explicitly normalized rational coordinates. Session/kind and bounded exact coefficient/expression work are checked.",true,(a:Handle<Operator>,bounds:Vec<usize>)->Maybe<usize>,|ctx,a,bounds|{let limits=limits(bounds)?;a.with(ctx,|a|Ok(Maybe(h::deciding_terms(&a.0,&limits))))});
object_binding!(recurrence,"num_analysis_holonomic_recurrence","crate::kernels::holonomic::recurrence","Native exact holonomic operation in explicitly normalized rational coordinates. Session/kind and bounded exact coefficient/expression work are checked.",true,(a:Handle<Operator>)->Many<Dynamic>,|ctx,a|a.with(ctx,|a|Ok(Many(h::recurrence(&a.0).into_iter().map(|(shift,p)|Ok(Dynamic(api::record(vec![("shift",Output::encode(shift)?),("polynomial",Output::encode(p)?)])))).collect::<Result<Vec<_>,OpError>>()?))));
object_binding!(operator_data,"num_analysis_holonomic_operator_data","","Native exact holonomic operation in explicitly normalized rational coordinates. Session/kind and bounded exact coefficient/expression work are checked.",true,(a:Handle<Operator>)->Many<Vec<Q>>,|ctx,a|a.with(ctx,|a|Ok(Many(a.0.clone()))));
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &rational_new::BINDING,
            &rational_constant::BINDING,
            &rational_poly::BINDING,
            &rational_zero::BINDING,
            &rational_add::BINDING,
            &rational_sub::BINDING,
            &rational_mul::BINDING,
            &rational_div::BINDING,
            &rational_neg::BINDING,
            &rational_deriv::BINDING,
            &rational_is_zero::BINDING,
            &rational_series::BINDING,
            &atom_new::BINDING,
            &atom_is_valid::BINDING,
            &atom_series::BINDING,
            &atom_annihilator::BINDING,
            &expression_atom::BINDING,
            &expression_binary::BINDING,
            &expression_unary::BINDING,
            &expression_series::BINDING,
            &annihilator::BINDING,
            &is_zero::BINDING,
            &deciding_terms::BINDING,
            &recurrence::BINDING,
            &operator_data::BINDING,
        ],
        "numerical_vector",
    )
}
