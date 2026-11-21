#![allow(clippy::wildcard_imports)] // Enclosures retain physical dimensions; polynomial data has explicit normalized-coordinate references.
use super::*;
use crate::backend::tape::Interval;
use crate::backend::tape::Scalar;
use crate::kernels::ball::CTerm;
use crate::kernels::ball::Cball;
use crate::kernels::ball::{
    self,
};
use crate::kernels::rootfind::PolyTerm;
use num_complex::Complex64;
use num_rational::BigRational;
use num_traits::ToPrimitive;
pub(super) struct PolynomialSystem {
    pub(super) real: Vec<Vec<PolyTerm>>,
    pub(super) complex: Vec<Vec<CTerm>>,
    pub(super) coordinate: f64,
    pub(super) residual: f64,
    pub(super) n: usize,
}
impl NativeType for PolynomialSystem {
    const KIND: &'static str = "rssn/native/discrete/polynomial-system/v1";
    const MAX_BYTES: usize = 1_048_576;

    fn validate(&self) -> Result<(), OpError> {
        limit(self.n, 6)?;
        if self.n == 0
            || self.complex.len() != self.n
            || self.coordinate <= 0.0
            || self.residual <= 0.0
            || !self.coordinate.is_finite()
            || !self.residual.is_finite()
        {
            return Err(error("invalid square polynomial system references"));
        }
        let mut count = 0usize;
        let mut paths = 1usize;
        for row in &self.complex {
            count += row.len();
            let mut degree = 0u32;
            for t in row {
                if t.exps.len() != self.n || !t.coef.re.is_finite() || !t.coef.im.is_finite() {
                    return Err(error("invalid polynomial term"));
                }
                degree = degree.max(
                    t.exps
                        .iter()
                        .try_fold(0u32, |s, e| s.checked_add(*e))
                        .ok_or_else(|| error("polynomial degree overflow"))?,
                );
            }
            if degree == 0 || degree > 16 {
                return Err(error("polynomial row degree must be1..16"));
            }
            paths = paths
                .checked_mul(degree as usize)
                .filter(|n| *n <= 256)
                .ok_or_else(|| error("polynomial total degree path budget exceeds256"))?;
        }
        limit(count, 4096)?;
        Ok(())
    }
}
impl NativeType for Cball {
    const KIND: &'static str = "rssn/native/discrete/complex-ball/v1";
    const MAX_BYTES: usize = 64;

    fn validate(&self) -> Result<(), OpError> {
        if !self.is_finite() || self.re.lo > self.re.hi || self.im.lo > self.im.hi {
            return Err(error("complex ball must be finite with ordered endpoints"));
        }
        Ok(())
    }
}
impl NativeType for CTerm {
    const KIND: &'static str = "rssn/native/discrete/complex-term/v1";
    const MAX_BYTES: usize = 128;

    fn validate(&self) -> Result<(), OpError> {
        limit(self.exps.len(), 7)?;
        if !self.coef.re.is_finite()
            || !self.coef.im.is_finite()
            || self.exps.iter().any(|e| *e > 16)
        {
            return Err(error("invalid normalized complex polynomial term"));
        }
        Ok(())
    }
}
struct Balls(Vec<Handle<Cball>>);
impl NativeArg for Balls {
    fn ty() -> TypeExpr {
        api::tensor_type(1, TypeExpr::Opaque(Cball::KIND.into()))
    }

    fn decode(
        v: &Value,
        _: usize,
        c: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        let Value::Tensor(t) = v else {
            return Err(error("expected vector of ball handles"));
        };
        if t.shape != vec![t.data.len()] || t.data.len() > 6 {
            return Err(error("complex box dimension budget is6"));
        }
        Ok(Self(
            t.data
                .iter()
                .map(|v| <Handle<Cball> as NativeArg>::decode(v, 0, c))
                .collect::<Result<_, _>>()?,
        ))
    }
}
fn system_dims(
    s: &Handle<PolynomialSystem>,
    ctx: &CallCtx<'_>,
) -> Result<(DimExpr, DimExpr), OpError> {
    let d = s.dimensions(ctx)?;
    if d.len() != 2 {
        return Err(error("polynomial reference dimension witnesses missing"));
    }
    Ok((d[0].clone(), d[1].clone()))
}
fn complex_input(v: &Input) -> Result<(Complex64, DimExpr), OpError> {
    let d = api::dimension(&v.0)?;
    let z = <Complex64 as Arg>::decode(&v.0, 0)?;
    if !z.re.is_finite() || !z.im.is_finite() {
        return Err(error("complex point must be finite"));
    }
    Ok((z, d))
}
fn box_values(
    v: &Balls,
    s: &PolynomialSystem,
    d: &DimExpr,
    ctx: &CallCtx<'_>,
) -> Result<Vec<Cball>, OpError> {
    if v.0.len() != s.n {
        return Err(error("polynomial box dimension mismatch"));
    }
    v.0.iter()
        .map(|h| {
            api::same_dimension(&object_dimension(h, ctx)?, d)?;
            h.with(ctx, |b| {
                b.div(&Cball::point(s.coordinate, 0.0))
                    .ok_or_else(|| error("coordinate normalization failed"))
            })
        })
        .collect()
}
fn point_values(
    v: &Input,
    s: &PolynomialSystem,
    d: &DimExpr,
) -> Result<Vec<Complex64>, OpError> {
    api::same_dimension(&api::dimension(&v.0)?, d)?;
    let x = <Vec<Complex64> as Arg>::decode(&v.0, 0)?;
    if x.len() != s.n {
        return Err(error("polynomial point dimension mismatch"));
    }
    let x = x.into_iter().map(|v| v / s.coordinate).collect::<Vec<_>>();
    if x.iter().any(|z| !z.re.is_finite() || !z.im.is_finite()) {
        return Err(error("normalized polynomial point must be finite"));
    }
    Ok(x)
}
fn ball_output(
    v: Cball,
    scale: f64,
    d: DimExpr,
) -> NewObject<Cball> {
    NewObject::new(v * Cball::point(scale, 0.0), vec![d])
}
fn complex_output(
    v: Vec<Complex64>,
    scale: f64,
    d: DimExpr,
) -> Result<Dynamic, OpError> {
    let values = v
        .into_iter()
        .map(|z| {
            let value = Output::encode(z * scale)?;
            if d == DimExpr::Dimless {
                Ok(value)
            } else {
                Ok(Value::Dyn(crate::kernel::value::DynValue {
                    dim: d.clone(),
                    value: Box::new(value),
                }))
            }
        })
        .collect::<Result<Vec<_>, OpError>>()?;
    Ok(Dynamic(Value::Tensor(TensorValue {
        shape: vec![values.len()],
        data: values,
        elem: TypeExpr::Complex.with_dim(d),
    })))
}
object_binding!(cball_point,"num_native_cball_point","crate::kernels::ball::Cball::point","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(re:Input,im:Input)->NewObject<Cball>,|_ctx:&CallCtx<'_>,re:Input,im:Input|{let(r,d)=scalar(&re)?;let(i,e)=scalar(&im)?;api::same_dimension(&d,&e)?;Ok(NewObject::new(Cball::point(r,i),vec![d]))});
object_binding!(cball_of,"num_native_cball_of","crate::kernels::ball::Cball::of","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(point:Input)->NewObject<Cball>,|_ctx:&CallCtx<'_>,point:Input|{let(z,d)=complex_input(&point)?;Ok(NewObject::new(Cball::of(z),vec![d]))});
object_binding!(cball_square,"num_native_cball_square","crate::kernels::ball::Cball::square","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(re:Input,im:Input,radius:Input)->NewObject<Cball>,|_ctx:&CallCtx<'_>,re:Input,im:Input,radius:Input|{let(r,d)=scalar(&re)?;let(i,e)=scalar(&im)?;let(w,f)=scalar(&radius)?;api::same_dimension(&d,&e)?;api::same_dimension(&d,&f)?;if w<0.0{return Err(error("complex ball radius must be nonnegative"));}Ok(NewObject::new(Cball::square(r,i,w),vec![d]))});
object_binding!(cball_rectangle,"num_native_cball_rectangle","crate::kernels::ball::Cball::rectangle","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(re:Input,im:Input)->NewObject<Cball>,|_ctx:&CallCtx<'_>,re:Input,im:Input|{let d=api::dimension(&re.0)?;api::same_dimension(&d,&api::dimension(&im.0)?)?;let r=<[f64;2] as Arg>::decode(&re.0,0)?;let i=<[f64;2] as Arg>::decode(&im.0,1)?;if r[0]>r[1]||i[0]>i[1]{return Err(error("complex rectangle endpoints must increase"));}Ok(NewObject::new(Cball::rectangle((r[0],r[1]),(i[0],i[1])),vec![d]))});
object_binding!(cball_zero,"num_native_cball_zero","crate::kernels::ball::Cball::zero","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,()->NewObject<Cball>,|_ctx:&CallCtx<'_>|{Ok(NewObject::normalized(Cball::zero()))});
object_binding!(cball_one,"num_native_cball_one","crate::kernels::ball::Cball::one","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,()->NewObject<Cball>,|_ctx:&CallCtx<'_>|{Ok(NewObject::normalized(Cball::one()))});
object_binding!(cball_centre,"num_native_cball_centre","crate::kernels::ball::Cball::centre","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Cball>)->Dynamic,|ctx:&CallCtx<'_>,a:Handle<Cball>|{let d=object_dimension(&a,ctx)?;a.with(ctx,|b|{let value=Output::encode(b.centre())?;Ok(Dynamic(if d==DimExpr::Dimless{value}else{Value::Dyn(crate::kernel::value::DynValue{dim:d,value:Box::new(value)})}))})});
object_binding!(cball_half_width,"num_native_cball_half_width","crate::kernels::ball::Cball::half_width","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Cball>)->Dynamic,|ctx:&CallCtx<'_>,a:Handle<Cball>|{let d=object_dimension(&a,ctx)?;a.with(ctx,|b|extended_quantity(b.half_width(),d))});
object_binding!(cball_modulus_lower,"num_native_cball_modulus_lower","crate::kernels::ball::Cball::modulus_lower","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Cball>)->Dynamic,|ctx:&CallCtx<'_>,a:Handle<Cball>|{let d=object_dimension(&a,ctx)?;a.with(ctx,|b|extended_quantity(b.modulus_lower(),d))});
object_binding!(cball_modulus_upper,"num_native_cball_modulus_upper","crate::kernels::ball::Cball::modulus_upper","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Cball>)->Dynamic,|ctx:&CallCtx<'_>,a:Handle<Cball>|{let d=object_dimension(&a,ctx)?;a.with(ctx,|b|extended_quantity(b.modulus_upper(),d))});
object_binding!(cball_is_finite,"num_native_cball_is_finite","crate::kernels::ball::Cball::is_finite","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Cball>)->bool,|ctx:&CallCtx<'_>,a:Handle<Cball>|{a.with(ctx,|b|Ok(b.is_finite()))});
object_binding!(cball_contains,"num_native_cball_contains","crate::kernels::ball::Cball::contains","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Cball>,point:Input)->bool,|ctx:&CallCtx<'_>,a:Handle<Cball>,point:Input|{let(z,d)=complex_input(&point)?;api::same_dimension(&d,&object_dimension(&a,ctx)?)?;a.with(ctx,|b|Ok(b.contains(z)))});
object_binding!(cball_inside_strictly,"num_native_cball_inside_strictly","crate::kernels::ball::Cball::inside_strictly","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Cball>,b:Handle<Cball>)->bool,|ctx:&CallCtx<'_>,a:Handle<Cball>,b:Handle<Cball>|{same(&a,&b,ctx)?;a.with(ctx,|x|b.with(ctx,|y|Ok(x.inside_strictly(y))))});
object_binding!(cball_intersect,"num_native_cball_intersect","crate::kernels::ball::Cball::intersect","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Cball>,b:Handle<Cball>)->Maybe<NewObject<Cball>>,|ctx:&CallCtx<'_>,a:Handle<Cball>,b:Handle<Cball>|{let d=same(&a,&b,ctx)?;a.with(ctx,|x|b.with(ctx,|y|Ok(Maybe(x.intersect(y).map(|v|NewObject::new(v,vec![d]))))))});
object_binding!(cball_div,"num_native_cball_div","crate::kernels::ball::Cball::div","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Cball>,b:Handle<Cball>)->Maybe<NewObject<Cball>>,|ctx:&CallCtx<'_>,a:Handle<Cball>,b:Handle<Cball>|{let d=canonical(&(DimExpr::Mul(vec![object_dimension(&a,ctx)?,DimExpr::Inverse(Box::new(object_dimension(&b,ctx)?))])))?;a.with(ctx,|x|b.with(ctx,|y|Ok(Maybe(x.div(y).map(|v|NewObject::new(v,vec![d]))))))});
object_binding!(cball_powi,"num_native_cball_powi","crate::kernels::ball::Cball::powi","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Cball>,exponent:u32)->NewObject<Cball>,|ctx:&CallCtx<'_>,a:Handle<Cball>,exponent:u32|{limit(exponent as usize,1024)?;let d=canonical(&(object_dimension(&a,ctx)?.pow(Exp::Rat(Rational::int(i64::from(exponent))))))?;a.with(ctx,|b|Ok(NewObject::new(b.powi(exponent),vec![d])))});
object_binding!(cterm_real,"num_native_cterm_real","crate::kernels::ball::CTerm::real","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(coefficient:f64,exponents:Vec<u32>)->NewObject<CTerm>,|_ctx:&CallCtx<'_>,coefficient:f64,exponents:Vec<u32>|{limit(exponents.len(),6)?;if exponents.iter().any(|e|*e>16){return Err(error("normalized polynomial exponent budget is16"));}Ok(NewObject::normalized(CTerm::real(coefficient,exponents)))});
object_binding!(polynomial_system_new,"num_native_polynomial_system_new","","Normalized polynomial coefficients describe F(x/coordinate_reference)/residual_reference; two positive references preserve independent physical coordinate and residual dimensions. Maximum6variables/4096terms/256total-degree paths.",false,(coefficients:Vec<f64>,exponents:Vec<Vec<u32>>,offsets:Vec<usize>,coordinate_reference:Input,residual_reference:Input)->NewObject<PolynomialSystem>,|_ctx:&CallCtx<'_>,coefficients:Vec<f64>,exponents:Vec<Vec<u32>>,offsets:Vec<usize>,coordinate_reference:Input,residual_reference:Input|{let(cr,d)=api::reference(&coordinate_reference.0)?;let(vr,e)=api::reference(&residual_reference.0)?;let n=offsets.len().checked_sub(1).filter(|n|*n>0&&*n<=6).ok_or_else(||error("square polynomial descriptor needs1..6 equations"))?;if coefficients.len()!=exponents.len()||offsets.first()!=Some(&0)||offsets.last()!=Some(&coefficients.len())||offsets.windows(2).any(|w|w[0]>w[1]){return Err(error("invalid polynomial offsets"));}limit(coefficients.len(),4096)?;let terms=coefficients.into_iter().zip(exponents).map(|(coef,exps)|PolyTerm{coef,exps}).collect::<Vec<_>>();let real=offsets.windows(2).map(|w|terms[w[0]..w[1]].to_vec()).collect::<Vec<_>>();let complex=ball::real_terms(&real);let system=PolynomialSystem{real,complex,coordinate:cr,residual:vr,n};system.validate()?;Ok(NewObject::new(system,vec![d,e]))});
object_binding!(ball_real_terms,"num_native_ball_real_terms","crate::kernels::ball::real_terms","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(system:Handle<PolynomialSystem>)->NewObject<PolynomialSystem>,|ctx:&CallCtx<'_>,system:Handle<PolynomialSystem>|{let dims=system.dimensions(ctx)?;system.with(ctx,|s|Ok(NewObject::new(PolynomialSystem{real:s.real.clone(),complex:ball::real_terms(&s.real),coordinate:s.coordinate,residual:s.residual,n:s.n},dims)))});
object_binding!(ball_eval_terms,"num_native_ball_eval_terms","crate::kernels::ball::eval_terms","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(system:Handle<PolynomialSystem>,point:Balls)->Many<NewObject<Cball>>,|ctx:&CallCtx<'_>,system:Handle<PolynomialSystem>,point:Balls|{let(cd,rd)=system_dims(&system,ctx)?;system.with(ctx,|s|{let x=box_values(&point,s,&cd,ctx)?;let d=rd;let scale=s.residual;Ok(Many(ball::eval_terms(&s.complex,&x).into_iter().map(|v|ball_output(v,scale,d.clone())).collect()))})});
object_binding!(ball_eval_polys,"num_native_ball_eval_polys","crate::kernels::ball::eval_polys","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(system:Handle<PolynomialSystem>,point:Balls)->Many<NewObject<Cball>>,|ctx:&CallCtx<'_>,system:Handle<PolynomialSystem>,point:Balls|{let(cd,rd)=system_dims(&system,ctx)?;system.with(ctx,|s|{let x=box_values(&point,s,&cd,ctx)?;let d=rd;let scale=s.residual;Ok(Many(ball::eval_polys(&s.real,&x).into_iter().map(|v|ball_output(v,scale,d.clone())).collect()))})});
object_binding!(ball_eval_jacobian_terms,"num_native_ball_eval_jacobian_terms","crate::kernels::ball::eval_jacobian_terms","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(system:Handle<PolynomialSystem>,point:Balls)->Many<Many<NewObject<Cball>>>,|ctx:&CallCtx<'_>,system:Handle<PolynomialSystem>,point:Balls|{let(cd,rd)=system_dims(&system,ctx)?;system.with(ctx,|s|{let x=box_values(&point,s,&cd,ctx)?;let d=canonical(&(DimExpr::Mul(vec![rd,DimExpr::Inverse(Box::new(cd))])))?;let scale=s.residual/s.coordinate;Ok(Many(ball::eval_jacobian_terms(&s.complex,&x).into_iter().map(|row|Many(row.into_iter().map(|v|ball_output(v,scale,d.clone())).collect())).collect()))})});
object_binding!(ball_eval_jacobian,"num_native_ball_eval_jacobian","crate::kernels::ball::eval_jacobian","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(system:Handle<PolynomialSystem>,point:Balls)->Many<Many<NewObject<Cball>>>,|ctx:&CallCtx<'_>,system:Handle<PolynomialSystem>,point:Balls|{let(cd,rd)=system_dims(&system,ctx)?;system.with(ctx,|s|{let x=box_values(&point,s,&cd,ctx)?;let d=canonical(&(DimExpr::Mul(vec![rd,DimExpr::Inverse(Box::new(cd))])))?;let scale=s.residual/s.coordinate;Ok(Many(ball::eval_jacobian(&s.real,&x).into_iter().map(|row|Many(row.into_iter().map(|v|ball_output(v,scale,d.clone())).collect())).collect()))})});
object_binding!(ball_values_at_terms,"num_native_ball_values_at_terms","crate::kernels::ball::values_at_terms","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(system:Handle<PolynomialSystem>,point:Input)->Dynamic,|ctx:&CallCtx<'_>,system:Handle<PolynomialSystem>,point:Input|{let(cd,rd)=system_dims(&system,ctx)?;system.with(ctx,|s|{let x=point_values(&point,s,&cd)?;let d=rd;let scale=s.residual;complex_output(ball::values_at_terms(&s.complex,&x),scale,d)})});
object_binding!(ball_values_at,"num_native_ball_values_at","crate::kernels::ball::values_at","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(system:Handle<PolynomialSystem>,point:Input)->Dynamic,|ctx:&CallCtx<'_>,system:Handle<PolynomialSystem>,point:Input|{let(cd,rd)=system_dims(&system,ctx)?;system.with(ctx,|s|{let x=point_values(&point,s,&cd)?;let d=rd;let scale=s.residual;complex_output(ball::values_at(&s.real,&x),scale,d)})});
object_binding!(ball_jacobian_at_terms,"num_native_ball_jacobian_at_terms","crate::kernels::ball::jacobian_at_terms","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(system:Handle<PolynomialSystem>,point:Input)->Many<Dynamic>,|ctx:&CallCtx<'_>,system:Handle<PolynomialSystem>,point:Input|{let(cd,rd)=system_dims(&system,ctx)?;system.with(ctx,|s|{let x=point_values(&point,s,&cd)?;let d=canonical(&(DimExpr::Mul(vec![rd,DimExpr::Inverse(Box::new(cd))])))?;let scale=s.residual/s.coordinate;Ok(Many(ball::jacobian_at_terms(&s.complex,&x).into_iter().map(|v|complex_output(v,scale,d.clone())).collect::<Result<_,_>>()?))})});
object_binding!(ball_jacobian_at,"num_native_ball_jacobian_at","crate::kernels::ball::jacobian_at","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(system:Handle<PolynomialSystem>,point:Input)->Many<Dynamic>,|ctx:&CallCtx<'_>,system:Handle<PolynomialSystem>,point:Input|{let(cd,rd)=system_dims(&system,ctx)?;system.with(ctx,|s|{let x=point_values(&point,s,&cd)?;let d=canonical(&(DimExpr::Mul(vec![rd,DimExpr::Inverse(Box::new(cd))])))?;let scale=s.residual/s.coordinate;Ok(Many(ball::jacobian_at(&s.real,&x).into_iter().map(|v|complex_output(v,scale,d.clone())).collect::<Result<_,_>>()?))})});
struct BallEvaluation {
    values: Many<NewObject<Cball>>,
    jacobian: Vec<NewObject<Cball>>,
    size: usize,
}

impl NativeOutput for BallEvaluation {
    fn ty() -> TypeExpr {
        let ball = TypeExpr::Opaque(Cball::KIND.into());
        TypeExpr::Record(vec![
            ("values".into(), api::tensor_type(1, ball.clone())),
            ("jacobian".into(), api::tensor_type(2, ball)),
        ])
    }

    fn encode(
        self,
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        let scope = api::bank(ctx)?.scoped()?;
        let values = self.values.encode(ctx)?;
        let jacobian = self
            .jacobian
            .into_iter()
            .map(|cell| cell.encode(ctx))
            .collect::<Result<Vec<_>, _>>()?;
        let jacobian = Value::Tensor(TensorValue {
            shape: vec![self.size, self.size],
            data: jacobian,
            elem: TypeExpr::Opaque(Cball::KIND.into()),
        });
        let output = api::record(vec![("values", values), ("jacobian", jacobian)]);
        scope.commit();
        Ok(output)
    }
}

#[cold]
fn exact_scaled_endpoint(
    value: f64,
    ratio: &BigRational,
    upper: bool,
) -> Result<f64, OpError> {
    let exact = BigRational::from_float(value)
        .ok_or_else(|| error("nonfinite polynomial enclosure"))?
        * ratio;
    let rounded = exact
        .to_f64()
        .ok_or_else(|| error("polynomial enclosure exceeds Float range"))?;
    let encoded = BigRational::from_float(rounded)
        .ok_or_else(|| error("polynomial enclosure exceeds Float range"))?;
    Ok(if upper && encoded < exact {
        rounded.next_up()
    } else if !upper && encoded > exact {
        rounded.next_down()
    } else {
        rounded
    })
}

enum JacobianScale {
    Outward(Interval),
    Exact(BigRational),
}

impl JacobianScale {
    fn new(
        residual: f64,
        coordinate: f64,
    ) -> Result<Self, OpError> {
        let scale = residual / coordinate;
        let upper = scale.next_up();
        if scale.is_normal() && upper.is_finite() {
            return Ok(Self::Outward(Interval {
                lo: scale.next_down(),
                hi: upper,
            }));
        }
        // A reference ratio outside the normal range can still give
        // representable endpoints when applied to the normalized enclosure.
        let residual = BigRational::from_float(residual)
            .ok_or_else(|| error("invalid polynomial residual reference"))?;
        let coordinate = BigRational::from_float(coordinate)
            .ok_or_else(|| error("invalid polynomial coordinate reference"))?;
        Ok(Self::Exact(residual / coordinate))
    }

    fn restore(
        &self,
        cell: Cball,
        dimension: DimExpr,
    ) -> Result<NewObject<Cball>, OpError> {
        let interval = |part: Interval| -> Result<Interval, OpError> {
            match self {
                | Self::Outward(scale) => Ok(part.mul(*scale)),
                | Self::Exact(scale) => {
                    Ok(Interval {
                        lo: exact_scaled_endpoint(part.lo, scale, false)?,
                        hi: exact_scaled_endpoint(part.hi, scale, true)?,
                    })
                },
            }
        };
        let value = Cball {
            re: interval(cell.re)?,
            im: interval(cell.im)?,
        };
        Ok(NewObject::new(value, vec![dimension]))
    }
}

object_binding!(
    ball_eval_value_and_jacobian_terms,
    "num_native_ball_eval_value_and_jacobian_terms",
    "crate::kernels::ball::eval_value_and_jacobian_terms",
    "Evaluates a polynomial system and its Jacobian over one complex box, sharing term products. Returns {values, jacobian} as enclosure handles. Value units equal the residual unit. Jacobian units are residual divided by the coordinate unit.",
    false,
    (system: Handle<PolynomialSystem>, point: Balls) -> BallEvaluation,
    |ctx, system, point| {
        let (coordinate, residual) = system_dims(&system, ctx)?;
        let derivative = canonical(&DimExpr::Mul(vec![residual.clone(), coordinate.clone().pow(-1)]))?;
        system.with(ctx, |s| {
            let point = box_values(&point, s, &coordinate, ctx)?;
            let (values, jacobian) = ball::eval_value_and_jacobian_terms(&s.complex, &point);
            let scale = JacobianScale::new(s.residual, s.coordinate)?;
            Ok(BallEvaluation {
                values: Many(values.into_iter().map(|cell| ball_output(cell, s.residual, residual.clone())).collect()),
                jacobian: jacobian.into_iter().flatten().map(|cell| scale.restore(cell, derivative.clone())).collect::<Result<Vec<_>, _>>()?,
                size: s.n,
            })
        })
    }
);

object_binding!(
    ball_values_and_jacobian_at_terms,
    "num_native_ball_values_and_jacobian_at_terms",
    "crate::kernels::ball::values_and_jacobian_at_terms",
    "Evaluates a polynomial system and its Jacobian at one complex point. Returns {values, jacobian} with residual units and residual-per-coordinate units respectively. The Jacobian is a rank-two tensor with one row per equation.",
    false,
    (system: Handle<PolynomialSystem>, point: Input) -> Dynamic,
    |ctx, system, point| {
        let (coordinate, residual) = system_dims(&system, ctx)?;
        let derivative = canonical(&DimExpr::Mul(vec![residual.clone(), coordinate.clone().pow(-1)]))?;
        system.with(ctx, |s| {
            let point = point_values(&point, s, &coordinate)?;
            let (values, jacobian) = ball::values_and_jacobian_at_terms(&s.complex, &point);
            let Dynamic(values) = complex_output(values, s.residual, residual)?;
            let matrix = jacobian.into_iter().flatten().map(|value| {
                Output::encode(Complex64::new(
                    api::scale_product(value.re, &[s.residual], &[s.coordinate])?,
                    api::scale_product(value.im, &[s.residual], &[s.coordinate])?,
                ))
            }).collect::<Result<Vec<_>, OpError>>()?;
            let jacobian = Value::Tensor(TensorValue {
                shape: vec![s.n, s.n], data: matrix,
                elem: TypeExpr::Complex.with_dim(derivative),
            });
            Ok(Dynamic(api::record(vec![("values", values), ("jacobian", jacobian)])))
        })
    }
);

pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &cball_point::BINDING,
            &cball_of::BINDING,
            &cball_square::BINDING,
            &cball_rectangle::BINDING,
            &cball_zero::BINDING,
            &cball_one::BINDING,
            &cball_centre::BINDING,
            &cball_half_width::BINDING,
            &cball_modulus_lower::BINDING,
            &cball_modulus_upper::BINDING,
            &cball_is_finite::BINDING,
            &cball_contains::BINDING,
            &cball_inside_strictly::BINDING,
            &cball_intersect::BINDING,
            &cball_div::BINDING,
            &cball_powi::BINDING,
            &cterm_real::BINDING,
            &polynomial_system_new::BINDING,
            &ball_real_terms::BINDING,
            &ball_eval_terms::BINDING,
            &ball_eval_polys::BINDING,
            &ball_eval_jacobian_terms::BINDING,
            &ball_eval_jacobian::BINDING,
            &ball_values_at_terms::BINDING,
            &ball_values_at::BINDING,
            &ball_jacobian_at_terms::BINDING,
            &ball_jacobian_at::BINDING,
            &ball_eval_value_and_jacobian_terms::BINDING,
            &ball_values_and_jacobian_at_terms::BINDING,
        ],
        "numerical_vector",
    )
}
