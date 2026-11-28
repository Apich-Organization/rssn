#![allow(clippy::wildcard_imports)] // Explicitly scaled polynomial homotopy and native proof objects.
use super::ball::PolynomialSystem;
use super::*;
use crate::kernels::ball::CTerm;
use crate::kernels::ball::Cball;
use crate::kernels::homotopy::HomotopyOptions;
use crate::kernels::homotopy::PathEnd;
use crate::kernels::homotopy::Precision;
use crate::kernels::homotopy::{
    self,
};
use crate::kernels::homotopy_cert::Certificate;
use crate::kernels::homotopy_cert::CertificationReport;
use crate::kernels::homotopy_cert::RigorousBall;
use crate::kernels::homotopy_cert::RigorousReport;
use crate::kernels::homotopy_cert::{
    self as cert,
};
use num_complex::Complex64;
struct Ends {
    paths: Vec<PathEnd>,
    coordinate: f64,
    n: usize,
}
struct Proof {
    value: RigorousBall,
    coordinate: f64,
}
impl NativeType for Ends {
    const KIND: &'static str = "rssn/native/discrete/path-ends/v1";
    const MAX_BYTES: usize = 1_048_576;

    fn validate(&self) -> Result<(), OpError> {
        limit(self.paths.len(), 256)?;
        limit(self.n, 6)?;
        if self.coordinate <= 0.0
            || !self.coordinate.is_finite()
            || self.paths.iter().any(|p| {
                p.x.len() != self.n || p.x_lo.len() != self.n || p.homogeneous.len() != self.n + 1
            })
        {
            return Err(error("invalid native path endpoint record"));
        }
        Ok(())
    }
}
impl NativeType for Proof {
    const KIND: &'static str = "rssn/native/discrete/rigorous-zero/v1";
    const MAX_BYTES: usize = 1024;

    fn validate(&self) -> Result<(), OpError> {
        limit(self.value.enclosure.len(), 7)?;
        if self.value.enclosure.len() != self.value.centre.len()
            || self.coordinate <= 0.0
            || !self.coordinate.is_finite()
            || !self.value.radius.is_finite()
            || self.value.radius < 0.0
        {
            return Err(error("invalid rigorous zero proof"));
        }
        for b in &self.value.enclosure {
            b.validate()?;
        }
        Ok(())
    }
}
fn options(v: &Input) -> Result<HomotopyOptions, OpError> {
    api::normalized(&v.0)?;
    let f = |key| numeric(api::field(&v.0, key)?, 0);
    let precision = match <String as Arg>::decode(api::field(&v.0, "precision")?, 0)?.as_str() {
        | "Double" => Precision::Double,
        | "DoubleDouble" => Precision::DoubleDouble,
        | "Adaptive" => Precision::Adaptive,
        | _ => return Err(error("precision must be Double, DoubleDouble, or Adaptive")),
    };
    let o = HomotopyOptions {
        seed: <u64 as Arg>::decode(api::field(&v.0, "seed")?, 0)?,
        endgame_start: f("endgame_start")?,
        endgame_tol: f("endgame_tol")?,
        max_step: f("max_step")?,
        max_retries: <usize as Arg>::decode(api::field(&v.0, "max_retries")?, 0)?,
        infinity_tol: f("infinity_tol")?,
        precision,
        certify_steps: <bool as Arg>::decode(api::field(&v.0, "certify_steps")?, 0)?,
    };
    if !(0.000001..0.5).contains(&o.endgame_start)
        || o.endgame_tol <= 0.0
        || o.max_step < 0.00001
        || o.max_step > 0.5
        || o.max_retries > 8
        || o.infinity_tol <= 0.0
    {
        return Err(error(
            "homotopy options exceed positive tolerance/step/retry work bounds",
        ));
    }
    Ok(o)
}
fn dims(
    s: &Handle<PolynomialSystem>,
    ctx: &CallCtx<'_>,
) -> Result<DimExpr, OpError> {
    s.dimensions(ctx)?
        .first()
        .cloned()
        .ok_or_else(|| error("missing coordinate dimension"))
}
fn guess(
    v: &Input,
    s: &PolynomialSystem,
    d: &DimExpr,
) -> Result<Vec<(f64, f64)>, OpError> {
    api::same_dimension(&api::dimension(&v.0)?, d)?;
    let x = <Vec<Complex64> as Arg>::decode(&v.0, 0)?;
    if x.len() != s.n {
        return Err(error("homotopy guess dimension differs from system"));
    }
    x.into_iter()
        .map(|z| {
            let z = z / s.coordinate;
            if !z.re.is_finite() || !z.im.is_finite() {
                return Err(error("normalized homotopy point not finite"));
            }
            Ok((z.re, z.im))
        })
        .collect()
}
fn certificate(
    c: Certificate,
    r: f64,
    d: DimExpr,
) -> Result<Value, OpError> {
    Ok(api::record(vec![
        ("beta", extended_quantity(c.beta * r, d.clone())?.0),
        (
            "gamma",
            extended_quantity(
                c.gamma / r,
                canonical(&(DimExpr::Inverse(Box::new(d.clone()))))?,
            )?
            .0,
        ),
        ("alpha", Output::encode(c.alpha)?),
        ("certified", Value::Bool(c.certified)),
        ("radius", extended_quantity(c.radius * r, d)?.0),
    ]))
}
fn floating(
    report: CertificationReport,
    r: f64,
    d: &DimExpr,
    ctx: &CallCtx<'_>,
) -> Result<Value, OpError> {
    let certificates = Many(
        report
            .certificates
            .into_iter()
            .map(|c| {
                c.map(|c| certificate(c, r, d.clone()).map(Dynamic))
                    .transpose()
                    .map(Maybe)
            })
            .collect::<Result<Vec<_>, _>>()?,
    )
    .encode(ctx)?;
    Ok(api::record(vec![
        ("certificates", certificates),
        ("certified", Value::int(report.certified)),
        ("distinct", Value::int(report.distinct)),
        ("paths", Value::int(report.paths)),
        ("complete", Value::Bool(report.complete)),
    ]))
}
fn rigorous(
    report: RigorousReport,
    r: f64,
    d: &DimExpr,
    ctx: &CallCtx<'_>,
) -> Result<Value, OpError> {
    let certificates = Many(
        report
            .certificates
            .into_iter()
            .map(|p| {
                Maybe(
                    p.map(|value| NewObject::new(Proof { value, coordinate: r }, vec![d.clone()])),
                )
            })
            .collect(),
    )
    .encode(ctx)?;
    Ok(api::record(vec![
        ("certificates", certificates),
        ("certified", Value::int(report.certified)),
        ("distinct", Value::int(report.distinct)),
        ("paths", Value::int(report.paths)),
        ("complete", Value::Bool(report.complete)),
    ]))
}
fn compatible(
    s: &PolynomialSystem,
    e: &Ends,
) -> Result<(), OpError> {
    if s.n != e.n || !crate::kernels::exact_eq(s.coordinate, e.coordinate) {
        return Err(error(
            "endpoint coordinate reference differs from polynomial system",
        ));
    }
    Ok(())
}
fn homogenized(s: &PolynomialSystem) -> (Vec<Vec<crate::kernels::rootfind::PolyTerm>>, Vec<u32>) {
    let degree = s
        .real
        .iter()
        .map(|row| {
            row.iter()
                .map(|t| t.exps.iter().sum::<u32>())
                .max()
                .unwrap_or(0)
        })
        .collect::<Vec<_>>();
    let rows = s
        .real
        .iter()
        .zip(&degree)
        .map(|(row, d)| {
            row.iter()
                .map(|t| {
                    let mut exps = vec![d - t.exps.iter().sum::<u32>()];
                    exps.extend_from_slice(&t.exps);
                    crate::kernels::rootfind::PolyTerm { coef: t.coef, exps }
                })
                .collect()
        })
        .collect();
    (rows, degree)
}
object_binding!(homotopy_solve,"num_native_homotopy_solve","crate::kernels::homotopy::homotopy_solve","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(system:Handle<PolynomialSystem>,config:Input)->NewObject<Ends>,|ctx:&CallCtx<'_>,system:Handle<PolynomialSystem>,config:Input|{let d=dims(&system,ctx)?;let o=options(&config)?;system.with(ctx,|s|{let paths=homotopy::homotopy_solve(&s.real,&o).map_err(|e|error(format!("{e:?}")))?;Ok(NewObject::new(Ends{paths,coordinate:s.coordinate,n:s.n},vec![d]))})});
object_binding!(homotopy_ends_data,"num_native_homotopy_ends_data","","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(ends:Handle<Ends>)->Many<Dynamic>,|ctx:&CallCtx<'_>,ends:Handle<Ends>|{let d=object_dimension(&ends,ctx)?;ends.with(ctx,|e|{Ok(Many(e.paths.iter().map(|p|{let coordinates=p.x.iter().zip(&p.x_lo).map(|(&(re,im),&(lo_re,lo_im))|{let z=Complex64::new((re+lo_re)*e.coordinate,(im+lo_im)*e.coordinate);if z.re.is_finite()&&z.im.is_finite(){let v=Output::encode(z)?;Ok(Maybe(Some(Dynamic(if d==DimExpr::Dimless{v}else{Value::Dyn(crate::kernel::value::DynValue{dim:d.clone(),value:Box::new(v)})}))))}else{Ok(Maybe(None))}}).collect::<Result<Vec<Maybe<Dynamic>>,OpError>>()?;Ok(Dynamic(api::record(vec![("coordinates",Many(coordinates).encode(ctx)?),("at_infinity",Value::Bool(p.at_infinity)),("extended",Value::Bool(p.extended)),("singular",Value::Bool(p.singular)),("converged",Value::Bool(p.converged)),("suspect",Value::Bool(p.suspect)),("winding",Value::int(p.winding)),("retries",Value::int(p.retries)),("residual",Output::encode(p.residual.is_finite().then_some(p.residual))?),("verified_steps",Value::int(p.verified_steps)),("unverified_steps",Value::int(p.unverified_steps)),("steps_certified",Value::Bool(p.steps_certified))]))) }).collect::<Result<_,OpError>>()?))})});
object_binding!(homotopy_smale_alpha,"num_native_homotopy_smale_alpha","crate::kernels::homotopy_cert::smale_alpha","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(system:Handle<PolynomialSystem>,point:Input)->Maybe<Dynamic>,|ctx:&CallCtx<'_>,system:Handle<PolynomialSystem>,point:Input|{let d=dims(&system,ctx)?;system.with(ctx,|s|{let p=guess(&point,s,&d)?;Ok(Maybe(cert::smale_alpha(&s.real,&p).map(|c|certificate(c,s.coordinate,d).map(Dynamic)).transpose()?))})});
object_binding!(homotopy_verify_zero,"num_native_homotopy_verify_zero","crate::kernels::homotopy_cert::verify_zero","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(system:Handle<PolynomialSystem>,point:Input,relative:f64)->Maybe<NewObject<Proof>>,|ctx:&CallCtx<'_>,system:Handle<PolynomialSystem>,point:Input,relative:f64|{let d=dims(&system,ctx)?;if relative<=0.0||relative>1.0{return Err(error("verification relative radius must be in(0,1]"));}system.with(ctx,|s|{let p=guess(&point,s,&d)?;Ok(Maybe(cert::verify_zero(&s.real,&p,relative).map(|value|NewObject::new(Proof{value,coordinate:s.coordinate},vec![d]))))})});
object_binding!(homotopy_certify_zero,"num_native_homotopy_certify_zero","crate::kernels::homotopy_cert::certify_zero","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(system:Handle<PolynomialSystem>,point:Input)->Maybe<NewObject<Proof>>,|ctx:&CallCtx<'_>,system:Handle<PolynomialSystem>,point:Input|{let d=dims(&system,ctx)?;system.with(ctx,|s|{let p=guess(&point,s,&d)?;Ok(Maybe(cert::certify_zero(&s.real,&p).map(|value|NewObject::new(Proof{value,coordinate:s.coordinate},vec![d]))))})});
object_binding!(homotopy_verify_zero_terms,"num_native_homotopy_verify_zero_terms","crate::kernels::homotopy_cert::verify_zero_terms","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(system:Handle<PolynomialSystem>,point:Input,relative:f64)->Maybe<NewObject<Proof>>,|ctx:&CallCtx<'_>,system:Handle<PolynomialSystem>,point:Input,relative:f64|{let d=dims(&system,ctx)?;if relative<=0.0||relative>1.0{return Err(error("verification relative radius must be in(0,1]"));}system.with(ctx,|s|{let p=guess(&point,s,&d)?;Ok(Maybe(cert::verify_zero_terms(&s.complex,&p,relative).map(|value|NewObject::new(Proof{value,coordinate:s.coordinate},vec![d]))))})});
object_binding!(homotopy_certify_zero_terms,"num_native_homotopy_certify_zero_terms","crate::kernels::homotopy_cert::certify_zero_terms","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(system:Handle<PolynomialSystem>,point:Input)->Maybe<NewObject<Proof>>,|ctx:&CallCtx<'_>,system:Handle<PolynomialSystem>,point:Input|{let d=dims(&system,ctx)?;system.with(ctx,|s|{let p=guess(&point,s,&d)?;Ok(Maybe(cert::certify_zero_terms(&s.complex,&p).map(|value|NewObject::new(Proof{value,coordinate:s.coordinate},vec![d]))))})});
object_binding!(homotopy_certify_endpoints,"num_native_homotopy_certify_endpoints","crate::kernels::homotopy_cert::certify_endpoints","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(system:Handle<PolynomialSystem>,ends:Handle<Ends>,safety:f64)->Dynamic,|ctx:&CallCtx<'_>,system:Handle<PolynomialSystem>,ends:Handle<Ends>,safety:f64|{let d=dims(&system,ctx)?;api::same_dimension(&d,&object_dimension(&ends,ctx)?)?;if safety<1.0{return Err(error("certification safety must be>=1"));}system.with(ctx,|s|ends.with(ctx,|e|{compatible(s,e)?;let report=cert::certify_endpoints(&s.real,&e.paths,safety);Ok(Dynamic(floating(report,s.coordinate,&(d),ctx)?))}))});
object_binding!(homotopy_certify_endpoints_rigorous,"num_native_homotopy_certify_endpoints_rigorous","crate::kernels::homotopy_cert::certify_endpoints_rigorous","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(system:Handle<PolynomialSystem>,ends:Handle<Ends>)->Dynamic,|ctx:&CallCtx<'_>,system:Handle<PolynomialSystem>,ends:Handle<Ends>|{let d=dims(&system,ctx)?;api::same_dimension(&d,&object_dimension(&ends,ctx)?)?;system.with(ctx,|s|ends.with(ctx,|e|{compatible(s,e)?;let report=cert::certify_endpoints_rigorous(&s.real,&e.paths);Ok(Dynamic(rigorous(report,s.coordinate,&(d),ctx)?))}))});
object_binding!(homotopy_homotopy_solve_certified,"num_native_homotopy_homotopy_solve_certified","crate::kernels::homotopy_cert::homotopy_solve_certified","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(system:Handle<PolynomialSystem>,config:Input)->Dynamic,|ctx:&CallCtx<'_>,system:Handle<PolynomialSystem>,config:Input|{let d=dims(&system,ctx)?;let o=options(&config)?;system.with(ctx,|s|{let out=cert::homotopy_solve_certified(&s.real,&o).map_err(|e|error(format!("{e:?}")))?;let ends=NewObject::new(Ends{paths:out.ends,coordinate:s.coordinate,n:s.n},vec![d.clone()]).encode(ctx)?;let report=floating(out.report,s.coordinate,&(d),ctx)?;Ok(Dynamic(api::record(vec![("ends",ends),("report",report)])))})});
object_binding!(homotopy_homotopy_solve_verified,"num_native_homotopy_homotopy_solve_verified","crate::kernels::homotopy_cert::homotopy_solve_verified","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(system:Handle<PolynomialSystem>,config:Input)->Dynamic,|ctx:&CallCtx<'_>,system:Handle<PolynomialSystem>,config:Input|{let d=dims(&system,ctx)?;let o=options(&config)?;system.with(ctx,|s|{let out=cert::homotopy_solve_verified(&s.real,&o).map_err(|e|error(format!("{e:?}")))?;let ends=NewObject::new(Ends{paths:out.ends,coordinate:s.coordinate,n:s.n},vec![d.clone()]).encode(ctx)?;let report=floating(out.floating,s.coordinate,&d,ctx)?;let rigorous=rigorous(out.rigorous,s.coordinate,&(d),ctx)?;Ok(Dynamic(api::record(vec![("ends",ends),("floating",report),("rigorous",rigorous)])))})});
object_binding!(rigorous_ball_contains,"num_native_rigorous_ball_contains","crate::kernels::homotopy_cert::RigorousBall::contains","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(proof:Handle<Proof>,re:Input,im:Input)->bool,|ctx:&CallCtx<'_>,proof:Handle<Proof>,re:Input,im:Input|{let(r,d)=scalar(&re)?;let(i,e)=scalar(&im)?;api::same_dimension(&d,&e)?;api::same_dimension(&d,&object_dimension(&proof,ctx)?)?;proof.with(ctx,|p|Ok(p.value.contains(r/p.coordinate,i/p.coordinate)))});
object_binding!(rigorous_ball_distance,"num_native_rigorous_ball_distance","crate::kernels::homotopy_cert::RigorousBall::distance","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(proof:Handle<Proof>,re:Input,im:Input)->Dynamic,|ctx:&CallCtx<'_>,proof:Handle<Proof>,re:Input,im:Input|{let(r,d)=scalar(&re)?;let(i,e)=scalar(&im)?;api::same_dimension(&d,&e)?;api::same_dimension(&d,&object_dimension(&proof,ctx)?)?;proof.with(ctx,|p|extended_quantity(p.value.distance(r/p.coordinate,i/p.coordinate)*p.coordinate,d))});
object_binding!(homotopy_disjoint,"num_native_homotopy_disjoint","crate::kernels::homotopy_cert::disjoint","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Proof>,b:Handle<Proof>)->bool,|ctx:&CallCtx<'_>,a:Handle<Proof>,b:Handle<Proof>|{same(&a,&b,ctx)?;a.with(ctx,|a|b.with(ctx,|b|{if !crate::kernels::exact_eq(a.coordinate,b.coordinate)||a.value.enclosure.len()!=b.value.enclosure.len(){return Err(error("proofs require identical coordinate reference and dimension"));}Ok(cert::disjoint(&a.value,&b.value))}))});
object_binding!(rigorous_ball_enclosure,"num_native_rigorous_ball_enclosure","","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(proof:Handle<Proof>)->Many<NewObject<Cball>>,|ctx:&CallCtx<'_>,proof:Handle<Proof>|{let d=object_dimension(&proof,ctx)?;proof.with(ctx,|p|Ok(Many(p.value.enclosure.iter().map(|b|NewObject::new(*b*Cball::point(p.coordinate,0.0),vec![d.clone()])).collect())))});
object_binding!(homotopy_node_terms,"num_native_homotopy_node_terms","crate::kernels::homotopy_cert::homotopy_node_terms","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(system:Handle<PolynomialSystem>,gamma:Complex64,patch:Vec<Complex64>,time:Complex64)->Many<Many<NewObject<CTerm>>>,|ctx:&CallCtx<'_>,system:Handle<PolynomialSystem>,gamma:Complex64,patch:Vec<Complex64>,time:Complex64|{system.with(ctx,|s|{if patch.len()!=s.n+1{return Err(error("homogeneous patch needs n+1 coefficients"));}let(h,d)=homogenized(s);let rows=cert::homotopy_node_terms(&h,&d,gamma,&patch,time).ok_or_else(||error("invalid homotopy node"))?;Ok(Many(rows.into_iter().map(|row|Many(row.into_iter().map(NewObject::normalized).collect())).collect()))})});
object_binding!(homotopy_certify_homotopy_node,"num_native_homotopy_certify_homotopy_node","crate::kernels::homotopy_cert::certify_homotopy_node","Native rigorous homotopy node in explicitly normalized homogeneous projective coordinates. Patch coefficients, gamma and complex continuation time are dimensionless; this proof object retains dimensionless projective coordinates.",true,(system:Handle<PolynomialSystem>,gamma:Complex64,patch:Vec<Complex64>,point:Vec<Complex64>,time:Complex64)->Maybe<NewObject<Proof>>,|ctx:&CallCtx<'_>,system:Handle<PolynomialSystem>,gamma:Complex64,patch:Vec<Complex64>,point:Vec<Complex64>,time:Complex64|{system.with(ctx,|s|{if patch.len()!=s.n+1||point.len()!=s.n+1{return Err(error("homogeneous patch/point needs n+1 coordinates"));}let(h,d)=homogenized(s);let x=point.into_iter().map(|z|(z.re,z.im)).collect::<Vec<_>>();Ok(Maybe(cert::certify_homotopy_node(&h,&d,gamma,&patch,&x,(time.re,time.im)).map(|value|NewObject::normalized(Proof{value,coordinate:1.0}))))})});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &homotopy_solve::BINDING,
            &homotopy_ends_data::BINDING,
            &homotopy_smale_alpha::BINDING,
            &homotopy_verify_zero::BINDING,
            &homotopy_certify_zero::BINDING,
            &homotopy_verify_zero_terms::BINDING,
            &homotopy_certify_zero_terms::BINDING,
            &homotopy_certify_endpoints::BINDING,
            &homotopy_certify_endpoints_rigorous::BINDING,
            &homotopy_homotopy_solve_certified::BINDING,
            &homotopy_homotopy_solve_verified::BINDING,
            &rigorous_ball_contains::BINDING,
            &rigorous_ball_distance::BINDING,
            &homotopy_disjoint::BINDING,
            &rigorous_ball_enclosure::BINDING,
            &homotopy_node_terms::BINDING,
            &homotopy_certify_homotopy_node::BINDING,
        ],
        "numerical_vector",
    )
}
