#![allow(clippy::wildcard_imports)] // Checked family method declarations.
use super::*;
use crate::kernels::dd::Cdd;
use crate::kernels::dd::Dd;
use crate::kernels::fractal_geometry_and_chaos::AffineTransform2D;
use crate::kernels::fractal_geometry_and_chaos::FractalData;
use crate::kernels::geometric_algebra::Multivector3D;
impl NativeType for Dd {
    const KIND: &'static str = "rssn/native/discrete/dd/v1";
    const MAX_BYTES: usize = 64;

    fn validate(&self) -> Result<(), OpError> {
        if !self.is_finite() {
            return Err(error("non-finite double-double object"));
        }
        Ok(())
    }
}
impl NativeType for Cdd {
    const KIND: &'static str = "rssn/native/discrete/cdd/v1";
    const MAX_BYTES: usize = 128;

    fn validate(&self) -> Result<(), OpError> {
        if !self.re.is_finite() || !self.im.is_finite() {
            return Err(error("non-finite complex double-double object"));
        }
        Ok(())
    }
}
impl NativeType for Multivector3D {
    const KIND: &'static str = "rssn/native/discrete/multivector3/v1";
    const MAX_BYTES: usize = 128;

    fn validate(&self) -> Result<(), OpError> {
        if [
            self.s, self.v1, self.v2, self.v3, self.b12, self.b23, self.b31, self.pss,
        ]
        .iter()
        .any(|v| !v.is_finite())
        {
            return Err(error("multivector coefficients must be finite"));
        }
        Ok(())
    }
}
impl NativeType for AffineTransform2D {
    const KIND: &'static str = "rssn/native/discrete/affine2/v1";
    const MAX_BYTES: usize = 128;

    fn validate(&self) -> Result<(), OpError> {
        if [self.a, self.b, self.c, self.d, self.e, self.f]
            .iter()
            .any(|v| !v.is_finite())
        {
            return Err(error("affine coefficients must be finite"));
        }
        Ok(())
    }
}
impl NativeType for FractalData {
    const KIND: &'static str = "rssn/native/discrete/fractal-data/v1";
    const MAX_BYTES: usize = 1048576;

    fn validate(&self) -> Result<(), OpError> {
        if self.width.checked_mul(self.height) != Some(self.data.len()) || self.data.len() > 262144
        {
            return Err(error("fractal grid exceeds retained shape budget"));
        }
        Ok(())
    }
}
object_binding!(dd_from_f64,"num_native_dd_from_f64","crate::kernels::dd::Dd::from_f64","Construct native double-double from a finite physical scalar; coherent values and physical dimensions are retained.",false,(x:Input)->NewObject<Dd>,|_ctx:&CallCtx<'_>,x:Input|{let(x,d)=scalar(&x)?;Ok(NewObject::new(Dd::from_f64(x),vec![d]))});
object_binding!(dd_from_parts,"num_native_dd_from_parts","crate::kernels::dd::Dd::from_parts","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(hi:Input,lo:Input)->NewObject<Dd>,|_ctx:&CallCtx<'_>,hi:Input,lo:Input|{let(h,d)=scalar(&hi)?;let(l,e)=scalar(&lo)?;api::same_dimension(&d,&e)?;Ok(NewObject::new(Dd::from_parts(h,l),vec![d]))});
object_binding!(dd_abs,"num_native_dd_abs","crate::kernels::dd::Dd::abs","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Dd>)->NewObject<Dd>,|ctx:&CallCtx<'_>,a:Handle<Dd>|{let d=object_dimension(&a,ctx)?;a.with(ctx,|o|Ok(NewObject::new(o.abs(),vec![d])))});
object_binding!(dd_sqrt,"num_native_dd_sqrt","crate::kernels::dd::Dd::sqrt","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Dd>)->NewObject<Dd>,|ctx:&CallCtx<'_>,a:Handle<Dd>|{let d=canonical(&(object_dimension(&a,ctx)?.pow(Exp::Rat(Rational::new(1,2)))))?;a.with(ctx,|o|Ok(NewObject::new(o.sqrt(),vec![d])))});
object_binding!(dd_is_finite,"num_native_dd_is_finite","crate::kernels::dd::Dd::is_finite","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Dd>)->bool,|ctx:&CallCtx<'_>,a:Handle<Dd>|{a.with(ctx,|o|Ok(o.is_finite()))});
object_binding!(dd_to_f64,"num_native_dd_to_f64","crate::kernels::dd::Dd::to_f64","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Dd>)->Dynamic,|ctx:&CallCtx<'_>,a:Handle<Dd>|{let d=object_dimension(&a,ctx)?;a.with(ctx,|o|Ok(Dynamic(api::quantity(o.to_f64(),d)?)))});
object_binding!(cdd_new,"num_native_cdd_new","crate::kernels::dd::Cdd::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(re:Input,im:Input)->NewObject<Cdd>,|_ctx:&CallCtx<'_>,re:Input,im:Input|{let(r,d)=scalar(&re)?;let(i,e)=scalar(&im)?;api::same_dimension(&d,&e)?;Ok(NewObject::new(Cdd::new(r,i),vec![d]))});
object_binding!(cdd_norm,"num_native_cdd_norm","crate::kernels::dd::Cdd::norm","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Cdd>)->Dynamic,|ctx:&CallCtx<'_>,a:Handle<Cdd>|{let d=object_dimension(&a,ctx)?;a.with(ctx,|o|Ok(Dynamic(api::quantity(o.norm(),d)?)))});
object_binding!(multivector_new,"num_native_multivector_new","crate::kernels::geometric_algebra::Multivector3D::new","Construct a homogeneous-coefficient Cl(3,0) multivector; the eight basis blades are dimensionless and all coefficients share physical units.",false,(s:Input,v1:Input,v2:Input,v3:Input,b12:Input,b23:Input,b31:Input,pss:Input)->NewObject<Multivector3D>,|_ctx:&CallCtx<'_>,s:Input,v1:Input,v2:Input,v3:Input,b12:Input,b23:Input,b31:Input,pss:Input|{let values=[scalar(&s)?,scalar(&v1)?,scalar(&v2)?,scalar(&v3)?,scalar(&b12)?,scalar(&b23)?,scalar(&b31)?,scalar(&pss)?];let d=values[0].1.clone();for(_,e)in &values{api::same_dimension(&d,e)?;}Ok(NewObject::new(Multivector3D::new(values[0].0,values[1].0,values[2].0,values[3].0,values[4].0,values[5].0,values[6].0,values[7].0),vec![d]))});
object_binding!(multivector_reverse,"num_native_multivector_reverse","crate::kernels::geometric_algebra::Multivector3D::reverse","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Multivector3D>)->NewObject<Multivector3D>,|ctx:&CallCtx<'_>,a:Handle<Multivector3D>|{let d=object_dimension(&a,ctx)?;a.with(ctx,|o|Ok(NewObject::new(o.reverse(),vec![d])))});
object_binding!(multivector_conjugate,"num_native_multivector_conjugate","crate::kernels::geometric_algebra::Multivector3D::conjugate","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Multivector3D>)->NewObject<Multivector3D>,|ctx:&CallCtx<'_>,a:Handle<Multivector3D>|{let d=object_dimension(&a,ctx)?;a.with(ctx,|o|Ok(NewObject::new(o.conjugate(),vec![d])))});
object_binding!(multivector_norm,"num_native_multivector_norm","crate::kernels::geometric_algebra::Multivector3D::norm","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Multivector3D>)->Dynamic,|ctx:&CallCtx<'_>,a:Handle<Multivector3D>|{let d=object_dimension(&a,ctx)?;a.with(ctx,|o|Ok(Dynamic(api::quantity(o.norm(),d)?)))});
object_binding!(multivector_norm_sq,"num_native_multivector_norm_sq","crate::kernels::geometric_algebra::Multivector3D::norm_sq","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Multivector3D>)->Dynamic,|ctx:&CallCtx<'_>,a:Handle<Multivector3D>|{let d=canonical(&(object_dimension(&a,ctx)?.pow(Exp::Rat(Rational::int(2)))))?;a.with(ctx,|o|Ok(Dynamic(api::quantity(o.norm_sq(),d)?)))});
object_binding!(multivector_inv,"num_native_multivector_inv","crate::kernels::geometric_algebra::Multivector3D::inv","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Multivector3D>)->Maybe<NewObject<Multivector3D>>,|ctx:&CallCtx<'_>,a:Handle<Multivector3D>|{let d=canonical(&(DimExpr::Inverse(Box::new(object_dimension(&a,ctx)?))))?;a.with(ctx,|o|Ok(Maybe(o.inv().map(|v|NewObject::new(v,vec![d])))))});
object_binding!(multivector_dot,"num_native_multivector_dot","crate::kernels::geometric_algebra::Multivector3D::dot","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Multivector3D>,b:Handle<Multivector3D>)->NewObject<Multivector3D>,|ctx:&CallCtx<'_>,a:Handle<Multivector3D>,b:Handle<Multivector3D>|{let d=canonical(&(DimExpr::Mul(vec![object_dimension(&a,ctx)?,object_dimension(&b,ctx)?])))?;a.with(ctx,|o|b.with(ctx,|v|Ok(NewObject::new(o.dot(*v),vec![d]))))});
object_binding!(multivector_wedge,"num_native_multivector_wedge","crate::kernels::geometric_algebra::Multivector3D::wedge","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Multivector3D>,b:Handle<Multivector3D>)->NewObject<Multivector3D>,|ctx:&CallCtx<'_>,a:Handle<Multivector3D>,b:Handle<Multivector3D>|{let d=canonical(&(DimExpr::Mul(vec![object_dimension(&a,ctx)?,object_dimension(&b,ctx)?])))?;a.with(ctx,|o|b.with(ctx,|v|Ok(NewObject::new(o.wedge(*v),vec![d]))))});
object_binding!(affine2_new,"num_native_affine2_new","crate::kernels::fractal_geometry_and_chaos::AffineTransform2D::new","Construct an affine map in explicitly nondimensionalized coordinate space; coefficients represent a,b,c,d,e,f in the native map.",true,(a:f64,b:f64,c:f64,d:f64,e:f64,f:f64)->NewObject<AffineTransform2D>,|_ctx:&CallCtx<'_>,a:f64,b:f64,c:f64,d:f64,e:f64,f:f64|{Ok(NewObject::normalized(AffineTransform2D::new(a,b,c,d,e,f)))});
object_binding!(affine2_apply,"num_native_affine2_apply","crate::kernels::fractal_geometry_and_chaos::AffineTransform2D::apply","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<AffineTransform2D>,point:(f64,f64))->(f64,f64),|ctx:&CallCtx<'_>,a:Handle<AffineTransform2D>,point:(f64,f64)|{a.with(ctx,|o|Ok(o.apply(point)))});
object_binding!(fractal_data_new,"num_native_fractal_data_new","crate::kernels::fractal_geometry_and_chaos::FractalData::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(width:usize,height:usize,iterations:u32)->NewObject<FractalData>,|_ctx:&CallCtx<'_>,width:usize,height:usize,iterations:u32|{limit(width,262144)?;limit(height,262144)?;if width.checked_mul(height).is_none_or(|n|n>262144){return Err(error("fractal image allocation exceeds budget"));}Ok(NewObject::normalized(FractalData::new(width,height,iterations)))});
object_binding!(fractal_data_get,"num_native_fractal_data_get","crate::kernels::fractal_geometry_and_chaos::FractalData::get","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<FractalData>,x:usize,y:usize)->Option<u32>,|ctx:&CallCtx<'_>,a:Handle<FractalData>,x:usize,y:usize|{a.with(ctx,|o|Ok(o.get(x,y)))});
object_binding!(fractal_data_set,"num_native_fractal_data_set","crate::kernels::fractal_geometry_and_chaos::FractalData::set","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<FractalData>,x:usize,y:usize,value:u32)->bool,|ctx:&CallCtx<'_>,a:Handle<FractalData>,x:usize,y:usize,value:u32|{a.with_mut(ctx,|o|{if x>=o.width||y>=o.height{return Err(error("fractal pixel lies outside image bounds"));}o.set(x,y,value);Ok(true)})});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &dd_from_f64::BINDING,
            &dd_from_parts::BINDING,
            &dd_abs::BINDING,
            &dd_sqrt::BINDING,
            &dd_is_finite::BINDING,
            &dd_to_f64::BINDING,
            &cdd_new::BINDING,
            &cdd_norm::BINDING,
            &multivector_new::BINDING,
            &multivector_reverse::BINDING,
            &multivector_conjugate::BINDING,
            &multivector_norm::BINDING,
            &multivector_norm_sq::BINDING,
            &multivector_inv::BINDING,
            &multivector_dot::BINDING,
            &multivector_wedge::BINDING,
            &affine2_new::BINDING,
            &affine2_apply::BINDING,
            &fractal_data_new::BINDING,
            &fractal_data_get::BINDING,
            &fractal_data_set::BINDING,
        ],
        "numerical_vector",
    )
}

impl NativeType for crate::kernels::polynomial::Polynomial {
    const KIND: &'static str = "rssn/native/polynomial/v1";
    const MAX_BYTES: usize = 4096;

    fn validate(&self) -> Result<(), OpError> {
        if self.coeffs.len() > 129 || self.coeffs.iter().any(|v| !v.is_finite()) {
            return Err(error(
                "native polynomial exceeds degree/finite coefficient budget",
            ));
        }
        Ok(())
    }
}
