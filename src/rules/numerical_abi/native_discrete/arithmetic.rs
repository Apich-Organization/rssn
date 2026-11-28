#![allow(clippy::wildcard_imports)] // Checked family method declarations.
use super::super::codecs::IfsPoints;
use super::super::codecs::ifs_work_budget;
use super::*;
use crate::kernels::dd::Cdd;
use crate::kernels::dd::Dd;
use crate::kernels::finite_field::Gf256Multiplier;
use crate::kernels::fractal_geometry_and_chaos::AffineIfs;
use crate::kernels::fractal_geometry_and_chaos::AffineTransform2D;
use crate::kernels::fractal_geometry_and_chaos::FractalData;
use crate::kernels::geometric_algebra::Multivector3D;
use crate::kernels::reuse::PlanCache;
use num_complex::Complex64;
use std::sync::Arc;

fn gf256_multiplier(coefficient: u8) -> Arc<Gf256Multiplier> {
    static MULTIPLIERS: PlanCache<u8, Gf256Multiplier> = PlanCache::new(16);
    MULTIPLIERS.get_or_insert_with(coefficient, || Arc::new(Gf256Multiplier::new(coefficient)))
}

object_binding!(
    gf256_scale,
    "num_finite_field_gf256_scale",
    "crate::kernels::finite_field::Gf256Multiplier::scale_assign",
    "Multiplies every byte by the supplied coefficient in GF(256) with polynomial 0x11d. Inputs and results are exact integers from 0 through 255. Returns the products in input order.",
    true,
    (coefficient: u8, values: Vec<u8>) -> Vec<u8>,
    |_, coefficient, mut values| {
        gf256_multiplier(coefficient).scale_assign(&mut values);
        Ok(values)
    }
);

object_binding!(
    gf256_mul_add,
    "num_finite_field_gf256_mul_add",
    "crate::kernels::finite_field::Gf256Multiplier::mul_add_assign",
    "Returns destination XOR coefficient times source over GF(256) with polynomial 0x11d. Supply equal-length byte vectors and an exact byte coefficient.",
    true,
    (coefficient: u8, source: Vec<u8>, destination: Vec<u8>) -> Vec<u8>,
    |_, coefficient, source, mut destination| {
        gf256_multiplier(coefficient).mul_add_assign(&source, &mut destination).map_err(error)?;
        Ok(destination)
    }
);

struct WindingVertices(Vec<Complex64>);

impl NativeArg for WindingVertices {
    fn ty() -> TypeExpr {
        api::tensor_type(1, TypeExpr::Complex.with_dim(DimExpr::Var("D".into())))
    }

    fn decode(
        value: &Value,
        index: usize,
        _: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        api::dimension(value)?;
        <Vec<Complex64> as Arg>::decode(value, index).map(Self)
    }
}

object_binding!(
    complex_winding_number,
    "num_complex_winding_number",
    "crate::kernels::complex::winding_number",
    "Returns the winding number about zero of the polygon formed by a rank-one tensor of complex vertices with a common physical dimension. The last vertex connects to the first. Counterclockwise winding is positive. Empty paths and paths through zero return a domain error.",
    false,
    (vertices: WindingVertices) -> i64,
    |_, vertices| {
        crate::kernels::complex::winding_number(&vertices.0)
            .ok_or_else(|| error("complex polygon has no defined winding number"))
    }
);

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
impl NativeType for AffineIfs {
    const KIND: &'static str = "rssn/native/discrete/affine-ifs/v1";
    const MAX_BYTES: usize = 262_144;

    fn validate(&self) -> Result<(), OpError> {
        Ok(())
    }

    fn retained_bytes(&self) -> usize {
        std::mem::size_of::<Self>().saturating_add(self.heap_bytes())
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
object_binding!(
    dd_hypot,
    "num_native_dd_hypot",
    "crate::kernels::dd::Dd::hypot",
    "Computes the two-coordinate Euclidean norm at double-double precision. Both session-owned inputs must have the same physical dimension, which the retained result preserves.",
    false,
    (x: Handle<Dd>, y: Handle<Dd>) -> NewObject<Dd>,
    |ctx, x, y| {
        let dimension = object_dimension(&x, ctx)?;
        api::same_dimension(&dimension, &object_dimension(&y, ctx)?)?;
        x.with(ctx, |a| {
            y.with(ctx, |b| Ok(NewObject::new(a.hypot(*b), vec![dimension])))
        })
    }
);
struct DdComponents {
    value: Dd,
    dimension: DimExpr,
}

impl NativeOutput for DdComponents {
    fn ty() -> TypeExpr {
        TypeExpr::Record(vec![
            ("item0".into(), api::dynamic_type()),
            ("item1".into(), api::dynamic_type()),
        ])
    }

    fn encode(
        self,
        _: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        let (high, low) = self.value.parts();
        Ok(api::record(vec![
            ("item0", api::quantity(high, self.dimension.clone())?),
            ("item1", api::quantity(low, self.dimension)?),
        ]))
    }
}

object_binding!(
    dd_parts,
    "num_native_dd_parts",
    "crate::kernels::dd::Dd::parts",
    "Returns the high and low components as item0 and item1. Both have the physical dimension of the double-double input, and the low component is kept separately from the high one.",
    false,
    (x: Handle<Dd>) -> DdComponents,
    |ctx, x| {
        let dimension = object_dimension(&x, ctx)?;
        x.with(ctx, |value| Ok(DdComponents { value: *value, dimension }))
    }
);
object_binding!(dd_is_finite,"num_native_dd_is_finite","crate::kernels::dd::Dd::is_finite","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Dd>)->bool,|ctx:&CallCtx<'_>,a:Handle<Dd>|{a.with(ctx,|o|Ok(o.is_finite()))});
object_binding!(dd_to_f64,"num_native_dd_to_f64","crate::kernels::dd::Dd::to_f64","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Dd>)->Dynamic,|ctx:&CallCtx<'_>,a:Handle<Dd>|{let d=object_dimension(&a,ctx)?;a.with(ctx,|o|Ok(Dynamic(api::quantity(o.to_f64(),d)?)))});
object_binding!(cdd_new,"num_native_cdd_new","crate::kernels::dd::Cdd::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(re:Input,im:Input)->NewObject<Cdd>,|_ctx:&CallCtx<'_>,re:Input,im:Input|{let(r,d)=scalar(&re)?;let(i,e)=scalar(&im)?;api::same_dimension(&d,&e)?;Ok(NewObject::new(Cdd::new(r,i),vec![d]))});
object_binding!(cdd_norm,"num_native_cdd_norm","crate::kernels::dd::Cdd::norm","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Cdd>)->Dynamic,|ctx:&CallCtx<'_>,a:Handle<Cdd>|{let d=object_dimension(&a,ctx)?;a.with(ctx,|o|Ok(Dynamic(api::quantity(o.norm(),d)?)))});
object_binding!(
    cdd_norm_dd,
    "num_native_cdd_norm_dd",
    "crate::kernels::dd::Cdd::norm_dd",
    "Computes the modulus of a complex double-double value and keeps both components of the real double-double result. The result preserves the input's physical dimension.",
    false,
    (z: Handle<Cdd>) -> NewObject<Dd>,
    |ctx, z| {
        let dimension = object_dimension(&z, ctx)?;
        z.with(ctx, |value| Ok(NewObject::new(value.norm_dd(), vec![dimension])))
    }
);
object_binding!(
    cdd_recip,
    "num_native_cdd_recip",
    "crate::kernels::dd::Cdd::recip",
    "Computes the reciprocal of a complex double-double value. The retained result has the inverse physical dimension. Zero inputs and nonfinite results are rejected.",
    false,
    (z: Handle<Cdd>) -> NewObject<Cdd>,
    |ctx, z| {
        let dimension = canonical(&DimExpr::Inverse(Box::new(object_dimension(&z, ctx)?)))?;
        z.with(ctx, |value| Ok(NewObject::new(value.recip(), vec![dimension])))
    }
);
object_binding!(multivector_new,"num_native_multivector_new","crate::kernels::geometric_algebra::Multivector3D::new","Construct a homogeneous-coefficient Cl(3,0) multivector; the eight basis blades are dimensionless and all coefficients share physical units.",false,(s:Input,v1:Input,v2:Input,v3:Input,b12:Input,b23:Input,b31:Input,pss:Input)->NewObject<Multivector3D>,|_ctx:&CallCtx<'_>,s:Input,v1:Input,v2:Input,v3:Input,b12:Input,b23:Input,b31:Input,pss:Input|{let values=[scalar(&s)?,scalar(&v1)?,scalar(&v2)?,scalar(&v3)?,scalar(&b12)?,scalar(&b23)?,scalar(&b31)?,scalar(&pss)?];let d=values[0].1.clone();for(_,e)in &values{api::same_dimension(&d,e)?;}Ok(NewObject::new(Multivector3D::new(values[0].0,values[1].0,values[2].0,values[3].0,values[4].0,values[5].0,values[6].0,values[7].0),vec![d]))});
object_binding!(multivector_reverse,"num_native_multivector_reverse","crate::kernels::geometric_algebra::Multivector3D::reverse","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Multivector3D>)->NewObject<Multivector3D>,|ctx:&CallCtx<'_>,a:Handle<Multivector3D>|{let d=object_dimension(&a,ctx)?;a.with(ctx,|o|Ok(NewObject::new(o.reverse(),vec![d])))});
object_binding!(multivector_conjugate,"num_native_multivector_conjugate","crate::kernels::geometric_algebra::Multivector3D::conjugate","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Multivector3D>)->NewObject<Multivector3D>,|ctx:&CallCtx<'_>,a:Handle<Multivector3D>|{let d=object_dimension(&a,ctx)?;a.with(ctx,|o|Ok(NewObject::new(o.conjugate(),vec![d])))});
object_binding!(
    multivector_mul_central,
    "num_native_multivector_mul_central",
    "crate::kernels::geometric_algebra::Multivector3D::mul_central",
    "Multiplies a Cl(3,0) multivector by scalar + pseudoscalar times e123. The two factor coefficients share one physical dimension. The product coefficients have that dimension times the input coefficient dimension.",
    false,
    (a: Handle<Multivector3D>, central_scalar: Input, central_pseudoscalar: Input) -> NewObject<Multivector3D>,
    |ctx, a, central_scalar, central_pseudoscalar| {
        let (real, factor_dimension) = scalar(&central_scalar)?;
        let (imaginary, second_dimension) = scalar(&central_pseudoscalar)?;
        api::same_dimension(&factor_dimension, &second_dimension)?;
        let dimension = canonical(&DimExpr::Mul(vec![object_dimension(&a, ctx)?, factor_dimension]))?;
        a.with(ctx, |value| Ok(NewObject::new(value.mul_central(real, imaginary), vec![dimension])))
    }
);
object_binding!(multivector_norm,"num_native_multivector_norm","crate::kernels::geometric_algebra::Multivector3D::norm","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Multivector3D>)->Dynamic,|ctx:&CallCtx<'_>,a:Handle<Multivector3D>|{let d=object_dimension(&a,ctx)?;a.with(ctx,|o|Ok(Dynamic(api::quantity(o.norm(),d)?)))});
object_binding!(multivector_norm_sq,"num_native_multivector_norm_sq","crate::kernels::geometric_algebra::Multivector3D::norm_sq","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Multivector3D>)->Dynamic,|ctx:&CallCtx<'_>,a:Handle<Multivector3D>|{let d=canonical(&(object_dimension(&a,ctx)?.pow(Exp::Rat(Rational::int(2)))))?;a.with(ctx,|o|Ok(Dynamic(api::quantity(o.norm_sq(),d)?)))});
object_binding!(multivector_inv,"num_native_multivector_inv","crate::kernels::geometric_algebra::Multivector3D::inv","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Multivector3D>)->Maybe<NewObject<Multivector3D>>,|ctx:&CallCtx<'_>,a:Handle<Multivector3D>|{let d=canonical(&(DimExpr::Inverse(Box::new(object_dimension(&a,ctx)?))))?;a.with(ctx,|o|Ok(Maybe(o.inv().map(|v|NewObject::new(v,vec![d])))))});
object_binding!(multivector_dot,"num_native_multivector_dot","crate::kernels::geometric_algebra::Multivector3D::dot","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Multivector3D>,b:Handle<Multivector3D>)->NewObject<Multivector3D>,|ctx:&CallCtx<'_>,a:Handle<Multivector3D>,b:Handle<Multivector3D>|{let d=canonical(&(DimExpr::Mul(vec![object_dimension(&a,ctx)?,object_dimension(&b,ctx)?])))?;a.with(ctx,|o|b.with(ctx,|v|Ok(NewObject::new(o.dot(*v),vec![d]))))});
object_binding!(multivector_wedge,"num_native_multivector_wedge","crate::kernels::geometric_algebra::Multivector3D::wedge","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Multivector3D>,b:Handle<Multivector3D>)->NewObject<Multivector3D>,|ctx:&CallCtx<'_>,a:Handle<Multivector3D>,b:Handle<Multivector3D>|{let d=canonical(&(DimExpr::Mul(vec![object_dimension(&a,ctx)?,object_dimension(&b,ctx)?])))?;a.with(ctx,|o|b.with(ctx,|v|Ok(NewObject::new(o.wedge(*v),vec![d]))))});
object_binding!(affine2_new,"num_native_affine2_new","crate::kernels::fractal_geometry_and_chaos::AffineTransform2D::new","Construct an affine map in explicitly nondimensionalized coordinate space; coefficients represent a,b,c,d,e,f in the native map.",true,(a:f64,b:f64,c:f64,d:f64,e:f64,f:f64)->NewObject<AffineTransform2D>,|_ctx:&CallCtx<'_>,a:f64,b:f64,c:f64,d:f64,e:f64,f:f64|{Ok(NewObject::normalized(AffineTransform2D::new(a,b,c,d,e,f)))});
object_binding!(affine2_apply,"num_native_affine2_apply","crate::kernels::fractal_geometry_and_chaos::AffineTransform2D::apply","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<AffineTransform2D>,point:(f64,f64))->(f64,f64),|ctx:&CallCtx<'_>,a:Handle<AffineTransform2D>,point:(f64,f64)|{a.with(ctx,|o|Ok(o.apply(point)))});
object_binding!(
    affine_ifs_new,
    "num_native_affine_ifs_new",
    "crate::kernels::fractal_geometry_and_chaos::AffineIfs::new",
    "Prepares a session-owned affine iterated function system from up to 4096 dimensionless rows (a,b,c,d,e,f) and matching finite nonnegative relative weights. A positive weight is required. Reuse the plan across independent starting points and seeds.",
    true,
    (maps: Vec<AffineTransform2D>, weights: Vec<f64>) -> NewObject<AffineIfs>,
    |_, maps: Vec<AffineTransform2D>, weights: Vec<f64>| {
        limit(maps.len(), 4096)?;
        Ok(NewObject::normalized(AffineIfs::new(&maps, &weights).map_err(error)?))
    }
);
object_binding!(
    affine_ifs_generate,
    "num_native_affine_ifs_generate",
    "crate::kernels::fractal_geometry_and_chaos::AffineIfs::generate",
    "Generates a dimensionless two-column point tensor after the requested burn-in, using an exact unsigned 64-bit seed. The starting point is excluded. A zero point count returns shape [0, 2]. Nonfinite trajectories and excessive work return errors.",
    true,
    (plan: Handle<AffineIfs>, start: (f64, f64), points: usize, skip: usize, seed: u64) -> IfsPoints,
    |ctx, plan, start, points, skip, seed| {
        ifs_work_budget(points, skip)?;
        plan.with(ctx, |prepared| prepared.generate(start, points, skip, seed).map(IfsPoints).map_err(error))
    }
);
object_binding!(fractal_data_new,"num_native_fractal_data_new","crate::kernels::fractal_geometry_and_chaos::FractalData::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(width:usize,height:usize,iterations:u32)->NewObject<FractalData>,|_ctx:&CallCtx<'_>,width:usize,height:usize,iterations:u32|{limit(width,262144)?;limit(height,262144)?;if width.checked_mul(height).is_none_or(|n|n>262144){return Err(error("fractal image allocation exceeds budget"));}Ok(NewObject::normalized(FractalData::new(width,height,iterations)))});
object_binding!(fractal_data_get,"num_native_fractal_data_get","crate::kernels::fractal_geometry_and_chaos::FractalData::get","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<FractalData>,x:usize,y:usize)->Option<u32>,|ctx:&CallCtx<'_>,a:Handle<FractalData>,x:usize,y:usize|{a.with(ctx,|o|Ok(o.get(x,y)))});
object_binding!(fractal_data_set,"num_native_fractal_data_set","crate::kernels::fractal_geometry_and_chaos::FractalData::set","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<FractalData>,x:usize,y:usize,value:u32)->bool,|ctx:&CallCtx<'_>,a:Handle<FractalData>,x:usize,y:usize,value:u32|{a.with_mut(ctx,|o|{if x>=o.width||y>=o.height{return Err(error("fractal pixel lies outside image bounds"));}o.set(x,y,value);Ok(true)})});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &gf256_scale::BINDING,
            &gf256_mul_add::BINDING,
            &complex_winding_number::BINDING,
            &dd_from_f64::BINDING,
            &dd_from_parts::BINDING,
            &dd_abs::BINDING,
            &dd_sqrt::BINDING,
            &dd_hypot::BINDING,
            &dd_parts::BINDING,
            &dd_is_finite::BINDING,
            &dd_to_f64::BINDING,
            &cdd_new::BINDING,
            &cdd_norm::BINDING,
            &cdd_norm_dd::BINDING,
            &cdd_recip::BINDING,
            &multivector_new::BINDING,
            &multivector_reverse::BINDING,
            &multivector_conjugate::BINDING,
            &multivector_norm::BINDING,
            &multivector_mul_central::BINDING,
            &multivector_norm_sq::BINDING,
            &multivector_inv::BINDING,
            &multivector_dot::BINDING,
            &multivector_wedge::BINDING,
            &affine2_new::BINDING,
            &affine2_apply::BINDING,
            &affine_ifs_new::BINDING,
            &affine_ifs_generate::BINDING,
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
