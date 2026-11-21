#![allow(clippy::wildcard_imports)] // Native graphics with physical coordinates, direction contracts and explicitly normalized frames.
use super::*;
use crate::kernels::computer_graphics::Color;
use crate::kernels::computer_graphics::InnerProduct;
use crate::kernels::computer_graphics::Intersection;
use crate::kernels::computer_graphics::Plane;
use crate::kernels::computer_graphics::Point2D;
use crate::kernels::computer_graphics::Point3D;
use crate::kernels::computer_graphics::PreparedTriangle;
use crate::kernels::computer_graphics::Quaternion;
use crate::kernels::computer_graphics::Ray;
use crate::kernels::computer_graphics::Sphere;
use crate::kernels::computer_graphics::Vector2D;
use crate::kernels::computer_graphics::Vector3D;
use crate::kernels::computer_graphics::{
    self as g,
};
use crate::kernels::matrix::Matrix;
impl NativeType for Point2D {
    const KIND: &'static str = "rssn/native/discrete/graphics-point2d/v1";
    const MAX_BYTES: usize = 128;

    fn validate(&self) -> Result<(), OpError> {
        if [self.x, self.y].iter().any(|v| !v.is_finite()) {
            return Err(error("graphics components must be finite"));
        }
        Ok(())
    }
}
impl NativeType for Point3D {
    const KIND: &'static str = "rssn/native/discrete/graphics-point3d/v1";
    const MAX_BYTES: usize = 128;

    fn validate(&self) -> Result<(), OpError> {
        if [self.x, self.y, self.z].iter().any(|v| !v.is_finite()) {
            return Err(error("graphics components must be finite"));
        }
        Ok(())
    }
}
impl NativeType for Vector2D {
    const KIND: &'static str = "rssn/native/discrete/graphics-vector2d/v1";
    const MAX_BYTES: usize = 128;

    fn validate(&self) -> Result<(), OpError> {
        if [self.x, self.y].iter().any(|v| !v.is_finite()) {
            return Err(error("graphics components must be finite"));
        }
        Ok(())
    }
}
impl NativeType for Vector3D {
    const KIND: &'static str = "rssn/native/discrete/graphics-vector3d/v1";
    const MAX_BYTES: usize = 128;

    fn validate(&self) -> Result<(), OpError> {
        if [self.x, self.y, self.z].iter().any(|v| !v.is_finite()) {
            return Err(error("graphics components must be finite"));
        }
        Ok(())
    }
}
impl NativeType for Color {
    const KIND: &'static str = "rssn/native/discrete/graphics-color/v1";
    const MAX_BYTES: usize = 128;

    fn validate(&self) -> Result<(), OpError> {
        if [self.r, self.g, self.b, self.a]
            .iter()
            .any(|v| !v.is_finite())
        {
            return Err(error("graphics components must be finite"));
        }
        Ok(())
    }
}
impl NativeType for Quaternion {
    const KIND: &'static str = "rssn/native/discrete/graphics-quaternion/v1";
    const MAX_BYTES: usize = 128;

    fn validate(&self) -> Result<(), OpError> {
        if [self.w, self.x, self.y, self.z]
            .iter()
            .any(|v| !v.is_finite())
        {
            return Err(error("graphics components must be finite"));
        }
        Ok(())
    }
}
struct Frame {
    matrix: Matrix<f64>,
    coordinate: f64,
    output: f64,
    framed: bool,
    projective: bool,
}
struct GfxRay {
    ray: Ray,
    reference: f64,
}
struct GfxPreparedTriangle {
    triangle: PreparedTriangle,
    reference: f64,
}

const MAX_TRIANGLE_RAYS: usize = 4096;

struct TriangleOrigins(Vec<Handle<Point3D>>);

impl NativeArg for TriangleOrigins {
    fn ty() -> TypeExpr {
        api::tensor_type(1, TypeExpr::Opaque(Point3D::KIND.into()))
    }

    fn decode(
        value: &Value,
        index: usize,
        ctx: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        let Value::Tensor(tensor) = value else {
            return Err(error("expected a rank-one tensor of ray origins"));
        };
        if tensor.shape != [tensor.data.len()] || tensor.data.len() > MAX_TRIANGLE_RAYS {
            return Err(error(
                "ray origins must be rank one with at most 4096 points",
            ));
        }
        let origins = tensor
            .data
            .iter()
            .map(|v| Handle::<Point3D>::decode(v, index, ctx))
            .collect::<Result<_, _>>()?;
        Ok(Self(origins))
    }
}

impl NativeType for GfxPreparedTriangle {
    const KIND: &'static str = "rssn/native/discrete/graphics-prepared-triangle/v1";
    const MAX_BYTES: usize = 512;

    fn validate(&self) -> Result<(), OpError> {
        if self.reference <= 0.0 || !self.reference.is_finite() {
            return Err(error(
                "triangle coordinate reference must be positive and finite",
            ));
        }
        Ok(())
    }
}
impl NativeType for Frame {
    const KIND: &'static str = "rssn/native/discrete/graphics-frame/v1";
    const MAX_BYTES: usize = 512;

    fn validate(&self) -> Result<(), OpError> {
        if self.matrix.rows() != 4
            || self.matrix.cols() != 4
            || self.matrix.data().iter().any(|v| !v.is_finite())
            || self.coordinate <= 0.0
            || self.output <= 0.0
            || !self.coordinate.is_finite()
            || !self.output.is_finite()
        {
            return Err(error("invalid finite4x4normalized graphics frame"));
        }
        Ok(())
    }
}
impl NativeType for GfxRay {
    const KIND: &'static str = "rssn/native/discrete/graphics-ray/v1";
    const MAX_BYTES: usize = 128;

    fn validate(&self) -> Result<(), OpError> {
        self.ray.origin.validate()?;
        unit(&self.ray.direction)?;
        if self.reference <= 0.0 || !self.reference.is_finite() {
            return Err(error("ray coordinate reference must be positive"));
        }
        Ok(())
    }
}
impl NativeType for Sphere {
    const KIND: &'static str = "rssn/native/discrete/graphics-sphere/v1";
    const MAX_BYTES: usize = 128;

    fn validate(&self) -> Result<(), OpError> {
        self.center.validate()?;
        if self.radius <= 0.0 || !self.radius.is_finite() {
            return Err(error("sphere radius must be positive finite"));
        }
        Ok(())
    }
}
impl NativeType for Plane {
    const KIND: &'static str = "rssn/native/discrete/graphics-plane/v1";
    const MAX_BYTES: usize = 128;

    fn validate(&self) -> Result<(), OpError> {
        self.point.validate()?;
        unit(&self.normal)
    }
}
fn unit(v: &Vector3D) -> Result<(), OpError> {
    v.validate()?;
    if (v.magnitude() - 1.0).abs() > 1e-10 {
        return Err(error(
            "native direction/normal must be a dimensionless unit vector",
        ));
    }
    Ok(())
}
fn normalized_direction(
    h: &Handle<Vector3D>,
    ctx: &CallCtx<'_>,
) -> Result<(), OpError> {
    api::same_dimension(&object_dimension(h, ctx)?, &DimExpr::Dimless)?;
    h.with(ctx, unit)
}
const fn frame(matrix: Matrix<f64>) -> NewObject<Frame> {
    NewObject::normalized(Frame {
        matrix,
        coordinate: 1.0,
        output: 1.0,
        framed: false,
        projective: false,
    })
}
fn framed(
    matrix: Matrix<f64>,
    r: f64,
    out: f64,
    input: DimExpr,
    output: DimExpr,
    projective: bool,
) -> NewObject<Frame> {
    NewObject::new(
        Frame {
            matrix,
            coordinate: r,
            output: out,
            framed: true,
            projective,
        },
        vec![input, output],
    )
}
fn scaled_point(
    p: Point3D,
    scale: f64,
) -> Point3D {
    Point3D::new(p.x * scale, p.y * scale, p.z * scale)
}
fn scaled_vector(
    p: Vector3D,
    scale: f64,
) -> Vector3D {
    Vector3D::new(p.x * scale, p.y * scale, p.z * scale)
}
fn hit(
    v: Option<Intersection>,
    r: f64,
    d: DimExpr,
    ctx: &CallCtx<'_>,
) -> Result<Maybe<Dynamic>, OpError> {
    Ok(Maybe(
        v.map(|h| {
            let point = NewObject::new(scaled_point(h.point, r), vec![d.clone()]).encode(ctx)?;
            let normal = NewObject::normalized(h.normal).encode(ctx)?;
            Ok::<_, OpError>(Dynamic(api::record(vec![
                ("t", api::quantity(h.t * r, d)?),
                ("point", point),
                ("normal", normal),
            ])))
        })
        .transpose()?,
    ))
}

fn rescale_triangle_point(
    p: Point3D,
    numerator: &[f64],
    denominator: &[f64],
) -> Result<Point3D, OpError> {
    let scale = |value| api::scale_without_underflow(value, numerator, denominator);
    Ok(Point3D::new(scale(p.x)?, scale(p.y)?, scale(p.z)?))
}

fn triangle_point(
    point: &Handle<Point3D>,
    dimension: &DimExpr,
    reference: f64,
    ctx: &CallCtx<'_>,
) -> Result<Point3D, OpError> {
    api::same_dimension(dimension, &object_dimension(point, ctx)?)?;
    point.with(ctx, |p| rescale_triangle_point(*p, &[], &[reference]))
}

fn prepared_triangle(
    vertices: [&Handle<Point3D>; 3],
    coordinate_reference: &Input,
    ctx: &CallCtx<'_>,
) -> Result<Maybe<NewObject<GfxPreparedTriangle>>, OpError> {
    let (reference, dimension) = api::reference(&coordinate_reference.0)?;
    let [a, b, c] = vertices;
    let a = triangle_point(a, &dimension, reference, ctx)?;
    let b = triangle_point(b, &dimension, reference, ctx)?;
    let c = triangle_point(c, &dimension, reference, ctx)?;
    Ok(Maybe(PreparedTriangle::new(&a, &b, &c).map(|triangle| {
        NewObject::new(GfxPreparedTriangle { triangle, reference }, vec![dimension])
    })))
}

fn triangle_ray_bounds(
    min_t: &Input,
    max_t: Option<&api::Input>,
    dimension: &DimExpr,
    reference: f64,
) -> Result<(f64, f64), OpError> {
    let read = |value: &Value| {
        api::same_dimension(dimension, &api::dimension(value)?)?;
        numeric(value, 0)
    };
    let min = read(&min_t.0)?;
    let max = max_t.map(|input| read(&input.0)).transpose()?;
    if min < 0.0 || max.is_some_and(|max| max < min) {
        return Err(error("ray bounds require 0 <= min_t <= max_t"));
    }
    let scale = |value| api::scale_without_underflow(value, &[], &[reference]);
    Ok((
        scale(min)?,
        max.map(scale).transpose()?.unwrap_or(f64::INFINITY),
    ))
}

fn physical_triangle_hit(
    h: Intersection,
    reference: f64,
) -> Result<Intersection, OpError> {
    unit(&h.normal)?;
    Ok(Intersection {
        t: api::scale_without_underflow(h.t, &[reference], &[])?,
        point: rescale_triangle_point(h.point, &[reference], &[])?,
        normal: h.normal,
    })
}

fn prepared_triangle_hit(
    triangle: &Handle<GfxPreparedTriangle>,
    ray: &Handle<GfxRay>,
    min_t: &Input,
    max_t: Option<&api::Input>,
    ctx: &CallCtx<'_>,
) -> Result<Maybe<Dynamic>, OpError> {
    let dimension = same(triangle, ray, ctx)?;
    let physical = triangle.with(ctx, |t| {
        ray.with(ctx, |r| {
            let (min, max) = triangle_ray_bounds(min_t, max_t, &dimension, t.reference)?;
            let origin = rescale_triangle_point(r.ray.origin, &[r.reference], &[t.reference])?;
            let ray = Ray::new(origin, r.ray.direction);
            t.triangle
                .intersect(&ray, min, max)
                .map(|h| physical_triangle_hit(h, t.reference))
                .transpose()
        })
    })?;
    hit(physical, 1.0, dimension, ctx)
}

fn batched_triangle_record(
    hits: Vec<Option<Intersection>>,
    reference: f64,
    dimension: &DimExpr,
) -> Result<Dynamic, OpError> {
    let count = hits.len();
    let mut mask = Vec::with_capacity(count);
    let mut parameters = Vec::with_capacity(count);
    let mut points = Vec::with_capacity(3 * count);
    let mut normals = Vec::with_capacity(3 * count);
    let zero_hit = Intersection {
        t: 0.0,
        point: Point3D::new(0.0, 0.0, 0.0),
        normal: Vector3D::new(0.0, 0.0, 0.0),
    };
    for hit in hits {
        mask.push(Value::Bool(hit.is_some()));
        let hit = hit
            .map(|h| physical_triangle_hit(h, reference))
            .transpose()?
            .unwrap_or(zero_hit);
        parameters.push(api::quantity(hit.t, dimension.clone())?);
        for coordinate in [hit.point.x, hit.point.y, hit.point.z] {
            points.push(api::quantity(coordinate, dimension.clone())?);
        }
        normals.extend([hit.normal.x, hit.normal.y, hit.normal.z].map(Value::Float));
    }
    let measured = TypeExpr::Float.with_dim(dimension.clone());
    let tensor = |shape, elem, data| Value::Tensor(TensorValue { shape, data, elem });
    Ok(Dynamic(api::record(vec![
        ("hit_mask", tensor(vec![count], TypeExpr::Bool, mask)),
        ("t", tensor(vec![count], measured.clone(), parameters)),
        ("points", tensor(vec![count, 3], measured, points)),
        ("normals", tensor(vec![count, 3], TypeExpr::Float, normals)),
    ])))
}

fn prepared_triangle_parallel_hits(
    triangle: &Handle<GfxPreparedTriangle>,
    origins: &TriangleOrigins,
    direction: &Handle<Vector3D>,
    min_t: &Input,
    max_t: Option<&api::Input>,
    ctx: &CallCtx<'_>,
) -> Result<Dynamic, OpError> {
    normalized_direction(direction, ctx)?;
    let direction = direction.with(ctx, |v| Ok(*v))?;
    let dimension = object_dimension(triangle, ctx)?;
    triangle.with(ctx, |t| {
        let (min, max) = triangle_ray_bounds(min_t, max_t, &dimension, t.reference)?;
        let inputs = origins
            .0
            .iter()
            .map(|p| triangle_point(p, &dimension, t.reference, ctx))
            .collect::<Result<Vec<_>, _>>()?;
        let hits = t.triangle.intersect_parallel(&inputs, direction, min, max);
        batched_triangle_record(hits, t.reference, &dimension)
    })
}
object_binding!(graphics_point2_new,"num_native_graphics_point2_new","crate::kernels::computer_graphics::Point2D::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(x:Input,y:Input)->NewObject<Point2D>,|_ctx:&CallCtx<'_>,x:Input,y:Input|{let(x,d)=scalar(&x)?;let(y,e)=scalar(&y)?;api::same_dimension(&d,&e)?;Ok(NewObject::new(Point2D::new(x,y),vec![d]))});
object_binding!(graphics_point2_components,"num_native_graphics_point2_components","","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Point2D>)->Dynamic,|ctx:&CallCtx<'_>,a:Handle<Point2D>|{let d=object_dimension(&a,ctx)?;a.with(ctx,|v|quantities(vec![v.x,v.y],&(d)))});
object_binding!(graphics_point3_new,"num_native_graphics_point3_new","crate::kernels::computer_graphics::Point3D::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(x:Input,y:Input,z:Input)->NewObject<Point3D>,|_ctx:&CallCtx<'_>,x:Input,y:Input,z:Input|{let(x,d)=scalar(&x)?;let(y,e)=scalar(&y)?;api::same_dimension(&d,&e)?;let(z,e)=scalar(&z)?;api::same_dimension(&d,&e)?;Ok(NewObject::new(Point3D::new(x,y,z),vec![d]))});
object_binding!(graphics_point3_components,"num_native_graphics_point3_components","","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Point3D>)->Dynamic,|ctx:&CallCtx<'_>,a:Handle<Point3D>|{let d=object_dimension(&a,ctx)?;a.with(ctx,|v|quantities(vec![v.x,v.y,v.z],&(d)))});
object_binding!(graphics_vector2_new,"num_native_graphics_vector2_new","crate::kernels::computer_graphics::Vector2D::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(x:Input,y:Input)->NewObject<Vector2D>,|_ctx:&CallCtx<'_>,x:Input,y:Input|{let(x,d)=scalar(&x)?;let(y,e)=scalar(&y)?;api::same_dimension(&d,&e)?;Ok(NewObject::new(Vector2D::new(x,y),vec![d]))});
object_binding!(graphics_vector2_components,"num_native_graphics_vector2_components","","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Vector2D>)->Dynamic,|ctx:&CallCtx<'_>,a:Handle<Vector2D>|{let d=object_dimension(&a,ctx)?;a.with(ctx,|v|quantities(vec![v.x,v.y],&(d)))});
object_binding!(graphics_vector3_new,"num_native_graphics_vector3_new","crate::kernels::computer_graphics::Vector3D::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(x:Input,y:Input,z:Input)->NewObject<Vector3D>,|_ctx:&CallCtx<'_>,x:Input,y:Input,z:Input|{let(x,d)=scalar(&x)?;let(y,e)=scalar(&y)?;api::same_dimension(&d,&e)?;let(z,e)=scalar(&z)?;api::same_dimension(&d,&e)?;Ok(NewObject::new(Vector3D::new(x,y,z),vec![d]))});
object_binding!(graphics_vector3_components,"num_native_graphics_vector3_components","","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Vector3D>)->Dynamic,|ctx:&CallCtx<'_>,a:Handle<Vector3D>|{let d=object_dimension(&a,ctx)?;a.with(ctx,|v|quantities(vec![v.x,v.y,v.z],&(d)))});
object_binding!(graphics_point2_distance_to,"num_native_graphics_point2_distance_to","crate::kernels::computer_graphics::Point2D::distance_to","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Point2D>,b:Handle<Point2D>)->Dynamic,|ctx:&CallCtx<'_>,a:Handle<Point2D>,b:Handle<Point2D>|{let d=same(&a,&b,ctx)?;a.with(ctx,|a|b.with(ctx,|b|extended_quantity(a.distance_to(b),d)))});
object_binding!(graphics_point3_distance_to,"num_native_graphics_point3_distance_to","crate::kernels::computer_graphics::Point3D::distance_to","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Point3D>,b:Handle<Point3D>)->Dynamic,|ctx:&CallCtx<'_>,a:Handle<Point3D>,b:Handle<Point3D>|{let d=same(&a,&b,ctx)?;a.with(ctx,|a|b.with(ctx,|b|extended_quantity(a.distance_to(b),d)))});
object_binding!(graphics_point3_to_vector,"num_native_graphics_point3_to_vector","crate::kernels::computer_graphics::Point3D::to_vector","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Point3D>)->NewObject<Vector3D>,|ctx:&CallCtx<'_>,a:Handle<Point3D>|{let d=object_dimension(&a,ctx)?;a.with(ctx,|a|Ok(NewObject::new(a.to_vector(),vec![d])))});
object_binding!(graphics_vector2_magnitude,"num_native_graphics_vector2_magnitude","crate::kernels::computer_graphics::Vector2D::magnitude","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Vector2D>)->Dynamic,|ctx:&CallCtx<'_>,a:Handle<Vector2D>|{let d=object_dimension(&a,ctx)?;a.with(ctx,|a|extended_quantity(a.magnitude(),d))});
object_binding!(graphics_vector2_normalize,"num_native_graphics_vector2_normalize","crate::kernels::computer_graphics::Vector2D::normalize","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Vector2D>)->NewObject<Vector2D>,|ctx:&CallCtx<'_>,a:Handle<Vector2D>|{a.with(ctx,|a|Ok(NewObject::normalized(a.normalize())))});
object_binding!(graphics_vector3_magnitude,"num_native_graphics_vector3_magnitude","crate::kernels::computer_graphics::Vector3D::magnitude","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Vector3D>)->Dynamic,|ctx:&CallCtx<'_>,a:Handle<Vector3D>|{let d=object_dimension(&a,ctx)?;a.with(ctx,|a|extended_quantity(a.magnitude(),d))});
object_binding!(graphics_vector3_normalize,"num_native_graphics_vector3_normalize","crate::kernels::computer_graphics::Vector3D::normalize","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Vector3D>)->NewObject<Vector3D>,|ctx:&CallCtx<'_>,a:Handle<Vector3D>|{a.with(ctx,|a|Ok(NewObject::normalized(a.normalize())))});
object_binding!(graphics_vector3_magnitude_squared,"num_native_graphics_vector3_magnitude_squared","crate::kernels::computer_graphics::Vector3D::magnitude_squared","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Vector3D>)->Dynamic,|ctx:&CallCtx<'_>,a:Handle<Vector3D>|{let d=canonical(&(object_dimension(&a,ctx)?.pow(Exp::Rat(Rational::int(2)))))?;a.with(ctx,|a|extended_quantity(a.magnitude_squared(),d))});
object_binding!(graphics_vector2_rotate,"num_native_graphics_vector2_rotate","crate::kernels::computer_graphics::Vector2D::rotate","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Vector2D>,angle:f64)->NewObject<Vector2D>,|ctx:&CallCtx<'_>,a:Handle<Vector2D>,angle:f64|{let d=object_dimension(&a,ctx)?;a.with(ctx,|a|Ok(NewObject::new(a.rotate(angle),vec![d])))});
object_binding!(graphics_vector2_perpendicular,"num_native_graphics_vector2_perpendicular","crate::kernels::computer_graphics::Vector2D::perpendicular","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Vector2D>)->NewObject<Vector2D>,|ctx:&CallCtx<'_>,a:Handle<Vector2D>|{let d=object_dimension(&a,ctx)?;a.with(ctx,|a|Ok(NewObject::new(a.perpendicular(),vec![d])))});
object_binding!(graphics_color_new,"num_native_graphics_color_new","crate::kernels::computer_graphics::Color::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(r:f64,g:f64,b:f64,a:f64)->NewObject<Color>,|_ctx:&CallCtx<'_>,r:f64,g:f64,b:f64,a:f64|{Ok(NewObject::normalized(Color::new(r,g,b,a)))});
object_binding!(graphics_color_rgb,"num_native_graphics_color_rgb","crate::kernels::computer_graphics::Color::rgb","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(r:f64,g:f64,b:f64)->NewObject<Color>,|_ctx:&CallCtx<'_>,r:f64,g:f64,b:f64|{Ok(NewObject::normalized(Color::rgb(r,g,b)))});
object_binding!(graphics_color_clamp,"num_native_graphics_color_clamp","crate::kernels::computer_graphics::Color::clamp","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Color>)->NewObject<Color>,|ctx:&CallCtx<'_>,a:Handle<Color>|{a.with(ctx,|a|Ok(NewObject::normalized(a.clamp())))});
object_binding!(graphics_color_lerp,"num_native_graphics_color_lerp","crate::kernels::computer_graphics::Color::lerp","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Color>,b:Handle<Color>,t:f64)->NewObject<Color>,|ctx:&CallCtx<'_>,a:Handle<Color>,b:Handle<Color>,t:f64|{a.with(ctx,|a|b.with(ctx,|b|Ok(NewObject::normalized(a.lerp(b,t)))))});
object_binding!(graphics_color_components,"num_native_graphics_color_components","","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Color>)->[f64;4],|ctx:&CallCtx<'_>,a:Handle<Color>|{a.with(ctx,|a|Ok([a.r,a.g,a.b,a.a]))});
object_binding!(graphics_dot_product_2d,"num_native_graphics_dot_product_2d","crate::kernels::computer_graphics::dot_product_2d","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Vector2D>,b:Handle<Vector2D>)->Dynamic,|ctx:&CallCtx<'_>,a:Handle<Vector2D>,b:Handle<Vector2D>|{let d=canonical(&(DimExpr::Mul(vec![object_dimension(&a,ctx)?,object_dimension(&b,ctx)?])))?;a.with(ctx,|a|b.with(ctx,|b|extended_quantity(g::dot_product_2d(a,b),d)))});
object_binding!(graphics_dot_product,"num_native_graphics_dot_product","crate::kernels::computer_graphics::dot_product","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Vector3D>,b:Handle<Vector3D>)->Dynamic,|ctx:&CallCtx<'_>,a:Handle<Vector3D>,b:Handle<Vector3D>|{let d=canonical(&(DimExpr::Mul(vec![object_dimension(&a,ctx)?,object_dimension(&b,ctx)?])))?;a.with(ctx,|a|b.with(ctx,|b|extended_quantity(g::dot_product(a,b),d)))});
object_binding!(graphics_dot_product_3d,"num_native_graphics_dot_product_3d","crate::kernels::computer_graphics::dot_product_3d","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Vector3D>,b:Handle<Vector3D>)->Dynamic,|ctx:&CallCtx<'_>,a:Handle<Vector3D>,b:Handle<Vector3D>|{let d=canonical(&(DimExpr::Mul(vec![object_dimension(&a,ctx)?,object_dimension(&b,ctx)?])))?;a.with(ctx,|a|b.with(ctx,|b|extended_quantity(g::dot_product_3d(a,b),d)))});
object_binding!(graphics_dot,"num_native_graphics_dot","crate::kernels::computer_graphics::dot","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Vector3D>,b:Handle<Vector3D>)->Dynamic,|ctx:&CallCtx<'_>,a:Handle<Vector3D>,b:Handle<Vector3D>|{let d=canonical(&(DimExpr::Mul(vec![object_dimension(&a,ctx)?,object_dimension(&b,ctx)?])))?;a.with(ctx,|a|b.with(ctx,|b|extended_quantity(g::dot(a,b),d)))});
object_binding!(graphics_cross_product_2d,"num_native_graphics_cross_product_2d","crate::kernels::computer_graphics::cross_product_2d","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Vector2D>,b:Handle<Vector2D>)->Dynamic,|ctx:&CallCtx<'_>,a:Handle<Vector2D>,b:Handle<Vector2D>|{let d=canonical(&(DimExpr::Mul(vec![object_dimension(&a,ctx)?,object_dimension(&b,ctx)?])))?;a.with(ctx,|a|b.with(ctx,|b|extended_quantity(g::cross_product_2d(a,b),d)))});
object_binding!(graphics_inner,"num_native_graphics_inner","crate::kernels::computer_graphics::InnerProduct::inner","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Vector3D>,b:Handle<Vector3D>)->Dynamic,|ctx:&CallCtx<'_>,a:Handle<Vector3D>,b:Handle<Vector3D>|{let d=canonical(&(DimExpr::Mul(vec![object_dimension(&a,ctx)?,object_dimension(&b,ctx)?])))?;a.with(ctx,|a|b.with(ctx,|b|extended_quantity(a.inner(b),d)))});
object_binding!(graphics_cross_product,"num_native_graphics_cross_product","crate::kernels::computer_graphics::cross_product","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Vector3D>,b:Handle<Vector3D>)->NewObject<Vector3D>,|ctx:&CallCtx<'_>,a:Handle<Vector3D>,b:Handle<Vector3D>|{let d=canonical(&(DimExpr::Mul(vec![object_dimension(&a,ctx)?,object_dimension(&b,ctx)?])))?;a.with(ctx,|a|b.with(ctx,|b|Ok(NewObject::new(g::cross_product(a,b),vec![d]))))});
object_binding!(graphics_angle_between,"num_native_graphics_angle_between","crate::kernels::computer_graphics::angle_between","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Vector3D>,b:Handle<Vector3D>)->f64,|ctx:&CallCtx<'_>,a:Handle<Vector3D>,b:Handle<Vector3D>|{a.with(ctx,|a|b.with(ctx,|b|Ok(g::angle_between(a,b))))});
object_binding!(graphics_project,"num_native_graphics_project","crate::kernels::computer_graphics::project","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Vector3D>,b:Handle<Vector3D>)->NewObject<Vector3D>,|ctx:&CallCtx<'_>,a:Handle<Vector3D>,b:Handle<Vector3D>|{let d=object_dimension(&a,ctx)?;a.with(ctx,|a|b.with(ctx,|b|Ok(NewObject::new(g::project(a,b),vec![d]))))});
object_binding!(graphics_reflect,"num_native_graphics_reflect","crate::kernels::computer_graphics::reflect","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Vector3D>,normal:Handle<Vector3D>)->NewObject<Vector3D>,|ctx:&CallCtx<'_>,a:Handle<Vector3D>,normal:Handle<Vector3D>|{normalized_direction(&normal,ctx)?;let d=object_dimension(&a,ctx)?;a.with(ctx,|a|normal.with(ctx,|n|Ok(NewObject::new(g::reflect(a,n),vec![d]))))});
object_binding!(graphics_refract,"num_native_graphics_refract","crate::kernels::computer_graphics::refract","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Vector3D>,normal:Handle<Vector3D>,eta:f64)->Maybe<NewObject<Vector3D>>,|ctx:&CallCtx<'_>,a:Handle<Vector3D>,normal:Handle<Vector3D>,eta:f64|{normalized_direction(&normal,ctx)?;if eta<=0.0{return Err(error("refractive index ratio must be positive"));}a.with(ctx,|a|normal.with(ctx,|n|{if a.magnitude()==0.0{return Err(error("incident direction must be nonzero"));}Ok(Maybe(g::refract(a,n,eta).map(NewObject::normalized)))}))});
object_binding!(graphics_lerp,"num_native_graphics_lerp","crate::kernels::computer_graphics::lerp","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Vector3D>,b:Handle<Vector3D>,t:f64)->NewObject<Vector3D>,|ctx:&CallCtx<'_>,a:Handle<Vector3D>,b:Handle<Vector3D>,t:f64|{let d=same(&a,&b,ctx)?;a.with(ctx,|a|b.with(ctx,|b|Ok(NewObject::new(g::lerp(a,b,t),vec![d]))))});
object_binding!(graphics_slerp,"num_native_graphics_slerp","crate::kernels::computer_graphics::slerp","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Vector3D>,b:Handle<Vector3D>,t:f64)->NewObject<Vector3D>,|ctx:&CallCtx<'_>,a:Handle<Vector3D>,b:Handle<Vector3D>,t:f64|{normalized_direction(&a,ctx)?;normalized_direction(&b,ctx)?;a.with(ctx,|a|b.with(ctx,|b|Ok(NewObject::normalized(g::slerp(a,b,t)))))});
object_binding!(graphics_quaternion_new,"num_native_graphics_quaternion_new","crate::kernels::computer_graphics::Quaternion::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(w:f64,x:f64,y:f64,z:f64)->NewObject<Quaternion>,|_ctx:&CallCtx<'_>,w:f64,x:f64,y:f64,z:f64|{Ok(NewObject::normalized(Quaternion::new(w,x,y,z)))});
object_binding!(graphics_quaternion_identity,"num_native_graphics_quaternion_identity","crate::kernels::computer_graphics::Quaternion::identity","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,()->NewObject<Quaternion>,|_ctx:&CallCtx<'_>|{Ok(NewObject::normalized(Quaternion::identity()))});
object_binding!(graphics_quaternion_from_axis_angle,"num_native_graphics_quaternion_from_axis_angle","crate::kernels::computer_graphics::Quaternion::from_axis_angle","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(axis:Handle<Vector3D>,angle:f64)->NewObject<Quaternion>,|ctx:&CallCtx<'_>,axis:Handle<Vector3D>,angle:f64|{normalized_direction(&axis,ctx)?;axis.with(ctx,|a|Ok(NewObject::normalized(Quaternion::from_axis_angle(a,angle))))});
object_binding!(graphics_quaternion_from_euler,"num_native_graphics_quaternion_from_euler","crate::kernels::computer_graphics::Quaternion::from_euler","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(roll:f64,pitch:f64,yaw:f64)->NewObject<Quaternion>,|_ctx:&CallCtx<'_>,roll:f64,pitch:f64,yaw:f64|{Ok(NewObject::normalized(Quaternion::from_euler(roll,pitch,yaw)))});
object_binding!(graphics_quaternion_magnitude,"num_native_graphics_quaternion_magnitude","crate::kernels::computer_graphics::Quaternion::magnitude","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Quaternion>)->f64,|ctx:&CallCtx<'_>,a:Handle<Quaternion>|{a.with(ctx,|a|{Ok(a.magnitude())})});
object_binding!(graphics_quaternion_normalize,"num_native_graphics_quaternion_normalize","crate::kernels::computer_graphics::Quaternion::normalize","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Quaternion>)->NewObject<Quaternion>,|ctx:&CallCtx<'_>,a:Handle<Quaternion>|{a.with(ctx,|a|{Ok(NewObject::normalized(a.normalize()))})});
object_binding!(graphics_quaternion_conjugate,"num_native_graphics_quaternion_conjugate","crate::kernels::computer_graphics::Quaternion::conjugate","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Quaternion>)->NewObject<Quaternion>,|ctx:&CallCtx<'_>,a:Handle<Quaternion>|{a.with(ctx,|a|{Ok(NewObject::normalized(a.conjugate()))})});
object_binding!(graphics_quaternion_inverse,"num_native_graphics_quaternion_inverse","crate::kernels::computer_graphics::Quaternion::inverse","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Quaternion>)->NewObject<Quaternion>,|ctx:&CallCtx<'_>,a:Handle<Quaternion>|{a.with(ctx,|a|{if a.magnitude()==0.0{return Err(error("quaternion inverse requires nonzero norm"));}Ok(NewObject::normalized(a.inverse()))})});
object_binding!(graphics_quaternion_multiply,"num_native_graphics_quaternion_multiply","crate::kernels::computer_graphics::Quaternion::multiply","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Quaternion>,b:Handle<Quaternion>)->NewObject<Quaternion>,|ctx:&CallCtx<'_>,a:Handle<Quaternion>,b:Handle<Quaternion>|{a.with(ctx,|a|b.with(ctx,|b|Ok(NewObject::normalized(a.multiply(b)))))});
object_binding!(graphics_quaternion_slerp,"num_native_graphics_quaternion_slerp","crate::kernels::computer_graphics::Quaternion::slerp","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Quaternion>,b:Handle<Quaternion>,t:f64)->NewObject<Quaternion>,|ctx:&CallCtx<'_>,a:Handle<Quaternion>,b:Handle<Quaternion>,t:f64|{a.with(ctx,|a|b.with(ctx,|b|Ok(NewObject::normalized(a.slerp(b,t)))))});
object_binding!(graphics_quaternion_rotate_vector,"num_native_graphics_quaternion_rotate_vector","crate::kernels::computer_graphics::Quaternion::rotate_vector","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(a:Handle<Quaternion>,vector:Handle<Vector3D>)->NewObject<Vector3D>,|ctx:&CallCtx<'_>,a:Handle<Quaternion>,vector:Handle<Vector3D>|{let d=object_dimension(&vector,ctx)?;a.with(ctx,|q|vector.with(ctx,|v|Ok(NewObject::new(q.rotate_vector(v),vec![d]))))});
object_binding!(graphics_quaternion_to_matrix,"num_native_graphics_quaternion_to_matrix","crate::kernels::computer_graphics::Quaternion::to_matrix","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Quaternion>)->NewObject<Frame>,|ctx:&CallCtx<'_>,a:Handle<Quaternion>|{a.with(ctx,|q|Ok(frame(q.to_matrix())))});
object_binding!(graphics_identity_matrix,"num_native_graphics_identity_matrix","crate::kernels::computer_graphics::identity_matrix","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,()->NewObject<Frame>,|_ctx:&CallCtx<'_>|{Ok(frame(g::identity_matrix()))});
object_binding!(graphics_scaling_matrix,"num_native_graphics_scaling_matrix","crate::kernels::computer_graphics::scaling_matrix","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(sx:f64,sy:f64,sz:f64)->NewObject<Frame>,|_ctx:&CallCtx<'_>,sx:f64,sy:f64,sz:f64|{Ok(frame(g::scaling_matrix(sx,sy,sz)))});
object_binding!(graphics_uniform_scaling_matrix,"num_native_graphics_uniform_scaling_matrix","crate::kernels::computer_graphics::uniform_scaling_matrix","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(scale:f64)->NewObject<Frame>,|_ctx:&CallCtx<'_>,scale:f64|{Ok(frame(g::uniform_scaling_matrix(scale)))});
object_binding!(graphics_rotation_matrix_x,"num_native_graphics_rotation_matrix_x","crate::kernels::computer_graphics::rotation_matrix_x","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(angle_rad:f64)->NewObject<Frame>,|_ctx:&CallCtx<'_>,angle_rad:f64|{Ok(frame(g::rotation_matrix_x(angle_rad)))});
object_binding!(graphics_rotation_matrix_y,"num_native_graphics_rotation_matrix_y","crate::kernels::computer_graphics::rotation_matrix_y","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(angle_rad:f64)->NewObject<Frame>,|_ctx:&CallCtx<'_>,angle_rad:f64|{Ok(frame(g::rotation_matrix_y(angle_rad)))});
object_binding!(graphics_rotation_matrix_z,"num_native_graphics_rotation_matrix_z","crate::kernels::computer_graphics::rotation_matrix_z","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(angle_rad:f64)->NewObject<Frame>,|_ctx:&CallCtx<'_>,angle_rad:f64|{Ok(frame(g::rotation_matrix_z(angle_rad)))});
object_binding!(graphics_shearing_matrix,"num_native_graphics_shearing_matrix","crate::kernels::computer_graphics::shearing_matrix","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(xy:f64,xz:f64,yx:f64,yz:f64,zx:f64,zy:f64)->NewObject<Frame>,|_ctx:&CallCtx<'_>,xy:f64,xz:f64,yx:f64,yz:f64,zx:f64,zy:f64|{Ok(frame(g::shearing_matrix(xy,xz,yx,yz,zx,zy)))});
object_binding!(graphics_rotation_matrix_axis,"num_native_graphics_rotation_matrix_axis","crate::kernels::computer_graphics::rotation_matrix_axis","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(axis:Handle<Vector3D>,angle_rad:f64)->NewObject<Frame>,|ctx:&CallCtx<'_>,axis:Handle<Vector3D>,angle_rad:f64|{normalized_direction(&axis,ctx)?;axis.with(ctx,|a|Ok(frame(g::rotation_matrix_axis(a,angle_rad))))});
object_binding!(graphics_translation_matrix,"num_native_graphics_translation_matrix","crate::kernels::computer_graphics::translation_matrix","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(dx:Input,dy:Input,dz:Input,coordinate_reference:Input)->NewObject<Frame>,|_ctx:&CallCtx<'_>,dx:Input,dy:Input,dz:Input,coordinate_reference:Input|{let(r,d)=api::reference(&coordinate_reference.0)?;let(x,xd)=scalar(&dx)?;let(y,yd)=scalar(&dy)?;let(z,zd)=scalar(&dz)?;for dim in [xd,yd,zd]{api::same_dimension(&d,&dim)?;}Ok(framed(g::translation_matrix(x/r,y/r,z/r),r,r,d.clone(),d,false))});
object_binding!(graphics_perspective_matrix,"num_native_graphics_perspective_matrix","crate::kernels::computer_graphics::perspective_matrix","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(fov_y:f64,aspect:f64,near:Input,far:Input,coordinate_reference:Input)->NewObject<Frame>,|_ctx:&CallCtx<'_>,fov_y:f64,aspect:f64,near:Input,far:Input,coordinate_reference:Input|{let(r,d)=api::reference(&coordinate_reference.0)?;let(n,nd)=scalar(&near)?;let(f,fd)=scalar(&far)?;api::same_dimension(&d,&nd)?;api::same_dimension(&d,&fd)?;if fov_y<=0.0||fov_y>=std::f64::consts::PI||aspect<=0.0||n<=0.0||f<=n{return Err(error("perspective requires0<fov<pi,aspect>0,0<near<far"));}Ok(framed(g::perspective_matrix(fov_y,aspect,n/r,f/r),r,1.0,d,DimExpr::Dimless,true))});
object_binding!(graphics_orthographic_matrix,"num_native_graphics_orthographic_matrix","crate::kernels::computer_graphics::orthographic_matrix","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(left:Input,right:Input,bottom:Input,top:Input,near:Input,far:Input,coordinate_reference:Input)->NewObject<Frame>,|_ctx:&CallCtx<'_>,left:Input,right:Input,bottom:Input,top:Input,near:Input,far:Input,coordinate_reference:Input|{let(r,d)=api::reference(&coordinate_reference.0)?;let read=|v:&Input|{let(x,e)=scalar(v)?;api::same_dimension(&d,&e)?;Ok::<_,OpError>(x/r)};let(l,rr,b,t,n,f)=(read(&left)?,read(&right)?,read(&bottom)?,read(&top)?,read(&near)?,read(&far)?);if l>=rr||b>=t||n>=f{return Err(error("orthographic bounds must increase"));}Ok(framed(g::orthographic_matrix(l,rr,b,t,n,f),r,1.0,d,DimExpr::Dimless,false))});
object_binding!(graphics_look_at_matrix,"num_native_graphics_look_at_matrix","crate::kernels::computer_graphics::look_at_matrix","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(eye:Handle<Vector3D>,center:Handle<Vector3D>,up:Handle<Vector3D>,coordinate_reference:Input)->NewObject<Frame>,|ctx:&CallCtx<'_>,eye:Handle<Vector3D>,center:Handle<Vector3D>,up:Handle<Vector3D>,coordinate_reference:Input|{let(r,d)=api::reference(&coordinate_reference.0)?;api::same_dimension(&d,&same(&eye,&center,ctx)?)?;normalized_direction(&up,ctx)?;eye.with(ctx,|e|center.with(ctx,|c|up.with(ctx,|u|{let e=scaled_vector(*e,1.0/r);let c=scaled_vector(*c,1.0/r);let delta=c-e;if delta.magnitude()==0.0||g::cross_product(&delta,u).magnitude()==0.0{return Err(error("look-at needs distinct eye/center and nonparallel up"));}Ok(framed(g::look_at_matrix(&e,&c,u),r,r,d.clone(),d,false))})))});
object_binding!(graphics_transform_point,"num_native_graphics_transform_point","crate::kernels::computer_graphics::transform_point","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(matrix:Handle<Frame>,value:Handle<Point3D>)->NewObject<Point3D>,|ctx:&CallCtx<'_>,matrix:Handle<Frame>,value:Handle<Point3D>|{let vd=object_dimension(&value,ctx)?;let dims=matrix.dimensions(ctx)?;matrix.with(ctx,|f|{let d=if f.framed{if dims.len()!=2{return Err(error("graphics frame dimension witnesses missing"));}api::same_dimension(&vd,&dims[0])?;dims[1].clone()}else{vd};value.with(ctx,|v|{let v=scaled_point(*v,1.0/f.coordinate);if f.projective{let w=f.matrix.get(3,0)*v.x+f.matrix.get(3,1)*v.y+f.matrix.get(3,2)*v.z+f.matrix.get(3,3);if w.abs()<=1e-10{return Err(error("projective point lies at/near infinite homogeneous coordinate"));}}let out=g::transform_point(&f.matrix,&v);Ok(NewObject::new(scaled_point(out,f.output),vec![d]))})})});
object_binding!(graphics_transform_vector,"num_native_graphics_transform_vector","crate::kernels::computer_graphics::transform_vector","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(matrix:Handle<Frame>,value:Handle<Vector3D>)->NewObject<Vector3D>,|ctx:&CallCtx<'_>,matrix:Handle<Frame>,value:Handle<Vector3D>|{let vd=object_dimension(&value,ctx)?;let dims=matrix.dimensions(ctx)?;matrix.with(ctx,|f|{let d=if f.framed{if dims.len()!=2{return Err(error("graphics frame dimension witnesses missing"));}api::same_dimension(&vd,&dims[0])?;dims[1].clone()}else{vd};value.with(ctx,|v|{let v=scaled_vector(*v,1.0/f.coordinate);let out=g::transform_vector(&f.matrix,&v);Ok(NewObject::new(scaled_vector(out,f.output),vec![d]))})})});
object_binding!(graphics_ray_new,"num_native_graphics_ray_new","crate::kernels::computer_graphics::Ray::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(origin:Handle<Point3D>,direction:Handle<Vector3D>,coordinate_reference:Input)->NewObject<GfxRay>,|ctx:&CallCtx<'_>,origin:Handle<Point3D>,direction:Handle<Vector3D>,coordinate_reference:Input|{let(r,d)=api::reference(&coordinate_reference.0)?;api::same_dimension(&d,&object_dimension(&origin,ctx)?)?;normalized_direction(&direction,ctx)?;origin.with(ctx,|o|direction.with(ctx,|v|Ok(NewObject::new(GfxRay{ray:Ray::new(scaled_point(*o,1.0/r),*v),reference:r},vec![d]))))});
object_binding!(graphics_ray_at,"num_native_graphics_ray_at","crate::kernels::computer_graphics::Ray::at","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(ray:Handle<GfxRay>,t:Input)->NewObject<Point3D>,|ctx:&CallCtx<'_>,ray:Handle<GfxRay>,t:Input|{let(t,d)=scalar(&t)?;api::same_dimension(&d,&object_dimension(&ray,ctx)?)?;ray.with(ctx,|r|Ok(NewObject::new(scaled_point(r.ray.at(t/r.reference),r.reference),vec![d])))});
object_binding!(graphics_sphere_new,"num_native_graphics_sphere_new","crate::kernels::computer_graphics::Sphere::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(center:Handle<Point3D>,radius:Input)->NewObject<Sphere>,|ctx:&CallCtx<'_>,center:Handle<Point3D>,radius:Input|{let(r,d)=scalar(&radius)?;api::same_dimension(&d,&object_dimension(&center,ctx)?)?;if r<=0.0{return Err(error("sphere radius must be positive"));}center.with(ctx,|c|Ok(NewObject::new(Sphere::new(*c,r),vec![d])))});
object_binding!(graphics_plane_new,"num_native_graphics_plane_new","crate::kernels::computer_graphics::Plane::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(point:Handle<Point3D>,normal:Handle<Vector3D>)->NewObject<Plane>,|ctx:&CallCtx<'_>,point:Handle<Point3D>,normal:Handle<Vector3D>|{normalized_direction(&normal,ctx)?;let d=object_dimension(&point,ctx)?;point.with(ctx,|p|normal.with(ctx,|n|Ok(NewObject::new(Plane::new(*p,*n),vec![d]))))});
object_binding!(graphics_ray_sphere_intersection,"num_native_graphics_ray_sphere_intersection","crate::kernels::computer_graphics::ray_sphere_intersection","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(ray:Handle<GfxRay>,sphere:Handle<Sphere>)->Maybe<Dynamic>,|ctx:&CallCtx<'_>,ray:Handle<GfxRay>,sphere:Handle<Sphere>|{let d=same(&ray,&sphere,ctx)?;ray.with(ctx,|r|sphere.with(ctx,|s|{let s=Sphere::new(scaled_point(s.center,1.0/r.reference),s.radius/r.reference);hit(g::ray_sphere_intersection(&r.ray,&s),r.reference,d,ctx)}))});
object_binding!(graphics_ray_plane_intersection,"num_native_graphics_ray_plane_intersection","crate::kernels::computer_graphics::ray_plane_intersection","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(ray:Handle<GfxRay>,plane:Handle<Plane>)->Maybe<Dynamic>,|ctx:&CallCtx<'_>,ray:Handle<GfxRay>,plane:Handle<Plane>|{let d=same(&ray,&plane,ctx)?;ray.with(ctx,|r|plane.with(ctx,|s|{let s=Plane::new(scaled_point(s.point,1.0/r.reference),s.normal);hit(g::ray_plane_intersection(&r.ray,&s),r.reference,d,ctx)}))});
object_binding!(graphics_ray_triangle_intersection,"num_native_graphics_ray_triangle_intersection","crate::kernels::computer_graphics::ray_triangle_intersection","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(ray:Handle<GfxRay>,v0:Handle<Point3D>,v1:Handle<Point3D>,v2:Handle<Point3D>)->Maybe<Dynamic>,|ctx:&CallCtx<'_>,ray:Handle<GfxRay>,v0:Handle<Point3D>,v1:Handle<Point3D>,v2:Handle<Point3D>|{let d=object_dimension(&ray,ctx)?;for h in [&v0,&v1,&v2]{api::same_dimension(&d,&object_dimension(h,ctx)?)?;}ray.with(ctx,|r|v0.with(ctx,|a|v1.with(ctx,|b|v2.with(ctx,|c|hit(g::ray_triangle_intersection(&r.ray,&scaled_point(*a,1.0/r.reference),&scaled_point(*b,1.0/r.reference),&scaled_point(*c,1.0/r.reference)),r.reference,d,ctx)))))});
object_binding!(graphics_prepared_triangle_new,"num_native_graphics_prepared_triangle_new","crate::kernels::computer_graphics::PreparedTriangle::new","Prepare three point handles with one coordinate dimension and a positive finite reference in that dimension. Degenerate triangles return None. The resulting handle owns its geometry.",false,(v0:Handle<Point3D>,v1:Handle<Point3D>,v2:Handle<Point3D>,coordinate_reference:Input)->Maybe<NewObject<GfxPreparedTriangle>>,|ctx:&CallCtx<'_>,v0:Handle<Point3D>,v1:Handle<Point3D>,v2:Handle<Point3D>,coordinate_reference:Input|{prepared_triangle([&v0,&v1,&v2],&coordinate_reference,ctx)});
object_binding!(graphics_prepared_triangle_intersect,"num_native_graphics_prepared_triangle_intersect","crate::kernels::computer_graphics::PreparedTriangle::intersect","Intersect a prepared triangle and a ray with matching coordinate dimensions. Inclusive bounds use the same dimension as the coordinates. The minimum is finite and nonnegative. A Some upper bound is finite and None is unbounded. The hit record contains physical t and point plus a dimensionless unit normal.",false,(triangle:Handle<GfxPreparedTriangle>,ray:Handle<GfxRay>,min_t:Input,max_t:api::Maybe<api::Input>)->Maybe<Dynamic>,|ctx:&CallCtx<'_>,triangle:Handle<GfxPreparedTriangle>,ray:Handle<GfxRay>,min_t:Input,max_t:api::Maybe<api::Input>|{prepared_triangle_hit(&triangle,&ray,&min_t,max_t.0.as_ref(),ctx)});
object_binding!(graphics_prepared_triangle_intersect_parallel,"num_native_graphics_prepared_triangle_intersect_parallel","crate::kernels::computer_graphics::PreparedTriangle::intersect_parallel","Intersect at most 4096 origin handles with a common dimensionless unit direction. Origins and inclusive bounds share the triangle coordinate dimension. The minimum is finite and nonnegative. A Some upper bound is finite and None is unbounded. The record has hit_mask and physical t of shape [n], physical points and dimensionless normals of shape [n,3]. Misses have zero placeholders and no result handles are retained.",false,(triangle:Handle<GfxPreparedTriangle>,origins:TriangleOrigins,direction:Handle<Vector3D>,min_t:Input,max_t:api::Maybe<api::Input>)->Dynamic,|ctx:&CallCtx<'_>,triangle:Handle<GfxPreparedTriangle>,origins:TriangleOrigins,direction:Handle<Vector3D>,min_t:Input,max_t:api::Maybe<api::Input>|{prepared_triangle_parallel_hits(&triangle,&origins,&direction,&min_t,max_t.0.as_ref(),ctx)});
object_binding!(graphics_bezier_quadratic,"num_native_graphics_bezier_quadratic","crate::kernels::computer_graphics::bezier_quadratic","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(p0:Handle<Point3D>,p1:Handle<Point3D>,p2:Handle<Point3D>,t:f64)->NewObject<Point3D>,|ctx:&CallCtx<'_>,p0:Handle<Point3D>,p1:Handle<Point3D>,p2:Handle<Point3D>,t:f64|{let d=object_dimension(&p0,ctx)?;api::same_dimension(&d,&object_dimension(&p1,ctx)?)?;api::same_dimension(&d,&object_dimension(&p2,ctx)?)?;p0.with(ctx,|a0|p1.with(ctx,|a1|p2.with(ctx,|a2|Ok(NewObject::new(g::bezier_quadratic(a0,a1,a2,t),vec![d])))))});
object_binding!(graphics_bezier_cubic,"num_native_graphics_bezier_cubic","crate::kernels::computer_graphics::bezier_cubic","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(p0:Handle<Point3D>,p1:Handle<Point3D>,p2:Handle<Point3D>,p3:Handle<Point3D>,t:f64)->NewObject<Point3D>,|ctx:&CallCtx<'_>,p0:Handle<Point3D>,p1:Handle<Point3D>,p2:Handle<Point3D>,p3:Handle<Point3D>,t:f64|{let d=object_dimension(&p0,ctx)?;api::same_dimension(&d,&object_dimension(&p1,ctx)?)?;api::same_dimension(&d,&object_dimension(&p2,ctx)?)?;api::same_dimension(&d,&object_dimension(&p3,ctx)?)?;p0.with(ctx,|a0|p1.with(ctx,|a1|p2.with(ctx,|a2|p3.with(ctx,|a3|Ok(NewObject::new(g::bezier_cubic(a0,a1,a2,a3,t),vec![d]))))))});
object_binding!(graphics_catmull_rom,"num_native_graphics_catmull_rom","crate::kernels::computer_graphics::catmull_rom","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(p0:Handle<Point3D>,p1:Handle<Point3D>,p2:Handle<Point3D>,p3:Handle<Point3D>,t:f64)->NewObject<Point3D>,|ctx:&CallCtx<'_>,p0:Handle<Point3D>,p1:Handle<Point3D>,p2:Handle<Point3D>,p3:Handle<Point3D>,t:f64|{let d=object_dimension(&p0,ctx)?;api::same_dimension(&d,&object_dimension(&p1,ctx)?)?;api::same_dimension(&d,&object_dimension(&p2,ctx)?)?;api::same_dimension(&d,&object_dimension(&p3,ctx)?)?;p0.with(ctx,|a0|p1.with(ctx,|a1|p2.with(ctx,|a2|p3.with(ctx,|a3|Ok(NewObject::new(g::catmull_rom(a0,a1,a2,a3,t),vec![d]))))))});
object_binding!(graphics_barycentric_coordinates,"num_native_graphics_barycentric_coordinates","crate::kernels::computer_graphics::barycentric_coordinates","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(p0:Handle<Point3D>,p1:Handle<Point3D>,p2:Handle<Point3D>,p3:Handle<Point3D>)->(f64,f64,f64),|ctx:&CallCtx<'_>,p0:Handle<Point3D>,p1:Handle<Point3D>,p2:Handle<Point3D>,p3:Handle<Point3D>|{let d=object_dimension(&p0,ctx)?;api::same_dimension(&d,&object_dimension(&p1,ctx)?)?;api::same_dimension(&d,&object_dimension(&p2,ctx)?)?;api::same_dimension(&d,&object_dimension(&p3,ctx)?)?;p0.with(ctx,|a0|p1.with(ctx,|a1|p2.with(ctx,|a2|p3.with(ctx,|a3|Ok(g::barycentric_coordinates(a0,a1,a2,a3))))))});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &graphics_point2_new::BINDING,
            &graphics_point2_components::BINDING,
            &graphics_point3_new::BINDING,
            &graphics_point3_components::BINDING,
            &graphics_vector2_new::BINDING,
            &graphics_vector2_components::BINDING,
            &graphics_vector3_new::BINDING,
            &graphics_vector3_components::BINDING,
            &graphics_point2_distance_to::BINDING,
            &graphics_point3_distance_to::BINDING,
            &graphics_point3_to_vector::BINDING,
            &graphics_vector2_magnitude::BINDING,
            &graphics_vector2_normalize::BINDING,
            &graphics_vector3_magnitude::BINDING,
            &graphics_vector3_normalize::BINDING,
            &graphics_vector3_magnitude_squared::BINDING,
            &graphics_vector2_rotate::BINDING,
            &graphics_vector2_perpendicular::BINDING,
            &graphics_color_new::BINDING,
            &graphics_color_rgb::BINDING,
            &graphics_color_clamp::BINDING,
            &graphics_color_lerp::BINDING,
            &graphics_color_components::BINDING,
            &graphics_dot_product_2d::BINDING,
            &graphics_dot_product::BINDING,
            &graphics_dot_product_3d::BINDING,
            &graphics_dot::BINDING,
            &graphics_cross_product_2d::BINDING,
            &graphics_inner::BINDING,
            &graphics_cross_product::BINDING,
            &graphics_angle_between::BINDING,
            &graphics_project::BINDING,
            &graphics_reflect::BINDING,
            &graphics_refract::BINDING,
            &graphics_lerp::BINDING,
            &graphics_slerp::BINDING,
            &graphics_quaternion_new::BINDING,
            &graphics_quaternion_identity::BINDING,
            &graphics_quaternion_from_axis_angle::BINDING,
            &graphics_quaternion_from_euler::BINDING,
            &graphics_quaternion_magnitude::BINDING,
            &graphics_quaternion_normalize::BINDING,
            &graphics_quaternion_conjugate::BINDING,
            &graphics_quaternion_inverse::BINDING,
            &graphics_quaternion_multiply::BINDING,
            &graphics_quaternion_slerp::BINDING,
            &graphics_quaternion_rotate_vector::BINDING,
            &graphics_quaternion_to_matrix::BINDING,
            &graphics_identity_matrix::BINDING,
            &graphics_scaling_matrix::BINDING,
            &graphics_uniform_scaling_matrix::BINDING,
            &graphics_rotation_matrix_x::BINDING,
            &graphics_rotation_matrix_y::BINDING,
            &graphics_rotation_matrix_z::BINDING,
            &graphics_shearing_matrix::BINDING,
            &graphics_rotation_matrix_axis::BINDING,
            &graphics_translation_matrix::BINDING,
            &graphics_perspective_matrix::BINDING,
            &graphics_orthographic_matrix::BINDING,
            &graphics_look_at_matrix::BINDING,
            &graphics_transform_point::BINDING,
            &graphics_transform_vector::BINDING,
            &graphics_ray_new::BINDING,
            &graphics_ray_at::BINDING,
            &graphics_sphere_new::BINDING,
            &graphics_plane_new::BINDING,
            &graphics_ray_sphere_intersection::BINDING,
            &graphics_ray_plane_intersection::BINDING,
            &graphics_ray_triangle_intersection::BINDING,
            &graphics_prepared_triangle_new::BINDING,
            &graphics_prepared_triangle_intersect::BINDING,
            &graphics_prepared_triangle_intersect_parallel::BINDING,
            &graphics_bezier_quadratic::BINDING,
            &graphics_bezier_cubic::BINDING,
            &graphics_catmull_rom::BINDING,
            &graphics_barycentric_coordinates::BINDING,
        ],
        "numerical_vector",
    )
}
