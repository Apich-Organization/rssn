//! Native interpolation closures retain independent coordinate and ordinate references.
use super::Dynamic;
use super::Handle;
use super::Input;
use super::Installer;
use super::LIMIT;
use super::NativeType;
use super::NewObject;
use super::Normalized;
use super::OpError;
use super::RuleError;
use super::Value;
use super::api;
use super::err;
use super::object_binding;
use crate::kernels::interpolate;
use crate::kernels::polynomial::Polynomial;
struct Curve {
    data: std::sync::Arc<dyn Fn(f64) -> f64>,
    xscale: f64,
    yscale: f64,
}
impl NativeType for Curve {
    const KIND: &'static str = "rssn/native/interpolation-curve/v1";
    const MAX_BYTES: usize = 262144;

    fn validate(&self) -> Result<(), OpError> {
        if !self.xscale.is_finite()
            || !self.yscale.is_finite()
            || self.xscale <= 0.0
            || self.yscale <= 0.0
        {
            return Err(err("interpolation reference must be positive finite"));
        }
        Ok(())
    }
}
fn points(data: &[(f64, f64)]) -> Result<(), OpError> {
    if data.len() < 2 || data.len() > 129 || data.windows(2).any(|p| p[0].0 >= p[1].0) {
        return Err(err(
            "interpolation needs 2..=129 strictly increasing normalized points",
        ));
    }
    Ok(())
}
object_binding!(interp_linear_interpolation,"num_linear_linear_interpolation","crate::kernels::interpolate::linear_interpolation","Native interpolation over normalized samples with explicit coordinate and ordinate references. Physical queries and outputs are checked against retained domains; point/order/allocation bounds are enforced.",false,(data:Normalized<Vec<(f64,f64)>>,coordinate_reference:Input,value_reference:Input)->NewObject<Curve>,|ctx,data,coordinate_reference,value_reference|{let _=ctx;points(&data.0)?;let(x,dx)=api::reference(&coordinate_reference.0)?;let(y,dy)=api::reference(&value_reference.0)?;Ok(NewObject::new(Curve{data:interpolate::linear_interpolation(&data.0).map_err(err)?,xscale:x,yscale:y},vec![dx,dy]))});
object_binding!(interp_cubic_spline_interpolation,"num_linear_cubic_spline_interpolation","crate::kernels::interpolate::cubic_spline_interpolation","Native interpolation over normalized samples with explicit coordinate and ordinate references. Physical queries and outputs are checked against retained domains; point/order/allocation bounds are enforced.",false,(data:Normalized<Vec<(f64,f64)>>,coordinate_reference:Input,value_reference:Input)->NewObject<Curve>,|ctx,data,coordinate_reference,value_reference|{let _=ctx;points(&data.0)?;let(x,dx)=api::reference(&coordinate_reference.0)?;let(y,dy)=api::reference(&value_reference.0)?;Ok(NewObject::new(Curve{data:interpolate::cubic_spline_interpolation(&data.0).map_err(err)?,xscale:x,yscale:y},vec![dx,dy]))});
object_binding!(interp_pchip_interpolation,"num_linear_pchip_interpolation","crate::kernels::interpolate::pchip_interpolation","Native interpolation over normalized samples with explicit coordinate and ordinate references. Physical queries and outputs are checked against retained domains; point/order/allocation bounds are enforced.",false,(data:Normalized<Vec<(f64,f64)>>,coordinate_reference:Input,value_reference:Input)->NewObject<Curve>,|ctx,data,coordinate_reference,value_reference|{let _=ctx;points(&data.0)?;let(x,dx)=api::reference(&coordinate_reference.0)?;let(y,dy)=api::reference(&value_reference.0)?;Ok(NewObject::new(Curve{data:interpolate::pchip_interpolation(&data.0).map_err(err)?,xscale:x,yscale:y},vec![dx,dy]))});
object_binding!(interp_cubic_spline_with,"num_linear_cubic_spline_with","crate::kernels::interpolate::cubic_spline_with","Native interpolation over normalized samples with explicit coordinate and ordinate references. Physical queries and outputs are checked against retained domains; point/order/allocation bounds are enforced.",false,(data:Normalized<Vec<(f64,f64)>>,coordinate_reference:Input,value_reference:Input,end:String,left:Normalized<f64>,right:Normalized<f64>)->NewObject<Curve>,|ctx,data,coordinate_reference,value_reference,end,left,right|{let _=ctx;points(&data.0)?;let(x,dx)=api::reference(&coordinate_reference.0)?;let(y,dy)=api::reference(&value_reference.0)?;let ends=match end.as_str(){"Natural"=>interpolate::SplineEnd::Natural,"NotAKnot"=>interpolate::SplineEnd::NotAKnot,"Clamped"=>interpolate::SplineEnd::Clamped(left.0,right.0),_=>return Err(err("spline end must be Natural, NotAKnot or Clamped"))};Ok(NewObject::new(Curve{data:interpolate::cubic_spline_with(&data.0,ends).map_err(err)?,xscale:x,yscale:y},vec![dx,dy]))});
object_binding!(interp_curve_eval,"num_linear_curve_eval","","Native interpolation over normalized samples with explicit coordinate and ordinate references. Physical queries and outputs are checked against retained domains; point/order/allocation bounds are enforced.",false,(curve:Handle<Curve>,point:Input)->Dynamic,|ctx,curve,point|{let _=ctx;let domains=curve.dimensions(ctx)?;if domains.len()!=2{return Err(err("curve domains missing"));}api::same_dimension(&api::dimension(&point.0)?,&domains[0])?;let x=crate::rules::numerical_abi::numeric(&point.0,1)?;curve.with(ctx,|c|Ok(Dynamic(api::quantity((c.data)(x/c.xscale)*c.yscale,domains[1].clone())?)))});
object_binding!(interp_lagrange_interpolation,"num_linear_lagrange_interpolation","crate::kernels::interpolate::lagrange_interpolation","Native normalized-only convenience API: coordinates and ordinates are strictly dimensionless. Corresponding physical adapters accept explicit independent coordinate and ordinate references. Point/order/allocation bounds are enforced.",true,(data:Vec<(f64,f64)>)->NewObject<Polynomial>,|ctx,data|{let _=ctx;points(&data)?;Ok(NewObject::normalized(interpolate::lagrange_interpolation(&data).map_err(err)?))});
object_binding!(interp_interpolate_nd,"num_linear_interpolate_nd","crate::kernels::interpolate::interpolate_nd","Native normalized-only convenience API: coordinates and ordinates are strictly dimensionless. Corresponding physical adapters accept explicit independent coordinate and ordinate references. Point/order/allocation bounds are enforced.",true,(axes:Vec<Vec<f64>>,values:Vec<f64>,point:Vec<f64>,order:String)->f64,|ctx,axes,values,point,order|{let _=ctx;if axes.is_empty()||axes.len()>6||axes.iter().any(|a|a.len()>129)||axes.iter().try_fold(1usize,|n,a|n.checked_mul(a.len()))!=Some(values.len())||values.len()>LIMIT{return Err(err("ND interpolation shape/allocation budget exceeded"));}let order=match order.as_str(){"Linear"=>interpolate::TensorOrder::Linear,"CubicSpline"=>interpolate::TensorOrder::CubicSpline,_=>return Err(err("interpolation order must be Linear or CubicSpline"))};interpolate::interpolate_nd(&axes,&values,&point,order).map_err(err)});
object_binding!(interp_delaunay_linear_2d,"num_linear_delaunay_linear_2d","crate::kernels::interp_ext::delaunay_linear_2d","Native normalized-only convenience API: coordinates and ordinates are strictly dimensionless. Corresponding physical adapters accept explicit independent coordinate and ordinate references. Point/order/allocation bounds are enforced.",true,(data:Vec<(f64,f64)>,values:Vec<f64>,x:f64,y:f64)->f64,|ctx,data,values,x,y|{let _=ctx;if data.len()>1024{return Err(err("scattered interpolation size/work budget exceeded"));}crate::kernels::interp_ext::delaunay_linear_2d(&data,&values,x,y).map_err(err)});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(i, BINDINGS, "numerical_matrix")
}
static BINDINGS: &[&api::Binding] = &[
    &lagrange_physical::BINDING,
    &interpolate_nd_physical::BINDING,
    &delaunay_physical::BINDING,
    &interp_linear_interpolation::BINDING,
    &interp_cubic_spline_interpolation::BINDING,
    &interp_pchip_interpolation::BINDING,
    &interp_cubic_spline_with::BINDING,
    &interp_curve_eval::BINDING,
    &interp_lagrange_interpolation::BINDING,
    &interp_interpolate_nd::BINDING,
    &interp_delaunay_linear_2d::BINDING,
];
object_binding!(lagrange_physical,"num_linear_lagrange_physical","crate::kernels::interpolate::lagrange_interpolation","Native Lagrange interpolation with explicit independent X/Y references; normalized data define a polynomial in x/X and y/Y. Retained query closure checks physical X inputs and restores physical Y outputs.",false,(data:Normalized<Vec<(f64,f64)>>,coordinate_reference:Input,value_reference:Input)->NewObject<Curve>,|ctx,data,coordinate_reference,value_reference|{let _=ctx;points(&data.0)?;let(x,dx)=api::reference(&coordinate_reference.0)?;let(y,dy)=api::reference(&value_reference.0)?;let p=interpolate::lagrange_interpolation(&data.0).map_err(err)?;if p.coeffs.iter().any(|v|!v.is_finite()){return Err(err("Lagrange coefficient arithmetic overflow"));}Ok(NewObject::new(Curve{data:std::sync::Arc::new(move|x|p.eval(x)),xscale:x,yscale:y},vec![dx,dy]))});
object_binding!(interpolate_nd_physical,"num_linear_interpolate_nd_physical","crate::kernels::interpolate::interpolate_nd","Native ND interpolation with separate positive physical references for every axis and the ordinate. Axis references and the query are records named axis0..axisN; this permits heterogeneous axis dimensions without mixing units in a numeric tensor. Sample axes/values are normalized by those explicit references.",false,(axes:Normalized<Vec<Vec<f64>>>,values:Normalized<Vec<f64>>,point:Input,coordinate_references:Input,value_reference:Input,order:String)->Dynamic,|ctx,axes,values,point,coordinate_references,value_reference,order|{let _=ctx;let n=axes.0.len();if n==0||n>6||axes.0.iter().any(|a|a.len()>129)||axes.0.iter().try_fold(1usize,|n,a|n.checked_mul(a.len()))!=Some(values.0.len())||values.0.len()>LIMIT{return Err(err("ND interpolation shape/allocation budget exceeded"));}let Value::Record(points)=&point.0 else{return Err(err("physical ND point must be a record"));};let Value::Record(references)=&coordinate_references.0 else{return Err(err("physical ND references must be a record"));};
if points.fields.len()!=n||references.fields.len()!=n{return Err(err("ND axis reference/query field count mismatch"));}let mut query=Vec::new();for axis in 0..n{let name=format!("axis{axis}");let reference=api::field(&coordinate_references.0,&name)?;let(s,d)=api::reference(reference)?;let q=api::field(&point.0,&name)?;api::same_dimension(&d,&api::dimension(q)?)?;query.push(crate::rules::numerical_abi::numeric(q,2)?/s);}let(s,d)=api::reference(&value_reference.0)?;let order=match order.as_str(){"Linear"=>interpolate::TensorOrder::Linear,"CubicSpline"=>interpolate::TensorOrder::CubicSpline,_=>return Err(err("interpolation order must be Linear or CubicSpline"))};Ok(Dynamic(api::quantity(interpolate::interpolate_nd(&axes.0,&values.0,&query,order).map_err(err)?*s,d)?))});
object_binding!(delaunay_physical,"num_linear_delaunay_linear_2d_physical","crate::kernels::interp_ext::delaunay_linear_2d","Native Delaunay interpolation in the explicitly scaled coordinate metric: samples are normalized by independent positive X/Y references, queries carry the corresponding axis units, and values use a positive ordinate reference. Changing unequal coordinate reference scales intentionally changes the triangulation metric.",false,(data:Normalized<Vec<(f64,f64)>>,values:Normalized<Vec<f64>>,x:Input,y:Input,x_reference:Input,y_reference:Input,value_reference:Input)->Dynamic,|ctx,data,values,x,y,x_reference,y_reference,value_reference|{let _=ctx;if data.0.len()>1024{return Err(err("Delaunay interpolation work budget exceeded"));}let(sx,dx)=api::reference(&x_reference.0)?;let(sy,dy)=api::reference(&y_reference.0)?;api::same_dimension(&dx,&api::dimension(&x.0)?)?;api::same_dimension(&dy,&api::dimension(&y.0)?)?;let(sv,dv)=api::reference(&value_reference.0)?;let x=crate::rules::numerical_abi::numeric(&x.0,2)?/sx;let y=crate::rules::numerical_abi::numeric(&y.0,3)?/sy;Ok(Dynamic(api::quantity(crate::kernels::interp_ext::delaunay_linear_2d(&data.0,&values.0,x,y).map_err(err)?*sv,dv)?))});
