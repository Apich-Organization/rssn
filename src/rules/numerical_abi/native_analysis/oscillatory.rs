//! Physical multidimensional and oscillatory native integrators on compiled functions.
use crate::graph::RuleError;
use crate::graph::rule::Installer;
use crate::kernel::CallCtx;
use crate::kernel::OpError;
use crate::kernel::Value;
use crate::kernel::signature::DimExpr;
use crate::kernel::signature::Exp;
use crate::kernel::signature::Rational;
use crate::kernel::value::ComplexValue;
use crate::kernel::value::DynValue;
use crate::kernel::value::VariantValue;
use crate::kernels::quadrature as q;
use crate::rules::numerical_abi::codecs::Arg;
use crate::rules::numerical_abi::codecs::Output;
use crate::rules::numerical_abi::native_analysis_parallel::bounds;
use crate::rules::numerical_abi::native_analysis_parallel::limit;
use crate::rules::numerical_abi::native_analysis_parallel::output;
use crate::rules::numerical_abi::native_backend::compiled_snapshot;
use crate::rules::numerical_abi::native_object_api::Dynamic;
use crate::rules::numerical_abi::native_object_api::Handle;
use crate::rules::numerical_abi::native_object_api::Input;
use crate::rules::numerical_abi::native_object_api::object_binding;
use crate::rules::numerical_abi::native_object_api::{
    self as object,
};
fn complex(
    value: num_complex::Complex64,
    d: DimExpr,
) -> Result<Value, OpError> {
    if !value.re.is_finite() || !value.im.is_finite() {
        return Err(object::error("oscillatory native result is not finite"));
    }
    let value = Value::Complex(ComplexValue::Float {
        re: value.re,
        im: value.im,
    });
    Ok(Value::Dyn(DynValue {
        value: Box::new(value),
        dim: d,
    }))
}
fn optional_error(
    error: f64,
    d: DimExpr,
) -> Result<Value, OpError> {
    Ok(Value::Variant(VariantValue {
        tag: if error.is_finite() {
            "Some"
        } else {
            "None"
        }
        .into(),
        payload: if error.is_finite() {
            vec![object::quantity(error, d)?]
        } else {
            vec![]
        },
    }))
}
fn complex_result(
    r: q::OscillatoryIntegral,
    d: DimExpr,
) -> Result<Dynamic, OpError> {
    Ok(Dynamic(object::record(vec![
        ("value", complex(r.value, d.clone())?),
        ("error", optional_error(r.error, d)?),
        ("evaluations", Output::encode(r.evaluations)?),
        ("converged", Output::encode(r.converged)?),
    ])))
}
object_binding!(triangle,"num_analysis_integrate_triangle_compiled","crate::kernels::quadrature::integrate_triangle","Native adaptive triangle quadrature on three physical planar vertices; compiled function has matching coordinate references and absolute tolerance has area times ordinate units.",false,(vertices:Input,absolute:Input,relative:f64,budget:usize,function:Input)->Dynamic,|ctx:&CallCtx<'_>,vertices:Input,absolute:Input,relative:f64,budget:usize,function:Input|{
 limit(budget)?;let coordinate=object::dimension(&vertices.0)?;let vertices=<Vec<Vec<f64>> as Arg>::decode(&vertices.0,0)?;
 if vertices.len()!=3||vertices.iter().any(|r|r.len()!=2){return Err(object::error("triangle needs a3by2 vertex tensor"));}
 let f=compiled_snapshot(ctx,&function.0,budget)?;if f.input_dimensions.len()!=2{return Err(object::error("triangle integrand needs two inputs"));}for d in&f.input_dimensions{object::same_dimension(&coordinate,d)?;}
 let d=DimExpr::Mul(vec![DimExpr::Pow(Box::new(coordinate),Exp::Rat(Rational::new(2,1))),f.output_dimension.clone()]);object::same_dimension(&d,&object::dimension(&absolute.0)?)?;let abs=super::super::numeric(&absolute.0,1)?;
 if abs<0.0||relative<0.0||abs==0.0&&relative==0.0{return Err(object::error("triangle tolerances must be nonnegative and one positive"));}
 let v=[[vertices[0][0],vertices[0][1]],[vertices[1][0],vertices[1][1]],[vertices[2][0],vertices[2][1]]];let r=q::integrate_triangle(&|x,y|f.sample(&[x,y]),v,abs,relative,budget);f.finish()?;output(r,d)
});
object_binding!(importance,"num_analysis_importance_sampling_gaussian_compiled","crate::kernels::quadrature::importance_sampling_gaussian","Native Gaussian importance sampling using measured target/proposal models; integrand references match their coordinate domain. Output has ordinate units, not ordinate times length; ESS remains dimensionless.",false,(target:Handle<q::Gaussian>,proposal:Handle<q::Gaussian>,samples:usize,seed:u64,budget:usize,function:Input)->Dynamic,|ctx:&CallCtx<'_>,target:Handle<q::Gaussian>,proposal:Handle<q::Gaussian>,samples:usize,seed:u64,budget:usize,function:Input|{
 limit(samples)?;let td=target.dimensions(ctx)?.into_iter().next().ok_or_else(||object::error("Gaussian lacks coordinate domain"))?;let pd=proposal.dimensions(ctx)?.into_iter().next().ok_or_else(||object::error("Gaussian lacks coordinate domain"))?;object::same_dimension(&td,&pd)?;
 let f=compiled_snapshot(ctx,&function.0,budget)?;if f.input_dimensions.len()!=1{return Err(object::error("importance integrand must have one input"));}object::same_dimension(&td,&f.input_dimensions[0])?;
 let target=target.with(ctx,|v|Ok(*v))?;let proposal=proposal.with(ctx,|v|Ok(*v))?;let r=q::importance_sampling_gaussian(|x|f.sample(&[x]),target,proposal,samples,seed);f.finish()?;
 let Dynamic(integral)=output(r.integral,f.output_dimension)?;Ok(Dynamic(object::record(vec![("integral",integral),("ess",Output::encode(r.ess)?)])))
});
object_binding!(filon_real,"num_analysis_filon2d_compiled","crate::kernels::quadrature::filon2d","Native2D Filon cosine/sine quadrature: phase coefficients have one dimension P,omega has1/(P*coordinate),amplitude references preserve arbitrary ordinate units.",false,(lo:Input,hi:Input,omega:Input,phase:Input,cosine:bool,tolerance:Input,budget:usize,function:Input)->Dynamic,|ctx:&CallCtx<'_>,lo:Input,hi:Input,omega:Input,phase:Input,cosine:bool,tolerance:Input,budget:usize,function:Input|{
 let coordinate=object::dimension(&lo.0)?;let phase_dim=object::dimension(&phase.0)?;let(lo,hi,f,d)=bounds(ctx,&lo,&hi,&function,budget)?;if lo.len()!=2{return Err(object::error("Filon2D needs two coordinates"));}
 object::same_dimension(&DimExpr::Inverse(Box::new(DimExpr::Mul(vec![coordinate,phase_dim]))),&object::dimension(&omega.0)?)?;object::same_dimension(&d,&object::dimension(&tolerance.0)?)?;
 let phase=<Vec<f64> as Arg>::decode(&phase.0,3)?;if phase.len()!=2{return Err(object::error("linear phase needs two coefficients"));}let w=super::super::numeric(&omega.0,2)?;let tol=super::super::numeric(&tolerance.0,5)?;if tol<=0.0{return Err(object::error("Filon absolute tolerance must be positive"));}
 let r=q::filon2d(|x,y|f.sample(&[x,y]),&lo,&hi,w,&phase,cosine,tol);f.finish()?;output(r,d)
});
object_binding!(filon_complex,"num_analysis_filon2d_complex_compiled","crate::kernels::quadrature::filon2d_complex","Native complex2D Filon quadrature with physical coordinate,amplitude,linear-phase and frequency contracts.",false,(lo:Input,hi:Input,omega:Input,phase:Input,tolerance:Input,budget:usize,function:Input)->Dynamic,|ctx:&CallCtx<'_>,lo:Input,hi:Input,omega:Input,phase:Input,tolerance:Input,budget:usize,function:Input|{
 let coordinate=object::dimension(&lo.0)?;let phase_dim=object::dimension(&phase.0)?;let(lo,hi,f,d)=bounds(ctx,&lo,&hi,&function,budget)?;if lo.len()!=2{return Err(object::error("Filon2D needs two coordinates"));}
 object::same_dimension(&DimExpr::Inverse(Box::new(DimExpr::Mul(vec![coordinate,phase_dim]))),&object::dimension(&omega.0)?)?;object::same_dimension(&d,&object::dimension(&tolerance.0)?)?;
 let phase=<Vec<f64> as Arg>::decode(&phase.0,3)?;if phase.len()!=2{return Err(object::error("linear phase needs two coefficients"));}let w=super::super::numeric(&omega.0,2)?;let tol=super::super::numeric(&tolerance.0,4)?;if tol<=0.0{return Err(object::error("Filon absolute tolerance must be positive"));}
 let r=q::filon2d_complex(|x,y|f.sample(&[x,y]),&lo,&hi,w,&phase,tol);f.finish()?;complex_result(r,d)
});
object_binding!(oscillatory,"num_analysis_oscillatory2d_compiled","crate::kernels::quadrature::oscillatory2d","Native general-phase2D oscillatory cubature; owned phase program references match physical coordinates and omega has inverse phase dimension. Evaluation limit applies separately to each compiled program.",false,(lo:Input,hi:Input,omega:Input,tolerance:Input,budget:usize,function:Input,phase:Input)->Dynamic,|ctx:&CallCtx<'_>,lo:Input,hi:Input,omega:Input,tolerance:Input,budget:usize,function:Input,phase:Input|{
 let coordinate=object::dimension(&lo.0)?;let(lo,hi,f,d)=bounds(ctx,&lo,&hi,&function,budget)?;if lo.len()!=2{return Err(object::error("oscillatory2D requires two coordinates"));}let phase=compiled_snapshot(ctx,&phase.0,budget)?;if phase.input_dimensions.len()!=2{return Err(object::error("phase program must have two inputs"));}for p in&phase.input_dimensions{object::same_dimension(&coordinate,p)?;}
 object::same_dimension(&DimExpr::Inverse(Box::new(phase.output_dimension.clone())),&object::dimension(&omega.0)?)?;object::same_dimension(&d,&object::dimension(&tolerance.0)?)?;let w=super::super::numeric(&omega.0,2)?;let tol=super::super::numeric(&tolerance.0,3)?;if tol<=0.0{return Err(object::error("oscillatory absolute tolerance must be positive"));}
 let r=q::oscillatory2d(|x,y|f.sample(&[x,y]),|x,y|phase.sample(&[x,y]),&lo,&hi,w,tol);f.finish()?;phase.finish()?;complex_result(r,d)
});
#[allow(
    clippy::needless_pass_by_value,
    clippy::too_many_arguments,
    reason = "The typed native ABI transfers decoded arguments by value; this shared dispatcher mirrors its explicit source signature."
)]
fn levin(
    ctx: &CallCtx<'_>,
    lo: Input,
    hi: Input,
    omega: Input,
    tolerance: Input,
    points: usize,
    budget: usize,
    function: Input,
    phase: Input,
    derivative: Input,
    adaptive: bool,
) -> Result<Dynamic, OpError> {
    if !(2..=512).contains(&points) {
        return Err(object::error("Levin collocation points must be2..512"));
    }
    let coordinate = object::dimension(&lo.0)?;
    object::same_dimension(&coordinate, &object::dimension(&hi.0)?)?;
    let f = compiled_snapshot(ctx, &function.0, budget)?;
    let phase = compiled_snapshot(ctx, &phase.0, budget)?;
    let derivative = compiled_snapshot(ctx, &derivative.0, budget)?;
    for p in [&f, &phase, &derivative] {
        if p.input_dimensions.len() != 1 {
            return Err(object::error("Levin programs must have one input"));
        }
        object::same_dimension(&coordinate, &p.input_dimensions[0])?;
    }
    object::same_dimension(
        &DimExpr::Mul(vec![
            phase.output_dimension.clone(),
            DimExpr::Inverse(Box::new(coordinate.clone())),
        ]),
        &derivative.output_dimension,
    )?;
    object::same_dimension(
        &DimExpr::Inverse(Box::new(phase.output_dimension.clone())),
        &object::dimension(&omega.0)?,
    )?;
    let d = DimExpr::Mul(vec![coordinate, f.output_dimension.clone()]);
    object::same_dimension(&d, &object::dimension(&tolerance.0)?)?;
    let a = super::numeric(&lo.0, 0)?;
    let b = super::numeric(&hi.0, 1)?;
    let w = super::numeric(&omega.0, 2)?;
    let tol = super::numeric(&tolerance.0, 3)?;
    if tol <= 0.0 {
        return Err(object::error("Levin tolerance must be positive"));
    }
    let af = |x| f.sample(&[x]);
    let gp = |x| phase.sample(&[x]);
    let dg = |x| derivative.sample(&[x]);
    let r = if adaptive {
        crate::kernels::oscillatory::levin_adaptive(&af, &gp, &dg, a, b, w, tol, points)
    } else {
        crate::kernels::oscillatory::levin(&af, &gp, &dg, a, b, w, points)
    };
    f.finish()?;
    phase.finish()?;
    derivative.finish()?;
    let r = r.map_err(object::error)?;
    Ok(Dynamic(object::record(vec![
        (
            "value",
            complex(num_complex::Complex64::new(r.re, r.im), d.clone())?,
        ),
        ("error", optional_error(r.error, d)?),
    ])))
}
object_binding!(levin_fixed,"num_analysis_levin_compiled","crate::kernels::oscillatory::levin","Actual native Levin collocation with independent physical amplitude,phase,and derivative programs. Absolute reference tolerance is checked for interface consistency; fixed collocation supplies no error estimate.",false,(lo:Input,hi:Input,omega:Input,tolerance:Input,points:usize,budget:usize,function:Input,phase:Input,derivative:Input)->Dynamic,|ctx:&CallCtx<'_>,lo:Input,hi:Input,omega:Input,tolerance:Input,points:usize,budget:usize,function:Input,phase:Input,derivative:Input|levin(ctx,lo,hi,omega,tolerance,points,budget,function,phase,derivative,false));
object_binding!(levin_adaptive,"num_analysis_levin_adaptive_compiled","crate::kernels::oscillatory::levin_adaptive","Actual native adaptive Levin collocation with positive absolute physical tolerance and bounded collocation/evaluation work.",false,(lo:Input,hi:Input,omega:Input,tolerance:Input,points:usize,budget:usize,function:Input,phase:Input,derivative:Input)->Dynamic,|ctx:&CallCtx<'_>,lo:Input,hi:Input,omega:Input,tolerance:Input,points:usize,budget:usize,function:Input,phase:Input,derivative:Input|levin(ctx,lo,hi,omega,tolerance,points,budget,function,phase,derivative,true));
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    object::install(
        i,
        &[
            &triangle::BINDING,
            &importance::BINDING,
            &filon_real::BINDING,
            &filon_complex::BINDING,
            &oscillatory::BINDING,
            &levin_fixed::BINDING,
            &levin_adaptive::BINDING,
        ],
        "numerical_vector",
    )
}
