//! Native tape and multiple-output execution with explicit heterogeneous scales.
use super::native_object_api::Dynamic;
use super::native_object_api::Handle;
use super::native_object_api::Input;
use super::native_object_api::NativeType;
use super::native_object_api::NewObject;
use super::native_object_api::object_binding;
use super::native_object_api::{
    self as api,
};
use crate::api::Session;
use crate::backend::CompiledMulti;
use crate::backend::TapeFunction;
use crate::backend::tape::Interval;
use crate::backend::tape::Scalar;
use crate::backend::tape::Tape;
use crate::backend::tape::{
    self,
};
use crate::backend::{
    self,
};
use crate::graph::RuleError;
use crate::graph::rule::Installer;
use crate::kernel::CallCtx;
use crate::kernel::OpError;
use crate::kernel::Value;
use crate::kernel::signature::DimExpr;
use crate::kernel::signature::TypeExpr;
use crate::kernel::value::TensorValue;
use std::sync::Arc;
const LIMIT: usize = 262_144;
#[derive(Clone)]
struct Scale {
    value: f64,
    dim: DimExpr,
}
#[derive(Clone)]
struct Frame {
    inputs: Vec<Scale>,
    params: Vec<Scale>,
    outputs: Vec<Scale>,
}
struct NativeTape {
    tape: Tape,
    frame: Frame,
}
struct Multi {
    native: Arc<dyn CompiledMulti>,
    frame: Frame,
    work_per_row: usize,
}
impl NativeType for NativeTape {
    const KIND: &'static str = "rssn/native/backend/tape/v1";
    const MAX_BYTES: usize = 4_194_304;

    fn validate(&self) -> Result<(), OpError> {
        if self.tape.insts.len() > 4096
            || self.tape.args.len() > 65536
            || self.tape.inputs != self.frame.inputs.len()
            || self.tape.params != self.frame.params.len()
            || self.tape.outputs.len() != self.frame.outputs.len()
        {
            return Err(api::error("native tape shape/work budget exceeded"));
        }
        self.tape.validate().map_err(api::error)
    }
}
impl NativeType for Multi {
    const KIND: &'static str = "rssn/native/backend/multiple/v1";
    const MAX_BYTES: usize = 4_194_304;

    fn validate(&self) -> Result<(), OpError> {
        if self.native.inputs() != self.frame.inputs.len()
            || self.native.params() != self.frame.params.len()
            || self.native.outputs() != self.frame.outputs.len()
        {
            Err(api::error("compiled multiple arity mismatch"))
        } else {
            Ok(())
        }
    }
}
fn vector(value: &Value) -> Result<&[Value], OpError> {
    let Value::Tensor(t) = value else {
        return Err(api::error("expected a vector"));
    };
    if t.shape != [t.data.len()] || t.data.len() > LIMIT {
        return Err(api::error("invalid native vector shape/budget"));
    }
    Ok(&t.data)
}
fn scale(
    ctx: &CallCtx<'_>,
    v: &Value,
) -> Result<Scale, OpError> {
    let (value, dim) = if matches!(v, Value::Opaque(_)) {
        super::native_backend::scalar_reference(ctx, v)?
    } else {
        api::reference(v)?
    };
    if value <= 0.0 || !value.is_finite() {
        return Err(api::error("physical scales must be positive"));
    }
    Ok(Scale { value, dim })
}
fn scales(
    ctx: &CallCtx<'_>,
    v: &Value,
) -> Result<Vec<Scale>, OpError> {
    vector(v)?.iter().map(|v| scale(ctx, v)).collect()
}
fn frame(
    ctx: &CallCtx<'_>,
    i: &Value,
    p: &Value,
    o: &Value,
) -> Result<Frame, OpError> {
    let frame = Frame {
        inputs: scales(ctx, i)?,
        params: scales(ctx, p)?,
        outputs: scales(ctx, o)?,
    };
    if frame.inputs.len() + frame.params.len() + frame.outputs.len() > 64
        || frame.outputs.is_empty()
    {
        return Err(api::error(
            "native program requires1..64 total coordinate/parameter/output scales",
        ));
    }
    Ok(frame)
}
fn domains(f: &Frame) -> Vec<DimExpr> {
    f.inputs
        .iter()
        .chain(&f.params)
        .chain(&f.outputs)
        .map(|s| s.dim.clone())
        .collect()
}
fn points(
    ctx: &CallCtx<'_>,
    v: &Value,
    scales: &[Scale],
) -> Result<Vec<f64>, OpError> {
    let values = vector(v)?;
    if values.len() != scales.len() {
        return Err(api::error("native program point arity mismatch"));
    }
    values
        .iter()
        .zip(scales)
        .map(|(v, s)| {
            let (value, dim) = if matches!(v, Value::Opaque(_)) {
                super::native_backend::scalar_value(ctx, v)?
            } else {
                (super::numeric(v, 0)?, api::dimension(v)?)
            };
            api::same_dimension(&dim, &s.dim)?;
            let value = value / s.value;
            if !value.is_finite() {
                return Err(api::error("normalized coordinate is nonfinite"));
            }
            Ok(value)
        })
        .collect()
}
fn quantity_any(
    value: f64,
    dim: DimExpr,
) -> Result<Value, OpError> {
    if value.is_nan() {
        return Err(api::error("native interval result is NaN"));
    }
    Ok(if dim == DimExpr::Dimless {
        Value::Float(value)
    } else {
        Value::Dyn(crate::kernel::value::DynValue {
            value: Box::new(Value::Float(value)),
            dim,
        })
    })
}
fn values(
    data: Vec<f64>,
    scales: &[Scale],
) -> Result<Dynamic, OpError> {
    if data.len() != scales.len() {
        return Err(api::error("native output count mismatch"));
    }
    Ok(Dynamic(Value::Tensor(TensorValue {
        shape: vec![data.len()],
        elem: api::dynamic_type(),
        data: data
            .into_iter()
            .zip(scales)
            .map(|(v, s)| api::quantity(api::scale_product(v, &[s.value], &[])?, s.dim.clone()))
            .collect::<Result<_, _>>()?,
    })))
}
fn ratio(
    e: &DimExpr,
    a: &DimExpr,
    b: Option<&DimExpr>,
) -> Result<DimExpr, OpError> {
    let mut dims = vec![e.clone(), DimExpr::Inverse(Box::new(a.clone()))];
    if let Some(b) = b {
        dims.push(DimExpr::Inverse(Box::new(b.clone())));
    }
    super::runtime_dimensions::expression(&super::runtime_dimensions::canonical(&DimExpr::Mul(
        dims,
    ))?)
}
fn derivative(
    data: Vec<f64>,
    f: &Frame,
    k: usize,
    second: bool,
) -> Result<Dynamic, OpError> {
    let output = f
        .outputs
        .get(k)
        .ok_or_else(|| api::error("output index out of range"))?;
    let n = f.inputs.len();
    if data.len() != if second { n * n } else { n } {
        return Err(api::error("native derivative shape mismatch"));
    }
    let mut result = Vec::with_capacity(data.len());
    for (i, v) in data.into_iter().enumerate() {
        let a = &f.inputs[if second { i / n } else { i }];
        let b = second.then(|| &f.inputs[i % n]);
        let factor = if let Some(b) = b {
            api::scale_product(v, &[output.value], &[a.value, b.value])?
        } else {
            api::scale_product(v, &[output.value], &[a.value])?
        };
        result.push(api::quantity(
            factor,
            ratio(&output.dim, &a.dim, b.map(|s| &s.dim))?,
        )?);
    }
    Ok(Dynamic(Value::Tensor(TensorValue {
        shape: if second {
            vec![n, n]
        } else {
            vec![n]
        },
        elem: api::dynamic_type(),
        data: result,
    })))
}
fn source(
    expressions: &[String],
    inputs: &[String],
    params: &[String],
) -> Result<(), OpError> {
    if expressions.is_empty()
        || expressions.len() > 64
        || inputs.len() + params.len() > 63
        || expressions
            .iter()
            .any(|s| s.is_empty() || s.len() > 16384 || !s.is_ascii())
        || expressions.iter().map(String::len).sum::<usize>() > 65536
    {
        return Err(api::error("native source size/arity budget exceeded"));
    }
    let mut names = std::collections::BTreeSet::new();
    for n in inputs.iter().chain(params) {
        if n.is_empty()
            || n.len() > 128
            || !n.chars().all(|c| c.is_ascii_alphanumeric() || c == '_')
            || !names.insert(n)
        {
            return Err(api::error(
                "native symbols must be unique ASCII identifiers",
            ));
        }
    }
    Ok(())
}
object_binding!(new_tape,"num_native_tape_lower","crate::backend::tape::lower","Lower elementary expressions in declared symbol/parameter order. Every coordinate, parameter and output has an explicit positive physical reference; mixed units stay in checked Value cells.",false,(expressions:Vec<String>,inputs:Vec<String>,params:Vec<String>,input_scales:Input,param_scales:Input,output_scales:Input)->NewObject<NativeTape>,|ctx,expressions,inputs,params,input_scales,param_scales,output_scales|{source(&expressions,&inputs,&params)?;let f=frame(ctx,&input_scales.0,&param_scales.0,&output_scales.0)?;if f.inputs.len()!=inputs.len()||f.params.len()!=params.len()||f.outputs.len()!=expressions.len(){return Err(api::error("program references do not match source arities"));}let session=Session::with_rules(&[crate::rules::elementary()]).map_err(api::error)?;let tape=session.with_native_compilation(&expressions,&inputs,&params,tape::lower).map_err(api::error)?;let dims=domains(&f);Ok(NewObject::new(NativeTape{tape,frame:f},dims))});
object_binding!(compile_multi,"num_native_backend_compile_multi","crate::backend::Backend::compile_multi","Compile multiple elementary expressions with physical coordinate, parameter and output references. Uses the actual selected native backend and bounded source work.",false,(expressions:Vec<String>,inputs:Vec<String>,params:Vec<String>,input_scales:Input,param_scales:Input,output_scales:Input)->NewObject<Multi>,|ctx,expressions,inputs,params,input_scales,param_scales,output_scales|{source(&expressions,&inputs,&params)?;let f=frame(ctx,&input_scales.0,&param_scales.0,&output_scales.0)?;if f.inputs.len()!=inputs.len()||f.params.len()!=params.len()||f.outputs.len()!=expressions.len(){return Err(api::error("program references do not match source arities"));}let session=Session::with_rules(&[crate::rules::elementary()]).map_err(api::error)?;let chosen=backend::current();if chosen.name()!="interpreter" {return Err(api::error("bounded multi-output compilation requires interpreter backend"));}let mut work_per_row=0usize;let native=session.with_native_compilation(&expressions,&inputs,&params,|g,r,i,p|{let tape=tape::lower(g,r,i,p)?;work_per_row=tape.insts.len()+tape.args.len();if tape.insts.len()>4096{return Err(backend::BackendError::Codegen("program exceeds4096instructions".into()));}chosen.compile_multi(g,r,i,p)}).map_err(api::error)?;let dims=domains(&f);Ok(NewObject::new(Multi{native:Arc::from(native),frame:f,work_per_row},dims))});
object_binding!(optimize,"num_native_tape_optimise","crate::backend::tape::Tape::optimise","Return the native optimized tape retaining its physical frame.",true,(a:Handle<NativeTape>)->NewObject<NativeTape>,|ctx,a|a.with(ctx,|t|Ok(NewObject::new(NativeTape{tape:t.tape.optimise(),frame:t.frame.clone()},domains(&t.frame)))));
object_binding!(hash,"num_native_tape_hash","crate::backend::tape::Tape::structural_hash","Native deterministic structural tape hash.",true,(a:Handle<NativeTape>)->u64,|ctx,a|a.with(ctx,|t|Ok(t.tape.structural_hash())));
object_binding!(validate,"num_native_tape_validate","crate::backend::tape::Tape::validate","Validate native tape register and output invariants.",true,(a:Handle<NativeTape>)->bool,|ctx,a|a.with(ctx,|t|{t.tape.validate().map_err(api::error)?;Ok(true)}));
object_binding!(eval,"num_native_tape_eval","crate::backend::tape::Tape::eval","Evaluate native tape using measured heterogeneous coordinates/parameters; returns checked measured Value cells.",false,(a:Handle<NativeTape>,x:Input,p:Input)->Dynamic,|ctx,a,x,p|a.with(ctx,|t|{let x=points(ctx,&x.0,&t.frame.inputs)?;let p=points(ctx,&p.0,&t.frame.params)?;let mut out=vec![0.0;t.frame.outputs.len()];t.tape.eval(&x,&p,&mut out);values(out,&t.frame.outputs)}));
object_binding!(run,"num_native_tape_run_float","crate::backend::tape::Tape::run","Explicit Float specialization of generic native tape execution, preserving physical coordinate/output references.",false,(a:Handle<NativeTape>,x:Input,p:Input)->Dynamic,|ctx,a,x,p|a.with(ctx,|t|{let x=points(ctx,&x.0,&t.frame.inputs)?;let p=points(ctx,&p.0,&t.frame.params)?;let mut out=vec![0.0;t.frame.outputs.len()];t.tape.run::<f64>(&x,&p,&mut Vec::new(),&mut out);values(out,&t.frame.outputs)}));
object_binding!(gradient,"num_native_tape_gradient","crate::backend::tape::Tape::gradient","Native forward tape gradient; each returned Value cell has its own output/input unit.",false,(a:Handle<NativeTape>,k:usize,x:Input,p:Input)->Dynamic,|ctx,a,k,x,p|a.with(ctx,|t|{if k>=t.frame.outputs.len(){return Err(api::error("output index out of range"));}let x=points(ctx,&x.0,&t.frame.inputs)?;let p=points(ctx,&p.0,&t.frame.params)?;derivative(t.tape.gradient(k,&x,&p),&t.frame,k,false)}));
object_binding!(reverse,"num_native_tape_gradient_reverse","crate::backend::tape::Tape::gradient_reverse","Native reverse tape gradient with heterogeneous coordinate units.",false,(a:Handle<NativeTape>,k:usize,x:Input,p:Input)->Dynamic,|ctx,a,k,x,p|a.with(ctx,|t|{if k>=t.frame.outputs.len(){return Err(api::error("output index out of range"));}let x=points(ctx,&x.0,&t.frame.inputs)?;let p=points(ctx,&p.0,&t.frame.params)?;derivative(t.tape.gradient_reverse(k,&x,&p),&t.frame,k,false)}));
object_binding!(hessian,"num_native_tape_hessian","crate::backend::tape::Tape::hessian","Native hyperdual tape Hessian; output(i,j) has units E/(D_i D_j).",false,(a:Handle<NativeTape>,k:usize,x:Input,p:Input)->Dynamic,|ctx,a,k,x,p|a.with(ctx,|t|{if k>=t.frame.outputs.len(){return Err(api::error("output index out of range"));}
if t.frame.inputs.len()>32{return Err(api::error("Hessian work budget exceeds32coordinates"));}let x=points(ctx,&x.0,&t.frame.inputs)?;let p=points(ctx,&p.0,&t.frame.params)?;derivative(t.tape.hessian(k,&x,&p),&t.frame,k,true)}));
object_binding!(hvp,"num_native_tape_hessian_vector_product","crate::backend::tape::Tape::hessian_vector_product","Native Hessian-vector product with a measured coordinate displacement direction; result units E/D_i.",false,(a:Handle<NativeTape>,k:usize,x:Input,p:Input,direction:Input)->Dynamic,|ctx,a,k,x,p,direction|a.with(ctx,|t|{if k>=t.frame.outputs.len(){return Err(api::error("output index out of range"));}let x=points(ctx,&x.0,&t.frame.inputs)?;let p=points(ctx,&p.0,&t.frame.params)?;let d=points(ctx,&direction.0,&t.frame.inputs)?;derivative(t.tape.hessian_vector_product(k,&x,&p,&d),&t.frame,k,false)}));
object_binding!(multi_from_tape,"num_native_tape_function_from_tape","crate::backend::TapeFunction::from_tape","Wrap a real native tape in the multiple-output executor.",true,(a:Handle<NativeTape>)->NewObject<Multi>,|ctx,a|a.with(ctx,|t|Ok(NewObject::new(Multi{native:Arc::new(TapeFunction::from_tape(t.tape.clone())),frame:t.frame.clone(),work_per_row:t.tape.insts.len()+t.tape.args.len()},domains(&t.frame)))));
object_binding!(multi_eval,"num_native_compiled_multi_eval","crate::backend::CompiledMulti::eval","Evaluate a native multiple-output program with physical references.",false,(a:Handle<Multi>,x:Input,p:Input)->Dynamic,|ctx,a,x,p|a.with(ctx,|t|{let x=points(ctx,&x.0,&t.frame.inputs)?;let p=points(ctx,&p.0,&t.frame.params)?;let mut out=vec![0.0;t.frame.outputs.len()];t.native.eval(&x,&p,&mut out);values(out,&t.frame.outputs)}));
object_binding!(multi_inputs,"num_native_compiled_multi_inputs","crate::backend::CompiledMulti::inputs","Native program inputs count.",true,(a:Handle<Multi>)->usize,|ctx,a|a.with(ctx,|t|Ok(t.native.inputs())));
object_binding!(multi_params,"num_native_compiled_multi_params","crate::backend::CompiledMulti::params","Native program params count.",true,(a:Handle<Multi>)->usize,|ctx,a|a.with(ctx,|t|Ok(t.native.params())));
object_binding!(multi_outputs,"num_native_compiled_multi_outputs","crate::backend::CompiledMulti::outputs","Native program outputs count.",true,(a:Handle<Multi>)->usize,|ctx,a|a.with(ctx,|t|Ok(t.native.outputs())));

fn physical_columns(
    ctx: &CallCtx<'_>,
    v: &Value,
    scales: &[Scale],
) -> Result<Vec<Vec<f64>>, OpError> {
    let columns = vector(v)?;
    if columns.len() != scales.len() {
        return Err(api::error("native column count mismatch"));
    }
    let mut rows = None;
    let mut cells = 0usize;
    columns
        .iter()
        .zip(scales)
        .map(|(column, s)| {
            let n = vector(column)?.len();
            if rows.is_some_and(|r| r != n) {
                return Err(api::error("native columns have unequal lengths"));
            }
            rows = Some(n);
            cells = cells
                .checked_add(n)
                .filter(|n| *n <= LIMIT)
                .ok_or_else(|| api::error("column data budget exceeded"))?;
            points(ctx, column, &vec![s.clone(); n])
        })
        .collect()
}
fn column_values(
    data: Vec<Vec<f64>>,
    scales: &[Scale],
) -> Result<Dynamic, OpError> {
    if data.len() != scales.len() || data.iter().map(Vec::len).sum::<usize>() > LIMIT {
        return Err(api::error("native output column budget exceeded"));
    }
    let data = data
        .into_iter()
        .zip(scales)
        .map(|(v, s)| {
            let n = v.len();
            Ok(Value::Tensor(TensorValue {
                shape: vec![n],
                elem: TypeExpr::Float.with_dim(s.dim.clone()),
                data: v
                    .into_iter()
                    .map(|v| api::quantity(api::scale_product(v, &[s.value], &[])?, s.dim.clone()))
                    .collect::<Result<_, _>>()?,
            }))
        })
        .collect::<Result<Vec<_>, OpError>>()?;
    Ok(Dynamic(Value::Tensor(TensorValue {
        shape: vec![data.len()],
        elem: api::dynamic_type(),
        data,
    })))
}
fn coherent_value(
    ctx: &CallCtx<'_>,
    v: &Value,
) -> Result<(f64, DimExpr), OpError> {
    if matches!(v, Value::Opaque(_)) {
        super::native_backend::scalar_value(ctx, v)
    } else {
        Ok((super::numeric(v, 0)?, api::dimension(v)?))
    }
}
fn bound_value(v: &Value) -> Result<f64, OpError> {
    match v {
        | Value::Float(x) if x.is_infinite() => Ok(*x),
        | Value::Dyn(d) => bound_value(&d.value),
        | Value::Quantity(q) => {
            let factor = super::numeric(
                &Value::Quantity(crate::kernel::value::QuantityValue {
                    value: Box::new(Value::Float(1.0)),
                    unit: q.unit.clone(),
                }),
                0,
            )?;
            let x = bound_value(&q.value)?;
            if x.is_infinite() {
                Ok(x * factor)
            } else {
                Ok(super::numeric(v, 0)?)
            }
        },
        | _ => super::numeric(v, 0),
    }
}
fn intervals(
    v: &Value,
    scales: &[Scale],
) -> Result<Vec<Interval>, OpError> {
    let data = vector(v)?;
    if data.len() != scales.len() {
        return Err(api::error("native interval count mismatch"));
    }
    data.iter()
        .zip(scales)
        .map(|(v, s)| {
            let lo = api::field(v, "lo")?;
            let hi = api::field(v, "hi")?;
            api::same_dimension(&api::dimension(lo)?, &s.dim)?;
            api::same_dimension(&api::dimension(hi)?, &s.dim)?;
            let lo = bound_value(lo)?;
            let hi = bound_value(hi)?;
            if lo > hi || lo.is_nan() || hi.is_nan() {
                return Err(api::error("native interval bounds invalid"));
            }
            Ok(<Interval as Scalar>::div(
                Interval { lo, hi },
                Interval::point(s.value),
            ))
        })
        .collect()
}
struct FunctionObject {
    native: TapeFunction,
    frame: Frame,
}
impl NativeType for FunctionObject {
    const KIND: &'static str = "rssn/native/backend/tape-function/v1";
    const MAX_BYTES: usize = 4_194_304;

    fn validate(&self) -> Result<(), OpError> {
        self.native.tape().validate().map_err(api::error)
    }
}
object_binding!(tape_eval_columns,"num_native_tape_eval_columns","crate::backend::tape::Tape::eval_columns","Native column tape evaluation; each coordinate/output column has independent retained physical references.",false,(a:Handle<NativeTape>,rows:usize,columns:Input,p:Input)->Dynamic,|ctx,a,rows,columns,p|{a.with(ctx,|t|{let cols=physical_columns(ctx,&columns.0,&t.frame.inputs)?;if cols.first().is_some_and(|c|c.len()!=rows){return Err(api::error("explicit batch row extent disagrees with columns"));}
if rows.checked_mul(t.tape.insts.len()+t.tape.args.len()).is_none_or(|n|n>100_000_000)||rows.checked_mul(t.frame.outputs.len()).is_none_or(|n|n>LIMIT){return Err(api::error("native column work/allocation budget exceeded"));}let params=points(ctx,&p.0,&t.frame.params)?;let slices=cols.iter().map(Vec::as_slice).collect::<Vec<_>>();let mut output=vec![vec![0.0;rows];t.frame.outputs.len()];let mut out=output.iter_mut().map(Vec::as_mut_slice).collect::<Vec<_>>();t.tape.eval_columns(&slices,&params,&mut out);column_values(output,&t.frame.outputs)})});
object_binding!(compiled_multi_eval_batch,"num_native_compiled_multi_eval_batch","crate::backend::CompiledMulti::eval_batch","Actual multiple-output batch executor with per-column physical domains and bounded rows.",false,(a:Handle<Multi>,rows:usize,columns:Input,p:Input)->Dynamic,|ctx,a,rows,columns,p|{a.with(ctx,|t|{let cols=physical_columns(ctx,&columns.0,&t.frame.inputs)?;if cols.first().is_some_and(|c|c.len()!=rows){return Err(api::error("explicit batch row extent disagrees with columns"));}
if rows.checked_mul(t.frame.outputs.len()).is_none_or(|n|n>LIMIT)||rows>16384||rows.checked_mul(t.work_per_row).is_none_or(|n|n>100_000_000){return Err(api::error("native compiled batch budget exceeded"));}let params=points(ctx,&p.0,&t.frame.params)?;let slices=cols.iter().map(Vec::as_slice).collect::<Vec<_>>();let mut output=vec![vec![0.0;rows];t.frame.outputs.len()];let mut out=output.iter_mut().map(Vec::as_mut_slice).collect::<Vec<_>>();t.native.eval_batch(&slices,&params,&mut out);column_values(output,&t.frame.outputs)})});
object_binding!(tape_vjp,"num_native_tape_vjp","crate::backend::tape::Tape::vjp","Native VJP; weights have functional/output units and each answer has functional/input units. Positive functional reference specifies the normalized objective measure.",false,(a:Handle<NativeTape>,x:Input,p:Input,weights:Input,functional_scale:Input)->Dynamic,|ctx,a,x,p,weights,functional_scale|{a.with(ctx,|t|{let x=points(ctx,&x.0,&t.frame.inputs)?;let p=points(ctx,&p.0,&t.frame.params)?;let functional=scale(ctx,&functional_scale.0)?;let w=vector(&weights.0)?;if w.len()!=t.frame.outputs.len(){return Err(api::error("VJP weight count mismatch"));}let weights=w.iter().zip(&t.frame.outputs).map(|(v,o)|{let(value,dim)=coherent_value(ctx,v)?;api::same_dimension(&dim,&ratio(&functional.dim,&o.dim,None)?)?;api::scale_product(value,&[o.value],&[functional.value])}).collect::<Result<Vec<_>,_>>()?;let mut frame=t.frame.clone();frame.outputs=vec![functional];derivative(t.tape.vjp(&x,&p,&weights),&frame,0,false)})});
object_binding!(tape_enclose,"num_native_tape_enclose","crate::backend::tape::Tape::enclose","Native outward interval enclosure; measured bounds and coordinate references are checked. Infinite endpoints retain their sound interval meaning.",false,(a:Handle<NativeTape>,x:Input,p:Input)->Dynamic,|ctx,a,x,p|{a.with(ctx,|t|{let x=intervals(&x.0,&t.frame.inputs)?;let p=intervals(&p.0,&t.frame.params)?;let bounds=t.tape.enclose(&x,&p);let data=bounds.into_iter().zip(&t.frame.outputs).map(|(b,s)|{let b=<Interval as Scalar>::mul(b,Interval::point(s.value));Ok(api::record(vec![("lo",quantity_any(b.lo,s.dim.clone())?),("hi",quantity_any(b.hi,s.dim.clone())?)]))}).collect::<Result<Vec<_>,OpError>>()?;Ok(Dynamic(Value::Tensor(TensorValue{shape:vec![data.len()],elem:api::dynamic_type(),data})))})});
object_binding!(tape_function_object,"num_native_tape_function_object","crate::backend::TapeFunction::from_tape","Retain the exact native TapeFunction type, including its physical frame.",false,(a:Handle<NativeTape>)->NewObject<FunctionObject>,|ctx,a|{a.with(ctx,|t|Ok(NewObject::new(FunctionObject{native:TapeFunction::from_tape(t.tape.clone()),frame:t.frame.clone()},domains(&t.frame))))});
object_binding!(tape_function_tape,"num_native_tape_function_tape","crate::backend::TapeFunction::tape","Recover the defining native tape from a retained TapeFunction without serializing internal pointers.",false,(a:Handle<FunctionObject>)->NewObject<NativeTape>,|ctx,a|{a.with(ctx,|t|Ok(NewObject::new(NativeTape{tape:t.native.tape().clone(),frame:t.frame.clone()},domains(&t.frame))))});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &new_tape::BINDING,
            &compile_multi::BINDING,
            &optimize::BINDING,
            &hash::BINDING,
            &validate::BINDING,
            &eval::BINDING,
            &run::BINDING,
            &gradient::BINDING,
            &reverse::BINDING,
            &hessian::BINDING,
            &hvp::BINDING,
            &multi_from_tape::BINDING,
            &multi_eval::BINDING,
            &multi_inputs::BINDING,
            &multi_params::BINDING,
            &multi_outputs::BINDING,
            &tape_eval_columns::BINDING,
            &compiled_multi_eval_batch::BINDING,
            &tape_vjp::BINDING,
            &tape_enclose::BINDING,
            &tape_function_object::BINDING,
            &tape_function_tape::BINDING,
        ],
        "numerical_vector",
    )
}

#[derive(Clone)]
pub(in crate::rules) struct TapeSnapshot {
    pub tape: Tape,
    pub input_scales: Vec<f64>,
    pub input_dimensions: Vec<DimExpr>,
    pub parameter_scales: Vec<f64>,
    pub parameter_dimensions: Vec<DimExpr>,
    pub output_scales: Vec<f64>,
    pub output_dimensions: Vec<DimExpr>,
}
pub(in crate::rules) fn tape_snapshot(
    ctx: &CallCtx<'_>,
    handle: &Value,
) -> Result<TapeSnapshot, OpError> {
    api::get::<NativeTape, _>(ctx, handle, NativeTape::KIND, |t| {
        Ok(TapeSnapshot {
            tape: t.tape.clone(),
            input_scales: t.frame.inputs.iter().map(|s| s.value).collect(),
            input_dimensions: t.frame.inputs.iter().map(|s| s.dim.clone()).collect(),
            parameter_scales: t.frame.params.iter().map(|s| s.value).collect(),
            parameter_dimensions: t.frame.params.iter().map(|s| s.dim.clone()).collect(),
            output_scales: t.frame.outputs.iter().map(|s| s.value).collect(),
            output_dimensions: t.frame.outputs.iter().map(|s| s.dim.clone()).collect(),
        })
    })
}

pub(super) fn normalized_snapshot(tape: Tape) -> TapeSnapshot {
    TapeSnapshot {
        input_scales: vec![1.; tape.inputs],
        input_dimensions: vec![DimExpr::Dimless; tape.inputs],
        parameter_scales: vec![1.; tape.params],
        parameter_dimensions: vec![DimExpr::Dimless; tape.params],
        output_scales: vec![1.; tape.outputs.len()],
        output_dimensions: vec![DimExpr::Dimless; tape.outputs.len()],
        tape,
    }
}
pub(super) fn tape_value(
    ctx: &CallCtx<'_>,
    s: TapeSnapshot,
) -> Result<Value, OpError> {
    use super::native_object_api::NativeOutput;
    let frame = Frame {
        inputs: s
            .input_scales
            .into_iter()
            .zip(s.input_dimensions)
            .map(|(value, dim)| Scale { value, dim })
            .collect(),
        params: s
            .parameter_scales
            .into_iter()
            .zip(s.parameter_dimensions)
            .map(|(value, dim)| Scale { value, dim })
            .collect(),
        outputs: s
            .output_scales
            .into_iter()
            .zip(s.output_dimensions)
            .map(|(value, dim)| Scale { value, dim })
            .collect(),
    };
    let d = domains(&frame);
    NewObject::new(NativeTape { tape: s.tape, frame }, d).encode(ctx)
}

#[cfg(feature = "jit")]
pub(super) fn snapshot_dimensions(s: &TapeSnapshot) -> Vec<DimExpr> {
    s.input_dimensions
        .iter()
        .chain(&s.parameter_dimensions)
        .chain(&s.output_dimensions)
        .cloned()
        .collect()
}
#[cfg(feature = "jit")]
pub(super) fn empty_arguments() -> Value {
    Value::Tensor(TensorValue {
        shape: vec![0],
        elem: api::dynamic_type(),
        data: vec![],
    })
}
#[cfg(feature = "jit")]
pub(super) fn snapshot_arguments(
    ctx: &CallCtx<'_>,
    s: &TapeSnapshot,
    x: &Value,
    p: &Value,
) -> Result<(Vec<f64>, Vec<f64>), OpError> {
    let scales = |v: &[f64], d: &[DimExpr]| {
        v.iter()
            .copied()
            .zip(d.iter().cloned())
            .map(|(value, dim)| Scale { value, dim })
            .collect::<Vec<_>>()
    };
    Ok((
        points(ctx, x, &scales(&s.input_scales, &s.input_dimensions))?,
        points(
            ctx,
            p,
            &scales(&s.parameter_scales, &s.parameter_dimensions),
        )?,
    ))
}
#[cfg(feature = "jit")]
pub(super) fn snapshot_outputs(
    v: Vec<f64>,
    s: &TapeSnapshot,
) -> Result<Dynamic, OpError> {
    values(
        v,
        &s.output_scales
            .iter()
            .copied()
            .zip(s.output_dimensions.iter().cloned())
            .map(|(value, dim)| Scale { value, dim })
            .collect::<Vec<_>>(),
    )
}
