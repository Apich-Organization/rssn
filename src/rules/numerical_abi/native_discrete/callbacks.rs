#![allow(clippy::wildcard_imports)] // Typed host callbacks are invoked serially with preserved dimensions and first failure.
use super::*;
use crate::kernel::value::CallbackRef;
use crate::kernels::differential_geometry as dg;
use crate::kernels::indefinite_sum as sums;
use crate::kernels::robust;
use crate::kernels::stats_ext;
use std::cell::Cell;
use std::cell::RefCell;
fn scalar_type(t: &TypeExpr) -> Result<DimExpr, OpError> {
    match t {
        | TypeExpr::Float => Ok(DimExpr::Dimless),
        | TypeExpr::Dim(inner, d) if **inner == TypeExpr::Float => canonical(d),
        | _ => {
            Err(error(
                "callback scalar must be Float with a concrete dimension",
            ))
        },
    }
}
fn tensor_dimension(
    t: &TypeExpr,
    rank: u64,
) -> Result<DimExpr, OpError> {
    match t {
        | TypeExpr::Tensor { rank: r, elem, .. } if *r == NatExpr::Lit(rank) => scalar_type(elem),
        | _ => Err(error("callback tensor rank/scalar schema mismatch")),
    }
}
struct ScalarFn {
    handle: CallbackRef,
    input: DimExpr,
    output: DimExpr,
}
struct StatisticFn {
    handle: CallbackRef,
    input: DimExpr,
    output: DimExpr,
}
struct MetricFn {
    handle: CallbackRef,
    input: DimExpr,
    output: DimExpr,
}

impl MetricFn {
    fn normalized_point(
        &self,
        point: &Input,
        coordinate_reference: &Input,
        metric_reference: &Input,
    ) -> Result<(Vec<f64>, f64, f64), OpError> {
        let (coordinate, input) = api::reference(&coordinate_reference.0)?;
        let (metric, output) = api::reference(&metric_reference.0)?;
        api::same_dimension(&input, &self.input)?;
        api::same_dimension(&output, &self.output)?;
        api::same_dimension(&api::dimension(&point.0)?, &input)?;
        let mut point = <Vec<f64> as Arg>::decode(&point.0, 0)?;
        if !(1..=6).contains(&point.len()) {
            return Err(error("metric dimension must be1..6"));
        }
        for coordinate_value in &mut point {
            *coordinate_value /= coordinate;
            if !coordinate_value.is_finite() {
                return Err(error("normalized metric coordinates not finite"));
            }
        }
        Ok((point, coordinate, metric))
    }
}
struct CdfFn(ScalarFn);
fn fn_type(
    input: TypeExpr,
    output: TypeExpr,
) -> TypeExpr {
    TypeExpr::Function(vec![input], Box::new(output))
}
fn var(n: &str) -> TypeExpr {
    TypeExpr::Float.with_dim(DimExpr::Var(n.into()))
}
fn signature<'a>(
    v: &Value,
    c: &'a CallCtx<'_>,
) -> Result<(CallbackRef, &'a crate::kernel::signature::Signature), OpError> {
    let Value::Function(h) = v else {
        return Err(error("expected a typed callback"));
    };
    let s = c.callback_signature(*h)?;
    if s.params.len() != 1 {
        return Err(error("callback must have one argument"));
    }
    Ok((*h, s))
}
impl NativeArg for ScalarFn {
    fn ty() -> TypeExpr {
        fn_type(var("X"), var("Y"))
    }

    fn decode(
        v: &Value,
        _: usize,
        c: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        let (h, s) = signature(v, c)?;
        Ok(Self {
            handle: h,
            input: scalar_type(&s.params[0])?,
            output: scalar_type(&s.ret)?,
        })
    }
}
impl NativeArg for CdfFn {
    fn ty() -> TypeExpr {
        fn_type(var("X"), TypeExpr::Float)
    }

    fn decode(
        v: &Value,
        i: usize,
        c: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        let f = <ScalarFn as NativeArg>::decode(v, i, c)?;
        api::same_dimension(&f.output, &DimExpr::Dimless)?;
        Ok(Self(f))
    }
}
impl NativeArg for StatisticFn {
    fn ty() -> TypeExpr {
        fn_type(api::tensor_type(1, var("X")), var("Y"))
    }

    fn decode(
        v: &Value,
        _: usize,
        c: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        let (h, s) = signature(v, c)?;
        Ok(Self {
            handle: h,
            input: tensor_dimension(&s.params[0], 1)?,
            output: scalar_type(&s.ret)?,
        })
    }
}
impl NativeArg for MetricFn {
    fn ty() -> TypeExpr {
        fn_type(api::tensor_type(1, var("X")), api::tensor_type(2, var("Y")))
    }

    fn decode(
        v: &Value,
        _: usize,
        c: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        let (h, s) = signature(v, c)?;
        Ok(Self {
            handle: h,
            input: tensor_dimension(&s.params[0], 1)?,
            output: tensor_dimension(&s.ret, 2)?,
        })
    }
}
struct Gate<'a, 'b> {
    ctx: &'a CallCtx<'b>,
    remaining: Cell<usize>,
    failure: RefCell<Option<OpError>>,
}
impl<'a, 'b> Gate<'a, 'b> {
    fn new(
        ctx: &'a CallCtx<'b>,
        budget: usize,
    ) -> Result<Self, OpError> {
        if budget == 0 || budget > 1_000_000 {
            return Err(error("callback budget must be1..1000000"));
        }
        Ok(Self {
            ctx,
            remaining: Cell::new(budget),
            failure: RefCell::new(None),
        })
    }

    fn call<T>(
        &self,
        f: impl FnOnce() -> Result<T, OpError>,
    ) -> Result<T, ()> {
        if self.failure.borrow().is_some() {
            return Err(());
        }
        let result = self.remaining.get().checked_sub(1).map_or_else(
            || Err(error("native callback budget exhausted")),
            |n| {
                self.remaining.set(n);
                f()
            },
        );
        match result {
            | Ok(v) => Ok(v),
            | Err(e) => {
                *self.failure.borrow_mut() = Some(e);
                Err(())
            },
        }
    }

    fn done(&self) -> Result<(), OpError> {
        self.failure.borrow_mut().take().map_or(Ok(()), Err)
    }

    fn scalar(
        &self,
        f: &ScalarFn,
        x: f64,
    ) -> f64 {
        self.call(|| {
            let a = api::quantity(x, f.input.clone())?;
            let out = self.ctx.invoke_callback(f.handle, &[a])?;
            api::same_dimension(&api::dimension(&out)?, &f.output)?;
            numeric(&out, 0)
        })
        .unwrap_or(f64::NAN)
    }

    fn statistic(
        &self,
        f: &StatisticFn,
        x: &[f64],
    ) -> f64 {
        self.call(|| {
            let a = quantities(x.to_vec(), &f.input)?.0;
            let out = self.ctx.invoke_callback(f.handle, &[a])?;
            api::same_dimension(&api::dimension(&out)?, &f.output)?;
            numeric(&out, 0)
        })
        .unwrap_or(f64::NAN)
    }

    fn metric(
        &self,
        f: &MetricFn,
        x: &[f64],
        coordinate: f64,
        reference: f64,
    ) -> Vec<Vec<f64>> {
        let n = x.len();
        self.call(|| {
            let a = quantities(x.iter().map(|x| x * coordinate).collect(), &f.input)?.0;
            let out = self.ctx.invoke_callback(f.handle, &[a])?;
            api::same_dimension(&api::dimension(&out)?, &f.output)?;
            let rows = <Vec<Vec<f64>> as Arg>::decode(&out, 0)?;
            if rows.len() != n || rows.iter().any(|r| r.len() != n) {
                return Err(error("metric callback must return n-by-n values"));
            }
            if rows.iter().enumerate().any(|(i, r)| {
                r.iter().enumerate().any(|(j, x)| {
                    (x - rows[j][i]).abs() > 1e-10 * x.abs().max(rows[j][i].abs()).max(reference)
                })
            }) {
                return Err(error("metric callback must be symmetric"));
            }
            let rows = rows
                .into_iter()
                .map(|r| r.into_iter().map(|v| v / reference).collect::<Vec<_>>())
                .collect::<Vec<_>>();
            if rows.iter().flatten().any(|x| !x.is_finite()) {
                return Err(error("normalized metric not finite"));
            }
            Ok(rows)
        })
        .unwrap_or_else(|()| vec![vec![f64::NAN; n]; n])
    }
}
fn sum_config(
    v: &Input,
    coordinate: f64,
    d: &DimExpr,
) -> Result<sums::IndefiniteSumConfig, OpError> {
    let h = api::field(&v.0, "h")?;
    let step = api::field(&v.0, "step")?;
    api::same_dimension(&api::dimension(h)?, d)?;
    api::same_dimension(&api::dimension(step)?, d)?;
    let o = sums::IndefiniteSumConfig {
        h: numeric(h, 0)? / coordinate,
        step: numeric(step, 0)? / coordinate,
        tolerance: numeric(api::field(&v.0, "tolerance")?, 0)?,
        shift: <usize as Arg>::decode(api::field(&v.0, "shift")?, 0)?,
        taylor_terms: <usize as Arg>::decode(api::field(&v.0, "taylor_terms")?, 0)?,
    };
    if o.step == 0.0
        || !o.step.is_finite()
        || !o.h.is_finite()
        || o.tolerance <= 0.0
        || o.shift > 10000
        || o.taylor_terms == 0
        || o.taylor_terms > 32
    {
        return Err(error(
            "invalid sum configuration or shift/Taylor work bound",
        ));
    }
    api::normalized(api::field(&v.0, "tolerance")?)?;
    Ok(o)
}

struct RecurrenceInput {
    point: f64,
    coordinate_reference: f64,
    value_reference: f64,
    config: sums::IndefiniteSumConfig,
}

fn normalize_recurrence(
    function: &ScalarFn,
    point: &Input,
    config: &Input,
    coordinate_reference: &Input,
    value_reference: &Input,
) -> Result<RecurrenceInput, OpError> {
    let (coordinate_reference, coordinate_dim) = api::reference(&coordinate_reference.0)?;
    let (value_reference, value_dim) = api::reference(&value_reference.0)?;
    api::same_dimension(&coordinate_dim, &function.input)?;
    api::same_dimension(&value_dim, &function.output)?;
    let (point, point_dim) = scalar(point)?;
    api::same_dimension(&coordinate_dim, &point_dim)?;
    let config = sum_config(config, coordinate_reference, &coordinate_dim)?;
    let point = point / coordinate_reference;
    if !point.is_finite() || ((point - config.h) / config.step).abs() > 100_000.0 {
        return Err(error("indefinite recurrence step budget exceeds 100000"));
    }
    Ok(RecurrenceInput {
        point,
        coordinate_reference,
        value_reference,
        config,
    })
}

fn tensor3(v: Vec<Vec<Vec<f64>>>) -> Result<Value, OpError> {
    let n = v.len();
    let values = v
        .into_iter()
        .flatten()
        .flatten()
        .map(Output::encode)
        .collect::<Result<Vec<_>, _>>()?;
    if n.checked_pow(3) != Some(values.len()) {
        return Err(error("native tensor3 shape mismatch"));
    }
    Ok(Value::Tensor(TensorValue {
        shape: vec![n; 3],
        data: values,
        elem: TypeExpr::Float,
    }))
}
fn tensor4(v: Vec<Vec<Vec<Vec<f64>>>>) -> Result<Value, OpError> {
    let n = v.len();
    let values = v
        .into_iter()
        .flatten()
        .flatten()
        .flatten()
        .map(Output::encode)
        .collect::<Result<Vec<_>, _>>()?;
    if n.checked_pow(4) != Some(values.len()) {
        return Err(error("native tensor4 shape mismatch"));
    }
    Ok(Value::Tensor(TensorValue {
        shape: vec![n; 4],
        data: values,
        elem: TypeExpr::Float,
    }))
}
fn scaled_tensor(
    value: Value,
    scale: f64,
    d: &DimExpr,
) -> Result<Dynamic, OpError> {
    restore_tensor(value, d, |value| Ok(value * scale))
}

fn restore_tensor(
    mut value: Value,
    d: &DimExpr,
    restore: impl Fn(f64) -> Result<f64, OpError>,
) -> Result<Dynamic, OpError> {
    let Value::Tensor(t) = &mut value else {
        return Err(error("expected tensor output"));
    };
    t.elem = TypeExpr::Float.with_dim(d.clone());
    for v in &mut t.data {
        *v = api::quantity(restore(numeric(v, 0)?)?, d.clone())?;
    }
    Ok(Dynamic(value))
}

struct RicciResult {
    tensor: Dynamic,
    scalar: Dynamic,
}

impl api::NativeOutput for RicciResult {
    fn ty() -> TypeExpr {
        TypeExpr::Record(vec![
            ("tensor".into(), api::dynamic_type()),
            ("scalar".into(), api::dynamic_type()),
        ])
    }

    fn encode(
        self,
        _: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        Ok(api::record(vec![
            ("tensor", self.tensor.0),
            ("scalar", self.scalar.0),
        ]))
    }
}

object_binding!(
    geometry_ricci_curvature,
    "num_native_geometry_ricci_curvature",
    "crate::kernels::differential_geometry::ricci_curvature",
    "Computes the covariant Ricci tensor and scalar curvature with one finite-difference metric stencil and factorization. Coordinates use a common dimension in one to six dimensions. With coordinate dimension X and metric dimension G, tensor entries have X^-2 and the scalar has X^-2/G. Use positive coordinate and metric references and supply a callback evaluation budget.",
    false,
    (
        metric: MetricFn,
        point: Input,
        coordinate_reference: Input,
        metric_reference: Input,
        callback_budget: usize,
    ) -> RicciResult,
    |ctx, metric, point, coordinate_reference, metric_reference, callback_budget| {
        let (x, coordinate, reference) =
            metric.normalized_point(&point, &coordinate_reference, &metric_reference)?;
        let gate = Gate::new(ctx, callback_budget)?;
        let curvature = dg::ricci_curvature(|x| gate.metric(&metric, x, coordinate, reference), &x);
        gate.done()?;
        let curvature = curvature.map_err(error)?;
        let tensor_dimension = metric.input.pow(Exp::Rat(Rational::int(-2)));
        let scalar_dimension = canonical(&DimExpr::Mul(vec![
            tensor_dimension.clone(),
            DimExpr::Inverse(Box::new(metric.output)),
        ]))?;
        let tensor = restore_tensor(
            Output::encode(curvature.tensor)?,
            &canonical(&tensor_dimension)?,
            |value| api::scale_product(value, &[], &[coordinate, coordinate]),
        )?;
        let scalar = api::scale_product(curvature.scalar, &[], &[coordinate, coordinate, reference])?;
        Ok(RicciResult {
            tensor,
            scalar: Dynamic(api::quantity(scalar, scalar_dimension)?),
        })
    }
);
object_binding!(stats_ks_test_one_sample,"num_native_stats_ks_test_one_sample","crate::kernels::stats_ext::ks_test_one_sample","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(sample:Input,cdf:CdfFn,callback_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,sample:Input,cdf:CdfFn,callback_budget:usize|{api::same_dimension(&api::dimension(&sample.0)?,&cdf.0.input)?;let x=<Vec<f64> as Arg>::decode(&sample.0,0)?;if x.is_empty(){return Err(error("KS sample is empty"));}limit(x.len(),262144)?;let gate=Gate::new(ctx,callback_budget)?;let out=stats_ext::ks_test_one_sample(&x,&|x|gate.scalar(&cdf.0,x));gate.done()?;let out=out.map_err(error)?;Ok(Dynamic(api::record(vec![("statistic",Output::encode(out.statistic)?),("p_value",Output::encode(out.p_value)?),("df",Output::encode(out.df)?)])))});
object_binding!(robust_bootstrap_ci,"num_native_robust_bootstrap_ci","crate::kernels::robust::bootstrap_ci","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(data:Input,statistic:StatisticFn,resamples:usize,alpha:f64,seed:u64,callback_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,data:Input,statistic:StatisticFn,resamples:usize,alpha:f64,seed:u64,callback_budget:usize|{api::same_dimension(&api::dimension(&data.0)?,&statistic.input)?;let x=<Vec<f64> as Arg>::decode(&data.0,0)?;if x.is_empty()||!(2..=100000).contains(&resamples)||x.len().checked_mul(resamples).is_none_or(|n|n>10_000_000)||!(0.0..1.0).contains(&alpha){return Err(error("bootstrap requires nonempty data,2..100000resamples,0<alpha<1 and bounded work"));}let gate=Gate::new(ctx,callback_budget)?;let(a,b)=robust::bootstrap_ci(&x,|x|gate.statistic(&statistic,x),resamples,alpha,seed);gate.done()?;quantities(vec![a,b],&(statistic.output))});
object_binding!(sum_numeric_taylor_coefficients,"num_native_sum_numeric_taylor_coefficients","crate::kernels::indefinite_sum::numeric_taylor_coefficients","Actual Taylor coefficients in the explicitly normalized coordinate u=x/coordinate_reference. Every coefficient retains output units; divide coefficient k by coordinate_reference^k for derivatives in physical coordinates.",false,(function:ScalarFn,center:Input,radius:Input,coordinate_reference:Input,value_reference:Input,terms:usize,callback_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,function:ScalarFn,center:Input,radius:Input,coordinate_reference:Input,value_reference:Input,terms:usize,callback_budget:usize|{let(cr,d)=api::reference(&coordinate_reference.0)?;let(vr,e)=api::reference(&value_reference.0)?;api::same_dimension(&d,&function.input)?;api::same_dimension(&e,&function.output)?;let(c,cd)=scalar(&center)?;let(r,rd)=scalar(&radius)?;api::same_dimension(&d,&cd)?;api::same_dimension(&d,&rd)?;if terms==0||terms>32||r<=0.0{return Err(error("Taylor terms must be1..32 and radius positive"));}let gate=Gate::new(ctx,callback_budget)?;let values=sums::numeric_taylor_coefficients(|u|gate.scalar(&function,u*cr)/vr,c/cr,r/cr,terms);gate.done()?;let coefficients=quantities(values.into_iter().map(|x|x*vr).collect(),&(e))?.0;Ok(Dynamic(api::record(vec![("coefficients",coefficients),("coordinate_reference",coordinate_reference.0),("center",center.0)])))});
object_binding!(sum_series_antidiff,"num_native_sum_series_antidiff","crate::kernels::indefinite_sum::series_antidiff","Native antidifference of coefficients in normalized coordinate u=x/coordinate_reference. Coefficients and result share physical output dimension; this is a discrete sum with unit step in u.",false,(coefficients:Input,center:Input,x:Input,coordinate_reference:Input)->Dynamic,|_ctx:&CallCtx<'_>,coefficients:Input,center:Input,x:Input,coordinate_reference:Input|{let(cr,d)=api::reference(&coordinate_reference.0)?;let(c,cd)=scalar(&center)?;let(x,xd)=scalar(&x)?;api::same_dimension(&d,&cd)?;api::same_dimension(&d,&xd)?;let e=api::dimension(&coefficients.0)?;let coefficients=<Vec<f64> as Arg>::decode(&coefficients.0,0)?;limit(coefficients.len(),32)?;extended_quantity(sums::series_antidiff(&coefficients,c/cr,x/cr),e)});
object_binding!(
    sum_indefinite_sum,
    "num_native_sum_indefinite_sum",
    "crate::kernels::indefinite_sum::indefinite_sum_with_budget",
    "Evaluates the discrete sum from config.h to x with steps of config.step. Fractional offsets use analytic continuation in normalized coordinates. The sum has the callback's output dimension. An integral interval that exceeds the callback budget is rejected before the callback runs.",
    false,
    (function: ScalarFn, x: Input, config: Input, coordinate_reference: Input, value_reference: Input, callback_budget: usize) -> Dynamic,
    |ctx, function, x, config, coordinate_reference, value_reference, callback_budget| {
        let input = normalize_recurrence(&function, &x, &config, &coordinate_reference, &value_reference)?;
        let gate = Gate::new(ctx, callback_budget)?;
        let out = sums::indefinite_sum_with_budget(
            |u| gate.scalar(&function, u * input.coordinate_reference) / input.value_reference,
            input.point,
            &input.config,
            callback_budget,
        );
        gate.done()?;
        extended_quantity(out.map_err(error)? * input.value_reference, function.output)
    }
);
object_binding!(
    sum_indefinite_product,
    "num_native_sum_indefinite_product",
    "crate::kernels::indefinite_sum::indefinite_product_with_budget",
    "Evaluates the dimensionless product of positive f(x)/value_reference factors from config.h with steps of config.step. Fractional offsets use analytic continuation of the sum of logarithms. An integral interval that exceeds the callback budget is rejected before the callback runs.",
    false,
    (function: ScalarFn, x: Input, config: Input, coordinate_reference: Input, value_reference: Input, callback_budget: usize) -> Dynamic,
    |ctx, function, x, config, coordinate_reference, value_reference, callback_budget| {
        let input = normalize_recurrence(&function, &x, &config, &coordinate_reference, &value_reference)?;
        let gate = Gate::new(ctx, callback_budget)?;
        let out = sums::indefinite_product_with_budget(
            |u| gate.scalar(&function, u * input.coordinate_reference) / input.value_reference,
            input.point,
            &input.config,
            callback_budget,
        );
        gate.done()?;
        extended_quantity(out.map_err(error)?, DimExpr::Dimless)
    }
);
object_binding!(geometry_metric_tensor_at_point,"num_native_geometry_metric_tensor_at_point","crate::kernels::differential_geometry::metric_tensor_at_point","Native finite-difference differential geometry over homogeneous physical coordinates and a typed symmetric metric callback. Explicit positive coordinate/metric references normalize derivative step scales; metric/connection/curvature result dimensions are restored. Dimension<=6 with first-failure and callback budget checks.",false,(metric:MetricFn,point:Input,coordinate_reference:Input,metric_reference:Input,callback_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,metric:MetricFn,point:Input,coordinate_reference:Input,metric_reference:Input,callback_budget:usize|{let(x,cr,mr)=metric.normalized_point(&point,&coordinate_reference,&metric_reference)?;let gate=Gate::new(ctx,callback_budget)?;let out=dg::metric_tensor_at_point(|x|gate.metric(&metric,x,cr,mr),&x);gate.done()?;let out=out.map_err(error)?;let d=metric.output;scaled_tensor(Output::encode(out)?,mr,&(d))});
object_binding!(geometry_christoffel_symbols,"num_native_geometry_christoffel_symbols","crate::kernels::differential_geometry::christoffel_symbols","Native finite-difference differential geometry over homogeneous physical coordinates and a typed symmetric metric callback. Explicit positive coordinate/metric references normalize derivative step scales; metric/connection/curvature result dimensions are restored. Dimension<=6 with first-failure and callback budget checks.",false,(metric:MetricFn,point:Input,coordinate_reference:Input,metric_reference:Input,callback_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,metric:MetricFn,point:Input,coordinate_reference:Input,metric_reference:Input,callback_budget:usize|{let(x,cr,mr)=metric.normalized_point(&point,&coordinate_reference,&metric_reference)?;let gate=Gate::new(ctx,callback_budget)?;let out=dg::christoffel_symbols(|x|gate.metric(&metric,x,cr,mr),&x);gate.done()?;let out=out.map_err(error)?;let d=canonical(&(metric.input.pow(Exp::Rat(Rational::int(-1)))))?;scaled_tensor(tensor3(out)?,1.0/cr,&(d))});
object_binding!(geometry_riemann_tensor,"num_native_geometry_riemann_tensor","crate::kernels::differential_geometry::riemann_tensor","Native finite-difference differential geometry over homogeneous physical coordinates and a typed symmetric metric callback. Explicit positive coordinate/metric references normalize derivative step scales; metric/connection/curvature result dimensions are restored. Dimension<=6 with first-failure and callback budget checks.",false,(metric:MetricFn,point:Input,coordinate_reference:Input,metric_reference:Input,callback_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,metric:MetricFn,point:Input,coordinate_reference:Input,metric_reference:Input,callback_budget:usize|{let(x,cr,mr)=metric.normalized_point(&point,&coordinate_reference,&metric_reference)?;let gate=Gate::new(ctx,callback_budget)?;let out=dg::riemann_tensor(|x|gate.metric(&metric,x,cr,mr),&x);gate.done()?;let out=out.map_err(error)?;let d=canonical(&(metric.input.pow(Exp::Rat(Rational::int(-2)))))?;scaled_tensor(tensor4(out)?,1.0/(cr*cr),&(d))});
object_binding!(geometry_ricci_tensor,"num_native_geometry_ricci_tensor","crate::kernels::differential_geometry::ricci_tensor","Native finite-difference differential geometry over homogeneous physical coordinates and a typed symmetric metric callback. Explicit positive coordinate/metric references normalize derivative step scales; metric/connection/curvature result dimensions are restored. Dimension<=6 with first-failure and callback budget checks.",false,(metric:MetricFn,point:Input,coordinate_reference:Input,metric_reference:Input,callback_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,metric:MetricFn,point:Input,coordinate_reference:Input,metric_reference:Input,callback_budget:usize|{let(x,cr,mr)=metric.normalized_point(&point,&coordinate_reference,&metric_reference)?;let gate=Gate::new(ctx,callback_budget)?;let out=dg::ricci_tensor(|x|gate.metric(&metric,x,cr,mr),&x);gate.done()?;let out=out.map_err(error)?;let d=canonical(&(metric.input.pow(Exp::Rat(Rational::int(-2)))))?;scaled_tensor(Output::encode(out)?,1.0/(cr*cr),&(d))});
object_binding!(geometry_ricci_scalar,"num_native_geometry_ricci_scalar","crate::kernels::differential_geometry::ricci_scalar","Native finite-difference differential geometry over homogeneous physical coordinates and a typed symmetric metric callback. Explicit positive coordinate/metric references normalize derivative step scales; metric/connection/curvature result dimensions are restored. Dimension<=6 with first-failure and callback budget checks.",false,(metric:MetricFn,point:Input,coordinate_reference:Input,metric_reference:Input,callback_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,metric:MetricFn,point:Input,coordinate_reference:Input,metric_reference:Input,callback_budget:usize|{let(x,cr,mr)=metric.normalized_point(&point,&coordinate_reference,&metric_reference)?;let gate=Gate::new(ctx,callback_budget)?;let out=dg::ricci_scalar(|x|gate.metric(&metric,x,cr,mr),&x);gate.done()?;let out=out.map_err(error)?;let d=canonical(&(DimExpr::Mul(vec![metric.input.clone().pow(Exp::Rat(Rational::int(-2))),DimExpr::Inverse(Box::new(metric.output))])))?;extended_quantity(out*(1.0/(cr*cr*mr)),d)});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &stats_ks_test_one_sample::BINDING,
            &robust_bootstrap_ci::BINDING,
            &sum_numeric_taylor_coefficients::BINDING,
            &sum_series_antidiff::BINDING,
            &sum_indefinite_sum::BINDING,
            &sum_indefinite_product::BINDING,
            &geometry_metric_tensor_at_point::BINDING,
            &geometry_christoffel_symbols::BINDING,
            &geometry_riemann_tensor::BINDING,
            &geometry_ricci_tensor::BINDING,
            &geometry_ricci_scalar::BINDING,
            &geometry_ricci_curvature::BINDING,
        ],
        "numerical_vector",
    )
}
