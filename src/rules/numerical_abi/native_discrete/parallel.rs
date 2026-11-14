#![allow(clippy::wildcard_imports)] // Native Rayon work with owner-thread checked host-callback dispatch.
use super::*;
use crate::kernel::value::CallbackRef;
use crate::kernels::par;
use std::sync::mpsc::Sender;
use std::sync::mpsc::{
    self,
};
struct IndexedFn {
    h: CallbackRef,
    output: DimExpr,
}
struct MappedFn {
    h: CallbackRef,
    input: DimExpr,
    output: DimExpr,
}
struct MutableFn {
    h: CallbackRef,
    dimension: DimExpr,
}
struct ChunkFn {
    h: CallbackRef,
    dimension: DimExpr,
}
fn d(t: &TypeExpr) -> Result<DimExpr, OpError> {
    match t {
        | TypeExpr::Float => Ok(DimExpr::Dimless),
        | TypeExpr::Dim(t, d) if **t == TypeExpr::Float => canonical(d),
        | _ => {
            Err(error(
                "parallel callback must use concrete Float dimensions",
            ))
        },
    }
}
fn td(t: &TypeExpr) -> Result<DimExpr, OpError> {
    match t {
        | TypeExpr::Tensor {
            rank: NatExpr::Lit(1),
            elem,
            ..
        } => d(elem),
        | _ => Err(error("parallel chunk callback must use rank1 Float tensor")),
    }
}
fn sig<'a>(
    v: &Value,
    c: &'a CallCtx<'_>,
) -> Result<(CallbackRef, &'a crate::kernel::signature::Signature), OpError> {
    let Value::Function(h) = v else {
        return Err(error("parallel operation requires typed callback"));
    };
    Ok((*h, c.callback_signature(*h)?))
}
fn ft(
    args: Vec<TypeExpr>,
    ret: TypeExpr,
) -> TypeExpr {
    TypeExpr::Function(args, Box::new(ret))
}
fn v(n: &str) -> TypeExpr {
    TypeExpr::Float.with_dim(DimExpr::Var(n.into()))
}
impl NativeArg for IndexedFn {
    fn ty() -> TypeExpr {
        ft(vec![TypeExpr::Int], v("Y"))
    }

    fn decode(
        v: &Value,
        _: usize,
        c: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        let (h, s) = sig(v, c)?;
        if s.params != vec![TypeExpr::Int] {
            return Err(error("indexed callback expects Int argument"));
        }
        Ok(Self {
            h,
            output: d(&s.ret)?,
        })
    }
}
impl NativeArg for MappedFn {
    fn ty() -> TypeExpr {
        ft(vec![v("X")], v("Y"))
    }

    fn decode(
        v: &Value,
        _: usize,
        c: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        let (h, s) = sig(v, c)?;
        if s.params.len() != 1 {
            return Err(error("map callback requires one argument"));
        }
        Ok(Self {
            h,
            input: d(&s.params[0])?,
            output: d(&s.ret)?,
        })
    }
}
impl NativeArg for MutableFn {
    fn ty() -> TypeExpr {
        ft(vec![TypeExpr::Int, v("X")], v("X"))
    }

    fn decode(
        v: &Value,
        _: usize,
        c: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        let (h, s) = sig(v, c)?;
        if s.params.len() != 2 || s.params[0] != TypeExpr::Int {
            return Err(error("mutable map callback requires index and scalar"));
        }
        let dimension = d(&s.params[1])?;
        api::same_dimension(&dimension, &d(&s.ret)?)?;
        Ok(Self { h, dimension })
    }
}
impl NativeArg for ChunkFn {
    fn ty() -> TypeExpr {
        ft(
            vec![TypeExpr::Int, api::tensor_type(1, v("X"))],
            api::tensor_type(1, v("X")),
        )
    }

    fn decode(
        v: &Value,
        _: usize,
        c: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        let (h, s) = sig(v, c)?;
        if s.params.len() != 2 || s.params[0] != TypeExpr::Int {
            return Err(error("chunk callback requires index and vector"));
        }
        let dimension = td(&s.params[1])?;
        api::same_dimension(&dimension, &td(&s.ret)?)?;
        Ok(Self { h, dimension })
    }
}
struct Request {
    index: usize,
    values: Vec<f64>,
    reply: Sender<Vec<f64>>,
    size: usize,
}
enum Event<R> {
    Request(Request),
    Done(Result<R, String>),
}
fn request<R>(
    tx: &Sender<Event<R>>,
    index: usize,
    values: Vec<f64>,
    size: usize,
) -> Vec<f64> {
    let (a, b) = mpsc::channel();
    if tx
        .send(Event::Request(Request {
            index,
            values,
            reply: a,
            size,
        }))
        .is_err()
    {
        return vec![f64::NAN; size];
    }
    b.recv().unwrap_or_else(|_| vec![f64::NAN; size])
}
fn run<R: Send>(
    budget: usize,
    work: impl FnOnce(Sender<Event<R>>) -> R + Send,
    invoke: impl Fn(usize, &[f64]) -> Result<Vec<f64>, OpError>,
) -> Result<R, OpError> {
    if budget == 0 || budget > 1_000_000 {
        return Err(error("parallel callback budget must be1..1000000"));
    }
    std::thread::scope(|scope| {
        let (tx, rx) = mpsc::channel();
        let thread = scope.spawn(move || {
            let result =
                std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| work(tx.clone())))
                    .map_err(|_| "native parallel worker panicked".to_string());
            let _ = tx.send(Event::Done(result));
        });
        let mut first = None;
        let mut remaining = budget;
        let result = loop {
            match rx.recv() {
                | Ok(Event::Request(r)) => {
                    let out = if first.is_some() {
                        Err(error("parallel callback already failed"))
                    } else if remaining == 0 {
                        Err(error("parallel callback budget exhausted"))
                    } else {
                        remaining -= 1;
                        invoke(r.index, &r.values)
                    };
                    let values = match out {
                        | Ok(values)
                            if values.len() == r.size && values.iter().all(|x| x.is_finite()) =>
                        {
                            values
                        },
                        | Ok(_) => {
                            if first.is_none() {
                                first = Some(error(
                                    "parallel callback result has invalid size or nonfinite values",
                                ));
                            }
                            vec![f64::NAN; r.size]
                        },
                        | Err(e) => {
                            if first.is_none() {
                                first = Some(e);
                            }
                            vec![f64::NAN; r.size]
                        },
                    };
                    let _ = r.reply.send(values);
                },
                | Ok(Event::Done(result)) => break result.map_err(error),
                | Err(_) => break Err(error("parallel worker terminated without completion")),
            }
        };
        thread
            .join()
            .map_err(|_| error("parallel worker panicked during completion"))?;
        if let Some(e) = first {
            return Err(e);
        }
        result
    })
}
fn sizes(
    n: usize,
    grain: usize,
) -> Result<(), OpError> {
    limit(n, 262144)?;
    if grain == 0 || grain > 262144 {
        return Err(error("parallel grain must be1..262144"));
    }
    Ok(())
}
fn scalar_result(
    ctx: &CallCtx<'_>,
    h: CallbackRef,
    args: &[Value],
    dimension: &DimExpr,
) -> Result<Vec<f64>, OpError> {
    let out = ctx.invoke_callback(h, args)?;
    api::same_dimension(&api::dimension(&out)?, dimension)?;
    Ok(vec![numeric(&out, 0)?])
}
object_binding!(par_map_indexed,"num_native_par_map_indexed","crate::kernels::par::map_indexed","Actual native parallel map with host callbacks serviced serially on their owning thread. Index order and native Rayon scheduling are retained; every worker is joined and callback failures/budgets are propagated.",false,(n:usize,grain:usize,function:IndexedFn,callback_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,n:usize,grain:usize,function:IndexedFn,callback_budget:usize|{sizes(n,grain)?;let out=run(callback_budget,move|tx|par::map_indexed(n,grain,|i|request(&tx,i,vec![],1)[0]),|i,_|scalar_result(ctx,function.h,&[Value::int(i)],&function.output))?;quantities(out,&(function.output))});
object_binding!(par_map_slice,"num_native_par_map_slice","crate::kernels::par::map_slice","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(items:Input,grain:usize,function:MappedFn,callback_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,items:Input,grain:usize,function:MappedFn,callback_budget:usize|{api::same_dimension(&api::dimension(&items.0)?,&function.input)?;let data=<Vec<f64> as Arg>::decode(&items.0,0)?;sizes(data.len(),grain)?;let out=run(callback_budget,move|tx|par::map_slice(&data,grain,|x|request(&tx,0,vec![*x],1)[0]),|_,x|scalar_result(ctx,function.h,&[api::quantity(x[0],function.input.clone())?],&function.output))?;quantities(out,&(function.output))});
object_binding!(par_sum_fixed,"num_native_par_sum_fixed","crate::kernels::par::sum_fixed","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(n:usize,block:usize,grain:usize,function:IndexedFn,callback_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,n:usize,block:usize,grain:usize,function:IndexedFn,callback_budget:usize|{sizes(n,grain)?;if block==0||block>262144{return Err(error("parallel block must be1..262144"));}let out=run(callback_budget,move|tx|par::sum_fixed(n,block,grain,|i|request(&tx,i,vec![],1)[0]),|i,_|scalar_result(ctx,function.h,&[Value::int(i)],&function.output))?;extended_quantity(out,function.output)});
object_binding!(par_for_each_mut,"num_native_par_for_each_mut","crate::kernels::par::for_each_mut","Native mutable parallel map specializes mutation to returning a replacement value. The input tensor is owned/copied; outputs retain identical physical units and original index semantics.",false,(items:Input,grain:usize,function:MutableFn,callback_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,items:Input,grain:usize,function:MutableFn,callback_budget:usize|{api::same_dimension(&api::dimension(&items.0)?,&function.dimension)?;let mut data=<Vec<f64> as Arg>::decode(&items.0,0)?;sizes(data.len(),grain)?;let out=run(callback_budget,move|tx|{par::for_each_mut(&mut data,grain,|i,x|*x=request(&tx,i,vec![*x],1)[0]);data},|i,x|scalar_result(ctx,function.h,&[Value::int(i),api::quantity(x[0],function.dimension.clone())?],&function.dimension))?;quantities(out,&(function.dimension))});
object_binding!(par_for_chunks_mut,"num_native_par_for_chunks_mut","crate::kernels::par::for_chunks_mut","Native disjoint mutable chunk scheduling specializes mutation to returning a replacement vector of identical chunk length and units, including the short final chunk. Checked owner-thread callbacks preserve failure and budget semantics.",false,(items:Input,chunk:usize,grain:usize,function:ChunkFn,callback_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,items:Input,chunk:usize,grain:usize,function:ChunkFn,callback_budget:usize|{api::same_dimension(&api::dimension(&items.0)?,&function.dimension)?;let mut data=<Vec<f64> as Arg>::decode(&items.0,0)?;sizes(data.len(),grain)?;if chunk==0||chunk>262144{return Err(error("parallel chunk must be1..262144"));}let out=run(callback_budget,move|tx|{par::for_chunks_mut(&mut data,chunk,grain,|i,x|x.copy_from_slice(&request(&tx,i,x.to_vec(),x.len())));data},|i,x|{let a=quantities(x.to_vec(),&function.dimension)?.0;let out=ctx.invoke_callback(function.h,&[Value::int(i),a])?;api::same_dimension(&api::dimension(&out)?,&function.dimension)?;<Vec<f64> as Arg>::decode(&out,0)})?;quantities(out,&(function.dimension))});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &par_map_indexed::BINDING,
            &par_map_slice::BINDING,
            &par_sum_fixed::BINDING,
            &par_for_each_mut::BINDING,
            &par_for_chunks_mut::BINDING,
        ],
        "numerical_vector",
    )
}
