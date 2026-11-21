//! Shared checked object and method adapters. Numerical work stays in native types.
use super::codecs::Arg;
use super::codecs::Output;
use super::numeric;
use super::output_scalar;
use super::runtime_dimensions;
use crate::graph::Arity;
use crate::graph::OpDescriptor;
use crate::graph::RuleError;
use crate::graph::op::EvalDyn;
use crate::graph::rule::Installer;
use crate::kernel::call::CallCtx;
use crate::kernel::call::OpError;
use crate::kernel::native_objects::NativeObjects;
use crate::kernel::signature::ArityInfo;
use crate::kernel::signature::BackendSupport;
use crate::kernel::signature::Determinism;
use crate::kernel::signature::DimExpr;
use crate::kernel::signature::Effect;
use crate::kernel::signature::EffectRow;
use crate::kernel::signature::NatExpr;
use crate::kernel::signature::Purity;
use crate::kernel::signature::RepExpr;
use crate::kernel::signature::Signature;
use crate::kernel::signature::Stability;
use crate::kernel::signature::ThreadSafety;
use crate::kernel::signature::TypeExpr;
use crate::kernel::value::DynValue;
use crate::kernel::value::RecordValue;
use crate::kernel::value::Value;
use std::any::Any;
use std::marker::PhantomData;
use std::sync::Arc;

pub(in crate::rules) const VALUE_LIMIT: usize = 262_144;
pub(in crate::rules) fn error(message: impl std::fmt::Display) -> OpError {
    OpError::Domain(message.to_string())
}
pub(in crate::rules) fn bank<'a>(ctx: &'a CallCtx<'_>) -> Result<&'a NativeObjects, OpError> {
    ctx.native_objects
        .as_ref()
        .ok_or_else(|| error("native objects require a session context"))
}
pub(in crate::rules) fn get<T: Any, R>(
    ctx: &CallCtx<'_>,
    handle: &Value,
    kind: &str,
    f: impl FnOnce(&T) -> Result<R, OpError>,
) -> Result<R, OpError> {
    bank(ctx)?.with(handle, kind, f)
}
pub(in crate::rules) fn get_mut<T: Any, R>(
    ctx: &CallCtx<'_>,
    handle: &Value,
    kind: &str,
    f: impl FnOnce(&mut T) -> Result<R, OpError>,
) -> Result<R, OpError> {
    bank(ctx)?.with_mut(handle, kind, f)
}
pub(in crate::rules) fn store<T: Any>(
    ctx: &CallCtx<'_>,
    kind: &str,
    dimensions: Vec<DimExpr>,
    value: T,
    bytes: usize,
) -> Result<Value, OpError> {
    bank(ctx)?.insert(kind, dimensions, value, bytes)
}
pub(in crate::rules) fn dimensions(
    ctx: &CallCtx<'_>,
    handle: &Value,
) -> Result<Vec<DimExpr>, OpError> {
    bank(ctx)?.dimensions(handle)
}
pub(in crate::rules) fn dimension(value: &Value) -> Result<DimExpr, OpError> {
    let dim = match value {
        | Value::Quantity(q) => {
            if matches!(&*q.value, Value::Quantity(_) | Value::Dyn(_)) {
                same_dimension(&dimension(&q.value)?, &DimExpr::Dimless)?;
            }
            q.unit.dimension.clone()
        },
        | Value::Dyn(d) => {
            if matches!(
                &*d.value,
                Value::Quantity(_) | Value::Dyn(_) | Value::Tensor(_)
            ) {
                same_dimension(&dimension(&d.value)?, &d.dim)?;
            }
            d.dim.clone()
        },
        | Value::Tensor(t) => {
            let declared = type_dimension(&t.elem)?;
            for v in &t.data {
                if matches!(v, Value::Quantity(_) | Value::Dyn(_) | Value::Tensor(_)) {
                    same_dimension(&dimension(v)?, &declared)?;
                }
            }
            return Ok(declared);
        },
        | _ => DimExpr::Dimless,
    };
    runtime_dimensions::expression(&runtime_dimensions::canonical(&dim)?)
}
fn type_dimension(ty: &TypeExpr) -> Result<DimExpr, OpError> {
    match ty {
        | TypeExpr::Dim(_, d) => runtime_dimensions::expression(&runtime_dimensions::canonical(d)?),
        | TypeExpr::Tensor { elem, .. } => type_dimension(elem),
        | TypeExpr::Point(_) => Err(error("native numerical values cannot be affine points")),
        | _ => Ok(DimExpr::Dimless),
    }
}
pub(in crate::rules) fn same_dimension(
    a: &DimExpr,
    b: &DimExpr,
) -> Result<(), OpError> {
    if runtime_dimensions::canonical(a)? == runtime_dimensions::canonical(b)? {
        Ok(())
    } else {
        Err(error("native argument physical dimensions disagree"))
    }
}
pub(in crate::rules) fn quantity(
    value: f64,
    dim: DimExpr,
) -> Result<Value, OpError> {
    let value = output_scalar(value)?;
    if dim == DimExpr::Dimless {
        Ok(value)
    } else {
        Ok(Value::Dyn(DynValue {
            dim,
            value: Box::new(value),
        }))
    }
}
pub(in crate::rules) fn reference(value: &Value) -> Result<(f64, DimExpr), OpError> {
    let magnitude = numeric(value, 0)?;
    if magnitude <= 0.0 {
        return Err(error("reference scale must be finite and positive"));
    }
    Ok((magnitude, dimension(value)?))
}
pub(in crate::rules) fn dynamic_type() -> TypeExpr {
    TypeExpr::Dyn(RepExpr {
        repr: "term".into(),
        dim: None,
    })
}
pub(in crate::rules) fn tensor_type(
    rank: u64,
    elem: TypeExpr,
) -> TypeExpr {
    TypeExpr::Tensor {
        rank: NatExpr::Lit(rank),
        shape: None,
        elem: Box::new(elem),
    }
}
pub(in crate::rules) fn record(fields: Vec<(&str, Value)>) -> Value {
    Value::Record(RecordValue {
        fields: fields.into_iter().map(|(k, v)| (k.into(), v)).collect(),
    })
}
pub(in crate::rules) fn field<'a>(
    v: &'a Value,
    name: &str,
) -> Result<&'a Value, OpError> {
    let Value::Record(r) = v else {
        return Err(error("expected a record"));
    };
    let mut found = r.fields.iter().filter(|(k, _)| k == name);
    let result = &found
        .next()
        .ok_or_else(|| error(format!("missing field {name}")))?
        .1;
    if found.next().is_some() {
        return Err(error("duplicate record field"));
    }
    Ok(result)
}
/// Bound shape, nested data, text and retained handle encodings before native decoding.
pub(in crate::rules) fn bounded(value: &Value) -> Result<(), OpError> {
    fn visit_children(
        values: &[Value],
        depth: usize,
        nodes: &mut usize,
    ) -> Result<(), OpError> {
        if values.iter().all(|value| {
            matches!(
                value,
                Value::Bool(_) | Value::Float(_) | Value::Int(_) | Value::Complex(_)
            )
        }) {
            if !values.is_empty() && (depth > 32 || values.len() > VALUE_LIMIT - *nodes) {
                return Err(error("native value exceeds element/depth budget"));
            }
            *nodes += values.len();
        } else {
            for child in values.iter().rev() {
                visit(child, depth, nodes)?;
            }
        }
        Ok(())
    }

    // Descend one child at a time so a wide tensor needs no traversal buffer.
    fn visit(
        value: &Value,
        depth: usize,
        nodes: &mut usize,
    ) -> Result<(), OpError> {
        *nodes += 1;
        if *nodes > VALUE_LIMIT || depth > 32 {
            return Err(error("native value exceeds element/depth budget"));
        }
        match value {
            | Value::Tensor(t) => {
                if t.shape.len() > 8
                    || t.shape.iter().any(|d| *d > VALUE_LIMIT)
                    || t.shape.iter().try_fold(1usize, |n, d| n.checked_mul(*d))
                        != Some(t.data.len())
                    || t.data.len() > VALUE_LIMIT
                {
                    return Err(error("native tensor has invalid shape or exceeds budget"));
                }
                let mut prefix = 1usize;
                let mut slots = 0usize;
                for extent in &t.shape {
                    prefix = prefix
                        .checked_mul(*extent)
                        .ok_or_else(|| error("native tensor nesting overflows"))?;
                    slots = slots
                        .checked_add(prefix)
                        .filter(|n| *n <= VALUE_LIMIT)
                        .ok_or_else(|| error("native nested tensor allocation budget exceeded"))?;
                }
                return visit_children(&t.data, depth + 1, nodes);
            },
            | Value::Record(r) => {
                if r.fields.len() > 256 {
                    return Err(error("native record field budget exceeded"));
                }
                let mut names = std::collections::BTreeSet::new();
                for (k, _) in &r.fields {
                    if k.len() > 256 || !names.insert(k) {
                        return Err(error(
                            "native record has duplicate or excessive field names",
                        ));
                    }
                }
                drop(names);
                for (_, child) in r.fields.iter().rev() {
                    visit(child, depth + 1, nodes)?;
                }
            },
            | Value::Variant(v) => return visit_children(&v.payload, depth + 1, nodes),
            | Value::Quantity(q) => return visit(&q.value, depth + 1, nodes),
            | Value::Dyn(d) => return visit(&d.value, depth + 1, nodes),
            | Value::String(s) if s.len() > 1_048_576 => {
                return Err(error("native text budget exceeded"));
            },
            | Value::Opaque(o) if o.bytes.len() > 16_777_216 => {
                return Err(error("native opaque payload budget exceeded"));
            },
            | _ => {},
        }
        Ok(())
    }
    let mut nodes = 0;
    visit(value, 0, &mut nodes)
}
pub(in crate::rules) fn normalized(value: &Value) -> Result<(), OpError> {
    fn ty(t: &TypeExpr) -> Result<(), OpError> {
        match t {
            | TypeExpr::Dim(inner, d) => {
                if !runtime_dimensions::canonical(d)?.is_empty() {
                    return Err(error(
                        "normalized native coordinates require explicit nondimensionalization",
                    ));
                }
                ty(inner)
            },
            | TypeExpr::Point(_) => {
                Err(error(
                    "normalized native coordinates cannot be affine points",
                ))
            },
            | TypeExpr::Tensor { elem, .. } | TypeExpr::Sparse { elem, .. } => ty(elem),
            | _ => Ok(()),
        }
    }
    bounded(value)?;
    ty(&value.type_expr())?;
    match value {
        | Value::Quantity(q) => {
            if !runtime_dimensions::canonical(&q.unit.dimension)?.is_empty() {
                return Err(error(
                    "normalized native coordinates require explicit reference scales",
                ));
            }
            normalized(&q.value)
        },
        | Value::Dyn(d) => {
            if !runtime_dimensions::canonical(&d.dim)?.is_empty() {
                return Err(error(
                    "normalized native coordinates require explicit reference scales",
                ));
            }
            normalized(&d.value)
        },
        | Value::Tensor(t) => {
            for v in &t.data {
                normalized(v)?;
            }
            Ok(())
        },
        | Value::Record(r) => {
            for (_, v) in &r.fields {
                normalized(v)?;
            }
            Ok(())
        },
        | Value::Variant(v) => {
            for v in &v.payload {
                normalized(v)?;
            }
            Ok(())
        },
        | _ => Ok(()),
    }
}
/// Native object contract. Mutable values reserve their maximum retained allocation
/// at creation. Immutable values can report their allocated capacity once.
pub(in crate::rules) trait NativeType: Any {
    const KIND: &'static str;
    const MAX_BYTES: usize;
    fn validate(&self) -> Result<(), OpError>;

    fn retained_bytes(&self) -> usize {
        Self::MAX_BYTES
    }
}
/// A typed, authenticated capability, not a decoded native layout.
pub(in crate::rules) struct Handle<T: NativeType> {
    pub value: Value,
    _marker: PhantomData<T>,
}
impl<T: NativeType> Clone for Handle<T> {
    fn clone(&self) -> Self {
        Self {
            value: self.value.clone(),
            _marker: PhantomData,
        }
    }
}
impl<T: NativeType> Handle<T> {
    pub const fn heap_bytes(&self) -> usize {
        match &self.value {
            | Value::Opaque(opaque) => opaque.kind.capacity() + opaque.bytes.capacity(),
            | _ => 0,
        }
    }

    pub fn with<R>(
        &self,
        ctx: &CallCtx<'_>,
        f: impl FnOnce(&T) -> Result<R, OpError>,
    ) -> Result<R, OpError> {
        get(ctx, &self.value, T::KIND, f)
    }

    pub fn with_mut<R>(
        &self,
        ctx: &CallCtx<'_>,
        f: impl FnOnce(&mut T) -> Result<R, OpError>,
    ) -> Result<R, OpError> {
        get_mut(ctx, &self.value, T::KIND, f)
    }

    pub fn dimensions(
        &self,
        ctx: &CallCtx<'_>,
    ) -> Result<Vec<DimExpr>, OpError> {
        dimensions(ctx, &self.value)
    }
}
pub(in crate::rules) struct NewObject<T: NativeType> {
    pub value: T,
    pub dimensions: Vec<DimExpr>,
    pub parents: Vec<Value>,
}
impl<T: NativeType> NewObject<T> {
    pub const fn new(
        value: T,
        dimensions: Vec<DimExpr>,
    ) -> Self {
        Self {
            value,
            dimensions,
            parents: vec![],
        }
    }

    pub const fn dependent(
        value: T,
        dimensions: Vec<DimExpr>,
        parents: Vec<Value>,
    ) -> Self {
        Self {
            value,
            dimensions,
            parents,
        }
    }

    pub const fn normalized(value: T) -> Self {
        Self::new(value, vec![])
    }
}
pub(in crate::rules) trait NativeArg: Sized {
    fn ty() -> TypeExpr;
    fn decode(
        value: &Value,
        index: usize,
        ctx: &CallCtx<'_>,
    ) -> Result<Self, OpError>;
}
impl<T: Arg> NativeArg for T {
    fn ty() -> TypeExpr {
        T::ty()
    }

    fn decode(
        v: &Value,
        index: usize,
        _: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        normalized(v)?;
        T::decode(v, index)
    }
}
impl<T: NativeType> NativeArg for Handle<T> {
    fn ty() -> TypeExpr {
        TypeExpr::Opaque(T::KIND.into())
    }

    fn decode(
        v: &Value,
        _: usize,
        ctx: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        get::<T, _>(ctx, v, T::KIND, NativeType::validate)?;
        Ok(Self {
            value: v.clone(),
            _marker: PhantomData,
        })
    }
}
pub(in crate::rules) trait NativeOutput {
    fn ty() -> TypeExpr;
    fn encode(
        self,
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError>;
}
impl<T: Output> NativeOutput for T {
    fn ty() -> TypeExpr {
        T::ty()
    }

    fn encode(
        self,
        _: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        self.encode()
    }
}
impl<T: NativeType> NativeOutput for NewObject<T> {
    fn ty() -> TypeExpr {
        TypeExpr::Opaque(T::KIND.into())
    }

    fn encode(
        self,
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        self.value.validate()?;
        let retained_bytes = self.value.retained_bytes();
        if retained_bytes > T::MAX_BYTES {
            return Err(error("native object exceeds its retained allocation bound"));
        }
        bank(ctx)?.insert_dependent(
            T::KIND,
            self.dimensions,
            self.value,
            retained_bytes,
            &self.parents,
        )
    }
}
impl<T: NativeType> NativeOutput for Handle<T> {
    fn ty() -> TypeExpr {
        TypeExpr::Opaque(T::KIND.into())
    }

    fn encode(
        self,
        _: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        Ok(self.value)
    }
}
/// Escape hatch for a checked heterogeneous/dimensioned result. The signature is
/// deliberately Value; it never fabricates a compile-time-known numeric domain.
pub(in crate::rules) struct Dynamic(pub Value);
impl NativeOutput for Dynamic {
    fn ty() -> TypeExpr {
        dynamic_type()
    }

    fn encode(
        self,
        _: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        bounded(&self.0)?;
        Ok(self.0)
    }
}
pub(in crate::rules) struct Binding {
    pub name: &'static str,
    pub source: &'static str,
    pub doc: &'static str,
    pub signature: fn() -> Signature,
    pub arity: u8,
    pub call: fn(&[Value], &CallCtx<'_>) -> Result<Value, OpError>,
    pub normalized: bool,
}
fn dimension_variables(
    ty: &TypeExpr,
    names: &mut std::collections::BTreeSet<String>,
) {
    fn dim(
        d: &DimExpr,
        n: &mut std::collections::BTreeSet<String>,
    ) {
        match d {
            | DimExpr::Var(v) => {
                n.insert(v.clone());
            },
            | DimExpr::Mul(parts) => {
                for p in parts {
                    dim(p, n);
                }
            },
            | DimExpr::Pow(inner, _) | DimExpr::Inverse(inner) => dim(inner, n),
            | _ => {},
        }
    }
    match ty {
        | TypeExpr::Dim(t, d) => {
            dim(d, names);
            dimension_variables(t, names);
        },
        | TypeExpr::Point(d) | TypeExpr::Unit(d) => dim(d, names),
        | TypeExpr::Tensor { elem, .. } | TypeExpr::Sparse { elem, .. } => {
            dimension_variables(elem, names);
        },
        | TypeExpr::Record(fields) => {
            for (_, t) in fields {
                dimension_variables(t, names);
            }
        },
        | TypeExpr::Variant(cases) => {
            for (_, types) in cases {
                for t in types {
                    dimension_variables(t, names);
                }
            }
        },
        | TypeExpr::Function(args, ret) => {
            for a in args {
                dimension_variables(a, names);
            }
            dimension_variables(ret, names);
        },
        | TypeExpr::Dyn(rep) => {
            if let Some(d) = &rep.dim {
                dim(d, names);
            }
        },
        | _ => {},
    }
}
impl EvalDyn for &'static Binding {
    fn signature(&self) -> Signature {
        let mut s = (self.signature)();
        let mut names = std::collections::BTreeSet::new();
        for t in s.params.iter().chain(std::iter::once(&s.ret)) {
            dimension_variables(t, &mut names);
        }
        for name in names {
            if !s.type_params.iter().any(|p| p.name == name) {
                s.type_params.push(crate::kernel::signature::TypeParam {
                    name,
                    kind: crate::kernel::signature::TypeParamKind::Dim,
                    bounds: vec![],
                });
            }
        }
        s.effects = EffectRow {
            effects: vec![Effect::Stateful, Effect::Diverges],
        };
        if self.source.contains("simulate_and_save")
            || self.source
                == "crate::sim::models_linear_elasticity::simulate_cantilever_beam_scenario"
            || self.source
                == "crate::sim::models_navier_stokes_fluid::simulate_lid_driven_cavity_scenario"
            || self.source
                == "crate::sim::models_gpe_superfluidity::simulate_bose_einstein_vortex_scenario"
            || self.source
                == "crate::sim::models_schrodinger_quantum::simulate_double_slit_scenario"
            || self.source
                == "crate::sim::models_ising_statistical::simulate_ising_phase_transition_scenario"
            || self.source
                == "crate::sim::models_geodesic_relativity::simulate_black_hole_orbits_scenario"
        {
            s.effects.effects.push(Effect::Io);
        }
        s
    }

    fn eval(
        &self,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        if args.len() != usize::from(self.arity) {
            return Err(OpError::Arity {
                expected: ArityInfo::fixed(self.arity),
                got: args.len(),
            });
        }
        for v in args {
            bounded(v)?;
            if self.normalized {
                normalized(v)?;
            }
        }
        let scope = bank(ctx)?.scoped()?;
        let result = (self.call)(args, ctx);
        if result.is_ok() {
            scope.commit();
        }
        result
    }
}
fn has_callback(ty: &TypeExpr) -> bool {
    match ty {
        | TypeExpr::Function(..) => true,
        | TypeExpr::Tensor { elem, .. }
        | TypeExpr::Sparse { elem, .. }
        | TypeExpr::Dim(elem, _) => has_callback(elem),
        | TypeExpr::Record(fields) => fields.iter().any(|(_, t)| has_callback(t)),
        | TypeExpr::Variant(cases) => {
            cases
                .iter()
                .any(|(_, fields)| fields.iter().any(has_callback))
        },
        | _ => false,
    }
}
pub(in crate::rules) fn install(
    i: &mut Installer<'_>,
    bindings: &[&'static Binding],
    namespace: &str,
) -> Result<(), RuleError> {
    for b in bindings {
        let s = b.signature();
        let higher = s.params.iter().any(has_callback);
        let mut d = OpDescriptor::new(b.name, Arity::Fixed(b.arity))
            .category(namespace)
            .doc(b.doc)
            .stability(Stability::Experimental)
            .purity(Purity::Effectful)
            .determinism(Determinism::EnvDependent)
            .thread_safety(ThreadSafety::SessionBound)
            .backend(BackendSupport {
                interpreter: true,
                jit: false,
                ffi: false,
                gpu: false,
            })
            .signature(s)
            .eval_dyn(Arc::new(*b));
        if higher {
            d = d.higher_order();
        }
        if !b.source.is_empty() {
            d = d.rust_path(b.source);
        }
        let op = i.op(d)?;
        i.dim(op, crate::graph::dim::DimRule::Container);
    }
    Ok(())
}
/// Declare actual typed calls once, deriving arity, schemas and checked codecs.
macro_rules! object_binding {
 ($module:ident,$name:literal,$source:literal,$doc:literal,$normalized:expr,($($arg:ident:$ty:ty),*$(,)?) -> $ret:ty,$body:expr)=>{
 pub(super) mod $module {
 #[allow(unused_imports)]
 use super::*;
 const ARITY:u8={let arguments:&[&str]=&[$(stringify!($arg)),*];assert!(arguments.len() <= u8::MAX as usize);arguments.len()as u8};
 fn signature()->crate::kernel::signature::Signature {crate::kernel::signature::Signature::new(vec![$(<$ty as crate::rules::numerical_abi::native_object_api::NativeArg>::ty()),*],<$ret as crate::rules::numerical_abi::native_object_api::NativeOutput>::ty())}
 fn call(args:&[crate::kernel::Value],ctx:&crate::kernel::CallCtx<'_>)->Result<crate::kernel::Value,crate::kernel::OpError>{#[allow(unused_mut)]let mut index=0usize;let _=args;$(let $arg=<$ty as crate::rules::numerical_abi::native_object_api::NativeArg>::decode(&args[index],index,ctx)?;index+=1;)*let _=index;let body:fn(&crate::kernel::CallCtx<'_>,$($ty),*)->Result<$ret,crate::kernel::OpError>=$body;let result:$ret=body(ctx,$($arg),*)?;<$ret as crate::rules::numerical_abi::native_object_api::NativeOutput>::encode(result,ctx)}
 pub(super) static BINDING:crate::rules::numerical_abi::native_object_api::Binding=crate::rules::numerical_abi::native_object_api::Binding{name:$name,source:$source,doc:$doc,signature,arity:ARITY,call,normalized:$normalized};
 }
 };
}
pub(in crate::rules) use object_binding;

/// A real rank-one tensor with a checked uniform dimension and bounded length.
/// `values` are decoded into coherent units before the numerical method runs.
pub(in crate::rules) struct MeasuredVector<const MAX: usize> {
    pub values: Vec<f64>,
    pub dimension: DimExpr,
}

impl<const MAX: usize> NativeArg for MeasuredVector<MAX> {
    fn ty() -> TypeExpr {
        tensor_type(1, TypeExpr::Float.with_dim(DimExpr::Var("Samples".into())))
    }

    fn decode(
        value: &Value,
        index: usize,
        _: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        let Value::Tensor(tensor) = value else {
            return Err(error("expected a measured rank-one tensor"));
        };
        if tensor.shape != [tensor.data.len()] || tensor.data.len() > MAX {
            return Err(error("measured samples exceed the shape or work budget"));
        }
        Ok(Self {
            dimension: dimension(value)?,
            values: <Vec<f64> as Arg>::decode(value, index)?,
        })
    }
}

/// Explicit dynamically typed argument, still bounded before use.
pub(in crate::rules) struct Input(pub Value);
impl NativeArg for Input {
    fn ty() -> TypeExpr {
        dynamic_type()
    }

    fn decode(
        v: &Value,
        _: usize,
        _: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        bounded(v)?;
        Ok(Self(v.clone()))
    }
}
/// Optional native results preserve an actual None/Some envelope.
pub(in crate::rules) struct Maybe<T>(pub Option<T>);
impl<T: NativeArg> NativeArg for Maybe<T> {
    fn ty() -> TypeExpr {
        TypeExpr::Variant(vec![
            ("None".into(), vec![]),
            ("Some".into(), vec![T::ty()]),
        ])
    }

    fn decode(
        v: &Value,
        index: usize,
        ctx: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        bounded(v)?;
        let Value::Variant(v) = v else {
            return Err(error("expected native None/Some argument"));
        };
        match (v.tag.as_str(), v.payload.as_slice()) {
            | ("None", []) => Ok(Self(None)),
            | ("Some", [value]) => Ok(Self(Some(T::decode(value, index, ctx)?))),
            | _ => Err(error("invalid native Option tag/payload")),
        }
    }
}
impl<T: NativeOutput> NativeOutput for Maybe<T> {
    fn ty() -> TypeExpr {
        TypeExpr::Variant(vec![
            ("None".into(), vec![]),
            ("Some".into(), vec![T::ty()]),
        ])
    }

    fn encode(
        self,
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        Ok(Value::Variant(crate::kernel::value::VariantValue {
            tag: if self.0.is_some() {
                "Some"
            } else {
                "None"
            }
            .into(),
            payload: self
                .0
                .map(|v| v.encode(ctx))
                .transpose()?
                .into_iter()
                .collect(),
        }))
    }
}
/// A homogeneous list of native results with exact element schema.
pub(in crate::rules) struct Many<T>(pub Vec<T>);
impl<T: NativeArg> NativeArg for Many<T> {
    fn ty() -> TypeExpr {
        tensor_type(1, T::ty())
    }

    fn decode(
        v: &Value,
        index: usize,
        ctx: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        bounded(v)?;
        let Value::Tensor(t) = v else {
            return Err(error("expected native argument vector"));
        };
        if t.shape != [t.data.len()] {
            return Err(error("native argument vector shape mismatch"));
        }
        Ok(Self(
            t.data
                .iter()
                .map(|v| T::decode(v, index, ctx))
                .collect::<Result<_, _>>()?,
        ))
    }
}
impl<T: NativeOutput> NativeOutput for Many<T> {
    fn ty() -> TypeExpr {
        tensor_type(1, T::ty())
    }

    fn encode(
        self,
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        if self.0.len() > VALUE_LIMIT {
            return Err(error("native result count budget exceeded"));
        }
        let n = self.0.len();
        Ok(Value::Tensor(crate::kernel::value::TensorValue {
            shape: vec![n],
            elem: T::ty(),
            data: self
                .0
                .into_iter()
                .map(|v| v.encode(ctx))
                .collect::<Result<_, _>>()?,
        }))
    }
}

pub(in crate::rules) fn scale_product(
    value: f64,
    numerator: &[f64],
    denominator: &[f64],
) -> Result<f64, OpError> {
    crate::kernels::numeric_common::scale_product(value, numerator, denominator).map_err(error)
}

/// Scales an input whose nonzero value must survive normalization.
pub(in crate::rules) fn scale_without_underflow(
    value: f64,
    numerator: &[f64],
    denominator: &[f64],
) -> Result<f64, OpError> {
    let scaled = scale_product(value, numerator, denominator)?;
    if value != 0.0 && scaled == 0.0 {
        return Err(error(
            "nonzero input underflows in the selected reference frame",
        ));
    }
    Ok(scaled)
}

/// Reads a scalar while preserving infinite endpoints.
pub(in crate::rules) fn extended_numeric(value: &Value) -> Result<f64, OpError> {
    match value {
        | Value::Float(v) if v.is_infinite() => Ok(*v),
        | Value::Dyn(d) => extended_numeric(&d.value),
        | Value::Quantity(q) => {
            let unit = Value::Quantity(crate::kernel::value::QuantityValue {
                value: Box::new(Value::Float(1.0)),
                unit: q.unit.clone(),
            });
            Ok(extended_numeric(&q.value)? * numeric(&unit, 0)?)
        },
        | _ => numeric(value, 0),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kernel::value::TensorValue;
    use crate::kernel::value::VariantValue;

    object_binding!(
        zero_args, "test_zero_args", "", "Test only.", true,
        () -> bool,
        |_ctx| Ok(true)
    );
    object_binding!(
        two_args, "test_two_args", "", "Test only.", true,
        (a: bool, b: bool) -> bool,
        |_ctx, a, b| Ok(a && b)
    );

    #[test]
    fn declared_arity_is_shared_by_signature_and_dispatch() {
        let ctx = CallCtx::new();
        for (binding, n) in [(&zero_args::BINDING, 0), (&two_args::BINDING, 2)] {
            assert_eq!(usize::from(binding.arity), n);
            assert_eq!(binding.signature().params.len(), n);
            assert!(matches!(
                binding.eval(&vec![Value::Bool(true); n + 1], &ctx),
                Err(OpError::Arity { got, .. }) if got == n + 1
            ));
        }
    }

    fn tensor(
        shape: Vec<usize>,
        data: Vec<Value>,
    ) -> Value {
        Value::Tensor(TensorValue {
            shape,
            elem: TypeExpr::Bool,
            data,
        })
    }

    #[test]
    fn value_limits_are_checked_in_original_child_order() {
        let wide = tensor(
            vec![VALUE_LIMIT - 1],
            vec![Value::Bool(true); VALUE_LIMIT - 1],
        );
        assert!(bounded(&wide).is_ok());
        drop(wide);
        let too_wide = tensor(vec![VALUE_LIMIT], vec![Value::Bool(true); VALUE_LIMIT]);
        assert!(
            matches!(bounded(&too_wide), Err(OpError::Domain(msg)) if msg.contains("element/depth"))
        );

        let never_allocates = tensor(vec![VALUE_LIMIT, 2, 0], Vec::new());
        assert!(
            matches!(bounded(&never_allocates), Err(OpError::Domain(msg)) if msg.contains("nested tensor"))
        );
        let bad_shape = tensor(vec![2], Vec::new());
        let long = Value::String("a".repeat(1_048_577));
        let first_error = Value::Variant(VariantValue {
            tag: "Test".into(),
            payload: vec![bad_shape, long],
        });
        assert!(
            matches!(bounded(&first_error), Err(OpError::Domain(msg)) if msg.contains("text budget"))
        );
    }

    #[test]
    fn depth_is_bounded_even_for_wrappers() {
        let nest = |levels| {
            (0..levels).fold(Value::Bool(true), |value, _| {
                Value::Dyn(DynValue {
                    dim: DimExpr::Dimless,
                    value: Box::new(value),
                })
            })
        };
        assert!(bounded(&nest(32)).is_ok());
        assert!(
            matches!(bounded(&nest(33)), Err(OpError::Domain(msg)) if msg.contains("element/depth"))
        );
    }
}
