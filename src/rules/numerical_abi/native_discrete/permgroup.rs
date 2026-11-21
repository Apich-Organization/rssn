#![allow(clippy::wildcard_imports)] // Checked native permutation-group methods.
use super::*;
use crate::kernel::value::CallbackRef;
use crate::kernels::permgroup::Group;
use crate::kernels::permgroup::{
    self,
};
use num_bigint::BigUint;
use rand_v09::SeedableRng;
use rand_v09::rngs::StdRng;
use std::cell::Cell;
use std::cell::RefCell;
impl NativeType for Group {
    const KIND: &'static str = "rssn/native/discrete/permutation-group/v1";
    const MAX_BYTES: usize = 8388608;

    fn validate(&self) -> Result<(), OpError> {
        limit(self.degree(), 16)?;
        limit(self.generators().len(), 256)?;
        limit(self.strong_generators().len(), 4096)?;
        Ok(())
    }
}
impl NativeType for StdRng {
    const KIND: &'static str = "rssn/native/discrete/std-rng/v1";
    const MAX_BYTES: usize = 512;

    fn validate(&self) -> Result<(), OpError> {
        Ok(())
    }
}
fn permutation(
    p: &[u32],
    n: usize,
) -> Result<(), OpError> {
    if p.len() != n || !permgroup::is_permutation(p) {
        return Err(error("group requires a valid permutation of its degree"));
    }
    Ok(())
}
fn generators(
    p: &[Vec<u32>],
    n: usize,
) -> Result<(), OpError> {
    limit(n, 16)?;
    limit(p.len(), 128)?;
    for a in p {
        permutation(a, n)?;
    }
    Ok(())
}
fn work(g: &Group) -> Result<(), OpError> {
    if g.order() > BigUint::from(10000u32) {
        return Err(error(
            "native group backtracking order exceeds10000 work bound",
        ));
    }
    Ok(())
}
fn points(
    p: &[u32],
    n: usize,
) -> Result<(), OpError> {
    if p.iter().any(|x| *x as usize >= n) {
        return Err(error("group point outside degree"));
    }
    let mut q = p.to_vec();
    q.sort_unstable();
    q.dedup();
    if q.len() != p.len() {
        return Err(error("group point prefix must contain distinct points"));
    }
    Ok(())
}
struct GroupProperty(CallbackRef);
struct GroupPrune(CallbackRef);
fn callback(
    value: &Value,
    ctx: &CallCtx<'_>,
    arity: usize,
) -> Result<CallbackRef, OpError> {
    let Value::Function(h) = value else {
        return Err(error("expected group callback"));
    };
    let s = ctx.callback_signature(*h)?;
    if s.params.len() != arity || s.ret != TypeExpr::Bool {
        return Err(error(
            "group callback must return Bool at its declared arity",
        ));
    }
    Ok(*h)
}
impl NativeArg for GroupProperty {
    fn ty() -> TypeExpr {
        TypeExpr::Function(
            vec![api::tensor_type(1, TypeExpr::Int)],
            Box::new(TypeExpr::Bool),
        )
    }

    fn decode(
        v: &Value,
        _: usize,
        c: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        Ok(Self(callback(v, c, 1)?))
    }
}
impl NativeArg for GroupPrune {
    fn ty() -> TypeExpr {
        TypeExpr::Function(vec![TypeExpr::Int, TypeExpr::Int], Box::new(TypeExpr::Bool))
    }

    fn decode(
        v: &Value,
        _: usize,
        c: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        Ok(Self(callback(v, c, 2)?))
    }
}
object_binding!(std_rng_new,"num_native_std_rng_new","","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(seed:u64)->NewObject<StdRng>,|_ctx:&CallCtx<'_>,seed:u64|{Ok(NewObject::normalized(StdRng::seed_from_u64(seed)))});
object_binding!(group_new,"num_native_group_new","crate::kernels::permgroup::Group::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(n:usize,gens:Vec<Vec<u32>>)->NewObject<Group>,|_ctx:&CallCtx<'_>,n:usize,gens:Vec<Vec<u32>>|{generators(&gens,n)?;Ok(NewObject::normalized(Group::new(n,&gens)))});
object_binding!(group_trivial,"num_native_group_trivial","crate::kernels::permgroup::Group::trivial","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(n:usize)->NewObject<Group>,|_ctx:&CallCtx<'_>,n:usize|{limit(n,16)?;Ok(NewObject::normalized(Group::trivial(n)))});
object_binding!(group_degree,"num_native_group_degree","crate::kernels::permgroup::Group::degree","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Group>)->usize,|ctx:&CallCtx<'_>,a:Handle<Group>|{a.with(ctx,|g|Ok(g.degree()))});
object_binding!(group_base,"num_native_group_base","crate::kernels::permgroup::Group::base","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Group>)->Vec<u32>,|ctx:&CallCtx<'_>,a:Handle<Group>|{a.with(ctx,|g|Ok(g.base()))});
object_binding!(group_basic_orbit_lengths,"num_native_group_basic_orbit_lengths","crate::kernels::permgroup::Group::basic_orbit_lengths","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Group>)->Vec<usize>,|ctx:&CallCtx<'_>,a:Handle<Group>|{a.with(ctx,|g|Ok(g.basic_orbit_lengths()))});
object_binding!(group_order,"num_native_group_order","crate::kernels::permgroup::Group::order","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Group>)->BigUint,|ctx:&CallCtx<'_>,a:Handle<Group>|{a.with(ctx,|g|Ok(g.order()))});
object_binding!(group_strong_generators,"num_native_group_strong_generators","crate::kernels::permgroup::Group::strong_generators","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Group>)->Vec<Vec<u32>>,|ctx:&CallCtx<'_>,a:Handle<Group>|{a.with(ctx,|g|Ok(g.strong_generators()))});
object_binding!(group_is_transitive,"num_native_group_is_transitive","crate::kernels::permgroup::Group::is_transitive","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Group>)->bool,|ctx:&CallCtx<'_>,a:Handle<Group>|{a.with(ctx,|g|Ok(g.is_transitive()))});
object_binding!(group_is_abelian,"num_native_group_is_abelian","crate::kernels::permgroup::Group::is_abelian","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Group>)->bool,|ctx:&CallCtx<'_>,a:Handle<Group>|{a.with(ctx,|g|Ok(g.is_abelian()))});
object_binding!(group_generators,"num_native_group_generators","crate::kernels::permgroup::Group::generators","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Group>)->Vec<Vec<u32>>,|ctx:&CallCtx<'_>,a:Handle<Group>|{a.with(ctx,|g|Ok(g.generators().to_vec()))});
object_binding!(group_orbits,"num_native_group_orbits","crate::kernels::permgroup::Group::orbits","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Group>)->Dynamic,|ctx:&CallCtx<'_>,a:Handle<Group>|{a.with(ctx,|g|{let mut offsets=vec![0usize];let mut values=vec![];for o in g.orbits(){values.extend(o);offsets.push(values.len());}Ok(Dynamic(api::record(vec![("offsets",Output::encode(offsets)?),("points",Output::encode(values)?)])))})});
object_binding!(group_contains,"num_native_group_contains","crate::kernels::permgroup::Group::contains","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Group>,permutation_value:Vec<u32>)->bool,|ctx:&CallCtx<'_>,a:Handle<Group>,permutation_value:Vec<u32>|{a.with(ctx,|g|{permutation(&permutation_value,g.degree())?;Ok(g.contains(&permutation_value))})});
object_binding!(group_extend,"num_native_group_extend","crate::kernels::permgroup::Group::extend","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Group>,permutation_value:Vec<u32>)->bool,|ctx:&CallCtx<'_>,a:Handle<Group>,permutation_value:Vec<u32>|{a.with_mut(ctx,|g|{permutation(&permutation_value,g.degree())?;if g.generators().len()>=128{return Err(error("group generator growth budget exceeded"));}g.extend(permutation_value);g.validate()?;Ok(true)})});
object_binding!(group_random_element,"num_native_group_random_element","crate::kernels::permgroup::Group::random_element","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Group>,rng:Handle<StdRng>)->Vec<u32>,|ctx:&CallCtx<'_>,a:Handle<Group>,rng:Handle<StdRng>|{a.with(ctx,|g|rng.with_mut(ctx,|r|Ok(g.random_element(r))))});
object_binding!(group_elements,"num_native_group_elements","crate::kernels::permgroup::Group::elements","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Group>,max_elements:usize)->Option<Vec<Vec<u32>>>,|ctx:&CallCtx<'_>,a:Handle<Group>,max_elements:usize|{limit(max_elements,10000)?;a.with(ctx,|g|Ok(g.elements(max_elements)))});
object_binding!(group_center,"num_native_group_center","crate::kernels::permgroup::Group::center","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Group>)->NewObject<Group>,|ctx:&CallCtx<'_>,a:Handle<Group>|{a.with(ctx,|g|{work(g)?;Ok(NewObject::normalized(g.center()))})});
object_binding!(group_derived_subgroup,"num_native_group_derived_subgroup","crate::kernels::permgroup::Group::derived_subgroup","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Group>)->NewObject<Group>,|ctx:&CallCtx<'_>,a:Handle<Group>|{a.with(ctx,|g|{work(g)?;Ok(NewObject::normalized(g.derived_subgroup()))})});
object_binding!(group_is_solvable,"num_native_group_is_solvable","crate::kernels::permgroup::Group::is_solvable","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Group>)->bool,|ctx:&CallCtx<'_>,a:Handle<Group>|{a.with(ctx,|g|{work(g)?;Ok(g.is_solvable())})});
object_binding!(group_solvable_composition_factors,"num_native_group_solvable_composition_factors","crate::kernels::permgroup::Group::solvable_composition_factors","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Group>)->Option<Vec<u64>>,|ctx:&CallCtx<'_>,a:Handle<Group>|{a.with(ctx,|g|{work(g)?;Ok(g.solvable_composition_factors())})});
object_binding!(group_contains_group,"num_native_group_contains_group","crate::kernels::permgroup::Group::contains_group","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Group>,b:Handle<Group>)->bool,|ctx:&CallCtx<'_>,a:Handle<Group>,b:Handle<Group>|{a.with(ctx,|g|b.with(ctx,|h|{if g.degree()!=h.degree(){return Err(error("group degrees differ"));}work(g)?;work(h)?;Ok(g.contains_group(h))}))});
object_binding!(group_normalizer,"num_native_group_normalizer","crate::kernels::permgroup::Group::normalizer","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Group>,b:Handle<Group>)->NewObject<Group>,|ctx:&CallCtx<'_>,a:Handle<Group>,b:Handle<Group>|{a.with(ctx,|g|b.with(ctx,|h|{if g.degree()!=h.degree(){return Err(error("group degrees differ"));}work(g)?;work(h)?;Ok(NewObject::normalized(g.normalizer(h)))}))});
object_binding!(group_intersection,"num_native_group_intersection","crate::kernels::permgroup::Group::intersection","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Group>,b:Handle<Group>)->NewObject<Group>,|ctx:&CallCtx<'_>,a:Handle<Group>,b:Handle<Group>|{a.with(ctx,|g|b.with(ctx,|h|{if g.degree()!=h.degree(){return Err(error("group degrees differ"));}work(g)?;work(h)?;Ok(NewObject::normalized(g.intersection(h)))}))});
object_binding!(group_is_normal,"num_native_group_is_normal","crate::kernels::permgroup::Group::is_normal","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Group>,b:Handle<Group>)->bool,|ctx:&CallCtx<'_>,a:Handle<Group>,b:Handle<Group>|{a.with(ctx,|g|b.with(ctx,|h|{if g.degree()!=h.degree(){return Err(error("group degrees differ"));}work(g)?;work(h)?;Ok(g.is_normal(h))}))});
object_binding!(group_pointwise_stabilizer,"num_native_group_pointwise_stabilizer","crate::kernels::permgroup::Group::pointwise_stabilizer","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Group>,selected:Vec<u32>)->NewObject<Group>,|ctx:&CallCtx<'_>,a:Handle<Group>,selected:Vec<u32>|{a.with(ctx,|g|{points(&selected,g.degree())?;work(g)?;Ok(NewObject::normalized(g.pointwise_stabilizer(&selected)))})});
object_binding!(group_with_base_prefix,"num_native_group_with_base_prefix","crate::kernels::permgroup::Group::with_base_prefix","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Group>,selected:Vec<u32>)->NewObject<Group>,|ctx:&CallCtx<'_>,a:Handle<Group>,selected:Vec<u32>|{a.with(ctx,|g|{points(&selected,g.degree())?;work(g)?;Ok(NewObject::normalized(g.with_base_prefix(&selected)))})});
object_binding!(group_set_stabilizer,"num_native_group_set_stabilizer","crate::kernels::permgroup::Group::set_stabilizer","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Group>,selected:Vec<u32>)->NewObject<Group>,|ctx:&CallCtx<'_>,a:Handle<Group>,selected:Vec<u32>|{a.with(ctx,|g|{points(&selected,g.degree())?;work(g)?;Ok(NewObject::normalized(g.set_stabilizer(&selected)))})});
object_binding!(group_centralizer,"num_native_group_centralizer","crate::kernels::permgroup::Group::centralizer","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Group>,x:Vec<u32>)->NewObject<Group>,|ctx:&CallCtx<'_>,a:Handle<Group>,x:Vec<u32>|{a.with(ctx,|g|{permutation(&x,g.degree())?;work(g)?;Ok(NewObject::normalized(g.centralizer(&x)))})});
object_binding!(group_conjugating_element,"num_native_group_conjugating_element","crate::kernels::permgroup::Group::conjugating_element","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Group>,x:Vec<u32>,y:Vec<u32>)->Option<Vec<u32>>,|ctx:&CallCtx<'_>,a:Handle<Group>,x:Vec<u32>,y:Vec<u32>|{a.with(ctx,|g|{permutation(&x,g.degree())?;permutation(&y,g.degree())?;work(g)?;Ok(g.conjugating_element(&x,&y))})});
object_binding!(group_are_conjugate,"num_native_group_are_conjugate","crate::kernels::permgroup::Group::are_conjugate","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Group>,x:Vec<u32>,y:Vec<u32>)->bool,|ctx:&CallCtx<'_>,a:Handle<Group>,x:Vec<u32>,y:Vec<u32>|{a.with(ctx,|g|{permutation(&x,g.degree())?;permutation(&y,g.degree())?;work(g)?;Ok(g.are_conjugate(&x,&y))})});
object_binding!(group_normal_closure,"num_native_group_normal_closure","crate::kernels::permgroup::Group::normal_closure","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Group>,gens:Vec<Vec<u32>>)->NewObject<Group>,|ctx:&CallCtx<'_>,a:Handle<Group>,gens:Vec<Vec<u32>>|{a.with(ctx,|g|{generators(&gens,g.degree())?;work(g)?;Ok(NewObject::normalized(g.normal_closure(&gens)))})});
object_binding!(group_sylow_subgroup,"num_native_group_sylow_subgroup","crate::kernels::permgroup::Group::sylow_subgroup","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Group>,prime:u64,seed:u64,max_tries:usize)->Maybe<NewObject<Group>>,|ctx:&CallCtx<'_>,a:Handle<Group>,prime:u64,seed:u64,max_tries:usize|{limit(max_tries,10000)?;if !crate::kernels::homology::is_prime(prime){return Err(error("Sylow prime must be prime >=2"));}a.with(ctx,|g|{work(g)?;Ok(Maybe(g.sylow_subgroup(prime,seed,max_tries).map(NewObject::normalized)))})});
object_binding!(group_exponent,"num_native_group_exponent","crate::kernels::permgroup::Group::exponent","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Group>,seed:u64)->Option<BigUint>,|ctx:&CallCtx<'_>,a:Handle<Group>,seed:u64|{a.with(ctx,|g|{if g.order()>BigUint::from(1000u32){return Err(error("group exponent class-search order budget is1000"));}Ok(g.exponent(seed))})});
object_binding!(group_conjugacy_classes,"num_native_group_conjugacy_classes","crate::kernels::permgroup::Group::conjugacy_classes","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Group>,seed:u64,max_tries:usize)->Maybe<Dynamic>,|ctx:&CallCtx<'_>,a:Handle<Group>,seed:u64,max_tries:usize|{limit(max_tries,10000)?;a.with(ctx,|g|{work(g)?;let out=g.conjugacy_classes(seed,max_tries).map(|classes|{let(rep,size):(Vec<_>,Vec<_>)=classes.into_iter().unzip();Ok::<_,OpError>(Dynamic(api::record(vec![("representatives",Output::encode(rep)?),("sizes",Output::encode(size.into_iter().map(num_bigint::BigInt::from).collect::<Vec<_>>())?)])))}).transpose()?;Ok(Maybe(out))})});
object_binding!(group_subgroup,"num_native_group_subgroup","crate::kernels::permgroup::Group::subgroup","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Group>,property:GroupProperty,prune:GroupPrune,callback_budget:usize)->NewObject<Group>,|ctx:&CallCtx<'_>,a:Handle<Group>,property:GroupProperty,prune:GroupPrune,callback_budget:usize|{limit(callback_budget,1_000_000)?;a.with(ctx,|g|{work(g)?;let remaining=Cell::new(callback_budget);let failure=RefCell::new(None);let invoke=|h,args:Vec<Value>|{if failure.borrow().is_some(){return false;}let Some(n)=remaining.get().checked_sub(1)else{*failure.borrow_mut()=Some(error("group callback budget exhausted"));return false;};remaining.set(n);match ctx.invoke_callback(h,&args){Ok(Value::Bool(v))=>v,Ok(_)=>{*failure.borrow_mut()=Some(error("group callback returned a non-Bool"));false},Err(e)=>{*failure.borrow_mut()=Some(e);false}}};let result=g.subgroup(&|p|{match Output::encode(p.clone()){Ok(v)=>invoke(property.0,vec![v]),Err(e)=>{*failure.borrow_mut()=Some(e);false}}},&|x,y|invoke(prune.0,vec![Value::int(x),Value::int(y)]));if let Some(e)=failure.into_inner(){return Err(e);}Ok(NewObject::normalized(result))})});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &std_rng_new::BINDING,
            &group_new::BINDING,
            &group_trivial::BINDING,
            &group_degree::BINDING,
            &group_base::BINDING,
            &group_basic_orbit_lengths::BINDING,
            &group_order::BINDING,
            &group_strong_generators::BINDING,
            &group_is_transitive::BINDING,
            &group_is_abelian::BINDING,
            &group_generators::BINDING,
            &group_orbits::BINDING,
            &group_contains::BINDING,
            &group_extend::BINDING,
            &group_random_element::BINDING,
            &group_elements::BINDING,
            &group_center::BINDING,
            &group_derived_subgroup::BINDING,
            &group_is_solvable::BINDING,
            &group_solvable_composition_factors::BINDING,
            &group_contains_group::BINDING,
            &group_normalizer::BINDING,
            &group_intersection::BINDING,
            &group_is_normal::BINDING,
            &group_pointwise_stabilizer::BINDING,
            &group_with_base_prefix::BINDING,
            &group_set_stabilizer::BINDING,
            &group_centralizer::BINDING,
            &group_conjugating_element::BINDING,
            &group_are_conjugate::BINDING,
            &group_normal_closure::BINDING,
            &group_sylow_subgroup::BINDING,
            &group_exponent::BINDING,
            &group_conjugacy_classes::BINDING,
            &group_subgroup::BINDING,
        ],
        "numerical_vector",
    )
}
