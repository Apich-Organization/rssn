#![allow(clippy::wildcard_imports)] // Checked long-range and fluid native particle algorithms.
use super::engine::Field;
use super::engine::System;
use super::engine::{
    self,
};
mod fluids;
mod rigid;
use super::*;
use crate::rules::numerical_abi::codecs::Arg;
use crate::sim::ewald;
use crate::sim::fmm;
use crate::sim::md_engine as e;
use crate::sim::nbody;
struct Tree<const D: usize> {
    native: nbody::BhTree<D>,
    frame: Frame,
    source_dimension: DimExpr,
    source_reference: f64,
    n: usize,
}
impl<const D: usize> NativeType for Tree<D> {
    const KIND: &'static str = match D {
        | 1 => "rssn/native/sim/bh-tree/d1/v1",
        | 2 => "rssn/native/sim/bh-tree/d2/v1",
        | _ => "rssn/native/sim/bh-tree/d3/v1",
    };
    const MAX_BYTES: usize = 4_194_304;

    fn validate(&self) -> Result<(), OpError> {
        engine::bound(self.n, 512)?;
        self.frame.validate()
    }
}
fn source_values(
    f: &Frame,
    v: &Input,
    reference: &Input,
) -> Result<(Vec<f64>, f64, DimExpr), OpError> {
    let (r, d) = api::reference(&reference.0)?;
    api::same_dimension(&api::dimension(&v.0)?, &d)?;
    let values = <Vec<f64> as Arg>::decode(&v.0, 0)?
        .into_iter()
        .map(|x| x / r)
        .collect::<Vec<_>>();
    engine::finite(values.iter().copied())?;
    let _ = f;
    Ok((values, r, d))
}
fn source_field<const D: usize>(
    f: &Frame,
    data: Vec<[f64; D]>,
    r: f64,
    d: DimExpr,
    power: usize,
) -> Result<Value, OpError> {
    let dim = DimExpr::Mul(vec![d, f.dim([0, -(power as i32), 0, 0])?]);
    let dim = crate::rules::numerical_abi::runtime_dimensions::expression(
        &crate::rules::numerical_abi::runtime_dimensions::canonical(&dim)?,
    )?;
    let den = vec![f.values[1]; power];
    let n = data.len();
    let data = data
        .into_iter()
        .flatten()
        .map(|x| api::quantity(api::scale_product(x, &[r], &den)?, dim.clone()))
        .collect::<Result<Vec<_>, _>>()?;
    Ok(Value::Tensor(TensorValue {
        shape: vec![n, D],
        data,
        elem: TypeExpr::Float.with_dim(dim),
    }))
}
fn source_potential(
    f: &Frame,
    data: Vec<f64>,
    r: f64,
    d: DimExpr,
) -> Result<Value, OpError> {
    let dim = DimExpr::Mul(vec![d, DimExpr::Inverse(Box::new(f.dim(LENGTH)?))]);
    let dim = crate::rules::numerical_abi::runtime_dimensions::expression(
        &crate::rules::numerical_abi::runtime_dimensions::canonical(&dim)?,
    )?;
    Ok(Value::Tensor(TensorValue {
        shape: vec![data.len()],
        data: data
            .into_iter()
            .map(|x| api::quantity(api::scale_product(x, &[r], &[f.values[1]])?, dim.clone()))
            .collect::<Result<_, _>>()?,
        elem: TypeExpr::Float.with_dim(dim),
    }))
}
fn source_check<const D: usize>(
    p: &[[f64; D]],
    s: &[f64],
    eps: f64,
) -> Result<(), OpError> {
    if p.is_empty() || p.len() != s.len() || p.len() > 512 || eps < 0.0 {
        return Err(api::error(
            "long-range sources require equal nonempty<=512data and nonnegative softening",
        ));
    }
    engine::finite(p.iter().flatten().copied())?;
    engine::finite(s.iter().copied())?;
    if eps == 0.0 {
        let mut seen = std::collections::BTreeSet::new();
        for point in p {
            let key = point
                .iter()
                .map(|v| if *v == 0.0 { 0 } else { v.to_bits() })
                .collect::<Vec<_>>();
            if !seen.insert(key) {
                return Err(api::error("unsoftened source positions must be distinct"));
            }
        }
    }
    Ok(())
}
fn fmm_options(
    f: &Frame,
    v: &Input,
) -> Result<fmm::FmmOptions, OpError> {
    let order = <usize as Arg>::decode(api::field(&v.0, "order")?, 0)?;
    let leaf_size = <usize as Arg>::decode(api::field(&v.0, "leaf_size")?, 0)?;
    api::normalized(api::field(&v.0, "theta")?)?;
    let theta = crate::rules::numerical_abi::numeric(api::field(&v.0, "theta")?, 0)?;
    let eps = f.scalar(api::field(&v.0, "softening")?, LENGTH, false)?;
    if order == 0
        || order > 12
        || leaf_size == 0
        || leaf_size > 512
        || theta <= 0.0
        || theta >= 1.0
        || eps < 0.0
    {
        return Err(api::error(
            "FMM expansion/leaf/angle/softening controls exceedwork/domain bounds",
        ));
    }
    Ok(fmm::FmmOptions {
        order,
        leaf_size,
        theta,
        eps,
    })
}
fn recip(
    f: &Frame,
    out: ewald::Recip,
) -> Result<Dynamic, OpError> {
    Ok(Dynamic(api::record(vec![
        ("energy", f.output(out.energy, ENERGY, false)?),
        ("forces", engine::packed(f, out.forces, FORCE)?),
        ("virial", engine::packed(f, out.virial.to_vec(), ENERGY)?),
    ])))
}
fn wave_budget(
    n: usize,
    k: i32,
) -> Result<usize, OpError> {
    if !(0..=16).contains(&k) {
        return Err(api::error("reciprocal truncation kmax must be0..16"));
    }
    let work = (2 * k as usize + 1).pow(3) * n;
    if work > 50_000_000 {
        return Err(api::error("reciprocal wave-particle work budget exceeded"));
    }
    Ok(work)
}
fn mesh_budget(
    grid: [usize; 3],
    order: usize,
) -> Result<usize, OpError> {
    if !(2..=12).contains(&order)
        || !order.is_multiple_of(2)
        || grid.iter().any(|n| *n < order || *n > 128)
    {
        return Err(api::error(
            "PME/P3M need even order2..12 and gridorder..128",
        ));
    }
    let total = grid
        .iter()
        .try_fold(1usize, |n, g| n.checked_mul(*g))
        .filter(|n| *n <= 262144)
        .ok_or_else(|| api::error("reciprocal mesh allocation budget exceeded"))?;
    Ok(total * 125)
}
fn electrostatic(
    f: &Frame,
    alpha: &Input,
    constant: &Input,
) -> Result<(f64, f64), OpError> {
    let a = f.scalar(&alpha.0, [0, -1, 0, 0], false)?;
    let k = f.scalar(&constant.0, [1, 3, -2, -2], false)?;
    if a <= 0.0 {
        return Err(api::error("Ewald splitting alpha must bepositive"));
    }
    Ok((a, k))
}
struct Extra {
    native: Option<e::ExtraForce<3>>,
    frame: Frame,
    work: usize,
}
impl NativeType for Extra {
    const KIND: &'static str = "rssn/native/sim/electrostatic-extra/v1";
    const MAX_BYTES: usize = 4096;

    fn validate(&self) -> Result<(), OpError> {
        if self.native.is_none() {
            return Err(api::error("extra force ownership has been transferred"));
        }
        self.frame.validate()
    }
}
object_binding!(nbody_tree_new_d1,"num_native_nbody_tree_new_d1","crate::sim::nbody::BhTree<D>::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,positions:Input,sources:Input,source_reference:Input,theta:f64,softening:Input)->NewObject<Tree<1>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,positions:Input,sources:Input,source_reference:Input,theta:f64,softening:Input|{let f=frame(ctx,&f)?;let p=engine::rows::<1>(&f,&positions.0,LENGTH)?;let(src,r,d)=source_values(&f,&sources,&source_reference)?;let eps=f.scalar(&softening.0,LENGTH,false)?;source_check(&p,&src,eps)?;if theta<=0.0||theta>1.0{return Err(api::error("opening angle must be in(0,1]"));}let native=nbody::BhTree::new(&p,&src,theta,eps).map_err(api::error)?;let mut dims=f.dims.to_vec();dims.push(d.clone());Ok(NewObject::new(Tree{native,frame:f,source_reference:r,source_dimension:d,n:p.len()},dims))});
object_binding!(nbody_tree_field_d1,"num_native_nbody_tree_field_d1","crate::sim::nbody::BhTree<D>::field","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(tree:Handle<Tree<1>>,point:Input,skip:Option<usize>)->Dynamic,|ctx:&CallCtx<'_>,tree:Handle<Tree<1>>,point:Input,skip:Option<usize>|{tree.with(ctx,|t|{let x=engine::array::<1>(&t.frame,&point.0,LENGTH)?;if skip.is_some_and(|i|i>=t.n){return Err(api::error("skipped source index outside tree"));}Ok(Dynamic(source_field(&t.frame,vec![t.native.field(&x,skip.unwrap_or(usize::MAX))],t.source_reference,t.source_dimension.clone(),0)?))})});
object_binding!(nbody_direct_field_d1,"num_native_nbody_direct_field_d1","crate::sim::nbody::direct_field","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,positions:Input,sources:Input,source_reference:Input,softening:Input)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,positions:Input,sources:Input,source_reference:Input,softening:Input|{let f=frame(ctx,&f)?;let p=engine::rows::<1>(&f,&positions.0,LENGTH)?;let(src,r,d)=source_values(&f,&sources,&source_reference)?;let eps=f.scalar(&softening.0,LENGTH,false)?;source_check(&p,&src,eps)?;Ok(Dynamic(source_field(&f,nbody::direct_field(&p,&src,eps),r,d,0)?))});
object_binding!(nbody_gravity_accelerations_d1,"num_native_nbody_gravity_accelerations_d1","crate::sim::nbody::gravity_accelerations","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,positions:Input,masses:Input,constant:Input,softening:Input,theta:f64)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,positions:Input,masses:Input,constant:Input,softening:Input,theta:f64|{let f=frame(ctx,&f)?;let p=engine::rows::<1>(&f,&positions.0,LENGTH)?;let q=f.vector(&masses.0,MASS)?;let eps=f.scalar(&softening.0,LENGTH,false)?;let k=f.scalar(&constant.0,[-1,1,-2,0],false)?;source_check(&p,&q,eps)?;if theta<=0.0||theta>1.0||q.iter().any(|m|*m<=0.0){return Err(api::error("opening angle/masses outside domain"));}let out=nbody::gravity_accelerations(&p,&q,k,eps,theta).map_err(api::error)?;Ok(Dynamic(engine::packed(&f,out,[0,1,-2,0])?))});
object_binding!(nbody_coulomb_forces_d1,"num_native_nbody_coulomb_forces_d1","crate::sim::nbody::coulomb_forces","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,positions:Input,charges:Input,constant:Input,softening:Input,theta:f64)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,positions:Input,charges:Input,constant:Input,softening:Input,theta:f64|{let f=frame(ctx,&f)?;let p=engine::rows::<1>(&f,&positions.0,LENGTH)?;let q=f.vector(&charges.0,CHARGE)?;let eps=f.scalar(&softening.0,LENGTH,false)?;let k=f.scalar(&constant.0,[1,1,-2,-2],false)?;source_check(&p,&q,eps)?;if theta<=0.0||theta>1.0{return Err(api::error("opening angle/masses outside domain"));}let out=nbody::coulomb_forces(&p,&q,k,eps,theta).map_err(api::error)?;Ok(Dynamic(engine::packed(&f,out,FORCE)?))});
object_binding!(nbody_tree_new_d2,"num_native_nbody_tree_new_d2","crate::sim::nbody::BhTree<D>::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,positions:Input,sources:Input,source_reference:Input,theta:f64,softening:Input)->NewObject<Tree<2>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,positions:Input,sources:Input,source_reference:Input,theta:f64,softening:Input|{let f=frame(ctx,&f)?;let p=engine::rows::<2>(&f,&positions.0,LENGTH)?;let(src,r,d)=source_values(&f,&sources,&source_reference)?;let eps=f.scalar(&softening.0,LENGTH,false)?;source_check(&p,&src,eps)?;if theta<=0.0||theta>1.0{return Err(api::error("opening angle must be in(0,1]"));}let native=nbody::BhTree::new(&p,&src,theta,eps).map_err(api::error)?;let mut dims=f.dims.to_vec();dims.push(d.clone());Ok(NewObject::new(Tree{native,frame:f,source_reference:r,source_dimension:d,n:p.len()},dims))});
object_binding!(nbody_tree_field_d2,"num_native_nbody_tree_field_d2","crate::sim::nbody::BhTree<D>::field","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(tree:Handle<Tree<2>>,point:Input,skip:Option<usize>)->Dynamic,|ctx:&CallCtx<'_>,tree:Handle<Tree<2>>,point:Input,skip:Option<usize>|{tree.with(ctx,|t|{let x=engine::array::<2>(&t.frame,&point.0,LENGTH)?;if skip.is_some_and(|i|i>=t.n){return Err(api::error("skipped source index outside tree"));}Ok(Dynamic(source_field(&t.frame,vec![t.native.field(&x,skip.unwrap_or(usize::MAX))],t.source_reference,t.source_dimension.clone(),1)?))})});
object_binding!(nbody_direct_field_d2,"num_native_nbody_direct_field_d2","crate::sim::nbody::direct_field","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,positions:Input,sources:Input,source_reference:Input,softening:Input)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,positions:Input,sources:Input,source_reference:Input,softening:Input|{let f=frame(ctx,&f)?;let p=engine::rows::<2>(&f,&positions.0,LENGTH)?;let(src,r,d)=source_values(&f,&sources,&source_reference)?;let eps=f.scalar(&softening.0,LENGTH,false)?;source_check(&p,&src,eps)?;Ok(Dynamic(source_field(&f,nbody::direct_field(&p,&src,eps),r,d,1)?))});
object_binding!(nbody_gravity_accelerations_d2,"num_native_nbody_gravity_accelerations_d2","crate::sim::nbody::gravity_accelerations","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,positions:Input,masses:Input,constant:Input,softening:Input,theta:f64)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,positions:Input,masses:Input,constant:Input,softening:Input,theta:f64|{let f=frame(ctx,&f)?;let p=engine::rows::<2>(&f,&positions.0,LENGTH)?;let q=f.vector(&masses.0,MASS)?;let eps=f.scalar(&softening.0,LENGTH,false)?;let k=f.scalar(&constant.0,[-1,2,-2,0],false)?;source_check(&p,&q,eps)?;if theta<=0.0||theta>1.0||q.iter().any(|m|*m<=0.0){return Err(api::error("opening angle/masses outside domain"));}let out=nbody::gravity_accelerations(&p,&q,k,eps,theta).map_err(api::error)?;Ok(Dynamic(engine::packed(&f,out,[0,1,-2,0])?))});
object_binding!(nbody_coulomb_forces_d2,"num_native_nbody_coulomb_forces_d2","crate::sim::nbody::coulomb_forces","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,positions:Input,charges:Input,constant:Input,softening:Input,theta:f64)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,positions:Input,charges:Input,constant:Input,softening:Input,theta:f64|{let f=frame(ctx,&f)?;let p=engine::rows::<2>(&f,&positions.0,LENGTH)?;let q=f.vector(&charges.0,CHARGE)?;let eps=f.scalar(&softening.0,LENGTH,false)?;let k=f.scalar(&constant.0,[1,2,-2,-2],false)?;source_check(&p,&q,eps)?;if theta<=0.0||theta>1.0{return Err(api::error("opening angle/masses outside domain"));}let out=nbody::coulomb_forces(&p,&q,k,eps,theta).map_err(api::error)?;Ok(Dynamic(engine::packed(&f,out,FORCE)?))});
object_binding!(nbody_tree_new_d3,"num_native_nbody_tree_new_d3","crate::sim::nbody::BhTree<D>::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,positions:Input,sources:Input,source_reference:Input,theta:f64,softening:Input)->NewObject<Tree<3>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,positions:Input,sources:Input,source_reference:Input,theta:f64,softening:Input|{let f=frame(ctx,&f)?;let p=engine::rows::<3>(&f,&positions.0,LENGTH)?;let(src,r,d)=source_values(&f,&sources,&source_reference)?;let eps=f.scalar(&softening.0,LENGTH,false)?;source_check(&p,&src,eps)?;if theta<=0.0||theta>1.0{return Err(api::error("opening angle must be in(0,1]"));}let native=nbody::BhTree::new(&p,&src,theta,eps).map_err(api::error)?;let mut dims=f.dims.to_vec();dims.push(d.clone());Ok(NewObject::new(Tree{native,frame:f,source_reference:r,source_dimension:d,n:p.len()},dims))});
object_binding!(nbody_tree_field_d3,"num_native_nbody_tree_field_d3","crate::sim::nbody::BhTree<D>::field","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(tree:Handle<Tree<3>>,point:Input,skip:Option<usize>)->Dynamic,|ctx:&CallCtx<'_>,tree:Handle<Tree<3>>,point:Input,skip:Option<usize>|{tree.with(ctx,|t|{let x=engine::array::<3>(&t.frame,&point.0,LENGTH)?;if skip.is_some_and(|i|i>=t.n){return Err(api::error("skipped source index outside tree"));}Ok(Dynamic(source_field(&t.frame,vec![t.native.field(&x,skip.unwrap_or(usize::MAX))],t.source_reference,t.source_dimension.clone(),2)?))})});
object_binding!(nbody_direct_field_d3,"num_native_nbody_direct_field_d3","crate::sim::nbody::direct_field","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,positions:Input,sources:Input,source_reference:Input,softening:Input)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,positions:Input,sources:Input,source_reference:Input,softening:Input|{let f=frame(ctx,&f)?;let p=engine::rows::<3>(&f,&positions.0,LENGTH)?;let(src,r,d)=source_values(&f,&sources,&source_reference)?;let eps=f.scalar(&softening.0,LENGTH,false)?;source_check(&p,&src,eps)?;Ok(Dynamic(source_field(&f,nbody::direct_field(&p,&src,eps),r,d,2)?))});
object_binding!(nbody_gravity_accelerations_d3,"num_native_nbody_gravity_accelerations_d3","crate::sim::nbody::gravity_accelerations","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,positions:Input,masses:Input,constant:Input,softening:Input,theta:f64)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,positions:Input,masses:Input,constant:Input,softening:Input,theta:f64|{let f=frame(ctx,&f)?;let p=engine::rows::<3>(&f,&positions.0,LENGTH)?;let q=f.vector(&masses.0,MASS)?;let eps=f.scalar(&softening.0,LENGTH,false)?;let k=f.scalar(&constant.0,[-1,3,-2,0],false)?;source_check(&p,&q,eps)?;if theta<=0.0||theta>1.0||q.iter().any(|m|*m<=0.0){return Err(api::error("opening angle/masses outside domain"));}let out=nbody::gravity_accelerations(&p,&q,k,eps,theta).map_err(api::error)?;Ok(Dynamic(engine::packed(&f,out,[0,1,-2,0])?))});
object_binding!(nbody_coulomb_forces_d3,"num_native_nbody_coulomb_forces_d3","crate::sim::nbody::coulomb_forces","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,positions:Input,charges:Input,constant:Input,softening:Input,theta:f64)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,positions:Input,charges:Input,constant:Input,softening:Input,theta:f64|{let f=frame(ctx,&f)?;let p=engine::rows::<3>(&f,&positions.0,LENGTH)?;let q=f.vector(&charges.0,CHARGE)?;let eps=f.scalar(&softening.0,LENGTH,false)?;let k=f.scalar(&constant.0,[1,3,-2,-2],false)?;source_check(&p,&q,eps)?;if theta<=0.0||theta>1.0{return Err(api::error("opening angle/masses outside domain"));}let out=nbody::coulomb_forces(&p,&q,k,eps,theta).map_err(api::error)?;Ok(Dynamic(engine::packed(&f,out,FORCE)?))});
object_binding!(nbody_fmm2d_field,"num_native_nbody_fmm2d_field","crate::sim::nbody::fmm2d_field","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,positions:Input,sources:Input,source_reference:Input,order:usize,points_per_leaf:usize)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,positions:Input,sources:Input,source_reference:Input,order:usize,points_per_leaf:usize|{let f=frame(ctx,&f)?;let p=engine::rows::<2>(&f,&positions.0,LENGTH)?;let(src,r,d)=source_values(&f,&sources,&source_reference)?;source_check(&p,&src,0.0)?;if order==0||order>12||points_per_leaf==0||points_per_leaf>512{return Err(api::error("FMM order/leaf budget exceeded"));}let out=nbody::fmm2d_field(&p,&src,order,points_per_leaf).map_err(api::error)?;Ok(Dynamic(source_field(&f,out,r,d,1)?))});
object_binding!(nbody_fmm3d_field,"num_native_nbody_fmm3d_field","crate::sim::nbody::fmm3d_field","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,positions:Input,sources:Input,source_reference:Input,order:usize,points_per_leaf:usize)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,positions:Input,sources:Input,source_reference:Input,order:usize,points_per_leaf:usize|{let f=frame(ctx,&f)?;let p=engine::rows::<3>(&f,&positions.0,LENGTH)?;let(src,r,d)=source_values(&f,&sources,&source_reference)?;source_check(&p,&src,0.0)?;if order==0||order>12||points_per_leaf==0||points_per_leaf>512{return Err(api::error("FMM order/leaf budget exceeded"));}let out=nbody::fmm3d_field(&p,&src,order,points_per_leaf).map_err(api::error)?;Ok(Dynamic(source_field(&f,out,r,d,2)?))});
object_binding!(fmm_fmm2d,"num_native_fmm_fmm2d","crate::sim::fmm::fmm2d","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,positions:Input,sources:Input,source_reference:Input,options:Input)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,positions:Input,sources:Input,source_reference:Input,options:Input|{let f=frame(ctx,&f)?;let p=engine::rows::<2>(&f,&positions.0,LENGTH)?;let(src,r,d)=source_values(&f,&sources,&source_reference)?;let o=fmm_options(&f,&options)?;source_check(&p,&src,o.eps)?;let out=fmm::fmm2d(&p,&src,&o).map_err(api::error)?;Ok(Dynamic(source_field(&f,out,r,d,1)?))});
object_binding!(fmm_fmm3d,"num_native_fmm_fmm3d","crate::sim::fmm::fmm3d","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,positions:Input,sources:Input,source_reference:Input,options:Input)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,positions:Input,sources:Input,source_reference:Input,options:Input|{let f=frame(ctx,&f)?;let p=engine::rows::<3>(&f,&positions.0,LENGTH)?;let(src,r,d)=source_values(&f,&sources,&source_reference)?;let o=fmm_options(&f,&options)?;source_check(&p,&src,o.eps)?;let out=fmm::fmm3d(&p,&src,&o).map_err(api::error)?;Ok(Dynamic(api::record(vec![("potential",source_potential(&f,out.potential,r,d.clone())?),("field",source_field(&f,out.field,r,d,2)?)])))});
object_binding!(fmm_periodic_fmm3d,"num_native_fmm_periodic_fmm3d","crate::sim::fmm::periodic_fmm3d","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,positions:Input,sources:Input,source_reference:Input,options:Input,cell:Input,near_shells:usize,far_half_order:usize,soften_far:bool)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,positions:Input,sources:Input,source_reference:Input,options:Input,cell:Input,near_shells:usize,far_half_order:usize,soften_far:bool|{let f=frame(ctx,&f)?;let p=engine::rows::<3>(&f,&positions.0,LENGTH)?;let(src,r,d)=source_values(&f,&sources,&source_reference)?;let o=fmm_options(&f,&options)?;source_check(&p,&src,o.eps)?;if near_shells>1||far_half_order==0||far_half_order>8||p.len()>128{return Err(api::error("periodic FMM shell/order/source work bound exceeded"));}let h=engine::matrix::<3>(&f,&cell.0,LENGTH)?;e::SimBox::triclinic(h).map_err(api::error)?;let o=fmm::PeriodicFmmOptions{fmm:o,near_shells,far_half_order,soften_far};let out=fmm::periodic_fmm3d(&p,&src,&h,&o).map_err(api::error)?;Ok(Dynamic(api::record(vec![("potential",source_potential(&f,out.potential,r,d.clone())?),("field",source_field(&f,out.field,r,d,2)?)])))});
object_binding!(ewald_bspline,"num_native_ewald_bspline","crate::sim::ewald::bspline","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(order:usize,x:f64)->f64,|_ctx:&CallCtx<'_>,order:usize,x:f64|{if order==0||order>12{return Err(api::error("B-spline recursive order must be1..12"));}Ok(ewald::bspline(order,x))});
object_binding!(ewald_reciprocal_basis,"num_native_ewald_reciprocal_basis","crate::sim::ewald::reciprocal_basis","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,cell:Input)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,cell:Input|{let f=frame(ctx,&f)?;let h=engine::matrix::<3>(&f,&cell.0,LENGTH)?;let out=ewald::reciprocal_basis(&h).map_err(api::error)?;Ok(Dynamic(engine::packed(&f,out.to_vec(),[0,-1,0,0])?))});
object_binding!(ewald_real_space_pair,"num_native_ewald_real_space_pair","crate::sim::ewald::real_space_pair","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,distance:Input,charge_product:Input,alpha:Input,constant:Input)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,distance:Input,charge_product:Input,alpha:Input,constant:Input|{let f=frame(ctx,&f)?;let(a,k)=electrostatic(&f,&alpha,&constant)?;let r=f.scalar(&distance.0,LENGTH,false)?;let q=f.scalar(&charge_product.0,[0,0,0,2],false)?;if r<=0.0{return Err(api::error("pair distance must bepositive"));}Ok(Dynamic(f.output(ewald::real_space_pair(r,q,a,k),ENERGY,false)?))});
object_binding!(ewald_ewald_reciprocal,"num_native_ewald_ewald_reciprocal","crate::sim::ewald::ewald_reciprocal","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(system:Handle<System<3>>,alpha:Input,constant:Input,kmax:i32)->Dynamic,|ctx:&CallCtx<'_>,system:Handle<System<3>>,alpha:Input,constant:Input,kmax:i32|{system.with(ctx,|s|{let f=&s.frame;let(a,k)=electrostatic(f,&alpha,&constant)?;wave_budget(s.native.n(),kmax)?;if s.native.bx.cell.is_some(){return Err(api::error("orthorhombic API requires an orthorhombic box; use cell method"));}let out=ewald::ewald_reciprocal(&s.native.pos,&s.native.charge,s.native.bx.len,a,kmax,k);recip(f,out)})});
object_binding!(ewald_ewald_reciprocal_cell,"num_native_ewald_ewald_reciprocal_cell","crate::sim::ewald::ewald_reciprocal_cell","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(system:Handle<System<3>>,alpha:Input,constant:Input,kmax:i32)->Dynamic,|ctx:&CallCtx<'_>,system:Handle<System<3>>,alpha:Input,constant:Input,kmax:i32|{system.with(ctx,|s|{let f=&s.frame;let(a,k)=electrostatic(f,&alpha,&constant)?;wave_budget(s.native.n(),kmax)?;let out=ewald::ewald_reciprocal_cell(&s.native.pos,&s.native.charge,&s.native.bx.matrix(),a,kmax,k).map_err(api::error)?;recip(f,out)})});
object_binding!(ewald_ewald_energy_forces,"num_native_ewald_ewald_energy_forces","crate::sim::ewald::ewald_energy_forces","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(system:Handle<System<3>>,alpha:Input,constant:Input,kmax:i32,cutoff:Input)->Dynamic,|ctx:&CallCtx<'_>,system:Handle<System<3>>,alpha:Input,constant:Input,kmax:i32,cutoff:Input|{system.with(ctx,|s|{let f=&s.frame;let(a,k)=electrostatic(f,&alpha,&constant)?;wave_budget(s.native.n(),kmax)?;if s.native.bx.cell.is_some(){return Err(api::error("orthorhombic API requires an orthorhombic box; use cell method"));}let rcut=f.scalar(&cutoff.0,LENGTH,false)?;if rcut<=0.0{return Err(api::error("realspace cutoff must bepositive"));}
if s.native.bx.len.iter().any(|l|2.0*rcut>*l){return Err(api::error("minimum image Ewald cutoff must be<=halfbox"));}let out=ewald::ewald_energy_forces(&s.native.pos,&s.native.charge,s.native.bx.len,a,kmax,rcut,k);let(energy,forces)=out;Ok(Dynamic(api::record(vec![("energy",f.output(energy,ENERGY,false)?),("forces",engine::packed(f,forces,FORCE)?)])))})});
object_binding!(ewald_ewald_energy_forces_cell,"num_native_ewald_ewald_energy_forces_cell","crate::sim::ewald::ewald_energy_forces_cell","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(system:Handle<System<3>>,alpha:Input,constant:Input,kmax:i32,cutoff:Input)->Dynamic,|ctx:&CallCtx<'_>,system:Handle<System<3>>,alpha:Input,constant:Input,kmax:i32,cutoff:Input|{system.with(ctx,|s|{let f=&s.frame;let(a,k)=electrostatic(f,&alpha,&constant)?;wave_budget(s.native.n(),kmax)?;let rcut=f.scalar(&cutoff.0,LENGTH,false)?;if rcut<=0.0{return Err(api::error("realspace cutoff must bepositive"));}let images=(0..3).try_fold(1usize,|n,axis|{let layers=(rcut/s.native.bx.height(axis)).ceil()+1.0;if layers>16.0{return None;}n.checked_mul(2*layers as usize+1)}).ok_or_else(||api::error("Ewald image budget exceeded"))?;if images.checked_mul(s.native.n()*s.native.n()).is_none_or(|n|n>50_000_000){return Err(api::error("Ewald image pair work bound exceeded"));}let out=ewald::ewald_energy_forces_cell(&s.native.pos,&s.native.charge,&s.native.bx.matrix(),a,kmax,rcut,k).map_err(api::error)?;let(energy,forces)=out;Ok(Dynamic(api::record(vec![("energy",f.output(energy,ENERGY,false)?),("forces",engine::packed(f,forces,FORCE)?)])))})});
object_binding!(ewald_pme_reciprocal,"num_native_ewald_pme_reciprocal","crate::sim::ewald::pme_reciprocal","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(system:Handle<System<3>>,alpha:Input,constant:Input,grid:[usize;3],order:usize)->Dynamic,|ctx:&CallCtx<'_>,system:Handle<System<3>>,alpha:Input,constant:Input,grid:[usize;3],order:usize|{system.with(ctx,|s|{let f=&s.frame;let(a,k)=electrostatic(f,&alpha,&constant)?;mesh_budget(grid,order)?;if s.native.bx.cell.is_some(){return Err(api::error("orthorhombic API requires an orthorhombic box; use cell method"));}let out=ewald::pme_reciprocal(&s.native.pos,&s.native.charge,s.native.bx.len,a,grid,order,k);recip(f,out)})});
object_binding!(ewald_pme_reciprocal_cell,"num_native_ewald_pme_reciprocal_cell","crate::sim::ewald::pme_reciprocal_cell","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(system:Handle<System<3>>,alpha:Input,constant:Input,grid:[usize;3],order:usize)->Dynamic,|ctx:&CallCtx<'_>,system:Handle<System<3>>,alpha:Input,constant:Input,grid:[usize;3],order:usize|{system.with(ctx,|s|{let f=&s.frame;let(a,k)=electrostatic(f,&alpha,&constant)?;mesh_budget(grid,order)?;let out=ewald::pme_reciprocal_cell(&s.native.pos,&s.native.charge,&s.native.bx.matrix(),a,grid,order,k).map_err(api::error)?;recip(f,out)})});
object_binding!(ewald_p3m_reciprocal_cell,"num_native_ewald_p3m_reciprocal_cell","crate::sim::ewald::p3m_reciprocal_cell","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(system:Handle<System<3>>,alpha:Input,constant:Input,grid:[usize;3],order:usize)->Dynamic,|ctx:&CallCtx<'_>,system:Handle<System<3>>,alpha:Input,constant:Input,grid:[usize;3],order:usize|{system.with(ctx,|s|{let f=&s.frame;let(a,k)=electrostatic(f,&alpha,&constant)?;mesh_budget(grid,order)?;let out=ewald::p3m_reciprocal_cell(&s.native.pos,&s.native.charge,&s.native.bx.matrix(),a,grid,order,k).map_err(api::error)?;recip(f,out)})});
object_binding!(ewald_ewald_extra,"num_native_ewald_ewald_extra","crate::sim::ewald::ewald_extra","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,alpha:Input,constant:Input,kmax:i32)->NewObject<Extra>,|ctx:&CallCtx<'_>,f:Handle<Frame>,alpha:Input,constant:Input,kmax:i32|{let f=frame(ctx,&f)?;let(a,k)=electrostatic(&f,&alpha,&constant)?;let work=wave_budget(512,kmax)?;let native=ewald::ewald_extra(a,kmax,k);let dims=f.dims.to_vec();Ok(NewObject::new(Extra{native:Some(native),frame:f,work},dims))});
object_binding!(ewald_pme_extra,"num_native_ewald_pme_extra","crate::sim::ewald::pme_extra","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,alpha:Input,constant:Input,grid:[usize;3],order:usize)->NewObject<Extra>,|ctx:&CallCtx<'_>,f:Handle<Frame>,alpha:Input,constant:Input,grid:[usize;3],order:usize|{let f=frame(ctx,&f)?;let(a,k)=electrostatic(&f,&alpha,&constant)?;let work=mesh_budget(grid,order)?;let native=ewald::pme_extra(a,grid,order,k);let dims=f.dims.to_vec();Ok(NewObject::new(Extra{native:Some(native),frame:f,work},dims))});
object_binding!(ewald_p3m_extra,"num_native_ewald_p3m_extra","crate::sim::ewald::p3m_extra","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,alpha:Input,constant:Input,grid:[usize;3],order:usize)->NewObject<Extra>,|ctx:&CallCtx<'_>,f:Handle<Frame>,alpha:Input,constant:Input,grid:[usize;3],order:usize|{let f=frame(ctx,&f)?;let(a,k)=electrostatic(&f,&alpha,&constant)?;let work=mesh_budget(grid,order)?;let native=ewald::p3m_extra(a,grid,order,k);let dims=f.dims.to_vec();Ok(NewObject::new(Extra{native:Some(native),frame:f,work},dims))});
object_binding!(engine_forcefield_add_extra,"num_native_engine_forcefield_add_extra","","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(field:Handle<Field<3>>,extra:Handle<Extra>)->bool,|ctx:&CallCtx<'_>,field:Handle<Field<3>>,extra:Handle<Extra>|{field.with_mut(ctx,|f|extra.with_mut(ctx,|x|{f.frame.same(&x.frame)?;let ff=f.native.as_mut().ok_or_else(||api::error("force field consumed"))?;if ff.extra.len()>=8||f.extra_work.checked_add(x.work).is_none_or(|n|n>50_000_000){return Err(api::error("extra force callback count/work budget exceeded"));}let native=x.native.take().ok_or_else(||api::error("extra force already consumed"))?;ff.extra.push(native);f.extra_work+=x.work;Ok(true)}))});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    fluids::install(i)?;
    rigid::install(i)?;
    api::install(
        i,
        &[
            &nbody_tree_new_d1::BINDING,
            &nbody_tree_field_d1::BINDING,
            &nbody_direct_field_d1::BINDING,
            &nbody_gravity_accelerations_d1::BINDING,
            &nbody_coulomb_forces_d1::BINDING,
            &nbody_tree_new_d2::BINDING,
            &nbody_tree_field_d2::BINDING,
            &nbody_direct_field_d2::BINDING,
            &nbody_gravity_accelerations_d2::BINDING,
            &nbody_coulomb_forces_d2::BINDING,
            &nbody_tree_new_d3::BINDING,
            &nbody_tree_field_d3::BINDING,
            &nbody_direct_field_d3::BINDING,
            &nbody_gravity_accelerations_d3::BINDING,
            &nbody_coulomb_forces_d3::BINDING,
            &nbody_fmm2d_field::BINDING,
            &nbody_fmm3d_field::BINDING,
            &fmm_fmm2d::BINDING,
            &fmm_fmm3d::BINDING,
            &fmm_periodic_fmm3d::BINDING,
            &ewald_bspline::BINDING,
            &ewald_reciprocal_basis::BINDING,
            &ewald_real_space_pair::BINDING,
            &ewald_ewald_reciprocal::BINDING,
            &ewald_ewald_reciprocal_cell::BINDING,
            &ewald_ewald_energy_forces::BINDING,
            &ewald_ewald_energy_forces_cell::BINDING,
            &ewald_pme_reciprocal::BINDING,
            &ewald_pme_reciprocal_cell::BINDING,
            &ewald_p3m_reciprocal_cell::BINDING,
            &ewald_ewald_extra::BINDING,
            &ewald_pme_extra::BINDING,
            &ewald_p3m_extra::BINDING,
            &engine_forcefield_add_extra::BINDING,
        ],
        "numerical_vector",
    )
}
