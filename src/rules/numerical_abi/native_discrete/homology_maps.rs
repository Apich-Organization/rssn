#![allow(clippy::wildcard_imports)] // Shared checked exact homology map adapters.
use super::qmatrix::QMatrix;
use super::*;
use crate::kernels::homology;
use crate::kernels::topology::Simplex;
impl NativeType for Vec<Simplex> {
    const KIND: &'static str = "rssn/native/homology/complex/v1";
    const MAX_BYTES: usize = 1_048_576;

    fn validate(&self) -> Result<(), OpError> {
        limit(self.len(), 4096)?;
        if self
            .iter()
            .any(|s| s.len() > 12 || s.iter().any(|v| *v > 1_000_000))
        {
            return Err(error("native complex exceeds simplex/vertex-label budget"));
        }
        Ok(())
    }
}
fn map_values(values: Vec<[usize; 2]>) -> Result<std::collections::HashMap<usize, usize>, OpError> {
    limit(values.len(), 4096)?;
    let mut result = std::collections::HashMap::new();
    for [a, b] in values {
        if result.insert(a, b).is_some() {
            return Err(error(
                "simplicial vertex map contains duplicate source vertices",
            ));
        }
    }
    Ok(result)
}
fn map_budget(
    k: &[Simplex],
    l: &[Simplex],
    dim: usize,
) -> Result<(), OpError> {
    limit(dim, 12)?;
    for complex in [k, l] {
        let mut counts = [0usize; 13];
        for s in complex {
            if !s.is_empty() {
                counts[s.len() - 1] += 1;
            }
        }
        if counts.iter().any(|n| *n > 128) {
            return Err(error(
                "induced homology exact chain rank exceeds elimination budget",
            ));
        }
    }
    Ok(())
}
object_binding!(homology_is_simplicial_map,"num_native_homology_is_simplicial_map","crate::kernels::homology::is_simplicial_map","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(k:Handle<Vec<Simplex>>,l:Handle<Vec<Simplex>>,map:Vec<[usize;2]>)->bool,|ctx:&CallCtx<'_>,k:Handle<Vec<Simplex>>,l:Handle<Vec<Simplex>>,map:Vec<[usize;2]>|{let map=map_values(map)?;k.with(ctx,|k|l.with(ctx,|l|Ok(homology::is_simplicial_map(k,l,&map))))});
object_binding!(homology_induced_chain_map,"num_native_homology_induced_chain_map","crate::kernels::homology::induced_chain_map","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(k:Handle<Vec<Simplex>>,l:Handle<Vec<Simplex>>,map:Vec<[usize;2]>,dimension:usize)->Option<Vec<Vec<i64>>>,|ctx:&CallCtx<'_>,k:Handle<Vec<Simplex>>,l:Handle<Vec<Simplex>>,map:Vec<[usize;2]>,dimension:usize|{let map=map_values(map)?;k.with(ctx,|k|l.with(ctx,|l|{map_budget(k,l,dimension)?;Ok(homology::induced_chain_map(k,l,&map,dimension))}))});
object_binding!(homology_induced_homology,"num_native_homology_induced_homology","crate::kernels::homology::induced_homology","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(k:Handle<Vec<Simplex>>,l:Handle<Vec<Simplex>>,map:Vec<[usize;2]>,dimension:usize)->Maybe<NewObject<QMatrix>>,|ctx:&CallCtx<'_>,k:Handle<Vec<Simplex>>,l:Handle<Vec<Simplex>>,map:Vec<[usize;2]>,dimension:usize|{let map=map_values(map)?;k.with(ctx,|k|l.with(ctx,|l|{map_budget(k,l,dimension)?;Ok(Maybe(homology::induced_homology(k,l,&map,dimension).map(|m|{let cols=m.first().map_or(0,Vec::len);NewObject::normalized(QMatrix::new(m,cols))})))}))});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &homology_is_simplicial_map::BINDING,
            &homology_induced_chain_map::BINDING,
            &homology_induced_homology::BINDING,
        ],
        "numerical_vector",
    )
}
