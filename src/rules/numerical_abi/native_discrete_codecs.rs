// Discrete group adapters validate permutation domains before native routines
// whose direct slice interfaces assume bijections and common degrees.
integer_output!(num_bigint::BigUint);
native_vector_codec!(num_bigint::BigInt);

fn permutation_rows(rows: Vec<Vec<u32>>) -> Vec<Vec<usize>> {
    rows.into_iter()
        .map(|row| row.into_iter().map(|v| v as usize).collect())
        .collect()
}
#[allow(clippy::indexing_slicing)] // Called only after exact CodecExport arity validation.
fn native_discrete_preflight(
    name: &str,
    args: &[Value],
) -> Result<(), OpError> {
    let Some(operation) = name.strip_prefix("num_permgroup_") else {
        return Ok(());
    };
    if operation == "identity" {
        native_budget(&[native_count(args, 0)?], 4096)?;
        return Ok(());
    }
    if matches!(
        operation,
        "orbits" | "minimal_blocks" | "is_primitive" | "minimal_block_systems"
    ) {
        let n = native_count(args, 0)?;
        native_budget(&[n], 4096)?;
        let generators = Vec::<Vec<u32>>::decode(&args[1], 1)?;
        if generators
            .iter()
            .any(|g| g.len() != n || !crate::kernels::permgroup::is_permutation(g))
        {
            return Err(OpError::Domain(
                "group generators must be permutations of the declared degree".into(),
            ));
        }
        native_budget(&[n, n, n, generators.len().max(1)], 50_000_000)?;
        if operation == "minimal_blocks"
            && (native_count(args, 2)? >= n || native_count(args, 3)? >= n)
        {
            return Err(OpError::Domain(
                "block seeds must lie inside the permutation degree".into(),
            ));
        }
        return Ok(());
    }
    let permutation = Vec::<u32>::decode(&args[0], 0)?;
    native_budget(&[permutation.len()], 4096)?;
    if matches!(operation, "is_permutation" | "power_exact") {
        return Ok(());
    }
    if !crate::kernels::permgroup::is_permutation(&permutation) {
        return Err(OpError::Domain(
            "native group operation requires a bijective permutation".into(),
        ));
    }
    if matches!(operation, "compose" | "conjugate") {
        let other = Vec::<u32>::decode(&args[1], 1)?;
        if other.len() != permutation.len() || !crate::kernels::permgroup::is_permutation(&other) {
            return Err(OpError::Domain(
                "permutation operands must be bijections of the same degree".into(),
            ));
        }
    }
    Ok(())
}

struct NativeBlockSystems {
    offsets: Vec<usize>,
    blocks: Vec<Vec<usize>>,
}
record!(NativeBlockSystems{offsets:Vec<usize>,blocks:Vec<Vec<usize>>});
fn block_systems(systems: Vec<Vec<Vec<u32>>>) -> Result<NativeBlockSystems, OpError> {
    let mut offsets = vec![0usize];
    let mut blocks = Vec::new();
    for system in systems {
        if blocks
            .len()
            .checked_add(system.len())
            .is_none_or(|n| n > NATIVE_VALUE_LIMIT)
        {
            return Err(OpError::Domain(
                "native block-system output exceeds its allocation budget".into(),
            ));
        }
        blocks.extend(permutation_rows(system));
        offsets.push(blocks.len());
    }
    Ok(NativeBlockSystems { offsets, blocks })
}
binding!(native_perm_cycles,"num_permgroup_cycles","numerical_vector","Cycles of a checked bijective permutation (fixed points included). CSR offsets/indices preserve independently sized cycles; degree is bounded by 4096.","crate::kernels::permgroup::cycles",DimensionPolicy::Dimensionless,
(p:Vec<u32>)->Vec<Vec<usize>>,|p:Vec<u32>|super::permutation_rows(crate::kernels::permgroup::cycles(&p)));
binding!(native_perm_orbits,"num_permgroup_orbits","numerical_vector","Orbits of generators encoded as rows of a common-degree permutation matrix. Output uses CSR offsets/indices with finite degree and group-work budgets.","crate::kernels::permgroup::orbits",DimensionPolicy::Dimensionless,
(n:usize,generators:Vec<Vec<u32>>)->Vec<Vec<usize>>,|n:usize,generators:Vec<Vec<u32>>|super::permutation_rows(crate::kernels::permgroup::orbits(n,&generators)));
binding!(
    native_perm_power_exact, "num_permgroup_power_exact", "numerical_vector",
    "Raise a dimensionless permutation of degree up to 4096 to an exact signed integer power. A negative exponent rotates cycles in reverse, and the magnitude can have at most 65536 bits. Invalid images return an error.",
    "crate::kernels::permgroup::power_exact", DimensionPolicy::Dimensionless,
    (p: Vec<u32>, exponent: num_bigint::BigInt) -> Result<crate::kernels::permgroup::Perm, OpError>,
    |p: Vec<u32>, exponent: num_bigint::BigInt| {
        if exponent.bits() > 65_536 {
            return Err(OpError::Domain("exact permutation exponent exceeds 65536 bits".into()));
        }
        crate::kernels::permgroup::power_exact(&p, &exponent)
            .ok_or_else(|| OpError::Domain("native group operation requires a bijective permutation".into()))
    }
);
binding!(native_perm_blocks,"num_permgroup_minimal_blocks","numerical_vector","Minimal block system identifying two checked points. Generators must be bijections of the declared degree, and seed points must be in range. Output uses CSR offsets/indices.","crate::kernels::permgroup::minimal_blocks",DimensionPolicy::Dimensionless,
(n:usize,generators:Vec<Vec<u32>>,a:u32,b:u32)->Vec<Vec<usize>>,|n:usize,generators:Vec<Vec<u32>>,a:u32,b:u32|super::permutation_rows(crate::kernels::permgroup::minimal_blocks(n,&generators,a,b)));
binding!(native_perm_systems,"num_permgroup_minimal_block_systems","numerical_vector","All native minimal block systems under a finite degree/work budget. Top-level offsets group CSR block rows into systems; no rectangular padding changes block membership.","crate::kernels::permgroup::minimal_block_systems",DimensionPolicy::Dimensionless,
(n:usize,generators:Vec<Vec<u32>>)->Result<crate::rules::numerical_abi::codecs::NativeBlockSystems,OpError>,|n:usize,generators:Vec<Vec<u32>>|super::block_systems(crate::kernels::permgroup::minimal_block_systems(n,&generators)));

/// Discrete outputs that require CSR variable-length rows.
pub const NATIVE_DISCRETE_EXPORTS: &[&CodecExport] = &[
    &native_perm_cycles::EXPORT,
    &native_perm_orbits::EXPORT,
    &native_perm_power_exact::EXPORT,
    &native_perm_blocks::EXPORT,
    &native_perm_systems::EXPORT,
];
