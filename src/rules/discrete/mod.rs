//! Discrete mathematics: cryptography, coding theory, finite fields, graphs,
//! algebraic topology, fractals and chaos, computer graphics and finite
//! groups.
//!
//! These peripheral domains are utility wrappers around the graph engine:
//! their inputs and outputs are ordinary terms (numbers, symbols, lists), so
//! they compose with the rest of the system. A graph's adjacency matrix is a
//! `list(list(...))` that `det`, `charpoly` or `eigenvals` of the linear
//! algebra rules accept directly; a transformation matrix of the graphics
//! module multiplies with `matmul`.
//!
//! All operators are *pure functions of their arguments*: nothing is
//! random. Where the legacy code drew a random number (key pairs, ECDSA
//! nonces) the value is now an explicit argument.
//!
//! Every operator is reduced by a kernel that fires when its arguments have
//! the expected literal shape; anything else (a symbolic modulus, a
//! malformed graph, a failing precondition) stays unreduced.
//!
//! # Conventions
//!
//! * Polynomials over a finite field are coefficient lists, **highest degree
//!   first**, without leading zeros (the zero polynomial is `list()`).
//! * "No result" answers are the boolean `false` (no bipartition, no
//!   topological order, negative cycle).
//!
//! # Operators
//!
//! | module | operators |
//! |---|---|
//! | [`crypto`] | `ec_curve`, `ec_on_curve`, `ec_neg`, `ec_double`, `ec_add`, `ec_mul`, `ec_order`, `ec_is_infinity`, `ec_x`, `ec_y`, `ecdh_public`, `ecdh_shared`, `ec_compress`, `ec_decompress`, `ecdsa_sign`, `ecdsa_verify`, `rsa_keygen`, `rsa_encrypt`, `rsa_decrypt` |
//! | [`coding`] | `hamming_distance`, `hamming_weight`, `hamming_encode`, `hamming_check`, `hamming_decode`, `rs_encode`, `rs_check`, `rs_decode`, `rs_error_count`, `bch_encode`, `bch_decode`, `crc32`, `crc32_verify`, `crc32_update`, `crc32_finalize`, `crc16`, `crc8`, `gf256_*`, `gf256_poly_*` |
//! | [`finite_field`] | `gf_*` (prime field), `gfp_*` (polynomials over GF(p)), `gfx_*` (extension fields GF(p)\[x\]/(m)) |
//! | [`gf_factor`] | factorisation over GF(p): `gfp_squarefree`, `gfp_ddf`, `gfp_edf`, `gfp_factor`, `gfp_berlekamp`, `factor_mod`, `gfp_powmod`, `gfp_invmod` |
//! | [`graphs`] | `graph`/`digraph` terms, `graph_*` queries, traversals, shortest paths, spanning trees, flows, matchings, colouring, isomorphism, products |
//! | [`topology`] | `sc`, `sc_complex`, `sc_dimension`, `sc_simplices`, `sc_euler_characteristic`, `sc_boundary`, `sc_boundary_matrix`, `sc_coboundary_matrix`, `sc_chain_boundary`, `sc_betti`, `sc_cohomology_betti`, `sc_verify_boundary`, `sc_verify_coboundary`, `sc_components`, `sc_grid`, `sc_torus`, `vietoris_rips`, `vietoris_rips_filtration`, `betti_at_radius`, `persistence`, `euclidean_distance`, `sc_homology`, `sc_homology_all`, `sc_reduced_homology`, `sc_cohomology`, `sc_homology_mod`, `sc_relative_homology`, `sc_union`, `sc_intersection`, `sc_is_subcomplex`, `sc_skeleton`, `sc_mayer_vietoris_check`, `sc_mayer_vietoris_ranks`, `sc_euler_from_betti`, `sc_verify_euler`, `sc_simplex`, `sc_sphere`, `sc_klein_bottle`, `sc_projective_plane`, `sc_dunce_cap`, `sc_cone`, `sc_suspension`, `sc_join`, `sc_wedge`, `sc_product`, `sc_fundamental_group`, `sc_fundamental_group_raw`, `sc_pi1_abelianization`, `presentation_abelianization`, `sc_is_simplicial_map`, `sc_induced_homology`, `sc_induced_homology_rank`, `persistent_homology`, `rips_filtration`, `rips_persistence`, `rips_persistence_diagram`, `bottleneck_distance`, `cubical_betti`, `cubical_euler_characteristic` |
//! | [`fractal`] | `mandelbrot_escape`, `julia_escape`, `burning_ship_escape`, `multibrot_escape`, `newton_fractal_root`, `mandelbrot_iterate`, `mandelbrot_orbit`, `mandelbrot_fixed_points`, `mandelbrot_stability`, `complex_map_fixed_points`, `complex_map_stability`, `map_fixed_points`, `map_stability`, `lyapunov_exponent`, `logistic_iterate`, `logistic_bifurcation`, `logistic_lyapunov`, `lorenz`, `lorenz_orbit`, `lorenz_lyapunov`, `rossler_orbit`, `henon_orbit`, `tinkerbell_orbit`, `ifs_apply`, `ifs_generate`, `similarity_dimension`, `moran_dimension`, `box_counting`, `correlation_dimension`, `orbit_density`, `orbit_entropy` |
//! | [`graphics`] | `translation_2d/3d`, `scaling_2d/3d`, `shear_2d`, `rotation_2d`, `rotation_3d_x/y/z`, `rotation_axis_angle`, `reflection_2d/3d`, `perspective`, `orthographic`, `look_at`, `apply_transform`, `apply_transform_vector`, `bezier`, `bezier_derivative`, `bezier_split`, `bspline`, `catmull_rom`, `quat_mul`, `quat_conj`, `quat_inverse`, `quat_norm`, `quat_normalize`, `quat_from_axis_angle`, `quat_rotate`, `quat_to_matrix`, `quat_slerp`, `mesh_transform`, `mesh_normals`, `mesh_triangulate`, `ray_sphere`, `ray_plane`, `ray_triangle`, `reflect`, `refract`, `barycentric` |
//! | [`groups`] | `group` term; `cyclic_group`, `dihedral_group`, `symmetric_group`, `klein_four_group`, `group_from_table`; `group_elements`, `group_order`, `group_identity`, `group_mul`, `group_inverse`, `group_is_abelian`, `group_element_order`, `group_conjugacy_classes`, `group_center`, `group_is_valid`, `group_subgroups`, `group_cosets`, `group_is_normal`; `representation_is_valid`, `group_character`; `perm_compose`, `perm_inverse`, `perm_order`, `perm_cycles`, `perm_sign` |
//! | [`perm_groups`] | Schreier–Sims: `perm_group_order`, `perm_group_contains`, `perm_group_base`, `perm_group_strong_generators`, `perm_group_basic_orbits`, `perm_group_orbits`, `perm_group_is_transitive`, `perm_group_stabilizer_order`, `perm_group_elements`, `group_from_perms`, `perm_group_is_abelian`, `perm_group_derived_series`, `perm_group_is_solvable`, `perm_group_lower_central_series`, `perm_group_is_nilpotent`, `perm_group_derived_subgroup`; Todd–Coxeter: `todd_coxeter`, `todd_coxeter_index`, `fp_group_order`, `fp_group` |
//! | [`group_theory`] | `group_generate`, `group_subgroup`, `group_derived_subgroup`, `group_derived_series`, `group_is_solvable`, `group_lower_central_series`, `group_is_nilpotent`, `group_nilpotency_class`, `group_normal_subgroups`, `group_is_simple`, `group_sylow_subgroup`, `group_sylow_count`, `group_sylow_subgroups`, `group_quotient`, `group_direct_product`, `group_isomorphism`, `group_is_isomorphic`, `group_automorphism_count`, `group_inner_automorphism_count`, `group_outer_automorphism_count`, `group_exponent`, `group_element_orders`, `group_order_statistics` |
//! | [`representations`] | Burnside–Dixon character tables: `group_class_count`, `group_class_sizes`, `group_class_representatives`, `group_class_index`, `group_character_table`, `group_character_degrees`, `character_table_is_orthogonal`, `character_inner_product`, `character_is_irreducible`, `character_decompose`, `character_tensor`, `character_sym_square`, `character_alt_square`, `character_adams`, `character_conjugate`, `character_regular`, `character_of_matrices`, `representation_decompose`, `character_projection` |
//! | [`point_groups`] | `point_group`, `point_group_order`, `point_group_is_crystallographic`, `point_group_classes`, `point_group_class_sizes`, `point_group_irreps`, `point_group_irrep_dimensions`, `point_group_character_table`, `point_group_decompose`, `point_group_multiplicities`, `point_group_vector_character`, `point_group_rotation_character`, `point_group_ir_active`, `point_group_raman_active`, `point_group_function_irreps`, `point_group_hm`, `point_group_from_hm`, `point_group_crystal_system`, `crystallographic_point_groups`, `crystal_systems`, `bravais_lattices`, `crystallographic_restriction`, `crystallographic_min_dimension`, `molecule_symmetry_operations`, `molecule_point_group`, `molecule_decomposition`, `molecule_vibrations`, `molecule_vibrations_table`, `molecule_spectroscopy` |
//!
//! See the documentation of each module for the exact term formats.

pub mod coding;
pub mod crypto;
pub mod finite_field;
pub mod fractal;
pub mod gf_factor;
pub mod graphics;
pub mod graphs;
#[cfg(test)]
mod graphs_tests;
pub mod group_theory;
pub mod groups;
pub mod perm_groups;
pub mod point_groups;
pub mod representations;
pub mod topology;

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Signed;
use num_traits::ToPrimitive;
use num_traits::Zero;

use crate::graph::Arity;
use crate::graph::Cx;
use crate::graph::Graph;
use crate::graph::Kernel;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::OpDescriptor;
use crate::graph::OpFlags;
use crate::graph::OpId;
use crate::graph::Outcome;
use crate::graph::Payload;
use crate::graph::RuleError;
use crate::graph::RuleSet;
use crate::graph::Tier;
use crate::graph::dim::DimCtx;
use crate::graph::dim::DimError;
use crate::graph::dim::DimRule;
use crate::graph::dim::DimTerm;
use crate::graph::op::core;
use crate::graph::rule::Installer;

use super::complex::complex;
use super::elementary::elementary;
use super::linalg::linalg;
use super::number_theory::number_theory;

/// The discrete mathematics rule set.
#[must_use]
pub fn discrete() -> RuleSet {
    RuleSet::new("discrete", install)
        .needs(elementary())
        .needs(number_theory())
        .needs(linalg())
        .needs(complex())
}

/// A count, a field element or a coordinate ratio: dimensionless whatever
/// the dimensions of the objects it queries.
#[allow(clippy::missing_const_for_fn, clippy::unnecessary_wraps)]
fn dimensionless_result(
    _ctx: &DimCtx<'_>,
    _args: &[DimTerm],
) -> Result<DimTerm, DimError> {
    Ok(DimTerm::dimensionless())
}

/// The dimension signature of a discrete operator.
///
/// Discrete mathematics mostly produces integers, truth values, labels and
/// opaque structures: counts and finite-field elements are dimensionless,
/// homogeneous integer lists are `List`, and graphs, groups, simplicial
/// complexes, permutations and predicates are `NonMetric` containers. The
/// graphics/fractal operators that act on geometric points keep the point
/// dimension (`Same`) or multiply dimensions (`Product`).
#[allow(clippy::too_many_lines)]
fn dim_rule(name: &str) -> DimRule {
    use DimRule::Custom;
    use DimRule::Inverse;
    use DimRule::List;
    use DimRule::NonMetric;
    use DimRule::Product;
    use DimRule::Same;
    // Counts, invariants and finite-field elements.
    if matches!(
        name,
        "graph_node_count"
            | "graph_chromatic_number"
            | "algebraic_connectivity"
            | "ec_x"
            | "ec_y"
            | "ec_order"
            | "rsa_encrypt"
            | "rsa_decrypt"
            | "gf_add"
            | "gf_sub"
            | "gf_mul"
            | "gf_div"
            | "gf_neg"
            | "gf_inv"
            | "gf_pow"
            | "gfp_norm"
            | "gfp_degree"
            | "gfp_eval"
            | "gf256_add"
            | "gf256_mul"
            | "gf256_div"
            | "gf256_inv"
            | "gf256_pow"
            | "gf256_exp"
            | "gf256_log"
            | "gf256_poly_eval"
            | "hamming_distance"
            | "hamming_weight"
            | "rs_error_count"
            | "code_min_distance"
            | "code_rate"
            | "crc32"
            | "crc32_update"
            | "crc32_finalize"
            | "crc16"
            | "crc8"
            | "group_order"
            | "group_element_order"
            | "group_sylow_count"
            | "group_class_count"
            | "group_class_index"
            | "group_automorphism_count"
            | "group_inner_automorphism_count"
            | "group_outer_automorphism_count"
            | "group_exponent"
            | "group_nilpotency_class"
            | "perm_order"
            | "perm_sign"
            | "perm_group_order"
            | "perm_group_stabilizer_order"
            | "todd_coxeter_index"
            | "fp_group_order"
            | "character_inner_product"
            | "point_group_order"
            | "crystallographic_min_dimension"
            | "sc_dimension"
            | "sc_simplex_dimension"
            | "sc_euler_characteristic"
            | "sc_betti"
            | "sc_cohomology_betti"
            | "sc_components"
            | "sc_euler_from_betti"
            | "sc_induced_homology_rank"
            | "betti_at_radius"
            | "bottleneck_distance"
            | "cubical_euler_characteristic"
            | "mandelbrot_escape"
            | "julia_escape"
            | "burning_ship_escape"
            | "multibrot_escape"
            | "mandelbrot_stability"
            | "complex_map_stability"
            | "map_stability"
            | "lyapunov_exponent"
            | "logistic_lyapunov"
            | "lorenz_lyapunov"
            | "similarity_dimension"
            | "moran_dimension"
            | "box_counting"
            | "correlation_dimension"
            | "orbit_density"
            | "orbit_entropy"
    ) {
        return Custom(dimensionless_result);
    }
    // Homogeneous integer lists, tables and matrices.
    if matches!(
        name,
        "graph_adjacency"
            | "graph_incidence"
            | "graph_laplacian"
            | "spectral_analysis"
            | "graph_greedy_coloring"
            | "gfp_add"
            | "gfp_sub"
            | "gfp_mul"
            | "gfp_divmod"
            | "gfp_gcd"
            | "gfp_egcd"
            | "gfp_monic"
            | "gfp_derivative"
            | "gfp_powmod"
            | "gfp_invmod"
            | "gfp_squarefree"
            | "gfp_ddf"
            | "gfp_edf"
            | "gfp_factor"
            | "gfp_berlekamp"
            | "gfx_reduce"
            | "gfx_add"
            | "gfx_sub"
            | "gfx_mul"
            | "gfx_div"
            | "gfx_neg"
            | "gfx_inv"
            | "gfx_pow"
            | "hamming_encode"
            | "hamming_check"
            | "hamming_decode"
            | "rs_encode"
            | "rs_check"
            | "rs_decode"
            | "bch_encode"
            | "bch_decode"
            | "interleave"
            | "deinterleave"
            | "conv_encode"
            | "gf256_poly_add"
            | "gf256_poly_mul"
            | "gf256_poly_scale"
            | "gf256_poly_derivative"
            | "gf256_poly_mod"
            | "gf256_poly_divmod"
            | "gf256_poly_gcd"
            | "group_character"
            | "group_character_table"
            | "group_character_degrees"
            | "group_element_orders"
            | "group_order_statistics"
            | "group_class_sizes"
            | "perm_compose"
            | "perm_inverse"
            | "perm_group_base"
            | "perm_group_strong_generators"
            | "perm_group_basic_orbits"
            | "perm_group_orbits"
            | "character_tensor"
            | "character_sym_square"
            | "character_alt_square"
            | "character_adams"
            | "character_conjugate"
            | "character_regular"
            | "character_of_matrices"
            | "character_projection"
            | "point_group_class_sizes"
            | "point_group_irrep_dimensions"
            | "point_group_character_table"
            | "point_group_multiplicities"
            | "point_group_vector_character"
            | "point_group_rotation_character"
            | "molecule_symmetry_operations"
            | "molecule_vibrations"
            | "molecule_vibrations_table"
            | "molecule_spectroscopy"
            | "sc_simplices"
            | "sc_boundary"
            | "sc_boundary_matrix"
            | "sc_coboundary_matrix"
            | "sc_chain_boundary"
            | "sc_homology"
            | "sc_homology_all"
            | "sc_reduced_homology"
            | "sc_cohomology"
            | "sc_homology_mod"
            | "sc_relative_homology"
            | "sc_mayer_vietoris_ranks"
            | "sc_induced_homology"
            | "persistence"
            | "persistent_homology"
            | "rips_persistence"
            | "rips_persistence_diagram"
            | "vietoris_rips_filtration"
            | "cubical_betti"
            | "mandelbrot_iterate"
            | "mandelbrot_orbit"
            | "mandelbrot_fixed_points"
            | "complex_map_fixed_points"
            | "map_fixed_points"
            | "logistic_iterate"
            | "logistic_bifurcation"
            | "lorenz"
            | "lorenz_orbit"
            | "rossler_orbit"
            | "henon_orbit"
            | "tinkerbell_orbit"
            | "ifs_generate"
            | "translation_2d"
            | "translation_3d"
            | "scaling_2d"
            | "scaling_3d"
            | "shear_2d"
            | "rotation_2d"
            | "rotation_3d_x"
            | "rotation_3d_y"
            | "rotation_3d_z"
            | "rotation_axis_angle"
            | "reflection_2d"
            | "reflection_3d"
            | "perspective"
            | "orthographic"
            | "look_at"
            | "quat_from_axis_angle"
            | "quat_to_matrix"
            | "mesh_triangulate"
    ) {
        return List;
    }
    // Geometric values that keep or combine point dimensions.
    match name {
        | "euclidean_distance"
        | "newton_fractal_root"
        | "ifs_apply"
        | "apply_transform"
        | "apply_transform_vector"
        | "bezier"
        | "bezier_derivative"
        | "bezier_split"
        | "bspline"
        | "catmull_rom"
        | "quat_conj"
        | "quat_norm"
        | "quat_normalize"
        | "quat_slerp"
        | "mesh_transform"
        | "mesh_normals"
        | "reflect"
        | "refract" => return Same(0),
        | "quat_rotate" => return Same(1),
        | "quat_mul" => return Product,
        | "quat_inverse" => return Inverse,
        | "barycentric" => return Custom(dimensionless_result),
        | _ => {},
    }
    // Graphs, groups, complexes, permutations and predicates are opaque
    // containers of labels, not scalars.
    if name == "graph"
        || name == "digraph"
        || name == "sc"
        || name == "group"
        || name.starts_with("graph_")
        || name.starts_with("sc_")
        || name.starts_with("group_")
        || name.starts_with("perm_group_")
        || name.starts_with("point_group_")
        || name.starts_with("molecule_")
        || name.starts_with("character_")
        || name.starts_with("ec_")
        || name.starts_with("lie_")
        || matches!(
            name,
            "ecdh_public"
                | "ecdh_shared"
                | "ec_compress"
                | "ec_decompress"
                | "ecdsa_sign"
                | "ecdsa_verify"
                | "rsa_keygen"
                | "gate"
        )
    {
        return NonMetric;
    }
    // The remaining discrete operators are predicates, field objects and
    // symbolic factorisations: not scalars.
    NonMetric
}

fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    crypto::install(i)?;
    finite_field::install(i)?;
    gf_factor::install(i)?;
    coding::install(i)?;
    graphs::install(i)?;
    topology::install(i)?;
    fractal::install(i)?;
    graphics::install(i)?;
    groups::install(i)?;
    perm_groups::install(i)?;
    group_theory::install(i)?;
    representations::install(i)?;
    point_groups::install(i)?;
    Ok(())
}

// ----------------------------------------------------------------------
// Values and kernels
// ----------------------------------------------------------------------

/// What a discrete kernel computed.
pub(crate) enum V {
    /// An already built term.
    Node(NodeId),
    /// An exact integer.
    Int(BigInt),
    /// A truth value.
    Bool(bool),
    /// A float.
    Float(f64),
    /// An exact rational.
    Rat(BigRational),
    /// A string literal.
    Str(String),
    /// `list(...)`.
    List(Vec<Self>),
}

impl V {
    pub(crate) fn int(v: i64) -> Self {
        Self::Int(BigInt::from(v))
    }

    pub(crate) fn uint(v: usize) -> Self {
        Self::Int(BigInt::from(v))
    }

    pub(crate) fn ints<T: Into<BigInt>>(items: impl IntoIterator<Item = T>) -> Self {
        Self::List(items.into_iter().map(|v| Self::Int(v.into())).collect())
    }

    pub(crate) fn nodes(items: &[NodeId]) -> Self {
        Self::List(items.iter().map(|&n| Self::Node(n)).collect())
    }

    pub(crate) fn build(
        self,
        graph: &mut Graph,
    ) -> NodeId {
        match self {
            | Self::Node(n) => n,
            | Self::Int(v) => graph.num(Number::Int(v)),
            | Self::Bool(b) => graph.lit(Payload::Bool(b)),
            | Self::Float(x) => graph.float(x),
            | Self::Rat(r) => graph.num(Number::rat(r)),
            | Self::Str(s) => graph.lit(Payload::Str(s.into())),
            | Self::List(items) => {
                let nodes: Vec<NodeId> = items.into_iter().map(|v| v.build(graph)).collect();
                graph.node(core::LIST, &nodes)
            },
        }
    }
}

/// A computation on the argument nodes of an operator.
pub(crate) type Run = fn(&mut Cx<'_>, &[NodeId]) -> Option<V>;

/// Reduces `op` with a [`Run`] function.
struct Fun {
    op: OpId,
    run: Run,
    /// Lists are larger than the request; pin them so that extraction
    /// returns them.
    pin: bool,
}

impl Kernel for Fun {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let args = cx.graph.children(node).to_vec();
        let Some(value) = (self.run)(cx, &args) else {
            return Outcome::Pass;
        };
        let result = value.build(cx.graph);
        let head = cx.graph.op(result);
        let structured = head == core::LIST
            || matches!(
                &*cx.graph.ops().get(head).name,
                "graph" | "digraph" | "group" | "sc"
            );
        if self.pin && structured {
            Outcome::Pinned(result)
        } else {
            Outcome::Equal(result)
        }
    }

    fn revisit(&self) -> bool {
        true
    }
}

/// Registers an operator with a kernel; list results are pinned.
pub(crate) fn def(
    i: &mut Installer<'_>,
    name: &str,
    arity: Arity,
    run: Run,
) -> Result<OpId, RuleError> {
    let op = i.op(OpDescriptor::new(name, arity))?;
    i.dim(op, dim_rule(name));
    i.kernel(
        &format!("discrete/{name}"),
        Tier::Normalize,
        Fun { op, run, pin: true },
    );
    Ok(op)
}

/// Registers a *request* operator: it stands for a computation whose answer
/// is a term still to be simplified (a formula with `sin`, `ln`, ...), so
/// it is flagged heavy and the answer is not pinned.
pub(crate) fn def_request(
    i: &mut Installer<'_>,
    name: &str,
    arity: Arity,
    run: Run,
) -> Result<OpId, RuleError> {
    let op = i.op(OpDescriptor::new(name, arity)
        .flags(OpFlags::HEAVY)
        .cost(100))?;
    i.dim(op, dim_rule(name));
    i.kernel(
        &format!("discrete/{name}"),
        Tier::Reduce,
        Fun { op, run, pin: false },
    );
    Ok(op)
}

/// Registers a data constructor with no kernel (an inert term).
pub(crate) fn def_inert(
    i: &mut Installer<'_>,
    name: &str,
    arity: Arity,
) -> Result<OpId, RuleError> {
    let op = i.op(OpDescriptor::new(name, arity))?;
    i.dim(op, dim_rule(name));
    Ok(op)
}

// ----------------------------------------------------------------------
// Reading arguments
// ----------------------------------------------------------------------

/// A literal integer.
pub(crate) fn big(
    g: &Graph,
    n: NodeId,
) -> Option<BigInt> {
    match g.number_of(n)? {
        | Number::Int(v) => Some(v.clone()),
        | _ => None,
    }
}

/// A literal integer that fits `i64`.
pub(crate) fn small(
    g: &Graph,
    n: NodeId,
) -> Option<i64> {
    big(g, n)?.to_i64()
}

/// A literal non-negative integer.
pub(crate) fn idx(
    g: &Graph,
    n: NodeId,
) -> Option<usize> {
    big(g, n)?.to_usize()
}

/// A literal number as a float.
pub(crate) fn float(
    g: &Graph,
    n: NodeId,
) -> Option<f64> {
    g.number_of(n).map(Number::to_f64)
}

/// The children of a `list`.
pub(crate) fn items(
    g: &Graph,
    n: NodeId,
) -> Option<Vec<NodeId>> {
    if g.op(n) == core::LIST {
        return Some(g.children(n).to_vec());
    }
    // The argument may be a request that has since been reduced to a list.
    g.enodes(g.find(n))
        .find(|&e| g.op(e) == core::LIST)
        .map(|e| g.children(e).to_vec())
}

/// A `list` of literal integers.
pub(crate) fn bigs(
    g: &Graph,
    n: NodeId,
) -> Option<Vec<BigInt>> {
    items(g, n)?.into_iter().map(|c| big(g, c)).collect()
}

/// A `list` of literal non-negative integers.
pub(crate) fn idxs(
    g: &Graph,
    n: NodeId,
) -> Option<Vec<usize>> {
    items(g, n)?.into_iter().map(|c| idx(g, c)).collect()
}

/// A `list` of bytes (integers in `0..=255`).
pub(crate) fn bytes(
    g: &Graph,
    n: NodeId,
) -> Option<Vec<u8>> {
    items(g, n)?
        .into_iter()
        .map(|c| small(g, c).and_then(|v| u8::try_from(v).ok()))
        .collect()
}

/// A `list` of floats.
pub(crate) fn floats(
    g: &Graph,
    n: NodeId,
) -> Option<Vec<f64>> {
    items(g, n)?.into_iter().map(|c| float(g, c)).collect()
}

/// A list of lists.
pub(crate) fn rows(
    g: &Graph,
    n: NodeId,
) -> Option<Vec<Vec<NodeId>>> {
    items(g, n)?.into_iter().map(|r| items(g, r)).collect()
}

/// Bytes as a list value.
pub(crate) fn byte_list(data: &[u8]) -> V {
    V::ints(data.iter().map(|&b| i64::from(b)))
}

// ----------------------------------------------------------------------
// Building terms
// ----------------------------------------------------------------------

/// Whether `n` is a literal zero.
pub(crate) fn is_zero(
    g: &Graph,
    n: NodeId,
) -> bool {
    g.number_of(n).is_some_and(Number::is_zero)
}

/// The sum of `terms`; literal numbers are folded.
pub(crate) fn sum(
    g: &mut Graph,
    terms: &[NodeId],
) -> NodeId {
    let mut acc: Option<Number> = None;
    let mut rest = Vec::new();
    for &t in terms {
        match g.number_of(t).cloned() {
            | Some(v) => {
                acc = Some(match acc {
                    | Some(a) => a.add(&v),
                    | None => v,
                });
            },
            | None => rest.push(t),
        }
    }
    let mut all = Vec::new();
    if let Some(a) = acc
        && (!a.is_zero() || rest.is_empty())
    {
        all.push(g.num(a));
    }
    all.extend(rest);
    match all.as_slice() {
        | [] => g.int(0),
        | [only] => *only,
        | _ => g.node(core::ADD, &all),
    }
}

/// The product of `factors`; literal numbers are folded.
pub(crate) fn prod(
    g: &mut Graph,
    factors: &[NodeId],
) -> NodeId {
    let mut acc = Number::from(1_i64);
    let mut rest = Vec::new();
    for &f in factors {
        match g.number_of(f).cloned() {
            | Some(v) => acc = acc.mul(&v),
            | None => rest.push(f),
        }
    }
    if acc.is_zero() {
        return g.int(0);
    }
    let mut all = Vec::new();
    if !acc.is_one() || rest.is_empty() {
        all.push(g.num(acc));
    }
    all.extend(rest);
    match all.as_slice() {
        | [only] => *only,
        | _ => g.node(core::MUL, &all),
    }
}

/// `-n`.
pub(crate) fn neg(
    g: &mut Graph,
    n: NodeId,
) -> NodeId {
    let minus_one = g.int(-1);
    prod(g, &[minus_one, n])
}

/// `n^e`.
pub(crate) fn pow(
    g: &mut Graph,
    n: NodeId,
    e: NodeId,
) -> NodeId {
    g.node(core::POW, &[n, e])
}

/// `1 / n`.
pub(crate) fn inv(
    g: &mut Graph,
    n: NodeId,
) -> NodeId {
    let minus_one = g.int(-1);
    pow(g, n, minus_one)
}

/// `a / b`.
pub(crate) fn div(
    g: &mut Graph,
    a: NodeId,
    b: NodeId,
) -> NodeId {
    let b = inv(g, b);
    prod(g, &[a, b])
}

/// `name(args...)` for an operator registered by another rule set.
pub(crate) fn apply(
    g: &mut Graph,
    name: &str,
    args: &[NodeId],
) -> Option<NodeId> {
    let op = g.ops().lookup(name)?;
    g.try_node(op, args)
}

/// `list(items...)`.
pub(crate) fn list(
    g: &mut Graph,
    items: &[NodeId],
) -> NodeId {
    g.node(core::LIST, items)
}

/// A matrix term from rows of nodes.
pub(crate) fn matrix(
    g: &mut Graph,
    rows: &[Vec<NodeId>],
) -> NodeId {
    let rows: Vec<NodeId> = rows.iter().map(|r| list(g, r)).collect();
    list(g, &rows)
}

/// Simplifies every node with the running program.
pub(crate) fn simplified(
    cx: &mut Cx<'_>,
    nodes: &[NodeId],
) -> Vec<NodeId> {
    nodes.iter().map(|&n| cx.simplify(n)).collect()
}

/// The product of two matrices of nodes; entries are simplified.
pub(crate) fn matmul(
    cx: &mut Cx<'_>,
    a: &[Vec<NodeId>],
    b: &[Vec<NodeId>],
) -> Option<Vec<Vec<NodeId>>> {
    let inner = b.len();
    let cols = b.first().map_or(0, Vec::len);
    let mut out = Vec::with_capacity(a.len());
    for row in a {
        if row.len() != inner || b.iter().any(|r| r.len() != cols) {
            return None;
        }
        let mut new_row = Vec::with_capacity(cols);
        #[allow(clippy::needless_range_loop)]
        // index is used for more than one array / arithmetic; iterator form would not be clearer
        for j in 0..cols {
            let terms: Vec<NodeId> = (0..inner)
                .map(|k| prod(cx.graph, &[row[k], b[k][j]]))
                .collect();
            let s = sum(cx.graph, &terms);
            new_row.push(cx.simplify(s));
        }
        out.push(new_row);
    }
    Some(out)
}

// ----------------------------------------------------------------------
// Integer helpers shared by several modules
// ----------------------------------------------------------------------

/// Remainder in `[0, |m|)`.
pub(crate) fn modulo(
    a: &BigInt,
    m: &BigInt,
) -> BigInt {
    let m = m.abs();
    let r = a % &m;
    if r.is_negative() { r + m } else { r }
}

/// The inverse of `a` modulo `m`, if `gcd(a, m) = 1`.
pub(crate) fn mod_inverse(
    a: &BigInt,
    m: &BigInt,
) -> Option<BigInt> {
    if !m.is_positive() {
        return None;
    }
    let (mut r0, mut r1) = (modulo(a, m), m.clone());
    let (mut s0, mut s1) = (BigInt::one(), BigInt::zero());
    while !r1.is_zero() {
        let q = &r0 / &r1;
        let r2 = &r0 - &q * &r1;
        let s2 = &s0 - &q * &s1;
        (r0, r1) = (r1, r2);
        (s0, s1) = (s1, s2);
    }
    if r0.is_one() {
        Some(modulo(&s0, m))
    } else if m.is_one() {
        Some(BigInt::zero())
    } else {
        None
    }
}

/// An exact rational from a literal number; floats are converted exactly.
pub(crate) fn rational(n: &Number) -> Option<BigRational> {
    match n {
        | Number::Float(x) => BigRational::from_float(*x),
        | other => other.to_rational(),
    }
}

/// Helpers for the tests of the submodules.
#[cfg(test)]
pub(crate) mod test_util {
    use super::discrete;
    use crate::rules::testing::simplify;

    /// The closed form of `src`.
    pub fn s(src: &str) -> String {
        simplify(&[discrete()], src)
    }

    /// The numeric rows of the closed form of `src`, a list of lists.
    pub fn num_rows(src: &str) -> Vec<Vec<num_complex::Complex64>> {
        use std::collections::HashMap;

        use crate::graph::Budget;
        use crate::graph::ClosedForm;
        use crate::graph::Engine;
        use crate::graph::Env;
        use crate::graph::Extractor;
        use crate::graph::Graph;
        use crate::graph::Saturate;
        let mut g = Graph::new();
        let engine = Engine::install(&mut g, &[discrete()]).unwrap_or_else(|e| panic!("{e}"));
        let root = g
            .parse(src)
            .unwrap_or_else(|e| panic!("cannot parse `{src}`: {e}"));
        engine.run(
            &mut g,
            &[root],
            &Env::symbolic(),
            &Saturate,
            &Budget::default(),
        );
        let node = Extractor::new(&g, &[root], &ClosedForm)
            .build(&mut g, root)
            .unwrap_or_else(|| panic!("`{src}` not reduced"));
        let bindings = HashMap::new();
        let rows = super::items(&g, node).unwrap_or_else(|| panic!("`{src}` is not a list"));
        rows.into_iter()
            .map(|r| {
                super::items(&g, r)
                    .unwrap_or_default()
                    .into_iter()
                    .map(|x| {
                        g.eval_complex(x, &bindings)
                            .unwrap_or_else(|| panic!("cannot evaluate {}", g.display(x)))
                    })
                    .collect()
            })
            .collect()
    }

    /// A small deterministic generator.
    pub struct Lcg(pub u64);

    impl Lcg {
        pub(crate) fn next(
            &mut self,
            bound: u64,
        ) -> u64 {
            self.0 = self
                .0
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (self.0 >> 33) % bound
        }
    }

    /// Every integer in `text`, in order.
    pub fn nums(text: &str) -> Vec<i64> {
        text.split(|c: char| !(c.is_ascii_digit() || c == '-'))
            .filter(|t| !t.is_empty() && *t != "-")
            .filter_map(|t| t.parse().ok())
            .collect()
    }

    /// The rows of a `list(list(..), list(..))`.
    pub fn nested(text: &str) -> Vec<Vec<i64>> {
        let inner = text
            .strip_prefix("list(")
            .and_then(|t| t.strip_suffix(')'))
            .unwrap_or("");
        inner.split("list(").skip(1).map(nums).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rules::testing::reduce_with;
    use crate::rules::testing::simplify;

    #[test]
    fn discrete_installs_with_the_standard_sets() {
        let s = simplify(
            &[discrete()],
            "ec_add(list(2, 3, 97), list(3, 6), list(3, 6))",
        );
        assert_eq!(s, "list(80, 10)");
        let (text, reduced) = reduce_with(&[discrete()], "hamming_weight(list(1, 0, 1, 1))", &[]);
        assert!(reduced);
        assert_eq!(text, "3");
    }
}
