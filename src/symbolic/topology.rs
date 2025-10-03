//! # Computational Algebraic Topology
//!
//! This module provides data structures and algorithms for computational algebraic
//! topology. It includes implementations for simplices, simplicial complexes, chain
//! complexes, and the computation of homology groups and Betti numbers.

use crate::numerical::sparse::rank;
use sprs::{CsMat, TriMat};
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};

// =====================================================================================
// region: Core Data Structures
// =====================================================================================

/// Represents a k-simplex as a set of its vertex indices.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Simplex(pub BTreeSet<usize>);

impl Simplex {
    /// Creates a new `Simplex` from a slice of vertex indices.
    ///
    /// The vertex indices are stored in a `BTreeSet` to ensure uniqueness and ordering.
    ///
    /// # Arguments
    /// * `vertices` - A slice of `usize` representing the vertex indices of the simplex.
    ///
    /// # Returns
    /// A new `Simplex` instance.
    pub fn new(vertices: &[usize]) -> Self {
        Simplex(vertices.iter().cloned().collect())
    }

    /// Returns the dimension of the simplex.
    ///
    /// The dimension of a k-simplex is `k`. For example, a 0-simplex (point) has dimension 0,
    /// a 1-simplex (edge) has dimension 1, a 2-simplex (triangle) has dimension 2, etc.
    /// It is calculated as `number_of_vertices - 1`.
    ///
    /// # Returns
    /// The dimension as a `usize`.
    pub fn dimension(&self) -> usize {
        self.0.len().saturating_sub(1)
    }

    /// Computes the boundary of the simplex.
    ///
    /// The boundary of a k-simplex is a (k-1)-chain, which is a formal sum of its (k-1)-dimensional faces.
    /// The coefficients are `+1` or `-1` depending on the orientation.
    /// For a 0-simplex (point), the boundary is empty.
    ///
    /// # Returns
    /// A tuple `(faces, coeffs)` where `faces` is a `Vec<Simplex>` of the boundary faces
    /// and `coeffs` is a `Vec<f64>` of their corresponding coefficients.
    pub fn boundary(&self) -> (Vec<Simplex>, Vec<f64>) {
        let mut faces = Vec::new();
        let mut coeffs = Vec::new();
        let vertices: Vec<_> = self.0.iter().cloned().collect();
        if self.dimension() == 0 {
            return (faces, coeffs);
        }

        for i in 0..vertices.len() {
            let face_vertices: BTreeSet<_> = vertices
                .iter()
                .enumerate()
                .filter(|(j, _)| *j != i)
                .map(|(_, &v)| v)
                .collect();
            faces.push(Simplex(face_vertices));
            coeffs.push(if i % 2 == 0 { 1.0 } else { -1.0 });
        }
        (faces, coeffs)
    }
}

/// Represents a k-chain as a formal linear combination of k-simplices.
#[derive(Debug, Clone)]
pub struct Chain {
    pub terms: HashMap<Simplex, f64>,
    pub dimension: usize,
}

/// A k-cochain is an element of the dual space to the k-chains. It can be represented similarly.
pub type Cochain = Chain;

impl Chain {
    /// Creates a new, empty k-chain of a specified dimension.
    ///
    /// # Arguments
    /// * `dimension` - The dimension `k` of the simplices in the chain.
    ///
    /// # Returns
    /// A new `Chain` instance.
    pub fn new(dimension: usize) -> Self {
        Chain {
            terms: HashMap::new(),
            dimension,
        }
    }

    /// Adds a simplex with a given coefficient to the chain.
    ///
    /// If the simplex already exists in the chain, its coefficient is updated.
    ///
    /// # Arguments
    /// * `simplex` - The `Simplex` to add.
    /// * `coeff` - The coefficient for the simplex.
    ///
    /// # Panics
    /// Panics if the dimension of the simplex does not match the chain's dimension.
    pub fn add_term(&mut self, simplex: Simplex, coeff: f64) {
        if simplex.dimension() != self.dimension {
            panic!("Cannot add simplex of wrong dimension to chain.");
        }
        *self.terms.entry(simplex).or_insert(0.0) += coeff;
    }
}

/// Represents a simplicial complex.
#[derive(Debug, Clone, Default)]
pub struct SimplicialComplex {
    simplices: HashSet<Simplex>,
    simplices_by_dim: BTreeMap<usize, Vec<Simplex>>,
}

impl SimplicialComplex {
    /// Creates a new, empty `SimplicialComplex`.
    pub fn new() -> Self {
        Default::default()
    }

    /// Adds a simplex and all its faces (sub-simplices) to the complex.
    ///
    /// This ensures that the complex remains valid (i.e., every face of a simplex
    /// in the complex is also in the complex).
    ///
    /// # Arguments
    /// * `vertices` - A slice of `usize` representing the vertex indices of the simplex to add.
    pub fn add_simplex(&mut self, vertices: &[usize]) {
        let simplex = Simplex::new(vertices);
        pub(crate) fn add_faces(complex: &mut SimplicialComplex, s: Simplex) {
            if complex.simplices.insert(s.clone()) {
                let dim = s.dimension();
                complex
                    .simplices_by_dim
                    .entry(dim)
                    .or_default()
                    .push(s.clone());
                if dim > 0 {
                    let (boundary_faces, _) = s.boundary();
                    for face in boundary_faces {
                        add_faces(complex, face);
                    }
                }
            }
        }
        add_faces(self, simplex);
    }

    /// Returns the dimension of the simplicial complex.
    ///
    /// The dimension of a complex is the highest dimension of any simplex in the complex.
    ///
    /// # Returns
    /// An `Option<usize>` containing the dimension, or `None` if the complex is empty.
    pub fn dimension(&self) -> Option<usize> {
        self.simplices_by_dim.keys().max().cloned()
    }

    /// Returns a reference to the vector of simplices of a specific dimension.
    ///
    /// # Arguments
    /// * `dim` - The dimension of the simplices to retrieve.
    ///
    /// # Returns
    /// An `Option<&Vec<Simplex>>` containing the simplices, or `None` if no simplices of that dimension exist.
    pub fn get_simplices_by_dim(&self, dim: usize) -> Option<&Vec<Simplex>> {
        self.simplices_by_dim.get(&dim)
    }

    /// Constructs the k-th boundary matrix `∂_k` for the simplicial complex.
    ///
    /// The boundary matrix maps k-chains to (k-1)-chains. Its entries are `+1`, `-1`, or `0`,
    /// indicating the incidence and orientation of (k-1)-faces within k-simplices.
    ///
    /// # Arguments
    /// * `k` - The dimension `k` for which to construct the boundary matrix.
    ///
    /// # Returns
    /// An `Option<CsMat<f64>>` representing the sparse boundary matrix, or `None` if `k` is 0
    /// or if there are no simplices of dimension `k` or `k-1`.
    pub fn get_boundary_matrix(&self, k: usize) -> Option<CsMat<f64>> {
        if k == 0 {
            return None;
        }
        let k_simplices = self.get_simplices_by_dim(k)?;
        let k_minus_1_simplices = self.get_simplices_by_dim(k - 1)?;
        let k_minus_1_map: HashMap<&Simplex, usize> = k_minus_1_simplices
            .iter()
            .enumerate()
            .map(|(i, s)| (s, i))
            .collect();

        let mut triplets = Vec::new();
        for (j, simplex_k) in k_simplices.iter().enumerate() {
            let (boundary_faces, coeffs) = simplex_k.boundary();
            for (i, face) in boundary_faces.iter().enumerate() {
                if let Some(&row_idx) = k_minus_1_map.get(face) {
                    triplets.push((row_idx, j, coeffs[i]));
                }
            }
        }
        Some(csr_from_triplets(
            k_minus_1_simplices.len(),
            k_simplices.len(),
            &triplets,
        ))
    }

    /// Applies the k-th boundary operator `∂_k` to a k-chain.
    ///
    /// This effectively computes the boundary of the input chain.
    ///
    /// # Arguments
    /// * `chain` - The `Chain` to apply the boundary operator to.
    ///
    /// # Returns
    /// An `Option<Chain>` representing the resulting (k-1)-chain, or `None` if the boundary
    /// matrix for the given dimension `k` cannot be constructed.
    pub fn apply_boundary_operator(&self, chain: &Chain) -> Option<Chain> {
        let k = chain.dimension;
        if k == 0 {
            return Some(Chain::new(0));
        }

        let boundary_matrix = self.get_boundary_matrix(k)?;
        let k_simplices = self.get_simplices_by_dim(k)?;
        let k_minus_1_simplices = self.get_simplices_by_dim(k - 1)?;

        let mut input_vec = vec![0.0; k_simplices.len()];
        for (i, simplex) in k_simplices.iter().enumerate() {
            if let Some(&coeff) = chain.terms.get(simplex) {
                input_vec[i] = coeff;
            }
        }

        let output_vec = crate::numerical::sparse::sp_mat_vec_mul(&boundary_matrix, &input_vec);

        let mut result_chain = Chain::new(k - 1);
        for (i, &coeff) in output_vec.iter().enumerate() {
            if coeff.abs() > 1e-9 {
                result_chain.add_term(k_minus_1_simplices[i].clone(), coeff);
            }
        }
        Some(result_chain)
    }

    /// Computes the Euler characteristic `χ` of the simplicial complex.
    ///
    /// The Euler characteristic is an important topological invariant, defined as
    /// `χ = Σ_k (-1)^k * (number of k-simplices)`.
    ///
    /// # Returns
    /// An `isize` representing the Euler characteristic.
    pub fn compute_euler_characteristic(&self) -> isize {
        let mut ch = 0;
        for (dim, simplices) in &self.simplices_by_dim {
            let term = simplices.len() as isize;
            if dim % 2 == 0 {
                ch += term;
            } else {
                ch -= term;
            }
        }
        ch
    }
}

/// Represents the full chain complex and its dual, the cochain complex.
pub struct ChainComplex {
    pub complex: SimplicialComplex,
    pub boundary_operators: BTreeMap<usize, CsMat<f64>>,
    pub coboundary_operators: BTreeMap<usize, CsMat<f64>>,
}

impl ChainComplex {
    /// Creates a new `ChainComplex` from a `SimplicialComplex`.
    ///
    /// This constructor automatically computes and stores all relevant boundary
    /// and coboundary operators for the given simplicial complex.
    ///
    /// # Arguments
    /// * `complex` - The `SimplicialComplex` to build the chain complex from.
    ///
    /// # Returns
    /// A new `ChainComplex` instance.
    pub fn new(complex: SimplicialComplex) -> Self {
        let mut boundary_operators = BTreeMap::new();
        let mut coboundary_operators = BTreeMap::new();
        if let Some(max_dim) = complex.dimension() {
            for k in 1..=max_dim {
                if let Some(matrix) = complex.get_boundary_matrix(k) {
                    coboundary_operators.insert(k - 1, matrix.transpose_view().to_owned());
                    boundary_operators.insert(k, matrix);
                }
            }
        }
        ChainComplex {
            complex,
            boundary_operators,
            coboundary_operators,
        }
    }

    /// Verifies the fundamental property of boundary operators: `∂_k ∘ ∂_{k+1} = 0`.
    ///
    /// This means that the boundary of a boundary is always zero, which is a cornerstone
    /// of homology theory.
    ///
    /// # Returns
    /// `true` if the property holds for all relevant dimensions, `false` otherwise.
    pub fn verify_boundary_property(&self) -> bool {
        if let Some(max_dim) = self.complex.dimension() {
            for k in 1..max_dim {
                if let (Some(d_k), Some(d_k_plus_1)) = (
                    self.boundary_operators.get(&k),
                    self.boundary_operators.get(&(k + 1)),
                ) {
                    if (d_k * d_k_plus_1).nnz() != 0 {
                        return false;
                    }
                }
            }
        }
        true
    }

    /// Verifies the fundamental property of coboundary operators: `d^k ∘ d^{k-1} = 0`.
    ///
    /// This is the dual property to the boundary operator property and is fundamental
    /// to cohomology theory.
    ///
    /// # Returns
    /// `true` if the property holds for all relevant dimensions, `false` otherwise.
    pub fn verify_coboundary_property(&self) -> bool {
        if let Some(max_dim) = self.complex.dimension() {
            for k in 1..max_dim {
                if let (Some(d_k), Some(d_k_minus_1)) = (
                    self.coboundary_operators.get(&k),
                    self.coboundary_operators.get(&(k - 1)),
                ) {
                    if (d_k * d_k_minus_1).nnz() != 0 {
                        return false;
                    }
                }
            }
        }
        true
    }

    /// Computes the k-th Betti number, `β_k`, which is a topological invariant.
    ///
    /// The k-th Betti number measures the number of k-dimensional "holes" in a topological space.
    /// It is calculated as `β_k = dim(ker(∂_k)) - dim(im(∂_{k+1}))`.
    ///
    /// # Arguments
    /// * `k` - The dimension `k` for which to compute the Betti number.
    ///
    /// # Returns
    /// An `Option<usize>` containing the k-th Betti number, or `None` if the necessary
    /// boundary operators or simplex counts are not available.
    pub fn compute_homology_betti_number(&self, k: usize) -> Option<usize> {
        let num_k_simplices = self.complex.get_simplices_by_dim(k)?.len();
        let rank_dk = if k == 0 {
            0
        } else {
            self.boundary_operators.get(&k).map_or(0, rank)
        };
        let rank_dk_plus_1 = self.boundary_operators.get(&(k + 1)).map_or(0, rank);
        let dim_ker_k = num_k_simplices.saturating_sub(rank_dk);
        let dim_im_k_plus_1 = rank_dk_plus_1;
        Some(dim_ker_k.saturating_sub(dim_im_k_plus_1))
    }

    /// Computes the k-th cohomology Betti number, `β^k`, which is a topological invariant.
    ///
    /// The k-th cohomology Betti number is dual to the homology Betti number and
    /// is calculated as `β^k = dim(ker(d^k)) - dim(im(d^{k-1}))`.
    ///
    /// # Arguments
    /// * `k` - The dimension `k` for which to compute the cohomology Betti number.
    ///
    /// # Returns
    /// An `Option<usize>` containing the k-th cohomology Betti number, or `None` if the necessary
    /// coboundary operators or simplex counts are not available.
    pub fn compute_cohomology_betti_number(&self, k: usize) -> Option<usize> {
        let num_k_simplices = self.complex.get_simplices_by_dim(k)?.len();
        let rank_dk_t = self.coboundary_operators.get(&k).map_or(0, rank);
        let rank_dk_minus_1_t = if k == 0 {
            0
        } else {
            self.coboundary_operators.get(&(k - 1)).map_or(0, rank)
        };
        let dim_ker_dk = num_k_simplices.saturating_sub(rank_dk_t);
        let dim_im_dk_minus_1 = rank_dk_minus_1_t;
        Some(dim_ker_dk.saturating_sub(dim_im_dk_minus_1))
    }
}

/// Represents a filtration, a sequence of nested simplicial complexes.
pub struct Filtration {
    pub steps: Vec<(f64, SimplicialComplex)>,
}

// =====================================================================================
// endregion: Core Data Structures
// =====================================================================================

// =====================================================================================
// region: Constructors for Common Spaces
// =====================================================================================

// =====================================================================================
// endregion: Core Data Structures
// =====================================================================================

// =====================================================================================
// region: Constructors for Common Spaces
// =====================================================================================

/// Creates a 2D grid simplicial complex.
///
/// This function constructs a simplicial complex that represents a rectangular grid
/// by triangulating each square cell into two triangles.
///
/// # Arguments
/// * `width` - The number of cells in the horizontal direction.
/// * `height` - The number of cells in the vertical direction.
///
/// # Returns
/// A `SimplicialComplex` representing the grid.
pub fn create_grid_complex(width: usize, height: usize) -> SimplicialComplex {
    let mut complex = SimplicialComplex::new();
    for i in 0..height {
        for j in 0..width {
            let v0 = i * (width + 1) + j;
            let v1 = v0 + 1;
            let v2 = (i + 1) * (width + 1) + j;
            let v3 = v2 + 1;
            complex.add_simplex(&[v0, v1, v2]);
            complex.add_simplex(&[v1, v3, v2]);
        }
    }
    complex
}

/// Creates a 2D torus simplicial complex.
///
/// This function constructs a simplicial complex that represents a torus by
/// taking a rectangular grid and identifying opposite edges.
///
/// # Arguments
/// * `m` - The number of cells in one direction (e.g., around the major radius).
/// * `n` - The number of cells in the other direction (e.g., around the minor radius).
///
/// # Returns
/// A `SimplicialComplex` representing the torus.
pub fn create_torus_complex(m: usize, n: usize) -> SimplicialComplex {
    let mut complex = SimplicialComplex::new();
    for i in 0..m {
        for j in 0..n {
            let v0 = i * n + j;
            let v1 = i * n + (j + 1) % n;
            let v2 = ((i + 1) % m) * n + j;
            let v3 = ((i + 1) % m) * n + (j + 1) % n;
            complex.add_simplex(&[v0, v1, v2]);
            complex.add_simplex(&[v1, v3, v2]);
        }
    }
    complex
}

/// Creates a Vietoris-Rips filtration from a set of points in a metric space.
///
/// A Vietoris-Rips complex is a type of simplicial complex constructed from a set of points.
/// A filtration is a sequence of nested complexes, typically built by varying a distance parameter `epsilon`.
/// In a Vietoris-Rips complex, a simplex is included if all its vertices are within `epsilon` distance of each other.
///
/// # Arguments
/// * `points` - A slice of `Vec<f64>` representing the coordinates of the points.
/// * `max_epsilon` - The maximum distance parameter for the filtration.
/// * `steps` - The number of steps (complexes) in the filtration.
///
/// # Returns
/// A `Filtration` containing a sequence of `SimplicialComplex`es.
pub fn vietoris_rips_filtration(points: &[Vec<f64>], max_epsilon: f64, steps: usize) -> Filtration {
    let mut filtration = Filtration { steps: Vec::new() };
    let num_points = points.len();

    for i in 0..=steps {
        let epsilon = max_epsilon * (i as f64 / steps as f64);
        let mut complex = SimplicialComplex::new();
        // Add all points as 0-simplices
        for i in 0..num_points {
            complex.add_simplex(&[i]);
        }
        // Add edges for points within distance epsilon
        for i in 0..num_points {
            for j in (i + 1)..num_points {
                let dist_sq = points[i]
                    .iter()
                    .zip(&points[j])
                    .map(|(a, b)| (a - b).powi(2))
                    .sum::<f64>();
                if dist_sq.sqrt() <= epsilon {
                    complex.add_simplex(&[i, j]);
                }
            }
        }
        // Add all higher-dimensional simplices whose edges are all present (clique complex)
        // This is a simplified version; a full implementation is more complex.
        filtration.steps.push((epsilon, complex));
    }
    filtration
}

// =====================================================================================
// endregion: Constructors for Common Spaces
// =====================================================================================

// Helper function from the sparse module, included here for self-containment.
pub(crate) fn csr_from_triplets(
    rows: usize,
    cols: usize,
    triplets: &[(usize, usize, f64)],
) -> CsMat<f64> {
    let mut mat = TriMat::new((rows, cols));
    for &(r, c, v) in triplets {
        mat.add_triplet(r, c, v);
    }
    mat.to_csr()
}
