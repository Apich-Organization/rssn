//! Integration tests for `rssn::kernels`.

// Numeric tests index several arrays in step; index loops read closest to
// the formulas they check.
#![allow(clippy::needless_range_loop)]

#[path = "kernels/amg_general.rs"]
mod amg_general;
#[path = "kernels/andrews_squeezer.rs"]
mod andrews_squeezer;
#[path = "kernels/approx.rs"]
mod approx;
#[path = "kernels/calculus.rs"]
mod calculus;
#[path = "kernels/combinatorics.rs"]
mod combinatorics;
#[path = "kernels/computer_graphics.rs"]
mod computer_graphics;
#[path = "kernels/convergence.rs"]
mod convergence;
#[path = "kernels/coordinates.rs"]
mod coordinates;
#[path = "kernels/dd.rs"]
mod dd;
#[path = "kernels/dense.rs"]
mod dense;
#[path = "kernels/differential_geometry.rs"]
mod differential_geometry;
#[path = "kernels/eigs_ext.rs"]
mod eigs_ext;
#[path = "kernels/error_correction.rs"]
mod error_correction;
#[path = "kernels/fftconv.rs"]
mod fftconv;
#[path = "kernels/finite_field.rs"]
mod finite_field;
#[path = "kernels/fractal_geometry_and_chaos.rs"]
mod fractal_geometry_and_chaos;
#[path = "kernels/functional_analysis.rs"]
mod functional_analysis;
#[path = "kernels/geometric_algebra.rs"]
mod geometric_algebra;
#[path = "kernels/graph.rs"]
mod graph;
#[path = "kernels/homotopy.rs"]
mod homotopy;
#[path = "kernels/homotopy_cert.rs"]
mod homotopy_cert;
#[path = "kernels/hypergeometric.rs"]
mod hypergeometric;
#[path = "kernels/indefinite_sum.rs"]
mod indefinite_sum;
#[path = "kernels/integrate.rs"]
mod integrate;
#[path = "kernels/interpolate.rs"]
mod interpolate;
#[path = "kernels/krylov.rs"]
mod krylov;
#[path = "kernels/krylov_ext.rs"]
mod krylov_ext;
#[path = "kernels/matfun.rs"]
mod matfun;
#[path = "kernels/matrix.rs"]
mod matrix;
#[path = "kernels/multigrid.rs"]
mod multigrid;
#[path = "kernels/number_theory.rs"]
mod number_theory;
#[path = "kernels/ode.rs"]
mod ode;
#[path = "kernels/ode_adaptive.rs"]
mod ode_adaptive;
#[path = "kernels/ode_ext.rs"]
mod ode_ext;
#[path = "kernels/ode_stiff.rs"]
mod ode_stiff;
#[path = "kernels/ode_stiff_ext.rs"]
mod ode_stiff_ext;
#[path = "kernels/ode_stiff_struct.rs"]
mod ode_stiff_struct;
#[path = "kernels/optim.rs"]
mod optim;
#[path = "kernels/optimize.rs"]
mod optimize;
#[path = "kernels/ordering.rs"]
mod ordering;
#[path = "kernels/pde.rs"]
mod pde;
#[path = "kernels/polynomial.rs"]
mod polynomial;
#[path = "kernels/quadrature.rs"]
mod quadrature;
#[path = "kernels/random.rs"]
mod random;
#[path = "kernels/real_roots.rs"]
mod real_roots;
#[path = "kernels/robust.rs"]
mod robust;
#[path = "kernels/rootfind.rs"]
mod rootfind;
#[path = "kernels/sde.rs"]
mod sde;
#[path = "kernels/series.rs"]
mod series;
#[path = "kernels/signal.rs"]
mod signal;
#[path = "kernels/solve.rs"]
mod solve;
#[path = "kernels/sparse.rs"]
mod sparse;
#[path = "kernels/special.rs"]
mod special;
#[path = "kernels/special_ext.rs"]
mod special_ext;
#[path = "kernels/special_hyp.rs"]
mod special_hyp;
#[path = "kernels/special_limits.rs"]
mod special_limits;
#[path = "kernels/stats.rs"]
mod stats;
#[path = "kernels/tensor.rs"]
mod tensor;
#[path = "kernels/topology.rs"]
mod topology;
#[path = "kernels/transforms.rs"]
mod transforms;
#[path = "kernels/vector.rs"]
mod vector;
#[path = "kernels/wavelets.rs"]
mod wavelets;

#[path = "kernels/par_sites.rs"]
mod par_sites;
