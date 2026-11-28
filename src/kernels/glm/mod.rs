//! Generalised linear models fitted by iteratively reweighted least squares.
//!
//! The module separates the three ingredients of a generalised linear model:
//! a response [`Family`] fixing the mean–variance relationship, a [`Link`]
//! mapping the linear predictor to the mean, and the fitting loop in
//! [`fit`], which returns a [`GlmFit`]. [`predict`] applies a fitted model to
//! new covariates.
//!
//! [`Family`]: crate::kernels::glm::Family
//! [`Link`]: crate::kernels::glm::Link
//! [`fit`]: fn@crate::kernels::glm::fit
//! [`GlmFit`]: crate::kernels::glm::GlmFit
//! [`predict`]: fn@crate::kernels::glm::predict

pub mod family;
pub mod fit;
pub mod link;
pub mod predict;

pub use family::Family;
pub use family::variance;
pub use fit::GlmFit;
pub use fit::fit;
pub use link::Link;
pub use link::link;
pub use link::link_derivative;
pub use link::link_inverse;
pub use predict::predict;
