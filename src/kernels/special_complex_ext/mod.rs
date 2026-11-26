//! Extended complex special-function matrix covering Gamma/Beta, Error/Integrals,
//! Bessel/Airy/Kelvin/Struve, Zeta/Polylog, Elliptic/Modular, and Orthogonal
//! Polynomials / Whittaker / Wright Omega.

pub mod bessel_airy;
pub mod elliptic_modular;
pub mod erf_integrals;
pub mod gamma_beta;
pub mod orthopoly_hyper;
pub mod zeta_polylog;

pub use bessel_airy::*;
pub use elliptic_modular::*;
pub use erf_integrals::*;
pub use gamma_beta::*;
pub use orthopoly_hyper::*;
pub use zeta_polylog::*;
