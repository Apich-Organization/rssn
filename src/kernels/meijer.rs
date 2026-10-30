//! Meijer G-function evaluation.
//!
//! The kernel evaluates the Meijer G-function
//! `G^{m,n}_{p,q}(z | a_1..a_p; b_1..b_q)` for complex parameters and
//! argument, by Slater's contour theorem (residue series) and, where the
//! residue series does not apply, by numerical Mellin-Barnes integration.
//!
//! The symbolic layer lives in `rules::special::meijer`; this module is
//! plain numerics and is also usable directly.
