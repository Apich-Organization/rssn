//! Convex optimisation: dense and sparse quadratic programs, second-order cone and semidefinite programs.

pub mod lp;
pub mod qp;
pub mod sdp;
pub mod socp;

pub use lp::LpProblem;
pub use lp::LpSolution;
pub use lp::LpStatus;
pub use lp::Relation;
pub use lp::solve_lp;
pub use qp::QpProblem;
pub use qp::QpSolution;
pub use qp::solve_qp;
pub use sdp::SdpProblem;
pub use sdp::SdpSolution;
pub use sdp::solve_sdp;
pub use socp::SocCone;
pub use socp::SocpProblem;
pub use socp::SocpSolution;
pub use socp::solve_socp;
