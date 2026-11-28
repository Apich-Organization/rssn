//! Mesh generation: constrained Delaunay triangulation, quality refinement, tetrahedralisation and advancing front.

pub mod advancing_front;
pub mod delaunay;
pub mod refine;
pub mod tetra;
pub mod util;

pub use advancing_front::advancing_front_2d;
pub use delaunay::constrained_delaunay_2d;
pub use delaunay::delaunay_2d;
pub use delaunay::delaunay_2d_incremental;
pub use refine::refine_quality;
pub use tetra::TetMesh;
pub use tetra::boundary_faces;
pub use tetra::delaunay_3d;
pub use tetra::tetra_volumes;
pub use util::TriMesh2;
pub use util::min_angle_deg;
pub use util::triangle_areas;
