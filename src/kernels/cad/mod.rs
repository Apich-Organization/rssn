//! Geometry import and repair: STL and STEP readers, mesh repair and Boolean operations on triangle meshes.

pub mod csg;
pub mod repair;
pub mod step;
pub mod stl;
pub mod types;

pub use csg::BoolOp;
pub use csg::boolean;
pub use repair::RepairReport;
pub use repair::euler_characteristic;
pub use repair::is_closed;
pub use repair::is_manifold;
pub use repair::repair;
pub use repair::signed_volume;
pub use repair::weld;
pub use step::read_step;
pub use step::read_step_file;
pub use stl::read_stl;
pub use stl::write_stl_ascii;
pub use stl::write_stl_binary;
pub use types::Mesh;
