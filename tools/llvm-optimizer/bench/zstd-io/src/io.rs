// Compile the production modules and types directly, without unrelated RSSN kernels.
#[path = "../../../../../src/io/frame/column.rs"]
pub mod frame;
#[path = "../../../../../src/io/hdf5/mod.rs"]
pub mod hdf5;
#[path = "../../../../../src/io/parquet/mod.rs"]
pub mod parquet;
