//! Survival analysis: Kaplan-Meier, Nelson-Aalen, log-rank tests and Cox regression.
//!
//! The estimators work on paired slices of observation times and event flags,
//! where a false flag marks a censored observation. None of the routines modify
//! their inputs; ties are grouped internally on sorted copies.

mod common;

pub mod cox;
pub mod kaplan_meier;
pub mod log_rank;
pub mod nelson_aalen;

pub use cox::CoxFit;
pub use cox::Ties;
pub use cox::cox_fit;
pub use cox::cox_predict;
pub use kaplan_meier::KmEstimate;
pub use kaplan_meier::kaplan_meier;
pub use kaplan_meier::median_survival;
pub use log_rank::LogRank;
pub use log_rank::log_rank;
pub use nelson_aalen::NaEstimate;
pub use nelson_aalen::nelson_aalen;
