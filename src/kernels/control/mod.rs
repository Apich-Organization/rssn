//! Control systems: state space, transfer functions, Riccati equations, LQR/LQG/MPC, PID tuning, frequency response and system identification.

mod linalg;

pub mod ident;
pub mod mpc;
pub mod pid;
pub mod riccati;
pub mod state_space;
pub mod transfer;

pub use ident::ArxModel;
pub use ident::EraModel;
pub use ident::arx;
pub use ident::era;
pub use ident::impulse_response;
pub use mpc::MpcOptions;
pub use mpc::MpcStep;
pub use mpc::mpc_closed_loop;
pub use mpc::mpc_step;
pub use pid::Pid;
pub use pid::PidForm;
pub use pid::TuningRule;
pub use pid::pid_response;
pub use pid::simulate_pid;
pub use pid::tune;
pub use riccati::care;
pub use riccati::dare;
pub use riccati::dlqr;
pub use riccati::kalman_gain;
pub use riccati::lqg;
pub use riccati::lqr;
pub use state_space::StateSpace;
pub use transfer::TransferFunction;
pub use transfer::bode;
pub use transfer::feedback;
pub use transfer::gain_margin_phase_margin;
pub use transfer::nyquist;
pub use transfer::parallel;
pub use transfer::series;
pub use transfer::ss_to_tf;
pub use transfer::step_response;
