//! JIT Instructions for the logic grammar.

use serde::Deserialize;
use serde::Serialize;

#[derive(
    Debug,
    Clone,
    Copy,
    PartialEq,
    Eq,
    Serialize,
    Deserialize,
)]

pub enum JitType {
    I8,
    I16,
    I32,
    I64,
    F32,
    F64,
}

/// Basic instructions for the JIT engine.
/// The JIT operates on a stack-based model.
#[derive(
    Debug, Clone, Serialize, Deserialize,
)]

pub enum Instruction {
    /// Push a 64-bit integer constant onto the stack.
    ImmI(i64),
    /// Push a 64-bit float constant onto the stack.
    ImmF(f64),

    /// Pop address (I64), Load value from address. Push value.
    Load(JitType),
    /// Pop value, Pop address (I64). Store value to address.
    Store(JitType),

    /// Pop rhs, Pop lhs. Push lhs + rhs.
    Add(JitType),
    /// Pop rhs, Pop lhs. Push lhs - rhs.
    Sub(JitType),
    /// Pop rhs, Pop lhs. Push lhs * rhs.
    Mul(JitType),
    /// Pop rhs, Pop lhs. Push lhs / rhs.
    Div(JitType),

    /// Pop rhs, Pop lhs. Push lhs & rhs (Integer only).
    And,
    /// Pop rhs, Pop lhs. Push lhs | rhs (Integer only).
    Or,
    /// Pop rhs, Pop lhs. Push lhs ^ rhs (Integer only).
    Xor,
    /// Pop val. Push !val (Integer only).
    Not,

    /// Comparisons
    Eq(JitType),
    Ne(JitType),
    Lt(JitType),
    Gt(JitType),
    Le(JitType),
    Ge(JitType),

    /// Control Flow
    Label(u32),
    Jump(u32),
    BranchIfTrue(u32),
    BranchIfFalse(u32),

    /// Stack manipulation
    Dup,
    Swap,
    Drop,

    /// Call helper: Pop args_count, Pop function_ptr. Call(fn_ptr, args...).
    /// Note: Assumes signature (args...) -> f64. Arguments must be on stack.
    /// Used for calling helper C functions.
    Call(usize), // arg count

    /// Return the top value of the stack.
    Return,

    /// Custom instruction for user-defined interactions.
    /// `opcode`: Identifier for the custom operation.
    /// `payload`: Static data associated with the instruction.
    Custom {
        opcode: u32,
        payload: Vec<u64>,
    },
}
