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
/// The type of a JIT value.

pub enum JitType {
    /// 8-bit integer.
    I8,
    /// 16-bit integer.
    I16,
    /// 32-bit integer.
    I32,
    /// 64-bit integer.
    I64,
    /// 32-bit float.
    F32,
    /// 64-bit float.
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
    /// Pop rhs, Pop lhs. Push 1 if lhs == rhs else 0.
    Eq(JitType),
    /// Pop rhs, Pop lhs. Push 1 if lhs != rhs else 0.
    Ne(JitType),
    /// Pop rhs, Pop lhs. Push 1 if lhs < rhs else 0.
    Lt(JitType),
    /// Pop rhs, Pop lhs. Push 1 if lhs > rhs else 0.
    Gt(JitType),
    /// Pop rhs, Pop lhs. Push 1 if lhs <= rhs else 0.
    Le(JitType),
    /// Pop rhs, Pop lhs. Push 1 if lhs >= rhs else 0.
    Ge(JitType),

    /// Control Flow
    /// A label to jump to.
    Label(u32),
    /// Unconditionally jump to a label.
    Jump(u32),
    /// Pop value. Jump to label if value is not 0.
    BranchIfTrue(u32),
    /// Pop value. Jump to label if value is 0.
    BranchIfFalse(u32),

    /// Stack manipulation
    /// Duplicate the top value on the stack.
    Dup,
    /// Swap the top two values on the stack.
    Swap,
    /// Drop the top value from the stack.
    Drop,

    /// Call helper: Pop `args_count`, Pop `function_ptr`. `Call(fn_ptr`, args...).
    /// Note: Assumes signature (args...) -> f64. Arguments must be on stack.
    /// Used for calling helper C functions.
    Call(usize), // arg count

    /// Return the top value of the stack.
    Return,

    /// Custom instruction for user-defined interactions.
    Custom {
        /// Identifier for the custom operation.
        opcode: u32,
        /// Static data associated with the instruction.
        payload: Vec<u64>,
    },
}
