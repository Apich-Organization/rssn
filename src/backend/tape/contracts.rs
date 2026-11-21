//! Validation of public SSA tapes before native memory access, and exact
//! cache identity independent of the hash used to choose a cache bucket.
use super::BackendError;
use super::Inst;
#[cfg(feature = "jit")]
use super::ProgramKey;
use super::Tape;
use super::key_of;
use super::operands;

impl Tape {
    /// Validate SSA dependencies, input/parameter indices, calls and outputs.
    ///
    /// Every operand must refer to an earlier instruction; each input and
    /// parameter must lie within the declared arity. Call argument ranges
    /// must exist and all outputs must name instructions. Gather offsets are
    /// checked at runtime, since their values are data dependent.
    ///
    /// # Errors
    /// Returns [`BackendError::Codegen`] describing the invalid instruction or
    /// output. Native backends call this before generating any memory loads.
    pub fn validate(&self) -> Result<(), BackendError> {
        if self.insts.len() > u32::MAX as usize {
            return Err(BackendError::Codegen(
                "tape exceeds register address space".into(),
            ));
        }
        for (index, inst) in self.insts.iter().enumerate() {
            let invalid =
                |reason: &str| BackendError::Codegen(format!("instruction {index}: {reason}"));
            match *inst {
                | Inst::Input(i) if i as usize >= self.inputs => {
                    return Err(invalid("input index exceeds declared arity"));
                },
                | Inst::Param(i) if i as usize >= self.params => {
                    return Err(invalid("parameter index exceeds declared arity"));
                },
                | Inst::Call(_, start, len) => {
                    let end = (start as usize)
                        .checked_add(len as usize)
                        .ok_or_else(|| invalid("call argument range overflow"))?;
                    if self.args.get(start as usize..end).is_none() {
                        return Err(invalid("call argument range does not exist"));
                    }
                },
                | _ => {},
            }
            if operands(inst, &self.args)
                .iter()
                .any(|&r| r as usize >= index)
            {
                return Err(invalid("operand is not an earlier SSA register"));
            }
        }
        for (index, &reg) in self.outputs.iter().enumerate() {
            if reg as usize >= self.insts.len() {
                return Err(BackendError::Codegen(format!(
                    "output {index}: register does not exist"
                )));
            }
        }
        Ok(())
    }

    #[cfg(feature = "jit")]
    pub(in crate::backend) fn program_key(&self) -> ProgramKey {
        ProgramKey {
            inputs: self.inputs,
            params: self.params,
            insts: self.insts.iter().map(ordered_key_of).collect(),
            args: self.args.clone(),
            outputs: self.outputs.clone(),
        }
    }
}

// IEEE operand order and constant bit patterns are part of the program.
// Builder commoning may use canonical arithmetic keys; a native code cache
// must retain the exact instruction order instead.
pub(super) fn ordered_key_of(inst: &Inst) -> (u8, u64, u64, u64) {
    match *inst {
        | Inst::Add(a, b) => (3, u64::from(a), u64::from(b), 0),
        | Inst::Mul(a, b) => (5, u64::from(a), u64::from(b), 0),
        | _ => key_of(inst),
    }
}
