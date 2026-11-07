//! Native memory access, optimizer and cache contracts of public Tape IR.
#![cfg(feature = "jit")]
use rssn::backend::Compiled;
use rssn::backend::CompiledMulti;
use rssn::backend::jit::CraneliftBackend;
use rssn::backend::jit::JitFunction;
use rssn::backend::tape::Inst;
use rssn::backend::tape::Tape;
use std::sync::Arc;
use std::sync::Barrier;

fn same(
    a: f64,
    b: f64,
) {
    assert!(
        a.to_bits() == b.to_bits() || (a.is_nan() && b.is_nan()),
        "{a:?} != {b:?}"
    );
}

#[test]
fn constant_gather_is_not_constant_folded_and_still_reads_input_data() {
    let tape = Tape {
        insts: vec![Inst::Const(1.0), Inst::Gather(0, 0)],
        inputs: 3,
        outputs: vec![1],
        ..Tape::default()
    };
    tape.validate().unwrap();
    let optimized = tape.optimise();
    for program in [&tape, &optimized] {
        let native = JitFunction::compile(program).unwrap();
        for input in [[1.0, 7.0, 3.0], [1.0, 19.0, 3.0]] {
            let mut expected = [0.0];
            program.eval(&input, &[], &mut expected);
            assert_eq!(expected[0], input[1]);
            same(native.call(&input), expected[0]);
        }
    }
}

#[test]
fn gather_truncation_bounds_padding_and_extra_arguments_match_interpreter() {
    let tape = Tape {
        insts: vec![Inst::Input(0), Inst::Gather(1, 0)],
        inputs: 4,
        outputs: vec![1],
        ..Tape::default()
    };
    let native = JitFunction::compile(&tape).unwrap();
    for index in [
        -1.0,
        -0.9,
        -0.0,
        0.0,
        1.8,
        2.999,
        3.0,
        f64::NAN,
        f64::INFINITY,
        f64::NEG_INFINITY,
    ] {
        for input in [
            vec![index, 10.0, 20.0, 30.0],
            vec![index],
            vec![index, 10.0, 20.0, 30.0, 999.0],
        ] {
            let mut expected = [0.0];
            tape.eval(&input, &[], &mut expected);
            same(native.call(&input), expected[0]);
        }
    }
}

#[test]
fn malformed_tapes_are_rejected_before_any_native_memory_access() {
    fn opaque(x: &[f64]) -> f64 {
        x.iter().sum()
    }
    let malformed = [
        Tape {
            insts: vec![Inst::Input(1)],
            inputs: 1,
            outputs: vec![0],
            ..Tape::default()
        },
        Tape {
            insts: vec![Inst::Param(1)],
            params: 1,
            outputs: vec![0],
            ..Tape::default()
        },
        Tape {
            insts: vec![Inst::Neg(0)],
            outputs: vec![0],
            ..Tape::default()
        },
        Tape {
            insts: vec![Inst::Neg(1), Inst::Const(1.0)],
            outputs: vec![0],
            ..Tape::default()
        },
        Tape {
            insts: vec![Inst::Const(1.0)],
            outputs: vec![1],
            ..Tape::default()
        },
        Tape {
            insts: vec![Inst::Call(opaque, 0, 1)],
            outputs: vec![0],
            ..Tape::default()
        },
        Tape {
            insts: vec![Inst::Call(opaque, 0, 1)],
            args: vec![0],
            outputs: vec![0],
            ..Tape::default()
        },
    ];
    for tape in malformed {
        assert!(tape.validate().is_err());
        assert!(JitFunction::compile(&tape).is_err());
    }
}

#[test]
fn short_output_buffers_copy_the_computed_prefix_and_preserve_extra_slots() {
    let tape = Tape {
        insts: vec![Inst::Const(2.0), Inst::Const(3.0)],
        outputs: vec![0, 1],
        ..Tape::default()
    };
    let native = JitFunction::compile(&tape).unwrap();
    assert_eq!(native.call(&[]), 2.0);
    for n in [0, 1, 2, 3] {
        let mut expected = vec![99.0; n];
        let mut actual = expected.clone();
        tape.eval(&[], &[], &mut expected);
        native.eval(&[], &[], &mut actual);
        assert_eq!(actual, expected);
    }
    let empty = Tape::default();
    let native = JitFunction::compile(&empty).unwrap();
    assert!(native.call(&[]).is_nan());
}

#[test]
fn fingerprints_retain_ordered_operands_and_ieee_constant_bits() {
    let a = Tape {
        insts: vec![Inst::Input(0), Inst::Input(1), Inst::Add(0, 1)],
        inputs: 2,
        outputs: vec![2],
        ..Tape::default()
    };
    let mut b = a.clone();
    b.insts[2] = Inst::Add(1, 0);
    assert_ne!(a.structural_hash(), b.structural_hash());
    let backend = CraneliftBackend::new();
    let mut functions = Vec::new();
    for value in [
        0.0,
        -0.0,
        f64::from_bits(0x7ff8_0000_0000_0001),
        f64::from_bits(0x7ff8_0000_0000_0002),
    ] {
        let tape = Tape {
            insts: vec![Inst::Const(value)],
            outputs: vec![0],
            ..Tape::default()
        };
        let function = backend.compile_tape(&tape).unwrap();
        assert_eq!(function.call(&[]).to_bits(), value.to_bits());
        functions.push(function);
    }
    assert_eq!(backend.cached(), 4);
    for (i, a) in functions.iter().enumerate() {
        for b in &functions[i + 1..] {
            assert!(!Arc::ptr_eq(a, b));
        }
    }
}

#[test]
fn simultaneous_requests_share_one_native_code_owner() {
    let backend = Arc::new(CraneliftBackend::new());
    let barrier = Arc::new(Barrier::new(8));
    let tape = Tape {
        insts: vec![Inst::Input(0), Inst::Mul(0, 0)],
        inputs: 1,
        outputs: vec![1],
        ..Tape::default()
    };
    let handles: Vec<_> = (0..8)
        .map(|_| {
            let backend = Arc::clone(&backend);
            let barrier = Arc::clone(&barrier);
            let tape = tape.clone();
            std::thread::spawn(move || {
                barrier.wait();
                backend.compile_tape(&tape).unwrap()
            })
        })
        .collect();
    let compiled: Vec<_> = handles.into_iter().map(|h| h.join().unwrap()).collect();
    for function in &compiled {
        assert!(Arc::ptr_eq(&compiled[0], function));
        assert_eq!(function.call(&[3.0]), 9.0);
    }
    assert_eq!(backend.cached(), 1);
}

#[test]
fn failed_entries_do_not_count_as_compiled_and_context_changes_keep_old_code_alive() {
    let mut backend = CraneliftBackend::new();
    let bad = Tape {
        insts: vec![Inst::Input(4)],
        inputs: 1,
        outputs: vec![0],
        ..Tape::default()
    };
    for _ in 0..2 {
        assert!(backend.compile_tape(&bad).is_err());
    }
    assert_eq!(backend.cached(), 0);
    let good = Tape {
        insts: vec![Inst::Input(0)],
        inputs: 1,
        outputs: vec![0],
        ..Tape::default()
    };
    let old = backend.compile_tape(&good).unwrap();
    backend.set_key_extra(42);
    assert_eq!(backend.cached(), 0);
    let new = backend.compile_tape(&good).unwrap();
    assert!(!Arc::ptr_eq(&old, &new));
    assert_eq!(old.call(&[4.0]), 4.0);
    assert_eq!(new.call(&[5.0]), 5.0);
}
