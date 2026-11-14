//! Times batch evaluation of one compiled expression through each executor:
//! `cargo run --profile tools --features jit --example backend_timing`.
//!
//! * the reference interpreter (`backend::Interpreter`, the default backend);
//! * the tape interpreter row by row (`Tape::eval` in a loop) and column-wise
//!   (`Tape::eval_columns`, what `TapeFunction` batches use);
//! * Cranelift, one native call per row, and the native row loop that
//!   `JitFunction` compiles from LIR for large batches.

use std::time::Instant;

use rssn::api::Session;
use rssn::backend::Compiled;
use rssn::backend::CompiledMulti;

fn time<F: FnMut()>(
    name: &str,
    reps: u32,
    mut f: F,
) {
    f();
    let start = Instant::now();
    for _ in 0..reps {
        f();
    }
    println!("{name:<34} {:>10.3?}", start.elapsed() / reps);
}

fn main() {
    let s = Session::new();
    let source = "sin(x)*exp(-x^2/2) + x^3 - 2*x/(1 + x^2) + sqrt(x + 1)";
    let f = s.parse(source).expect("parses");
    let n = 100_000;
    let xs: Vec<f64> = (0..n).map(|i| f64::from(i) * 1e-4).collect();
    let mut out = vec![0.0; xs.len()];
    let reps = 20;
    println!("{source}, {n} rows");

    let interpreted = f
        .compile_with(&rssn::backend::Interpreter, &["x"])
        .expect("compiles");
    time("reference interpreter", reps, || {
        interpreted.call_batch(&[&xs], &mut out);
    });
    let reference = out.clone();

    let tape = f.tape(&["x"]).expect("lowers");
    let function = rssn::backend::TapeFunction::from_tape(tape.clone());
    time("tape, row by row", reps, || {
        let mut v = [0.0];
        for (row, slot) in out.iter_mut().enumerate() {
            tape.eval(&[xs[row]], &[], &mut v);
            *slot = v[0];
        }
    });
    time("tape, column blocks", reps, || {
        function.call_batch(&[&xs], &mut out);
    });
    let max_diff = out
        .iter()
        .zip(&reference)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max);
    println!("  (max |tape - reference| = {max_diff:.3e})");

    #[cfg(feature = "jit")]
    {
        let jit = rssn::backend::jit::JitFunction::compile(&tape).expect("compiles");
        let mut columns = vec![0.0; xs.len()];
        time("cranelift, row by row", reps, || {
            let mut v = [0.0];
            for (row, slot) in out.iter_mut().enumerate() {
                jit.eval(&[xs[row]], &[], &mut v);
                *slot = v[0];
            }
        });
        let rows_bits: Vec<u64> = out.iter().map(|x| x.to_bits()).collect();
        time("cranelift, native row loop (LIR)", reps, || {
            jit.call_batch(&[&xs], &mut columns);
        });
        let same = columns.iter().map(|x| x.to_bits()).eq(rows_bits);
        println!("  (native loop bit-identical to per-row: {same})");
    }
}
