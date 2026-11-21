//! Compiling a closed form through the `Backend` trait, calling it, and
//! taking its gradient from the tape.
//!
//! Run with `cargo run --release --example compile_expression`.

use std::sync::Mutex;

use rssn::api::*;
use rssn::backend::Backend;
use rssn::backend::BackendError;
use rssn::backend::Compiled;
use rssn::backend::tape::Tape;
use rssn::backend::tape::{
    self,
};
use rssn::graph::Graph;
use rssn::graph::NodeId;
use rssn::graph::SymbolId;

/// A backend that lowers to the reference tape and keeps a copy, so the
/// example can ask the same program for its gradient. `Interpreter` does the
/// same lowering; it just does not hand the tape back.
struct Taped(Mutex<Option<Tape>>);

struct Program(Tape);

impl Compiled for Program {
    fn arity(&self) -> usize {
        self.0.inputs
    }

    fn call(
        &self,
        args: &[f64],
    ) -> f64 {
        let mut out = [f64::NAN];
        self.0.eval(args, &[], &mut out);
        out[0]
    }
}

impl Backend for Taped {
    fn name(&self) -> &'static str {
        "tape"
    }

    fn compile(
        &self,
        graph: &Graph,
        root: NodeId,
        inputs: &[SymbolId],
    ) -> Result<Box<dyn Compiled>, BackendError> {
        let tape = tape::lower(graph, &[root], inputs, &[])?.optimise();
        let mut slot = self
            .0
            .lock()
            .map_err(|_| BackendError::Codegen("the tape slot is poisoned".to_owned()))?;
        *slot = Some(tape.clone());
        Ok(Box::new(Program(tape)))
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let s = Session::new();
    // Shared subterms survive the lowering: exp(x*y) is computed once a call.
    let f = s.parse("x^2*sin(y) + exp(x*y)")?;
    let (px, py) = (0.5, 1.2);

    // The default path: the interpreter behind `Compiled`. With the `jit`
    // feature `backend::set_default` swaps in Cranelift; call sites do not
    // change.
    let plain = f.compile(&["x", "y"])?;
    println!("interpreter f({px}, {py}) = {:.12}", plain.call(&[px, py]));

    let backend = Taped(Mutex::new(None));
    let compiled = f.compile_with(&backend, &["x", "y"])?;
    println!(
        "tape        f({px}, {py}) = {:.12}",
        compiled.call(&[px, py])
    );

    // Forward-mode gradient from the tape, against the closed form:
    // df/dx = 2x sin y + y e^(xy),  df/dy = x^2 cos y + x e^(xy).
    let gradient = {
        let guard = backend.0.lock().map_err(|_| "the tape slot is poisoned")?;
        let tape = guard.as_ref().ok_or("the backend kept no tape")?;
        tape.gradient(0, &[px, py], &[])
    };
    let want = [
        2.0 * px * py.sin() + py * (px * py).exp(),
        px * px * py.cos() + px * (px * py).exp(),
    ];
    println!("  df/dx = {:+.9}  (analytic {:+.9})", gradient[0], want[0]);
    println!("  df/dy = {:+.9}  (analytic {:+.9})", gradient[1], want[1]);

    // Compilation needs every free symbol: ask for f(x) alone and it says
    // which one is missing rather than returning a wrong number.
    match f.compile(&["x"]) {
        | Err(e) => println!("\nf as f(x) alone: {e}"),
        | Ok(_) => println!("\nf as f(x) alone: compiled (unexpected)"),
    }

    // One program, many points: the point of compiling.
    let xs: Vec<f64> = (0..5).map(|i| 0.1 * f64::from(i)).collect();
    let ys = vec![py; xs.len()];
    let mut out = vec![0.0; xs.len()];
    compiled.call_batch(&[&xs[..], &ys[..]], &mut out);
    println!("batch over x: {out:?}");
    Ok(())
}
