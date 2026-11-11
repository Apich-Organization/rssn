//! Panic and hang hunt over the dynamic ABI.
//!
//! Every operator the manifest advertises with a plain signature (numbers,
//! booleans, strings, rank 1/2 tensors) is called with degenerate and
//! adversarial arguments: zeros, NaN, infinities, huge magnitudes, empty and
//! 1x1 tensors, and tensors whose `shape` disagrees with their `data`.
//! A kernel may return an error or any value; it must neither panic nor run
//! for more than 20 seconds.
//!
//! A hung call cannot be interrupted, so the driver abandons its thread,
//! records the operator, and continues with a fresh session on a new thread.
//! Known offenders are listed in `KNOWN_PANICS` / `KNOWN_HANGS`; the lists
//! must only ever shrink, and the test fails when an entry stops
//! misbehaving so that it gets removed.

use std::collections::BTreeMap;
use std::collections::BTreeSet;
use std::panic::AssertUnwindSafe;
use std::panic::catch_unwind;
use std::sync::Arc;
use std::sync::Mutex;
use std::time::Duration;
use std::time::Instant;

use num_bigint::BigInt;
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::Value;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::TensorValue;

/// Operators that still panic on adversarial input (documented gaps).
const KNOWN_PANICS: &[&str] = &[];
/// Operators that still take longer than 20 s on adversarial input.
const KNOWN_HANGS: &[&str] = &[];

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }

    fn pick<'a, T>(
        &mut self,
        items: &'a [T],
    ) -> &'a T {
        &items[(self.next() % items.len() as u64) as usize]
    }
}

const FLOATS: &[f64] = &[
    0.0,
    -0.0,
    1.0,
    -1.0,
    0.5,
    2.0,
    3.0,
    10.0,
    -2.5,
    f64::NAN,
    f64::INFINITY,
    f64::NEG_INFINITY,
    1e308,
    -1e308,
    5e-324,
    1e-300,
    1e-8,
    1e8,
    std::f64::consts::PI,
];
const INTS: &[i64] = &[0, 1, -1, 2, 3, 7, 10, -5, 64, 100];

fn float_tensor(
    shape: Vec<usize>,
    rng: &mut Rng,
    consistent: bool,
) -> Value {
    let mut n: usize = shape.iter().product();
    if !consistent {
        n = n.saturating_sub(1);
    }
    let data = (0..n)
        .map(|_| Value::Float(*rng.pick(FLOATS)))
        .collect::<Vec<_>>();
    Value::Tensor(TensorValue {
        shape,
        data,
        elem: TypeExpr::Float,
    })
}

fn make_arg(
    ty: &TypeExpr,
    rng: &mut Rng,
) -> Option<Value> {
    Some(match ty {
        | TypeExpr::Float | TypeExpr::Real => Value::Float(*rng.pick(FLOATS)),
        | TypeExpr::Int => Value::Int(BigInt::from(*rng.pick(INTS))),
        | TypeExpr::Bool => Value::Bool(rng.next() % 2 == 0),
        | TypeExpr::String => {
            Value::String(String::from(*rng.pick(&["", "x", "sin(x)", "((", "1e999"])))
        },
        | TypeExpr::Symbol => Value::Symbol(String::from("x")),
        | TypeExpr::Tensor { rank, .. } => {
            let rank_text = format!("{rank:?}");
            let r = if rank_text.contains('2') {
                2
            } else if rank_text.contains('1') {
                1
            } else {
                (rng.next() % 2 + 1) as usize
            };
            let dims = [0usize, 1, 2, 3, 4];
            let shape: Vec<usize> = (0..r).map(|_| *rng.pick(&dims)).collect();
            let consistent = rng.next() % 6 != 0;
            float_tensor(shape, rng, consistent)
        },
        | _ => return None,
    })
}

#[derive(Default)]
struct Progress {
    generation: u32,
    done: BTreeSet<String>,
    hung: BTreeSet<String>,
    /// Operator -> first panic message and arguments.
    panicked: BTreeMap<String, String>,
    current: Option<(String, Instant)>,
    finished: bool,
    calls: usize,
}

fn spawn_worker(
    progress: &Arc<Mutex<Progress>>,
    generation: u32,
) {
    let progress = Arc::clone(progress);
    std::thread::Builder::new()
        .stack_size(64 * 1024 * 1024)
        .spawn(move || run(&progress, generation))
        .expect("spawn");
}

#[test]
fn abi_survives_adversarial_arguments() {
    std::panic::set_hook(Box::new(|_| {}));
    let progress = Arc::new(Mutex::new(Progress::default()));
    spawn_worker(&progress, 0);
    let mut generation = 0;
    loop {
        std::thread::sleep(Duration::from_millis(200));
        let mut p = progress.lock().unwrap();
        if p.finished {
            break;
        }
        if let Some((name, started)) = p.current.clone()
            && started.elapsed() > Duration::from_secs(20)
        {
            p.hung.insert(name.clone());
            p.done.insert(name);
            p.current = None;
            generation += 1;
            p.generation = generation;
            drop(p);
            spawn_worker(&progress, generation);
        }
    }
    let _ = std::panic::take_hook();
    let p = progress.lock().unwrap();
    eprintln!(
        "abi_fuzz: {} calls over {} operators, {} panicked, {} hung",
        p.calls,
        p.done.len(),
        p.panicked.len(),
        p.hung.len()
    );
    for (name, message) in &p.panicked {
        eprintln!("  panic {name}: {message}");
    }
    for name in &p.hung {
        eprintln!("  hang  {name}");
    }
    let unexpected_panics: Vec<_> = p
        .panicked
        .iter()
        .filter(|(name, _)| !KNOWN_PANICS.contains(&name.as_str()))
        .collect();
    let unexpected_hangs: Vec<_> = p
        .hung
        .iter()
        .filter(|name| !KNOWN_HANGS.contains(&name.as_str()))
        .collect();
    assert!(
        unexpected_panics.is_empty() && unexpected_hangs.is_empty(),
        "operators that panic: {unexpected_panics:#?}\noperators that hang: {unexpected_hangs:#?}"
    );
    let stale: Vec<_> = KNOWN_PANICS
        .iter()
        .filter(|name| !p.panicked.contains_key(**name))
        .chain(KNOWN_HANGS.iter().filter(|name| !p.hung.contains(**name)))
        .collect();
    assert!(stale.is_empty(), "stale known-offender entries: {stale:?}");
}

fn run(
    progress: &Arc<Mutex<Progress>>,
    generation: u32,
) {
    let session = Session::new();
    let manifest = session.manifest();
    for op in &manifest.ops {
        {
            let p = progress.lock().unwrap();
            if p.generation != generation {
                return;
            }
            if p.done.contains(&op.name) {
                continue;
            }
        }
        for (index, sig) in op.overloads.iter().enumerate() {
            if sig.variadic.is_some() {
                continue;
            }
            let mut rng = Rng(0x9E37_79B9_7F4A_7C15 ^ (op.name.len() as u64 + 1) * 7919);
            for _ in 0..40 {
                let args: Option<Vec<Value>> =
                    sig.params.iter().map(|t| make_arg(t, &mut rng)).collect();
                let Some(args) = args else { break };
                let Some(key) = op.keys.get(index).or_else(|| op.keys.first()) else {
                    break;
                };
                {
                    let mut p = progress.lock().unwrap();
                    if p.generation != generation {
                        return;
                    }
                    p.current = Some((op.name.clone(), Instant::now()));
                }
                let result = catch_unwind(AssertUnwindSafe(|| {
                    session.call_key(key, &args, &CallCtx::new())
                }));
                let mut p = progress.lock().unwrap();
                if p.generation != generation {
                    return;
                }
                p.calls += 1;
                if let Err(payload) = result {
                    let text = payload
                        .downcast_ref::<String>()
                        .cloned()
                        .or_else(|| payload.downcast_ref::<&str>().map(|s| (*s).to_string()))
                        .unwrap_or_default();
                    let shown = format!("{args:?}");
                    p.panicked
                        .entry(op.name.clone())
                        .or_insert_with(|| format!("{text} <- {}", &shown[..shown.len().min(200)]));
                }
            }
        }
        let mut p = progress.lock().unwrap();
        if p.generation != generation {
            return;
        }
        p.done.insert(op.name.clone());
        p.current = None;
    }
    let mut p = progress.lock().unwrap();
    if p.generation == generation {
        p.finished = true;
    }
}
