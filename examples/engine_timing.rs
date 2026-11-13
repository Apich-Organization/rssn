//! Times end-to-end symbolic requests (the `benches/engine.rs` set) without
//! criterion's LTO build: `cargo run --profile tools --example engine_timing`.

use rssn::api::Config;
use rssn::api::Session;
use std::time::Instant;

fn main() {
    let s = Session::new();
    let config = Config::new();
    let requests = [
        ("diff", "diff(sin(x)^3*exp(x^2), x)"),
        ("integral", "integral(x^2*exp(x), x)"),
        ("definite_integral", "defint(exp(-x^2), x, -oo, oo)"),
        ("gosper_sum", "sum(k*2^k, k, 0, n)"),
        ("solve_quadratic", "solve(x^2 - 5*x + 6, x)"),
        ("factor", "factor(x^4 - 1)"),
        ("simplify_trig", "sin(x)^2 + cos(x)^2"),
        ("dsolve", "dsolve(diff(diff(y(x), x), x) + y(x) = 0, y(x))"),
        ("expand", "expand((x+y+z)^6)"),
        ("integral_trig", "integral(sin(x)^3*cos(x)^2, x)"),
    ];
    let warm = Instant::now();
    let t = s.parse("x+1").expect("parses");
    let _ = s.compute(t, &config);
    println!("template build + first request: {:?}", warm.elapsed());
    let reps: u32 = std::env::args().nth(1).and_then(|a| a.parse().ok()).unwrap_or(20);
    s.set_cache_policy(rssn::api::CachePolicy::off());
    let no_cache = config.clone();
    let only = std::env::args().nth(2);
    for (name, src) in requests {
        if only.as_deref().is_some_and(|o| o != name) {
            continue;
        }
        let term = s.parse(src).expect("parses");
        let start = Instant::now();
        let mut answer = String::new();
        for _ in 0..reps {
            match s.compute(term, &no_cache) {
                | Ok(a) => answer = a.term.to_string(),
                | Err(e) => answer = format!("error: {e}"),
            }
        }
        let per = start.elapsed() / reps;
        println!("{name:<20} {per:>12.3?}  {}", answer.chars().take(70).collect::<String>());
    }
}
