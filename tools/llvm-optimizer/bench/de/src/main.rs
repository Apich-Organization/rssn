use std::hint::black_box;
use std::time::Instant;

#[cfg(any(target_os = "linux", windows))]
#[path = "../../device-lanes/clock.rs"]
mod clock;

#[cfg(not(any(target_os = "linux", windows)))]
mod clock {
    pub(super) fn now() -> u64 {
        panic!("thread CPU timing is available on Linux and Windows");
    }
}

use rssn_de_study::rssn_de_appell;
use rssn_de_study::rssn_de_negative;
use rssn_de_study::rssn_de_positive;
use rssn_de_study::rssn_de_small;

type Kernel = extern "C" fn(f64, f64, f64) -> f64;
type Case = (f64, f64, f64);

fn main() {
    let mut args = std::env::args().skip(1);
    let case = args.next().expect("small, positive, negative, or appell");
    let rounds: u64 = args.next().expect("round count").parse().unwrap();
    let thread_cpu = match args.next().as_deref() {
        | None => false,
        | Some("--thread-cpu") => true,
        | _ => panic!("expected --thread-cpu"),
    };
    if case == "appell" {
        measure(
            &[
                (0.8, 1.25, 0.5, 2.1, 0.96, 0.3),
                (1.2, 0.25, 2.0, 2.5, 0.4, 0.99),
                (0.05, 2.0, 0.25, 0.9, -1.2, 0.97),
                (1.55, 0.05, 0.25, 2.1, -2.1, -3.1),
            ],
            rounds,
            thread_cpu,
            |(a, b1, b2, c, x, y)| rssn_de_appell(a, b1, b2, c, x, y),
        );
        return;
    }

    let (kernel, data): (Kernel, &[Case]) = match case.as_str() {
        | "small" => {
            (
                rssn_de_small,
                &[
                    (-2.934, 0.09, 1.0),
                    (0.999, 0.001, 0.001),
                    (-1234.0, 0.1, 1.2),
                    (0.9999999, 0.12, 1.0),
                ],
            )
        },
        | "positive" => {
            (
                rssn_de_positive,
                &[
                    (-2.934, 0.805, 1.0),
                    (0.999, 0.805, 1.0),
                    (-6.5, 0.5, 1.7),
                    (0.999, 0.15, 0.001),
                ],
            )
        },
        | "negative" => {
            (
                rssn_de_negative,
                &[
                    (-1.0, -0.25, 1.0),
                    (0.95, -0.5, 1.2),
                    (-4.0, -3.5, 0.7),
                    (0.99, -2.25, 1.7),
                ],
            )
        },
        | _ => panic!("unknown case"),
    };
    let kernel = black_box(kernel);
    measure(data, rounds, thread_cpu, |(z, s, a)| kernel(z, s, a));
}

fn measure<T: Copy>(
    data: &[T],
    rounds: u64,
    thread_cpu: bool,
    kernel: impl Fn(T) -> f64,
) {
    for &args in data {
        assert!(kernel(args).is_finite());
    }
    let mut hash = 0_u64;
    let wall_start = Instant::now();
    let cpu_start = if thread_cpu {
        clock::now()
    } else {
        0
    };
    for _ in 0..rounds {
        for &args in black_box(data) {
            let bits = kernel(args).to_bits();
            hash = (hash.rotate_left(7) ^ bits).wrapping_mul(0x9e3779b97f4a7c15);
        }
    }
    let elapsed = if thread_cpu {
        u128::from(clock::now() - cpu_start)
    } else {
        wall_start.elapsed().as_nanos()
    };
    println!("{elapsed} {hash:016x}");
}
