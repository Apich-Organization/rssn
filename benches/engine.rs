//! Benchmarks of the compute engine: symbolic requests end to end, numeric
//! evaluation, and a compiled function evaluated in batch.

use criterion::BenchmarkId;
use criterion::Criterion;
use criterion::criterion_group;
use criterion::criterion_main;
use rssn::api::CachePolicy;
use rssn::api::Config;
use rssn::api::Session;
use std::hint::black_box;

fn symbolic(c: &mut Criterion) {
    let s = Session::new();
    s.set_cache_policy(CachePolicy::off());
    let config = Config::new();
    let mut group = c.benchmark_group("symbolic");
    for (name, request) in [
        ("diff", "diff(sin(x)^3*exp(x^2), x)"),
        ("integral", "integral(x^2*exp(x), x)"),
        // The Risch–Bronstein tower: a nested logarithm, and a repeated
        // pole in an exponential extension.
        (
            "integral_tower_nested",
            "integral(log(log(x))/(x*log(x)), x)",
        ),
        ("integral_tower_pole", "integral(1/(exp(x) + 1)^2, x)"),
        ("definite_integral", "defint(exp(-x^2), x, -oo, oo)"),
        ("gosper_sum", "sum(k*2^k, k, 0, n)"),
        ("solve_quadratic", "solve(x^2 - 5*x + 6, x)"),
        ("factor", "factor(x^4 - 1)"),
        ("simplify_trig", "sin(x)^2 + cos(x)^2"),
        ("dsolve", "dsolve(diff(diff(y(x), x), x) + y(x) = 0, y(x))"),
    ] {
        let term = s.parse(request).expect("benchmark request parses");
        group.bench_function(name, |b| b.iter(|| s.compute(black_box(term), &config)));
    }
    group.finish();
}

fn numeric(c: &mut Criterion) {
    let s = Session::new();
    s.set_cache_policy(CachePolicy::off());
    let term = s
        .parse("defint(exp(-x^2)*cos(x), x, 0, 2)")
        .expect("parses");
    let config = Config::new().numeric(1e-10);
    c.bench_function("numeric/quadrature", |b| {
        b.iter(|| s.compute(black_box(term), &config))
    });
    let term = s.parse("hann_window(8)").expect("parses");
    c.bench_function("numeric/bridged_window", |b| {
        b.iter(|| s.compute(black_box(term), &config))
    });
}

fn symbolic_batch(c: &mut Criterion) {
    let s = Session::new();
    s.set_cache_policy(CachePolicy::off());
    let config = Config::new();
    let terms: Vec<_> = (1..=32)
        .map(|power| {
            s.parse(&format!("diff(x^{power} * sin(x), x)"))
                .expect("benchmark request parses")
        })
        .collect();
    c.bench_function("symbolic_batch/ordered_loop_32", |b| {
        b.iter(|| {
            black_box(&terms)
                .iter()
                .map(|&term| s.compute(term, &config))
                .collect::<Vec<_>>()
        });
    });
    c.bench_function("symbolic_batch/compute_all_32", |b| {
        b.iter(|| s.compute_all(black_box(&terms), &config));
    });
}

fn compiled(c: &mut Criterion) {
    let s = Session::new();
    let f = s.parse("sin(x)*exp(-x^2/2) + x^3").expect("parses");
    let compiled = f.compile(&["x"]).expect("compiles");
    let xs: Vec<f64> = (0..10_000).map(|i| f64::from(i) * 1e-3).collect();
    let mut out = vec![0.0; xs.len()];
    c.bench_function("compiled/batch_10k", |b| {
        b.iter(|| compiled.call_batch(&[black_box(&xs)], &mut out));
    });
}

fn frobenius_norm(c: &mut Criterion) {
    use rssn::kernels::dense::Mat;

    let mut group = c.benchmark_group("frobenius_norm");
    for width in [8_u32, 64, 256] {
        let entries = (1..=width * width)
            .map(|index| f64::from(index).sin())
            .collect();
        let size = usize::try_from(width).expect("matrix width fits usize");
        let matrix = Mat {
            rows: size,
            cols: size,
            data: entries,
        };
        group.bench_function(BenchmarkId::from_parameter(width), |b| {
            b.iter(|| black_box(&matrix).norm_fro());
        });
    }
    group.finish();
}

fn circular_analysis(c: &mut Criterion) {
    use num_complex::Complex64;
    use rssn::kernels::analytic_continuation::laurent_coeffs_fft;
    use rssn::kernels::complex::cauchy_fft_derivatives;

    let center = Complex64::new(0.25, 0.125);
    let polynomial = |z: Complex64| z * z + 3.0 * z + 5.0;
    let meromorphic = |z: Complex64| 1.0 / (z - center) + polynomial(z);
    let mut group = c.benchmark_group("circular_analysis");
    for count in [32, 257, 1_024] {
        group.bench_function(BenchmarkId::new("cauchy", count), |b| {
            b.iter(|| cauchy_fft_derivatives(polynomial, black_box(center), 0.75, 5, count));
        });
        group.bench_function(BenchmarkId::new("laurent", count), |b| {
            b.iter(|| laurent_coeffs_fft(meromorphic, black_box(center), 0.75, 2, 5, count));
        });
    }
    group.finish();
}

fn polynomials(c: &mut Criterion) {
    use rssn::kernels::interpolate::lagrange_interpolation;
    use rssn::kernels::polynomial::Polynomial;

    let mut group = c.benchmark_group("polynomials");
    for count in [8, 32, 128] {
        // Binary-exact nodes and ordinates keep this quadratic well
        // conditioned at every size and give both input layouts the same data.
        let nodes: Vec<_> = (0..count)
            .map(|i| 2.0 * f64::from(i) / f64::from(count) - 1.0)
            .collect();
        let ordinates: Vec<_> = nodes.iter().map(|&x| x.mul_add(x, 1.0)).collect();
        let pairs: Vec<_> = nodes
            .iter()
            .copied()
            .zip(ordinates.iter().copied())
            .collect();
        group.bench_function(BenchmarkId::new("lagrange_columns", count), |b| {
            b.iter(|| Polynomial::interpolate_lagrange(black_box(&nodes), black_box(&ordinates)));
        });
        group.bench_function(BenchmarkId::new("lagrange_pairs", count), |b| {
            b.iter(|| lagrange_interpolation(black_box(&pairs)));
        });
    }

    let model = Polynomial::new((0..=64).map(|i| f64::from(i % 7 - 3) / 8.0).collect());
    for count in [16, 512, 32_768] {
        let points: Vec<_> = (0..count)
            .map(|i| 2.0 * f64::from(i) / f64::from(count) - 1.0)
            .collect();
        group.bench_function(BenchmarkId::new("scalar_evaluation", count), |b| {
            b.iter(|| {
                black_box(&points)
                    .iter()
                    .map(|&x| model.eval(x))
                    .collect::<Vec<_>>()
            });
        });
        group.bench_function(BenchmarkId::new("batch_evaluation", count), |b| {
            b.iter(|| model.eval_batch(black_box(&points)));
        });
    }
    group.finish();
}

fn fixed_prefix_scan(c: &mut Criterion) {
    use rssn::kernels::par::prefix_sum_f64_fixed;

    let mut group = c.benchmark_group("fixed_prefix_scan");
    for size in [256, 8_192, 131_072] {
        let items: Vec<_> = (0..size).map(|i| f64::from(i % 31 - 15) / 8.0).collect();
        group.bench_function(BenchmarkId::new("single_block", size), |b| {
            b.iter(|| prefix_sum_f64_fixed(black_box(&items), items.len(), 32_768));
        });
        group.bench_function(BenchmarkId::new("blocks_of_256", size), |b| {
            b.iter(|| prefix_sum_f64_fixed(black_box(&items), 256, 32_768));
        });
    }
    group.finish();
}

fn numerical_codec_dispatch(c: &mut Criterion) {
    use rssn::kernel::CallCtx;
    use rssn::kernel::Value;

    let session = Session::new();
    let ctx = CallCtx::new();
    let inputs = [Value::int(13), Value::int(97)];
    c.bench_function("numerical_codec_dispatch/gcd_u128", |b| {
        b.iter(|| {
            session
                .call_value("num_number_theory_gcd_u128", black_box(&inputs), &ctx)
                .expect("valid integer pair")
        });
    });
}

fn triangle_areas(c: &mut Criterion) {
    use rssn::execution::with_scoped_pool;
    use rssn::kernels::meshgen::triangle_areas;
    use std::time::Instant;

    let points = [[0.0, 0.0], [2.0, 0.0], [0.0, 3.0], [1.0, 0.0]];
    let mut group = c.benchmark_group("mesh_triangle_areas");
    for count in [1_024, 8_192, 65_536] {
        let cells: Vec<_> = [[0, 1, 2], [0, 3, 2], [0, 1, 3]]
            .into_iter()
            .cycle()
            .take(count)
            .collect();
        for threads in [1, 4] {
            group.bench_function(BenchmarkId::new(format!("{threads}_threads"), count), |b| {
                b.iter_custom(|iterations| {
                    with_scoped_pool(Some(threads), || {
                        let start = Instant::now();
                        for _ in 0..iterations {
                            black_box(triangle_areas(black_box(&points), black_box(&cells)));
                        }
                        start.elapsed()
                    })
                });
            });
        }
    }
    group.finish();
}

fn special_combinatorics(c: &mut Criterion) {
    use rssn::kernels::special::binomial;
    use rssn::kernels::special::ln_binomial;
    use rssn::kernels::special::multinomial;

    let mut group = c.benchmark_group("special_combinatorics");
    group.bench_function("log_binomial_unbalanced", |bench| {
        bench.iter(|| ln_binomial(black_box(u64::MAX), black_box(3)));
    });
    group.bench_function("log_binomial_balanced", |bench| {
        bench.iter(|| ln_binomial(black_box(1_u64 << 40), black_box(1_u64 << 39)));
    });
    group.bench_function("binomial_overflow", |bench| {
        bench.iter(|| binomial(black_box(2_048), black_box(1_024)));
    });
    group.bench_function("multinomial_partition", |bench| {
        bench.iter(|| multinomial(black_box(120), black_box(&[10, 20, 30, 60])));
    });
    group.finish();
}

fn beta_quantiles(c: &mut Criterion) {
    use rssn::kernels::special::inverse_beta_reg;

    let mut group = c.benchmark_group("beta_quantiles");
    for (name, a, b, probability) in [
        ("interior", 2.5, 7.0, 0.23),
        ("upper_tail", 3.0, 5.0, 0.9999),
        ("unit_shape_tiny_tail", 1.0, 2.0, 2.0_f64.powi(-54)),
    ] {
        group.bench_function(name, |bench| {
            bench.iter(|| inverse_beta_reg(black_box(a), black_box(b), black_box(probability)));
        });
    }
    group.finish();
}

fn zero_phase_filter(c: &mut Criterion) {
    use rssn::kernels::signal::checked_filtfilt;

    let mut group = c.benchmark_group("zero_phase_filter");
    for (taps, sample_count) in [(8_u16, 32_u32), (128, 32), (16, 16_384)] {
        let b = vec![1.0 / f64::from(taps); usize::from(taps)];
        let a = [1.0];
        let samples: Vec<_> = (0..sample_count).map(|i| f64::from(i % 13) - 6.0).collect();
        let id = BenchmarkId::new(format!("{taps}_taps"), sample_count);
        group.bench_function(id, |bench| {
            bench.iter(|| checked_filtfilt(black_box(&b), black_box(&a), black_box(&samples)));
        });
    }
    group.finish();
}

fn sharp_cauchy_radius(c: &mut Criterion) {
    use num_complex::Complex64;
    use rssn::kernels::complex_roots::cauchy_sharp_root_bound;

    let mut group = c.benchmark_group("sharp_cauchy_radius");
    for degree in [16, 128, 1024] {
        let coefficients: Vec<_> = (0..=degree)
            .map(|i| Complex64::new(1.0 + f64::from(i % 7), f64::from(i % 3)))
            .collect();
        group.bench_with_input(
            BenchmarkId::from_parameter(degree),
            &coefficients,
            |bench, coefficients| {
                bench.iter(|| cauchy_sharp_root_bound(black_box(coefficients)));
            },
        );
    }
    group.finish();
}

fn finite_field_ntt(c: &mut Criterion) {
    use criterion::BatchSize;
    use rssn::kernels::finite_field::poly_mul_ntt_mod;
    use rssn::kernels::finite_field::primitive_root_prime;
    use rssn::kernels::reuse::global_pool;

    let prime = 18_446_744_069_414_584_321;
    let mut group = c.benchmark_group("finite_field_ntt");
    group.bench_function("generator_cold", |bench| {
        bench.iter_batched(
            || global_pool().clear_all(),
            |()| primitive_root_prime(black_box(prime)),
            BatchSize::PerIteration,
        );
    });
    let _ = primitive_root_prime(prime);
    group.bench_function("generator_warm", |bench| {
        bench.iter(|| primitive_root_prime(black_box(prime)));
    });
    for size in [256, 8_192] {
        let coefficients: Vec<_> = (0..size).map(|i| prime - 1 - i).collect();
        let _ = poly_mul_ntt_mod(&coefficients, &coefficients, prime);
        group.bench_with_input(
            BenchmarkId::new("convolution_warm", size),
            &coefficients,
            |bench, coefficients| {
                bench.iter(|| {
                    poly_mul_ntt_mod(black_box(coefficients), black_box(coefficients), prime)
                });
            },
        );
    }
    group.finish();
}

fn economy_qr(c: &mut Criterion) {
    use rssn::kernels::dense::Mat;
    use rssn::kernels::dense::qr_thin;

    let mut group = c.benchmark_group("economy_qr");
    for (rows, cols) in [(128, 8), (1024, 16), (4096, 32), (128, 128)] {
        let values: Vec<Vec<f64>> = (0..rows)
            .map(|i| {
                (0..cols)
                    .map(|j| f64::from(i * 31 + j * 3).sin() + f64::from(i == j))
                    .collect()
            })
            .collect();
        let matrix = Mat::from_rows(&values).expect("rectangular benchmark input");
        group.bench_with_input(
            BenchmarkId::from_parameter(format!("{rows}x{cols}")),
            &matrix,
            |bench, matrix| {
                bench.iter(|| {
                    rssn::execution::with_scoped_pool(Some(1), || {
                        qr_thin(black_box(matrix)).expect("finite input has a QR")
                    })
                });
            },
        );
    }
    group.finish();
}

fn kalman_history(c: &mut Criterion) {
    use rssn::kernels::timeseries::kalman_filter;
    use rssn::kernels::timeseries::kalman_smooth;

    let f = [vec![1.0, 1.0], vec![0.0, 1.0]];
    let h = [vec![1.0, 0.0]];
    let q = [vec![0.01, 0.0], vec![0.0, 0.01]];
    let r = [vec![0.5]];
    let x0 = [0.0, 1.0];
    let p0 = [vec![1.0, 0.0], vec![0.0, 1.0]];
    let mut group = c.benchmark_group("kalman_history");
    for count in [16, 256, 4_096] {
        let observations: Vec<_> = (1..=count).map(|t| vec![f64::from(t)]).collect();
        group.bench_with_input(
            BenchmarkId::new("filter", count),
            &observations,
            |bench, y| {
                bench.iter(|| kalman_filter(black_box(y), &f, &h, &q, &r, &x0, &p0));
            },
        );
        group.bench_with_input(
            BenchmarkId::new("smooth", count),
            &observations,
            |bench, y| {
                bench.iter(|| kalman_smooth(black_box(y), &f, &h, &q, &r, &x0, &p0));
            },
        );
    }
    group.finish();
}

fn rolling_window(c: &mut Criterion) {
    use rssn::kernels::timeseries::rolling::rolling_moments;

    let mut group = c.benchmark_group("rolling_moments");
    for count in [65_536, 262_144] {
        let samples: Vec<_> = (0..count).map(|i| 1e12 + f64::from(i % 127)).collect();
        for window in [16, 64, 1_024] {
            group.bench_with_input(
                BenchmarkId::new(format!("width_{window}"), count),
                &samples,
                |bench, samples| {
                    bench.iter(|| rolling_moments(black_box(samples), window, 0));
                },
            );
        }
    }
    group.finish();
}

/// Interpreter against the Cranelift JIT on the same batch, plus the JIT's
/// compile latency (the tiering threshold is chosen from these numbers).
#[cfg(feature = "jit")]
fn jit(c: &mut Criterion) {
    use rssn::backend::jit::CraneliftBackend;
    let s = Session::new();
    let f = s
        .parse("sin(x)*exp(-x^2/2) + x^3 - 2*x/(1 + x^2) + sqrt(x + 1)")
        .expect("parses");
    let xs: Vec<f64> = (0..10_000).map(|i| f64::from(i) * 1e-3).collect();
    let mut out = vec![0.0; xs.len()];
    let interpreted = f
        .compile_with(&rssn::backend::Interpreter, &["x"])
        .expect("compiles");
    c.bench_function("jit/interpreter_batch_10k", |b| {
        b.iter(|| interpreted.call_batch(&[black_box(&xs)], &mut out));
    });
    let backend = CraneliftBackend::new();
    let compiled = f.compile_with(&backend, &["x"]).expect("compiles");
    c.bench_function("jit/cranelift_batch_10k", |b| {
        b.iter(|| compiled.call_batch(&[black_box(&xs)], &mut out));
    });
    c.bench_function("jit/compile_latency", |b| {
        b.iter(|| {
            f.compile_with(&CraneliftBackend::new(), &["x"])
                .expect("compiles")
        });
    });
}

#[cfg(not(feature = "jit"))]
fn jit(_: &mut Criterion) {}

criterion_group!(
    benches,
    symbolic,
    numeric,
    symbolic_batch,
    compiled,
    circular_analysis,
    frobenius_norm,
    polynomials,
    fixed_prefix_scan,
    numerical_codec_dispatch,
    triangle_areas,
    special_combinatorics,
    beta_quantiles,
    zero_phase_filter,
    sharp_cauchy_radius,
    finite_field_ntt,
    economy_qr,
    kalman_history,
    rolling_window,
    jit
);
criterion_main!(benches);
