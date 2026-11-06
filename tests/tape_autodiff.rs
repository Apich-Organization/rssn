//! Reverse AD uses independently seeded forward and hyper-dual runs as oracles.
use rssn::backend::Intrinsic;
use rssn::backend::tape::Cmp;
use rssn::backend::tape::Dual;
use rssn::backend::tape::Inst;
use rssn::backend::tape::Tape;

fn near(
    a: f64,
    b: f64,
) {
    assert!(
        (a - b).abs() <= 2e-10 * (1.0 + a.abs().max(b.abs())),
        "{a:e} != {b:e}"
    );
}

fn check(
    tape: &Tape,
    x: &[f64],
    p: &[f64],
    direction: &[f64],
) {
    let params: Vec<_> = p.iter().map(|&x| Dual(x, 0.0)).collect();
    for k in 0..tape.outputs.len() {
        let gradient = tape.gradient_reverse(k, x, p);
        let mut values = Vec::new();
        let mut output = vec![Dual(0.0, 0.0); tape.outputs.len()];
        for i in 0..x.len() {
            let input: Vec<_> = x
                .iter()
                .enumerate()
                .map(|(j, &x)| Dual(x, if j == i { 1.0 } else { 0.0 }))
                .collect();
            tape.run(&input, &params, &mut values, &mut output);
            near(gradient[i], output[k].1);
        }
        let hessian = tape.hessian(k, x, p);
        let product = tape.hessian_vector_product(k, x, p, direction);
        for i in 0..x.len() {
            let expected = (0..x.len())
                .map(|j| hessian[i * x.len() + j] * direction[j])
                .sum();
            near(product[i], expected);
        }
    }
}

#[test]
fn arithmetic_and_shared_subexpressions_match_independent_forward_ad() {
    let tape = Tape {
        insts: vec![
            Inst::Input(0),
            Inst::Input(1),
            Inst::Input(2),
            Inst::Param(0),
            Inst::Mul(0, 1),
            Inst::Add(4, 4),
            Inst::Div(5, 2),
            Inst::Neg(6),
            Inst::Sub(7, 3),
            Inst::PowI(8, 3),
            Inst::PowF(2, 1),
            Inst::Atan2(9, 10),
        ],
        inputs: 3,
        params: 1,
        outputs: vec![9, 10, 11],
        args: vec![],
    };
    check(&tape, &[0.3, 1.7, 2.1], &[0.2], &[0.7, -0.2, 1.1]);
    let a = tape.gradient_reverse(0, &[0.3, 1.7, 2.1], &[0.2]);
    let b = tape.gradient_reverse(2, &[0.3, 1.7, 2.1], &[0.2]);
    let vjp = tape.vjp(&[0.3, 1.7, 2.1], &[0.2], &[2.0, 0.0, -0.5]);
    for i in 0..3 {
        near(vjp[i], 2.0 * a[i] - 0.5 * b[i]);
    }
}

#[test]
fn every_analytic_unary_rule_matches_hyperdual_curvature() {
    for op in [
        Intrinsic::Exp,
        Intrinsic::Ln,
        Intrinsic::Sin,
        Intrinsic::Cos,
        Intrinsic::Tan,
        Intrinsic::Asin,
        Intrinsic::Acos,
        Intrinsic::Atan,
        Intrinsic::Sinh,
        Intrinsic::Cosh,
        Intrinsic::Tanh,
        Intrinsic::Asinh,
        Intrinsic::Atanh,
        Intrinsic::Sqrt,
        Intrinsic::Abs,
        Intrinsic::Floor,
        Intrinsic::Ceil,
        Intrinsic::Cbrt,
        Intrinsic::Acosh,
    ] {
        let x = if matches!(op, Intrinsic::Acosh) {
            1.4
        } else {
            0.4
        };
        let tape = Tape {
            insts: vec![Inst::Input(0), Inst::Unary(op, 0)],
            inputs: 1,
            outputs: vec![1],
            ..Tape::default()
        };
        check(&tape, &[x], &[], &[1.7]);
    }
}

#[test]
fn piecewise_selection_and_gather_scatter_only_the_active_dependencies() {
    let tape = Tape {
        insts: vec![
            Inst::Input(0),
            Inst::Const(0.0),
            Inst::Cmp(Cmp::Gt, 0, 1),
            Inst::Input(1),
            Inst::Unary(Intrinsic::Ln, 3),
            Inst::PowI(0, 2),
            Inst::Select(2, 5, 4),
            Inst::Const(1.9),
            Inst::Gather(0, 7),
            Inst::Mul(8, 8),
            Inst::Add(6, 9),
        ],
        inputs: 2,
        outputs: vec![10],
        ..Tape::default()
    };
    // Dead ln(-3) must not poison derivatives of the positive branch.
    assert_eq!(tape.gradient_reverse(0, &[2.0, -3.0], &[]), vec![4.0, -6.0]);
    assert_eq!(
        tape.hessian_vector_product(0, &[2.0, -3.0], &[], &[0.5, 2.0]),
        vec![1.0, 4.0]
    );
    check(&tape, &[-2.0, 3.0], &[], &[0.5, 2.0]);
    let repeated = Tape {
        outputs: vec![10, 10],
        ..tape
    };
    assert_eq!(
        repeated.vjp(&[2.0, -3.0], &[], &[2.0, 3.0]),
        vec![20.0, -30.0]
    );
}

#[test]
fn powers_at_zero_and_negative_bases_do_not_use_logarithmic_singularities() {
    for n in [0, 1, 2, 3] {
        let tape = Tape {
            insts: vec![Inst::Input(0), Inst::Const(f64::from(n)), Inst::PowF(0, 1)],
            inputs: 1,
            outputs: vec![2],
            ..Tape::default()
        };
        assert_eq!(
            tape.gradient_reverse(0, &[0.0], &[]),
            vec![if n == 1 { 1.0 } else { 0.0 }]
        );
        assert_eq!(
            tape.hessian_vector_product(0, &[0.0], &[], &[3.0]),
            vec![if n == 2 { 6.0 } else { 0.0 }]
        );
    }
    let tape = Tape {
        insts: vec![Inst::Input(0), Inst::Const(3.0), Inst::PowF(0, 1)],
        inputs: 1,
        outputs: vec![2],
        ..Tape::default()
    };
    assert_eq!(tape.gradient_reverse(0, &[-2.0], &[]), vec![12.0]);
    assert_eq!(
        tape.hessian_vector_product(0, &[-2.0], &[], &[2.0]),
        vec![-24.0]
    );
    let min_power = Tape {
        insts: vec![Inst::Input(0), Inst::PowI(0, i32::MIN)],
        inputs: 1,
        outputs: vec![1],
        ..Tape::default()
    };
    assert_eq!(
        min_power.gradient_reverse(0, &[1.0], &[]),
        vec![f64::from(i32::MIN)]
    );
    near(
        min_power.hessian_vector_product(0, &[1.0], &[], &[1.0])[0],
        f64::from(i32::MIN) * (f64::from(i32::MIN) - 1.0),
    );
}

fn opaque(x: &[f64]) -> f64 {
    x[0].sin() * x[1]
}
extern "C" fn external(x: f64) -> f64 {
    x * x
}
extern "C" fn external_two(
    x: f64,
    y: f64,
) -> f64 {
    x * y
}

#[test]
fn opaque_calls_keep_first_derivatives_but_do_not_invent_second_order_rules() {
    for inst in [
        Inst::Call(opaque, 0, 2),
        Inst::Extern1(external, 0),
        Inst::Extern2(external_two, 0, 1),
    ] {
        let tape = Tape {
            insts: vec![Inst::Input(0), Inst::Input(1), inst],
            args: vec![0, 1],
            inputs: 2,
            outputs: vec![2],
            ..Tape::default()
        };
        let reverse = tape.gradient_reverse(0, &[0.4, 1.2], &[]);
        assert!(reverse.iter().all(|x| x.is_finite()));
        let second = tape.hessian_vector_product(0, &[0.4, 1.2], &[], &[1.0, 0.5]);
        assert!(second[0].is_nan());
        if matches!(inst, Inst::Extern1(..)) {
            assert_eq!(second[1], 0.0);
        } else {
            assert!(second[1].is_nan());
        }
    }
}

#[test]
fn invalid_seeds_and_active_domains_are_reported_without_poisoning_unused_outputs() {
    let tape = Tape {
        insts: vec![
            Inst::Input(0),
            Inst::Input(1),
            Inst::Unary(Intrinsic::Ln, 1),
        ],
        inputs: 2,
        outputs: vec![0, 2],
        ..Tape::default()
    };
    assert_eq!(tape.vjp(&[2.0, -1.0], &[], &[3.0, 0.0]), vec![3.0, 0.0]);
    for gradient in [
        tape.gradient_reverse(5, &[2.0, -1.0], &[]),
        tape.gradient_reverse(1, &[2.0, -1.0], &[]),
        tape.vjp(&[2.0, -1.0], &[], &[1.0]),
        tape.hessian_vector_product(0, &[2.0, -1.0], &[], &[]),
    ] {
        assert!(gradient.iter().all(|x| x.is_nan()));
    }
}

#[test]
fn gradient_size_scales_with_inputs_without_one_tape_run_per_input() {
    let n = 1000;
    let mut tape = Tape {
        inputs: n,
        ..Tape::default()
    };
    for i in 0..n {
        tape.insts.push(Inst::Input(i as u32));
    }
    let mut sum = 0;
    for i in 1..n {
        tape.insts.push(Inst::Add(sum, i as u32));
        sum = (tape.insts.len() - 1) as u32;
    }
    tape.insts.push(Inst::PowI(sum, 2));
    tape.outputs.push((tape.insts.len() - 1) as u32);
    let input = vec![0.1; n];
    for x in tape.gradient_reverse(0, &input, &[]) {
        near(x, 200.0);
    }
    for x in tape.hessian_vector_product(0, &input, &[], &vec![1.0; n]) {
        near(x, 2000.0);
    }
}

static CALLS: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
fn counted_square(x: &[f64]) -> f64 {
    CALLS.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    x[0] * x[0]
}

#[test]
fn opaque_evaluation_count_is_independent_of_input_dimension() {
    let n = 1000;
    let mut tape = Tape {
        inputs: n,
        ..Tape::default()
    };
    for i in 0..n {
        tape.insts.push(Inst::Input(i as u32));
    }
    let mut sum = 0;
    for i in 1..n {
        tape.insts.push(Inst::Add(sum, i as u32));
        sum = (tape.insts.len() - 1) as u32;
    }
    tape.args.push(sum);
    tape.insts.push(Inst::Call(counted_square, 0, 1));
    tape.outputs.push((tape.insts.len() - 1) as u32);
    CALLS.store(0, std::sync::atomic::Ordering::Relaxed);
    let result = tape.gradient_reverse(0, &vec![1.0; n], &[]);
    assert!(CALLS.load(std::sync::atomic::Ordering::Relaxed) <= 4);
    for gradient in result {
        assert!((gradient - 2000.0).abs() < 1e-6);
    }
}

#[test]
fn graph_lowering_and_optimization_preserve_reverse_and_hvp_semantics() {
    use rssn::backend::tape::lower;
    use rssn::graph::Engine;
    use rssn::graph::Graph;
    let mut graph = Graph::new();
    Engine::install(&mut graph, &rssn::rules::standard()).unwrap();
    let node = graph.parse("exp(sin(x*y)) + x^3/(y+2)").unwrap();
    let inputs = [
        graph.interner_mut().symbol("x"),
        graph.interner_mut().symbol("y"),
    ];
    let tape = lower(&graph, &[node], &inputs, &[]).unwrap();
    for program in [&tape, &tape.optimise()] {
        check(program, &[0.3, 1.2], &[], &[0.4, -0.7]);
    }
}
