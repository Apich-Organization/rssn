#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::indexing_slicing,
    clippy::float_cmp,
    clippy::panic
)]

use std::sync::Arc;

use super::wgsl::KernelEmitter;
use super::wgsl::WgslEmitter;
use super::wgsl::WgslOptions;
use super::wgsl::emit_wgsl;
use super::*;
use crate::backend::Intrinsic;
use crate::backend::tape::Cmp;
use crate::backend::tape::Inst;
use crate::backend::tape::Tape;
use crate::graph::Engine;
use crate::rules;

fn same(
    a: f64,
    b: f64,
) -> bool {
    (a.is_nan() && b.is_nan()) || a.to_bits() == b.to_bits()
}

fn tape(
    insts: Vec<Inst>,
    inputs: usize,
    params: usize,
    outputs: Vec<u32>,
) -> Tape {
    Tape {
        insts,
        args: Vec::new(),
        inputs,
        params,
        outputs,
    }
}

/// Every non-host instruction, two outputs.
fn everything() -> KernelIr {
    let insts = vec![
        Inst::Input(0),                   // 0 x
        Inst::Input(1),                   // 1 y
        Inst::Param(0),                   // 2 p
        Inst::Const(0.5),                 // 3
        Inst::Add(0, 1),                  // 4
        Inst::Sub(4, 3),                  // 5
        Inst::Mul(5, 2),                  // 6
        Inst::Div(6, 1),                  // 7
        Inst::Neg(7),                     // 8
        Inst::PowI(0, 3),                 // 9
        Inst::PowI(1, -2),                // 10
        Inst::PowF(0, 1),                 // 11
        Inst::Unary(Intrinsic::Sin, 0),   // 12
        Inst::Unary(Intrinsic::Exp, 1),   // 13
        Inst::Unary(Intrinsic::Sqrt, 12), // 14
        Inst::Unary(Intrinsic::Cbrt, 0),  // 15
        Inst::Atan2(0, 1),                // 16
        Inst::Cmp(Cmp::Lt, 0, 1),         // 17
        Inst::Cmp(Cmp::Ne, 0, 1),         // 18
        Inst::Select(17, 8, 9),           // 19
        Inst::Const(1.0),                 // 20
        Inst::Gather(0, 20),              // 21 inputs[0 + 1] = y
        Inst::Add(19, 21),                // 22
        Inst::Add(22, 10),                // 23
        Inst::Add(23, 11),                // 24
        Inst::Add(24, 13),                // 25
        Inst::Add(25, 14),                // 26
        Inst::Add(26, 15),                // 27
        Inst::Add(27, 16),                // 28
        Inst::Add(28, 18),                // 29
    ];
    KernelIr::from_tape(&tape(insts, 2, 1, vec![29, 19])).unwrap()
}

/// A two-input, two-output elementwise program without `gather` (which
/// addresses columns and therefore cannot be spliced into another kernel).
fn smooth() -> KernelIr {
    let insts = vec![
        Inst::Input(0),                  // 0 x
        Inst::Input(1),                  // 1 y
        Inst::Param(0),                  // 2 p
        Inst::Mul(0, 1),                 // 3
        Inst::Unary(Intrinsic::Sin, 3),  // 4
        Inst::PowI(0, 2),                // 5
        Inst::Add(5, 2),                 // 6
        Inst::Unary(Intrinsic::Sqrt, 6), // 7
        Inst::Cmp(Cmp::Lt, 0, 1),        // 8
        Inst::Select(8, 4, 7),           // 9
        Inst::Sub(7, 4),                 // 10
    ];
    KernelIr::from_tape(&tape(insts, 2, 1, vec![9, 10])).unwrap()
}

fn sample(rows: usize) -> (Vec<f64>, Vec<f64>) {
    let mut s = 0x9E37_79B9_7F4A_7C15u64;
    let mut next = move || {
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        #[allow(clippy::cast_precision_loss)]
        let u = (s >> 11) as f64 / (1u64 << 53) as f64;
        u * 8.0 - 4.0
    };
    let mut xs: Vec<f64> = (0..rows).map(|_| next()).collect();
    let ys: Vec<f64> = (0..rows).map(|_| next()).collect();
    if rows > 4 {
        xs[0] = 0.0;
        xs[1] = f64::NAN;
        xs[2] = -2.0;
        xs[3] = f64::INFINITY;
    }
    (xs, ys)
}

#[test]
fn from_tape_rejects_host_calls_and_bad_programs() {
    extern "C" fn id(x: f64) -> f64 {
        x
    }
    let t = tape(vec![Inst::Input(0), Inst::Extern1(id, 0)], 1, 0, vec![1]);
    assert!(matches!(
        KernelIr::from_tape(&t),
        Err(DeviceError::Unsupported(_))
    ));
    let forward = tape(vec![Inst::Neg(1), Inst::Const(1.0)], 0, 0, vec![0]);
    assert!(matches!(
        KernelIr::from_tape(&forward),
        Err(DeviceError::InvalidKernel(_))
    ));
    let oob = tape(vec![Inst::Input(3)], 1, 0, vec![0]);
    assert!(KernelIr::from_tape(&oob).is_err());
    let no_out = tape(vec![Inst::Const(1.0)], 0, 0, vec![]);
    assert!(KernelIr::from_tape(&no_out).is_err());
}

#[test]
fn kernel_matches_the_tape_interpreter() {
    let ir = everything();
    let tape = ir.to_tape().unwrap();
    let (xs, ys) = sample(64);
    for r in 0..64 {
        let mut want = [0.0; 2];
        tape.eval(&[xs[r], ys[r]], &[1.25], &mut want);
        let got = ir.eval_row(&[xs[r], ys[r]], &[1.25]);
        assert!(same(got[0], want[0]) && same(got[1], want[1]), "row {r}");
    }
}

#[test]
fn cpu_devices_are_bit_identical_to_the_reference() {
    let ir = everything();
    // Spans several lane blocks and, for cpu-parallel, several splits.
    let rows = 2 * cpu::GRAIN + 777;
    let (xs, ys) = sample(rows);
    let seq = CpuDevice::sequential();
    let par = CpuDevice::parallel();
    let a = map_columns(&seq, &ir, &[&xs, &ys], &[1.25]).unwrap();
    let b = map_columns(&par, &ir, &[&xs, &ys], &[1.25]).unwrap();
    assert_eq!(a.len(), 2);
    for k in 0..2 {
        assert_eq!(a[k].len(), rows);
        for r in 0..rows {
            assert!(
                same(a[k][r], b[k][r]),
                "seq/par differ at output {k} row {r}"
            );
        }
    }
    for r in (0..rows).step_by(37) {
        let want = ir.eval_row(&[xs[r], ys[r]], &[1.25]);
        for k in 0..2 {
            assert!(
                same(a[k][r], want[k]),
                "output {k} row {r}: {} vs {}",
                a[k][r],
                want[k]
            );
        }
    }
}

#[test]
fn launch_checks_sizes_and_ownership() {
    let ir = everything();
    let d = CpuDevice::sequential();
    let other = CpuDevice::sequential();
    let k = d.compile(&ir).unwrap();
    let input = d.upload(&[0.0; 6], Dtype::F64).unwrap();
    let mut out = d.alloc(6, Dtype::F64).unwrap();
    // 2 inputs x 3 rows in, 2 outputs x 3 rows out: fine.
    d.launch(&*k, 3, &input, &[1.0], &mut out).unwrap();
    // Wrong row count.
    assert!(matches!(
        d.launch(&*k, 4, &input, &[1.0], &mut out),
        Err(DeviceError::Mismatch(_))
    ));
    // Buffers of another device are refused, never reinterpreted.
    let foreign = other.upload(&[0.0; 6], Dtype::F64).unwrap();
    assert!(matches!(
        d.launch(&*k, 3, &foreign, &[1.0], &mut out),
        Err(DeviceError::Mismatch(_))
    ));
    let mut host = [0.0; 6];
    assert!(d.download(&foreign, &mut host).is_err());
    assert!(d.download(&out, &mut [0.0; 5]).is_err());
    assert!(matches!(
        d.alloc(1, Dtype::F32),
        Err(DeviceError::Unsupported(_))
    ));
    // Missing parameter reads NaN rather than out of bounds.
    d.launch(&*k, 3, &input, &[], &mut out).unwrap();
}

#[test]
fn zero_rows_is_a_no_op() {
    let ir = everything();
    let d = CpuDevice::parallel();
    let r = map_columns(&d, &ir, &[&[], &[]], &[1.0]).unwrap();
    assert_eq!(r, vec![Vec::<f64>::new(), Vec::new()]);
}

#[test]
fn reductions_are_deterministic_across_devices() {
    let n = 3 * cpu::REDUCE_CHUNK + 123;
    let (xs, _) = sample(n);
    let xs: Vec<f64> = xs.into_iter().filter(|x| x.is_finite()).collect();
    let seq = CpuDevice::sequential();
    let par = CpuDevice::parallel();
    let bs = seq.upload(&xs, Dtype::F64).unwrap();
    let bp = par.upload(&xs, Dtype::F64).unwrap();
    for op in [ReduceOp::Sum, ReduceOp::Prod, ReduceOp::Min, ReduceOp::Max] {
        let a = seq.reduce(op, &bs).unwrap();
        let b = par.reduce(op, &bp).unwrap();
        assert!(same(a, b), "{op:?}");
    }
    let sum = seq.reduce(ReduceOp::Sum, &bs).unwrap();
    let naive: f64 = xs.iter().sum();
    assert!((sum - naive).abs() <= 1e-9 * naive.abs().max(1.0));
    assert_eq!(
        seq.reduce(ReduceOp::Min, &bs).unwrap(),
        xs.iter().copied().fold(f64::INFINITY, f64::min)
    );
    let empty = seq.upload(&[], Dtype::F64).unwrap();
    assert_eq!(seq.reduce(ReduceOp::Sum, &empty).unwrap(), 0.0);
    let with_nan = seq.upload(&[1.0, f64::NAN, 2.0], Dtype::F64).unwrap();
    assert!(seq.reduce(ReduceOp::Max, &with_nan).unwrap().is_nan());
}

#[test]
fn matmul_matches_the_naive_product() {
    let (m, k, n) = (37, 19, 23);
    let a: Vec<f64> = (0..m * k).map(|i| f64::from(i as u32 % 7) - 3.0).collect();
    let b: Vec<f64> = (0..k * n).map(|i| f64::from(i as u32 % 5) * 0.5).collect();
    let mut want = vec![0.0; m * n];
    for i in 0..m {
        for j in 0..n {
            for p in 0..k {
                want[i * n + j] += a[i * k + p] * b[p * n + j];
            }
        }
    }
    for dev in [CpuDevice::sequential(), CpuDevice::parallel()] {
        let (ba, bb) = (
            dev.upload(&a, Dtype::F64).unwrap(),
            dev.upload(&b, Dtype::F64).unwrap(),
        );
        let mut bc = dev.alloc(m * n, Dtype::F64).unwrap();
        dev.matmul(&ba, &bb, &mut bc, m, k, n).unwrap();
        let mut got = vec![0.0; m * n];
        dev.download(&bc, &mut got).unwrap();
        assert_eq!(got, want);
        let mut bad = dev.alloc(3, Dtype::F64).unwrap();
        assert!(dev.matmul(&ba, &bb, &mut bad, m, k, n).is_err());
    }
}

/// `dy/dt = -k y`, parameter `k`, state `y`.
fn decay() -> KernelIr {
    let insts = vec![
        Inst::Input(0),
        Inst::Param(0),
        Inst::Mul(0, 1),
        Inst::Neg(2),
    ];
    KernelIr::from_tape(&tape(insts, 1, 1, vec![3])).unwrap()
}

/// `y'' = -y` as `(y, v)' = (v, -y)`.
fn oscillator() -> KernelIr {
    let insts = vec![Inst::Input(0), Inst::Input(1), Inst::Neg(0)];
    KernelIr::from_tape(&tape(insts, 2, 0, vec![1, 2])).unwrap()
}

#[test]
fn rk4_ensemble_integrates_many_initial_conditions() {
    let y0: Vec<f64> = (0..5000).map(|i| 0.5 + f64::from(i) * 1e-3).collect();
    for dev in [CpuDevice::sequential(), CpuDevice::parallel()] {
        let r = ensemble_rk4(&dev, &decay(), &[2.0], &[&y0], 0.0, 1.0, 200).unwrap();
        let e = (-2.0f64).exp();
        for (got, y) in r[0].iter().zip(&y0) {
            assert!((got - y * e).abs() < 1e-9, "{got} vs {}", y * e);
        }
    }
}

#[test]
fn rk4_is_fourth_order_and_rows_are_independent() {
    let dev = CpuDevice::sequential();
    let err = |steps: usize| {
        let r = ensemble_rk4(&dev, &decay(), &[1.0], &[&[1.0]], 0.0, 1.0, steps).unwrap();
        (r[0][0] - (-1.0f64).exp()).abs()
    };
    let (e1, e2) = (err(10), err(20));
    let order = (e1 / e2).log2();
    assert!((3.7..4.3).contains(&order), "observed order {order}");
    // A row's result does not depend on its neighbours.
    let alone = ensemble_rk4(&dev, &decay(), &[1.0], &[&[0.7]], 0.0, 1.0, 20).unwrap();
    let batch = ensemble_rk4(&dev, &decay(), &[1.0], &[&[0.1, 0.7, 5.0]], 0.0, 1.0, 20).unwrap();
    assert!(same(alone[0][0], batch[0][1]));
}

#[test]
fn rk4_system_and_time_dependent_rhs() {
    let dev = CpuDevice::parallel();
    let y = [0.0, 1.0, 2.0];
    let v = [1.0, 0.0, -1.0];
    let r = ensemble_rk4(&dev, &oscillator(), &[], &[&y, &v], 0.0, 3.0, 300).unwrap();
    for i in 0..3 {
        let (s, c) = 3.0f64.sin_cos();
        let want_y = y[i] * c + v[i] * s;
        let want_v = v[i] * c - y[i] * s;
        assert!((r[0][i] - want_y).abs() < 1e-9);
        assert!((r[1][i] - want_v).abs() < 1e-9);
    }
    // y' = t  =>  y(2) = y0 + 2 (input 0 is time, input 1 the state).
    let timed = KernelIr::from_tape(&tape(vec![Inst::Input(0)], 2, 0, vec![0])).unwrap();
    let r = ensemble_rk4(&dev, &timed, &[], &[&[1.0, 2.0]], 0.0, 2.0, 7).unwrap();
    assert!((r[0][0] - 3.0).abs() < 1e-12 && (r[0][1] - 4.0).abs() < 1e-12);
}

#[test]
fn rk4_rejects_bad_shapes() {
    let dev = CpuDevice::sequential();
    assert!(ensemble_rk4(&dev, &decay(), &[], &[&[1.0]], 0.0, 1.0, 4).is_err());
    assert!(ensemble_rk4(&dev, &decay(), &[1.0], &[&[1.0]], 0.0, 1.0, 0).is_err());
    assert!(ensemble_rk4(&dev, &decay(), &[1.0], &[&[1.0], &[1.0]], 0.0, 1.0, 4).is_err());
    assert!(rk4_step(&everything(), 2).is_err());
}

#[test]
fn registry_lists_and_resolves_devices() {
    let reg = DeviceRegistry::new();
    let ids: Vec<String> = reg.list().into_iter().map(|i| i.id).collect();
    assert!(ids.contains(&"cpu".to_string()) && ids.contains(&"cpu-parallel".to_string()));
    assert!(reg.get("cpu").is_some() && reg.get("nope").is_none());
    assert!(reg.first_of(&DeviceKind::CpuParallel).is_some());
    assert!(reg.default_device().is_some());
    // Duplicate ids are refused.
    assert!(reg.register(Arc::new(CpuDevice::sequential())).is_err());
    assert!(DeviceRegistry::empty().default_device().is_none());
    let info = reg.get("cpu-parallel").unwrap().info().clone();
    assert!(info.ops.elementwise && info.ops.reduce && info.ops.matmul);
    assert_eq!(info.dtypes, vec![Dtype::F64]);
    let json = serde_json::to_string(&info).unwrap();
    assert!(json.contains("CpuParallel"));
    let build = build_devices();
    assert!(build.contains(&"cpu".to_string()) && build.contains(&"cpu-parallel".to_string()));
    assert_eq!(
        build.contains(&"gpu-wgpu".to_string()),
        cfg!(feature = "gpu-wgpu")
    );
    assert!(registry().get("cpu").is_some());
}

#[test]
fn a_custom_device_registers_through_the_trait() {
    /// Delegates to a CPU device under a different identity.
    struct Mirror {
        inner: CpuDevice,
        info: DeviceInfo,
    }
    impl Device for Mirror {
        fn info(&self) -> &DeviceInfo {
            &self.info
        }

        fn uid(&self) -> u64 {
            self.inner.uid()
        }

        fn alloc(
            &self,
            len: usize,
            d: Dtype,
        ) -> Result<Buffer, DeviceError> {
            self.inner.alloc(len, d)
        }

        fn upload(
            &self,
            h: &[f64],
            d: Dtype,
        ) -> Result<Buffer, DeviceError> {
            self.inner.upload(h, d)
        }

        fn download(
            &self,
            b: &Buffer,
            o: &mut [f64],
        ) -> Result<(), DeviceError> {
            self.inner.download(b, o)
        }

        fn compile(
            &self,
            ir: &KernelIr,
        ) -> Result<Arc<dyn CompiledKernel>, DeviceError> {
            self.inner.compile(ir)
        }

        fn launch(
            &self,
            k: &dyn CompiledKernel,
            rows: usize,
            i: &Buffer,
            p: &[f64],
            o: &mut Buffer,
        ) -> Result<(), DeviceError> {
            self.inner.launch(k, rows, i, p, o)
        }
    }
    let mut info = CpuDevice::sequential().info().clone();
    info.id = "custom:mirror".into();
    info.kind = DeviceKind::Custom("mirror".into());
    info.ops.reduce = false;
    let reg = DeviceRegistry::empty();
    reg.register(Arc::new(Mirror {
        inner: CpuDevice::sequential(),
        info,
    }))
    .unwrap();
    let dev = reg.get("custom:mirror").unwrap();
    // Unimplemented optional operations report Unsupported.
    let b = dev.upload(&[1.0], Dtype::F64).unwrap();
    assert!(matches!(
        dev.reduce(ReduceOp::Sum, &b),
        Err(DeviceError::Unsupported(_))
    ));
    let r = ensemble_rk4(&*dev, &decay(), &[1.0], &[&[1.0]], 0.0, 1.0, 50).unwrap();
    assert!((r[0][0] - (-1.0f64).exp()).abs() < 1e-8);
    assert_eq!(DeviceKind::Custom("x".into()).tag(), "custom-x");
}

#[test]
fn device_backend_compiles_terms_and_batches() {
    let mut g = crate::graph::Graph::new();
    assert!(Engine::install(&mut g, &rules::standard()).is_ok());
    let root = g.parse("sin(x)^2 + y/(1 + x^2)").unwrap();
    let inputs = vec![g.interner_mut().symbol("x"), g.interner_mut().symbol("y")];
    let reference = crate::backend::Interpreter
        .compile(&g, root, &inputs)
        .unwrap();
    for dev in [
        Arc::new(CpuDevice::sequential()) as Arc<dyn Device>,
        Arc::new(CpuDevice::parallel()),
    ] {
        let backend = DeviceBackend::new(dev);
        let f = crate::backend::Backend::compile(&backend, &g, root, &inputs).unwrap();
        assert_eq!(f.arity(), 2);
        let (xs, ys) = sample(5000);
        let mut out = vec![0.0; 5000];
        f.call_batch(&[&xs, &ys], &mut out);
        for r in (0..5000).step_by(11) {
            let want = reference.call(&[xs[r], ys[r]]);
            assert!(
                same(out[r], want) || (out[r] - want).abs() <= 1e-12 * want.abs().max(1.0),
                "row {r}: {} vs {want}",
                out[r]
            );
        }
        assert!(
            same(f.call(&[0.5, 2.0]), reference.call(&[0.5, 2.0])) || {
                (f.call(&[0.5, 2.0]) - reference.call(&[0.5, 2.0])).abs() < 1e-12
            }
        );
    }
}

// --- structured kernels ---------------------------------------------------------

fn cpus() -> [CpuDevice; 2] {
    [CpuDevice::sequential(), CpuDevice::parallel()]
}

#[test]
fn structured_programs_are_typed_and_phased() {
    // An index cannot be added to a number without a conversion.
    let mut b = KernelBuilder::new(1, 0);
    let row = b.push(KernelOp::Row);
    let x = b.push(KernelOp::Input(0));
    let mixed = b.push(KernelOp::Add(row, x));
    assert!(matches!(
        b.finish(vec![mixed]),
        Err(DeviceError::InvalidKernel(_))
    ));
    // An index is not an output.
    let mut b = KernelBuilder::new(0, 0);
    let row = b.push(KernelOp::Row);
    assert!(b.finish(vec![row]).is_err());
    // The inner index needs a fold, and a body value is not an output.
    let mut b = KernelBuilder::new(0, 0);
    let j = b.push(KernelOp::Inner);
    let jn = b.push(KernelOp::ToNum(j));
    assert!(b.finish(vec![jn]).is_err());
    let mut b = KernelBuilder::new(0, 0);
    let j = b.push(KernelOp::Inner);
    let jn = b.push(KernelOp::ToNum(j));
    let count = b.idx(4);
    let fold = Fold {
        count,
        accs: vec![(ReduceOp::Sum, jn)],
    };
    assert!(b.finish_fold(vec![jn], Some(fold)).is_err());
    // A value cannot depend on both the inner index and an accumulator.
    let mut b = KernelBuilder::new(0, 0);
    let j = b.push(KernelOp::Inner);
    let jn = b.push(KernelOp::ToNum(j));
    let acc = b.push(KernelOp::Acc(0));
    let both = b.push(KernelOp::Add(jn, acc));
    let count = b.idx(4);
    let fold = Fold {
        count,
        accs: vec![(ReduceOp::Sum, both)],
    };
    assert!(b.finish_fold(vec![acc], Some(fold)).is_err());
    // sum_{j<5} j, then doubled after the loop: pre, body and post phases.
    let mut b = KernelBuilder::new(0, 0);
    let j = b.push(KernelOp::Inner);
    let jn = b.push(KernelOp::ToNum(j));
    let acc = b.push(KernelOp::Acc(0));
    let two = b.num(2.0);
    let out = b.push(KernelOp::Mul(acc, two));
    let count = b.idx(5);
    let ir = b
        .finish_fold(
            vec![out],
            Some(Fold {
                count,
                accs: vec![(ReduceOp::Sum, jn)],
            }),
        )
        .unwrap();
    let analysis = ir.analyse().unwrap();
    assert_eq!(analysis.body.len(), 2);
    assert_eq!(analysis.post.len(), 2);
    assert!(!ir.is_elementwise() && ir.to_tape().is_none());
    assert_eq!(ir.eval_row(&[], &[]), vec![20.0]);
    for dev in cpus() {
        let got = array::run_cols(&dev, &ir, 3, &[], &[]).unwrap();
        assert_eq!(got, vec![vec![20.0; 3]]);
    }
}

#[test]
fn a_device_without_structured_kernels_says_so() {
    /// Implements only the elementwise entry point.
    struct Plain(CpuDevice, DeviceInfo);
    impl Device for Plain {
        fn info(&self) -> &DeviceInfo {
            &self.1
        }

        fn uid(&self) -> u64 {
            self.0.uid()
        }

        fn alloc(
            &self,
            len: usize,
            d: Dtype,
        ) -> Result<Buffer, DeviceError> {
            self.0.alloc(len, d)
        }

        fn upload(
            &self,
            h: &[f64],
            d: Dtype,
        ) -> Result<Buffer, DeviceError> {
            self.0.upload(h, d)
        }

        fn download(
            &self,
            b: &Buffer,
            o: &mut [f64],
        ) -> Result<(), DeviceError> {
            self.0.download(b, o)
        }

        fn compile(
            &self,
            ir: &KernelIr,
        ) -> Result<Arc<dyn CompiledKernel>, DeviceError> {
            self.0.compile(ir)
        }

        fn launch(
            &self,
            k: &dyn CompiledKernel,
            rows: usize,
            i: &Buffer,
            p: &[f64],
            o: &mut Buffer,
        ) -> Result<(), DeviceError> {
            self.0.launch(k, rows, i, p, o)
        }
    }
    let inner = CpuDevice::sequential();
    let mut info = inner.info().clone();
    info.id = "plain".into();
    info.ops.structured = false;
    let dev = Plain(inner, info);
    // The elementwise case of `launch_cols` is forwarded to `launch`...
    let got = array::run_cols(&dev, &decay(), 2, &[&[1.0, 2.0]], &[3.0]).unwrap();
    assert_eq!(got, vec![vec![-3.0, -6.0]]);
    // ...anything else is refused, not emulated.
    let b = dev.upload(&[1.0, 2.0, 3.0], Dtype::F64).unwrap();
    assert!(matches!(
        array::reduce(&dev, ReduceOp::Sum, &b),
        Err(DeviceError::Unsupported(_))
    ));
    // The default `copy` goes through the host and is still correct.
    let mut dst = dev.alloc(4, Dtype::F64).unwrap();
    dev.copy(&b, 1, &mut dst, 2, 2).unwrap();
    let mut host = [9.0; 4];
    dev.download(&dst, &mut host).unwrap();
    assert_eq!(host, [0.0, 0.0, 2.0, 3.0]);
    assert!(dev.copy(&b, 2, &mut dst, 0, 2).is_err());
}

#[test]
fn tree_reduction_is_exact_on_integers_and_identical_across_devices() {
    for n in [0_usize, 1, 2, 255, 256, 257, 65_536, 70_001] {
        #[allow(clippy::cast_precision_loss)]
        let data: Vec<f64> = (0..n)
            .map(|i| ((i * 7 + 3) % 1000) as f64 - 500.0)
            .collect();
        for op in [ReduceOp::Sum, ReduceOp::Min, ReduceOp::Max] {
            let want = data.iter().fold(op.identity(), |a, &x| op.combine(a, x));
            let mut results = Vec::new();
            for dev in cpus() {
                let b = dev.upload(&data, Dtype::F64).unwrap();
                results.push(array::reduce(&dev, op, &b).unwrap());
            }
            // Small integers: every order of summation is exact.
            assert!(
                same(results[0], want),
                "{op:?} of {n}: {} vs {want}",
                results[0]
            );
            assert!(same(results[0], results[1]), "{op:?} of {n}");
        }
    }
    // Non-associative data: the two CPU devices still agree bit for bit.
    let (xs, _) = sample(10_000);
    let finite: Vec<f64> = xs
        .iter()
        .map(|x| if x.is_finite() { *x } else { 0.25 })
        .collect();
    let [a, b] = cpus();
    let ba = a.upload(&finite, Dtype::F64).unwrap();
    let bb = b.upload(&finite, Dtype::F64).unwrap();
    let ra = array::reduce(&a, ReduceOp::Sum, &ba).unwrap();
    let rb = array::reduce(&b, ReduceOp::Sum, &bb).unwrap();
    assert!(same(ra, rb));
    assert!((ra - finite.iter().sum::<f64>()).abs() < 1e-9);
    // NaN poisons a minimum, as on the host.
    let poisoned = a.upload(&[1.0, f64::NAN, -3.0], Dtype::F64).unwrap();
    assert!(
        array::reduce(&a, ReduceOp::Min, &poisoned)
            .unwrap()
            .is_nan()
    );
    // One column of a wider buffer, without copying it out.
    let wide = a
        .upload(&[1.0, 2.0, 3.0, 10.0, 20.0, 30.0], Dtype::F64)
        .unwrap();
    assert_eq!(
        array::reduce_column(&a, ReduceOp::Sum, &wide, &[3, 3], 1).unwrap(),
        60.0
    );
    assert!(array::reduce_column(&a, ReduceOp::Sum, &wide, &[3, 3], 2).is_err());
}

#[test]
fn a_fused_map_reduce_equals_map_then_reduce() {
    let ir = smooth();
    let (xs, ys) = sample(5000);
    // Keep the rows on which every output is finite.
    let keep: Vec<usize> = (0..xs.len())
        .filter(|&r| {
            ir.eval_row(&[xs[r], ys[r]], &[1.25])
                .iter()
                .all(|v| v.is_finite())
        })
        .collect();
    let xs: Vec<f64> = keep.iter().map(|&r| xs[r]).collect();
    let ys: Vec<f64> = keep.iter().map(|&r| ys[r]).collect();
    assert!(xs.len() > 1000);
    for dev in cpus() {
        let fused = array::map_reduce(&dev, &ir, &[&xs, &ys], &[1.25], ReduceOp::Sum).unwrap();
        let mapped = map_columns(&dev, &ir, &[&xs, &ys], &[1.25]).unwrap();
        assert_eq!(fused.len(), 2);
        for (k, column) in mapped.iter().enumerate() {
            let b = dev.upload(column, Dtype::F64).unwrap();
            let want = array::reduce(&dev, ReduceOp::Sum, &b).unwrap();
            assert!(same(fused[k], want), "output {k}: {} vs {want}", fused[k]);
        }
        let max = array::map_reduce(&dev, &ir, &[&xs, &ys], &[1.25], ReduceOp::Max).unwrap();
        let want = mapped[1].iter().copied().fold(f64::NEG_INFINITY, f64::max);
        assert!(same(max[1], want));
    }
    let dev = CpuDevice::sequential();
    assert_eq!(
        array::map_reduce(&dev, &decay(), &[&[]], &[1.0], ReduceOp::Prod).unwrap(),
        vec![1.0]
    );
    assert!(array::map_reduce(&dev, &decay(), &[], &[1.0], ReduceOp::Sum).is_err());
    // A map that uses `gather` cannot be fused, and says so.
    assert!(matches!(
        array::map_reduce_pass(&everything(), ReduceOp::Sum),
        Err(DeviceError::InvalidKernel(_))
    ));
}

#[test]
fn matmul_as_a_kernel_matches_the_naive_product() {
    let (m, k, n) = (7, 5, 9);
    #[allow(clippy::cast_precision_loss)]
    let a: Vec<f64> = (0..m * k).map(|i| (i as f64 * 0.37).sin()).collect();
    #[allow(clippy::cast_precision_loss)]
    let b: Vec<f64> = (0..k * n).map(|i| (i as f64 * 0.11).cos()).collect();
    let mut want = vec![0.0; m * n];
    for i in 0..m {
        for j in 0..n {
            // The kernel adds in the order p = 0, 1, ...
            let mut s = 0.0;
            for p in 0..k {
                s += a[i * k + p] * b[p * n + j];
            }
            want[i * n + j] = s;
        }
    }
    for dev in cpus() {
        let ba = dev.upload(&a, Dtype::F64).unwrap();
        let bb = dev.upload(&b, Dtype::F64).unwrap();
        let mut bc = dev.alloc(m * n, Dtype::F64).unwrap();
        array::matmul(&dev, &ba, &bb, &mut bc, m, k, n).unwrap();
        let mut got = vec![0.0; m * n];
        dev.download(&bc, &mut got).unwrap();
        assert!(got.iter().zip(&want).all(|(g, w)| same(*g, *w)));
        assert!(array::matmul(&dev, &ba, &bb, &mut bc, m, k + 1, n).is_err());
        // An empty inner dimension is the zero matrix.
        let e1 = dev.alloc(0, Dtype::F64).unwrap();
        let e2 = dev.alloc(0, Dtype::F64).unwrap();
        let mut z = dev.upload(&[5.0; 6], Dtype::F64).unwrap();
        array::matmul(&dev, &e1, &e2, &mut z, 2, 0, 3).unwrap();
        let mut host = [1.0; 6];
        dev.download(&z, &mut host).unwrap();
        assert_eq!(host, [0.0; 6]);
    }
}

#[test]
fn csr_times_vector_is_a_kernel() {
    // [[2, 0, 1], [0, 0, 0], [0, 3, 4], [5, 0, 0]] * [1, 10, 100]
    let ptr = [0.0, 2.0, 2.0, 4.0, 5.0];
    let cols = [0.0, 2.0, 1.0, 2.0, 0.0];
    let vals = [2.0, 1.0, 3.0, 4.0, 5.0];
    let x = [1.0, 10.0, 100.0];
    let ir = array::spmv_csr_kernel().unwrap();
    for dev in cpus() {
        let y = array::run_cols(&dev, &ir, 4, &[&ptr, &cols, &vals, &x], &[]).unwrap();
        assert_eq!(y, vec![vec![102.0, 0.0, 430.0, 5.0]]);
    }
    // A column index outside the vector reads NaN instead of other memory.
    let bad_cols = [0.0, 7.0, 1.0, 2.0, 0.0];
    let dev = CpuDevice::sequential();
    let y = array::run_cols(&dev, &ir, 4, &[&ptr, &bad_cols, &vals, &x], &[]).unwrap();
    assert!(y[0][0].is_nan() && y[0][2] == 430.0);
}

#[test]
fn stencils_honour_their_boundaries() {
    use array::Boundary;
    let u = [1.0, 4.0, 9.0, 16.0, 25.0];
    let second = [(-1, 1.0), (0, -2.0), (1, 1.0)];
    let dev = CpuDevice::parallel();
    let run = |boundary| {
        let ir = array::stencil1d_kernel(&second, boundary).unwrap();
        array::run_cols(&dev, &ir, u.len(), &[&u], &[])
            .unwrap()
            .remove(0)
    };
    assert_eq!(run(Boundary::Zero), vec![2.0, 2.0, 2.0, 2.0, -34.0]);
    assert_eq!(run(Boundary::Clamp), vec![3.0, 2.0, 2.0, 2.0, -9.0]);
    assert_eq!(run(Boundary::Periodic), vec![27.0, 2.0, 2.0, 2.0, -33.0]);
    assert!(array::stencil1d_kernel(&[], Boundary::Zero).is_err());

    // Five-point Laplacian of u(x, y) = x^2 + y^2 is 4 in the interior.
    let (nx, ny) = (6_usize, 5_usize);
    #[allow(clippy::cast_precision_loss)]
    let grid: Vec<f64> = (0..nx * ny)
        .map(|i| {
            let (x, y) = ((i % nx) as f64, (i / nx) as f64);
            x * x + y * y
        })
        .collect();
    let five = [
        (0, 0, -4.0),
        (1, 0, 1.0),
        (-1, 0, 1.0),
        (0, 1, 1.0),
        (0, -1, 1.0),
    ];
    let ir = array::stencil2d_kernel(&five, Boundary::Clamp).unwrap();
    #[allow(clippy::cast_precision_loss)]
    let lap = array::run_cols(&dev, &ir, nx * ny, &[&grid], &[nx as f64])
        .unwrap()
        .remove(0);
    for y in 1..ny - 1 {
        for x in 1..nx - 1 {
            assert_eq!(lap[y * nx + x], 4.0, "({x}, {y})");
        }
    }
    // Periodic wrap-around stays on its own grid line.
    let shift = array::stencil2d_kernel(&[(1, 0, 1.0)], Boundary::Periodic).unwrap();
    #[allow(clippy::cast_precision_loss)]
    let moved = array::run_cols(&dev, &shift, nx * ny, &[&grid], &[nx as f64])
        .unwrap()
        .remove(0);
    assert_eq!(moved[nx - 1], grid[0]);
    assert_eq!(moved[2 * nx - 1], grid[nx]);
}

#[test]
fn all_pairs_interaction_matches_the_double_loop() {
    // 1-D "gravity": force on i from j is m_j * (x_j - x_i) / (|d|^3 + eps).
    // Pair inputs: x_i, m_i, x_j, m_j.
    let t = tape(
        vec![
            Inst::Input(0),                 // 0 x_i
            Inst::Input(2),                 // 1 x_j
            Inst::Input(3),                 // 2 m_j
            Inst::Sub(1, 0),                // 3 d
            Inst::Unary(Intrinsic::Abs, 3), // 4
            Inst::PowI(4, 3),               // 5
            Inst::Param(0),                 // 6 eps
            Inst::Add(5, 6),                // 7
            Inst::Mul(2, 3),                // 8
            Inst::Div(8, 7),                // 9 force
            Inst::Const(1.0),               // 10 a count
        ],
        4,
        1,
        vec![9, 10],
    );
    let pair = KernelIr::from_tape(&t).unwrap();
    let ir = array::pairwise_kernel(&pair, 2, ReduceOp::Sum, true).unwrap();
    let analysis = ir.analyse().unwrap();
    // Body i's attributes are loaded before the loop, body j's inside it.
    let loads_in = |regs: &[usize]| {
        regs.iter()
            .filter(|&&d| matches!(ir.ops[d], KernelOp::Load(..)))
            .count()
    };
    assert_eq!(loads_in(&analysis.pre), 2);
    assert_eq!(loads_in(&analysis.body), 2);
    let n = 300_usize;
    #[allow(clippy::cast_precision_loss)]
    let xs: Vec<f64> = (0..n).map(|i| (i as f64 * 0.731).sin() * 10.0).collect();
    #[allow(clippy::cast_precision_loss)]
    let ms: Vec<f64> = (0..n).map(|i| 1.0 + (i % 5) as f64).collect();
    let eps = 1e-3;
    for dev in cpus() {
        let out = array::run_cols(&dev, &ir, n, &[&xs, &ms], &[eps]).unwrap();
        for i in (0..n).step_by(37) {
            let mut f = 0.0;
            for j in 0..n {
                if j != i {
                    let d = xs[j] - xs[i];
                    f += ms[j] * d / (d.abs().powi(3) + eps);
                }
            }
            assert!(
                (out[0][i] - f).abs() <= 1e-12 * f.abs().max(1.0),
                "body {i}: {} vs {f}",
                out[0][i]
            );
            #[allow(clippy::cast_precision_loss)]
            let partners = (n - 1) as f64;
            assert_eq!(out[1][i], partners);
        }
    }
    assert!(array::pairwise_kernel(&pair, 3, ReduceOp::Sum, true).is_err());
}

#[test]
fn adaptive_ensemble_gives_each_row_its_own_step() {
    // y' = -k y with the rate k carried as a second state (k' = 0), so one
    // launch integrates gentle and steep rows together.
    let t = tape(
        vec![
            Inst::Input(0),
            Inst::Input(1),
            Inst::Mul(1, 0),
            Inst::Neg(2),
            Inst::Const(0.0),
        ],
        2,
        0,
        vec![3, 4],
    );
    let rhs = KernelIr::from_tape(&t).unwrap();
    let rates: Vec<f64> = vec![0.1, 1.0, 10.0, 100.0, 400.0];
    let y0 = vec![1.0; rates.len()];
    let tol = array::Tolerance {
        rtol: 1e-9,
        atol: 1e-12,
        ..array::Tolerance::default()
    };
    let mut launches = Vec::new();
    for dev in cpus() {
        let sol = array::ensemble_dopri5(&dev, &rhs, &[], &[&y0, &rates], 0.0, 1.0, &tol).unwrap();
        for (row, k) in rates.iter().enumerate() {
            let want = (-k).exp();
            assert!(
                (sol.y[0][row] - want).abs() <= 1e-6 * want + 1e-9,
                "k = {k}: {} vs {want}",
                sol.y[0][row]
            );
            assert_eq!(sol.y[1][row], *k);
        }
        launches.push(sol.launches);
    }
    assert_eq!(launches[0], launches[1]);
    // The steep rows dictate the number of launches; a gentle batch alone
    // needs far fewer, which is what per-row step control buys.
    let gentle = array::ensemble_dopri5(
        &CpuDevice::sequential(),
        &rhs,
        &[],
        &[&[1.0], &[0.1]],
        0.0,
        1.0,
        &tol,
    )
    .unwrap();
    assert!(
        gentle.launches * 4 < launches[0],
        "{} vs {}",
        gentle.launches,
        launches[0]
    );
    // A row's answer does not depend on its batch.
    let alone = array::ensemble_dopri5(
        &CpuDevice::sequential(),
        &rhs,
        &[],
        &[&[1.0], &[10.0]],
        0.0,
        1.0,
        &tol,
    )
    .unwrap();
    assert!((alone.y[0][0] - (-10.0f64).exp()).abs() < 1e-9);
    // Degenerate requests.
    let dev = CpuDevice::sequential();
    let same_time =
        array::ensemble_dopri5(&dev, &rhs, &[], &[&y0, &rates], 2.0, 2.0, &tol).unwrap();
    assert_eq!(
        (same_time.launches, same_time.y[0].clone()),
        (0, y0.clone())
    );
    assert!(array::ensemble_dopri5(&dev, &rhs, &[], &[&y0, &rates], 1.0, 0.0, &tol).is_err());
    // A blow-up (y' = y^2, y(0) = 1 escapes at t = 1) is reported, not looped on.
    let blow =
        KernelIr::from_tape(&tape(vec![Inst::Input(0), Inst::Mul(0, 0)], 1, 0, vec![1])).unwrap();
    let limited = array::Tolerance {
        max_launches: 2000,
        ..tol
    };
    assert!(matches!(
        array::ensemble_dopri5(&dev, &blow, &[], &[&[1.0]], 0.0, 2.0, &limited),
        Err(DeviceError::Backend(_))
    ));
}

// --- golden tests of the emitters -----------------------------------------------

fn golden(
    file: &str,
    got: &str,
) {
    let path = format!(
        "{}/src/backend/device/golden/{file}",
        env!("CARGO_MANIFEST_DIR")
    );
    // `RSSN_BLESS_GOLDEN=1 cargo test --lib backend::device` rewrites the files.
    if std::env::var_os("RSSN_BLESS_GOLDEN").is_some() {
        std::fs::write(path, got).unwrap();
        return;
    }
    let want = std::fs::read_to_string(&path)
        .unwrap_or_default()
        .replace("\r\n", "\n");
    assert!(
        got == want,
        "golden `{file}` differs.\n--- got ---\n{got}\n--- want ---\n{want}"
    );
}

fn axpy() -> KernelIr {
    // out0 = a * x + y  with  a = param 0
    let t = tape(
        vec![
            Inst::Input(0),
            Inst::Input(1),
            Inst::Param(0),
            Inst::Mul(2, 0),
            Inst::Add(3, 1),
        ],
        2,
        1,
        vec![4],
    );
    KernelIr::from_tape(&t).unwrap()
}

/// Pulls in every elementwise helper: powi, powf, cbrt, gather; plus
/// select/cmp/consts.
fn special() -> KernelIr {
    let t = tape(
        vec![
            Inst::Input(0),                  // 0
            Inst::Const(f64::NAN),           // 1
            Inst::Const(-2.5e-8),            // 2
            Inst::Const(f64::INFINITY),      // 3
            Inst::PowI(0, -3),               // 4
            Inst::PowF(0, 2),                // 5
            Inst::Unary(Intrinsic::Cbrt, 0), // 6
            Inst::Unary(Intrinsic::Ln, 0),   // 7
            Inst::Cmp(Cmp::Ge, 0, 2),        // 8
            Inst::Select(8, 4, 5),           // 9
            Inst::Gather(0, 3),              // 10
            Inst::Atan2(6, 7),               // 11
            Inst::Add(9, 10),                // 12
            Inst::Add(12, 11),               // 13
            Inst::Add(13, 1),                // 14
        ],
        1,
        0,
        vec![14, 6],
    );
    KernelIr::from_tape(&t).unwrap()
}

#[test]
fn wgsl_goldens() {
    let got = emit_wgsl(&axpy(), &WgslOptions::default()).unwrap();
    golden("axpy.wgsl", &got);
    let got = WgslEmitter {
        options: WgslOptions {
            workgroup_size: 128,
            entry: "kernel_main".into(),
        },
    }
    .emit(&special())
    .unwrap();
    golden("special.wgsl", &got);
    let pass = array::reduce_pass(ReduceOp::Min, 1, 0).unwrap();
    golden(
        "reduce_min.wgsl",
        &emit_wgsl(&pass, &WgslOptions::default()).unwrap(),
    );
    let spmv = array::spmv_csr_kernel().unwrap();
    golden(
        "spmv.wgsl",
        &emit_wgsl(&spmv, &WgslOptions::default()).unwrap(),
    );
}

#[test]
fn every_target_has_goldens_for_an_elementwise_and_a_folded_kernel() {
    use emit::EmitOptions;
    use emit::Precision;
    use emit::Target;
    let matmul = array::matmul_kernel().unwrap();
    for target in [Target::Glsl, Target::OpenCl, Target::Cuda, Target::Metal] {
        let options = EmitOptions::new(target);
        let ext = target.name();
        golden(
            &format!("special.{ext}"),
            &emit::emit(&special(), &options).unwrap(),
        );
        golden(
            &format!("matmul.{ext}"),
            &emit::emit(&matmul, &options).unwrap(),
        );
        let double = options.clone().precision(Precision::F64);
        if target.has_f64() {
            golden(
                &format!("matmul_f64.{ext}"),
                &emit::emit(&matmul, &double).unwrap(),
            );
        } else {
            assert!(matches!(
                emit::emit(&matmul, &double),
                Err(DeviceError::Unsupported(_))
            ));
        }
        let emitter = wgsl::SourceEmitter {
            options: options.clone(),
        };
        assert_eq!(emitter.target(), ext);
        assert!(
            emitter
                .emit(&axpy())
                .unwrap()
                .contains("rssn_out[0u * rssn_rows + row]")
        );
    }
    let wgsl_f64 = EmitOptions::new(Target::Wgsl).precision(Precision::F64);
    assert!(emit::emit(&matmul, &wgsl_f64).is_err());
}

#[test]
fn wgsl_emitter_rejects_bad_input_and_literals_round() {
    let bad = KernelIr {
        ops: vec![KernelOp::Neg(5)],
        inputs: 0,
        params: 0,
        outputs: vec![0],
        fold: None,
    };
    assert!(emit_wgsl(&bad, &WgslOptions::default()).is_err());
    let ok = decay();
    let zero = WgslOptions {
        workgroup_size: 0,
        entry: "m".into(),
    };
    assert!(emit_wgsl(&ok, &zero).is_err());
    assert_eq!(wgsl::f32_literal(0.5), "0.5");
    assert_eq!(wgsl::f32_literal(1e300), "bitcast<f32>(0x7f800000u)");
    assert_eq!(
        wgsl::f32_literal(f64::NEG_INFINITY),
        "bitcast<f32>(0xff800000u)"
    );
    assert_eq!(wgsl::f32_literal(f64::NAN), "bitcast<f32>(0x7fc00000u)");
}

#[test]
fn wgsl_of_the_rk4_step_has_the_expected_shape() {
    let step = rk4_step(&oscillator(), 2).unwrap();
    let w = emit_wgsl(&step, &WgslOptions::default()).unwrap();
    assert!(w.contains("rssn_out[1u * rssn_rows + row]"));
    assert!(w.contains("rssn_params[0u]") && w.contains("rssn_params[1u]"));
    assert!(w.contains("@compute @workgroup_size(64)"));
}

// --- wgpu (feature gpu-wgpu) ----------------------------------------------------

#[cfg(feature = "gpu-wgpu")]
mod gpu_tests {
    use super::*;

    fn structured() -> Vec<KernelIr> {
        use array::Boundary;
        vec![
            array::reduce_pass(ReduceOp::Sum, 1, 0).unwrap(),
            array::reduce_pass(ReduceOp::Max, 3, 2).unwrap(),
            array::map_reduce_pass(&smooth(), ReduceOp::Sum).unwrap(),
            array::matmul_kernel().unwrap(),
            array::spmv_csr_kernel().unwrap(),
            array::stencil1d_kernel(&[(-1, 1.0), (0, -2.0), (1, 1.0)], Boundary::Periodic).unwrap(),
            array::stencil2d_kernel(&[(0, 0, -4.0), (1, 0, 1.0), (0, -1, 1.0)], Boundary::Zero)
                .unwrap(),
            array::pairwise_kernel(&axpy_pair(), 1, ReduceOp::Sum, true).unwrap(),
            array::dopri5_step(&oscillator(), 2).unwrap(),
        ]
    }

    /// `a * x_i + x_j` as a pair kernel over one column.
    fn axpy_pair() -> KernelIr {
        axpy()
    }

    /// The emitted WGSL parses and validates under naga, the compiler wgpu
    /// itself uses. No adapter needed.
    #[test]
    fn emitted_wgsl_validates_under_naga() {
        use wgpu::naga;
        let step = rk4_step(&oscillator(), 2).unwrap();
        let mut programs = vec![everything(), decay(), step, special()];
        programs.extend(structured());
        for ir in programs {
            let src = emit_wgsl(&ir, &WgslOptions::default()).unwrap();
            let module = naga::front::wgsl::parse_str(&src)
                .unwrap_or_else(|e| panic!("{}\n{src}", e.emit_to_string(&src)));
            naga::valid::Validator::new(
                naga::valid::ValidationFlags::all(),
                naga::valid::Capabilities::all(),
            )
            .validate(&module)
            .unwrap_or_else(|e| panic!("{e:?}\n{src}"));
        }
    }

    /// Runs on whatever adapter wgpu finds (a software one counts); skipped
    /// when there is none.
    #[test]
    fn gpu_matches_cpu_within_f32_when_an_adapter_exists() {
        let devices = gpu::probe();
        let Some(dev) = devices.first() else {
            eprintln!("no wgpu adapter; skipping");
            return;
        };
        eprintln!("gpu test device: {}", dev.info().name);
        let ir = decay();
        let y0: Vec<f64> = (0..3000).map(|i| 0.5 + f64::from(i) * 1e-3).collect();
        let r = ensemble_rk4(&**dev, &ir, &[2.0], &[&y0], 0.0, 1.0, 100).unwrap();
        let e = (-2.0f64).exp();
        for (got, y) in r[0].iter().zip(&y0) {
            assert!(
                (got - y * e).abs() < 1e-4 * y.max(1.0),
                "{got} vs {}",
                y * e
            );
        }
        // Reductions and products are kernels, so the GPU has them.
        assert!(dev.info().ops.reduce && dev.info().ops.structured);
        let data: Vec<f64> = (0..70_001).map(|i| f64::from(i % 17) - 8.0).collect();
        let b = dev.upload(&data, Dtype::F32).unwrap();
        let want: f64 = data.iter().sum();
        assert_eq!(dev.reduce(ReduceOp::Sum, &b).unwrap(), want);
        assert_eq!(dev.reduce(ReduceOp::Max, &b).unwrap(), 8.0);
        assert_eq!(dev.reduce(ReduceOp::Min, &b).unwrap(), -8.0);
        // Matrix product against the CPU kernel.
        let (m, k, n) = (6, 4, 5);
        let a: Vec<f64> = (0..m * k).map(|i| f64::from(i as u8 % 7) - 3.0).collect();
        let bm: Vec<f64> = (0..k * n).map(|i| f64::from(i as u8 % 5) + 1.0).collect();
        let cpu = CpuDevice::sequential();
        let product = |d: &dyn Device, dtype| {
            let ba = d.upload(&a, dtype).unwrap();
            let bb = d.upload(&bm, dtype).unwrap();
            let mut bc = d.alloc(m * n, dtype).unwrap();
            d.matmul(&ba, &bb, &mut bc, m, k, n).unwrap();
            let mut host = vec![0.0; m * n];
            d.download(&bc, &mut host).unwrap();
            host
        };
        assert_eq!(product(&**dev, Dtype::F32), product(&cpu, Dtype::F64));
        // CSR times vector, stencil.
        let ptr = [0.0, 2.0, 2.0, 4.0, 5.0];
        let cols = [0.0, 2.0, 1.0, 2.0, 0.0];
        let vals = [2.0, 1.0, 3.0, 4.0, 5.0];
        let x = [1.0, 10.0, 100.0];
        let spmv = array::spmv_csr_kernel().unwrap();
        let y = array::run_cols(&**dev, &spmv, 4, &[&ptr, &cols, &vals, &x], &[]).unwrap();
        assert_eq!(y, vec![vec![102.0, 0.0, 430.0, 5.0]]);
        let u = [1.0, 4.0, 9.0, 16.0, 25.0];
        let st = array::stencil1d_kernel(&[(-1, 1.0), (0, -2.0), (1, 1.0)], array::Boundary::Clamp)
            .unwrap();
        let d2 = array::run_cols(&**dev, &st, 5, &[&u], &[]).unwrap();
        assert_eq!(d2, vec![vec![3.0, 2.0, 2.0, 2.0, -9.0]]);
        // Adaptive ensemble on the device, to f32 accuracy.
        let rhs = KernelIr::from_tape(&tape(
            vec![
                Inst::Input(0),
                Inst::Input(1),
                Inst::Mul(1, 0),
                Inst::Neg(2),
                Inst::Const(0.0),
            ],
            2,
            0,
            vec![3, 4],
        ))
        .unwrap();
        let rates = [0.5, 2.0, 20.0];
        let tol = array::Tolerance {
            rtol: 1e-4,
            atol: 1e-6,
            ..array::Tolerance::default()
        };
        let sol = array::ensemble_dopri5(&**dev, &rhs, &[], &[&[1.0; 3], &rates], 0.0, 1.0, &tol)
            .unwrap();
        for (row, k) in rates.iter().enumerate() {
            assert!(
                (sol.y[0][row] - (-k).exp()).abs() < 2e-4,
                "k = {k}: {}",
                sol.y[0][row]
            );
        }
    }
}
