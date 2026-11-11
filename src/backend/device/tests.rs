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
    let tape = ir.to_tape();
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

// --- WGSL golden tests --------------------------------------------------------

fn golden(
    name: &str,
    got: &str,
    want: &str,
) {
    // `RSSN_BLESS_GOLDEN=1 cargo test --lib backend::device` rewrites the files.
    if std::env::var_os("RSSN_BLESS_GOLDEN").is_some() {
        let path = format!(
            "{}/src/backend/device/golden/{name}.wgsl",
            env!("CARGO_MANIFEST_DIR")
        );
        std::fs::write(path, got).unwrap();
        return;
    }
    assert!(
        got == want,
        "WGSL golden `{name}` differs.\n--- got ---\n{got}\n--- want ---\n{want}"
    );
}

#[test]
fn wgsl_axpy_golden() {
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
    let ir = KernelIr::from_tape(&t).unwrap();
    let got = emit_wgsl(&ir, &WgslOptions::default()).unwrap();
    golden("axpy", &got, include_str!("golden/axpy.wgsl"));
}

#[test]
fn wgsl_special_ops_golden() {
    // Pulls in every helper: powi, powf, cbrt, gather; plus select/cmp/consts.
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
    let ir = KernelIr::from_tape(&t).unwrap();
    let got = WgslEmitter {
        options: WgslOptions {
            workgroup_size: 128,
            entry: "kernel_main".into(),
        },
    }
    .emit(&ir)
    .unwrap();
    golden("special", &got, include_str!("golden/special.wgsl"));
}

#[test]
fn wgsl_emitter_rejects_bad_input_and_literals_round() {
    let bad = KernelIr {
        ops: vec![KernelOp::Neg(5)],
        inputs: 0,
        params: 0,
        outputs: vec![0],
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
    assert!(w.contains("rssn_out[1u * rssn_dims.rows + row]"));
    assert!(w.contains("rssn_params[0u]") && w.contains("rssn_params[1u]"));
    assert!(w.contains("@compute @workgroup_size(64)"));
}

// --- wgpu (feature gpu-wgpu) ----------------------------------------------------

#[cfg(feature = "gpu-wgpu")]
mod gpu_tests {
    use super::*;

    /// The emitted WGSL parses and validates under naga, the compiler wgpu
    /// itself uses. No adapter needed.
    #[test]
    fn emitted_wgsl_validates_under_naga() {
        use wgpu::naga;
        let step = rk4_step(&oscillator(), 2).unwrap();
        for ir in [everything(), decay(), step] {
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
        let b = dev.upload(&[1.0], Dtype::F32).unwrap();
        assert!(dev.reduce(ReduceOp::Sum, &b).is_err());
    }
}
