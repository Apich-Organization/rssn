//! Checked CPU compilation keeps validated kernels across repeat and concurrent requests.
use std::sync::Arc;
use std::sync::Barrier;

use rssn::backend::device::CpuDevice;
use rssn::backend::device::Device;
use rssn::backend::device::Dtype;
use rssn::backend::device::KernelIr;
use rssn::backend::device::KernelOp;

fn multiply_kernel(multiplier: f64) -> KernelIr {
    KernelIr {
        ops: vec![
            KernelOp::Input(0),
            KernelOp::Const(multiplier),
            KernelOp::Mul(0, 1),
        ],
        inputs: 1,
        params: 0,
        outputs: vec![2],
        fold: None,
    }
}

#[test]
fn validated_cache_hits_preserve_the_kernel_and_reject_invalid_programs() {
    let device = CpuDevice::sequential();
    let program = multiply_kernel(2.0);
    let first = device.compile(&program).unwrap();
    let mut invalid = program.clone();
    invalid.outputs = vec![u32::MAX];
    assert!(device.compile(&invalid).is_err());
    assert!(Arc::ptr_eq(&first, &device.compile(&program).unwrap()));

    let other = device.compile(&multiply_kernel(3.0)).unwrap();
    assert!(!Arc::ptr_eq(&first, &other));
    let inputs = device.upload(&[1.0, -4.0], Dtype::F64).unwrap();
    let mut output = device.alloc(2, Dtype::F64).unwrap();
    let mut host = [0.0; 2];
    for (kernel, expected) in [(first, [2.0, -8.0]), (other, [3.0, -12.0])] {
        device
            .launch(&*kernel, 2, &inputs, &[], &mut output)
            .unwrap();
        device.download(&output, &mut host).unwrap();
        assert_eq!(host, expected);
    }
}

#[test]
fn concurrent_cold_compiles_publish_one_shared_kernel() {
    const CALLERS: usize = 6;
    let device = CpuDevice::parallel();
    let barrier = Barrier::new(CALLERS);
    let program = multiply_kernel(5.0);
    std::thread::scope(|scope| {
        let requests: Vec<_> = (0..CALLERS)
            .map(|_| {
                scope.spawn(|| {
                    barrier.wait();
                    device.compile(&program).unwrap()
                })
            })
            .collect();
        let installed = device.compile(&program).unwrap();
        for request in requests {
            assert!(Arc::ptr_eq(&installed, &request.join().unwrap()));
        }
    });
}
