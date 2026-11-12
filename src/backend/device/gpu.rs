//! The wgpu device (feature `gpu-wgpu`).
//!
//! Kernels are emitted as WGSL by [`super::wgsl`] and run as compute shaders.
//! Computation is `f32` ([`Dtype::F32`] is the only supported type), so
//! results differ from the CPU devices by single-precision rounding; reduce
//! and matmul are not implemented (`DeviceOps` says so and the calls return
//! [`DeviceError::Unsupported`]).
//!
//! Adapters are discovered by [`probe`], which never fails: a machine with no
//! usable adapter simply registers no GPU device. Validation errors raised
//! while creating a pipeline are caught in an error scope and returned as
//! [`DeviceError::Backend`]; errors that wgpu reports later (uncaptured) are
//! remembered and returned by the next [`Device::synchronize`] /
//! [`Device::download`].
//!
//! There is no `unsafe` here: byte views of `f32` slices go through
//! `bytemuck`.

use std::any::Any;
use std::collections::HashMap;
use std::sync::Arc;
use std::sync::Mutex;

use super::Buffer;
use super::CompiledKernel;
use super::Device;
use super::DeviceError;
use super::DeviceInfo;
use super::DeviceKind;
use super::DeviceOps;
use super::Dtype;
use super::GpuApi;
use super::KernelIr;
use super::Storage;
use super::check_launch_cols;
use super::next_uid;
use super::wgsl::WgslOptions;
use super::wgsl::emit_wgsl;

const WORKGROUP: u32 = 64;
const MAX_GROUPS_X: u32 = 65_535;

struct GpuStorage {
    buffer: wgpu::Buffer,
}

impl Storage for GpuStorage {
    fn as_any(&self) -> &dyn Any {
        self
    }

    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }
}

struct GpuKernel {
    ir: KernelIr,
    pipeline: wgpu::ComputePipeline,
}

impl CompiledKernel for GpuKernel {
    fn ir(&self) -> &KernelIr {
        &self.ir
    }

    fn as_any(&self) -> &dyn Any {
        self
    }
}

/// One wgpu adapter as a [`Device`].
pub struct WgpuDevice {
    info: DeviceInfo,
    uid: u64,
    device: wgpu::Device,
    queue: wgpu::Queue,
    layout: wgpu::BindGroupLayout,
    pipeline_layout: wgpu::PipelineLayout,
    cache: Mutex<HashMap<u64, Arc<GpuKernel>>>,
    late_error: Arc<Mutex<Option<String>>>,
}

/// Every usable wgpu adapter, as devices `gpu:wgpu:0`, `gpu:wgpu:1`, ...
/// Software rasterisers (CPU adapters) are included and marked in the name.
#[must_use]
pub fn probe() -> Vec<Arc<dyn Device>> {
    let instance = wgpu::Instance::new(&wgpu::InstanceDescriptor::default());
    let mut out: Vec<Arc<dyn Device>> = Vec::new();
    for (n, adapter) in instance
        .enumerate_adapters(wgpu::Backends::all())
        .into_iter()
        .enumerate()
    {
        if let Ok(device) = WgpuDevice::new(&adapter, n) {
            out.push(Arc::new(device));
        }
    }
    out
}

fn storage_entry(
    binding: u32,
    read_only: bool,
) -> wgpu::BindGroupLayoutEntry {
    wgpu::BindGroupLayoutEntry {
        binding,
        visibility: wgpu::ShaderStages::COMPUTE,
        ty: wgpu::BindingType::Buffer {
            ty: wgpu::BufferBindingType::Storage { read_only },
            has_dynamic_offset: false,
            min_binding_size: None,
        },
        count: None,
    }
}

impl WgpuDevice {
    /// Opens `adapter`.
    ///
    /// # Errors
    /// [`DeviceError::Backend`] when the device cannot be created.
    pub fn new(
        adapter: &wgpu::Adapter,
        ordinal: usize,
    ) -> Result<Self, DeviceError> {
        let ainfo = adapter.get_info();
        let limits = adapter.limits();
        let (device, queue) =
            pollster::block_on(adapter.request_device(&wgpu::DeviceDescriptor::default()))
                .map_err(|e| DeviceError::Backend(e.to_string()))?;
        let late_error = Arc::new(Mutex::new(None));
        {
            let sink = late_error.clone();
            device.on_uncaptured_error(Arc::new(move |e: wgpu::Error| {
                if let Ok(mut slot) = sink.lock() {
                    *slot = Some(e.to_string());
                }
            }));
        }
        let layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("rssn"),
            entries: &[
                storage_entry(0, true),
                storage_entry(1, true),
                storage_entry(2, false),
                wgpu::BindGroupLayoutEntry {
                    binding: 3,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Uniform,
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
                // The column table: (offset, length) per input column.
                storage_entry(4, true),
            ],
        });
        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("rssn"),
            bind_group_layouts: &[&layout],
            push_constant_ranges: &[],
        });
        let software = matches!(ainfo.device_type, wgpu::DeviceType::Cpu);
        let info = DeviceInfo {
            id: format!("gpu:wgpu:{ordinal}"),
            name: format!(
                "{}{} ({:?})",
                ainfo.name,
                if software {
                    " [software]"
                } else {
                    ""
                },
                ainfo.backend
            ),
            kind: DeviceKind::Gpu(GpuApi::Wgpu),
            dtypes: vec![Dtype::F32],
            ops: DeviceOps {
                elementwise: true,
                // Both are kernels (`super::array`), not hand-written
                // shaders: a device that runs structured kernels has them.
                reduce: true,
                matmul: true,
                structured: true,
            },
            parallelism: limits.max_compute_invocations_per_workgroup as usize,
            max_buffer_elems: usize::try_from(limits.max_storage_buffer_binding_size / 4)
                .unwrap_or(usize::MAX),
            unified_memory: matches!(ainfo.device_type, wgpu::DeviceType::Cpu),
            deterministic: true,
        };
        Ok(Self {
            info,
            uid: next_uid(),
            device,
            queue,
            layout,
            pipeline_layout,
            cache: Mutex::new(HashMap::new()),
            late_error,
        })
    }

    fn take_late_error(&self) -> Result<(), DeviceError> {
        let e = self.late_error.lock().ok().and_then(|mut s| s.take());
        e.map_or(Ok(()), |m| Err(DeviceError::Backend(m)))
    }

    fn gpu_buffer<'a>(
        &self,
        b: &'a Buffer,
    ) -> Result<&'a wgpu::Buffer, DeviceError> {
        if b.owner() != self.uid {
            return Err(DeviceError::Mismatch(
                "buffer belongs to another device".into(),
            ));
        }
        b.storage()
            .as_any()
            .downcast_ref::<GpuStorage>()
            .map(|s| &s.buffer)
            .ok_or_else(|| DeviceError::Mismatch("not a wgpu buffer".into()))
    }

    fn new_storage(
        &self,
        len: usize,
    ) -> Result<Buffer, DeviceError> {
        let bytes = (len.max(1) as u64).saturating_mul(4);
        let buffer = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("rssn"),
            size: bytes,
            usage: wgpu::BufferUsages::STORAGE
                | wgpu::BufferUsages::COPY_SRC
                | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        Ok(Buffer::new(
            self.uid,
            Dtype::F32,
            len,
            Box::new(GpuStorage { buffer }),
        ))
    }

    fn small(
        &self,
        usage: wgpu::BufferUsages,
        bytes: &[u8],
    ) -> wgpu::Buffer {
        let size = (bytes.len().max(4) as u64).next_multiple_of(4);
        let buffer = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("rssn"),
            size,
            usage: usage | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        if !bytes.is_empty() {
            self.queue.write_buffer(&buffer, 0, bytes);
        }
        buffer
    }
}

#[allow(clippy::cast_possible_truncation)]
fn to_f32(v: &[f64]) -> Vec<f32> {
    v.iter().map(|&x| x as f32).collect()
}

impl Device for WgpuDevice {
    fn info(&self) -> &DeviceInfo {
        &self.info
    }

    fn uid(&self) -> u64 {
        self.uid
    }

    fn alloc(
        &self,
        len: usize,
        dtype: Dtype,
    ) -> Result<Buffer, DeviceError> {
        if dtype != Dtype::F32 {
            return Err(DeviceError::Unsupported(
                "wgpu devices compute in f32".into(),
            ));
        }
        self.new_storage(len)
    }

    fn upload(
        &self,
        host: &[f64],
        dtype: Dtype,
    ) -> Result<Buffer, DeviceError> {
        let b = self.alloc(host.len(), dtype)?;
        if !host.is_empty() {
            let data = to_f32(host);
            self.queue
                .write_buffer(self.gpu_buffer(&b)?, 0, bytemuck::cast_slice(&data));
        }
        Ok(b)
    }

    fn download(
        &self,
        buffer: &Buffer,
        out: &mut [f64],
    ) -> Result<(), DeviceError> {
        let src = self.gpu_buffer(buffer)?;
        if buffer.len() != out.len() {
            return Err(DeviceError::Mismatch(format!(
                "buffer has {} elements, destination {}",
                buffer.len(),
                out.len()
            )));
        }
        if out.is_empty() {
            return Ok(());
        }
        let bytes = (out.len() as u64) * 4;
        let staging = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("rssn-readback"),
            size: bytes,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let mut enc = self
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor { label: None });
        enc.copy_buffer_to_buffer(src, 0, &staging, 0, bytes);
        self.queue.submit([enc.finish()]);
        let (tx, rx) = std::sync::mpsc::channel();
        staging.slice(..).map_async(wgpu::MapMode::Read, move |r| {
            let _ = tx.send(r);
        });
        self.device
            .poll(wgpu::PollType::wait_indefinitely())
            .map_err(|e| DeviceError::Backend(format!("{e:?}")))?;
        rx.recv()
            .map_err(|e| DeviceError::Backend(e.to_string()))?
            .map_err(|e| DeviceError::Backend(e.to_string()))?;
        {
            let view = staging.slice(..).get_mapped_range();
            let floats: &[f32] = bytemuck::cast_slice(&view);
            for (o, &f) in out.iter_mut().zip(floats) {
                *o = f64::from(f);
            }
        }
        staging.unmap();
        self.take_late_error()
    }

    fn copy(
        &self,
        src: &Buffer,
        src_at: usize,
        dst: &mut Buffer,
        dst_at: usize,
        len: usize,
    ) -> Result<(), DeviceError> {
        super::check_copy(self.uid, src, src_at, dst, dst_at, len)?;
        if len == 0 {
            return Ok(());
        }
        let mut enc = self
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor { label: None });
        enc.copy_buffer_to_buffer(
            self.gpu_buffer(src)?,
            src_at as u64 * 4,
            self.gpu_buffer(dst)?,
            dst_at as u64 * 4,
            len as u64 * 4,
        );
        self.queue.submit([enc.finish()]);
        Ok(())
    }

    fn reduce(
        &self,
        op: super::ReduceOp,
        input: &Buffer,
    ) -> Result<f64, DeviceError> {
        super::array::reduce(self, op, input)
    }

    fn matmul(
        &self,
        a: &Buffer,
        b: &Buffer,
        c: &mut Buffer,
        m: usize,
        k: usize,
        n: usize,
    ) -> Result<(), DeviceError> {
        super::array::matmul(self, a, b, c, m, k, n)
    }

    fn compile(
        &self,
        ir: &KernelIr,
    ) -> Result<Arc<dyn CompiledKernel>, DeviceError> {
        let key = ir.structural_hash();
        {
            let cache = self
                .cache
                .lock()
                .map_err(|_| DeviceError::Backend("kernel cache poisoned".into()))?;
            if let Some(k) = cache.get(&key).filter(|k| k.ir == *ir) {
                return Ok(k.clone());
            }
        }
        let options = WgslOptions {
            workgroup_size: WORKGROUP,
            entry: "main".into(),
        };
        let source = emit_wgsl(ir, &options)?;
        self.device.push_error_scope(wgpu::ErrorFilter::Validation);
        let module = self
            .device
            .create_shader_module(wgpu::ShaderModuleDescriptor {
                label: Some("rssn"),
                source: wgpu::ShaderSource::Wgsl(source.into()),
            });
        let pipeline = self
            .device
            .create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                label: Some("rssn"),
                layout: Some(&self.pipeline_layout),
                module: &module,
                entry_point: Some("main"),
                compilation_options: wgpu::PipelineCompilationOptions::default(),
                cache: None,
            });
        if let Some(e) = pollster::block_on(self.device.pop_error_scope()) {
            return Err(DeviceError::Backend(e.to_string()));
        }
        let k = Arc::new(GpuKernel {
            ir: ir.clone(),
            pipeline,
        });
        self.cache
            .lock()
            .map_err(|_| DeviceError::Backend("kernel cache poisoned".into()))?
            .insert(key, k.clone());
        Ok(k)
    }

    fn launch(
        &self,
        kernel: &dyn CompiledKernel,
        rows: usize,
        input: &Buffer,
        params: &[f64],
        output: &mut Buffer,
    ) -> Result<(), DeviceError> {
        let lens = vec![rows; kernel.ir().inputs];
        self.launch_cols(kernel, rows, input, &lens, params, output)
    }

    fn launch_cols(
        &self,
        kernel: &dyn CompiledKernel,
        rows: usize,
        input: &Buffer,
        lens: &[usize],
        params: &[f64],
        output: &mut Buffer,
    ) -> Result<(), DeviceError> {
        let ir = kernel.ir();
        check_launch_cols(self.uid, ir, rows, input, lens, output)?;
        let k = kernel
            .as_any()
            .downcast_ref::<GpuKernel>()
            .ok_or_else(|| DeviceError::Mismatch("kernel compiled by another device".into()))?;
        if rows == 0 {
            return Ok(());
        }
        let rows32 = u32::try_from(rows)
            .map_err(|_| DeviceError::Unsupported("more than 2^32 rows".into()))?;
        let p = to_f32(params);
        let pbuf = self.small(wgpu::BufferUsages::STORAGE, bytemuck::cast_slice(&p));
        let dims: [u32; 4] = [
            rows32,
            u32::try_from(ir.inputs).unwrap_or(u32::MAX),
            u32::try_from(ir.outputs.len()).unwrap_or(u32::MAX),
            0,
        ];
        let dbuf = self.small(wgpu::BufferUsages::UNIFORM, bytemuck::cast_slice(&dims));
        // `check_launch_cols` bounded every length by 2^31 - 1.
        let mut table: Vec<u32> = Vec::with_capacity(lens.len() * 2);
        let mut at = 0_u32;
        for &len in lens {
            let len = u32::try_from(len).unwrap_or(u32::MAX);
            table.extend([at, len]);
            at = at.saturating_add(len);
        }
        let cbuf = self.small(wgpu::BufferUsages::STORAGE, bytemuck::cast_slice(&table));
        let bind = self.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("rssn"),
            layout: &self.layout,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: self.gpu_buffer(input)?.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: pbuf.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: self.gpu_buffer(output)?.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 3,
                    resource: dbuf.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 4,
                    resource: cbuf.as_entire_binding(),
                },
            ],
        });
        let groups = rows32.div_ceil(WORKGROUP);
        let gx = groups.min(MAX_GROUPS_X).max(1);
        let gy = groups.div_ceil(gx).max(1);
        let mut enc = self
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor { label: None });
        {
            let mut pass = enc.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: None,
                timestamp_writes: None,
            });
            pass.set_pipeline(&k.pipeline);
            pass.set_bind_group(0, &bind, &[]);
            pass.dispatch_workgroups(gx, gy, 1);
        }
        self.queue.submit([enc.finish()]);
        Ok(())
    }

    fn synchronize(&self) -> Result<(), DeviceError> {
        self.device
            .poll(wgpu::PollType::wait_indefinitely())
            .map_err(|e| DeviceError::Backend(format!("{e:?}")))?;
        self.take_late_error()
    }
}
