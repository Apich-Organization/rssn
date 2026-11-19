# RSSN 执行后端与设备

[English](BACKENDS.md)

RSSN 将数值表达式的编译与批量内核执行分开。选择 `Backend` 将表达式编译成可调用函数。选择 `Device` 决定数据存放位置、计算精度和批量任务的执行设备。一行数据可以代表一个数据点、一个网格单元，或一条独立轨迹的初始状态。

| 接口 | 实现与用途 |
|---|---|
| `backend::Backend` | `Interpreter` 是默认标量后端。共享指令带通过 `TapeFunction` 实现 `Compiled` 和 `CompiledMulti`。启用 `jit` 后，`CraneliftBackend` 生成 CPU 机器码，`TieredBackend` 在调用次数或批量行数达到设定阈值时在后台编译标量函数。 |
| `backend::device::Device` | `cpu` 和 `cpu-parallel` 在宿主运行。启用 `gpu-wgpu` 后，每个可用的 wgpu 适配器可注册为一个 `gpu:wgpu:N` 设备。 |
| `backend::device::DeviceBackend` | 连接两套接口。批量调用会上载输入、在设备上启动已编译程序，再下载输出。单点调用启动一个仅处理一行的内核。 |

`Backend::compile` 根据声明的输入生成一个标量 `Compiled` 函数。`compile_multi` 生成 `CompiledMulti` 程序，支持共享子表达式、批内保持不变的参数和多个输出。`backend::current()` 依次选用当前线程最内层 `with_backend` 中的后端、`set_default` 设置的进程默认后端，以及解释器。
[后端 API](../src/backend/mod.rs) · [Cranelift 与分层执行](../src/backend/jit.rs)

表达式被降低为直线指令带（tape），其中共享子表达式并折叠常量。指令带支持宿主算子调用、微分和区间求值。`KernelIr::from_tape` 转换可移植的算术、内建函数、比较、选择和索引选列操作。遇到持有宿主函数指针的 `Call`、`Extern1` 或 `Extern2` 时会返回 `DeviceError::Unsupported`。如果表达式已成功生成指令带，而内核转换或设备编译失败，`DeviceBackend` 会用 `TapeFunction` 在宿主执行。设备运行期间，`Compiled` 系列接口用 NaN 表示设备错误。需要处理具体错误时，可直接调用 `Device` API。
[指令带](../src/backend/tape.rs) · [内核转换](../src/backend/device/kernel.rs) · [设备适配后端](../src/backend/device/mod.rs)

## 设备选择与缓冲区

`registry().list()` 列出运行时可用设备。可通过 `get(id)` 指定设备，或使用 `default_device()` 依次选择首个已注册的 GPU、`cpu-parallel`、首个可用设备。wgpu 会在名称中标明软件适配器。即使缺少可用适配器，两种 CPU 设备仍然可用。内核清单的 `features.devices` 通过 `build_devices()` 描述本次构建包含的设备支持，运行时探测单独通过注册表完成。
[设备注册表](../src/backend/device/mod.rs) · [适配器探测](../src/backend/device/gpu.rs) · [内核清单](../src/kernel/manifest.rs)

`DeviceInfo` 记录设备类型、支持的 `Dtype` 和操作、并行度、缓冲区容量、统一内存及确定性。`ops.structured` 表示设备支持独立列长、按地址读取和折叠的程序。设备未实现的可选操作返回 `DeviceError::Unsupported`。可通过 `emulated_f64` 识别 `F64` 的实现方式。`preferred_dtype` 在设备支持 `F64` 时选择它，否则选择 `F32`。

`Buffer` 归一个设备实例所有。`upload` 将宿主 `f64` 切片转换为指定的设备类型，`download` 将数据转换回 `f64`。`copy` 在同一设备、同一类型的缓冲区之间复制指定范围。默认实现通过宿主中转，CPU 和 wgpu 设备使用各自的内存复制。缓冲区归属错误、复制类型不一致或形状不匹配时，返回 `DeviceError::Mismatch`。

输入列连接在同一缓冲区中。等长输入的第 `i` 列从 `i * rows` 开始。第 `k` 个输出列从 `k * rows` 开始。`launch_cols` 通过 `lens` 指定每个输入列的长度，以长度的前缀和计算各列偏移。所有输出列仍有 `rows` 个元素。`launch` 用于等长输入。设备即使仅实现等长输入，仍可支持这种调用。共用启动检查将行数和输入元素总数限制为 `2^31 - 1`，适配器的缓冲区限制可能更低。每次启动会通过单独的宿主切片提供参数。
[设备和缓冲区约定](../src/backend/device/mod.rs)

对于可移植指令带和等长的宿主列，`map_columns` 完成数据传输和内核启动。

```rust
use rssn::backend::{device, tape::Tape};
use rssn::backend::device::{DeviceError, KernelIr};

fn evaluate(
    tape: &Tape,
    columns: &[&[f64]],
    params: &[f64],
) -> Result<Vec<Vec<f64>>, DeviceError> {
    let dev = device::registry()
        .get("cpu-parallel")
        .ok_or_else(|| DeviceError::Unsupported("cpu-parallel is unavailable".into()))?;
    let ir = KernelIr::from_tape(tape)?;
    device::map_columns(dev.as_ref(), &ir, columns, params)
}
```

重复启动时可以保留已编译内核和缓冲区。输出空间需要 `ir.outputs.len() * rows` 个元素。通过 `compile`、`upload` 和 `alloc` 准备后，用 `launch` 或 `launch_cols` 运行内核，在宿主需要结果时调用 `download`。wgpu 可异步提交。`download` 和 `synchronize` 会等待设备工作完成，也可报告延后收到的错误。
[wgpu 设备](../src/backend/device/gpu.rs)

## 内核中间表示与数组程序

`KernelIr` 在逐元素指令带基础上，增加整数索引寄存器、按地址读取和一个可选折叠。`Row`、`Rows`、`Len` 和整数运算使用 `Ty::Idx`，数值使用 `Ty::Num`。`ToIdx` 和 `ToNum` 在两种类型间转换。`Load(col, idx)` 读取指定列的单个元素，`Gather` 为当前行选取一个输入列。`Fold` 用 `Inner` 作为内部索引反复执行循环体，以 `Sum`、`Prod`、`Min` 或 `Max` 累积数值，之后通过 `Acc` 读取累积量。

`KernelIr::analyse` 校验操作数类型，并按数据依赖确定三个执行阶段。`Pre` 每行执行一次，包含循环不变量和折叠次数。`Body` 对每个内部索引执行一次。`Post` 读取累积值。寄存器不能同时依赖内部索引和累积量。`eval_row` 和 `eval_at` 定义宿主参考语义。
[内核中间表示与校验](../src/backend/device/kernel.rs)

数组库通过同一种设备启动接口实现下列算法。

| 计算 | 构造函数 | 行与每行工作 |
|---|---|---|
| 归约 | `reduce_pass` | 每行对应一个至多有 256 个元素的块，在块内折叠。 |
| 融合映射与归约 | `map_reduce_pass` | 在折叠内部映射各元素，保存各块的部分结果。 |
| 行主序矩阵乘法 | `matmul_kernel` | 每行对应一个结果元素，折叠得出内积。 |
| CSR 矩阵乘向量 | `spmv_csr_kernel` | 每行对应一个矩阵行，折叠该行存储的条目。 |
| 一维和二维模板运算 | `stencil1d_kernel`、`stencil2d_kernel` | 每行对应一个网格点，执行固定的读取操作，支持补零、夹紧或周期边界。 |
| 全对相互作用 | `pairwise_kernel` | 每行对应一个对象，折叠与其相互作用的对象。 |
| 自适应常微分方程轨迹组 | `dopri5_step` | 每行对应一条轨迹，用该行的时间与步长尝试一次嵌入式步进。 |

`array::run_cols` 在宿主提供的列上运行结构化内核。`reduce` 和 `reduce_column` 按每块最多 256 个元素逐层合并，只下载一个数值。`map_reduce` 将映射与第一轮归约合入一次启动。`array::matmul` 用 `Device::copy` 将两个行主序操作数打包到同一输入缓冲区。wgpu 设备通过这些驱动执行归约和矩阵乘法。CPU 设备的专用 `Device::reduce` 按固定的 4,096 个元素分块并依序合并，结果不受线程数影响。不同设备和调用入口的归约顺序与浮点行为可能不同。
[数组内核与驱动](../src/backend/device/array.rs) · [CPU 归约](../src/backend/device/cpu.rs)

`ensemble_rk4` 运行固定步数的经典四阶 Runge–Kutta。`rk4_step` 将四次右端函数求值内联到同一内核，每步启动一次。`ensemble_dopri5` 使用自适应 Dormand–Prince 5(4) 向前积分，每条轨迹将各自的时间和步长保留在设备上。每启动四个步进内核，宿主对时间列取一次最小值，最后下载各状态列。返回值还包括步进内核的启动次数。轨迹组在启动次数预算内未完成时，返回 `DeviceError::Backend`。
[固定步长驱动](../src/backend/device/mod.rs) · [自适应驱动](../src/backend/device/array.rs)

## CPU 执行与 LIR

串行和并行 CPU 设备使用原生 `f64`。逐元素内核的解释器每次由一条指令处理 128 行，让编译器可对访问连续数据的算术循环做向量化。超越函数使用标量数学运算。结构化内核使用参考求值器。并行设备通过 Rayon 递归划分互不重叠的输出范围。逐元素解释器和原生路径以 2,048 行为分组上限，结构化解释器以 128 行为分组上限。

启用 `jit` 时，达到 256 行的 CPU 设备任务会在首次使用时尝试编译行循环，并缓存结果。含 `Gather` 的内核使用解释器。原生路径因索引溢出等原因退出时，解释器会重跑整次启动。通常两条路径采用相同的索引和数值规则。`Min` 或 `Max` 折叠在比较数值相等、符号相反的零时，可能选中不同符号的零。
[CPU 执行](../src/backend/device/cpu.rs) · [原生降低](../src/backend/device/native.rs)

LIR 是指令带批处理、CPU 内核和语言前端共用的 CPU 中间表示。它支持基本块、跳转、分支和对模块内函数的直接调用。`F64`、`I64`、`Bool` 和 `Slice` 类型的寄存器可修改，初始为零。跳转时并行赋值给目标块参数。对宿主切片的读写会检查边界，越界读返回 NaN，越界写不产生效果。校验会沿调用关系分析写入，因此能在进入写入函数前拒绝只读切片参数。

整数运算与浮点转整数都带检查。检查失败、执行 `Term::Bail` 或模块调用深度超过 `MAX_CALL_DEPTH`（512）时，可返回 `Outcome::Bailed`。当计算读取的数据仍有效时，宿主可以重跑。`Module::interpret` 提供参考解释器。启用 `jit` 时，`lir_jit` 通过 Cranelift 编译模块。
[LIR](../src/backend/lir.rs) · [LIR 代码生成](../src/backend/lir_jit.rs)

| 前端 | LIR 用途 |
|---|---|
| 指令带 | `lir::batch_from_tape` 将符合条件的指令放入对列切片循环的函数中。达到 64 行时，Cranelift 批处理路径会尝试使用该入口。 |
| `KernelIr` | `device::native::lower` 把行循环、折叠和内核索引语义转换为 CPU 设备执行的 LIR。 |
| Aevia Core IR | 符合条件的纯标量 `@jit` 函数用 LIR 执行分支、带检查的 `Int` 运算、递归，以及受步数预算约束的自尾调用。 |

[后端计时示例](../examples/backend_timing.rs) 在 100,000 行数据上比较参考解释器、指令带逐行和分列执行，以及 Cranelift 的两种模式。可用 `cargo run --profile tools --features jit --example backend_timing` 测量当前机器上的执行时间。

## GPU 源代码目标与精度

`emit::emit(&KernelIr, &EmitOptions)` 生成指定目标的源代码。各目标共用校验、阶段划分、折叠结构和边界检查规则。`wgsl::emit_wgsl` 是生成单精度 WGSL 的便捷入口。

| 目标 | 算术 | 集成 |
|---|---|---|
| WGSL | `f32` 或双单精度 | 通过 `gpu-wgpu` 使用 wgpu 执行，有 naga 校验测试和依赖适配器的执行测试。 |
| GLSL 4.50 计算着色器 | `f32` 或原生 `f64` 源代码 | 源代码导出和基准文本测试。 |
| OpenCL C 1.2 | `f32` 或原生 `f64` 源代码 | 源代码导出和基准文本测试。 |
| CUDA C++ | `f32` 或原生 `f64` 源代码 | 源代码导出和基准文本测试。 |
| Metal Shading Language | `f32` | 源代码导出和基准文本测试。 |

后四种集成为外部运行时和编译器提供源代码，基准测试比较生成文本。WGSL 校验测试使用 naga 解析和验证生成模块，无需适配器。
[源代码生成器](../src/backend/device/emit.rs) · [生成器与设备测试](../src/backend/device/tests.rs)

每个着色器调用处理一行。WGSL 的第 0 组和 GLSL 的第 0 个描述符集都使用五个绑定，依次为输入、参数、输出、尺寸（`rows, inputs, outputs, pad`）和列元数据（`offset, length` 对）。Metal 使用对应的缓冲区编号。OpenCL 和 CUDA 依次接收输入、参数、输出、列元数据、行数和输入列数等内核参数。wgpu 用存储缓冲区传递参数，用 uniform 传递尺寸。每个工作组处理 64 行。首轴工作组超过 65,535 个时，分派扩展至二维。

wgpu 按结构哈希缓存已编译内核，命中时检查原始中间表示是否一致。内核编译时生成单精度管线，首次用 `F64` 启动时再生成双单精度管线。wgpu 中，非空任务的输入和输出应使用同一种元素类型。报告的 `max_buffer_elems` 以 8 字节元素为单位保守计算。
[wgpu 设备](../src/backend/device/gpu.rs)

### 双单精度算术

wgpu 上的 `Dtype::F64` 使用两个 `f32` 分量表示 `hi + lo`。上传时拆分宿主数据，下载时将分量相加。着色器输入、参数和输出是由 8 字节数对构成的数组。这种表示约有 48 个有效位，指数范围与 `f32` 相同。`DeviceInfo::emulated_f64` 标记这种实现。`dtype_epsilon` 对 `F32`、原生 `F64` 和数对实现分别返回 `2^-23`、`2^-52` 和 `2^-44`。自适应轨迹组据此设置终点容差。

数对算术支持基本运算、比较、平方根、幂、取整和初等函数。超越函数以数对精度通过区间缩减、级数或 Newton 修正求值。依赖适配器的测试在代表性输入上与 CPU 对照，绝对误差界为 `5e-12 * max(abs(reference), 1e-3)`。实际精度随计算任务和驱动程序而变。使用双单精度时，非零数据和中间值的量级宜控制在约 `1e-30` 至 `1e30`。三角函数的绝对误差随参数增大而增长，`sin`、`cos` 和 `tan` 在 `|x| > 1e6` 时返回 NaN。GPU 对 NaN 和无穷大的处理取决于实现。

数对库通过 `ds_hide` 保护中间舍入点。它与尺寸 uniform 提供的零执行按位运算，用于防止着色器优化器重排运算或合并运算，保留无误差变换所需的舍入行为。修改数对算术时应保留这道屏障。

着色器的内核索引寄存器使用有符号 32 位整数。CSR 行指针和列索引等通过数值列传入的索引还要满足缓冲区的数值表示范围。`F32` 在 `2^24` 以上无法区分所有相邻整数。双单精度能精确表示更大范围内的整数，同时仍受着色器索引转换和设备边界的约束。
[数对算术](../src/backend/device/ds.rs) · [精度测试](../src/backend/device/tests.rs)

## 内核启动融合

`device::plan::Plan` 由等长外部列和先前内核的输出构建内核图。用 `result` 指定所需输出，调用 `fuse` 合并符合条件的启动，再通过 `run` 在选定设备上执行。对于 `y = f(x); z = g(y, w)` 这样的链，融合可用寄存器保存中间值 `y`。其他消费者仍需使用的列，以及显式请求的最终结果，可以保留为合并内核的输出。

融合会代入同一行的 `Input` 读取，重新编号参数和寄存器，并复用重复的输入源。两个内核都有折叠、消费者执行 `Gather`，或消费者在生产者的列上执行按地址读取的 `Load` 或 `Len` 时，保留独立的内核启动。合并内核须通过校验，并维持执行顺序。若还有较早执行的内核读取生产者的输出，则保留生产者原来的执行位置。对于等长列，计划会处理启动之间的设备复制。输入列长度不同的程序使用 `launch_cols`。
[计划器与融合测试](../src/backend/device/plan.rs)

## 实现设备

1. 用唯一 ID 向 `DeviceRegistry` 注册设备，并通过 `DeviceInfo` 报告已实现的精度、操作和限制。内置且受功能开关控制的 GPU 使用 `gpu:<api>:<ordinal>`，通过 `build_devices()` 发布构建能力。找不到可初始化设备时，探测返回空列表。
2. 用 `Storage` 封装存储，用 `next_uid()` 标识设备实例。在数据传输和内核启动入口检查归属、形状和类型。设备能直接移动数据时，可重写 `copy` 以省去宿主传输。
3. 校验并编译 `KernelIr`。目标匹配时可复用源代码生成器。按结构哈希缓存，复用前检查中间表示是否相同。通过 `DeviceError` 报告校验失败和后端故障。异步运行时应在 `download` 和 `synchronize` 收集延后收到的错误。
4. 实现并声明支持的内核启动与归约能力。用 `check_launch` 或 `check_launch_cols` 检查共用的缓冲区约定。结构化设备可用数组驱动实现归约和矩阵乘法。以声明的精度和参考结果对照，并覆盖外部设备缓冲区和不匹配的缓冲区大小。
5. 让生成器的源代码测试独立于硬件。基准输出位于设备目录，可通过 `RSSN_BLESS_GOLDEN=1` 重新生成。目标编译器可用时应验证生成代码，探测到合适的适配器时应执行设备测试。新增初等函数 `Intrinsic` 时需定义指令带和 CPU 语义，并为各源代码目标添加映射或辅助函数。
