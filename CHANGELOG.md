# RSSN changelog · 变更日志

Release changes to the RSSN kernel, algorithms and embedding APIs. 本文记录 RSSN 内核、算法和嵌入 API 的版本变化。

## [3.0.0-beta.1]

### Engine and embedding

- The scientific kernel exposes a deterministic manifest, type and dimension contracts, structured results and typed numerical callbacks. Embedded hosts generate checked interfaces from the operator registry that executes the calls.
- Indexed rewrites and incremental matching reduce repeated search. Independent requests use private computation graphs with ordered result import. Optional JIT and device backends define their supported precision and fallback behavior.
- Binary manifest, plugin and compute records use native bincode-next 3.2 traits. Signature and manifest bytes remain deterministic. Manifest verification reads the existing operator and unit tables in place and streams the canonical encoding into the hash. Compatibility tests compare the supported previous wire layouts with the native encoding. JSON and CBOR expose metadata to self-describing consumers.

### Numerical and symbolic algorithms

- Added weighted empirical quantiles and CDFs, discrete logarithms for composite moduli with non-invertible bases, and checked exact integer polynomial products. Improved LU factor reuse, eigensolver scaling, least-squares projections and parallel fractal correlation counts.
- The Rust LOBPCG API solves matrix-free generalized eigenproblems with a positive-definite mass operator and soft-locks converged vectors.
- Fixed-step ODE solvers reject invalid time grids and unrepresentable trajectory sizes before calling user functions. Event location reuses its state buffers.
- Binomial masses and beta densities handle support endpoints. A bounded binomial recurrence preserves low-count probabilities when the success probability is tiny. Binomial sampling uses discrete inversion with one uniform draw for nondegenerate parameters. Uniform density, CDF and sampling also support finite endpoints whose difference exceeds the floating-point range. Seeded binomial sample sequences change with the new algorithm.
- Improved extreme-tail probabilities for statistical tests and used thin SVD for ridge and OLS regression. Kendall correlation handles ties in O(n log n). Bootstrap procedures check sample and output sizes before allocating.
- Extended symbolic Mellin integrals, limits and infinite sums. Large-argument complex Faddeeva, error and Fresnel functions use scaled computations to avoid intermediate overflow.
- Poisson, heat and wave grid solvers accept per-face Dirichlet and Neumann conditions. Pure-Neumann Poisson checks the source–flux balance and returns a zero-mean solution. Structured kernel diagnostics identify rejected shapes, domains and convergence criteria.
- Three-dimensional shallow-water Euler and AMR runs accept a barotropic state in `(h, hu, hv, hw)`. The kernel manifest describes boundary conditions and tensor ranks for typed hosts.

### Build and compatibility

- Source builds require Rust 1.100 or newer. The default Rust artifact is an rlib. Build C libraries explicitly with their generated headers. C archives include English and Chinese instructions for linking directly from the extracted directory. [Build guide](docs/BUILDING.md)
- C batch computation accepts uninitialized output buffers. The `ffi` build profile supports converting unwinding Rust panics to C error codes at guarded entry points.
- The manifest API version stays at 2. Inspect the manifest hash before sharing serialized contracts across independently built kernels. [Kernel embedding](docs/KERNEL.md)

### 内核与嵌入

- 科学计算内核提供确定性的清单、类型与量纲契约、结构化结果以及类型化数值回调。宿主从实际执行调用的同一算子注册表生成受检查的接口。
- 重写索引和增量匹配减少重复搜索。独立请求使用各自的计算图，结果按输入顺序导入。可选 JIT 和设备后端明确支持的精度与回退行为。
- 二进制清单、插件消息和计算记录采用 bincode-next 3.2 原生编码。签名和清单保持确定性编码，校验清单时直接读取现有的算子表和单位表，逐段将规范编码传给哈希器。兼容性测试将支持的旧线格式与原生编码进行比较。JSON 与 CBOR 继续提供自描述元数据。

### 数值与符号算法

- 增加加权经验分位数和累积分布、含不可逆底数的复合模数离散对数，以及带溢出检查的精确整数多项式乘法。改进 LU 分解复用、特征值计算的缩放、最小二乘投影和并行分形关联计数。
- Rust LOBPCG API 支持带正定质量算子的无矩阵广义特征值问题，并软锁定已收敛的向量。
- 固定步长 ODE 求解器在调用用户函数前检查时间网格和轨迹容量。事件定位复用状态缓冲区。
- 二项分布概率质量与 Beta 密度正确处理支撑集端点。成功概率极小时，有界递推可保持小计数的二项概率。非退化二项分布使用一次均匀随机数反演采样。即使有限端点之差超出浮点范围，均匀分布的密度、累积分布与采样仍可求值。新算法会改变给定种子的二项采样序列。
- 改进统计检验的极端尾部概率，岭回归与普通最小二乘采用薄 SVD。Kendall 相关系数以 O(n log n) 处理并列值。自助法在分配内存前检查样本与输出大小。
- 扩展符号 Mellin 积分、极限和无穷级数。复数 Faddeeva 函数、误差函数和菲涅耳积分在大参数下采用缩放计算，避免中间溢出。
- 泊松、热传导和波动网格求解器支持按面指定 Dirichlet 或 Neumann 条件。全 Neumann 泊松问题检查源项与通量平衡并返回零均值解。结构化内核诊断指出不符合要求的形状、定义域或收敛准则。
- 三维浅水 Euler 和自适应网格细化求解器接收 `(h, hu, hv, hw)` 形式的正压状态。内核清单声明宿主所需的边界条件和张量秩约束。

### 构建与兼容性

- 源码构建需要 Rust 1.100 或更新版本。默认生成 Rust rlib。C 库按独立流程构建并携带生成的头文件。C 归档提供中英文说明，可从解压目录直接链接库和运行示例。[构建指南](docs/BUILDING.zh.md)
- C 批量计算允许传入未初始化的输出缓冲区。`ffi` 构建配置允许受保护的入口将可展开的 Rust panic 转为 C 错误码。
- 清单 API 版本保持为 2。不同构建的内核共享序列化契约前，应核对清单哈希。[内核嵌入](docs/KERNEL.md)

## [3.0.0-rc.1]

The release candidate establishes the graph rewriting engine, numerical ABI and standard operator registry on the 3.0 branch. The release notes above cover changes to this architecture through the beta series.

该候选版本在 3.0 分支建立图重写引擎、数值 ABI 与标准算子注册表。上方发布说明覆盖该架构在 beta 系列中的变更。
