# RSSN changelog · 变更日志

Release changes to the RSSN kernel, algorithms and embedding APIs. 本文记录 RSSN 内核、算法和嵌入 API 的版本变化。

## [3.0.0-beta.1]

### Engine and embedding

- The scientific kernel exposes a deterministic manifest, type and dimension contracts, structured results and typed numerical callbacks. Embedded hosts generate checked interfaces from the operator registry that executes the calls.
- Indexed rewrites and incremental matching reduce repeated search. Independent requests use private computation graphs with ordered result import. Optional JIT and device backends define their supported precision and fallback behavior.
- Workers prepare and evaluate symbolic batches on private graphs with bounded in-flight jobs. CPU dispatch reads the active Rayon pool size, and FFT routines share plans and scratch storage through `kernels::reuse`.
- Symbolic windows reuse extraction work across overlapping terms. Compiled CPU batches lease sized scratch buffers and skip parallel dispatch for small work. Concurrent FFT requests share a plan build for the same length and direction.
- Column-wise tapes read the selected source column for each gather. A prepared gather tape can take measured source columns with a common reference and return values in their original dimension.
- Prepared scalar programs can take one measured column per input, each with its own unit. Large batches run on the active CPU pool in fixed, ordered blocks and honor the embedding cancellation scope.
- Standard sessions share installed rules and their manifests until the global unit table changes. Per-session units, term graphs, native objects, and computation caches remain session-owned.
- The C ABI captures the shared standard manifest and reuses its JSON encoding. Each call still returns a separately owned string.
- Native object bindings derive runtime arity from the signature's argument list. Dense value preflight checks bounds without allocating a traversal buffer.
- Binary manifest, plugin and compute records use native bincode-next 3.2 traits. Signature and manifest bytes remain deterministic. Manifest verification reads the existing operator and unit tables in place and streams the canonical encoding into the hash. Compatibility tests compare the supported previous wire layouts with the native encoding. JSON and CBOR expose metadata to self-describing consumers.

### Numerical and symbolic algorithms

- Add stiff one-dimensional method-of-lines integration with time-dependent Dirichlet and Neumann boundaries, including measured callback programs for embedded hosts.
- Correct finite-temperature XXZ stochastic series expansion sampling. The corrected update is reproducible for a fixed seed and changes the trajectories produced by earlier builds.
- Add one-sided Welch spectral density, inverse STFT with overlap normalization, and polyphase FIR resampling. Add checked topological ordering and a nonrecursive strongly connected component API shared with symbolic graph operators.
- Add analytic signals, Hilbert transforms, zoom FFT, spectral coherence and wavelet packets. DCT and DST use FFT-based algorithms and fractional Fourier transforms can reuse their spectral basis.
- Add overlap-save convolution, autocorrelation, FIR filter design, frequency response, second-order-section filtering, polyphase resampling and stationary wavelet transforms with coefficient thresholding.
- Add batched quantiles and empirical CDF, Benjamini–Hochberg, Benjamini–Yekutieli and Holm adjustment, batched multivariate normal log densities, and atomic covariance batch updates with correlation output.
- Add signed shortest paths, weak components, maximum flow with its minimum cut, and column interpolative decomposition. The kernel manifest carries the numerical inputs, structured results and dimensional contracts.
- Add multi-source BFS, strong-component condensation, articulation points and bridges, and three sample time-series lag estimators. Large Welch workloads reduce fixed frame blocks in parallel. Dormand–Prince, conjugate gradients and restarted GMRES reuse inner-loop workspaces.
- Optimized dense multiplication uses row-major input and output views. Large products use the active Rayon pool.
- Euclidean norms, distance and vector directions scale intermediate squares at extreme ranges. Distance no longer materializes the difference vector, and Rust callers can normalize into their own buffers.
- Vector reflection writes one output buffer and preserves finite components when an intermediate correction overflows. Arnoldi scales its basis construction and reuses projection storage. Native dense eigensolvers write matrix-vector products into the iteration workspace.
- Bilinear interpolation uses the same scaled interval coordinates for single queries and prepared grids, including finite endpoints whose difference exceeds Float range.
- Prepared Newton–CG, polynomial antidifferences and geometric grid transfers retain reusable numerical state. Implicit Gauss–Legendre steps use a separate Jacobian at each collocation stage.
- The library adds logarithmic Bessel K for real orders and exact Chebyshev residues with full-width unsigned moduli. Sparse exact elimination reuses row storage and indexes pivot candidates on larger systems.
- Block sparse multiplication reuses row workspaces and accepts batches of dense right-hand sides. Repeated moderate Zoom FFT calls share the transformed chirp for the same grid.
- Drift-implicit SDE methods reuse Newton workspaces and provide a sampled path. PCHIP and Akima batch queries reuse segment coefficients. Online covariance accepts exact integer frequencies, and integer kernels provide a Möbius sieve and linear congruence solver.
- Sparse conjugate gradients rescale extreme right-hand sides. Fixed-step delay solvers reject time grids that cannot advance, and the wave stencil combines speed and step before squaring.

- Added weighted empirical quantiles and CDFs, discrete logarithms for composite moduli with non-invertible bases, and checked exact integer polynomial products. Improved LU factor reuse, eigensolver scaling, least-squares projections and parallel fractal correlation counts.
- The Rust LOBPCG API solves matrix-free generalized eigenproblems with a positive-definite mass operator and soft-locks converged vectors.
- Fixed-step ODE solvers reject invalid time grids and unrepresentable trajectory sizes before calling user functions. Event location reuses its state buffers.
- The single-pole Meijer G reduction computes the power and exponential on a common scale and preserves subnormal results. The commutative Platen solver reuses its correction workspace and rescales Brownian corrections that would overflow or underflow in intermediate products.
- Complex winding uses signed ray crossings with an exact comparison when product rounding obscures the sign. The numerical registry accepts measured closed polygons, and the argument-principle sampler evaluates each contour vertex once.
- Prepared ray–triangle intersections support measured coordinate references, finite or unbounded ray intervals, and batches with a shared direction. Near-edge tests use an exact sign check. The quadratic sieve validates its work limits and can return an ordered exact factor pair.
- Prepared Reed–Solomon encoders reuse feedback tables. Encoding, decoding and error counting share the 255-symbol codeword bound, and Rust callers can encode into an existing byte buffer.
- Exact factorial, binomial and ordered-selection counts use balanced integer products for large selections. Homotopy continuation checks degree sums and the Bézout path count before constructing start paths.
- Prepared inverse real FFTs reuse unpacking factors and report separate sample and bin counts. The measured interface preserves the spectrum dimension, including for subnormal samples.
- Linear regression borrows sample columns and preserves finite fits across extreme moment scales. Checked polynomial long division can recover a finite remainder after a product overflows. Adaptive Simpson quadrature reuses endpoint and midpoint samples across subdivisions.
- Matrix access and deserialization check both extents. Matrix transposition can reuse a correctly shaped destination. Nelder–Mead sorts vertex indices and reuses the storage for trial points.
- Persistent homology keeps the earliest filtration value for duplicate simplices, treats signed zeros as the same filtration level, and reuses the boundary reduction workspace.
- Double-double square root retains precision at extreme scales. The numerical registry provides double-double Euclidean norm and complex modulus, with a separate readout of the two components in their original physical units. Complex double-double division and reciprocal scale the intermediate products and preserve representable components.
- CSV and JSON array output use buffered writes. The CSV reader reuses its line storage and parses directly into the result array. Raster writers stream row-major pixels, and 16-bit PGM conversion uses a bounded reusable buffer.
- Binomial masses and beta densities handle support endpoints. A bounded binomial recurrence preserves low-count probabilities when the success probability is tiny. Binomial sampling uses discrete inversion with one uniform draw for nondegenerate parameters. Uniform density, CDF and sampling also support finite endpoints whose difference exceeds the floating-point range. Seeded binomial sample sequences change with the new algorithm.
- Improved extreme-tail probabilities for statistical tests and used thin SVD for ridge and OLS regression. Kendall correlation handles ties in O(n log n). Bootstrap procedures check sample and output sizes before allocating.
- Extended symbolic Mellin integrals, limits and infinite sums. Large-argument complex Faddeeva, error and Fresnel functions use scaled computations to avoid intermediate overflow.
- Poisson, heat and wave grid solvers accept per-face Dirichlet and Neumann conditions. Pure-Neumann Poisson checks the source–flux balance and returns a zero-mean solution. Structured kernel diagnostics identify rejected shapes, domains and convergence criteria.
- Three-dimensional shallow-water Euler and AMR runs accept a barotropic state in `(h, hu, hv, hw)`. The kernel manifest describes boundary conditions and tensor ranks for typed hosts.

- Fixed-step ODEs can return only the final state. Diagonal ETDRK4 shares a real-weight stepper across fixed and adaptive routes. Kerr–Newman trajectories can sample by stride and initialize either supported future-directed energy branch within the ergoregion.
- Finite-field homology computes all requested degrees in one pass over the boundary ranks. Binary elimination packs row bits, and odd-prime arithmetic covers the full unsigned word range. Permutation rotations accept signed arbitrary-precision exponents.
- Seeded Ising runs support online observables, and Potts sweeps reuse neighborhood weight tables. Complex Lambert W distinguishes both sides of its negative real cut, and inverse erfc retains accuracy for subnormal probabilities.
- The transverse magnetic FDTD solver reuses its field buffers and can return only the requested final electric field. The physical operator returns the field and elapsed time in their respective dimensions.
- The two-dimensional Schrödinger solver prepares splitting phases for repeated initial waves. Its measured terminal operator supports signed time increments and returns complex amplitudes with inverse-length dimensions.
- Cartesian finite-volume Euler and shallow-water solvers reuse their Runge–Kutta workspace and support the two-wave HLLE flux. Prepared triangular P1/P2 meshes reuse geometry and quadrature for independent diffusion–reaction problems. Extreme finite loads use scaled conjugate gradients.
- Prepared Cartesian Q1 Poisson solves reuse the sine-basis stiffness spectrum on uniform two- and three-dimensional grids. Prepared two-dimensional fast multipole evaluations reuse source geometry across strength vectors. Chain dynamics exposes measured accelerations in the original body order.
- AMR computes independent flux lines in parallel on large patches and validates refinement extents and final times. Prepared ideal-gas Riemann solutions reuse the wave pattern across queries, including vacuum and extreme common pressure and density scales.
- Checked explicit Euler, midpoint and Heun end at the requested time in either direction. Periodic one-dimensional DG advection can reuse a fixed mesh and velocity across time steps.
- Sparse MD neighbour searches scale with occupied cells. The one-dimensional method of lines avoids copying the full field at each right-hand-side evaluation. Truncated holonomic series can compute the first coefficients of large powers directly.
- Channel flow reuses the pressure system between steps. Closed fluid regions require compatible net flux and receive a pressure gauge. Rigid molecular dynamics precomputes its rotor composition before updating the state, and the cell-flow exponential uses closed forms for diagonal matrices and simple shear.
- Three-dimensional Laplace boundary-element solves with fixed boundary kinds can reuse the dense collocation operator and its factorization across boundary values.
- Boundary-element queries reuse already validated mesh and coefficient data and reject mixing interior and exterior solutions. Cut-cell callbacks retain the signs of subnormal distances. One-dimensional advection leaves its retained mesh intact if the measured output cannot be represented.
- The parallel CPU device splits large one-row matrix products over output columns while preserving the reduction order. Real algebraic field embeddings generate successive powers in the quotient basis and eliminate independent rows in large exact solves.

### Build and compatibility

- Source builds require Rust 1.100 or newer. The default Rust artifact is an rlib. Build C libraries explicitly with their generated headers. C archives include English and Chinese instructions for linking directly from the extracted directory. [Build guide](docs/BUILDING.md)
- C batch computation accepts uninitialized output buffers. The `ffi` build profile supports converting unwinding Rust panics to C error codes at guarded entry points.
- The manifest API version stays at 2. Inspect the manifest hash before sharing serialized contracts across independently built kernels. [Kernel embedding](docs/KERNEL.md)

### 内核与嵌入

- 原生对象绑定从同一参数列表生成签名与运行时参数个数。稠密值检查无需分配遍历缓冲区。
- 科学计算内核提供确定性的清单、类型与量纲契约、结构化结果以及类型化数值回调。宿主从实际执行调用的同一算子注册表生成受检查的接口。
- 重写索引和增量匹配减少重复搜索。独立请求使用各自的计算图，结果按输入顺序导入。可选 JIT 和设备后端明确支持的精度与回退行为。
- 工作线程在独立计算图中准备和执行符号批处理，并限制同时保留的任务数。CPU 调度读取当前 Rayon 池的大小，FFT 通过 `kernels::reuse` 共享计划并复用临时空间。
- 符号化简的树窗口复用相交子图的抽取结果。编译后的 CPU 批处理按容量借用临时缓冲区，小任务顺序执行。同尺寸、同方向的并发 FFT 请求共享计划构建。
- 列式计算带逐行读取被选中的来源列。预备的选列计算带可接收共用参考量的带量纲来源列，并保留结果的原始量纲。
- 标准会话共享已安装规则和清单，直至全局单位表变动。会话的自定义单位、项图、原生对象及计算缓存各自独立。
- C ABI 捕获共享的标准清单快照并复用 JSON 编码，每次调用仍返回独立所有权的字符串。
- 二进制清单、插件消息和计算记录采用 bincode-next 3.2 原生编码。签名和清单保持确定性编码，校验清单时直接读取现有的算子表和单位表，逐段将规范编码传给哈希器。兼容性测试将支持的旧线格式与原生编码进行比较。JSON 与 CBOR 继续提供自描述元数据。

### 数值与符号算法

- 横磁 FDTD 求解器复用场缓冲区，并可直接返回指定步数的最终电场。物理算子分别保留电场和经过时间的量纲。
- 二维薛定谔求解器为多个初始波函数复用分步相位。带量纲的末态算子支持带符号时间步，并返回逆长度量纲的复振幅。
- 笛卡尔有限体积欧拉和浅水求解器复用 Runge–Kutta 工作区，并支持双波 HLLE 通量。预备的 P1/P2 三角网格为独立扩散反应问题复用几何与求积数据。极端有限载荷使用缩放后的共轭梯度。
- 预备笛卡尔 Q1 Poisson 求解在二维和三维均匀网格上复用正弦基刚度谱。二维快速多极子法在多组源强度间复用源点几何。链式动力学按输入物体顺序返回带量纲的加速度。
- AMR 在较大的网格块上并行计算独立通量线，并校验细化范围和终止时间。预备理想气体 Riemann 求解在多次查询中复用波系，支持真空区及密度和压力的极端同比缩放。
- 受检的显式 Euler、中点法及 Heun 方法可沿任意时间方向到达请求终点。周期一维 DG 平流可在多个时间步之间复用固定网格与速度。
- 稀疏分子动力学近邻搜索的网格存储随占用单元数增加。一维线法在求值右端项时复用积分器状态。`holonomic_series` 可直接求高次幂的前若干系数。
- 通道流在各时间步间复用压力系统。封闭流体区域须满足净通量相容性，并固定压力基准。刚体分子动力学在更新状态前预计算转子组合，晶胞流的矩阵指数对对角矩阵和简单剪切使用闭式解。
- 边界类型固定时，三维 Laplace 边界元求解器可对不同边界值复用稠密配点算子与分解。
- 边界元求值复用已验证的网格与系数，并拒绝内域解和外域解混用。切割单元回调保留次正规距离的符号。一维平流无法恢复带量纲的输出时保留原网格状态。
- 并行 CPU 设备可按输出列划分大型单行矩阵乘积，并保持每个元素的累加顺序。实代数域嵌入在商代数的基上逐次生成原始元的幂，在大型精确求解中并行消去独立行。
- 新增多源 BFS、强连通分量缩点图、割点与桥，以及三种时间序列延迟统计量。大型 Welch 任务并行计算固定帧块并依序合并。Dormand–Prince、共轭梯度和重启 GMRES 复用迭代工作区。
- 优化后的稠密矩阵乘法直接使用行主序输入和输出。大型乘积使用当前 Rayon 线程池。
- 欧几里得范数、距离和向量方向计算会在极端范围缩放中间平方项。距离计算省去差向量的分配，Rust 调用方可将单位向量写入自备缓冲区。
- 向量反射直接写入输出缓冲区，并在中间校正项溢出时保留可表示的分量。Arnoldi 缩放基向量构造并复用投影存储。原生稠密特征求解器将矩阵向量乘积直接写入迭代工作区。
- 双线性插值在单点查询与预备网格中共用缩放后的区间坐标，包括端点有限但间距超过 Float 范围的情况。
- 预备的 Newton–CG、多项式原差分和几何网格传递保留可复用的计算状态。隐式 Gauss–Legendre 在各配点分别计算雅可比。
- 支持实阶 Bessel K 的对数求值，以及完整无符号模数范围内的精确 Chebyshev 余数。精确稀疏消元复用行存储，并为较大系统建立主元索引。

- 增加单边 Welch 功率谱密度、重叠归一化逆 STFT 和多相 FIR 重采样。增加受检查的拓扑排序，并让 Rust 和符号图算子共用非递归强连通分量实现。
- 增加解析信号、Hilbert 变换、zoom FFT、谱相干度和小波包。DCT、DST 采用 FFT 算法，分数阶 Fourier 变换支持复用谱基。
- 增加重叠保留卷积、自相关、FIR 滤波器设计、频率响应、二阶节滤波、多相重采样，以及支持系数阈值处理的平稳小波变换。
- 增加批量分位数和经验分布函数、Benjamini–Hochberg、Benjamini–Yekutieli 与 Holm 多重检验校正、多元正态批量对数密度，以及可输出相关矩阵的原子协方差批量更新。
- 增加含负边权的最短路径、弱连通分量、最大流及最小割、列插值分解。内核清单声明输入、结构化输出和量纲契约。
- 稀疏共轭梯度对极端尺度的右端项做缩放。固定步长时滞求解拒绝无法前进的时间网格，波动方程离散格式先合并波速和步长再计算平方。

- 增加加权经验分位数和累积分布、含不可逆底数的复合模数离散对数，以及带溢出检查的精确整数多项式乘法。改进 LU 分解复用、特征值计算的缩放、最小二乘投影和并行分形关联计数。
- Rust LOBPCG API 支持带正定质量算子的无矩阵广义特征值问题，并软锁定已收敛的向量。
- 固定步长 ODE 求解器在调用用户函数前检查时间网格和轨迹容量。事件定位复用状态缓冲区。
- Meijer G 的单极点公式在共同尺度上计算幂和指数，保留次正规结果。可交换 Platen 求解器复用校正工作区，并缩放会在中间乘积中上溢或下溢的布朗校正项。
- 复数绕数使用带符号的射线穿越计数，浮点乘积舍入无法区分符号时使用精确比较。数值注册表支持带量纲的闭合多边形，辩值原理的采样在每个轮廓顶点只求值一次。
- 预备的射线与三角形求交支持带量纲的坐标参考值、有限或无上界的射线区间，以及方向相同的批量射线。贴近边缘时检查精确符号。二次筛校验工作预算，并可返回有序的精确因子对。
- 预备的 Reed–Solomon 编码器复用反馈表。编码、解码和错误计数共用 255 个符号的码字上限，Rust 调用方可将码字写入自备字节缓冲区。
- 精确阶乘、二项式系数和有序选法数在大规模选择时使用平衡整数乘积。同伦延拓在构造起始路径之前检查次数和与 Bézout 路径数。
- 预备的逆实数 FFT 复用解包系数，分别报告样本数和频率分量数。带量纲的接口保留频谱量纲，支持次正规样本。
- 线性回归借用样本列，并在中间矩尺度极端时保留可表示的拟合结果。受检查的多项式长除法可在中间乘积溢出后恢复有限余式。自适应 Simpson 求积在各次细分之间复用端点和中点样本。
- 矩阵访问和反序列化检查两个轴的尺寸。矩阵转置可复用形状匹配的目标矩阵。Nelder–Mead 对顶点索引排序，并复用试探点的缓冲区。
- 持久同调对重复单形保留最早的过滤值，将带符号的零视为同一过滤层，并复用边界约化工作区。
- 双双精度平方根在极端尺度下保留精度。数值注册表提供双双精度的欧几里得模长和复数模长，支持按原量纲分别读取高、低分量。复双双精度除法与倒数缩放中间乘积，保留可表示的分量。
- CSV 和 JSON 数组采用缓冲写入。CSV 读取复用行缓冲区，并直接将解析值写入结果数组。光栅文件按行主序流式写入像素，16 位 PGM 转换复用有容量上限的缓冲区。
- 二项分布概率质量与 Beta 密度正确处理支撑集端点。成功概率极小时，有界递推可保持小计数的二项概率。非退化二项分布使用一次均匀随机数反演采样。即使有限端点之差超出浮点范围，均匀分布的密度、累积分布与采样仍可求值。新算法会改变给定种子的二项采样序列。
- 改进统计检验的极端尾部概率，岭回归与普通最小二乘采用薄 SVD。Kendall 相关系数以 O(n log n) 处理并列值。自助法在分配内存前检查样本与输出大小。
- 扩展符号 Mellin 积分、极限和无穷级数。复数 Faddeeva 函数、误差函数和菲涅耳积分在大参数下采用缩放计算，避免中间溢出。
- 泊松、热传导和波动网格求解器支持按面指定 Dirichlet 或 Neumann 条件。全 Neumann 泊松问题检查源项与通量平衡并返回零均值解。结构化内核诊断指出不符合要求的形状、定义域或收敛准则。
- 三维浅水 Euler 和自适应网格细化求解器接收 `(h, hu, hv, hw)` 形式的正压状态。内核清单声明宿主所需的边界条件和张量秩约束。

- 定步长 ODE 可以只返回末态。对角 ETDRK4 的定步与自适应路径共用实权重步进器。Kerr–Newman 轨道可按步距采样，并在能层内选择相容的未来指向能量分支。
- 有限域同调通过一次边界秩计算得到各次数结果。二元消元打包行位，奇素数运算覆盖完整无符号机器字范围。置换循环支持有符号任意精度指数。
- 带种子的 Ising 模拟支持在线观测，Potts 扫描复用邻域权重表。复数 Lambert W 区分负实割线两侧，反余误差函数保留次正规概率的尾部精度。

### 构建与兼容性

- 源码构建需要 Rust 1.100 或更新版本。默认生成 Rust rlib。C 库按独立流程构建并携带生成的头文件。C 归档提供中英文说明，可从解压目录直接链接库和运行示例。[构建指南](docs/BUILDING.zh.md)
- C 批量计算允许传入未初始化的输出缓冲区。`ffi` 构建配置允许受保护的入口将可展开的 Rust panic 转为 C 错误码。
- 清单 API 版本保持为 2。不同构建的内核共享序列化契约前，应核对清单哈希。[内核嵌入](docs/KERNEL.md)

## [3.0.0-rc.1]

The release candidate establishes the graph rewriting engine, numerical ABI and standard operator registry on the 3.0 branch. The release notes above cover changes to this architecture through the beta series.

该候选版本在 3.0 分支建立图重写引擎、数值 ABI 与标准算子注册表。上方发布说明覆盖该架构在 beta 系列中的变更。
