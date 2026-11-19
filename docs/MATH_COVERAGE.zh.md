# 数学覆盖与输入约定

[English](MATH_COVERAGE.md) · [内核嵌入](KERNEL.zh.md)

按数学适用范围、数据表示和误差契约选择算法。运行时清单列出当前构建所公开的算子。Rust 函数和动态算子可以使用不同接口调用同一实现。

## 数值输入

受检查的 ABI 对无效形状和不受支持的定义域返回错误或算子规定的 IEEE 值。覆盖空张量、形状和数据长度不一致、非有限数值以及超出范围的控制参数。底层 Rust 访问器要求索引有效，宿主回调遵循自身 Rust panic 契约。

| 领域 | 输入和结果约定 |
|---|---|
| 稠密线性代数 | 受检查的接口区分非有限数据、形状不匹配和奇异分解。`Lu::solve_into` 复用大小精确匹配的输出存储，要求右端项有限。秩亏问题可选 `pinv`、`lstsq_svd`、`null_space` 或 `orth`，并按问题设定奇异值阈值。 |
| 矩阵函数 | 对称矩阵的 `sqrtm` 和 `logm` 使用谱分解。对称半正定矩阵的平方根可以奇异。一般实矩阵使用 Denman–Beavers 或逆缩放平方方法，无法计算实主值时返回错误。 |
| 统计 | 各样本统计和检验遵循各自的缺失值与空样本约定。`percentile` 接收 `[0,100]` 内的百分位。生存函数可保留极小尾概率，避免从一减去已被舍入到一的 CDF。加权经验统计要求有限的非负权重和正的总权重，零权重不计质量。 |
| 实数与复数特殊函数 | 极点与分支切口采用对应函数规定的值。无效数值可能产生 NaN。复 `ln_gamma` 使用解析分支，Lambert W 接收分支索引，复 Bessel J 和 I 接收实阶。菲涅耳和复误差函数在浮点范围边缘使用缩放尾项。`bessel_i0e` 消去实数修正 Bessel I₀ 的主导指数因子，Kaiser 窗使用该形式保持归一化有限。 |
| 精确整数算法 | 64 位模运算使用更宽的精确中间整数。`is_prime` 对 `u64` 是确定性的，因式分解采用 Pollard–Brent。`checked_lcm` 报告溢出。离散对数先去除公因子，再在互素情形采用大步小步法，其表空间可能与约化模数的平方根成正比。符号规则集提供任意精度算法。 |
| 数值求积 | 可积端点奇点适合双曲正切–双曲正弦求积。平滑被积函数可用带误差估计的自适应 Gauss–Kronrod。提供相位导数时，可用 `levin_adaptive` 处理振荡积分。Simpson 和 Romberg 有工作量上限。 |
| ODE、SDE 和求根 | 求解器有各自的求值或步数限制，并返回类型化错误。刚性常微分方程可用隐式方法。SDE 可用漂移隐式方法，并须选择匹配的噪声结构和 Itô 或 Stratonovich 约定。Brent 需要两端异号或某端点恰为根。容差和诊断见具体签名。 |

## 结构化及精确算法

- `kernels::dense` 和 `kernels::schur` 提供实 Schur、广义 Schur、Sylvester 与 Lyapunov 求解。`kernels::lowrank` 提供列主元 QR、随机 SVD、LOBPCG 及 Anderson 加速。`lobpcg` 求解标准对称特征值问题。`lobpcg_generalized` 通过回调接收刚度算子与正定质量算子，返回关于质量内积正交归一的特征向量，并可跳过已收敛特征对的搜索更新，同时保留它们参与 Rayleigh–Ritz 投影。
- `kernels::sparse_exact` 和 `sparse_symbolic` 提供在有理域、有限域和符号域上的消元、使用 Hadamard 界的 CRT 行列式、带有理重建的 Dixon 提升，以及标量 Wiedemann 方法。模最小多项式使用随机投影。有理输入的秩用有理消元。需要 Smith 或 Hermite 标准形时，使用对应的稠密精确算法。
- `kernels::holonomic` 与 `holonomic_*` 支持湮没算子、在加法、乘法、乘方和微分下的闭包、系数递推以及常点处的精确零测试。这里的系数路径使用单变量及原点处的有理 Taylor 数据。定和使用独立的望远镜规则。闭包湮没算子满足所需方程，但阶数可能高于函数的全局最小算子。
- 数值有限域包括素域 `GF(p)` 和 `GF(256)`。符号 `gf_*`、`gfx_*` 和 `gfp_factor` 提供扩域和符号因式分解。构造和运算检查各自适用的域不变量。
- 适当区域中的整数阶 Bessel J 和 Y 使用递推。复数渐近展开在零点附近或相对参数很大的阶数下可能损失相对精度。对残差或单个分量有误差要求时，应在目标参数区域与独立参考值比较。


## 信号处理

- `kernels::signal` 使用可复用的 FFT 计划实现正交归一的 DCT-II、DST-II 及其逆变换。既有的 `signal::fft` 和 `signal::ifft`、FFT 卷积与 `kernels::fft` 共用计划和工作缓冲区。多维实 FFT 在计算其余各轴前利用 Hermitian 对称性缩减一轴。规模较大的独立 FFT 批次和 STFT 帧可并行计算，结果顺序固定。
- `welch_psd` 对完整加窗帧取平均，得到单边功率谱密度，窗口可以重叠。`welch_cross_spectrum` 还计算 `conj(X) * Y` 的互谱密度和幅值平方相干度。任一信号在某个频点功率为零时，该频点的相干度为零。`inverse_stft` 用完整的短时傅里叶频谱重建各帧覆盖的区间，某个样本缺少非零窗口贡献时返回错误。`upfirdn` 在重采样时应用给定的 FIR 系数，无需生成上采样信号。
- 带类型的 `num_signal_welch_psd` 和 `num_signal_welch_cross_spectrum` 为各返回字段记录各自的量纲。频率沿用采样率单位，谱密度按单位频率计，相干度无量纲。
- `analytic_signal`、`hilbert` 和 `analytic_envelope` 使用周期边界。`zoom_fft` 通过线性调频卷积计算等间隔频点，频率单位为每输入样本的周期数。带类型的 `num_signal_zoom_fft` 算子接收和返回复数向量。
- `kernels::wavelets::wavelet_packet` 对低频和高频分支继续分解，得到等宽的叶子频带，每个低频子带排列在相应高频子带之前。`inverse_wavelet_packet` 从叶子矩阵重建信号。对长度相同的信号计算多个阶数的离散分数傅里叶变换时，`FrftPlan` 可以复用特征基。
- 一维实 FFT 共用奇偶长度处理和旋转因子。宽跨步轴通过连续缓冲区分批并行计算。`kernels::transforms` 在二次幂长度上保留原有标量基二 FFT 的舍入次序。
- `fftconv::convolve_overlap_save` 用独立的重叠保留块计算完整线性卷积，`OverlapSavePlan` 可为多段信号复用滤波器频谱。`autocorrelate` 返回非负延迟的线性自相关，零延迟等于平方和。
- `firwin_lowpass` 和 `firwin_bandpass` 设计 Kaiser 窗 FIR，截止频率以 Nyquist 频率为 1。`frequency_response` 则使用每采样周期数（Nyquist 为 0.5）。`sos_filter` 接收二阶节系数，`resample_poly` 用居中的多相 FIR 完成有理倍率重采样。
- `swt` 和 `iswt` 提供周期平稳小波变换，各频带与输入等长，支持奇数长度。`wavelet_soft_threshold` 和 `wavelet_hard_threshold` 对细节频带去噪，阈值与样本量纲一致。

## 加权图

`kernels::graph::Graph` 存储有向边，支持自环和平行边。`Graph::from_edges` 先验证整批输入，并为每个邻接表预分配空间。无向加权图的每条边用两条权值相等、方向相反的有向边表示。

| 算法 | 输入与结果 |
|---|---|
| Dijkstra | 计算非负权图的单源最短路。 |
| Bellman–Ford | 支持带负权的单源最短路，拒绝从源点可达的负环。 |
| DAG 最短路 | 在有向无环图上以线性时间计算单源最短路，支持负权。 |
| Floyd–Warshall | 返回按行排列的全点对最短路矩阵。负对角元表明图中存在负环。 |
| 拓扑排序、强连通和弱连通分量 | 忽略边权。强连通分量按拓扑顺序编号，弱连通分量按各自的最小顶点编号。符号强连通分量算子复用同一迭代遍历。 |
| Prim 最小生成森林 | 输入为以对称有向边表示的无向图，结果保留各分量及孤立顶点。 |
| PageRank | 平行边分别计入出度，出度为零的顶点将其权重均匀分配。 |
| Dinic 最大流与最小割 | 容量必须有限且非负，源点与汇点必须不同。返回流量值和最小割源侧顶点的布尔掩码。 |

`num_native_graph_*` 句柄保留边权的物理量纲。批量构造接收两列端点矩阵和一阶权值张量。路径记录的 `distances` 使用边权的一贯单位，`previous` 为可选的前驱顶点编号。不可达顶点的距离为正无穷。最大流的 `value` 使用容量的一贯单位，`source_side` 为无量纲掩码。原生加权图最多支持 256 个顶点和 65,536 条有向边，最大流还要求 `V² E` 不超过 5,000 万。可调用的算子名称及类型签名见清单。

## 低秩矩阵

`kernels::lowrank::interpolative_decomposition` 通过截断列主元 QR 分解，用原矩阵的部分列和无量纲系数矩阵近似输入矩阵。输入秩不足时，实际秩可能低于指定的上限。`num_linear_interpolative_decomposition` 接收原生矩阵句柄，返回所选列的索引以及按原始列顺序排列的插值系数。

## 空间维数与网格

热传导、波动、Laplace 和 Poisson 符号方程每个系列有一个算子。问题本身决定空间维数，并保留 `_1d`、`_2d` 和 `_3d` 兼容别名。`heat_fd`、`wave_fd`、`poisson_fd` 和 `advect_fd` 等有限差分算子从输入形状推导单位盒上的空间秩。边界版本允许 Dirichlet 和 Neumann 面。全 Neumann Poisson 检查源项与外向通量的相容性，返回零均值解。线方法的辅助接口用于一维网格，具体方法规定其支持的边界条件。

`sim::dimensional` 使用 x 轴变化最快的平坦字段，轴大小为 `[nx]`、`[nx, ny]` 或 `[nx, ny, nz]`。共享输运接口支持平流、扩散及两者组合。Poisson 版本提供松弛及多重网格。其他算法包括波动、Euler 方程、几何物理场及专用场景。场景清单描述空间维数、字段布局和方法选项。`tests/sim/dimension_coverage.rs` 和 `tests/kernels/dimension_coverage.rs` 使用解析解、低维退化、对称性、收敛阶或守恒律验证结果。

## 选择可调用接口

`rules::numeric` 下的字面项算子接收数字、列表和表示方法的符号。输入未求值或不满足前置条件时，可以保留为符号请求。可用 `Session::apply_diag` 查看所报告的拒绝原因。表达式算子覆盖稳健统计、傅里叶与小波变换、稀疏迭代、精确整数卷积、张量分解、采样函数、多项式求根和网格求解。接口中的索引从零计数。

类型化数值 ABI 还公开受检查的记录、张量和会话所属模型，支持复用 Cholesky、LU 与多重网格状态、张量列车代数、复数变换数据、受检查的同伦认证，以及具备专用求解控制的 ODE 和 SDE 对象。可从 `KernelManifest` 或 Aevia 的 `aevia rssn` 输出查询确切名称与签名。需要 Rust 回调或借用宿主存储的原生方法也可直接调用。[内核嵌入契约](KERNEL.zh.md) · [数值回调](NUMERICAL_CALLBACKS.zh.md)

## 统计

`num_descriptive_quantiles` 使用第 7 型插值。`num_robust_quantiles` 支持九种 Hyndman–Fan 定义，`num_robust_weighted_quantiles` 求加权经验分布函数的逆。每批查询共享排序后的样本，按查询顺序返回结果并保留样本量纲。`num_robust_ecdf_batch` 要求样本与查询量纲相同，返回概率。

`num_stats_ext_adjust_p_bh`、`num_stats_ext_adjust_p_by` 和 `num_stats_ext_adjust_p_holm` 按原顺序返回校正后的无量纲 p 值。Benjamini–Hochberg 在检验独立或零假设满足正回归依赖时控制错误发现率，Benjamini–Yekutieli 适用于任意依赖，Holm 控制族错误率。

`num_native_online_covariance_push_batch` 接收与累加器量纲相符的行主序观测矩阵。整个批次成功后才提交，并始终按固定顺序合并大型输入的分块。`num_native_online_covariance_correlation` 返回无量纲 Pearson 相关矩阵。任一列为常量或观测不足两条时，返回 `None`。`num_distributions_multivariate_normal_log_pdfs` 对所有查询行复用一次协方差分解，输入须按同一参考尺度归一化为无量纲值。

## 验证位置

| 领域 | 代表测试 |
|---|---|
| 稠密线性代数、矩阵函数与零空间 | `tests/math_linear_algebra.rs`、`tests/kernels/dense.rs` |
| 复函数与实函数尾部 | `tests/math_special_complex.rs`、`tests/math_special_edges.rs`、`kernels::fresnel` 的回归测试 |
| 素性、因子、模根及离散对数 | `tests/math_number_theory.rs` |
| 统计、求积、插值与逆问题 | `tests/math_numeric_edges.rs`、`tests/math_robustness.rs`、`tests/kernels/robust.rs` |
| 批量分位数、p 值校正、协方差和正态密度 | `tests/statistics_batch.rs`、`tests/native_discrete.rs`、`tests/source_policy_math.rs` |
| ABI 定义域、形状及工作量控制 | `tests/abi_fuzz.rs`、`tests/analysis_native.rs`、`tests/source_policy_math.rs` |
| 积分、极限、求和、级数和符号求解 | `tests/math_symbolic_classics.rs`、`tests/symbolic_battery.rs` |
| 科学计算布局与不变量 | `tests/native_*.rs`、`tests/dimensioned_coverage.rs` |
| 信号分析、FIR 重采样与有向图遍历 | `tests/numerical_signal_graph.rs` |
| 实数张量 FFT、余弦与正弦变换、Hilbert 变换、局部 FFT、相干度、小波包、平稳小波、重叠保留卷积、数字滤波 | `tests/signal_transforms.rs` |
| 带权路径、最小割、图同构与原生图算子类型 | `tests/kernels/graph.rs`、`tests/native_discrete.rs` |
| 低秩近似及插值系数的量纲 | `src/kernels/lowrank.rs`、`tests/native_linear.rs` |

符号测试区分可检查的精确值与保持为符号请求的已记录样例。具体期望以测试用例为准。数值结果还受精度、条件数和相应误差阈值影响。
