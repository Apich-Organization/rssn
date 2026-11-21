# 数学覆盖与输入约定

[English](MATH_COVERAGE.md) · [内核嵌入](KERNEL.zh.md)

按数学适用范围、数据表示和误差契约选择算法。运行时清单列出当前构建所公开的算子。Rust 函数和动态算子可以使用不同接口调用同一实现。

## 数值输入

受检查的 ABI 对无效形状和不受支持的定义域返回错误或算子规定的 IEEE 值。覆盖空张量、形状和数据长度不一致、非有限数值以及超出范围的控制参数。底层 Rust 访问器要求索引有效，宿主回调遵循自身 Rust panic 契约。

| 领域 | 输入和结果约定 |
|---|---|
| 稠密线性代数 | 受检查的接口区分非有限数据、形状不匹配和奇异分解。`Lu::solve_into` 复用大小精确匹配的输出存储，要求右端项有限。秩亏问题可选 `pinv`、`lstsq_svd`、`null_space` 或 `orth`，并按问题设定奇异值阈值。 |
| 矩阵函数 | 对称矩阵的 `sqrtm`、`logm` 和 `powm_psd` 使用谱分解，对角矩阵直接逐项计算。半正定矩阵即使奇异，也可以求平方根和正指数幂。对称矩阵的对数和负指数幂要求特征值均为正。一般实矩阵使用 Denman–Beavers 或逆缩放平方方法。 |
| 统计 | 各样本统计和检验遵循各自的缺失值与空样本约定。`percentile` 接收 `[0,100]` 内的百分位。生存函数可保留极小尾概率，避免从一减去已被舍入到一的 CDF。加权经验统计要求有限的非负权重和正的总权重，零权重不计质量。 |
| 实数与复数特殊函数 | 极点与分支切口采用对应函数规定的值。无效数值可能产生 NaN。复 `ln_gamma` 使用解析分支，Lambert W 接收分支索引，复 Bessel J 和 I 接收实阶。菲涅耳和复误差函数在浮点范围边缘使用缩放尾项。`bessel_i0e` 消去实数修正 Bessel I₀ 的主导指数因子，Kaiser 窗使用该形式保持归一化有限。 |
| 精确整数算法 | 64 位模运算使用更宽的精确中间整数。`is_prime` 对 `u64` 是确定性的，因式分解采用 Pollard–Brent。`checked_lcm` 报告溢出。离散对数先去除公因子，再在互素情形采用大步小步法，其表空间可能与约化模数的平方根成正比。符号规则集提供任意精度算法。 |
| 数值求积 | 可积端点奇点适合双曲正切–双曲正弦求积。自适应 Gauss–Kronrod 提供误差估计，`integrate_gk_breakpoints` 可在已知内部间断点或奇点分段。提供相位导数时，可用 `levin_adaptive` 处理振荡积分。Simpson 和 Romberg 有工作量上限。 |
| ODE、SDE 和求根 | 求解器有各自的求值或步数限制，并返回类型化错误。刚性常微分方程可用隐式方法。SDE 可用漂移隐式方法，并须选择匹配的噪声结构和 Itô 或 Stratonovich 约定。Brent 需要两端异号或某端点恰为根。容差和诊断见具体签名。 |

## 结构化及精确算法

- `Complex::sqrt` 返回主平方根，在负实轴上根据虚部符号选择分支。它保留可表示的微小分量，输入的模长溢出时仍可得到有限根。带量纲的原生算子对输入量纲取平方根。
- `AffineIfs` 一次预备无量纲仿射映射和非负权重，再按指定无符号种子重复生成轨迹。预备采样器支持总和溢出的有限权重，并在轨迹超出有限范围时返回类型化错误。
- `IrfftPlan` 针对固定的实数输出长度复用解包系数。受检查的原生算子返回与半频谱同量纲的样本。两个长度查询分别返回实数样本数和频率分量数。
- `Matrix<T>::transpose_into` 复用调用方的可变目标矩阵，并在写入前检查转置后的形状。可移植 Float 和精确矩阵算子返回各自拥有的转置结果。
- 多项式长除法逐步移动系数窗口，避免在每步移位剩余输入。`Polynomial::try_long_division` 跳过除式前导零，拒绝零除式，并在中间乘积溢出时通过融合运算保留可表示的系数。
- 稀疏符号展开在一次转换中复用已展开的非多项式子项，并保留各生成元首次出现的顺序。遍历在堆上保存待处理项，跨越非多项式运算时仍检查深度。多项式乘法拒绝超出整数表示范围的指数。
- `Categorical` 为有限非负相对权重保留可复用的采样表。抽样返回原顺序的零基索引，并跳过零权重项。所有权重相等或仅一个类别有正权重时无需建表。
- `permgroup::power_exact` 按各循环的长度分别对指数取模，计算置换的有符号整数次幂。原生算子接收零基编号，要求所有映射构成双射、置换度数不超过 4,096、指数绝对值不超过 65,536 位。
- `ode::solve_fixed_final` 与完整定步长求解器共用 Euler、Heun 和经典 RK4 的阶段计算，只保留当前状态和工作缓冲区。原生算子接收已存在的右端项程序，积分前检查求值预算，并恢复末态各分量的物理量纲。
- 对角 ETDRK4 的定步和自适应步均使用实权重，复用阶段缓冲区，并在参数接近零或刚度极大且为负时选用稳定公式。带量纲的原生接口接收逆时间量纲的对角项和非线性右端项程序。
- `ChainComplex::betti_mod_p` 对每个边界矩阵只消元一次，返回各次数的 Betti 数。二元矩阵采用位打包的行异或。奇素数消元接受带符号矩阵元素，模数覆盖 64 位无符号素数。语言接口在消元前检查模数。
- 同调、稀疏精确消元、机器字因数分解与符号素性检查共用同一组确定性 `u64` 见证。普通乘积将溢出时，机器字模乘使用精确的宽乘积。
- `crt_general` 为模数为正且相容的同余方程组返回最小非负解与精确周期。即使模数有公因子，周期也是它们的最小公倍数。适用于两两互素模数的 `crt` 返回最小非负解。
- 当当前 CPU 池有多个工作线程，且两两互素的模数规模较大、位数分布较均衡时，`crt` 使用乘积树与余数树。每次求逆使用原始模数，合并按输入顺序执行。
- `inverse_erfc` 使用对数残差细化极端尾部，包括次正规概率。复数 Lambert W 按虚部符号区分负实轴割线的两侧，并在主分支保留微小分量。
- Kerr–Newman 赤道测地线初始化可选择未来指向的能量分支，并支持位于静止界线内、外视界外的相容轨道。采样求解器保留首末点，以及指定步距上的位置。原生物理参考系将两列直角坐标恢复为长度。
- `kerr::trace_shadow` 用径向方程判定光线是否被捕获，并与完整射线追踪共用径向力和势。执行池可同时追踪独立方向，按图像平面角度顺序返回。`num_native_relativity_trace_shadow` 和 `num_native_relativity_shadow_analytic` 都按所选物理参考系恢复图像平面上的长度。解析算子对曲线两半各采样 2 至 2,048 点。
- 二维周期 Ising 求解器复用棋盘缓冲区。正温下，若多格点系统有一条边的长度为一或二，则采用带停留概率的逐点更新。在线测量保留最终自旋，按固定步距统计热化后的可观测量。物理参考系用玻尔兹曼常数将开尔文温度转换为交换能量单位。全部磁化样本为零时没有 Binder 值。
- Potts 热浴扫描每个邻域只统计一次，并按颜色计数复用权重。平移后的权重在极端耦合下保持条件概率有限。Rust 调用方可以复用工作区，原生模型也会为重复调用保留工作区。只读格点算子返回各颜色的无量纲概率。
- Gross–Pitaevskii 定步求解器复用 FFT 计划与缓冲区，步数或时间增量为零时返回归一化初态。原生输出为按行存储的逆面积量纲概率密度。
- Q4 平面应力装配共用同一个受检查的单元刚度算法。求积前平移单元坐标，并可按公共尺度归一化。装配按网格单元顺序处理，并在求解前检查连接关系和载荷。
- `PreparedTriangles` 为相互独立的扩散反应问题复用三角形 P1/P2 网格的几何、自由度布局和求积数据。求解会检查实际残差，极端载荷使用带缩放的共轭梯度。原生计划和问题共用同一物理参考系。
- `PreparedCartesianPoisson` 在齐次 Dirichlet 边界的均匀正方形或立方体上复用可分离的 Q1 刚度谱。二维和三维 Poisson 入口共用同一正弦基求解。预备接口每个求积点求值一次新源项，原生接口为按 x 轴最快顺序排列的节点数据提供物理参考量。
- `PreparedLaplaceBem3d` 为固定的边界类型保留稠密配点算子及分解。每次求解按面片顺序接收新值。纯 Neumann 问题检查净外向通量，并将面积加权平均势固定为零。原生解可按网格参考量纲求取内部势。
- `SchrodingerPlan2D` 为不同初始波函数复用带符号的时间步及势能、动能相位。带量纲的末态算子保留逆长度复振幅和带符号时间。势能均匀时，多步末态请求可只需一对傅里叶变换。采样接口保留第 1、11、21 步等结果。
- `run_fdtd_final_field` 使用两个场缓冲区推进原有的横磁模型，并返回指定时间步的电场。受检查的物理接口允许单单元轴，按时间步 `L/(2c)` 换算源频率，并返回经过时间。采样接口保留第 1、6、11 步等结果。
- 笛卡尔网格上的欧拉和浅水求解器在各 Runge–Kutta 阶段复用填充状态和残差。不同工作线程数下，每个网格单元按同一轴顺序累加通量。`RiemannSolver::Hlle` 提供采用 Einfeldt 波速估计的双波通量。
- AMR 在较大的网格块上通过当前 CPU 线程池并行计算独立通量线。时间子循环、通量寄存与回流修正保持原有顺序。`ExactRiemann` 一次准备理想气体波系，支持单点或按输入顺序批量采样。带量纲接口接收速度 `x/t`，分别返回密度、速度及压力。
- 大盒子中的分子动力学近邻搜索对占用的网格单元建立稀疏索引，并保持与稠密遍历相同的顺序。
- `physics_md::morse_radial` 求取平移后的 Morse 势能及带符号的径向力，覆盖接近平衡和零间距的情形。`num_native_md_morse_radial` 按输入的参考量制返回带能量和力量纲的结果。
- `md_engine::Rdf` 按各帧自身的理想气体参考量汇总粒子对计数。`num_native_engine_rdf_merge_d1`、`num_native_engine_rdf_merge_d2` 和 `num_native_engine_rdf_merge_d3` 合并截断半径、分箱和参考量制相同的累积器。
- `physics_sm::advance_advection_diffusion_3d_modes` 直接推进可复用的傅里叶场，按 x 方向变化最快的顺序排列，并采用不归一化的正向 FFT 约定。`num_grid_sm_advance_advection_diffusion_3d_modes` 保留输入的场量参考。时间推进使用谱法前向欧拉。
- 一维线法直接从 ODE 状态读取差分模板。`method_of_lines_1d` 对固定边界数据只求值一次。`method_of_lines_1d_timed` 在求解器求值时刻采样 Dirichlet 场值及 Neumann 外法向导数，并恢复末时刻的 Dirichlet 场值。带类型的原生算子为两端程序传入带量纲的时间，并保留输出场的量纲。
- `SseXxz` 在算子串非空时按与当前状态无关的次数执行闭环更新，在算子串为空时采样全部自由自旋。相同种子在此版本内可复现，但与旧版本的轨迹不同。原生扫掠的预算涵盖算子和格点遍历、起点采样及环上的步数。推进开始后若耗尽预算，对应可变句柄将失效。
- 结构化连续介质多重网格在 V、W 和全多重网格循环中共用工作区。相对残差在极端有限值处使用缩放计算，并保留普通路径的累加顺序。全多重网格在校正循环前保留初始残差和全网格阶段的残差。
- 符号稀疏矩阵构造和 `sp_nnz` 在逻辑行数很大时对实际占用行建立索引。结果保留逻辑尺寸与规范行序，重复条目仍按精确值合并和相消。
- 符号 `det` 对精确矩阵条目调用共用的有理数行列式。满足工作量要求时，各待消元行在当前 CPU 池上独立执行。小型符号行列式仍按子式展开缓存中间结果，计算过程无需除法。
- 线性常微分方程的特征环求解器在有理系数微分算子的复合与右除中复用右操作数的高阶导数。满足工作量要求时，各系数的导数可在当前 CPU 池上独立计算。右除按算子次数递减更新精确余项。
- `ode_asymptotic` 将形式级数的系数转换到首个被舍弃项为止，并将该项的绝对值乘以解的前因子模长。若已用完所提供的项，则用最后一项估计误差。`ode_borel` 检查所提供级数中的指数间隔，并在有阶数上限的 Padé 方程组中复用定点系数。独立的 Padé 阶数区间可在当前 CPU 池上并行计算，最后按请求的阶数依次比较求和值。
- `ode_frobenius` 的 Newton 多边形解在不同的共轭域嵌入之间共享精确系数。报告复用已构造的两个解的前因子，数值求值按相同顺序选择分支。
- Gosper 的 `antidifference` 比较最高次系数，得出至多两个候选多项式次数，包括两个最高次项均相消时的特殊次数。较大的移位公因子搜索可在当前 CPU 池上并行执行，并选取最早匹配的移位。
- 符号有限和的备用路径拟合一阶递推，并在三个额外的下标验证所得表达式。有理函数按额外参数进行精确比较。含有仍可求值的特殊函数时使用数值检查。
- `taylor` 与 `laurent` 的直接展开路径递推生成含阶乘、反正弦及二项式的精确权重。级数组合利用参数的最低幂限制所需权重个数。`asinh` 与 `atanh` 在参数常数项为零时使用有理 Taylor 权重。
- `root_branch` 在正整数根次数下先对精确整数分支下标取模，再构造相位。`sqrt_branch` 与二次根使用同一分支奇偶性。直接复数求值保留平方根割线的两侧，并对极端幅值进行缩放。
- 有限群同构与自同构搜索从各生成元前缀已验证的映射继续扩展。每次试探检查新增的右乘关系，并在返回时恢复原前缀。首个同构按候选像在运算表中的顺序选取。
- Todd–Coxeter 的指数与群阶请求检查完整的转移表，并返回仍有效的陪集数。`fp_group` 先验证关系式，再把生成作用的群阶与该陪集数比较，最后构造群阶至多为 720 的 Cayley 表。
- 置换群 Cayley 表构造在每行的乘积计算中复用置换缓冲区。较大的表可在当前 CPU 池上独立构造各行。行与列均按群元素的排序排列。
- F4 按 sugar 次数处理临界对，同一次数内按插入顺序处理。在发现矩阵的行与列时检查上限。较大的批次可并行准备 Macaulay 矩阵行，并按同一顺序交给消元。
- F5 在不同工作线程数下保留同一签名顺序和约化顺序。生成元归一化和每轮较大的预处理约化可在当前 CPU 池上独立执行。签名循环检查初始基的大小，并在当前多项式的尾部改变时检查系数大小。
- `apart` 在有理数域上分解有理函数。它使用多项式中国剩余定理计算各个分母因子幂对应的分子，并按因子顺序与幂次升序输出非零项。较大的分解可在当前 CPU 池上独立计算各因子块。
- `root_of(k, c0, ..., cn)` 的第一个参数是从零开始的下标。系数为精确有理数时，算法分离相互重叠的有理包围区间，先对不同的实根排序，再列出非实根。较大多项式的独立因子可使用当前 CPU 池。全为实根的因子在实根隔离后完成。
- `resultant` 与 `discriminant` 对有理标量系数使用精确多项式余式。系数为独立表达式时使用 Sylvester 矩阵和无分数消元。较大矩阵每次选定主元后，可在当前 CPU 池上独立更新各行。两个算子均拒绝依赖消去变量的系数。
- 实数域量词消去中的代数域范数在各插值点共用精确有理结式。约定为 Sylvester 结式，输入系数按次数升序排列。特化多项式为零时保留该零插值点，并使用 Newton 差商重构范数。较大的独立样本集合可使用当前 CPU 池。代数域嵌入在商代数的基上逐次生成原始元的幂。求解较大的精确线性系统时，工作线程并行更新互不依赖的行。
- `square_free_factors(poly, var)` 返回最高次项系数以及按重数分组的首一有理因子。各因子两两互素，列表按重数升序排列。非零常数返回空因子列表。
- `solve` 与 `solve_complex` 的六次和七次 Galois 预解式分析在候选群之间复用定点置换单项式，并在指数选择之间复用陪集成员。独立计算可使用当前 CPU 池。陪集代表、定点乘法及加法的顺序在不同工作线程数下保持一致。
- 分圆 Gauss 周期在每个导子内复用单位的定点相位，并在不同指数选择间复用各特征的已排序陪集。较大的相位表可使用当前 CPU 池。实周期只需计算余弦。
- Galois 定点根精化只准备一次多项式及其导数的系数，并按原顺序精化独立初值。`groebner` 算子默认在有理数域上使用字典序。指定 `grlex` 或 `grevlex` 可选择其他次序。
- 符号验证直接判断精确数值残差是否为零，无需转为浮点数。近似常数和未化简残差的采样值各自使用相应的数值容差。串行检查在首个反驳结论的有限样本处结束。较大的检查可以并行计算同一组采样点。
- `verify_limit(f, x, a, L, [dir])` 共用微积分的极限求解器，分别比较带符号的无穷大和有限结果，并返回布尔值。可选的 `plus` 或 `minus` 选择单侧，数值采样只在指定侧进行。
- 公共子表达式消除将声明为绑定范围的操作数保留在原作用域内。`cse_cost` 精确计算改写前展开表达式树的操作数，以及改写后各临时定义和各输出的操作数。临时量引用无需额外操作。
- 切割单元根据穿过零点的边计算单形开口比例。有限值的差溢出时，缩放后的边比例仍可保留结果。单元递归复用固定大小的局部采样存储。
- 仿射与射影同伦求解器共用次数检查和起始点顺序。仿射路径可在当前 CPU 线程池中独立跟踪。各适配器在构造起始点之前执行其文档中规定的次数和路径预算。
- `simple_linear_regression` 借用成对样本，并在需要时缩放中心矩。缩放路径遇到跨度极大的二进制指数时，使用精确累加器。原生斜率量纲为响应量纲除以预测量纲，截距保留响应量纲。
- `LinearRegression::value_and_gradient` 共用一次残差计算，并从数值复用池借用工作空间。物理接口中，损失量纲为响应量纲的平方，截距导数量纲为响应量纲，各斜率导数量纲为响应与预测量纲的乘积。
- `Dd::sqrt` 在输入接近 Float 范围边界时仍保留低分量。`num_native_dd_hypot` 和 `num_native_cdd_norm_dd` 保留原量纲，`num_native_dd_parts` 分别返回两个分量。`num_native_cdd_recip` 对复数商进行缩放，并取输入量纲的倒数。
- `persistent_homology` 在最早的过滤值处理每个单形，将带符号的零视为同一层。带量纲的数值接口分别返回出生值、死亡值和永续区间布尔掩码张量。
- 实值 Meijer G 的单极点公式在对数尺度合并 `z^b` 和 `exp(-z)`，保留可表示的次正规结果。
- 二次筛在工作预算不足以处理一个筛块时，仍检查可精确整除的小因子。`qs_factor_pair` 验证乘积后按升序返回两个非平凡精确因子，每个因子都可能仍可分解。
- `PreparedTriangle` 为单条射线和共享方向的批量射线复用三角形几何。带量纲的数值接口保留射线参数和交点的坐标量纲，法线无量纲。批量结果带有命中掩码。
- `Bsr::matmul_dense` 将块稀疏矩阵乘以按行主序存储的多列右端项，无需把稀疏矩阵展开为稠密矩阵。每个输出元素对精确点积做一次舍入。`num_linear_bsr_matmul_dense` 算子返回两个操作数参考尺度及量纲的乘积。
- `sparse_native` 使用一次稀疏 LU 或 Cholesky 分解求解多列右端项。原生批量接口返回解数组，其量纲为右端项量纲除以算子量纲。自动接口还会报告所选的分解方法。
- 一般 AMG 的 V/W 循环共用可复用工作区。原生批量预条件算子对每行从零初值应用该层次，并保留右端项除以算子的量纲。
- 张量列车内积根据键秩选择两阶段收缩顺序，并复用有上限的工作区。原生结果量纲为两个张量列车量纲的乘积。
- `SimplexMap` 在自适应积分和无符号测度之间复用仿射几何。行列式在舍入前对主元乘积进行缩放，因此极大与极小轴尺度并存时仍可得到可表示的测度。
- `real_roots::find_roots_in` 在有限闭区间中搜索实根，并保留恰好落在端点的根。`num_native_real_roots_find_roots_in` 要求端点、容差与正坐标参考量具有相同量纲，返回同量纲的有序实根。
- `brent_with_bracket` 返回根及最终的有序夹逼区间。带量纲的 `num_brent_with_bracket` 保留回调输入的坐标量纲，并报告实际求值次数。找到精确零点时，两端点均为该点。
- `levin_condition_numbers` 通过一次共用递推计算各个可用阶数的无量纲条件数。无穷值表示条件数未定义或无界。原生算子支持同量纲级数项以及 `U`、`T` 和 `V` 三种变换。
- `sim::symplectic::integrate_sampled` 保存可分离辛积分的初态、指定步距上的状态和末态。原生算子保留精确步数，并恢复每个位置与速度分量的物理量纲。
- `ricci_curvature` 在 Ricci 张量和标量之间复用有限差分采样和中心点的度规分解。原生算子检查坐标与度规参考量，并恢复两种结果各自的量纲。
- `SpectralGrid::helmholtz_plan` 为重复或批量右端项准备周期傅里叶符号。原生批量算子接收场量除以长度平方的源项及长度负二次方的移位量。相容的奇异傅里叶模态取零振幅。
- `solve_burgers_nd_with_boundaries` 用共享的 Lax–Friedrichs 界面通量更新单元平均值，各轴可选周期或零梯度外流边界。全周期网格在舍入误差范围内保持单元平均值之和。使用速度参考系的检查型算子逐轴校验 CFL 条件。
- `direct_field_into` 将带软化的 N 体场写入调用方提供的存储。原生算子的场量纲为源强量纲乘以长度量纲的 `1 - 维数` 次方。
- `solve_wave_dirichlet_sampled` 从初始位移和位移速率出发求解二维或三维齐次 Dirichlet 波动。所有样本从初态起即满足零边界，并包含末态。原生算子支持有符号的时间步，保留位移量纲。
- `solve_heat_equation_2d_cn_adi_fixed_boundary` 在 Crank–Nicolson ADI 的两个半步中均保持初始边缘值。二维的两个入口都跨时间步复用三对角分解。
- `NsSolver::divergence_summary` 计算流体单元中的最大绝对散度及体积加权均方根散度。原生算子保留散度的时间倒数量纲、流体面积或体积和有符号的外向通量。`divergence_into` 将逐单元结果写入可复用的切片。
- MAC 压力投影为每个未指定压力的独立分量选择一个压力基准。矩阵与右端项使用相同的基准行。密度与压力模板相同时可复用分解。
- 稀疏重排序包括逆 Cuthill–McKee、近似最小度及嵌套剖分。`ordering_metrics` 扫描一次稀疏结构即可计算半带宽和下三角轮廓。`num_native_ordering_metrics` 接收 CSR 邻接偏移、邻居编号和零基置换 `perm[new] = old`。每条边按无向边处理，重复条目不影响结果。
- `pchip_many` 和 `akima_many` 可在相邻查询落在同一分段时复用系数。它们保留查询顺序并支持重复坐标。带类型的适配器要求查询坐标与结点量纲相同，并保留纵坐标量纲。
- `tensordot_views` 可借用张量视图，只为原布局无法展平的操作数重排数据。`num_linear_tensordot` 使用同一缩并实现，保留非缩并轴的顺序及两个参考量的乘积。
- `GeometricMg` 在红黑平滑时只访问当前颜色的网格点。V/W 循环与完全多重网格共用残差工作区，并在普通平方和超出浮点范围时比较缩放后的范数。
- `MultilinearGrid` 一次验证直线网格的坐标轴与扁平采样布局，重复查询无需分配。会话算子 `num_native_multilinear_grid_new` 和 `num_native_multilinear_grid_eval` 支持各轴独立的物理量纲，并保留采样值量纲。
- `kernels::dense` 和 `kernels::schur` 提供实 Schur、广义 Schur、Sylvester 与 Lyapunov 求解。`kernels::lowrank` 提供列主元 QR、随机 SVD、LOBPCG 及 Anderson 加速。`lobpcg` 求解标准对称特征值问题。`lobpcg_generalized` 通过回调接收刚度算子与正定质量算子，返回关于质量内积正交归一的特征向量，并可跳过已收敛特征对的搜索更新，同时保留它们参与 Rayleigh–Ritz 投影。
- `kernels::sparse_exact` 和 `sparse_symbolic` 提供在有理域、有限域和符号域上的消元、使用 Hadamard 界的 CRT 行列式、带有理重建的 Dixon 提升，以及标量 Wiedemann 方法。模最小多项式使用随机投影。有理输入的秩用有理消元。需要 Smith 或 Hermite 标准形时，使用对应的稠密精确算法。
- `kernels::holonomic` 与 `holonomic_*` 支持湮没算子、在加法、乘法、乘方和微分下的闭包、系数递推以及常点处的精确零测试。这里的系数路径使用单变量及原点处的有理 Taylor 数据。定和使用独立的望远镜规则。闭包湮没算子满足所需方程，但阶数可能高于函数的全局最小算子。
- 数值有限域包括素域 `GF(p)` 和 `GF(256)`。`Gf256Multiplier` 跨字节批量复用乘积，支持缩放和异或乘加，并与 Reed–Solomon 使用相同的多项式。符号 `gf_*`、`gfx_*` 和 `gfp_factor` 提供扩域和符号因式分解。
- 适当区域中的整数阶 Bessel J 和 Y 使用递推。复数渐近展开在零点附近或相对参数很大的阶数下可能损失相对精度。对残差或单个分量有误差要求时，应在目标参数区域与独立参考值比较。
- `carlson_rg` 计算第二类 Carlson 对称椭圆积分。`jacobi_amplitude` 保留累积的周期信息。`jacobi_elliptic_batch` 对全部输入复用同参数的变换计划，返回按 `[sn, cn, dn]` 排列的各行。
- `fresnel_sc_interval` 在整个区间上计算正弦和余弦 Fresnel 积分，支持无限端点、相邻浮点端点和远场积分尾项。`EllipticPiPlan` 在计算不同振幅时复用第三类完全椭圆积分。
- `cartesian_to_spherical` 和 `spherical_to_cartesian` 针对固定三元向量，Rust 内核无需堆分配。正向转换保留近轴的微小极角，以及半径溢出时的方向。清单中的固定形状算子使用统一参考长度和弧度。
- `SylvesterPlan` 在不同右端项之间复用固定系数矩阵的 Schur 分解。`SylvesterPlan::for_lyapunov` 仅做一次分解，准备转置系数组。原生解的量纲和参考尺度分别为右端项除以系数的对应量。
- `bessel_in_log_abs` 计算整数阶第一类变形 Bessel 函数绝对值的自然对数，可覆盖原值或缩放值超出浮点范围的情况。
- `EwaldReciprocalPlan` 为固定周期晶胞、屏蔽参数和截断准备共轭波矢对。粒子构型可以复用同一计划，并各自计算相位、能量、力与维里张量。并行结果按固定波矢块顺序合并。
- `SphDensityPlan` 为固定坐标和最大平滑长度保留邻居与距离。密度计算可更换质量、平滑长度或内核。求和包括粒子自身和周期或反射镜像。原生密度量纲为质量除以长度的空间维数次幂。
- `Complex::recip` 与复数除法在极端分量下进行缩放，保留可表示的结果。倒数算子的清单契约将输出量纲设为输入量纲的倒数。
- `Yee::plane_wave_sample` 读取正入射 TF/SF 辅助线上的电场、磁场及各自时间。辅助线只推进因果可达的部分，直到抵达预留终点。原生适配器恢复两种场各自的单位。
- `ChebSeries::from_values` 用已有样本构造第二类 Chebyshev 插值模型。大型变换使用共享的 FFT 计划和缓冲区。自适应细分复用重合结点采样，要求同一坐标对应相同函数值。
- `pic::AssignmentPlan` 在固定粒子坐标下，为周期线性或三次沉积和插值共用权重。带量纲的网格密度为电荷除以长度，粒子采样保留网格场量纲。
- `LinearSplitPlan` 为固定实势、质量和带符号时间步保留 Strang 或四阶 Yoshida 相位。它在一到三维周期网格上对多组波函数复用 FFT 计划。带量纲适配器保留波函数振幅的体积负二分之一次幂量纲。
- `L2Grid` 为有限非递减坐标轴保留梯形权重的平方根。范数舍入值溢出或下溢时，预备归一化仍可使用舍入前的尺度。测度为零时，带量纲接口不返回归一化采样。
- `JumpDistribution::fill_marks` 按顺序向调用方缓冲区写入跳跃标记，消耗的随机数与逐次抽样一致。固定步长跳扩散复用同一标记缓冲区，并限制预期和实际累计的跳跃工作量。带量纲批量接口每行返回一个样本。
- `Sobol::point_at` 按索引读取单位立方体中的点，`Sobol::seek` 改变后续位置。拟蒙特卡洛使用固定分块归约，并在每块内递推 Sobol 坐标。
- `LevinPlan` 为固定相位、频率和区间复用带权结点。可传入振幅函数或同序采样，清单适配器支持带量纲振幅。缓慢变化的相位使用 Gauss–Legendre 求积。
- `eval_value_and_jacobian_terms` 和 `values_and_jacobian_at_terms` 在同一复数盒或点上对值与导数复用幂和单项式乘积。复数盒计算采用向外舍入。带量纲 Jacobian 的单位为残差除以坐标单位。
- `sum_terms_to_infinity` 对迭代器给出的部分和做加速，最多消耗请求的项数。采样前缀适配器报告带量纲估计、可用的启发式绝对误差及实际项数。`sum_range` 包含最终整数，即使该整数为 `i64::MAX`。
- CSC 矩阵与向量乘积只遍历有效元素，即使各列留有预分配空间。支持正向和转置方向，也可复用调用方的输出缓冲区。原生目标保留自身参考量，仅在完成检查与缩放后写入。
- `SimplexResidualEstimator` 保存连接关系和几何权重，可反复估计有限元状态。常系数和变系数共用这些几何数据。在 `d` 维中，指标量纲为场量纲乘以长度的 `(d + 2)/2` 次幂，再除以时间。
- `matrix_norms` 从有效 CSC 元素一次计算 Frobenius 范数、诱导 1-范数和诱导无穷范数。缩放累加处理输入或中间结果尺度极端时的有限范数。
- `kruskal_wallis` 和 `mann_whitney_u` 共用合并秩并处理并列。各组样本数可以不同。合并秩退化时返回定义域错误。
- `bell_batch` 通过有限数值表按输入顺序返回 Bell 数，原生结果用 `None` 标记单个索引溢出。`PrimeSqrtPlan` 为同一素数模数复用素性与约化数据，求较小的模平方根。
- `Multivector::mul_central` 根据几何代数的伪标量平方，乘上由标量和伪标量组成的中心元素。原生算子将两个操作数的量纲与参考尺度相乘。求逆使用缩放，输入尺度极端时仍可保留可表示的逆元。
- `huber_location_with_scale` 使用调用方传入的正尺度并复用残差缓冲区，求解截断分数方程。`num_robust_huber_location_with_scale` 要求样本与尺度共用量纲。
- `Dual::powi_with_second_derivative` 返回归一化整数幂及沿输入方向的导数。中间幂溢出或下溢时仍可恢复有限的导数。
- `hessian_diagonal` 和 `hessian_diagonal_with_steps` 只计算各坐标方向的二阶导数。`laplacian` 按坐标顺序对同一组对角项求和。带量纲回调接口的结果为输出量纲除以坐标量纲的平方。
- `NewtonCgWorkspace` 为重复求解保留 Newton 和 CG 缓冲区。原生计算带计划拥有选定目标的快照，分别保留坐标、参数和输出的物理参考量。
- `PolynomialAntidifference` 一次编译升次幂系数，并用差商保留邻近端点之间的微小差值。原生计划接收带量纲的原点和正步长，结果保留系数量纲。
- `indefinite_sum_with_budget` 与 `indefinite_product_with_budget` 在执行回调前检查格点和延拓所需的求值次数，并拒绝非有限端点与超出可表示范围的格点下标。带量纲求和保留回调输出量纲，乘积使用除以参考量后的各正项。
- `randomized_svd` 在正交化之前缩放数值尺度极端的投影列。原生稠密适配器将正向与转置的矩阵向量乘积写入分解过程提供的缓冲区，奇异值参考量保留输入矩阵量纲。
- `hyp2f1_cplx` 对收敛的正实参数计算 Gauss 超几何函数在 1 处的值，避免 Gamma 中间量溢出。[NIST DLMF §15.4.20](https://dlmf.nist.gov/15.4.E20)
- `ln_bessel_k` 返回实阶第二类修正 Bessel 函数的对数。积分围绕峰值平移并缩放，原值上溢或下溢时仍可求得对数。[NIST DLMF §10.32.9](https://dlmf.nist.gov/10.32.E9)
- 精确 Chebyshev 函数按次数做二进制倍增，并保留周期点的有理值。`chebyshev_tu_mod` 同时计算两类函数的余数，支持完整的 `u64` 次数和模数范围。[NIST DLMF §18.5](https://dlmf.nist.gov/18.5)
- `GeometricMg` 支持相邻层之间的全权重限制和齐次 Dirichlet 边界下的多线性延拓。V 周期、W 周期和全多重网格共用传递缓冲区。原生适配器在切换网格层时保留输入数组的量纲和参考尺度。
- 具有一至三个配点的 Gauss–Legendre 方法使用各配点阶段的雅可比构造 Newton 矩阵。`GaussLegendreWorkspace` 保留阶段缓冲区，按阶段导数的尺度判断收敛。
- `holonomic` 模块的 `Expr::series` 用系数递推计算精确截断幂级数。正的最低次数会平移级数，有时会使整个幂落在所需截断范围之外。
- `deciding_terms` 根据精确的递推首系数和保守的整数根界，确定需要多少个初始项。搜索包含预算给定的末端。
- `fd_jacobian_dense_par` 为每个任务复用扰动状态，并连续写入该任务的各列。每次求值的输出从零初始化，各线程数下每列逐位一致。
- `fd_jacobian_colored_par` 在不同颜色之间复用扰动状态。对颜色较多的稀疏模式，只保留 CSR 输出用到的采样项。每次右端项求值的输出从零初始化。
- `weibull_quantiles` 对下尾概率复用参数验证与尺度，`weibull_inverse_survival` 直接接收上尾概率。带量纲接口保留尺度单位，并对无穷结果报告定义域错误。

## 信号处理

- `DwtPlan` 为等长信号的单层小波变换保留滤波器和边界布局。Rust 正逆变换支持调用方缓冲区，清单适配器在两个频带中保留样本量纲。

- `kernels::signal` 使用可复用的 FFT 计划实现正交归一的 DCT-II、DST-II 及其逆变换。既有的 `signal::fft` 和 `signal::ifft`、FFT 卷积与 `kernels::fft` 共用计划和工作缓冲区。多维实 FFT 在计算其余各轴前利用 Hermitian 对称性缩减一轴。规模较大的独立 FFT 批次和 STFT 帧可并行计算，结果顺序固定。
- `welch_psd` 对完整加窗帧取平均，得到单边功率谱密度，窗口可以重叠。`welch_cross_spectrum` 还计算 `conj(X) * Y` 的互谱密度和幅值平方相干度。任一信号在某个频点功率为零时，该频点的相干度为零。`inverse_stft` 用完整的短时傅里叶频谱重建各帧覆盖的区间，某个样本缺少非零窗口贡献时返回错误。`upfirdn` 在重采样时应用给定的 FIR 系数，无需生成上采样信号。
- 带类型的 `num_signal_welch_psd` 和 `num_signal_welch_cross_spectrum` 为各返回字段记录各自的量纲。频率沿用采样率单位，谱密度按单位频率计，相干度无量纲。
- 大规模 Welch 计算按固定的帧数分块，并按块顺序合并频谱。`num_signal_autocovariance` 先将样本中心化，再把各延迟的结果除以总样本数。`num_signal_autocorrelation` 用零延迟自协方差归一化，`num_signal_partial_autocorrelation` 使用 Yule–Walker 方程。三者均返回延迟 `0..=max_lag`。自协方差的量纲为样本量纲的平方，两种相关系数均无量纲。
- `analytic_signal`、`hilbert` 和 `analytic_envelope` 使用周期边界。`zoom_fft` 通过线性调频卷积计算等间隔频点，频率单位为每输入样本的周期数。中等规模的变换只要输入长度、频点数和频率步长相同，就可以复用卷积计划，起始频率可以不同。带类型的 `num_signal_zoom_fft` 算子接收和返回复数向量。
- `kernels::wavelets::wavelet_packet` 对低频和高频分支继续分解，得到等宽的叶子频带，每个低频子带排列在相应高频子带之前。`inverse_wavelet_packet` 从叶子矩阵重建信号。对长度相同的信号计算多个阶数的离散分数傅里叶变换时，`FrftPlan` 可以复用特征基。
- 一维实 FFT 共用奇偶长度处理和旋转因子。宽跨步轴通过连续缓冲区分批并行计算。`kernels::transforms` 在二次幂长度上保留原有标量基二 FFT 的舍入次序。
- `fftconv::convolve_overlap_save` 用独立的重叠保留块计算完整线性卷积。`OverlapSavePlan` 为多段信号复用滤波器频谱，并在变换前缩放极端振幅。原生批量算子接收带量纲的展平采样及零基行偏移，保留空行，结果量纲为采样与滤波系数的量纲乘积。`autocorrelate` 返回非负延迟的线性自相关，零延迟等于平方和。
- `firwin_lowpass` 和 `firwin_bandpass` 设计 Kaiser 窗 FIR，截止频率以 Nyquist 频率为 1。`frequency_response` 则使用每采样周期数（Nyquist 为 0.5）。`sos_filter` 接收二阶节系数，`resample_poly` 用居中的多相 FIR 完成有理倍率重采样。
- `swt` 和 `iswt` 提供周期平稳小波变换，各频带与输入等长，支持奇数长度。不同层级复用近似系数缓冲区，不同线程数保持相同的系数顺序。`wavelet_soft_threshold` 和 `wavelet_hard_threshold` 对细节频带去噪，阈值与样本量纲一致。

## 加权图

`kernels::graph::Graph` 存储有向边，支持自环和平行边。`Graph::from_edges` 先验证整批输入，并为每个邻接表预分配空间。无向加权图的每条边用两条权值相等、方向相反的有向边表示。

| 算法 | 输入与结果 |
|---|---|
| Dijkstra | 计算非负权图的单源最短路。 |
| Bellman–Ford | 支持带负权的单源最短路，拒绝从源点可达的负环。 |
| DAG 最短路 | 在有向无环图上以线性时间计算单源最短路，支持负权。 |
| Floyd–Warshall | 返回按行排列的全点对最短路矩阵。负对角元表明图中存在负环。 |
| 拓扑排序、强连通和弱连通分量 | 忽略边权。强连通分量按拓扑顺序编号，弱连通分量按各自的最小顶点编号。符号强连通分量算子复用同一迭代遍历。 |
| 多源 BFS 与缩点图 | BFS 返回到最近源点的最少边数，不可达顶点无距离。缩点图按拓扑顺序为强连通分量编号，并去除分量间的重复边。 |
| 割点与桥 | 每条无向边以一对端点表示，支持平行边。返回割点与桥的端点。 |
| Prim 最小生成森林 | 输入为以对称有向边表示的无向图，结果保留各分量及孤立顶点。 |
| PageRank | 平行边分别计入出度，出度为零的顶点将其权重均匀分配。 |
| Dinic 最大流与最小割 | 容量必须有限且非负，源点与汇点必须不同。返回流量值和最小割源侧顶点的布尔掩码。 |

`num_native_graph_*` 句柄保留边权的物理量纲。批量构造接收两列端点矩阵和一阶权值张量。路径记录的 `distances` 使用边权的一贯单位，`previous` 为可选的前驱顶点编号。不可达顶点的距离为正无穷。最大流的 `value` 使用容量的一贯单位，`source_side` 为无量纲掩码。原生加权图最多支持 256 个顶点和 65,536 条有向边，最大流还要求 `V² E` 不超过 5,000 万。可调用的算子名称及类型签名见清单。

## 低秩矩阵

`kernels::lowrank::interpolative_decomposition` 通过截断列主元 QR 分解，用原矩阵的部分列和无量纲系数矩阵近似输入矩阵。输入秩不足时，实际秩可能低于指定的上限。`num_linear_interpolative_decomposition` 接收原生矩阵句柄，返回所选列的索引以及按原始列顺序排列的插值系数。

## 空间维数与网格

`solve_transport_sampled_nd` 保留初态、指定步距上的快照以及未落在步距上的终态。其带类型的适配器 `num_grid_dimensional_solve_transport_sampled_nd` 分别限制保存的快照数和推进全部时间步的工作量。

热传导、波动、Laplace 和 Poisson 符号方程每个系列有一个算子。问题本身决定空间维数，并保留 `_1d`、`_2d` 和 `_3d` 兼容别名。`heat_fd`、`wave_fd`、`poisson_fd` 和 `advect_fd` 等有限差分算子从输入形状推导单位盒上的空间秩。边界版本允许 Dirichlet 和 Neumann 面。全 Neumann Poisson 检查源项与外向通量的相容性，返回零均值解。线方法的辅助接口用于一维网格，具体方法规定其支持的边界条件。

`pde_poisson_integral(g, radius, dimension)` 将角向边界数据延拓至正半径的圆盘或球体内部。常量边界数据保持原值，圆盘上单个整数频率的正弦或余弦有闭式解。一般数据用定积分表示，积分变量与数据中的参数分别绑定。

`sim::dimensional` 使用 x 轴变化最快的平坦字段，轴大小为 `[nx]`、`[nx, ny]` 或 `[nx, ny, nz]`。共享输运接口支持平流、扩散及两者组合。Poisson 版本提供松弛及多重网格。其他算法包括波动、Euler 方程、几何物理场及专用场景。场景清单描述空间维数、字段布局和方法选项。`tests/sim/dimension_coverage.rs` 和 `tests/kernels/dimension_coverage.rs` 使用解析解、低维退化、对称性、收敛阶或守恒律验证结果。

标量有限体积的定步长与自适应求解器对四种格式共用同一重构和通量计算，在各时间步之间复用工作缓冲区，并按原网格顺序调用回调。自适应求解器在固定网格上控制时间离散误差。经典 Lax–Friedrichs 格式的数值黏性依赖时间步长，因此可能需要放宽误差容限。参见 `tests/fvm_facade_reuse.rs` 和 `tests/pde_adaptive.rs`。

## 选择可调用接口

`rules::numeric` 下的字面项算子接收数字、列表和表示方法的符号。输入未求值或不满足前置条件时，可以保留为符号请求。可用 `Session::apply_diag` 查看所报告的拒绝原因。表达式算子覆盖稳健统计、傅里叶与小波变换、稀疏迭代、精确整数卷积、张量分解、采样函数、多项式求根和网格求解。接口中的索引从零计数。

类型化数值 ABI 还公开受检查的记录、张量和会话所属模型，支持复用 Cholesky、LU 与多重网格状态、张量列车代数、复数变换数据、受检查的同伦认证，以及具备专用求解控制的 ODE 和 SDE 对象。可从 `KernelManifest` 或 Aevia 的 `aevia rssn` 输出查询确切名称与签名。需要 Rust 回调或借用宿主存储的原生方法也可直接调用。[内核嵌入契约](KERNEL.zh.md) · [数值回调](NUMERICAL_CALLBACKS.zh.md)

## 统计

`num_descriptive_quantiles` 使用第 7 型插值。`num_robust_quantiles` 支持九种 Hyndman–Fan 定义，`num_robust_weighted_quantiles` 求加权经验分布函数的逆。每批查询共享排序后的样本，按查询顺序返回结果并保留样本量纲。`num_robust_ecdf_batch` 要求样本与查询量纲相同，返回概率。

`bootstrap_ci` 返回百分位置信区间。`bootstrap_basic_ci` 以原始样本的统计量为中心反射百分位边界。两者对有限且归一化的观测数据进行带种子的重采样，返回下界和上界。中位数、MAD 和 IQR 的大量重采样批次在当前工作线程之间并行计算，并复用重采样缓冲区。Aevia 的 `Rssn_Numeric` 提供这两个算子。

`num_descriptive_histogram` 根据样本生成等宽区间。`num_descriptive_histogram_edges` 接收有限、严格递增且与样本量纲一致的边界，可在不同数据集中使用同一组区间。内部边界计入右侧区间，最高边界计入最后一个区间。显式边界必须覆盖所有观测值。样本较大且区间数适中时，各工作线程分别计数再合并，最终整数计数与线程数无关。

`num_stats_ext_adjust_p_bh`、`num_stats_ext_adjust_p_by` 和 `num_stats_ext_adjust_p_holm` 按原顺序返回校正后的无量纲 p 值。Benjamini–Hochberg 在检验独立或零假设满足正回归依赖时控制错误发现率，Benjamini–Yekutieli 适用于任意依赖，Holm 控制族错误率。

`num_native_online_covariance_push_batch` 接收与累加器量纲相符的行主序观测矩阵。整个批次成功后才提交，并始终按固定顺序合并大型输入的分块。`num_native_online_covariance_correlation` 返回无量纲 Pearson 相关矩阵。任一列为常量或观测不足两条时，返回 `None`。`num_distributions_multivariate_normal_log_pdfs` 对所有查询行复用一次协方差分解，输入须按同一参考尺度归一化为无量纲值。

## 验证位置

| 领域 | 代表测试 |
|---|---|
| 稠密线性代数、矩阵函数与零空间 | `tests/math_linear_algebra.rs`、`tests/kernels/dense.rs` |
| 复函数与实函数尾部 | `tests/math_special_complex.rs`、`tests/math_special_edges.rs`、`kernels::fresnel` 的回归测试 |
| 素性、因子、模根、离散对数及广义同余方程组 | `tests/math_number_theory.rs`、`tests/symbolic_modular_exact.rs` |
| 统计、求积、插值与逆问题 | `tests/math_numeric_edges.rs`、`tests/math_robustness.rs`、`tests/kernels/robust.rs` |
| 批量分位数、p 值校正、协方差和正态密度 | `tests/statistics_batch.rs`、`tests/native_discrete.rs`、`tests/source_policy_math.rs` |
| ABI 定义域、形状及工作量控制 | `tests/abi_fuzz.rs`、`tests/analysis_native.rs`、`tests/source_policy_math.rs` |
| 积分、极限、求和、级数和符号求解 | `tests/math_symbolic_classics.rs`、`tests/symbolic_battery.rs` |
| Gosper 反差分与移位公因子搜索顺序 | `src/rules/calculus/gosper.rs`、`tests/symbolic_gosper.rs` |
| 有限和递推拟合与额外下标验证 | `src/rules/calculus/recurrence.rs` |
| 直接 Taylor 和 Laurent 级数组合 | `src/rules/calculus/powerseries.rs`、`tests/symbolic_power_series.rs` |
| 复数根的分支与整数周期 | `src/rules/complex/branches.rs`、`tests/symbolic_complex_branches.rs` |
| 有限群同构见证与自同构 | `src/rules/discrete/group_theory.rs`、`tests/symbolic_group_isomorphisms.rs` |
| 群呈示与 Todd–Coxeter 陪集作用 | `src/rules/discrete/perm_groups.rs`、`tests/symbolic_presented_groups.rs` |
| 并行构造 Cayley 表 | `src/rules/discrete/perm_groups.rs`、`tests/symbolic_cayley_tables.rs` |
| F4 临界对调度与精确矩阵准备 | `src/rules/poly/f4.rs`、`tests/symbolic_groebner_f4.rs` |
| F5 签名约化、资源上限与系数域 | `src/rules/poly/f5.rs`、`tests/f4_f5.rs` |
| 精确部分分式与重复分母幂 | `src/rules/poly/apart.rs`、`tests/math_symbolic_classics.rs` |
| 精确实根下标与隔离 | `src/rules/poly/rootof.rs`、`tests/symbolic_rootof_order.rs` |
| 精确有理与含参数的结式 | `src/rules/poly/univariate.rs`、`src/rules/poly/algebra.rs`、`tests/symbolic_resultant.rs` |
| 实代数域的精确范数 | `src/rules/solve/qe_alg.rs`、`src/rules/solve/qe_field.rs` |
| 有理无平方因子分解与符号参数类型 | `src/rules/poly/algebra.rs`、`tests/symbolic_square_free_factors.rs` |
| 六次和七次 Galois 预解式与根的标号 | `src/rules/poly/galois_high.rs`、`src/rules/solve/mod.rs` |
| 定点 Newton 根精化与符号 Gröbner 默认次序 | `src/rules/poly/galois.rs`、`src/rules/poly/mod.rs`、`tests/f4_f5.rs` |
| 符号结论验证的精确残差、数值采样和单侧极限 | `tests/symbolic_verify_residual.rs`、`tests/symbolic_verify_limit.rs` |
| 有理系数微分算子的复合与特征环分解 | `src/rules/ode/ore.rs`、`tests/symbolic_ode_ore.rs` |
| 常微分方程形式解的截断与 Borel–Padé 求值 | `src/rules/ode/formal.rs`、`tests/symbolic_formal_ode.rs` |
| Borel–Padé 阶数选择与 Stokes 中值 | `src/rules/ode/resum.rs`、`tests/symbolic_borel_stokes.rs` |
| Frobenius 展开报告与共轭分支选择 | `src/rules/ode/frobenius.rs`、`tests/symbolic_frobenius_embeddings.rs` |
| 共享子项的符号多项式展开与精确指数限制 | `tests/symbolic_poly_repr.rs` |
| 科学计算布局与不变量 | `tests/native_*.rs`、`tests/dimensioned_coverage.rs` |
| 信号分析、FIR 重采样与有向图遍历 | `tests/numerical_signal_graph.rs` |
| 实数张量 FFT、余弦与正弦变换、Hilbert 变换、局部 FFT、相干度、小波包、平稳小波、重叠保留卷积、数字滤波 | `tests/signal_transforms.rs` |
| 带权路径、最小割、图同构与原生图算子类型 | `tests/kernels/graph.rs`、`tests/native_discrete.rs` |
| 低秩近似及插值系数的量纲 | `src/kernels/lowrank.rs`、`tests/native_linear.rs` |

符号测试区分可检查的精确值与保持为符号请求的已记录样例。具体期望以测试用例为准。数值结果还受精度、条件数和相应误差阈值影响。
