# zstd 解码中的相邻标量读取

交付的 `target/rssn.dll` 在 `.rdata` RVA `0x362b830`、`0x362ba40` 和 `0x362bb50` 各有一张完整的 FSE 表，分别与锁定的 zstd 1.5.7 源码中的 `LL_defaultDTable`、`OF_defaultDTable` 和 `ML_defaultDTable` 逐字节一致。PE 导入表没有单独的 zstd DLL。表的指纹能确认数据内容，尚无法把交付 DLL 中的具体函数边界归因到三个解码变体。以下指令和槽位来自可重建的 Windows 与 Linux 产物。 [zstd 1.5.7 解码源码](https://github.com/facebook/zstd/blob/v1.5.7/lib/decompress/zstd_decompress_block.c)

工作区原有的 `target/debug` GCC 16.2.0 对象中，`ZSTD_decompressSequences` 是调度函数，分别调用 `_default` 和 `_bmi2` 解码体。Clang 23.1.5 的 AVX2 release 重建直接把 BMI2 解码体编进三个函数。普通序列解码循环每次读取 LL、ML 和 OF 三个表项。每个表项的四次读取到达寄存器分配后仍独立。在 Linux 的普通序列热循环里，三次 `movzwl` 的结果立即存到帧槽，跨过附加位解码后才用于更新 FSE 状态。三个基值及两个位数也占用独立活跃区间。Windows 重建采用 MinGW 的 Win64 ABI，交付 DLL 的 C 编译配置没有随样本提供。

## 字节范围与代价

`rssn-scalar-load-pack` 查找同一个循环块中来自同一基地址的窄整数读取。每组至少三个字段，地址必须连续且恰好覆盖 4 或 8 字节。首个原始读取的位置必须已经支配最低字段的地址。内存模型检查跨越区间中的写入与别名，过滤原子、volatile、调用等边界。新读取沿用最低字段的对齐保证。每个字段单独投影到原类型。IR 用 `<4 x i8>` 或 `<8 x i8>` 持有字节，避免一个字段的 poison 污染其他字段。LLVM 23 的 x86 后端在本次对照中将两种字节向量读取分别降低成一次普通的 32 位或 64 位读取。 [LLVM load 与 bitcast 语义](https://llvm.org/docs/LangRef.html#load-instruction)

每个基本块须有至少两组合规读取，并且现有活跃度估计达到目标 GP 预算。每个 pack 的读窗口最多 32 条 IR 指令。这个预算用于筛选高压力位置，实际收益还取决于后端的提取指令、访存折叠、调度和寄存器分配。Linux 采用 4 字节组，可以独立读取后半的 32 位字段。Win64 采用 8 字节组，从一个载体投影四个字段。选择来自下表对同一解码器和实际 RSSN 帧的宽度消融。

| LLVM 23.1.5 AVX2 | 静态 Spill / Reload 基线 | 4 字节合并 | 8 字节合并 |
| --- | --- | --- | --- |
| Linux 普通序列 | 56 / 90 | 43 / 85 | 39 / 84 |
| Linux 分段字面量缓冲区 | 88 / 228 | 80 / 203 | 79 / 189 |
| Linux 长偏移序列 | 136 / 362 | 134 / 352 | 127 / 311 |
| Win64 普通序列 | 51 / 74 | 47 / 92 | 46 / 92 |
| Win64 分段字面量缓冲区 | 85 / 179 | 74 / 195 | 65 / 207 |
| Win64 长偏移序列 | 114 / 365 | 114 / 386 | 110 / 363 |

计时使用 RSSN 的 `codec_zstd.parquet` 和 `filters_plugins.h5` 中各四个原样 zstd 帧。大文件对照在八个帧之外还加入 `dense_nested.h5` 以 level 3 编码的内容，压缩工作位于计时区外，解压结果逐字节比较。所有页面帧的解码哈希为 `dfdebb2e332b6efe`，大文件组为 `683487f5beb4982e`。比较的是底层解码器耗时。

| 真实输入，线程 CPU | 4 字节相对基线 | 95% 区间 | 8 字节相对基线 | 95% 区间 |
| --- | --- | --- | --- | --- |
| Linux 页面 | −13.68% | −14.55% 至 −12.74% | −5.25% | −6.42% 至 −3.17% |
| Linux 加入大文件 | −9.18% | −10.32% 至 −8.58% | −3.04% | −5.07% 至 −1.78% |
| Wine Win64 页面 | +5.88% | +2.94% 至 +8.33% | −5.71% | −5.88% 至 −2.94% |
| Wine Win64 加入大文件 | +6.25% | +6.00% 至 +8.33% | −2.13% | −4.08% 至 −2.04% |

Linux 交错比较 23 对，页面和大文件分别每样本 25,000 和 260 轮。Wine 10.0 三个 PE 解码器在同一进程交错运行 21 组，分别每样本 5,000 和 180 轮。Wine 线程 CPU 读数约以 10 ms 为单位变化。两平台均将 Huffman 汇编入口关闭，使 C 解码路径可对照，并固定到单个允许的逻辑 CPU。基线和改写版共用同一份 O3 IR 和同一组其余 zstd 对象。

## 复现

在仓库根目录使用与 rustc 匹配的 LLVM SDK。Bazel 构建插件，从锁定的研究清单解析 `zstd-sys` 源码，并为两组链接同一份 zstd 对象库。已有源码时可用 `--zstd-lib` 指向同版本的 `zstd/lib`。

```bash
study=$(mktemp -d)
bazel run //:llvm_c_benchmark_build -- zstd --out-dir "$study"
bazel run //:llvm_benchmark_compare -- pair "$study/base" "$study/transform" --pairs 23 -- \
  100 tests/fixtures/parquet/codec_zstd.parquet tests/fixtures/hdf5/filters_plugins.h5
```

加上 `tests/fixtures/hdf5/dense_nested.h5` 参数即可纳入较大负载。Windows GNU 对照使用 `--define=rust_target=x86_64-pc-windows-gnu` 和相应的 Clang 链接环境。生成的文件为 `base.exe` 和 `transform.exe`，在 Wine 下比较时为 `pair` 加上 `--runner wine`。

正常 C 编译可用 `RSSN_SPILL_MODE=off RSSN_SCALAR_LOAD_PACK=transform` 与 `clang -O3 -fpass-plugin="$plugin"`。该插件入口同样支持报告模式。现有 `zstd-sys` 使用独立的 C 构建步骤，Rust 的 `-Zllvm-plugins` 用于 Rust 编译动作。下文的 `zstd-io` 用例为 C 依赖单独设置编译参数。

## RSSN 的 Parquet 和 HDF5 读取入口

RSSN 的 Parquet Zstd 页曾通过 `zstd::bulk::decompress` 各建一个解码上下文。HDF5 原生 Zstd filter 曾对每个 chunk 用 `zstd::stream::read::Decoder::new` 建立上下文并给切片再套一层 `BufReader`。Rust 读路径现在按一个 Parquet 文件、一个 HDF5 数据集惰性复用上下文。Parquet 保留每页的输出上限和长度校验。HDF5 从原切片流式读取，以第一帧的大小作预分配提示，保留未知大小的帧与串联帧支持。两条路径在各次独立解码之间保留 zstd 本身的会话重置语义。 [zstd 1.5.7 API](https://facebook.github.io/zstd/doc/api_manual_latest.html) [Rust zstd 0.14.4](https://docs.rs/crate/zstd/0.14.4)

对照程序直接编译 RSSN 的 Parquet、HDF5 和列实现，并锁定 `zstd-sys 2.1.3+zstd.1.5.7`。两个 C 构建都用 Clang 23.1.5、`-march=x86-64-v3` 及同一插件，保留 `zstd-sys` 自带的 C 构建 flags 和 Huffman 汇编。通过 `-Rpass-analysis=rssn-scalar-load-pack` 确认 C 编译命中三个序列解码函数共五个循环块。检查 Cargo 实际生成的 `zstd_decompress_block.o`，普通序列热循环原有三条 `movzwl` 和六条内存源 `movzbl` 用于这三张表的低 4 字节，启用 pass 后对应三条 `movl` 加寄存器内提取。

A/B 使用现有 PyArrow 和 h5py fixture，以及 RSSN writer 生成的两组分片密集文件。Parquet 文件有 48,000 行、两列、每页最多 96 行及三个行组。HDF5 文件有 250 个 `48×8` 的 F64 chunks，filter 为 shuffle 加 Zstd。计时读内存中的整个文件，包括页或块解码后的列和数据集组装。生成文件、磁盘读取及结果哈希均在计时区外。每个场景在同一逻辑 CPU 上轮换四种条件并采集 21 组线程 CPU 时间。报告组内时间比变化的中位数和 bootstrap 95% 区间。

| 完整读取入口 | 原源码开启 pass | 仅 Rust 修改 | 新源码开启 pass | Rust 修改与 pass 合计 |
| --- | --- | --- | --- | --- |
| PyArrow Parquet fixture | +0.89% [−0.26%, +4.49%] | −10.32% [−11.45%, −8.30%] | +0.08% [−1.78%, +3.18%] | −9.84% [−11.23%, −7.47%] |
| h5py HDF5 fixture | +1.27% [−1.18%, +4.02%] | −3.03% [−7.40%, −0.37%] | −0.82% [−2.05%, +3.89%] | −4.09% [−7.96%, −2.09%] |
| RSSN 多页 Parquet | +3.25% [−0.88%, +6.39%] | +0.55% [−1.36%, +4.17%] | +0.81% [−2.14%, +1.71%] | +0.05% [−2.40%, +4.55%] |
| RSSN 多块 HDF5 | +0.19% [−2.28%, +1.45%] | −11.97% [−14.64%, −9.82%] | −1.05% [−3.65%, +0.35%] | −13.20% [−15.05%, −12.43%] |

四组条件在各输入上的读取哈希分别一致，为 `2eef8a8137f63dd2`、`2ab02b4395b4fad5`、`e94ec1ce94a8047b` 和 `3317c7c4832c2ed8`。端到端可重复的改善来自 Rust 修改。底层解码的正向 A/B 来自独立解码循环及其编译配置。在这组 RSSN 入口及 Cargo C 配置下，pass 的区间均跨零，继续保持默认关闭。

### 复现完整入口 A/B

在源码修改前后各构建一次 pass 关闭版和打开版，保存为 `old-off`、`old-transform`、`new-off`、`new-transform`。基线应包含同一份 harness 和同一组生产模块，仅撤销 Parquet 与 HDF5 的上下文复用修改。用 `--out-dir` 为四次构建指定独立目录，因为 `cc` 无法感知插件的模式环境变量。

```bash
bazel run //:llvm_benchmark_build -- zstd-io --out-dir "$target_dir" \
  --env "CC=$clang" --rustflag=-Ctarget-cpu=x86-64-v3 \
  --env "CFLAGS=-march=x86-64-v3 -fpass-plugin=$plugin -Rpass-analysis=rssn-scalar-load-pack" \
  --env RSSN_SPILL_MODE=off --env RSSN_SPILL_FILTER=ZSTD_decompressSequences \
  --env "RSSN_SCALAR_LOAD_PACK=$mode"
cp "$target_dir/release/rssn-io-codec-study" "$study/$variant"
```

```bash
bazel run //:llvm_benchmark_compare -- run "$study/new-off" -- generate "$study/generated.parquet"
bazel run //:llvm_benchmark_compare -- run "$study/new-off" -- generate-hdf5 "$study/generated.h5"
bazel run //:llvm_benchmark_compare -- --cpu "$cpu" zstd-io "$study"
```

## 通用候选与分析开销

去掉函数名过滤后，编译 zstd 依赖的 C IR 还能在 `FSE_decompress_wksp_bmi2` 的两个循环块和 `HUF_decompress4X2_usingDTable_internal_fast_c_loop` 的一个循环块识别同一模式。RSSN 入口实验使用原来的序列解码函数过滤。单独针对新增模式的复核见下节。

相邻字段、读取窗口、支配关系和内存合法性检查独立于活跃度。pass 先用这些检查筛选循环块。只有同一个块找到至少两个合法 pack，才求解该函数的 SSA 活跃度并计算候选窗口的 GP 压力。一个函数最多求解一次，仅在合格块中计算逐指令压力。这样能省掉不可能改写的全函数求解。

用 LLVM 23.1.5 的 `opt -time-passes` 在 zstd 的 common、compress、decompress C IR，以及原生 Linux 和交叉 Windows 的 RSSN 求积研究模块上对照。每对交错执行新旧插件，每次调用在各函数重复 25 次报告 pass，采样 13 对。下表为各模块单次报告 pass 的 CPU 时间中位数。IR 解析时间独立计量，实验未运行后端。两版插件在三个模块以 transform 模式执行时输出相同的 IR 和候选诊断，通过 `verify`。

| 输入模块 | 先求活跃度 | 先筛合法候选 | 配对变化 |
| --- | --- | --- | --- |
| zstd C | 83.164 ms | 11.088 ms | −86.81% |
| RSSN 求积 Linux | 14.592 ms | 2.012 ms | −86.35% |
| RSSN 求积 Win64 | 5.112 ms | 1.004 ms | −80.19% |

## 四字节解码表的依赖链

在 RSSN 实际生成的 Linux Cargo 产物上给解码器函数入口设断点，只执行一次文件读取。四份输入进入 `FSE_decompress_wksp_bmi2` 的次数依次为 1、2、3、0，普通序列解码次数为 3、4、507、250。输入顺序为 PyArrow Parquet、h5py HDF5、RSSN 多页 Parquet、RSSN 多块 HDF5。四份输入均未进入 Huffman X2 的 C 或汇编快循环。把 `dense_nested.h5` 以 level 3 压成一个帧后，解码一次进入 FSE 13 次、进入 Huffman X2 快循环 8 次。运行时关闭 Huffman 汇编后，这 8 次调用同样进入 C 循环。 [zstd 1.5.7 Huffman 实现](https://github.com/facebook/zstd/blob/v1.5.7/lib/decompress/huf_decompress.c)

从原始 RSSN 读路径的 FSE 调用截取 5 种压缩后的 Huffman 权重头部，直接重放底层解码。用合并原型改写 FSE 后，23 对线程 CPU 样本的配对中位数增加 **7.05%**，bootstrap 95% 区间为 **+6.03% 至 +8.42%**。42 个大型文件头部的重放增加 **1.84%**，区间为 **−0.48% 至 +4.60%**。两组的输出哈希分别为 `c0373ecb8fff6714` 和 `00af8e9cf4839d84`。第二组来自 level 3 压缩的 `dense_nested.h5` 和 `bindings/python/rssn/_ops.py`。

两个解码器都把 16 位字段和两个相邻的 8 位字段存成 4 字节表项。读取次数可从每项三次减到一次，随后增加字段投影。若一个字段经过移位、加法等操作形成后续表项的索引，投影可能进入串行依赖链，吞吐就不能只用访存次数估算。pass 沿基本块内的 SSA 操作跟踪哪些字段影响后续候选的指针；内存操作、调用和 PHI 结束跟踪。出现这类依赖时，要求候选字段读取点的平均 GP 活跃度超过预算的两倍。这个保守门槛来自当前 A/B。FSE 的两个块平均为 20、22，预算为 12，因此保留报告并沿用原始 IR。Huffman X2 平均为 Linux 33、Win64 42，仍可合并各缓存表项的读取。

Win64 只在基址由一个动态 GEP 索引得到、元素分配大小恰好为 4 字节、且字段覆盖该元素偏移 0 到 3 时选择 4 字节组。原有序列解码表按 8 字节索引，保留原有组宽。精确字节覆盖、对齐、别名和可达性检查保持适用。工作区 GCC 16.2.0 的 debug 对象保留独立的 Huffman X2 C 快循环，另有 BMI2 与通用入口。Win64 Clang 23.1.5 重建用 `ZSTD_DISABLE_ASM`，符合锁定的 `zstd-sys` Windows 构建步骤，此时 C 快循环内联进 `HUF_decompress4X2_usingDTable_internal`。该函数的静态 spill/reload 标注从 131/280 降到 119/267。Linux 独立的 X2 C 快循环从 24/26 降到 15/16。交付 DLL 有 CodeView PDB 记录但缺少配套 PDB，这两个函数的机器码位置仍只在重建产物中确认。

Wine Win64 A/B 将四个重建解码器装入同一进程交错调用 21 组，基线和变体共享压缩器及其他 zstd 对象。保留在最终策略中的仅 Huffman X2 改写，在八个原始页和 chunk 上变化 **+1.85%**，95% 区间 **0.00% 至 +3.77%**。加入 level 3 编码的 `dense_nested.h5` 后变化 **−1.05%**，区间 **−4.30% 至 0.00%**。共享基线的输出哈希分别为 `dfdebb2e332b6efe` 和 `683487f5beb4982e`。Linux 在原生汇编启用和运行时禁用两种配置下，也未测出清晰的完整帧收益。当前证据支持保留分析和 opt-in 改写，尚无依据在默认 release 开启。

用于计数和回放的源码为 [断点脚本](analysis/zstd_entropy.gdb) 和 [FSE 回放程序](bench/zstd-fse.c)。用未剥离符号的 Linux `zstd-io` 入口及 `rounds=0` 捕获头部，每个场景恰好读取一次。将以 `FSE ` 开头的行去重后写入 `headers.txt`。基准以从同一份 IR 生成的两份 `fse_decompress.o` 和相同的 zstd 对象库链接回放程序。比较器逐组检查输出哈希和头部数。

CI 回放使用 [固定 FSE 头部](bench/fixtures/zstd-fse.hex)。它来自跟踪的
`tests/fixtures/parquet/codec_zstd.parquet`，从文件偏移 178 起取 70 字节。
断点捕获的输出容量为 255，最大表对数为 6。

```bash
bazel run //:llvm_c_benchmark_build -- zstd-fse --out-dir "$study/fse"
bazel run //:llvm_benchmark_compare -- pair "$study/fse/base" "$study/fse/transform" \
  --pairs 23 -- tools/llvm-optimizer/bench/fixtures/zstd-fse.hex 25000
```

在最终策略下 FSE 输出 IR 与直接通过 verifier 的基线一致，回放的回退由测试中的依赖链用例覆盖。

## doubleFast 短哈希与长哈希共用读取

交付 DLL 的 unwind 区间 `0x2f6ed76..0x2f6f9cb` 有与 zstd 1.5.7 哈希一致的四个 64 位乘数、32 位乘数和两套索引存储。在默认分支，RVA `0x2f6ee3d` 的短哈希乘法读 `(%rdi)` 的低 32 位，到汇合点 `0x2f6ee53` 又从同址读 64 位做长哈希。另一个同族区间 `0x2f68560..0x2f69064` 的 `0x2f68709` 已经用一次宽读取分别计算两种哈希。无 PDB 的 DLL 可定位这些指令与相邻控制流，符号归因以有名字的重建为准。工作区 GCC debug 对象的 `ZSTD_fillDoubleHashTable` 在节内 `0x115` 与 `0x126` 分别调用 `ZSTD_hash4Ptr` 和 `ZSTD_hash8Ptr`，两个调用使用同一指针，重定位可核对。Debug 的 helper 调用边界与 release 内联配置不同。 [zstd 1.5.7 双哈希实现](https://github.com/facebook/zstd/blob/v1.5.7/lib/compress/zstd_double_fast.c)

Clang 23.1.5 的 O3 IR 在 `ZSTD_fillDoubleHashTable` 的两个循环里保留了相同模式。五个 switch 臂都读取同一个 64 位值，默认臂还先读一次同址的 32 位值。对同一 Win64 IR 追加 `simplifycfg<hoist-common-insts>,gvn` 或单独使用 `newgvn`，两个候选仍在。

`rssn-switch-load-reuse` 先检查每个不同的后继块只来自该 switch。每个后继都必须在短前缀内从同一个支配 switch 的 SSA 指针加载同宽整数。前缀不能有调用、写入、特殊控制流、原子或 volatile 操作。至少一个臂还必须有同址的窄整数读取，并且目标认为这个宽度的低位投影可以免费实现。宽读取移到分派之前，使用原宽读取中最小的对齐，所有执行路径仍读原来的最大范围。原有 CFG 与哈希运算保持。这样仅跨分派增加一个宽值的活跃区间，收益通过实际机器码验证。

这个合并使用 LLVM 23 的 `b64` 等字节类型保留局部 poison。窄读取先从原始字节中投影，再转成整数。直接截断 `i64` 会把无关高位的 poison 传播到窄结果。早期 `<8 x i8>` 原型虽保留语义，却产生了 `vmovq` 和 `vmovd`。最终方案用 `bitcast b64 to <2 x b32>` 再取 lane 0，重建的 Linux 和 Win64 都用同一 GPR 的低 32 位计算窄哈希。两个候选各从五条静态宽读取和一条窄读取合并为一次宽读取。Win64 分配后的 MIR 由内存操作数乘法 `IMUL32rmi` 变为寄存器版本 `IMUL32rri`。整个函数的静态指令数在 Linux 从 243 降到 239，在 Win64 从 257 降到 253。Win64 整函数仍有一条 spill 和一条 reload，Linux 重建没有相应标注。报表模式与关闭模式输出相同 IR，IR 和机器指令 verifier 均通过。 [LLVM poison 语义](https://llvm.org/docs/LangRef.html#poison-values) · [LLVM 字节类型](https://blog.llvm.org/posts/2025-08-29-gsoc-byte-type/)

开关为 `RSSN_SWITCH_LOAD_REUSE=off|report|transform`，默认关闭。跟随 zstd-sys C 编译动作使用 Clang 23 的 `-fpass-plugin`，可用 `-Rpass-analysis=rssn-switch-load-reuse` 验证是否命中。显式对照从一份 O3 IR 切换 pass，另链接相同的库对象。选择 `double-fast` 后，构建入口只改写 `zstd_double_fast.c`。

```bash
bazel run //:llvm_c_benchmark_build -- zstd --study double-fast --out-dir "$study/dfast"
tail -c 65536 tests/fixtures/hdf5/dense_nested.h5 > "$study/sample.dict"
bazel run //:llvm_benchmark_compare -- pair "$study/dfast/base" "$study/dfast/transform" --pairs 17 -- \
  --compress-min4 2500 tests/fixtures/hdf5/filters_plugins.h5 "$study/sample.dict"
```

`--compress` 保留 level 3 的压缩参数，并可附带原始字典。无字典的大文件对照没有进入表初始化，默认参数的原始字典输入进入初始化时 `minMatch=5`。定向对照以真实 HDF5 fixture 为输入，用另一个 HDF5 fixture 的最后 64 KiB 作为 raw prefix，并显式设置 doubleFast 和 `minMatch=4`。GDB 确认这组进入目标分支。各样本在计时窗口外逐字节校验解压结果，A/B 进一步比较编码后帧的哈希。

固定单核交替执行基线和改写版。默认参数在 Linux 跑 13 对，在 Wine Win64 跑 11 对。定向参数分别跑 17 和 15 对。两平台各用例的编码哈希一致，所有逐字节回读校验通过。

| 输入 | Linux 线程 CPU 变化（95% 区间） | Wine Win64 线程 CPU 变化（95% 区间） | 编码哈希 |
| --- | ---: | ---: | --- |
| 大 HDF5，默认参数，无字典 | +0.33% [−7.89%, +5.22%] | +2.27% [−2.28%, +4.70%] | `ae08e240c81f3a34` |
| 大 HDF5，默认参数，原始字典 | −9.03% [−14.13%, +0.99%] | +0.42% [−2.82%, +4.10%] | `8f2d3002cbdcc3e5` |
| 小 HDF5，默认参数，原始字典 | +0.68% [−1.16%, +2.10%] | 0.00% [−9.09%, +5.00%] | `1ac413971e677e5f` |
| 小 HDF5，doubleFast，`minMatch=4` | +1.07% [−1.04%, +13.45%] | −3.23% [−3.33%, 0.00%] | `fdbf1150e3ff9631` |

Wine 的定向结果接近 10 ms 线程计时粒度，区间触及零。其余组合的波动也覆盖零。这轮确认了重复读取和生成代码，无默认发布性能结论。

## 相邻字段合并的字节载体消融

交付 DLL 的三张 FSE 表具有与锁定依赖一致的字节指纹。工作区 GCC debug 对象还能通过符号区分普通、分段字面量和预取解码器。DLL 对应循环的位置仍未由符号或重定位确认，本节的指令和计数取自同一份 Clang 23.1.5 AVX2 IR 重建。

在 `rssn-scalar-load-pack` 上只改载体和字段投影。对照沿用 `<4 x i8>` 或 `<8 x i8>`。实验用 `b32` 或 `b64` 作一次读取，按字段宽度转换为等宽字节向量并提取对应 lane。非自然对齐的字段从 `b8` 向量切片。两者均在投影后转成整数，保留逐字节的 poison 隔离。候选、合并宽度、取数位置、字段使用位置和盈利门槛相同。两版都通过 IR 和机器指令验证。 [LLVM 字节类型](https://llvm.org/docs/LangRef.html#byte-type)

在 Win64 `ZSTD_decompressSequencesLong`，原始字节载体相比向量载体将静态 spill/reload 从 110/363 降到 108/339。普通序列解码仍为 46/92。Linux 普通序列解码从 43/85 变成 44/84，计数近似，实际耗时却明显增加。两种表示经过目标合法化与寄存器分配得到不同的取数、移位、寄存器复制和循环布局。这里只定位到差异来自 IR 表示和后端代码生成，静态槽位数无法解释全部运行差异。

下表为单核线程 CPU 时间的配对变化，区间为 bootstrap 95%。默认页面和大文件输入沿用本研究前面的 RSSN 帧。Linux 每用例轮换六种顺序取 18 组，Wine 默认输入取 21 组。冷字典用例两平台各取 18 组。所有组的解码字节哈希均一致。

| 平台与输入 | 向量载体相对关闭 | 原始字节载体相对向量载体 |
| --- | ---: | ---: |
| Linux 原始页面 | −13.90% [−15.24%, −11.63%] | +15.10% [+11.74%, +17.34%] |
| Linux 页面与大文件 | −9.91% [−11.63%, −7.39%] | +13.51% [+9.51%, +15.43%] |
| Linux 冷字典预取解码 | +0.99% [−0.34%, +2.88%] | +1.40% [+0.49%, +4.08%] |
| Wine Win64 原始页面 | −3.33% [−6.45%, 0.00%] | 0.00% [0.00%, +3.57%] |
| Wine Win64 页面与大文件 | −3.23% [−3.70%, 0.00%] | 0.00% [0.00%, +4.00%] |
| Wine Win64 冷字典预取解码 | 0.00% [0.00%, 0.00%] | 0.00% [−5.56%, +11.11%] |

冷字典对照压缩仓库中的 `dense_order.h5`，使用 `dense_nested.h5` 的最后 64 KiB 作为原始字典。解码交替使用内容相同且底层地址不同的两份 `DDict`。GDB 确认它通过正常入口进入 `ZSTD_decompressSequencesLong`。这组刻意触发预取分支。Wine 子样本的中位耗时约 80 至 90 ms，线程计时的 10 ms 粒度限制了可解释的效应。新增的基准模式在计时之外构建字典、压缩并逐字节验证回读：

```bash
tail -c 65536 tests/fixtures/hdf5/dense_nested.h5 > "$study/sample.dict"
bazel run //:llvm_benchmark_compare -- pair "$study/base" "$study/transform" --pairs 18 -- \
  --decode-cold-dict 700 tests/fixtures/hdf5/dense_order.h5 "$study/sample.dict"
```

Linux 原始页面还试了把字节投影推迟到各字段的第一次使用。在同批 12 组数据中，这仍比向量载体慢 **20.69%**，区间 **+19.02% 至 +23.98%**。另在单独构建中撤除 FSE 依赖链门槛，强制处理两组已捕获 Huffman 权重头部。原始字节载体相对未合并在 RSSN 页头部上慢 **9.27%** [**+7.15%, +13.02%**]，在 42 个较大文件的头部上快 **5.54%** [**−9.64%, −2.26%**]。后者相对向量载体为 **−1.69%** [**−10.06%, +2.38%**]。生产选择沿用向量载体和依赖链门槛，原始字节载体继续用于上一节已经单独验证的分支共用读取。

共享压力模型现在用对应的整数机器类型对 LLVM 原始字节值计权。比如跨调用的 `b64` 会占一个 GP 估计单位。两个平台的 zstd 序列解码用最终插件改写出的 IR 与现有向量方案逐字节相同。
