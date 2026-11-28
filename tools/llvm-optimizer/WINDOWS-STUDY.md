# Windows release DLL 与 RSSN 双指数内核

## 交付 DLL 的反汇编

样本 `target/rssn.dll` 为 x86-64 PE32+，SHA-256 为 `c2fdbb8a7b2f606c20690b5393a4fd8ea92b7d85adff9593f9fd8209cdf9d645`，映射基址为 `0x180000000`。CodeView 记录引用 `rssn.pdb`，GUID `63A3538F-CDCF-1B3F-4C4C-44205044422E`，age 1。工作区缺少配套 PDB 和同一构建的 Windows debug DLL，因此以下地址按解码后的指令、导入和控制流识别。当前源码的 debug 与 release 重建用于跟踪具体 IR 值和 ABI。它们的链接方式及编译选项可以独立于交付 DLL。 [PE 格式](https://learn.microsoft.com/en-us/windows/win32/debug/pe-format)

| RVA | 观察 |
| --- | --- |
| `0x2b6a30` | 外形与 Appell F₁ 的 tanh-sinh 实例吻合。首个比较交换逆序的端点并递归，使用 32 字节结果缓冲区。内循环的一组四次 `pow` 分别接收 `t`、`1-t` 及两组 `1-k*t` 形态的基数，累乘顺序与 RSSN 源码一致。首个 `pow` 在 `0x2b6e63`，返回值经 XMM10 穿过后续调用。入口保存八个 GPR 和全部 XMM6–XMM15。 |
| `0x2bd7f0` | `0x2bd9ae` 和 `0x2bd9c1` 的 `exp` 构成一组。两次调用之间 XMM8 保持前一次返回值。循环同时对其他值执行帧访问，例如 `0x2bd990` 和 `0x2bd996` 的内存形式乘法。 |
| `0x176d4b0` | RSSN CPU lane 后端的元素级分发。循环 `0x176e3d0` 内展开八次 `pow`，每次从两个输入寄存器的数组各读一个标量并立即存储结果。`0x1770f28` 起的八次 `exp` 同样直接读写输入和输出数组。整数幂臂另有八条重复的指数循环，详见下文。 |
| `0xe41990` | Tape 列执行器。`0xe41a3f` 计算不超过 256 行的 stride，按 24 字节解析指令。首个 `atan2` 在 `0xe4242f`，整数幂组从 `0xe431d9` 起按行重复移位。 |

PE 扫描共得到 36,079 个 unwind 区间和 3,903 个候选 libm 调用点。候选统计按静态指令，不能当作执行频次。上表地址从完整函数边界反汇编中复核。Win64 的 XMM6–XMM15 低 128 位由被调用者保存，YMM 和 ZMM 的高位为易失状态。入口的保存、调用者栈传参和寄存器分配产生的槽位需要分别跟踪。 [Microsoft x64 调用约定](https://learn.microsoft.com/en-us/cpp/build/x64-calling-convention?view=msvc-170)

同一源码的 Win64 debug Appell 外壳用 XMM0–XMM3 接收前四个参数，从入口栈槽读入 `x,y`，再调用单独的积分函数。fat LTO release 将积分循环内联到外壳。它在四连 `pow` 中从 `96/80/64/48(%rsp)` 重读四个指数。后两项从 `368/376(%rsp)` 读取 `x,y`。后两者是调用者入参，前四者可在寄存器分配后的 MIR 中跟踪为内部栈槽。第一项结果跨过第二次调用保存在 XMM15 或 XMM10，随后在 XMM6 或 XMM8 内累乘。这两种编译策略产生的静态栈访问数不能直接相减来估算 release 性能。

在一份同源 Win64 Appell 重建中，`rssn-machine-stack-traffic` 以 `virtregrewriter` 后的机器内存操作数为单位，报告外壳内 25 个静态调用点。该函数有 12 次 spill 读取和 9 次 spill 写入，读中 5 次宽于 8 字节。调用之后的同一块内有 4 次 spill 读取、6 次不可变固定栈槽读取和 13 次常量池读取。这个分解用于区分寄存器分配产生的流量、ABI 固定栈对象和后端物化，计数范围为重建的外壳函数。

Appell 的 Win64 IR 试验保持四个 `pow` 及累乘顺序，将后两个基数的乘和减合并为两 lane 指令。局部构造向量时，两个栈入参仍逐项读取，随后增加 lane 组装和提取，函数帧从 312 增至 328 字节。将 `x,y` 载体提升到循环外把循环内两次栈读合并成一次 16 字节折叠读取，但整函数帧增至 344 字节，带注释的 spill/reload 从 21/33 增至 23/39。原先第四个基数在第三次 `pow` 之后生成，改写时包含第四个基数的向量跨越第三次 `pow`。这个试验留在分析阶段。

## 已实现的向量 libm 分支恢复

RSSN 小正阶 Lerch 积分的 `z > 0` 对预备点的所有 lane 相同。优化后的 IR 把条件值广播到向量，并用 `select` 在 `llvm.exp` 与已计算的另一结果之间选择。LLVM 23.1.5 x86 SelectionDAG 在本次无向量数学库的配置中把这些向量 `exp` 降为逐 lane 的标量调用。调用之后才做选择会延长其他值的跨调用区间，并执行本可跳过的一组 libm。控制流恢复将单用途 intrinsic 放进选中分支，由 PHI 选择结果。相邻且使用同一谓词、方向和权重的候选可以共享分支。如果某一结果在新汇合点之前已被消费，就分开处理。 [LLVM SelectionDAG](https://llvm.org/docs/CodeGenerator.html)

每次到达该分支的预期收益包括 `p_skip × N × C_libm` 和节省的跨调用保存开销，再扣除增加的分支、汇合点搬运和寄存器压力。这里 `N` 是标量化后调用数。有显式分支权重时，要求 `p_skip × N ≥ 0.5`，并把权重传给生成的分支。没有权重时，只在已验证的 x86 双精度 `exp`、`log` 及 2/4/8 lane 范围内试用。已配置向量数学库、可能采用近似实现、strict FP、非统一掩码、共享结果和有副作用的外部 libm 均排除。进入分支前用 LLVM `freeze` 固定谓词，避免把 `select` 原本可传播的 poison 变成分支上的立即未定义行为。收益门槛是启发式，最终仍以机器码和运行 A/B 判定。 [LLVM IR 语义](https://llvm.org/docs/LangRef.html)

| 同一小正阶 Lerch 内核 | 改写 | 静态 Spill / Reload |
| --- | --- | --- |
| Linux x86-64，`target-cpu=native` | 2 组共 10 个 `exp` 移入条件分支 | 138 / 216 → 133 / 211 |
| Win64，`x86-64-v4`，fat LTO IR | 3 组共 8 个 `exp` 共享一个分支，栈帧 1104 → 1056 字节 | 61 / 100 → 57 / 102 |

计数包含机器码生成器标注的折叠栈访问。两侧 Win64 版本仍各保存十个非易失 XMM。条件分支保留全部静态 `exp` 调用，在绕过时省掉执行。Linux A/B 用四组从定义域选取的真实输入，先对测试进程里的 `exp` 计数，另做不带插桩的计时。小正阶 Lerch 每轮四次求值减少 20 次 `exp`（7,588 → 7,568，扣除预热），全部结果逐位一致。

| Linux，21 对交错运行 | A / B 中位耗时 | 配对耗时变化 | 95% bootstrap 区间 | 位级哈希 |
| --- | --- | --- | --- | --- |
| 小正阶，650 轮 | 260.5 / 281.2 ms | +3.91% | −3.98% 至 +7.99% | `5a9c07fab0b78372` |
| 正阶对照，650 轮 | 330.2 / 320.1 ms | −5.27% | −13.48% 至 +4.16% | `6d2ea9a36fae5b97` |
| 负阶对照，1,600 轮 | 346.7 / 367.3 ms | +4.07% | −7.10% 至 +14.13% | `2c325b399be58a33` |

计时使用 i5-7300U 上同一个逻辑 CPU，未独占机器，两个对照也有可见波动。区间均跨过零，当前保留显式启用。Win64 只验证 IR、MIR、汇编和 COFF 机器码。SysV 标量算术锥移动的两个试验各减少两个 SSA 跨调用活跃值，实际寄存器分配计数维持 138 / 216，因此该改写仅在 `function(rssn-libm-cone-trial)` 下运行。

## 复现

从 RSSN 根目录运行。变量 `plugin` 指向 `RssnLLVMOptimizer.so`，`llvm_bin` 指向匹配 rustc 版本的 LLVM SDK 工具目录。

```bash
python3 tools/llvm-optimizer/analysis/pe_math_calls.py target/rssn.dll
"$llvm_bin/llvm-objdump" -d --start-address=0x1802b6a30 \
  --stop-address=0x1802b74fd target/rssn.dll

study=$(mktemp -d)
bazel run --define=rust_target=x86_64-pc-windows-msvc //:llvm_benchmark_build -- de \
  --profile dev --lib --out-dir "$study/win" --env CARGO_PROFILE_DEV_CODEGEN_UNITS=1 \
  --env CARGO_PROFILE_DEV_PANIC=abort --rustflag=-Ctarget-cpu=x86-64-v4 \
  --rustflag=-Zverbose-asm -- --emit=llvm-ir,asm,obj

bazel run --define=rust_target=x86_64-pc-windows-msvc //:llvm_benchmark_build -- de \
  --out-dir "$study/win" --rustflag=-Ctarget-cpu=x86-64-v4 --rustflag=-Zverbose-asm \
  -- "--emit=llvm-ir=$study/win-full.ll,asm,obj,link" -Zno-link

"$llvm_bin/llvm-extract" -S --recursive --func=rssn_de_small "$study/win-full.ll" -o "$study/win.ll"
"$llvm_bin/opt" -S -mcpu=x86-64-v4 -load-pass-plugin="$plugin" \
  -passes='function(rssn-libm-pressure)' -pass-remarks-analysis=rssn-libm-pressure \
  -verify-each "$study/win.ll" -o "$study/win-guard.ll"
for stem in win win-guard; do
  for phase in stop-before=greedy stop-after=virtregrewriter; do
    "$llvm_bin/llc" -O3 -mcpu=x86-64-v4 "-$phase" "$study/$stem.ll" -o "$study/$stem-${phase#*=}.mir"
  done
  "$llvm_bin/llc" -O3 -mcpu=x86-64-v4 "$study/$stem.ll" -o "$study/$stem.s"
  "$llvm_bin/llc" -O3 -mcpu=x86-64-v4 -filetype=obj "$study/$stem.ll" -o "$study/$stem.obj"
done

for mode in off transform; do
  bazel run //:llvm_benchmark_build -- de --out-dir "$study/$mode" \
    --env RSSN_SPILL_MODE=off --env RSSN_SPILL_ABI=off \
    --env RSSN_SPILL_FILTER=rssn_de_small --env "RSSN_SPILL_LIBM=$mode" \
    --rustflag=-Ctarget-cpu=native --rustflag=-Zverbose-asm --rustflag=-Zverify-llvm-ir \
    --rustflag="-Zllvm-plugins=$plugin" -- --emit=llvm-ir,asm,link
done

# 将 3 替换为当前进程可用的 CPU 编号。
for sample in "small 650" "positive 650" "negative 1600"; do
  read -r case rounds <<< "$sample"
  bazel run //:llvm_benchmark_compare -- --cpu 3 pair \
    "$study/off/release/rssn-de-study" "$study/transform/release/rssn-de-study" \
    --pairs 21 -- "$case" "$rounds"
done
```

可在计时之外编译调用计数探针。每组以 `small 1` 的计数减去 `small 0` 的计数，扣除两个进程都会执行的预热。

```bash
bazel run //:llvm_benchmark_build -- count-exp --out-dir "$study"
for mode in off transform; do
  for rounds in 0 1; do
    bazel run //:llvm_benchmark_compare -- run "$study/$mode/release/rssn-de-study" \
      --env "LD_PRELOAD=$study/count_exp.so" -- small "$rounds"
  done
done
```

Appell 案例从 RSSN 提取 F₁ 的积分分支和实际 gamma、tanh-sinh 内核。运行 `rssn-de-study appell 40` 可测试四组会走积分分支的参数。Windows IR 可用同一套命令抽取 `rssn_de_appell`。

## Appell F₁ 正负节点的标量复用

交付 DLL 在 `0x2b6da9` 调用第一个 `exp`，立即执行加法和除法。四连 `pow` 之后，`0x2b6ed0` 的 XOR 翻转本轮求积节点的符号，`0x2b6f29` 再次调用 `exp` 并执行相同的加法和除法。当前 RSSN 源码的 Win64 和 SysV fat LTO IR 在初始网格和细化网格中都复现该结构。两边都保留 `sinh`、`cosh` 调用，它们在 IR 中可写 `errnomem`。 [LLVM 普通数学 intrinsic 与严格浮点语义](https://llvm.org/docs/LangRef.html)

改写识别同一循环中以同一个节点正负值为输入的两条链。`sinh` 形态只用于筛选高概率复用机会。编译器在原 `llvm.exp` 调用位置比较两个输入的 64 位表示。首次结果支配第二次计算，第二次所在块后支配首次所在块。输入一致时复用首条链的结果，一旦有任何位不同便执行原操作。两个输入在调用和比较中共用各自的 `freeze` 后值。后续算术只匹配同一常量位模式、同一操作数位置与无 fast-math 的单用途链。遇到任意外部函数、规定浮点环境或者无法支配使用点时停止延长。运行回归用可带副作用且可返回非奇对称值的 `sinh` 覆盖两条分支，检查调用次数与输出位模式。

收益模型以每个配对的期望开销为单位。若命中率为 `p`，则可跳过的期望工作量为 `p × (C_exp + C_add + C_div)`，由比较和分支开销以及输入位模式、已算结果的额外活跃区间抵扣。后者在 SysV 与 Win64 的保存预算不同，可能使静态 spill 增加。仅复用 `exp` 的试验节约同样多的调用，计时变化 −2.56%，95% 区间 −10.43% 至 +7.59%。实际改写只处理共同算术链包含双精度除法的候选。分支权重的 15:1 用于布局，标记为启发式预期。

| 四组真实 RSSN Appell 输入 | 原 IR | 带按位守卫的复用 |
| --- | --- | --- |
| Linux `exp` 次数，每轮四次求值 | 14,144 | 7,072 |
| Linux 静态 Spill / Reload | 59 / 104 | 65 / 110 |
| Win64 静态 Spill / Reload | 21 / 33 | 22 / 34 |
| Win64 函数帧 | 312 字节 | 328 字节 |

Win64 两组都保存十个非易失 XMM。静态计数含慢路径。Linux 两项实验各交错运行 31 对、每次 38 轮，位级哈希均为 `e8872e9d92154ed5`。

| Linux A/B | 原版 / 改写中位耗时 | 配对耗时变化 | 95% bootstrap 区间 |
| --- | --- | --- | --- |
| 同一份 fat LTO IR 经 `opt` 与 `llc` 生成两组对象，链接相同基准外壳 | 100.2 / 95.6 ms | −6.47% | −8.62% 至 −1.96% |
| cargo rustc 正常 fat LTO 入口，插件环境变量切换 | 96.6 / 91.5 ms | −1.94% | −5.62% 至 +0.39% |

Rust 入口的三个 Lerch 对照输出也逐位一致。第二组计时区间跨过零，改写保留显式启用。Win64 同时通过了 `opt` 对照和 rustc fat LTO 的插件加载、IR 校验与 COFF 生成，两个慢路径出现在最终汇编中。Win64 证据为静态机器码和 MIR。

将上述 Linux 构建命令中的过滤器改为 `rssn_de_appell`，并用 `RSSN_SPILL_LIBM_REUSE="$mode" RSSN_SPILL_LIBM=off` 切换两组。计数探针与计时脚本均可直接用 `appell` 案例。对 Win64 可抽取同一函数，改用 `-passes='function(rssn-libm-reuse)'` 生成对照 IR。

## Win64 跨 pow 调用保留指数对

原始 Appell release MIR 的四个指数是四个独立的虚拟 SIMD 活跃区间。普通 Win64 CRT 调用保留 XMM6–XMM15 的低 128 位。实际热块的 `pow` 调用点有 16 个 SIMD 值跨越调用，四个指数分配到独立的槽位，每个节点各重读一遍。两个 `1-k*t` 系数仍由调用者栈传入。debug 产物保留单独积分函数，release 中的循环和栈槽归因来自重建的 fat LTO MIR。

pass 枚举同一基本块内的相邻标量 `pow`，只接受跨多个静态窗口重用的循环不变指数。一对指数的 IR 使用必须全部位于打包点之前或被本次改写消耗，否则标量区间和向量区间会同时占用寄存器。打包位于使用循环的前置块或它们的公共支配点。

每个窗口分别查看两次调用的活跃区间。至少一次调用需要同时保留两个指数，而且在该点活跃的 XMM 型值超过十个，才按块频率给这个窗口计一次可消除的标量读取。原生 YMM/ZMM 值的保存需要额外处理，不能占用同一个非易失 XMM 容量。成本用 TTI 的两次 lane 插入估计前置块组装，用高 lane shuffle 的吞吐代价计每个消费窗口的开销。一次预期消除的依赖读取得分为标量 load 吞吐代价的两倍。这个系数近似涵盖发射和数据依赖延迟，实际 spill 与节省时间仍由后端和 A/B 校验。最多十二个指数参与最大权匹配，每提交一对，就在它能消除重叠的所有已建模调用点扣减一个超限槽，再为余下配对计算边际得分。 [Microsoft x64 寄存器约定](https://learn.microsoft.com/en-us/cpp/build/x64-software-conventions?view=msvc-170)

单纯的 LLVM `insertelement` 与 `extractelement` 会被 SelectionDAG 拆成标量，实际四个指数的栈槽仍然分离。改写用两个 lane 的 `freeze` 保持每个已定义指数的位模式。低 lane 使用输入输出绑定在同一 XMM 的空汇编，高 lane 使用一条 `vunpckhpd`。每个 `pow` 的位置、签名、实参位模式与累乘顺序保持一致。生成的机器码在四连 `pow` 中只加载两个 16 字节指数槽，第二次取指数从跨调用保留的寄存器读取。高 lane 增加一条寄存器 shuffle，Win64 ABI 的寄存器搬运仍有吞吐成本。 [LLVM inline asm 约束](https://llvm.org/docs/LangRef.html#inline-assembler-expressions)

| Win64 同一份 Appell IR | 整函数帧 | 静态 Spill / Reload | 静态 `pow` 调用 |
| --- | --- | --- | --- |
| `x86-64-v3` | 296 → 264 字节 | 20 / 32 → 17 / 23 | 18 |
| `x86-64-v4` | 312 → 280 字节 | 21 / 33 → 18 / 24 | 18 |

运行验证使用本机可以执行的 AVX2 产物。以与交付 DLL 相同的 `api-ms-win-crt-math-l1-1-0.dll` 导入生成四个实验 DLL，在 Wine 10.0 的 UCRT 实现下，用同一个进程交错运行原版、指数打包、标量复用和两项同时开启。AVX-512 产物用于静态比较。

| Wine Win64，24 组交错运行 | 配对墙钟耗时变化 | 95% bootstrap 区间 | 配对线程 CPU 时间变化 | 95% bootstrap 区间 |
| --- | --- | --- | --- | --- |
| 仅指数打包 | −1.80% | −9.59% 至 +4.61% | 0.00% | −4.44% 至 0.00% |
| 仅标量复用 | −12.38% | −19.54% 至 −7.99% | −7.14% | −9.30% 至 −4.65% |
| 两项开启 | −9.75% | −16.44% 至 −1.98% | −4.88% | −7.14% 至 −3.63% |

每个样本计算 120 轮，每轮四组真实 RSSN Appell 输入。运行器逐项核对预热输出，并验证所有计时样本的位级哈希均为 `8871e0382fe46c61`。固定在同一个逻辑 CPU，机器未独占。墙钟来自 `QueryPerformanceCounter`，线程 CPU 时间来自 `GetThreadTimes`，本机观察到约 10 ms 的分辨率。较短的 55 轮交错实验中，单独打包的墙钟变化为 +2.57%，区间 −19.80% 至 +23.94%。指数打包减少了静态重读，整段计时没有确认单独收益，因此默认关闭。正常 rustc fat LTO 插件入口也生成了 264 字节栈帧和 17 / 23 的静态 Spill / Reload。

同一函数的 Gamma 系数乘除有独立的 [退出表达式实验](#appell-退出表达式与跨调用活跃区间)。

### 按调用复核打包收益

交付 DLL 的第一组循环内 `pow` 在 RVA `0x2b6e63/0x2b6e83/0x2b6ea4/0x2b6ec1`。前两个指数计算通过 `0x30(%rsp)` 和 `0x28(%rsp)` 重读捕获地址，再沿指针取系数。后两个指数从 R15 和 R13 指向的对象取值。同源码的 Win64 调试重建把四次 `powf` 保留在单节点 integrand 中，四个调用前也分别经指针取得系数。当前源码 fat LTO 的调用形态不同。分配后 MIR 在相应的四连 `pow` 中从 `%stack.0` 至 `%stack.3` 各载入 128 位至 XMM1，每处只消费低位标量。`%fixed-stack.1/.0` 分别提供两个基数系数，仍属 ABI 栈参数。

新版模型在当前 Appell AVX2 重建中选中两组配对，每组四个窗口、八个指数同时活跃的调用。新旧模型输出的 PE `.text` 逐字节相同。与关闭打包的 AVX2 基线比较，局部帧为 296 → 264 字节，静态 spill/reload 为 20/32 → 18/24。在受检四连 `pow` 块，打包后只重读 `%stack.0/.1`，能跨调用保留载体的 XMM8、XMM12 随后用 `vunpckhpd` 提取高位指数。两个 `%fixed-stack` 读取仍存在。AVX-512 重建同样选中两组，并通过机器验证。

另测试了允许内存操作数的 `vpermilpd` 提取。LLVM 为八个高 lane 使用点生成独立的栈临时槽和内存形式 shuffle，局部帧增至 392 字节。普通 spill/reload 注释计数不含这些内联汇编内的读取。该提取方式保留为实验对照。

四种 AVX2 DLL 用同一 Wine 进程交错运行 24 组，每组 240 轮，真实 Appell 输出哈希均为 `7b5ac432739520d1`。新版与无打包基线的配对线程 CPU 时间变化为 −1.61% [−4.69%, 0.00%]，与机器码相同的旧打包版为 0.00% [−4.02%, 0.00%]。内存形式 `vpermilpd` 相对新版为 +1.60% [−0.75%, +2.17%]。这些耗时保留了 Wine 时钟量化和共享机器上的变动，打包仍需显式启用。额外回归覆盖只有原生 YMM 溢出的情况，以及两组指数竞争一个超限 XMM 槽的情况。

### 复现 PE 运行验证

使用匹配的 LLVM SDK、Clang、LLD 和 Rust 的 MSVC 目标。Bazel 生成 fat LTO IR、提取内核，并链接四个 DLL 和运行器。已有同版本 Win64 IR 时，可通过 `--ir` 复用。Windows 可直接运行，在 Linux 上可使用已配置的 Wine。

```bash
bazel run --define=rust_target=x86_64-pc-windows-msvc //:llvm_c_benchmark_build -- \
  win-appell --out-dir "$study/pe"
bazel run //:llvm_benchmark_compare -- run "$study/pe/appell-suite.exe" --runner wine \
  --env RSSN_WIN_ROUNDS=120 --env RSSN_WIN_PAIRS=24
```

可用重复的 `--variant NAME=PASSES` 为同一内核指定两组或四组管线，其中固定使用 `base=verify`。对于 Lerch，使用 `--kernel lerch-small`、`lerch-positive` 或 `lerch-negative`。DLL 名称与运行器的实验开关对应。

输出每行依次为版本、墙钟纳秒、线程 CPU 纳秒和位级哈希。每四行形成一组，版本顺序随组次轮换。

### 用 Windows opt 加载插件

将插件交叉编译为 x64 PE DLL，使用 MSYS2 CLANG64 的 LLVM 23.1.5 头文件、C++ runtime 和 `libLLVM-23.dll.a`。插件导出 `llvmGetPassPluginInfo`，同一 `libLLVM-23.dll` 为插件与 `opt.exe` 提供 C++ 符号。以下是同一 SDK 在 Windows CLANG64 shell 中的构建配置。确保 `rustc -vV` 报告 LLVM 23.1.5，且 `python`、`cmake` 和 `ninja` 来自 CLANG64。 [MSYS2 LLVM SDK](https://packages.msys2.org/packages/mingw-w64-clang-x86_64-llvm?repo=clang64) [工具及导入库](https://packages.msys2.org/packages/mingw-w64-clang-x86_64-llvm-tools?repo=clang64)

```bash
sdk=$(cygpath -m /clang64)
bazel test --repo_env=CMAKE_GENERATOR=Ninja \
  --repo_env=CXX="$sdk/bin/clang++.exe" \
  --repo_env=RSSN_LLVM_CONFIG="$sdk/bin/llvm-config.exe" \
  --repo_env=RSSN_LLVM_IMPORT_LIBRARY="$sdk/lib/libLLVM-23.dll.a" \
  //:llvm_optimizer_tools_test
```

实际验证使用 Linux 上的 Clang 23.1.5 交叉构建，并在 Wine 10.0 中执行 Windows 版 `opt.exe` 和 `llc.exe`。对 RSSN 的 MSVC AVX2 Appell 使用相同 IR 做 A/B。窗口内跳过三个单位底数 `pow` 的保护覆盖四个循环块。两组均通过 IR verifier 与 `-verify-machineinstrs`，反汇编与对应 Linux `llc` 生成的指令逐字节一致。独立链接后对四组真实 RSSN 输入交错测量 20 对样本，每样本 150 轮，哈希均为 `08344d98c3ccfef2`。配对线程 CPU 时间变化为 **−10.26% [−10.53%, −9.05%]**，墙钟为 **−7.95% [−22.01%, +5.40%]**。

同一插件 DLL 也处理了 zstd 1.5.7 的 Win64 解码和压缩 IR。相邻标量读取分析在三个序列解码函数中选中五个循环块。多分支读取复用在 doubleFast 的长短哈希路径选中两个共享读组。两份 IR 均通过机器指令验证。这里验证的是 Windows 进程中的 LLVM 插件加载。官方 MSVC rustc 1.100.0 的驱动仍缺少插件需要的 LLVM C++ 导出。 [依赖构建与完整入口 A/B](ZSTD-STUDY.md)

## Appell 退出表达式与跨调用活跃区间

交付 DLL 的 Appell 外层代码位于 `0x1811c6a60`。它用 `c`、`a`、`c-a` 三个输入调用同一个浮点例程 `0x1812fa360`，随后在 `0x1811c6c72` 和 `0x1811c6c76` 做乘除。系数由 XMM6 保持，跨过 `0x1811c6c93` 对独立求积函数 `0x1802b6a30` 的调用。配合 Appell 源码的三个 `gamma_numerical` 调用，此序列能归因到归一化系数。缺少配套 PDB，地址是反汇编归因。交付 DLL 已在进入求积器之前完成这段算术。

当前源码的 fat LTO 构建把求积器内联。LLVM 23.1.5 的 pass trace 显示，内联后 Gamma 乘除仍在第一个求积循环之前，`IndVarSimplifyPass` 随后把两条指令移到唯一退出块。其 `sinkUnusedInvariants` 向下移动循环内未使用的前置块指令。这个实例由三个浮点操作数取代了一个结果跨越循环，增加了调用点活跃度。重建过程的函数内联与调度差异解释了为什么同一系数在交付 DLL 和当前 fat LTO 机器码里占用不同数量的活跃区间。 [LLVM IndVarSimplify 源码](https://www.llvm.org/doxygen/IndVarSimplify_8cpp_source.html)

新 pass 从顶层循环的退出块选取小型、可安全提前执行的纯浮点表达式，只移动原有指令。所有外部输入必须支配放置点，原有浮点操作顺序保持不变。候选先以试排指令顺序重算整个函数的边活跃度和调用点压力，再与原版逐点比较。Win64 计算超出十个非易失 XMM 槽的压力，SysV 按调用易失 SIMD 值的存取代价计权。块频率用于对调用点加权，并扣除提前执行可能增加的算术开销。候选至少改善两个调用点，且在任一建模调用点增加 spill 代价或增加前置块的超限时退出。选择阶段保持原 IR 和 CFG，最终只提交得分最高的调整。

| 同一份 Appell IR | 静态 Spill / Reload | 栈帧 |
| --- | --- | --- |
| Win64 `x86-64-v3` | 20 / 32 → 17 / 29 | 296 字节 |
| Win64 `x86-64-v4` | 21 / 33 → 18 / 30 | 312 字节 |
| SysV 本机 | 59 / 104 → 58 / 103 | 216 字节 |

最热的 `pow` 链少了两个活跃的 Gamma 值。AVX2 Win64 循环内的静态 spill/reload 在两版都为 1 / 18，减少的 3 / 3 位于循环之外。与指数打包合用的 AVX2 机器码为 264 字节栈帧和 15 / 21，组合只做了静态验证。正常 rustc fat LTO 入口也命中两条运算和 16 个建模调用点。

| 真实 Appell A/B，改写相对基线 | 配对耗时变化 | 95% bootstrap 区间 |
| --- | --- | --- |
| Win64 PE 在 Wine 下，墙钟 | −0.64% | −9.57% 至 +23.17% |
| Win64 PE 在 Wine 下，线程 CPU 时间 | 0.00% | −1.06% 至 +2.33% |
| Wine 下，已有 `exp` 复用再加退出表达式改写，线程 CPU 时间 | 0.00% | 0.00% 至 +4.55% |
| Linux SysV，墙钟 | −7.99% | −21.00% 至 +19.36% |

Wine 实验为 24 组交错、每组 120 轮，四个变体的每轮校验一致，样本哈希为 `8871e0382fe46c61`。Linux 以同一份源码抽取的 IR 生成两组对象并链接同一个外壳，比较 25 对运行，每次 38 轮，哈希为 `e8872e9d92154ed5`。两组都使用未独占机器上的固定逻辑 CPU。Gamma 系数改写保留为显式开启。

复现时将上述 PE 命令中的 `pack` 变体替换为 `compression`，管线使用 `function(rssn-loop-exit-compression)`。`both` 使用 `function(rssn-libm-reuse,rssn-loop-exit-compression)`。运行时额外设置 `RSSN_WIN_EXIT_COMPRESSION=1`。要查看 LLVM 的下沉，向同一 Win64 fat LTO Rust 编译追加 `-Cllvm-args=-print-before=indvars -Cllvm-args=-print-after=indvars -Cllvm-args=-print-loop-func-scope -Cllvm-args=-filter-print-funcs=rssn_de_appell`。

## 双指数求积的端点检查

交付 DLL 的独立求积器在 `0x1802b6da9` 调用 `exp`，随后在 `0x1802b6dbe` 检查权重的指数位，在 `0x1802b6dd6` 检查求积距离。下端点、权重为零和上端点检查分别位于 `0x1802b6e0c`、`0x1802b6e17` 和 `0x1802b6e24`。失败分支汇合于 `0x1802b6ed0`。重建的 debug 保留独立的求积器，当前源码的 fat LTO 重建则在 Appell 和小参数 Lerch 内核中得到下端点为零的常量区间。交付 DLL 的实现保留动态端点，提供调用顺序的直接证据。以下改写与量化使用当前真实 RSSN 内核。

设距离为 `d`、两个正有限常量分别为 `C` 和 `H`，节点取 `C*d` 或 `H-C*d`。原先的 `fcmp ugt d, 0` 允许 NaN 进入端点检查。它拒绝的所有数，包括负无穷和两个带符号零，令第一种节点小于等于零，或令第二种节点大于等于上限。沿两条路径既有检查都会去往同一失败块。pass 要求所有失败边在该块每个 PHI 上的输入值相同，检查之间仅能存在可安全推测执行的无内存指令。改写用 `freeze` 固定节点选择器，覆盖原先提前退出时选择器可能为 poison 的路径。受影响的三个分支都清除旧权重元数据。 [LLVM 比较与 freeze 指令](https://llvm.org/docs/LangRef.html#other-operations)

合并后的权重检查接受所有非零有限 binary64。将位模式左移一位丢掉符号，再减一并与 `0xffdfffffffffffff` 作无符号严格小于比较。两个零绕回无符号最大值，无穷与 NaN 落在门限外，两个符号的次正规数都保留。pass 限于 IEEE 次正规模式与普通浮点语义，检查所匹配的算术没有 fast-math 标志。若分母来自可推测执行的 `llvm.exp`，且移动 `exp → fadd → fdiv` 后所有剩余用途依然受支配，就把这三条指令移到权重检查之后。外部 CRT 调用和受约束的浮点运算保持原位。与标量复用合用时先执行复用，Appell 的距离来自 PHI，端点检查仍可合并，已有 `exp` 复用得以保留。 [LLVM 浮点环境](https://llvm.org/docs/LangRef.html#function-attributes) [LLVM exp intrinsic 语义](https://llvm.org/docs/LangRef.html#llvm-exp-intrinsic)

| 同一份 Appell fat LTO IR | 静态浮点比较 | 静态条件跳转 | 静态 Spill / Reload |
| --- | --- | --- | --- |
| Win64 `x86-64-v3` | 23 → 15 | 34 → 22 | 20 / 32 → 20 / 32 |
| Win64 `x86-64-v4` | 23 → 15 | 34 → 22 | 21 / 33 → 21 / 33 |
| Linux SysV | 23 → 15 | 34 → 22 | 59 / 104 → 57 / 95 |

正常的 Linux 和 Win64 rustc fat LTO 插件入口都在 Appell 命中四处，在小参数 Lerch 命中两处。Win64 同时通过 IR 验证并生成 COFF 对象。以下 Linux 计数以 `LD_PRELOAD` 探针读取 `exp` 次数，再扣除每个进程的一轮预热。所有四组源代码路径在改写前后得到相同的位级哈希。

| 真实 RSSN 路径 | 每轮 `exp` 次数 | 每轮输出位级哈希 |
| --- | --- | --- |
| Appell F₁，四组输入 | 14,144 → 13,302 | `6bf2268c6ac1ffeb` |
| Lerch 小参数 | 7,588 → 7,466 | `1ed032a1bc29086b` |
| Lerch 正参数 | 8,349 → 8,349 | `cc602d0d9be79c1d` |
| Lerch 负参数 | 4,164 → 4,164 | `5db0b2774e8993ad` |

Linux 取同一份 Appell IR 的两个版本，交错执行 25 对、每个样本 40 轮。Win64 用 AVX2 PE 与同名 CRT math 导入，在 Wine 10.0 里以四个变体交错执行 24 组、每个样本 120 轮。所有 Win64 计时样本的位级哈希均为 `8871e0382fe46c61`。

| Appell A/B | 配对耗时变化 | 95% bootstrap 区间 |
| --- | --- | --- |
| Linux，端点改写相对基线，墙钟 | +4.88% | −6.55% 至 +18.53% |
| Wine Win64，端点改写相对基线，墙钟 | −4.12% | −28.63% 至 +14.52% |
| Wine Win64，端点改写相对基线，线程 CPU | −2.33% | −2.38% 至 0.00% |
| Wine Win64，已有标量复用再加端点改写，线程 CPU | 0.00% | −2.38% 至 +2.56% |

机器未独占，Wine 线程 CPU 时间约以 10 ms 为单位变化。运行验证确认少执行的调用和位级一致性，计时仍未确认独立收益，该 pass 保留为显式开启。回归测试另对 234 组权重、指数和端点组合执行原版与改写版，覆盖带符号零、次正规数、无穷、安静与信号 NaN。

正常 Rust 构建可用 `RSSN_SPILL_FILTER=rssn_de_`，并在两组之间切换 `RSSN_FP_ENDPOINT_GUARDS=off` 与 `RSSN_FP_ENDPOINT_GUARDS=transform`。复现 PE 时把上文的 `pack` 变体替换为 `endpoint`，分别用 `function(rssn-fp-endpoint-guards)` 与 `function(rssn-libm-reuse,rssn-fp-endpoint-guards)` 生成 `endpoint` 与 `both`，并在运行时设置 `RSSN_WIN_FP_ENDPOINT=1`。在 `llc` 增加 `-stop-after=virtregrewriter` 可输出寄存器分配后的 MIR。

## 用非易失通用寄存器保存浮点位模式

交付 DLL 的求积器保存八个非易失 GPR 和十个 XMM。它在 `0x1802b6e49` 与 `0x1802b6e71` 从自己的栈帧读出捕获指针，随后读取指针指向的 `a,c`。另外四个捕获源从 `r14/r15/r12/r13` 间接读取。重建的同源码 debug 闭包也通过环境指针读取六个捕获，同时 MIR 另有 `t.dbg.spill` 和 `powf` 返回值临时栈槽。fat LTO release 将闭包内联，新的四个重复读取可以直接归因到 SSA 值 `a-1`、`c-a-1`、`-b1`、`-b2` 对应的分配器槽位。另两次折叠到 `vmulsd` 的 `x,y` 读取来自 `%fixed-stack.1/.0`。此时 Appell 入口只保存两个 GPR，出现了非易失寄存器组之间的余量差。交付 DLL 的捕获读取、debug 物化和当前 release 的 spill 来自不同的 lowering，下面的优化数据来自可重复构建的 release。 [Microsoft Win64 寄存器约定](https://learn.microsoft.com/en-us/cpp/build/x64-software-conventions?view=msvc-170)

pass 在两个以上的循环块中寻找重复参与标量 libm 的循环不变 `f64`。候选来自普通浮点算术或可按 Win64 ABI 确认为栈传递的固定签名参数。载体放在所有消费者的公共支配路径上。计算值在定义块中第一次跨调用之前仍可由原浮点值消费，循环内的消费者按需把 `i64` 位模式送回 XMM。每个消费者只生成一次位搬运。若 IR 无法证明值已定义，先用 `freeze` 固定其位模式。还原指令是带 GPR 输入约束的无副作用 `vmovq`，因为只有 `bitcast` 的来回转换会在寄存器分配之前消失。浮点算术和原有调用顺序保持原位。 [LLVM inline asm 约束](https://llvm.org/docs/LangRef.html#inline-assembler-expressions) [Intel 指令手册](https://cdrdv2-public.intel.com/858447/253667-088-sdm-vol-2b.pdf)

模型用 SSA 活跃集合分别计算调用点前后的寄存器组占用。跨调用 GPR 预算为 6，从 Win64 的 8 个非易失 GPR 中预留 2 个给帧指针和 lowering 产生的地址临时值。普通指令使用独立的局部 GP 预算。保留 XMM 的估计只对 128 位及以下占用 SIMD 寄存器的活跃值计数，宽 YMM/ZMM 活跃值仍进入普通 SIMD 压力统计。每个候选的跨调用和局部区间都必须装入预计的 GP 预算。每个消费者只对邻近的一个调用点计收益，且该点剩余的 XMM 活跃度必须大于 10。每次潜在的独立重读计 3 分，可能折叠进算术的内存读取计 2 分，`vmovq` 扣 1 分。按块频率加权后，再扣载体初始化及一个非易失 GPR 的保存恢复。每接受一个候选就更新两组占用和剩余收益，每个函数最多改写四个值。分值用于筛选，实际分配和延迟由 MIR 和运行 A/B 复核。 [LLVM register banks](https://llvm.org/docs/GlobalISel/GMIR.html)

| 同一份 fat LTO IR | 静态 Spill / Reload | `sub rsp` 字节数 | 非易失 GPR 保存数 |
| --- | --- | --- | --- |
| Appell Win64 `x86-64-v3` | 20 / 32 → 15 / 15 | 296 → 232 | 2 → 6 |
| Appell Win64 `x86-64-v4` | 21 / 33 → 16 / 16 | 312 → 248 | 2 → 6 |
| 小参数 Lerch Win64 `x86-64-v3` | 54 / 92 → 52 / 83 | 752 → 744 | 3 → 4 |

Appell 的四个热块各从四次指数重读降到零，改用四条从非易失 GPR 到 XMM1 的 `vmovq`。四个热块各保留两次 ABI 入参的折叠读取，整函数仍有 18 个静态 `pow` 调用。另一个原型把 `x,y` 也送入 GPR，消除了循环内这两次参数重读，入口保存全部八个非易失 GPR，计时没有显示额外收益。Lerch 则选择一个指数 `s-1`，在 IR 涉及 36 个跨调用位置和 15 个使用点。将入口的 GPR 保存纳入帧大小后，Appell 的帧减少 32 字节，小参数 Lerch 的帧大小相同。

| AVX2 PE 在 Wine 下的真实输入 A/B | 配对墙钟变化 | 95% bootstrap 区间 | 配对线程 CPU 时间变化 | 95% bootstrap 区间 |
| --- | --- | --- | --- | --- |
| Appell，跨寄存器组改写相对基线 | +0.69% | −2.36% 至 +9.86% | 0.00% | 0.00% 至 +2.41% |
| Appell，标量 `exp` 复用相对基线 | −4.84% | −7.67% 至 −2.60% | −4.76% | −7.06% 至 −4.76% |
| Appell，标量复用再加跨寄存器组改写 | +0.26% | −1.39% 至 +2.96% | 0.00% | 0.00% 至 +2.41% |
| 小参数 Lerch，跨寄存器组改写相对基线 | +1.11% | −3.83% 至 +4.13% | 0.00% | 0.00% 至 +7.69% |
| 小参数 Lerch，端点改写再加跨寄存器组改写 | +0.27% | −2.89% 至 +11.63% | 0.00% | 0.00% 至 0.00% |

Appell 交错执行 28 组、每样本 120 轮。Lerch 交错执行 24 组、每样本 190 轮。每轮四组真实输入。运行器逐项核对预热，所有相同内核的计时样本哈希分别为 `8871e0382fe46c61` 和 `d4e66064190eff31`。CPU 时间读数约以 10 ms 为单位变化。已确认热区的重读消除，时间收益仍待确认，pass 默认关闭。

复现时用 `RSSN_FP_REGISTER_BANK=transform RSSN_SPILL_FILTER=rssn_de_` 走正常 rustc fat LTO 构建。PE 对照把上面的 `pack` 变体替换为 `bank`，用 `function(rssn-fp-register-bank)` 生成，并设置 `RSSN_WIN_FP_BANK=1`。Appell 的 `both` 用 `function(rssn-libm-reuse,rssn-fp-register-bank)`。对小参数 Lerch，用 `llvm-extract --func=rssn_de_small` 并导出该符号，给运行器添加 `-DRSSN_WIN_LERCH_SMALL=1`，把 `reuse` 替换为 `endpoint`，`both` 用 `function(rssn-fp-endpoint-guards,rssn-fp-register-bank)`，math 导入表包含 `log`、`log1p`、`expm1` 和 `tanh`。

## 跳过舍入至单位底数的 pow

交付 DLL 在 `0x1802b6e6c` 计算 `1-t`，随后在 `0x1802b6e83` 调用 `pow`。接下来两组 `1-k*t` 分别进入 `0x1802b6ea4` 与 `0x1802b6ec1` 的调用。进入这四连调用后，三个减法的结果直接传给 CRT。这与重建 debug 的独立闭包吻合。fat LTO release 内联闭包后，同一结构在四个循环块里重复。用实际输入在 Linux 上计数，每轮 40,412 次 `pow` 中有 9,587 次底数已是精确的 `+1.0`。

Lerch 的另一条线索是标量 `select` 同时引用 `expm1` 与 `llvm.exp`。手动把单用途 `llvm.exp` 放进条件分支后，Win64 和 SysV 的指令都与基线相同。LLVM 23 的 CodeGenPrepare 已将该标量调用下沉到 `select.false.sink`。提高 unswitch 门槛会复制这段循环，试验的 Win64 静态重读从 92 增到 116，未纳入插件。 [LLVM CodeGenPrepare 下沉实现](https://llvm.org/doxygen/CodeGenPrepare_8cpp_source.html)

对 Appell，pass 在一个基本块里查找共享节点 `t`、以 `1-t` 开头的三至四个 `llvm.pow`，要求紧随其后的乘法构成单用途累乘链。后续底数为 `1-k*t`，其中系数和指数在循环中不变。至少两个循环块共享同一组参数，且块频率估计的总访问数至少是初始化次数的八倍。这样一次比较可省一组三连调用，也把创建门限的成本摊到重复节点上。块频率用于衡量摊销，输入分布决定实际命中率。

设 `E` 为包括隐式系数 1 在内的最大 binary64 偏置指数。所有有限系数满足 `|k| < 2^(E−1022)`。当 `0 ≤ t ≤ 2^(968−E)`，实数乘积满足 `|k*t| < 2^-54`。舍入后的乘积最多达到 `2^-54`，此时 `1−2^-54` 恰好在 `1.0` 与其下方相邻浮点数的中点，向最近偶数舍入到 `1.0`。上方的舍入区间更宽。pass 在循环外构造门限的位表示，热路径只需按无符号整数比较节点位模式。负值和非有限节点走原路径。非有限系数、会使门限变成次正规数的系数，以及 NaN 指数都禁用快路径。 [LLVM 默认浮点环境与舍入](https://llvm.org/docs/LangRef.html#floating-point-environment) [LLVM pow intrinsic 语义](https://llvm.org/docs/LangRef.html#llvm-pow-intrinsic)

快路径保留累乘顺序和操作数位置。为保持信号 NaN 的处理，一般累乘初值为 NaN 时走原路径。在 Appell 实际形态中，初值是同一非负节点的第一个 `pow`。其指数的 NaN 检查并入循环外门限，循环内省去额外浮点比较。新分支依赖的未定义值由 `freeze` 处理，检查和回退代码读取同一版本。strict FP、非 IEEE 次正规模式、带副作用的外部数学函数，以及有额外使用者或副作用的指令区域均退出该规则。

Win64 AVX2 的四个热块各得到一条 `vmovq`、一条 `cmpq` 和一条条件跳转。命中时保留首个 `pow` 并绕过其余三个及对应的栈读取。AVX2 和 AVX-512 的静态 `pow` 数均为 18，静态 spill / reload 分别保持 20 / 32 和 21 / 33。AVX2 新增一个 GPR 保存，`sub rsp` 从 296 降至 288 字节，两者计入后栈帧大小相同。Linux SysV 的静态 spill / reload 从 58 / 103 增至 59 / 110。最终产物通过 IR 与机器指令验证，另用 4,704 组输入运行原版与改写版，覆盖舍入边界、带符号零、次正规数和 NaN 回退。

Linux 计数确认，每轮四组真实 Appell 输入的 `pow` 调用从 40,412 降到 30,932，单位底数调用从 9,587 降到 107。两组输出逐位一致，`exp` 调用次数均为 14,144。

| Appell 真实输入 A/B | 配对耗时变化 | 95% bootstrap 区间 |
| --- | --- | --- |
| Linux SysV，线程 CPU，单位底数改写相对基线 | −7.81% | −8.48% 至 −6.94% |
| Wine Win64，线程 CPU，单位底数改写相对基线 | −6.67% | −7.50% 至 −5.53% |
| Wine Win64，墙钟，单位底数改写相对基线 | −4.13% | −14.58% 至 +4.21% |
| Wine Win64，线程 CPU，已有 `exp` 复用再加单位底数改写 | −6.90% | −6.98% 至 −4.76% |
| Wine Win64，线程 CPU，两项同时开启相对基线 | −11.11% | −13.04% 至 −10.87% |

Linux 在同一份 IR 上交错执行 27 对、每个样本 70 轮，线程时间来自 `CLOCK_THREAD_CPUTIME_ID`，样本哈希为 `097d449b72b480b3`。Wine 10.0 使用与交付 DLL 同名的 CRT math 导入，四个 AVX2 PE 变体交错执行 24 组、每个样本 130 轮，样本哈希为 `00830d95dc58007b`。机器未独占。Wine 线程计时约以 10 ms 为单位变化。这里验证的是重建的真实内核，AVX-512 完成了静态机器码对照。编译选项默认关闭。

在正常 rustc fat LTO 构建中设置 `RSSN_POW_UNIT_BASE=transform RSSN_SPILL_FILTER=rssn_de_appell`。PE 复现用上面的构建步骤，把 `pack` 替换为 `unit`，单独管线为 `function(rssn-pow-unit-base)`，组合管线为 `function(rssn-libm-reuse,rssn-pow-unit-base)`，运行时设置 `RSSN_WIN_POW_UNIT_BASE=1`。Linux 计数用 `bazel run //:llvm_benchmark_build -- count-exp --mode pow --out-dir "$study"` 构建探针，输出总 `pow` 次数和单位底数次数，用 `appell 1` 减去 `appell 0` 扣除预热。

### 共享单位底数与跨目标验证

当前 RSSN 的 Appell fat LTO 在 Win64 和 Linux ARM64 都生成四个热 `pow` 块。两块在块内计算 `1-t`，另两块复用支配块在做端点检查之前计算的 `1-t`。交付 DLL 的对应机器码依然把 `1-t` 直接送入第二个 `pow`，调试重建的 integrand 保留独立的 `powf` 调用。IR 新增的共用定义来自当前构建，不能从 DLL 地址反推其 SSA 结构。

守卫现在识别支配慢路径或有其他使用者的首个 `1-t`。原定义及其其他使用保持原样。只有在守卫引入冻结节点时，才在慢路径用同一冻结值重算首个底数。后续 `1-k*t` 的匹配和原累乘顺序保持不变。这个改写与端点守卫均依赖 LLVM 的 IEEE binary64 语义。跨目标测试覆盖 ELF、Mach-O 和 COFF，包括 ARM64、ARM 和 x86 的 32 位目标。

Linux ARM64 使用 rustc 1.100.0 的常规 fat LTO 插件入口和实际 Appell 提取源码。两个块被守卫时，四组输入运行一轮的动态 `pow` 仍为 40,412 次。识别共享 `1-t` 后，四个块均生成守卫，动态 `pow` 由 40,412 降至 30,932，其中底数精确为一的调用由 9,587 降至 107。计数在 QEMU 的 ARM64 用户态用 glibc `pow` 探针获得，并扣除了相同的预热。另固定单核交错运行 15 对，每样本 12 轮，输出哈希均为 `a252b4cc0103f596`，模拟环境中的线程 CPU 时间变化为 **−11.62% [−14.27%, −9.85%]**。

指数打包器在 Win64 增加了同一循环内的单前驱边。边上两端的 `pow` 之间允许纯计算，并以两个块的较小频率估计成对执行频率。首个调用有绕过第二个调用的路径时，对这些执行另计一次可能的参数复制。活跃分析同时参考块内最后使用与块出口集合，覆盖本块中没有局部使用的跨块活跃指数。 [LLVM 块频率定义](https://www.llvm.org/docs/BlockFrequencyTerminology.html)

在只守卫两个块的旧版 Appell IR 上，跨边识别恢复两组指数，包含两个跨边窗口。Win64 AVX2 局部帧从仅单块打包的 272 字节降至 256 字节，静态 spill/reload 从 19/28 降至 18/24。在新版守卫覆盖四个块之后，未提供 profile 的 BFI 估算使含绕过路径的一组配对退出，留下另一组。新版守卫单独启用的静态 spill/reload 为 20/32，另加打包为 19/28，四种代码均保留 18 个静态 `pow`。

对确实命中跨边配对的旧版守卫 IR，另做 20 组、每样本 190 轮的 Wine 交错实验。跨边打包相对仅单块打包的配对线程 CPU 时间为 0.00% [0.00%, +2.08%]。同一次实验中，新版守卫单独相对旧版守卫单独为 −10.42% [−12.24%, −10.10%]。所有输出哈希均为 `f36a3fecb2c104dd`。

Win64 AVX2 的四种构建用同一个 Wine 进程交错执行 24 组，每样本 150 轮，所有输出哈希均为 `08344d98c3ccfef2`。仅新版守卫相对原始内核的线程 CPU 时间为 **−10.53% [−12.82%, −10.26%]**。再启用指数打包为 **−7.89% [−10.26%, −7.89%]**。相对仅新版守卫，进一步打包为 **+1.43% [0.00%, +2.94%]**。这组对照支持单独使用守卫作为 Appell 的实验对照，打包保持可选。

### Appell 单位底数的源码对照

RSSN 的 Appell 积分直接掌握两个尺度系数、正节点范围，以及排除 NaN 指数的输入检查。在 Rust 源码中从 `max(1, abs(x), abs(y))` 的偏置指数构造相同的二进制门限。积分闭包用一次无符号位模式比较判断三个底数是否都舍入到一。命中时返回首个 `powf` 的结果。未命中时按原顺序计算三次幂和累乘。尺度过大或无穷时门限为零，正节点进入原有计算。 [RSSN 实现](../../src/kernels/special_ext.rs)

Linux 上以当前生产内核为输入，使用相同 rustc、fat LTO、AVX2 和四组 RSSN 用例。固定单核交替测量 19 组、每样本 150 轮。源码和旧源码加 pass 每轮都将动态 `pow` 由 40,412 次降到 30,932 次，输出哈希一致。源码改动跨过部分求积步进闭包的内联门槛。Linux 主入口与被保留的闭包合计 66/102 条静态 spill/reload，基线为 57/94。Win64 对应合计 27/30，基线为 20/32。这里的静态计数包含冷路径，CPU 时间比较包含完整调用路径。

| 真实 Appell A/B，线程 CPU 时间 | 相对基线 | 95% bootstrap 区间 |
| --- | --- | --- |
| Linux，只改 Rust | −10.10% | −11.67% 至 −9.26% |
| Linux，旧源码只启用 pass | −10.48% | −12.25% 至 −8.32% |
| Linux，源码相对 pass | −0.76% | −2.72% 至 +1.67% |
| Wine Win64，只改 Rust | −10.53% | −12.82% 至 −8.11% |

Win64 由相同 MSVC AVX2 配置重建并用 Windows `opt.exe` 和 `llc.exe` 验证。在 Wine 下交错运行 17 对，每样本 150 轮。配对墙钟变化为 −8.50% [−14.89%, −6.05%]。两平台计时样本哈希都是 `08344d98c3ccfef2`。Linux 另核对 1,190 组提取的真实积分路径输入，包含接近 1、负的大尺度、次正规门限、负零和无穷系数。1,160 个结果逐位一致，其余 30 个结果均为 NaN。

结果支持把 Appell 的特例留在 Rust。Win64 的原版命中 pass 的四个块，新源码在 Linux 和 Win64 的报告模式均没有生成第二重守卫候选。`rssn-pow-unit-base` 保留为独立、默认关闭的结构化变换，可用于能提供相同 IR 证明的其他调用方。

## 复用 exp-sinh 单位节点的 log/exp 结果

交付 DLL 的 exp-sinh 形态例程 `0x1802bd7f0` 中，首个循环的 `sinh → exp` 位于 `0x1802bd928` 和 `0x1802bd931`。它将原始 `exp` 结果放进 XMM7，保留到 `0x1802bd98b` 的 `log`。随后的 `0x1802bd990` 和 `0x1802bd996` 从帧内读取两项系数，`0x1802bd9ae` 和 `0x1802bd9c1` 分别计算指数分子和分母所需的指数。这里进入 `log` 的是原始 `exp` 值。当前源码重建的小参数 Lerch 尾积分将 exp-sinh 的下限设为 1。其 debug 保留独立的尾积分闭包和 `exp_sinh` 调用，fat LTO release 则在四个循环块中生成 `x = 1 + exp(u)` 和 `exp((s-1)*log(x)-x)`。交付 DLL 缺少匹配 PDB，地址归因到上述机器码形态，单位节点的实验使用当前源码重建。

在 LLVM 的默认浮点环境中，设循环不变的系数 `C` 为有限数。`x` 舍入为 `+1.0` 时，`log(x)` 为 `+0.0`，`C * +0.0` 为带符号零，减去 `x` 后精确为 `−1.0`。pass 在共享支配点用原浮点运算构造该实参，再调用一次同一 `llvm.exp`。热路径比较 `x` 的位模式，通过 PHI 复用该运行时结果。原 `log`、乘法、减法和 `exp` 都保留在回退块。分母上的 `expm1` 及其调用位置保持原状。 [LLVM 默认浮点环境](https://llvm.org/docs/LangRef.html#floating-point-environment) [LLVM log intrinsic](https://llvm.org/docs/LangRef.html#llvm-log-intrinsic)

当 LLVM 能证明节点已定义时，有限系数用 `+1.0` 作匹配键，非有限系数用 `+0.0`。匹配的 `1 + exp(u)` 只能是大于等于 1 的数或 NaN，因而后一个键不会命中。可能未定义的节点用 `freeze` 后再结合显式有限性检查。至少两个链共用同一系数，总访问频率须达到共享初始化频率的八倍。匹配要求链内三段依赖都是单用途且连续，跳过 strict FP、非 IEEE 次正规模式、近似数学及外部 libm。

| 同一份小参数 Lerch IR | 静态 Spill / Reload |
| --- | --- |
| Win64 `x86-64-v3` | 54 / 92 → 55 / 102 |
| Win64 `x86-64-v4` | 61 / 100 → 62 / 111 |
| Linux SysV | 138 / 216 → 137 / 223 |

在 Win64 AVX2 汇编里，寄存器分配将缓存结果放进帧槽，并在单位节点比较之前重读它。原分子系数的折叠重读只在回退块执行。Linux 的逐指令计数在 `rssn_de_small` 内给 LLVM 的 Spill / Reload 位置加零长度符号，以 QEMU 的单指令翻译块记录执行地址，再扣除一轮预热。每轮四组输入执行的 Spill 从 32,848 降到 32,027，Reload 从 61,953 增到 62,358。记录范围限于该内核，CRT 内部指令、GPR 保存和 LLVM 未标为 spill 的栈对象不计入。 [QEMU 指令跟踪选项](https://www.qemu.org/docs/master/user/main.html)

另用两种原型把缓存的恢复限制到命中边。第一种在 Win64 上把缓存保留在非易失 GPR，命中时用 `vmovq` 恢复。静态 Reload 从 102 减到 95，保存的 GPR 从 4 个增至 5 个。第二种将 XMM 载体的重读下沉到命中边，静态 Reload 为 99。两种原型与普通缓存做了 24 组 Wine 交错对照，线程 CPU 变化的中位数都为 0.00%，因此发布的 pass 保留普通缓存。单独把分母乘减移到 `log` 前能降低 Linux 静态 Spill / Reload 到 130 / 212，该尝试的线程 CPU 时间区间仍跨零。

| Linux 真实输入，每轮四组 | `log` | `exp` | `log(1.0)` |
| --- | --- | --- | --- |
| 基线 | 2,092 | 7,588 | 427 |
| 单位节点复用 | 1,667 | 7,167 | 2 |
| 端点检查改写 | 2,092 | 7,466 | 427 |
| 两项一起启用 | 1,667 | 7,045 | 2 |

| 真实小参数 Lerch A/B | 配对耗时变化 | 95% bootstrap 区间 |
| --- | --- | --- |
| Linux SysV，线程 CPU，单位节点复用相对基线 | −1.58% | −4.06% 至 +0.80% |
| Linux SysV，线程 CPU，端点改写再加单位节点复用 | −2.79% | −4.94% 至 −0.17% |
| Wine Win64，线程 CPU，单位节点复用相对基线 | 0.00% | −3.33% 至 0.00% |
| Wine Win64，墙钟，单位节点复用相对基线 | +2.19% | −6.73% 至 +15.33% |

Linux 交错运行 27 对、每个样本 430 轮，哈希为 `352298e232fa8ce9`。Wine 10.0 用与交付 DLL 同名的 CRT math 导入，四个 AVX2 PE 变体交错运行 28 组、每个样本 420 轮，哈希为 `4ef68b2189bcb06e`。机器未独占，Wine 线程计时约以 10 ms 为单位变化。另运行了 260 组浮点边界对照，覆盖舍入临界值、非有限系数、带符号零和信号 NaN。正常 Win64 rustc fat LTO 入口命中四条链，生成的 IR 和 COFF 通过验证。

复现时以 `RSSN_LIBM_UNIT_LOG=transform RSSN_SPILL_FILTER=rssn_de_small` 走正常 rustc。PE 运行器以 `-DRSSN_WIN_LERCH_SMALL=1` 构建，导出 `rssn_de_small`，生成 `base`、`unit-log`、`endpoint` 和 `both`。后两项分别使用 `function(rssn-fp-endpoint-guards)` 和 `function(rssn-libm-unit-log,rssn-fp-endpoint-guards)`，运行时设置 `RSSN_WIN_LIBM_UNIT_LOG=1`。Linux 调用计数用 `bazel run //:llvm_benchmark_build -- count-exp --mode log --out-dir "$study"` 构建探针，输出依次为 `exp`、`log`、`log(1.0)` 次数。

### ARM 目标与 Rust 单位节点复用

当前源码重建的 `rssn_de_small` 在 Linux ARM64 和 ARM32 均保留四条可识别的 `log`/乘减/`exp` 链。这一证明只依赖 IEEE binary64 和 LLVM 数学 intrinsic 的语义，单位节点 pass 现可用于其它目标。由直接前驱的浮点比较或针对同一节点位模式的整数比较已排除单位节点时，pass 用 LLVM ValueTracking 跳过冗余改写。Rust 小正阶 Lerch 尾积分另在源码复用 `x == 1.0` 的完整被积函数值，覆盖依赖 `z` 和 `a` 的分母。对源码改写后的 ARM32 和 Win64 IR 再运行单位节点 pass，bitcode 逐字节一致。

两种 ARM 目标在 QEMU 中对实际 Lerch 的四组输入逐位一致。单独的探针运行扣掉预热，ARM64 与 ARM32 得到相同的动态调用数。ARM32 的基准计时器按公开 `clock_gettime` 的 32 位 `timespec` 字段读取线程时间。

| 每轮四组输入 | `exp` | `log` | `log(1.0)` | `expm1` |
| --- | ---: | ---: | ---: | ---: |
| 基线 | 6,630 | 2,092 | 427 | 3,658 |
| LLVM 缓存分子 | 6,209 | 1,667 | 2 | 3,658 |
| Rust 缓存完整被积函数 | 5,956 | 1,667 | 2 | 3,490 |

ARM64 在 `llc` 验证的提取函数内，静态 spill/reload 分别为 30/44、28/52、18/48。ARM32 对应 21/40、15/42、20/43。每种模拟平台做 18 个交错三元组，所列区间为配对线程 CPU 时间的 bootstrap 95% 区间。ARM64 每样本 45 轮，ARM32 每样本 35 轮。Wine 的 MSVC AVX2 对照为 24 组，每样本 460 轮，并通过四组预热和所有计时样本的输出位比较。

| 平台 | LLVM 分子缓存相对基线 | Rust 完整缓存相对基线 |
| --- | ---: | ---: |
| QEMU Linux ARM64 | −4.70% [−8.42%, −1.25%] | −4.98% [−5.92%, −2.18%] |
| QEMU Linux ARM32 | −3.33% [−4.30%, −2.12%] | −2.09% [−4.55%, +0.005%] |
| Wine Win64 MSVC | 0.00% [−4.00%, 0.00%] | −3.85% [−4.35%, 0.00%] |

在两种 ARM 目标上另各检验 432 组有效输入，覆盖近零的阶数、很小的尺度、较大的负 `z` 及 `z` 接近单位端点的情况。所有有限输出位模式一致。

## 提前截断 exp-sinh 远端节点

交付 DLL 的 exp-sinh 例程在 `0x1802bd985` 检查节点，再于 `0x1802bd98b` 调用 `log`。随后两条折叠内存读取的 `vmulsd` 在 `0x1802bd990` 与 `0x1802bd996` 形成对数项与线性衰减项，减法后到 `0x1802bd9ac` 才依据 `−745` 跳过后面的两个 `exp`。负向节点在 `0x1802bda46` 至 `0x1802bda67` 重复这一结构。当前重建的正阶 Lerch debug 将该截断保留在独立闭包，fat LTO release 则把它内联到四个循环块。交付 DLL 的分母序列与重建版本不同，后者还有按 `z` 选择 `exp` 或 `expm1` 的分支。PE A/B 使用重建内核。

设原阈值为 `−K`，对数分数为 `P − D`，其中 `P = C*log(x)`、`D = A*x`。对有限 `x ≥ 1`，有 `0 ≤ log(x) < 1024`。所以 `P ≤ B`，其中 `B = 1024*max(C,0)`。pass 在循环外构造节点门限 `max(1, fl(fl(2*B + K + 1)/A))`，并要求系数有限、尺度 `A` 有限且为正。求和、除法和热路径乘法的三次舍入使 `D` 的下界至多乘上 `(1−2^-53)^3`。在 `8 ≤ K ≤ 1024` 上，`(1−2^-53)^3 * (2*B + K + 1) > B + K + 0.5`。因此命中后，原先舍入后的分数也严格小于 `−K`。门限溢出或参数无效时，门限为正无穷，已证明有限的节点不会命中。 [LLVM 默认浮点环境](https://llvm.org/docs/LangRef.html#floating-point-environment) [LLVM log intrinsic](https://llvm.org/docs/LangRef.html#llvm-log-intrinsic)

LLVM IR 中，接受有限非零权重 `exp(u)*w` 的支配边证明同一个已定义 `exp(u)` 为正且有限。正阶模式借此复用原有的节点正值分支，在 `log` 之前判定提前退出。原有分数比较继续处理其余节点。两条失败边必须在汇合块的每个 PHI 上选到相同值，且跳过的区域只含纯运算。小参数模式使用 `x = 1 + exp(u)`、`A = 1` 和 `K = 746`。由于 `746 > 1075 ln(2)`，命中时分子舍入为正零。分母、除法和乘权重保留原次序，因而零分母、负分母或非有限分母仍有原有行为。严格浮点、非 IEEE 次正规环境及有副作用的数学调用退出该规则。

| Linux 每轮四组真实输入 | `log` 调用 | `exp` 调用 |
| --- | --- | --- |
| 正阶，基线 → 提前截断 | 4,576 → 3,219 | 8,349 → 8,349 |
| 小参数，基线 → 提前截断 | 2,092 → 1,393 | 7,588 → 6,889 |
| 小参数，已有单位节点复用 → 加提前截断 | 1,667 → 968 | 7,167 → 6,468 |

分母的 `expm1` 调用次数一致。Linux 和 Win64 的真实输入验证比较浮点位模式。边界回归执行 2,600 组输入，覆盖临界值附近的舍入、次正规尺度、带符号零、正负非有限系数及分母。NaN 结果按类别比较，其他值逐位比较。正常 Win64 rustc fat LTO 入口在两条路径各命中四处。

| 同一份 IR | 静态 Spill / Reload | Linux 实际执行的 Spill / Reload，每轮 |
| --- | --- | --- |
| 正阶 Linux SysV | 32 / 67 → 35 / 73 | 21,624 / 57,742 → 20,275 / 56,883 |
| 小参数 Linux SysV | 138 / 216 → 135 / 221 | 32,848 / 61,953 → 31,454 / 61,241 |
| 正阶 Win64 AVX2 | 16 / 30 → 19 / 34 | — |
| 小参数 Win64 AVX2 | 54 / 92 → 55 / 101 | — |
| 正阶 Win64 AVX-512 | 17 / 31 → 18 / 38 | — |
| 小参数 Win64 AVX-512 | 61 / 100 → 62 / 111 | — |

动态计数沿用单位节点实验的带符号汇编和 QEMU 单指令跟踪，扣除预热。Linux 正阶内核自身的执行指令从 254,099 降到 242,550；小参数从 253,602 增至 260,671，后者仍省掉动态 libm 调用。计数范围没有包含 CRT 内部指令。

| Linux 实际输入 A/B，线程 CPU | 配对耗时变化 | 95% bootstrap 区间 |
| --- | --- | --- |
| 正阶，提前截断相对基线 | −4.60% | −5.02% 至 −3.85% |
| 小参数，提前截断相对基线 | −2.98% | −4.71% 至 −1.40% |
| 小参数，已有单位节点复用再加提前截断 | −2.50% | −5.32% 至 −1.24% |

Linux 每组交错运行 27 对，每个样本 440 轮。正阶和小参数的样本哈希分别为 `71128f6fd5911f5a` 与 `97eeb7bada656f60`。

| Wine Win64 AVX2，实际输入 A/B | 配对墙钟变化 | 95% bootstrap 区间 | 配对线程 CPU 变化 | 95% bootstrap 区间 |
| --- | --- | --- | --- | --- |
| 正阶，提前截断相对基线 | −8.53% | −18.76% 至 +11.89% | −3.12% | −6.67% 至 0.00% |
| 小参数，提前截断相对基线 | −0.72% | −21.74% 至 +8.93% | −3.12% | −3.18% 至 0.00% |
| 小参数，已有单位节点复用再加提前截断 | −12.34% | −15.96% 至 −3.52% | −3.12% | −5.56% 至 −2.95% |

Wine 10.0 每条路径运行 24 组交错样本，每个样本 450 轮，正阶和小参数的哈希分别为 `09bb6acd1588b034` 与 `b8f19a49c06be176`。机器未独占，线程 CPU 时间约以 10 ms 为单位变化。AVX-512 仅做静态机器码对照。

复现用 `RSSN_LIBM_TAIL_CUTOFF=transform RSSN_SPILL_FILTER=rssn_de_`。Win64 正阶运行器用 `-DRSSN_WIN_LERCH_POSITIVE=1`，导出 `rssn_de_positive` 并设置 `RSSN_WIN_VARIANTS=2`。小参数运行器用 `-DRSSN_WIN_LERCH_SMALL=1`，导出 `rssn_de_small`，四个 DLL 分别为 `base`、`tail`、`unit-log`、`both`。两种运行器均设置 `RSSN_WIN_LIBM_TAIL_CUTOFF=1`，组合使用 `function(rssn-libm-tail-cutoff,rssn-libm-unit-log)`。

## 负阶 Lerch 的偏移量与跨调用重读

当前 RSSN 负阶积分的分数为 `(order−1)*log(t) − a*t − log_gamma`。Win64 debug 保留 integrand、`exp_sinh` 和 Horner 折叠为独立调用。integrand 通过闭包环境读取 `order`、`a`、`log_gamma`，计算完分数后才调用分母函数和 Horner 折叠。fat LTO 将它们内联，在四个循环块 `.unr-lcssa30/39/48/57` 中把分母和 Horner 调度到 `log` 之前。每个 `log` 返回后，寄存器分配后的 MIR 依次读取 `%stack.4/.1/.3`，对应 `order−1`、入口从 XMM2 接收的 `a` 和 `log_gamma`。AVX2 汇编分别折叠到 `72/48/40(%rsp)` 的两条 `vmulsd` 和一条 `vsubsd`。三个槽都是分配器 spill。以 `.unr-lcssa30` 为例，XMM10 中的原始节点跨过 `log`，XMM13 中的倒数、XMM14 中的 Horner 结果、XMM11 中的求积权重跨过 `log` 和 `exp`。

交付 DLL 中的 `0x1802bd7f0` 例程确认了相近的 exp-sinh 结构和两次系数重读。以 unwind 边界扫描 `log` 后 100 字节内的 `exp`，得到 101 对候选。复核其中具有两乘两减的候选后，仍未找到可同时归因到固定 `log_gamma` 减法和逆向 Horner 的例程。负阶的 IR、MIR、机器码及运行对照来自同一份当前 RSSN 源码的重建。

设 `G` 为第二次减法的循环不变量。令 `B = 1024*max(C,0)`，`H = max(−G,0)`。对有限节点 `x ≥ 1`，舍入后的对数项 `P ≤ B`。当 `C`、`G` 有限且 `A` 有限为正时，把节点门限设为 `max(1, RN(RN(RN(2B + 2H) + (K+1))/A))`。令 `r = 1−2^-53`。门限有限且 `x` 命中时，舍入后的 `D = A*x` 至少为 `r^4*(2B+2H+K+1)`。四次可能向下的舍入来自合并两项余量、加入常量、除法和衰减乘法。若 `D` 有限，首次差的负值至少为 `r*(D−B) > H+K+0.5`，其中 `8 ≤ K ≤ 1024`。再减 `G` 后，最后一次舍入仍使原分数严格小于 `−K`。若 `D` 溢出为正无穷，则原分数为负无穷。非法参数或溢出的门限改用正无穷，走原算式。 [LLVM 浮点环境](https://llvm.org/docs/LangRef.html#floating-point-environment) [LLVM exp intrinsic](https://llvm.org/docs/LangRef.html#llvm-exp-intrinsic)

在负阶零分子规则里，`K=746` 足以使 `exp` 舍入至 `+0.0`。四条链共享循环外的门限。命中路径绕过 `log`、三项重读和 `exp`，再用同一顺序做倒数、Horner 值及求积权重相关的原有乘法。倒数和 Horner 本身仍照常计算。Win64 的分配器每个节点增加一次门限栈读取，未命中时还执行原来的三次读取。要求四则运算有连续的单用途依赖链，并保留原浮点环境。Linux 回归共执行 3,380 行边界对照，其中 780 行新增尺度与偏移组合，覆盖溢出、符号零、次正规数和信号 NaN。

| 同一份负阶 Lerch IR | 静态 Spill / Reload | `sub rsp` 字节数 | 每轮实际执行的 Spill / Reload |
| --- | --- | --- | --- |
| Linux SysV AVX2 | 43 / 100 → 47 / 111 | 328 → 344 | 10,110 / 29,201 → 9,004 / 27,577 |
| Win64 AVX2 | 18 / 36 → 21 / 40 | 392 → 424 | — |
| Win64 AVX-512 | 18 / 35 → 20 / 40 | 408 → 440 | — |

每轮四组真实输入在 Linux 上执行的 `log` 从 1,672 降至 1,111，`exp` 从 4,164 降至 3,603，返回正零的 `exp` 从 841 降至 280，`expm1` 均为 834。QEMU 跟踪只计 `rssn_de_negative` 内且由 LLVM 标记的 Spill / Reload，扣除一轮预热。尽管动态重读下降了 5.56%，该内核自己的指令数由 134,370 增至 140,244。此计数排除 CRT 内部指令。四组输入的结果逐位一致。额外 512 组 Linux 输入以相同产物执行对照，非 NaN 结果逐位一致，NaN 按分类比较。

运行探针还识别出 1,672 个 Horner 节点中有 632 个 `u` 精确等于 1，另有 336 个精确等于 `1/(1−z)`。缓存一个值和用二 lane 同时缓存两个值的原型都维持了输出。Linux 线程 CPU 的配对中位数分别变化 −0.27% 和 +0.47%，两个 95% 区间均跨零。该路径保留现有 Horner 执行方式。

| 负阶 Lerch 实际输入 A/B | 配对耗时变化 | 95% bootstrap 区间 |
| --- | --- | --- |
| Linux SysV，线程 CPU，27 对 × 600 轮 | −5.80% | −7.75% 至 −3.75% |
| Linux 正常 rustc fat LTO 运行器，墙钟，27 对 × 600 轮 | −9.70% | −20.29% 至 +5.20% |
| Wine Win64 AVX2，线程 CPU，20 对 × 2,100 轮 | −7.95% | −8.93% 至 −5.41% |
| Wine Win64 AVX2，墙钟，同一批样本 | −7.04% | −16.03% 至 +10.64% |

Linux 和 Wine 样本哈希分别为 `c68f6e29c95dc76a` 和 `1cb96381fbf0300d`，对应不同的轮数。Wine 10.0 用与交付 DLL 同名的 CRT math 导入，线程计时约以 10 ms 为单位变化。运行器均固定在单个允许的逻辑 CPU 上，机器未独占。Linux 和 Win64 的独立运行器都从各自同一份 fat LTO IR 链接基线与改写版，双方都把冷 `panic_bounds_check` 调用替换为 trap 以免引入 Rust panic 运行时。正常 Linux rustc fat LTO 构建保留原 panic 路径，复核了调用计数和输出。正常 Win64 rustc fat LTO 构建也保留原 panic 路径，在负阶、正阶和小参数上各命中四处，通过 IR 校验并生成 COFF。AVX-512 仅用于静态对照。

复现时设置 `RSSN_SPILL_FILTER=rssn_de_negative`，切换 `RSSN_LIBM_TAIL_CUTOFF=off/transform`，每种模式使用独立的 `--out-dir`，使用上文 Rust 运行器的 `negative` 子命令。独立 PE 对照用 `--kernel lerch-negative` 提取并导出 `rssn_de_negative`，同时保留积分器和实际 Gamma 实现。Linux 调用总数可用 `count-exp --mode log` 探针复核。

```bash
bazel run --define=rust_target=x86_64-pc-windows-msvc //:llvm_c_benchmark_build -- \
  win-appell --out-dir "$study/pe-negative" --kernel lerch-negative \
  --variant base=verify --variant 'tail=function(rssn-libm-tail-cutoff),verify'
bazel run //:llvm_benchmark_compare -- run "$study/pe-negative/appell-suite.exe" --runner wine \
  --env RSSN_WIN_LIBM_TAIL_CUTOFF=1 --env RSSN_WIN_VARIANTS=2
```

## 合并有限非零类别的区间检查

交付 DLL 的 exp-sinh 例程在 `0x1802bd947` 形成求积权重，随后用 `bzhiq` 去掉符号位，又在 `0x1802bd955` 至 `0x1802bd97a` 分别构造并合并类别条件。当前 Win64 fat LTO 在三个 Lerch 内核中各产生四个 `llvm.is.fpclass.f64(x, 615)` 热分支。带状 LU 的主元检查复现同类问题。重建的 debug 保留零检查与 `is_finite` 调用，release 得到掩码 `579`，因为主元来自取绝对值后的比较。基线 Win64 release 在检查主元的块中从 `vmovq` 到条件跳转用 14 条指令，改写后用 5 条。包括立即数物化和分支。

对 binary64 的位模式 `b`，无符号 `((b << 1) - 2) mod 2^64 < 0xffdffffffffffffe` 恰好表示有限非零。移位删除符号位。减 2 将最小次正规数移到区间起点，将两个零绕到区间上方。两个无穷、所有安静与信号 NaN 也在区间外。`llvm.is.fpclass` 的分支层语义允许用整数表示完成判断，x86 将移位和减法折叠为 `lea`。pass 处理该集合及其补集。若 LLVM 先删除了仅在某个符号出现的类别，则通过 `KnownFPClass` 和支配点事实验证缺失类别。对于循环携带的 `phi`、`select`，遍历最多 24 个生产者，对每个叶子使用与迭代无关的证明。这样可以处理带状 LU 的主元累积，同时保留无法证明符号的掩码。 [LLVM is.fpclass 语义](https://llvm.org/docs/LangRef.html#llvm-is-fpclass-intrinsic)

| 已开启尾端截断的同一份 Win64 AVX2 IR | 静态 Spill / Reload | `sub rsp` 字节数 |
| --- | --- | --- |
| 负阶 Lerch | 21 / 40 → 20 / 39 | 424 → 408 |
| 正阶 Lerch | 19 / 34 → 19 / 37 | 256 → 264 |
| 小参数 Lerch | 55 / 101 → 55 / 101 | 752 → 760 |

Linux 负阶 Lerch 在现有尾端截断之上，线程 CPU 时间下降 1.97%，95% bootstrap 区间为下降 3.27% 至 1.00%。27 对样本、每样本 700 轮，哈希为 `26725a1f2e8d0d1b`。该内核的静态 Spill / Reload 维持 47 / 111。同一真实内核的 Win64 AVX2 PE 在 Wine 10.0 中线程 CPU 时间下降 2.41%，区间为下降 4.55% 至 1.19%。墙钟变化为 +2.89%，区间为 −9.04% 至 +21.73%。18 对样本、每样本 1,600 轮，哈希为 `2c325b399be58a33`。带状 LU 用真实 `BandLu::factor` 运行三种带宽和包含行交换的输入，共 23 对、每样本 200 轮。Linux 墙钟的配对变化为 −0.47%，区间为 −8.66% 至 +4.67%，哈希为 `0e786d6310ad3319`。LU 同时完成了 Win64 debug 与 release 重建、IR 验证及 COFF 生成。交付 DLL 中的地址归因限于已经解码的求积例程。

对照原型将尾端门限从 `double` 改为精确排序的 `i64` 位模式。Win64 负阶四个站点从每次重读门限变成 `vmovq` 加寄存器间 `cmpq`，静态帧分配从 424 降到 408 字节。单独开启该原型以及和类别区间改写合用的 Wine 线程计时区间均含零。Linux 负阶原型的静态 Reload 从 111 增到 118。现有尾端门限保留原浮点比较。报告与改写入口分别是 `function(rssn-fp-class-interval-report)` 与 `function(rssn-fp-class-interval)`。正常 rustc 以 `RSSN_FP_CLASS_INTERVAL=report` 或 `transform` 启用。

## 负阶 Lerch 源码复用单位分母

负阶积分的尾端可利用分母公式本身。对 `0 < z < 1`，`t ≥ 38` 时 `expm1(-t)` 舍入到 `−1`。`z ≥ 1/2` 时 `1−z` 精确，`z < 1/2` 时第一次减法的绝对误差至多 `2^-54`，两种情形的最终分母都舍入到 `1.0`。对 `z ≤ 0`，将门限提高到 `38 + b`，其中 `b` 是 `|z|` 的非负二进制指数上界。因为 `2 < e`，实数 `|z|exp(-t)` 小于 `exp(-38)`，它和 `2^-54` 之间的余量覆盖 `exp` 与乘法舍入。若 `38+b < 746`，次正规 `exp` 的绝对舍入误差乘上 `|z|` 后也在该余量内。门限超过 746 时直接采用 746，`exp(-t)` 已舍入到正零。

这一门限和 `P(1)` 在 `lerch_negative_order` 的 Rust 源码中计算一次。命中后可直接使用单位分母和按原次序求出的 `P(1)`。四组真实输入每轮减少 315 次 `exp` 和 316 次 `expm1`，其中 `log` 次数一致。只采用 `t ≥ 746` 的源码原型使 Linux 线程 CPU 中位数下降 2.50%，区间为 −5.59% 至 +0.61%。采用依赖 `z` 的门限后，线程 CPU 下降 4.69%，区间为 −6.40% 至 −2.19%。同时开启尾端截断与类别区间两个 pass 后，源码改动仍额外节省 3.73%，区间为 −5.32% 至 −2.29%，且两个 pass 均命中四处。三组线程计时各为 27 对样本、每样本 700 轮，位级哈希均为 `26725a1f2e8d0d1b`。

Win64 AVX2 源码 A/B 采用重建的 MSVC release COFF 和 Wine 10.0。18 对样本、每样本 1,500 轮，线程 CPU 时间下降 4.55%，区间为 −7.32% 至 −3.44%。墙钟下降 4.30%，区间为 −7.75% 至 −0.24%。样本哈希为 `63892b75263788ba`。Linux 另对旧、新源码各自带和不带上述两个 pass 的四种产物运行 1,055 组输入，覆盖 `z` 在 `−1.0` 附近、极大负值、接近 1 的正值及次正规的 `a`。四种产物的 1,005 个非 NaN 结果逐位一致，其余 50 个结果均为 NaN。对照使用从当前源码抽取的 RSSN 内核及相同输入，独立链接时双方的冷 panic 路径都用 trap 代替。

## CPU lane 后端的整数幂

交付 DLL 的 `0x176d4b0` 入口用指令数乘以 128 分配双精度 scratch，外层按最多 128 行处理，指令记录间距为 16 字节。输入和输出寄存器按 1,024 字节间距寻址。八连 `pow` 和各数学函数的分发与 `src/backend/device/cpu.rs` 中的 `run_range_chunk`、`exec_block` 及其内联 slice helper 对应。整数幂臂从 `0x176edd9` 生成无符号指数幅值。自 `0x176ee30` 起的八行展开逐行复制该幅值、掩取最低位、右移，用标量 `vmulsd` 更新基数和累乘值。整组之外的短尾块同样执行标量算法。工作区较早的 Linux debug 路径归档 `librssn-e2433476288a0385.rlib` 仍有命名的 `run_range_chunk`，对应 ELF 对象也在连续行重复 `shr` 和标量乘法。该归档含优化过的函数体。另用同一份提取源码生成 Win64 O0 与 O3 对照，O0 保留独立的 `int_pow` 调用，O3 把递推内联进八行循环。交付 DLL 和归档无法按统一的链接、profile 或机器特性归因。

指数是 `KernelOp::PowI` 的一个立即数。同一条操作的每一行都有自己的基数与乘法结果，可以共享指数控制流。Rust 的 `powi_lanes` 原型（现在共用 `eval_powi_slice`）先求最低置位的位置，执行该位之前的平方，用该次幂直接初始化累乘值。后续置位按从低到高的次序累乘。负指数在最后才取倒数，零指数直接得到 `1.0`。这种行间调度保留每行的运算次序，第一项省去乘 `1.0` 的操作也保持原语义。凑不满八行的尾部调用原来的 `int_pow`。LLVM 23 能在这个形态中把各行的平方与累乘生成为 `vmulpd`。Win64 AVX-512 对照的整份分发器中，静态 spill/reload 标注从 133/533 变为 132/567，收益主要来自整数幂热臂里共享的整数控制流和 SIMD 运算。相关前提属于 CPU 后端的操作语义，因此由 Rust 表达。通用循环互换或 unswitch 要求满足不同的循环结构与成本条件。 [LLVM 循环 unswitch](https://llvm.org/doxygen/classllvm_1_1SimpleLoopUnswitchPass.html)

### 运行对照

基准直接提取生产 `exec_block`、算术 helper 及仓库设备测试中的 `smooth` 和 `everything`，以相同的行块和列布局执行。scratch 预分配，计时包括分发、数值运算和结果列复制。程序生成和结果哈希位于计时区外。4099 行包含完整 lane block 和尾部，11 行包含一个八行组与三行标量尾部。两份程序在同一逻辑 CPU 上交错计量 23 对线程 CPU 时间。所有配对输出哈希一致。

| Linux x86-64-v3 输入 | 行数 | 配对耗时变化 | 95% bootstrap 区间 |
| --- | ---: | ---: | ---: |
| 仓库测试 `smooth` | 4099 | −17.71% | −19.38% 至 −15.46% |
| 仓库测试 `everything` | 4099 | −10.28% | −11.64% 至 −9.33% |
| 仓库测试 `smooth` | 11 | −11.58% | −12.10% 至 −9.89% |
| 仓库测试 `everything` | 11 | −6.46% | −7.34% 至 −5.42% |
| 两次整数幂的隔离对照 | 4099 | −82.41% | −82.51% 至 −82.02% |
| 乘法与加法的对照 | 4099 | −2.51% | −3.88% 至 +0.91% |

Win64 A/B 使用同一份提取源码，以 `x86_64-pc-windows-gnullvm` 和 AVX2 构建，在 Wine 下按交替顺序比较 19 对线程 CPU 样本。交付 DLL 采用的 CRT、链接方式和目标 CPU 继续单独归因。每个输入的新旧版本输出哈希一致。

| Wine Win64 AVX2，4099 行 | 配对耗时变化 | 95% bootstrap 区间 |
| --- | ---: | ---: |
| 仓库测试 `smooth` | −8.16% | −9.80% 至 −7.84% |
| 仓库测试 `everything` | −3.03% | −3.12% 至 −2.78% |
| 两次整数幂的隔离对照 | −80.77% | −81.48% 至 −78.57% |

比较相同算法的四行、八行和十六行组时，十六行组在整数幂隔离对照中的耗时比八行组少 24.86%。在 11 行 `smooth` 中则多 8.42%，因为这个长度只走标量尾部。4099 行 `everything` 的宽度差异区间均跨零。采用八行沿用 CPU 后端原有的展开粒度，并覆盖短批次。

定向正确性检查比较 18 个指数、0 到 260 行及各行的 binary64 位模式，包括有符号零、次正规数、无穷、NaN 和 `i32::MIN`。同样的检查用 AVX2 Win64 可执行程序在 Wine 中通过。生产测试另覆盖 CPU 设备入口和尾块。基准范围是提取后的 lane 执行器，未测量从 FFI、设备注册、上传到完整计算的端到端时间。

在仓库根目录复现。将 `old_cpu` 设为修改前的 `src/backend/device/cpu.rs` 副本。

```bash
dir=$(mktemp -d)
bazel run //:llvm_benchmark_build -- lanes --source "$old_cpu" \
  --out-dir "$dir/base" --rustflag=-Ctarget-cpu=x86-64-v3
bazel run //:llvm_benchmark_build -- lanes --out-dir "$dir/final" \
  --rustflag=-Ctarget-cpu=x86-64-v3
for version in base final; do
  bazel run //:llvm_benchmark_compare -- run "$dir/$version/lanes" -- check
done
for sample in "smooth 460" "everything 180"; do
  read -r case rounds <<< "$sample"
  bazel run //:llvm_benchmark_compare -- --cpu "$cpu" pair \
    "$dir/base/lanes" "$dir/final/lanes" --pairs 23 -- "$case" 4099 "$rounds"
done
```

## Win64 跨调用压力诊断

`function(rssn-spill-pressure-report)` 现在可在 Win64 上运行。它用 SSA 活出减去调用结果，分别统计 GP、可以放进 XMM 的值、目标上原生合法的宽 SIMD 值及尚未分类的 SIMD 部件。Win64 ABI 的上界是 8 个非易失 GPR（不计 RSP）与 10 个非易失的低 128 位 XMM；宽 YMM/ZMM 的高位属于易失状态。报告按块频率估计乘以容量欠额选择站点，列出被调函数和 IR 结果名。显式 ABI 调用与可能降为 libm 的标量 intrinsic 分开排名。需拆分而无法确定保存方式的 SIMD 值保留为未分类。 [Microsoft x64 调用约定](https://learn.microsoft.com/en-us/cpp/build/x64-calling-convention?view=msvc-170)

| 重建的 Win64 IR | 选中的显式调用 | GP / XMM 大小 / 原生宽 SIMD | 超容量的显式调用点 |
| --- | --- | --- | --- |
| Appell F₁ 求积 | `cosh`，`%i187` | 1 / 15 / 0 | 12 / 12 |
| 负阶 Lerch | `sinh`，`%378` | 10 / 12 / 0 | 12 / 13 |
| CPU lane 分发 | `atan2`，`%_0.i.i1750` | 220 / 0 / 85 | 24 / 87 |

Appell `bb182` 的 Machine IR 使用 `csr_win64`，在选中的 `cosh` 返回后用 XMM13 携带结果穿过下一次 `cosh`。常量池读取和后端重物化也能缩短实际区间。原始 DLL 的 Appell 样本在入口保存 XMM6–XMM15，并有额外的函数内槽位访问。静态容量报告与槽位数仍需通过机器码逐点对应。

CPU 分发器的首个 `atan2` 在 DLL 的 `0x176db80`，同一 unwind 区间有 15 个静态 `atan2` 调用点。重建 IR 的 `%_0.i.i1750` 前活跃的 305 个受跟踪 SSA 值中，有 283 个定义于分发循环的 `bb4.lr.ph`。比较、地址表达式和广播占据大部分，后端可以重物化或重新安排。这个分解把大的 SSA 估计归因到循环前置块，也解释了为何它不能直接预测 305 个机器栈槽。

用同一输入分别运行 verifier 与报告模式，CPU lane、Appell 和 Lerch 三份 Win64 bitcode 逐字节一致。Linux SysV 的 RSSN 小参数求积 IR 经新旧插件改写后也逐字节一致。用 `llc -stop-after=virtregrewriter -verify-machineinstrs` 可把具体报告站点和分配后的 MIR 对齐。DLL 入口可用 `python3 tools/llvm-optimizer/analysis/pe_math_calls.py target/rssn.dll --symbol atan2` 找到导入候选，再按 unwind 区间复核完整反汇编。

## Tape 列执行器复用整数幂批处理

通过同一个 `atan2` 导入定位到另一个 15 调用点区间 `0xe41990..0xe48a93`。该区间按 `min(rows, 256)` 计算 stride，使用 24 字节指令间距，与 `Tape::eval_columns_range` 一致。整数幂臂从 `0xe431d9` 起的组内八条标量路径各自移位相同的指数。工作区 Linux debug 路径归档中有名为 `eval_columns_range` 的符号，同样按 256 行取上限，整数幂臂也逐行执行移位。CPU lane 已验证的八行递推现在由两条路径共用 `eval_powi_slice`。Tape 的短输入仍按 NaN 补齐，零指数结果仍为 `1.0`。使用 LLVM 23.1.5 重建 Win64 MSVC AVX-512 机器码后，该函数内的静态 `vmulsd` 从 44 条降到 28 条，`vmulpd` 从 8 条增到 28 条。两份重建均通过 IR 与机器指令验证。

A/B 提取生产 `eval_columns_range` 和算术实现，使用单线程、可复用的 scratch，在两端执行相同的设备测试程序与输出列复制。Linux 和 Wine Win64 AVX2 各做 17 对交错样本，每个用例取 4099 行。每个平台的新旧结果哈希相同。

| 输入 | Linux 线程 CPU 变化（95% 区间） | Wine Win64 线程 CPU 变化（95% 区间） |
| --- | ---: | ---: |
| 仓库测试 `smooth` | −17.08% [−19.65%, −16.21%] | −5.77% [−10.00%, −3.92%] |
| 仓库测试 `everything` | −6.06% [−7.76%, −4.60%] | 0.00% [−2.94%, 0.00%] |
| 两次整数幂的隔离对照 | −82.37% [−82.91%, −82.18%] | −85.42% [−86.36%, −84.78%] |
| 乘法与加法对照 | +0.05% [−1.57%, +1.71%] | — |

以同样的 CPU lane 基准比较共享 helper 前后的实现，四个用例的区间均跨零。独立正确性检查覆盖 0 到 520 行、12 个指数、两种 NaN、带符号零、次正规数、无穷和缺失输入行，Linux 与 Wine 均通过。真实设备测试补充了 `CpuDevice` 和 `Tape::eval_columns` 入口的比较。上述计时和运行检查使用提取的解释器，完整 RSSN 测试套件未运行。

复现时将 `old_tape` 设为修改前的 `src/backend/tape.rs` 副本。

```bash
dir=$(mktemp -d)
bazel run //:llvm_benchmark_build -- tape --source "$old_tape" \
  --out-dir "$dir/base" --rustflag=-Ctarget-cpu=x86-64-v3
bazel run //:llvm_benchmark_build -- tape --out-dir "$dir/current" \
  --rustflag=-Ctarget-cpu=x86-64-v3
for version in base current; do
  bazel run //:llvm_benchmark_compare -- run "$dir/$version/tape" -- check
done
for sample in "smooth 1750" "everything 540" "powi 2400"; do
  read -r case rounds <<< "$sample"
  bazel run //:llvm_benchmark_compare -- --cpu "$cpu" pair \
    "$dir/base/tape" "$dir/current/tape" --pairs 17 -- "$case" 4099 "$rounds"
done
```

## 循环前置块的分支专用计算

原始 Windows DLL 中，CPU 行批处理区间从 RVA `0x176d4b0` 起，入口执行 `sub rsp, 0x238` 并保存 XMM6–15。归档于 `target/debug` 的 Linux `run_range_chunk` 对应区间执行 `sub rsp, 0x278`。两者的 ABI、优化和内联条件有差异。DLL 的第一组八次 `atan2`（首调用 RVA `0x176db80`）直接从 lane 缓冲区读两个参数并写回结果，该组内没有显式的 RSP 相对访问。紧随其后的余数路径在调用后从常量池重建三个向量广播。这个站点说明 LLVM 已经消化了部分 IR 活跃区间，局部 A/B 必须分别检查整个函数的分配、热臂机器码和实际耗时。

将 LLVM 23.1.5 的现有 `sink` 追加到完成 O3 优化的 Win64 MSVC IR 上，选中 `atan2` 站点只减少三个 GP 值，宽 SIMD 数量相同。对已有后端运行 `-machine-sink-cycle-limit=0` 和 `512`，对照版本的帧和标注计数相同。LLVM 的通用 IR 下沉与机器循环下沉仍在完整构建中运行。[LLVM IR Sink](https://llvm.org/docs/doxygen/Sink_8cpp_source.html) · [LLVM MachineSink](https://llvm.org/doxygen/MachineSink_8cpp_source.html)

### 选择位置与评估

按前置块内从使用到定义的顺序规划，使地址、掩码和广播的依赖链可以一起下沉。每条指令只移动一次，候选位置支配所有使用点。PHI 的使用归到其前驱边。位置必须在同一直接循环内，避开循环头、嵌套循环和异常处理入口。操作从整数、地址计算和纯向量构造中挑选，使用 LLVM ValueTracking 证明输入无 undef。没有 CFG 改写、内存访问移动或 `freeze` 移动。数值算法和调用约定保持原状。

先要求目的块的 BFI 估计频率不超过原块。随后将整批位置作为虚拟调度，重新求解包括循环回边和 PHI 边在内的 SSA 活跃集合。对每个普通 Win64 ABI 调用或可能降到标量 libm 的 intrinsic，按块频率加权求和，成本来自 GP 和 XMM 保存容量欠额与原生宽向量的保存恢复估计。新调度只有在总分严格下降时才应用。这能扣除因下沉而新增的跨调用输入，也能避免把同一依赖链的收益重复归到每条指令。尚未分类的类型、常量池重物化和实际 spill 成本仍要通过 MIR 及运行实验核实。Tape 列解释器和重建的 Appell 求积在这次扫描中均无改写。

### 代码生成与运行对照

当前源码提取的 `exec_block` 在 MSVC AVX-512 上移动了 270 条指令。选中 `atan2` 的 GP／原生宽 SIMD 估计从 `235 / 85` 降到 `66 / 44`。用 DLL 调查时的较早源码独立重建，253 条位置变更将相同类别的估计从 `220 / 85` 降到 `67 / 44`。这个较早版本在 `virtregrewriter` 后，8 字节 spill slot 从 72 个降到 59 个，64 字节槽从 43 个降到 22 个。当前与较早两份重建的首组八次 `atan2` 都保留了直接读写 lane 缓冲区的形态。机器码差异主要分布在整个分发循环及其他分支。

| 输入与配置 | `sub rsp` 字节 | 整函数静态 Spill / Reload 注释 |
| --- | ---: | ---: |
| 当前源码，MSVC AVX-512 基线 | 3416 | 132 / 567 |
| 当前源码，MSVC AVX-512 加 LLVM `sink` | 3416 | 132 / 571 |
| 当前源码，MSVC AVX-512 加插件 | 1944 | 100 / 479 |
| 较早源码，MSVC AVX-512 基线 | 3400 | 133 / 533 |
| 较早源码，MSVC AVX-512 加插件 | 1944 | 98 / 443 |
| 当前源码，GNU Win64 AVX2 基线 | 472 | 41 / 395 |
| 当前源码，GNU Win64 AVX2 加插件 | 472 | 41 / 392 |

较早 AVX-512 重建还试过把计算推迟到块内首次使用前。相同的 253 条移动产生了 1960 字节帧和 100 / 466 条注释，收益低于从块头开始的方案。LLVM IR verifier 和机器指令 verifier 均通过，报告模式与基线 bitcode 逐字节相同。

当前 CPU 不支持 AVX-512。运行 A/B 使用交叉链接的 GNU Win64 AVX2 可执行文件，在 Wine 上固定单核，输入 4099 行，交替运行 17 对。所有用例新旧输出哈希一致，整数幂特殊值及尾行的验证也逐位通过。Wine 线程 CPU 时间以 10 ms 为单位变化，短对照的区间较宽。四组区间均包含零，故插件默认关闭。

| 输入 | 配对线程 CPU 时间变化（95% 区间） |
| --- | ---: |
| RSSN 设备测试 `everything` | +1.28% [−1.30%, +1.39%] |
| RSSN 设备测试 `smooth` | −2.22% [−4.17%, 0.00%] |
| 乘法与加法对照 | 0.00% [−12.50%, 0.00%] |
| 两组整数幂对照 | −8.33% [−8.33%, 0.00%] |

可用当前工作树提取真实执行器并复现三组 MSVC 静态对照。将 `opt`、`llc` 和 `plugin` 设为同一版本 LLVM 23 的工具及插件路径。

```bash
dir=$(mktemp -d)
bench=tools/llvm-optimizer/bench/device-lanes
python3 "$bench/extract.py" "$dir/device.rs"
rustc --target x86_64-pc-windows-msvc --crate-type lib --crate-name rssn_lane \
  --edition=2024 -Copt-level=3 -Ccodegen-units=1 -Cpanic=abort \
  -Ctarget-cpu=x86-64-v4 --emit="llvm-ir=$dir/base.ll" "$dir/device.rs"
"$opt" -passes='function(sink),verify' "$dir/base.ll" -o "$dir/llvm-sink.bc"
RSSN_SPILL_MODE=off RSSN_SPILL_FILTER=exec_block "$opt" \
  -load-pass-plugin="$plugin" \
  -passes='function(rssn-loop-preheader-pressure,rssn-spill-pressure-report),verify' \
  -pass-remarks-analysis='rssn-(loop-preheader-pressure|spill-pressure)' \
  "$dir/base.ll" -o "$dir/plugin.bc"
for input in base.ll llvm-sink.bc plugin.bc; do
  "$llc" -O3 -verify-machineinstrs "$dir/$input" -o "$dir/$input.s"
  "$llc" -O3 -verify-machineinstrs -stop-after=virtregrewriter \
    "$dir/$input" -o "$dir/$input.mir"
done
```

运行对照从同一份提取源码构建 `bench/device-lanes/main.rs`，设 `LANE_SOURCE` 为提取文件的绝对路径。使用 `x86_64-pc-windows-gnullvm` 的 Clang／MinGW 链接工具链及 `-Ctarget-cpu=x86-64-v3`。分别设置 `RSSN_SPILL_MODE=off RSSN_SPILL_FILTER=exec_block RSSN_LOOP_PREHEADER_PRESSURE=off` 或 `transform`，两组用同一 `-Zllvm-plugins` 和 `-Zverify-llvm-ir`。执行 `check`，然后按用例分别使用 500、1700、6500 和 6500 轮，记录每次程序输出的时间和哈希。

## 标量调用余数链上的提前重物化

原始 DLL 的 CPU 分发器调用 `atan2` 后，从 RVA `0x176dcb2` 起连续执行三条从常量池读取的 `vbroadcastsd`。继续处理第二行的路径在 RVA `0x176dd06` 又执行同组广播，中间没有使用这些 YMM 值。Linux debug 路径归档里，`run_range_chunk` 的 `atan2` 与 `pow` 余数调用分别位于节内 `0x2802` 和 `0x2af2`，GOT 重定位确认符号。两条路径从源缓冲区读取参数，调用后直接写回并跳到余数循环的条件检查，未出现 DLL 里这一串广播。两种产物的 ABI、调用形式和优化设置有差异。下面的 MIR 结论基于同一 MSVC 重建和同一 LLVM 23.1.5 后端阶段。

`virtregrewriter` 后，重建的 `atan2` 与 `pow` 余数路径各有七个静态调用点。每次调用后生成五条常量池广播、一个宽向量清零和一个传入栈参数读取（`%fixed-stack.1`）。继续处理下一行的一侧没有广播寄存器活入，下一次 `csr_win64` 调用前也没有使用它们。退出到共同分发器续点的一侧才需要这些值。七个块中的局部栈访问另有余数条件和索引，不能按广播一起移走。LLVM 的 post-RA MachineSink 主要下沉寄存器复制，这些常量池广播是重物化产生的目标机器指令。 [LLVM MachineSink 源码](https://llvm.org/docs/doxygen/MachineSink_8cpp_source.html)

`analysis/mir_call_remat.py` 对这种结构做静态筛查。它使用物理寄存器别名、后继活入、调用 regmask 和后继块下一次调用前的指令，并分开记录固定传入栈与局部栈的读取。对两条余数链各报告 30 个单侧广播和 6 个传入栈参数读取候选。最后一次调用直接前往退出块，所以它的五条广播计入调用后总量。在 worktree 中复现时，先按上一节生成基线 MIR：

```bash
python3 tools/llvm-optimizer/analysis/mir_call_remat.py "$dir/base.ll.mir" \
  --function exec_block --callee atan2
python3 tools/llvm-optimizer/analysis/mir_call_remat.py "$dir/base.ll.mir" \
  --function exec_block --callee pow
```

一个手工 MIR 对照把两臂各七个余数出口重定向到同一共享边块，在那里一次性恢复七个值。用 `llc -O3 -verify-machineinstrs -start-after=virtregrewriter` 对原始和调整后的 MIR 走相同后端。整函数静态 RIP 相对广播从 95 条降至 30 条，栈帧同为 3544 字节，spill/reload 注释也同为 132/567。观察到的是常量读取和重物化位置的变化。通过正常完整流水线编译的 3416 字节帧来自不同的后端起点，不能与这个对照直接比较。当前主机无法运行 AVX-512 实验产物。另在同一工作负载上构建了可运行的 AVX2 对照。

IR 插件仍可先调整使值跨调用的定义和使用。本轮把前置块下沉的旧寄存器个数目标换成保存代价。GP 和 128 位保存容量缺口按窄保存恢复对计价，原生宽 SIMD 使用 TTI 的 load/store 代价和按 16 字节分档的传输权重。计分在 BFI 加权后使用。寄存器容量报告保留单位为寄存器的字段，以便同 MIR 对照。在 AVX-512 测试中，以相同控制流分别保留一个 YMM 和一个 ZMM，对应的估计开销为 192 和 256；已复现的真实 MSVC 候选仍选择 270 条下沉。这个 IR 模型没有把目标指令级的广播计为实际 spill。

GNU Win64 AVX2 的 MIR 在同一模式里还有真实 spill reload。两条余数链各有 12 个仅退出边需要的常量池广播和 12 个仅退出边需要的局部 spill reload。用一个共享退出块恢复两条链共同需要的两个向量广播、一个 32 字节 spill 和一个 8 字节 spill。同一可执行文件的 MIR 基线与对照均从 `virtregrewriter` 后继续，且对两组都移除 LLVM 23 MIR 打印器留下的一个无法重新解析的 `!range` 地址。之后用相同的 LLVM 23、Clang 和 MinGW 链接。整函数的 RIP 相对广播为 54 → 28，spill/reload 注释为 41/395 → 41/369，帧分配均为 472 字节。

两组程序均通过整数幂和尾行的逐位验证，下表每个用例的输出哈希也一致。Wine 固定单核，以相反顺序交替跑 11 对。127 行包含 7 行标量余数，128 行无余数。4099 行来自常规设备实验，`smooth` 不执行这两条调用臂。

| 输入 | 配对线程 CPU 时间变化（95% 区间） |
| --- | ---: |
| `everything` 127 行 | 0.00% [−2.22%, +2.27%] |
| `everything` 128 行 | 0.00% [−2.00%, +2.33%] |
| `everything` 4099 行 | 0.00% [−2.08%, +1.89%] |
| `smooth` 4099 行 | 0.00% [0.00%, +2.63%] |

Rust 层也试过将 `powf` 和 `atan2` 批循环移到 `inline(never)` helper，以及进一步拆出保留标量 libm 的一元分发。只拆二元时，Win64 AVX2 分发器和两个 helper 的静态 reload 总数从 395 降到 237，但 4099 行 `everything` 的耗时区间仍跨零。再拆一元后，Linux 同尺寸 `smooth` 变化 −2.96% [−4.61%, −0.94%]，Wine Win64 同用例区间含零，且 Wine 127 行的 `everything` 中位数上升 2.70% [0.00%, +5.88%]。因此保留 Rust 试验为对照，本轮进入仓库的是诊断和 IR 保存成本模型。

### 在机器层合并退出边恢复

重新解码交付 DLL，首个 CPU lane `atan2` 余数调用结束于 RVA `0x176dcb2`，下一次调用在 `0x176dd01`，两者都跳往同一个 libm thunk。前后两组广播从同样三个 RIP 相对常量地址读取。Tape 列执行器的 `0xe42564` 和 `0xe425c5` 调用同一 thunk，每次返回后重复四条常量池广播和两次栈读取。归档于 `target/debug/build/rssn/e2433476288a0385/out` 的 CPU lane ELF 对象在余数调用后的 `0x2808` 直接储存结果并进入分支。两份产物保留各自的 ABI 与布局差异。基线和变换的因果对照使用从当前 CPU lane 和 Tape 源码抽取的同一份 MIR。

`rssn-edge-restore` 处理寄存器分配完成之后、栈帧下标消除之前的机器指令。对含调用的分支尾部，区分常量池读取、不可变传入栈参数、分配器 spill slot，以及独立向量清零。用 LLVM 物理寄存器单元跟踪尾部指令及另一后继的活跃值，检查恢复产生的值是否只流向被重定向的边。栈读取要有目标后端识别的直接帧下标及匹配的内存操作数，跨越的写入也要满足帧对象别名条件。带副作用、使用别名寄存器、带指令调试引用或可能修改同一栈槽的候选退出。多个入边的恢复指令及物理目的寄存器完全相同才共享一份，在已有屏障后的布局位置插入新机器块，并重算活入与 kill/dead 信息。 [LLVM 物理寄存器活跃集合](https://github.com/rust-lang/llvm-project/blob/0a5579eb649d3634a9c30ac873300cc2944fa070/llvm/include/llvm/CodeGen/LivePhysRegs.h)

设一个入边原在频率为 `f` 的块中恢复一组值，进入共享块的概率为 `p`。以成本 `c` 表示这组恢复，预计节省 `f × ((1−p) × c − p)`，其中最后一项按一次新增跳转计费。16 字节以内的 spill 恢复、常量池读取和参数读取各计 2，更宽的 spill 恢复计 3，清零计 1。分组后还要求至少节省两条静态恢复指令，以抵消共享块跳转。这个模型使用机器块频率和分支概率，最终取舍仍需查实际代码和运行时间。

在 LLVM 23.1.5 中，用 `llc -run-pass` 显式安排该机器 pass。基线与变换都从同一份分配后 MIR 恢复后端。常规 Rust 构建继续使用原 IR 插件入口。CPU lane 在 MSVC AVX-512 对照中识别了 16 个入边，每边 7 条恢复，包含一次传入栈读取。GNU Win64 AVX2 下，CPU lane 和 Tape 列执行器都识别了 15 个入边，每边 4 条恢复。在本次扫描中 Appell 和 zstd 序列解码无候选。Linux 和 Wine 下的 Windows `llc.exe` 都能装载插件并通过机器指令验证。

| 分配后 MIR 对照 | RIP 相对广播 | 静态 spill / reload | 局部帧 |
| --- | ---: | ---: | ---: |
| CPU lane，MSVC AVX-512 | 99 → 29 | 132 / 567 → 132 / 567 | 3544 字节 |
| CPU lane，GNU Win64 AVX2 | 54 → 28 | 41 / 395 → 41 / 369 | 472 字节 |
| Tape，GNU Win64 AVX2 | 41 → 41 | 81 / 457 → 81 / 416 | 904 字节 |
| CPU lane，Linux SysV AVX2 | 54 → 30 | 45 / 435 → 45 / 387 | 408 字节 |
| Tape，Linux SysV AVX2 | 41 → 41 | 112 / 578 → 112 / 537 | 984 字节 |

广播计数包含整函数整数和浮点 RIP 相对广播。沿 GNU Win64 AVX2 中七行余数的实际控制流，`everything` 的 `pow` 和 `atan2` 各从七次恢复序列降到一次。仅计算这两条链，CPU lane 每批省 24 次常量池读取和 24 次 spill reload，Tape 每批省 36 次 spill reload 和 12 条向量清零。这是逐边与指令核算，耗时统计仍覆盖整个执行器。

两种执行器的 AVX2 基线和变换都通过原有的逐位回归检查。在 Wine 下固定单核，CPU lane 做了 15 对交错样本，Tape 做了 13 对。Linux SysV 固定单核，两个执行器均做了 13 对。表中是线程 CPU 的配对中位变化和 95% bootstrap 区间。各用例的 A/B 输出哈希相同。

| 执行器和真实输入 | Win64 在 Wine 下 | Linux SysV |
| --- | ---: | ---: |
| CPU lane，`everything` 127 行 | 0.00% [−2.22%, +2.27%] | −0.74% [−1.06%, +1.11%] |
| CPU lane，`everything` 128 行 | 0.00% [0.00%, +2.27%] | — |
| CPU lane，`everything` 4099 行 | +2.13% [0.00%, +4.26%] | −0.68% [−2.52%, +0.02%] |
| CPU lane，`smooth` 4099 行 | +2.38% [0.00%, +2.70%] | +0.18% [−1.57%, +3.35%] |
| Tape，`everything` 255 行 | 0.00% [−1.18%, +5.00%] | +1.00% [−0.42%, +1.88%] |
| Tape，`everything` 4099 行 | −1.11% [−1.14%, 0.00%] | +0.07% [−2.43%, +1.67%] |
| Tape，`smooth` 4099 行 | 0.00% [−2.17%, 0.00%] | +3.02% [−0.45%, +4.39%] |

这些输入未确认可重复的整段提速。LLVM 的分支位置与栈帧布局继续参与后端选择，省下的局部恢复也可能被新增跳转抵消。该变换继续通过显式入口使用。

```bash
"$llc" -O3 -no-x86-call-frame-opt -stop-after=virtregrewriter \
  -verify-machineinstrs base.ll -o base.mir
RSSN_SPILL_FILTER=exec_block "$llc" -load-pass-plugin="$plugin" \
  -run-pass=rssn-edge-restore -pass-remarks-analysis=rssn-edge-restore \
  -verify-machineinstrs base.mir -o edge.mir
for mode in base edge; do
  "$llc" -O3 -start-after=virtregrewriter -verify-machineinstrs "$mode.mir" -o "$mode.s"
done
```

AVX2 完整运行的 MIR 打印器在一个冷路径的范围附件中输出了不能重解析的内存地址。两组都删除了同一个 `!range <0x…>` 附件后再从 MIR 生成对象。这个附件不在被改写的分发器中。

Linux SysV 的原始调用帧优化会把栈参数写改成 `push`。LLVM 23.1.5 的 MIR YAML 会丢失 X86 的 `HasPushSequences`。CPU lane 的原始 IR 直接生成对象后通过运行检查。同一输入仅做原始 MIR 往返就出现了函数栈调整失衡。Linux 两侧均在生成 MIR 前使用 `-no-x86-call-frame-opt`，并以 PIC 模式生成、重启和链接，随后用原有检查验证基线与变换。 [X86 YAML 映射](https://github.com/rust-lang/llvm-project/blob/0a5579eb649d3634a9c30ac873300cc2944fa070/llvm/lib/Target/X86/X86MachineFunctionInfo.h) [X86 栈帧选择](https://github.com/rust-lang/llvm-project/blob/0a5579eb649d3634a9c30ac873300cc2944fa070/llvm/lib/Target/X86/X86FrameLowering.cpp)

另作了一个 Win64 文本 MIR 对照，把共享组中的三个无条件前驱恢复为原来的块。该产物通过了机器 verifier。运行 `check` 时的访问异常使这组对照失效。仓库中的变换用 LLVM 的块和寄存器接口重定向整组，生成的两组可执行文件均通过上述运行检查。

### 退出边的布局和逐边代价

恢复序列覆盖续点的全部前驱时，直接放到续点块入口并重算物理寄存器活入集合。多个路径携带不同恢复序列时，仍需独立块。若续点的布局前驱以明确的控制屏障或可解析到 `noreturn nounwind` 函数的调用结束，可以在续点前放置恢复块并直接落入续点。调用支持直接符号和 RIP 相对的精确 GOT 重定位，任意内存间接调用退出该规则。分析核对前驱和续点的机器代码节与隐式落入关系。启用 basic block sections 时，独立恢复块保留显式跳转。

逐边计算 `(1−p) × c − p × j`。其中 `p` 为退出概率，`c` 为恢复代价，`j` 为新增跳转成本。落入续点或在续点内部恢复时 `j=0`。独立跳转模式剔除预计净损失的入边，再用机器块频率对保留项加权，并检查静态恢复指令数量。即使一组条件边预计获益，无条件退出边也需通过自己的代价检查。另一条入边携带不同恢复值的情形、续点原有布局落入、分节或退出概率过高的情形都有独立 MIR 回归。

相同的 CPU lane 和 Tape 输入在 Win64 AVX2 与 Linux SysV AVX2 都从可证明的冷退出块之后直接布局共享恢复。汇编显示 LLVM 的后续块布局已经在旧版中折叠了多处跳转。新版与旧版在四份内核上的最终静态跳转和重读计数相同。两份 Tape 对象逐字节一致。每份 CPU lane 的 `.text` 缩小 16 字节。 [LLVM 机器块布局](https://github.com/rust-lang/llvm-project/blob/0a5579eb649d3634a9c30ac873300cc2944fa070/llvm/lib/CodeGen/MachineBlockPlacement.cpp)

基线、旧版和新版各自通过提取内核的逐位回归。每种输入以固定单核交错运行 12 组，三版输出哈希一致。下表中的线程 CPU 时间以旧版为分母，区间为配对 bootstrap 95%。相对于关闭 pass 的基线，本轮这些输入同样没有确认稳定提速。

| `everything` 输入 | Wine Win64 新版相对旧版 | Linux SysV 新版相对旧版 |
| --- | ---: | ---: |
| CPU lane，127 行 | 0.00% [−2.20%, 0.00%] | +0.46% [−0.93%, +0.93%] |
| CPU lane，128 行 | 0.00% [0.00%, +2.25%] | — |
| CPU lane，4099 行 | +1.00% [−0.98%, +2.08%] | +0.43% [−1.01%, +5.75%] |
| Tape，255 行 | +0.65% [0.00%, +1.91%] | −0.99% [−1.52%, +0.18%] |
| Tape，4099 行 | +0.62% [−1.41%, +1.44%] | −1.53% [−2.37%, +0.24%] |

### 只在单一后继使用的恢复

交付 DLL 的 Tape 在 `0xe42564` 和 `0xe425c5` 两次 `atan2` 返回后恢复分发器状态。Linux debug 的 `eval_columns_range` 对象在节内 `0x1d3c` 通过 GOT 调用 `atan2`，返回后直接写回输出并检查循环边界。Win64 release 的展开与物理分配来自另一套构建。本轮在 AVX2 release 重建的 `eval_columns_range` 还找到另一种模式。块 `bb.148` 调用单行 `tape::apply<f64>` 后，无条件恢复多个物理寄存器，再决定继续逐行执行或退出。只在退出用到的配方是 SysV 的三次 spill 重读和两次清零，以及 GNU Win64 的一次 spill 重读和一次清零。机器分支给退出边的概率为 `0x04000000/0x80000000`，约 3.1%。相邻的 MSVC AVX-512 重建使用不同的配方，出口上需要一次重读和两次广播。

继续处理下一行的一侧没有这些寄存器的活入。需要恢复的后继只有一个前驱，所以直接把指令移到后继入口，分配器槽语义和静态指令数不变。这里的收益项只有 `f × (1−p) × c`，其中 `p` 是边频率估计。显式和落入边都可用这一证明。无条件边、目的块本身含帧调整以及另一边仍需要对应物理寄存器的情形不会用单边规则。上游的 PostRA MachineSink 主要处理物理寄存器复制，这些候选是分配器栈读取和目标重物化指令。 [LLVM MachineSink 源码](https://www.llvm.org/doxygen/MachineSink_8cpp_source.html)

同配置的 MSVC AVX2 Appell 重建也命中一次。计算归一化因子的最后一条调用返回后，原先在检查指数域之前把 quiet NaN 常量读入 XMM0。它实际只服务非有限结果的出口，正常路径随即拿 XMM0 计算循环参数。最终汇编把载入留在非有限出口，常量池读取少一条，`rssn_de_appell` 的带标注重读保持 22 条。这个发现归属于重建中归一化因子检查，交付 DLL 四连 `pow` 的栈参数和系数仍是不同的问题。

三个版本使用同一份分配后 MIR。Tape 的新旧版静态重读总数相同，Win64 为 416 条，SysV 为 537 条。两平台的基线、旧版和新版均通过整数幂及余数逐位检查。表中每组在相同逻辑 CPU 上交错运行 12 次，输出哈希全相同。时间取线程 CPU 的配对中位变化，区间为配对 bootstrap 95%。

| Tape 输入 | SysV 新版相对基线 | SysV 新版相对旧版 | Wine Win64 新版相对旧版 |
| --- | ---: | ---: | ---: |
| `everything`，7 行 | −2.85% [−4.40%, −1.01%] | +0.05% [−0.80%, +4.05%] | 0.00% [−4.17%, +4.08%] |
| `everything`，255 行 | −1.11% [−3.74%, +2.17%] | +1.94% [−1.98%, +4.35%] | 0.00% [−4.00%, +4.17%] |
| `everything`，4099 行 | −2.12% [−5.09%, +2.99%] | −0.72% [−5.91%, +3.19%] | 0.00% [−6.25%, 0.00%] |
| `smooth`，4099 行 | −0.71% [−1.91%, +1.43%] | −1.48% [−2.61%, +1.65%] | 0.00% [−8.33%, 0.00%] |

Appell 的四组原有输入逐位一致，Wine 下 16 对样本的线程 CPU 中位变化为 0.00% [−1.39%, +2.56%]。Tape 的新版对旧版各区间均包含零，当前证据只确认恢复从通用路径移到对应出口。

### FreeBSD amd64 独立 LLVM 宿主

在 KVM 的 FreeBSD 15.1 amd64 来宾系统中，用原生 Clang 19.1.7 和 FreeBSD 官方包提供的 LLVM 23.1.5 构建了插件。该 SDK 包以 FreeBSD 14.4 amd64 为基线，测试的包校验和为 `142592d8f9a56ea3896e3a7d334ded341c1f67c24767faa9101169b565ed2558`。构建使用原生 `llvm-config --cxxflags`、`-fPIC` 和 `-fno-rtti`。FreeBSD 的 `opt` 在真实 Appell IR 命中四个单位底数守卫。FreeBSD 的 `llc` 在 Win64 CPU lane MIR 的 15 条入边共享恢复，并通过机器指令验证生成 COFF。 [FreeBSD llvm23 port](https://cgit.freebsd.org/ports/tree/devel/llvm23)

官方 FreeBSD rustc 1.100.0 的 `-vV` 报告 LLVM 23.1.5。实测直接使用 `-Zllvm-plugins` 时驱动缺少 `llvm::InnerAnalysisManagerProxy<…>::Key` 等 C++ 导出。可运行流程由 rustc 生成原生 FreeBSD IR，LLVM 工具进程加载插件并生成对象，最后用原生 Rust 工具链链接。对从当前源码提取的 CPU lane，`llc` 在 14 条入边共享七条恢复。原生基线和变换均通过整数幂及余数用例的逐位检查。

原生受测函数 `exec_block` 的静态重读从 435 条减至 387 条，其中 `vbroadcastsd` 指令的总数从 52 条减至 28 条。汇编标注的 spill 及无条件跳转分别保持 45 和 160 条。基线和变换均以 PIC 生成，MIR 往返时使用相同的 `-no-x86-call-frame-opt`。每种输入在虚拟机的同一逻辑 CPU 交替运行 13 组，逐组输出哈希一致。下表对每组计算相对基线的线程 CPU 时间变化，再报告中位数和配对 bootstrap 95% 区间。来宾计时波动较大，三组区间均覆盖零。

| `everything` 输入 | 变换相对基线 |
| --- | ---: |
| 127 行 | +2.52% [−34.33%, +12.41%] |
| 128 行 | −8.42% [−26.73%, +47.24%] |
| 4099 行 | −10.33% [−15.37%, +5.26%] |

### 用独立 LLVM 工具处理 fat LTO bitcode

用官方 rustc 1.100.0 为 Windows GNU 目标构建提取的 Tape 执行器。`-Clto=fat -Ccodegen-units=1 --emit=llvm-bc,link -Zno-link` 同时输出完成 fat LTO 的模块 bitcode、单一 CGU 对象和 `.rlink`。独立 LLVM 23.1.5 `opt` 在这份 bitcode 中识别到 `eval_columns_range` 的 35 个 Win64 调用点和 4 个实际使用的栈参数。再由 `llc` 直接生成替代对象，原 rustc 的 `-Zlink-only` 使用原有链接参数生成 PE。原对象和替代对象生成的两个程序都通过 Wine 下的逐位回归。

GNU 试验的外部插件运行报告入口。把整份链接后模块序列化为 MIR 再恢复代码生成时，基线就产生空符号重复定义错误。独立发布管线直接让 `llc` 从 IR 生成对象。

MSVC 对照使用 `x86_64-pc-windows-msvc` 的 `no_std` 壳生成 `.rlink`。壳声明真实 RSSN Appell 的同名 C 入口。基线和变换都把同一份提取的 Appell release IR 交给 Windows 版 `opt.exe -passes=rssn-release,verify` 和 `llc.exe`，再以生成的 COFF 对象替换壳的单个 CGU，交给官方 rustc 的 `-Zlink-only`。一次链接生成 DLL、MSVC 导入库和静态库。变换只启用单位底数保护，在四个循环块各守卫三个 `pow`。静态客户端和通过导入库链接的客户端均验证四组输入，结果哈希为 `08344d98c3ccfef2`。

先用 Linux 宿主的目标编译器生成 `.rlink`，再独立用 Wine 中官方 Windows MSVC `rustc.exe` 和 Windows `lld-link.exe` 重做两组。两个客户端都返回相同哈希。小型 Cargo 探针确认 `cargo rustc` 可以通过 `-Zno-link` 成功结束，同时 `compiler-artifact` 的 `filenames` 为空。接入发布规则时应找到本次最终 crate 的链接记录，在发布前核对实际产物。

Wine 上固定单核做 16 对 DLL 样本，每组使用四个真实输入并运行 150 轮。哈希均与静态客户端一致。两份 DLL 都有 35 个静态调用，`.text` 从 2,822 增至 2,982 字节。线程 CPU 时间的配对中位变化为 **−10.26% [−11.38%, −7.89%]**，墙钟为 −2.98% [−27.12%, +20.99%]。区间使用配对 bootstrap。该实验的裁剪对象通过数学导入库调用 CRT，DLL 以无 CRT 入口的方式链接。整套 RSSN 的 MSVC 发布及 CRT 初始化仍需验证。 [rustc 链接实现](https://github.com/rust-lang/rust/blob/0ef1bc386059596dada8fb65b50425aab2ea637a/compiler/rustc_codegen_ssa/src/back/link.rs)

对这两份 Appell bitcode 重新生成 `virtregrewriter` MIR，外壳 `rssn_de_appell` 的 spill 读取、写入分别维持 22 和 10。调用之后同一块内的 spill 读取从 12 降到 8，不可变固定栈槽从 10 降到 8，常量池读取从 14 降到 12。单位底数检查改变了基本块边界，以上为静态布局计数。运行收益按跳过 `pow` 的实际输入单独测量。

GNU Tape 的链接后 IR 中，选中的 `atan2` 和 `pow` 各有 73 个整数比较结果跨调用活跃，其中 70 个在外层循环头 `bb.7` 定义。该块在进入内层指令分发前便计算多条分支的向量尾部条件和行数门槛。用现有运行频率限制和 `undef` 证明试排循环头指令时，该块没有可用候选。仅为估算潜在收益，在报告原型里放宽循环头候选的 `undef` 前提后，模型从该块选出 82 个放置方案，整函数估计的跨调用保存成本从约 146 万降到 81 万。该对照未改写 IR。这些是 SSA 活跃度和块频率的代理量，尚未证明可减少实际 spill。

### Tape 的循环头去重与物理重读

从实际 `Tape::eval_columns_range` 提取同一份 fat LTO、AVX2、Windows GNU bitcode。LLVM 23.1.5 的 `function(early-cse)` 和 `function(gvn)` 均把选定 `pow` 和 `atan2` 调用处的跨调用比较结果从 73 个减到 24 个，SSA 活跃值总数从 193 个减到 100 个。用 `-opt-bisect=0 -opt-bisect-verbose` 找到 `eval_columns_range` 的 pass 编号，再以 `-opt-bisect=<编号>` 将 EarlyCSE 限制到该函数。基线与变换均由 `llc -O3 -verify-machineinstrs` 生成对象，使用原 rustc `.rlink` 重链，Wine 下均通过逐位回归。 [LLVM EarlyCSE 实现](https://github.com/llvm/llvm-project/blob/main/llvm/lib/Transforms/Scalar/EarlyCSE.cpp)

| 变体 | `eval_columns_range` 静态 Spill | 静态 Reload |
| --- | ---: | ---: |
| 基线 | 68 | 383 |
| 仅 `bb.7` 内的完全相同指令去重 | 68 | 383 |
| 仅 `eval_columns_range` 跑 EarlyCSE | 74 | 393 |
| 全模块跑函数级 GVN | 74 | 393 |

同块去重按操作码、类型、标志和操作数的完全相等进行。`bb.7` 去掉 48 条比较、28 条 `and`、28 条 `add`、14 条 `insertelement` 和 14 条 `shufflevector`。去除注释后，函数的汇编与基线相同。EarlyCSE 单函数变体中，8 字节 `movq` Spill 增加 6 条。8 字节 `movq` 和 `cmpq` Reload 分别增加 4 和 7 条，`addq` Reload 减少 2 条。`vmovsd`、`vmovq`、16 和 32 字节的栈访问计数均相同。两类 libm 各 15 个静态调用点，各调用后 12 条指令内的重读计数逐调用点相同。

在 `virtregrewriter` 后运行 `rssn-machine-stack-traffic`，同一函数的 spill 读取为 377 → 387，spill 写入为 58 → 64。调用后同一块内的 spill 读取均为 22 次，宽于 8 字节的 spill 读取均为 15 次。不可变固定栈槽读取为 16 → 14，常量池读取均为 100 次。MIR 以机器内存操作数计数，汇编以最终带注释的指令计数，两个阶段各自支持 A/B。

Wine 下在同一逻辑 CPU 上交错运行两份可执行文件。两个用例各做 12 对，每样本 4,099 行、600 轮。`everything` 的线程 CPU 配对中位变化为 +1.14%，95% bootstrap 区间为 [0.00%, +7.78%]，输出哈希为 `68d9f12cf7c56c66`。`smooth` 的变化为 0.00% [0.00%, +13.39%]，输出哈希为 `5ae2e4eb225ca50b`。线程 CPU 时间的观测粒度约为 10 ms，主机同时运行构建任务。

该案例表明 LLVM 后端能在局部吸收重复 SSA 表达式。同块去重消除了中间表示的冗余，物理代码无需变化。整函数去重改变了比较和整数值的分配，新增重读也出现在调用返回窗口之外。循环头去重保持为对照试验。

## 源码复用双指数求积的对称节点

交付 DLL 的 Appell 区间用同一网格点先执行 `sinh`、两次 `cosh` 和一个 `exp`。返回并消费第一个权重后，RVA `0x2b6ed0` 翻转同一个网格点的符号，再从 `0x2b6ede` 调用另一组双曲函数。exp-sinh 例程同样在 `0x2bd928` 和 `0x2bd9e5` 处理正负网格点。原有的 Win64 MSVC debug 重建把求积驱动、映射闭包及双曲函数包装器分开，驱动把正负点交给同一映射。DLL 的 release 直接展开这些调用。跨调用的保存与栈访问不能在这两个不同的内联形态之间直接相减。

三种 RSSN 双指数映射只在非负的二进制网格上生成节点对，因而在 Rust 驱动中提供成对映射。tanh-sinh 用偶对称权重和同一个端点距离，分别保留 `b - distance` 与 `a + distance` 的浮点运算及端点检查。exp-sinh 共用 `cosh(t)` 和来自 `sinh(t)` 的尺度，仍独立求两个指数。sinh-sinh 对奇对称的节点取负并共用权重。中心点只传给用户回调一次。回调顺序、跳过无效权重的规则、成对累加顺序和求积预算均保持。Rust 的双曲函数由平台 libm 提供。这里的逐位验证覆盖 glibc 和 Wine CRT。 [Rust 浮点双曲函数](https://doc.rust-lang.org/stable/std/primitive.f64.html#method.sinh)

遍历了从 0 到 6.5、步长为 `2^-20` 的 6,815,745 个输入，覆盖 tanh-sinh 最大预算可以到达的网格。Linux 和 Wine 上复核了 `sinh` 的奇对称、`cosh` 的偶对称，以及嵌套双曲函数。另有 89 个求积对照验证节点序列、估计值、误差、求值次数和收敛标志，包含逆序区间、相邻浮点端点、超大或次正规的区间和非有限输入。Linux 另通过 1,024 个随机位模式区间。源码回归检查了中心点、正负次序和零细化预算。

Linux 用独立的 libm 调用计数探针测量四组实际 RSSN 输入每轮的调用，扣掉相同预热。对应的优化后 Win64 MSVC AVX-512 IR 和机器码均通过 verifier。以下负阶及小参数计数包括成对后被 LLVM 外提的热闭包。

| 实际内核 | 动态 `sinh` | 动态 `cosh` | 动态 `exp` | Win64 MSVC AVX-512 静态 Spill / Reload |
| --- | ---: | ---: | ---: | ---: |
| Appell | 14,144 → 7,072 | 28,288 → 14,144 | 14,144 → 7,072 | 21/33 → 21/33 |
| 小参数 Lerch | 4,008 → 2,028 | 5,936 → 3,016 | 7,588 → 7,030 | 61/100 → 41/68 |
| 正阶 Lerch | 4,576 → 2,288 | 4,576 → 2,288 | 8,349 → 8,349 | 17/31 → 17/32 |
| 负阶 Lerch | 1,664 → 832 | 1,664 → 832 | 3,849 → 3,849 | 21/45 → 21/21 |

外提改变了 ABI 压力报告所需覆盖的调用边。配对后的负阶内核存在 Win64 `fastcc` 热闭包调用，`virtregrewriter` 后确认使用 `csr_win64`，传参为 RCX、RDX 和 XMM2。此前 Win64 压力统计仅接受 C 与显式 `win64cc`，会漏过这些 Rust 边。模型现在纳入 `fastcc`，普通 Win64 跨调用报告在负阶样本由 1 个可识别点变为 12 个，其中 11 个超过模型的保存容量。回归同时确认机器调用的寄存器掩码及报告计数。 [LLVM x86 保存掩码](https://llvm.org/doxygen/X86RegisterInfo_8cpp_source.html)

运行使用从 RSSN 源码提取的生产内核、原有的四组输入、相同 rustc 1.100.0 fat LTO 和 AVX2。配对样本固定单核并交替 A/B 顺序，各平台和内核均为 15 对。线程 CPU 时间记录在主计时区间内，输出哈希全部逐位一致。Wine 线程 CPU 时间约以 10 ms 为单位。独立的墙钟实验波动更大。

| 实际内核 | Linux 线程 CPU 变化（95% 区间） | Wine Win64 线程 CPU 变化（95% 区间） |
| --- | ---: | ---: |
| Appell F₁ | −26.32% [−27.09%, −25.79%] | −18.42% [−21.43%, −17.95%] |
| 小参数 Lerch | −31.41% [−32.47%, −29.59%] | −22.67% [−24.00%, −19.28%] |
| 正阶 Lerch | −28.36% [−29.99%, −27.90%] | −19.12% [−20.83%, −18.31%] |
| 负阶 Lerch | −22.16% [−24.94%, −19.95%] | −16.67% [−18.37%, −14.58%] |

另在 Linux 两侧都启用单位底数、尾端截断、浮点类别和退出算式压缩，再做 11 对线程 CPU 对照。Appell 为 **−19.43% [−20.26%, −18.92%]**，正阶为 **−29.13% [−31.16%, −26.21%]**，负阶为 **−18.01% [−20.59%, −11.50%]**。各自输出哈希仍逐位一致。成对后的负阶原分子链被外提，所以要求调用与算式同处求积循环的尾端 pass 在该函数减少了命中。正阶尾端截断和 Appell 的单位底数改写仍命中，组合管线在 Win64 MSVC 和 GNU Win64 上均通过 IR 与机器验证。

从修改前后两个 RSSN 源码目录分别运行 `bench/de/extract.py`，再按前文的独立 `CARGO_TARGET_DIR` 构建。运行四个子命令 `appell`、`small`、`positive`、`negative`，每样本分别运行 120、630、630 和 1150 轮，并添加 `--thread-cpu`。直接使用 `x86_64-pc-windows-gnullvm` 生成的 Win64 AVX2 可执行文件在 Wine 上运行。此测量覆盖提取内核和独立 Rust 回归，全库测试套件未运行。
