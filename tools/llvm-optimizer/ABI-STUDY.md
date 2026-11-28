# 从 RSSN 汇编到 ABI 与私有栈读取优化

## 已有二进制中的三类栈访问

抽查的是本机留存的 RSSN 探针产物。符号与地址只对应这些构建。没有用这些文件推断各路径占整库运行时间的比例。

| 产物与符号 | 指令证据 | 归类 |
| --- | --- | --- |
| `/tmp/rssn-shift-bench.5YCw7C/shift`，`shift_axes_in_place::<Complex<f64>>`，入口 `0x4a0b70` | 6 次寄存器 push 加 `sub $0xf8,%rsp` 把入口栈指针下移 `0x128`。函数随后读 `0x130(%rsp)` 和 `0x138(%rsp)`，分别对应入口栈指针以上 8 和 16 字节。结合源码的三个切片及方向标志，可将这两处识别为调用者传来的轴数量和标志。调用者在 `0x4a3238` 与 `0x4a323a` 执行两次 push，调用结束在 `0x4a3246` 回收 16 字节，和被调端偏移吻合。 | ABI 栈传参 |
| 同一函数，`0x4a0b9b` 附近 | `lea 0xe0(%rsp),%rdi` 将该地址传给 `check_shape`，返回后读 `0xe0(%rsp)` 并检查返回对象的判别值。 | 显式栈返回对象，具备可能可证明的局部所有权 |
| `/tmp/rssn-release-numerics.fvMKBx/real_probe`，双指数求积闭包，`0x4938a3 → 0x4938af` | `xmm1` 写入 `0x40(%rsp)`，经过 `exp@GLIBC_2.29` 之后立刻从同一地址恢复。其后的 `pow` 附近也有类似的向量寄存器保存。 | 跨外部 libm 调用的寄存器保存 |

保留寄存器的入口 push 和尾部 pop 另行记录。以栈指针为底的访存总数会混合上述类别。特别是把一个外部调用的 ABI 溢出参数计入分配器 spill 会误判收益，优化前应同时查看调用者设置参数的指令和被调函数的入口。

同一 rustc 的 fat-LTO IR 中，313 个 RustFFT 定义有 309 个含 `noundef`、173 个含 `range`、306 个含 `captures(none)`、306 个含 `noalias`、267 个含 `dereferenceable`。比如，`&[u64]` 在实际 Rust ABI 微基准中被展开为有 `captures(none)` 的指针和带 `range(i64 0, 1152921504606846976)` 的长度。Rust 的 `fastcc` 在此次检查的 LLVM 23 x86-64 调用分派中仍走 SysV 的参数寄存器分配，因此压力扫描也应覆盖它。[LLVM 23 x86 调用约定定义](https://github.com/rust-lang/llvm-project/blob/0a5579eb649d3634a9c30ac873300cc2944fa070/llvm/lib/Target/X86/X86CallingConv.td)

## Win64 栈参数与跨调用保存

原始 Windows release DLL 中，CPU lane `run_range_chunk` 区间从 RVA `0x176d4b0` 起，八次 push 和 `sub rsp, 0x238` 将栈指针下移 `0x278`。`0x176d538` 的 `0x2c0(%rsp)` 和 `0x176d540` 的 `0x2b8(%rsp)` 分别读入口 RSP 以上 `0x48` 和 `0x40` 字节，即第九和第八个 Win64 参数位置。现有 Linux debug 目录下的同名 ELF 对象以六次 push 和 `sub rsp, 0x278` 建帧，同样两个寻址位移分别对应入口以上 `0x18` 和 `0x10` 字节。两个产物的 ABI 和内联边界不同，无法从位移核实它们读取了同一源参数。[Microsoft x64 调用约定](https://learn.microsoft.com/en-us/cpp/build/x64-calling-convention?view=msvc-170)

从同一生产路径提取并以 LLVM 23 重建的 CPU lane `exec_block` 有八个标量参数。Win64 `fastcc` 将头四个放入 `RCX/RDX/R8/R9`，其余四个在 MIR 为 `fixedStack` 参数，偏移 `32/40/48/56` 不含返回地址。循环里的 `%fixed-stack.1` 对应 `regs.0`，会反复取用。负阶 Lerch 的 `exp_sinh` 闭包同样走 `fastcc` 和 `csr_win64`，但只通过 `RCX/RDX/XMM2` 传入两个捕获环境指针和一个 `double`。捕获环境本身可能在内存中，调用两侧还各有独立的分配器 spill。常规 `cgscc(argpromotion)` 对这个闭包的签名没有产生改动。

`rssn-spill-pressure-report` 对 Win64 普通固定标量签名分别列出入参槽位数、其中有 IR 使用的参数数、带出参槽位的站点数，以及按块频率乘以槽位数选中的出参调用。槽位计数从第五个参数开始，32 字节 home space 另计。向量、聚合、可变参数及特殊 ABI 属性留给后端确认。槽位诊断可独立于函数级 SSA 活跃分析运行，读写次数、尾调用复用和槽位是否被再次分配仍以 MIR 或机器码为准。

以下每组实验只修改重建 IR 中目标内部函数及全部直接调用点。表中的 spill 和 reload 是带 LLVM 汇编注释的静态指令数，入口保存另计。Win64 MSVC x86-64-v4 的负阶 Lerch 调用方在三个变体均为 `sub rsp, 448`、5 次 push、10 条 spill 和 10 条 reload。其被调闭包的变化如下。

| 负阶 Lerch 闭包约定 | `sub rsp` 字节 | push | spill / reload |
| --- | ---: | ---: | ---: |
| `fastcc` | 208 | 3 | 11 / 11 |
| `preserve_mostcc` | 208 | 9 | 11 / 11 |
| `x86_regcallcc` | 176 | 3 | 9 / 9 |

在 CPU lane 的 Win64 GNU AVX2 重建中，将 `regs` 指针和长度放到第三和第四位后，`exec_block` 的局部分配从 472 字节减到 456 字节，静态 spill/reload 由 41/395 增至 42/420。`x86_regcallcc` 省去调用者的四个出参槽位，局部分配为 464 字节，静态 spill/reload 增至 44/447。参数优先级和可用寄存器组变化会影响整个函数的分配。

在 Wine 中用实际 RSSN 提取内核运行 Win64 GNU 版本，单核轮换三个变体，每组 12 个匹配样本，以线程 CPU 时间计算配对中位变化和 bootstrap 95% 区间。Wine 的线程 CPU 时钟在这次测量中以 10 ms 递增。三组内全部输出哈希逐位一致，CPU lane 的 `check` 用例也全部通过。

| 输入 | 变体 | 相对 `fastcc` 的线程 CPU 变化 |
| --- | --- | ---: |
| 负阶 Lerch 500 轮 | `preserve_mostcc` | 0.0% [−5.4%, +5.6%] |
| 负阶 Lerch 500 轮 | `x86_regcallcc` | 0.0% [−5.7%, +2.8%] |
| CPU lane `everything`，4099 行 × 260 轮 | `regs` 优先 | +1.2% [−2.3%, +2.4%] |
| CPU lane `everything`，4099 行 × 260 轮 | `x86_regcallcc` | +2.4% [−3.4%, +3.6%] |
| CPU lane `smooth`，4099 行 × 1300 轮 | `regs` 优先 | −2.8% [−2.8%, +1.4%] |
| CPU lane `smooth`，4099 行 × 1300 轮 | `x86_regcallcc` | 0.0% [−2.8%, 0.0%] |

目前证据支持保留诊断。仅按参数数量、局部帧大小或调用者保存容量选择调用约定，无法预测这两条路径的收益。

## 闭包环境中的间接标量读取

交付 DLL 的 Appell 例程在中心节点第一次 `pow` 前，于 RVA `0x2b6c49` 从环境读指针，再在 `0x2b6c50` 沿指针读系数。第一次调用后的 `0x2b6c7d` 沿另一个从环境取得的指针读第二个系数。Win64 MSVC 调试重建的负阶 Lerch 独立 integrand 在函数偏移 `0xb0` 调用 `f64::ln`。返回后从 `0x38(%rsp)` 恢复环境地址，从 `0x40(%rsp)` 恢复 `t`，在 `0xce` 和 `0xda` 分别读 `a` 和 `log_gamma` 的地址，随后执行内存操作数形式的乘、减。交付 DLL 和调试重建属于不同函数、不同构建。

同源码的 fat LTO 基线把 integrand 内联进两个环境指针和一个 `double` 的热闭包。IR 在正负节点的两段积分体各读取 `order`、`a`、`log_gamma`、`unit_cutoff`、`z` 和 `unit_polynomial`。每个系数经过一次环境槽读取和一次沿指针读取。Horner 系数数组的迭代读取单独计数。分配后的 MIR 中，`a` 保存在 XMM13，`log_gamma` 从环境加载后写进 `%stack.0`，并在 `log` 返回后折叠到 `VSUBSDZrm`。这里可以区分环境读入的依赖链和寄存器分配后来新增的 spill/reload。

模块入口 `rssn-call-abi-report` 现在能在 Linux SysV 和 Win64 报告这个模式。它识别本地定义中从只读、`captures(none)` 指针参数以常量偏移读取的指针，再沿该指针读取常规标量的链。报告区分字段数、实际 IR load 站点数、直接调用点和传入调用方静态 alloca 地址的调用点。最高站点频率以调用方入口频率为分母，由 BFI 估算。该分析仅用于定位内存访问和调用边，字段属性本身无法证明把标量搬到调用方的安全性，也无法预测改写后的寄存器分配。Win64 MSVC 当前负阶 Lerch 样本报告 6 个字段、12 个标量读取和 7 个栈环境调用点。小正阶 Lerch 的 tanh-sinh 闭包也命中，报告 7 个字段、11 个读取和 6 个栈环境调用点。

负阶 Lerch 的环境可直接用 Rust `move` 按值捕获系数。数组仍通过切片访问。LLVM 23 fat LTO 将新闭包内联，消除其七个直接调用和系数的指针追踪。该内联决策同时复制了部分 libm 路径，Win64 MSVC 下 `rssn_de_negative` 与闭包合计的静态 spill/reload 从 21/21 变为 22/46，`exp` 调用从 7 个增至 13 个。SysV 静态 spill/reload 从 53/62 变为 50/123。静态汇编包含循环之外及不同细化阶段的指令，无法按总数判断单次求积开销。[Rust `move` 捕获](https://doc.rust-lang.org/std/keyword.move.html)

使用实际负阶 Lerch 提取内核，依次输入基准内置的四组参数，固定单核交错 17 对样本。Linux AVX2 线程 CPU 时间变化为 **−6.39% [−11.18%, −4.22%]**，每样本 1200 轮。Wine Win64 GNU AVX2 的 3500 轮复测为 **−2.46% [−5.79%, +0.84%]**。各组输出指纹相同。另在 Linux 和 Wine 对 720 组输入逐位比较，覆盖接近允许的负阶下界、接近零的阶数、接近单位奇点的 `z` 和较大的负 `z`。Linux 与 Win64 MSVC 的提取代码通过 IR 和机器指令验证。对原版和按值版运行报告均保持 bitcode 逐字节一致。

小正阶的 tanh-sinh 闭包也做了独立按值捕获实验。Linux AVX2 每样本 630 轮、15 对的线程 CPU 时间为 −0.69% [−2.60%, +1.78%]，输出哈希一致。该路径保留原捕获方式。

## 调用约定的净收益模型

新的模块 pass 在 fat LTO 完成符号封闭后识别本地、直接调用、固定原型的候选。它可以同时把函数和全部调用点改为 `x86_regcallcc`，类型和参数属性原样保留，包括切片的范围、指针属性、以及符合模型的首参数 `sret`。它只处理确定的标量参数。外部函数、地址被取用的函数、可能尾调用或展开的边、安全性和寄存器预算不明的原型均退出。可疑的非叶函数和 SIMD 压力较高的叶函数也退出，避免引入额外的 XMM 保存。[System V AMD64 psABI](https://gitlab.com/x86-psABIs/x86-64-ABI/-/blob/master/x86-64-ABI/low-level-sys-info.tex)

SysV 的 6 个 GP 参数寄存器与本版 LLVM RegCall 的 7 个易失和 4 个非易失 GP 参数寄存器，形成两种边界资源配置。具体排列也影响调用者的复制。试验中，使用 `preserve_nonecc` 的八参数用例在 SysV 调用者中新增了 5 组保存/恢复，即使 SSA 值本身跨调用很少。SysV 调用者还要替它自己的调用者保存非易失寄存器。选用的 RegCall 仍保留这些寄存器，使用 `R12–R15` 传参时才对由此占用的寄存器另计代价。外部库函数的 SysV 接口保持原调用约定。[LLVM 23 callee-saved 定义](https://github.com/rust-lang/llvm-project/blob/0a5579eb649d3634a9c30ac873300cc2944fa070/llvm/lib/Target/X86/X86CallingConv.td)

对每条边，根据类型和位置计算：

- `U` 是被移出栈且有 SSA 使用的 GP 参数数。另算 `R`，即有普通使用所在块支配所有可达返回块的参数数。一个仅在冷分支上使用的参数先只计调用者的写，支配所有返回的使用才给被调端读的收益。
- `A = max(0, GP 参数数 − 7)` 是被参数占去的非易失寄存器数。调用者的跨调用 GP 活跃度为 `L`，潜在额外 spill 数估为 `clamp(L − (5 − A), 0, A)`。
- 每条边的静态分值是消除的栈写和读扣除潜在 spill/reload、以及直接转发同槽参数损失的零复制边。栈访存各按 2 分、溢出往返按 4 分、转发复制按 1 分。把分值按 `BFI(call block) / BFI(entry)` 加总到调用者，再按入口扣一次 `2 × A` 的保存/恢复代价。每条边和每个调用者分别检查，避免热调用掩盖另一调用者的预计回退。

这些分值是机器码资源和活跃范围的估计，对最终延迟仍需 A/B。已有类型范围直接保留在参数上，指针保持原表示。RustFFT IR 中通过签名过滤的 28 个函数均没有通过全量直接调用点检查，ABI report 对这批候选全部放弃改写。

## 用 Rust 捕获属性证明调用不可见的私有栈对象

对静态入口 alloca，从根追踪由 GEP 得到的指针。允许直接的 load/store 地址用途和生命周期标记。对传到调用点的指针，要求 LLVM 的 `captures(none)`，并记录具体哪些调用收到了该指针。指针被存为值、返回、转成地址整数、合流到未解析的 PHI/select、通过 `llvm.localescape` 隐式逃逸，或传给不具备完全捕获约束的调用时，整个对象退出这条规则。

一个普通 opaque 调用没有收到该对象或其派生指针，又没有泄露途径时，才可用这条证明替代原先对该调用的 `nosync`/`nofree`/ModRef 条件。收到对象的 noncapturing 调用仍可能写对象。生命周期、原子、volatile、inline asm、stackrestore、coroutine 和特殊控制流屏障继续检查。在正例中，Rust 常见的 `captures(none)` 初始化调用显式收到栈对象。之后另一个无法取得该对象的 opaque 调用成为可越过的边界，SSA 浮点值的跨调用 spill/reload 被原有对象的单次延后读取替代。[LLVM 长度、别名、捕获和指针语义](https://llvm.org/docs/LangRef.html)

## 复现

```bash
plugin="$PWD/tools/llvm-optimizer/build/RssnLLVMOptimizer.so"

# 独立查看已优化 IR 中的 ABI 候选
opt-23 -load-pass-plugin="$plugin" -passes=rssn-call-abi-report \
  -pass-remarks-analysis=rssn-call-abi -disable-output input.ll

# Rust 三切片、标志与步进参数的真实 fat-LTO ABI A/B
python3 tools/llvm-optimizer/bench/run.py "$plugin" --suite abi
```

rustc 自动流程用 `RSSN_SPILL_ABI=report` 或 `RSSN_SPILL_ABI=transform`，默认关闭。加上 `-Cremark=rssn-call-abi` 和 `-Cllvm-args=-pass-remarks-analysis=rssn-call-abi` 可查看诊断。ABI 模式独立于函数压力 pass，故 A/B 可设 `RSSN_SPILL_MODE=off`。完整结果见 [BENCHMARKS.md](BENCHMARKS.md)。
