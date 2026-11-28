# RSSN 源码绑定生成

[English](README.md) | [简体中文](README.zh-CN.md)

源码生成器无需链接 RSSN 库即可独立构建，即使已签入的代码片段缺失或无效。在 RSSN 仓库根目录运行：

```sh
bazel test //:api_generator_test //:generated_sources_check_test
bazel run //:refresh_contracts
```

使用 `bazel test //:contracts_test` 结合已链接的算子清单校验生成片段。底层生成器支持 `--bindings-report` 或 `--structured-report` 输出单个片段，并可通过可选末尾参数指定备用 RSSN 源码根目录。Aevia 侧的 `rssn-api-audit` 另行比对在线清单与公开 Rust 清单，生成过程本身不依赖该清单。

- **`contracts.tsv`**：遗留数值量纲策略的已审核源，声明参数/嵌套记录绑定、相等约束及标量或结构化结果公式（可按名称复用结构化策略）。缺少枚举分支对应的策略会导致编译失败。
- **`codec-policy.tsv`**：通过可组合编解码器注册具体内核函数（仅含源码路径、ABI 名称、命名空间与已审核量纲策略）。`--codecs-report`、`--codecs-write` 与 `--codecs-check` 从公开 Rust 函数签名自动推导完整调用、借用参数转换、元数与结果类型。`extent_pow(D,arg,axis,multiplier)` 声明由已校验输入张量维度决定的结果幂次；元组映射为记录（`item0`、`item1`…），可选值映射为 `None`/`Some`，固定宽度点向量使用二阶张量。结构化结果编解码器仍来自 `structured-policy.tsv`。
- **回调、单位与标量检查**：`numerical_callbacks` 中的回调算法使用带显式求值预算的串行适配器，保留宿主错误。固定 SI 定律可使用 `src/dim/units.rs` 中的 `base("Length")` 等七个规范基本量纲名。具体公开 `sim` 标量函数遵循与 `kernels`/`backend` 相同的源码检查（`unsafe`、泛型、条件编译、私有或 `&mut` 函数需单独适配器）。
- **采样泛函与特例**：采样泛函分析数据使用 `{item0: coordinates, item1: ordinates}`，坐标网格须有限且严格递增。`num_stats_percentile` 在调用原生可变切片 API 前显式拷贝输入。`classical::damped_harmonic_oscillator`（临界/过阻尼分支相位表示不一致）、有状态仿真模型、写文件场景、后端内部函数及通用回调不在标量策略覆盖范围内。
