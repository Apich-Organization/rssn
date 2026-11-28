# RSSN R 绑定

[English](README.md) | [简体中文](README.zh-CN.md)

基于 rssn C ABI 的 R 扩展包，无第三方 R 包依赖（内置 JSON 编解码器，运行时动态加载共享库）。

```r
library(rssn)

s <- rssn_session()
x <- sym("x", s)

compute(rssn_diff(sin(x) * cos(x), x))$text                    # "cos(2*x)"
compute(rssn_diff(sin(x) * cos(x), x),
        rssn_config() |> cfg_numeric(1e-12) |> cfg_bind("x", 0))$value   # 1
```

## 安装

```sh
R CMD INSTALL bindings/r
export RSSN_LIB=/path/to/extracted-rssn
```

库探测顺序：`options(rssn.library =)` → `RSSN_LIB` → 已安装包与 SDK 根目录 → Bazel FFI 输出目录 → Cargo 输出目录 → 系统加载路径。可通过 `rssn_library()` 查看路径，或在创建会话前调用 `rssn_load(path)`。

## API 概览

| 对象 | 职责 |
|---|---|
| `rssn_session()` | 拥有所有项；单线程；`close()` 或 GC 释放 |
| `sym()`、`parse_term()`、`term()`、`rational()` | 构建符号、解析文本、强制转换与精确分数 |
| `rssn_term` | 句柄；`term_text()`、`term_latex()`、`term_typst()`、`term_pretty()`、`term_value()`、`to_float()` |
| `rssn_config()` + `cfg_*()` | 目标阶段、绑定、假设、预算与量纲模式 |
| `rssn_answer` | `$term`、`$text`、`$value`、`$error`、`$reduced`、`$iterations`、`$dimension`、`$diagnostic_codes`、`tensor(answer)` |
| `rssn_diff(f, x)` | 生成的算子函数（带 `rssn_` 前缀，不遮蔽基础 R 函数） |
| `OPS`、`OPERATORS`、`op_help()`、`op_table()` | 算子元数与可检索文档表 |
| `call_value()`、`rssn_value`、`rssn_callback()` | 动态值 ABI 与 R 函数宿主回调 |
| `sim_run()`、`SCENARIOS`、`UNITS`、`CONSTANTS` | 仿真调度、单位与常数表 |
| `manifest()`、`manifest_json()`、`manifest_hash()` | 内核清单 |

## 常见用法与注意事项

- **张量与 NumPy**：引擎行主序与 R 列主序在进出时自动转置；一维单元素向量请使用 `array(x, dim = 1L)` 保留向量形状；`rssn_encode_npy`/`rssn_decode_npy` 读写 `.npy` 字节。
- **整数与回调**：R 数值字面量默认为双精度，`call_value()` 会按引擎提示自动翻转整数/浮点类型；超精度整数使用 `value_integer_decimal()`。189 个高阶算子可直接传入 R `function`（或长生命周期 `rssn_callback()`）。
- **错误处理与后端**：引擎错误抛出携带 `$status` 与 `$status_name` 的 `rssn_error` 条件；`rssn_available_backends()`/`rssn_current_backend()`/`rssn_set_default_backend()` 管理后端；`checked = TRUE` 直接按声明签名调用动态算子。

## 测试

```sh
export RSSN_LIB=/path/to/rssn.dll
Rscript bindings/tests/run_conformance.R
python bindings/tests/run_readme_examples.py r
Rscript bindings/r/tests/test-rssn.R
```
