# 语言绑定依赖的 ABI 契约

[English](CONTRACT.md) | [简体中文](CONTRACT.zh-CN.md)

`python/`、`julia/`、`cpp/` 与 `r/` 均基于 `src/ffi/` 生成的 `rssn.h` C ABI 与内核清单编写。规范性定义见 `src/ffi/mod.rs`、[docs/KERNEL.zh.md](../../docs/KERNEL.zh.md) 与 [docs/DIMENSIONS.zh.md](../../docs/DIMENSIONS.zh.md)。

## 句柄与生命周期

| 规则 | 绑定层行为 |
|---|---|
| `rssn_session_new` / `rssn_config_new` 返回堆对象，由 `rssn_*_free` 释放 | RAII：Python `__del__`/上下文管理器、Julia `finalizer`/`close`、C++ 析构函数、R 外部指针终结器 |
| 项为所属会话内的 `uint32_t` 索引 | 项保持对其会话的引用；跨会话使用会在调用前被拒绝 |
| `RSSN_NO_TERM`（`UINT32_MAX`）表示构造失败 | 转换为携带对应状态提示（如 `Parse`、`UnknownOperator`）的异常 |
| 字符串与字节缓冲区由调用方通过 `rssn_string_free` / `rssn_bytes_free` 释放 | 统一读取辅助函数在拷贝后立即释放原指针 |
| `RssnTensor` 持有 `shape` 与 `data`，由 `rssn_tensor_free` 释放 | 拷贝数据后在 `finally`/`defer` 等价块中释放 |

## 状态码与线程模型

- `RssnStatus` 为统一错误枚举；绑定层在异常中携带状态码并立即读取线程局部 `rssn_last_error()`（除 R 的 `call_value` 对数值字面量做类型重试外，不解析错误文本）。
- `rssn_compute_all_with_metadata` 按请求顺序填充元数据数组；`rssn_compute` 与 `rssn_dimension` 将量纲诊断写入线程局部报告（通过 `rssn_dimension_count` / `rssn_dimension_code` 读取）。
- `Session` 为单线程（内部 `RefCell`）：绑定层为每个线程/任务提供独立默认会话，并通过 `compute_all` 在单会话内并行执行。

## 两套 ABI 与可达性

| ABI | 入口 | 载体 |
|---|---|---|
| term | `rssn_term_apply`、`rssn_compute`、`rssn_compute_all`、渲染器 | 指向图节点的 `uint32_t` 句柄 |
| value | `rssn_value_*` | 会话拥有的 `Value` 细胞（`0` 无效，由 `rssn_value_free` 释放） |

`KernelManifest.backend_support.ffi` 对所有注册算子均为 `true`（其中 189 个高阶算子通过 value ABI 与 `rssn_callback_new*` 路由）。超过有符号 64 位的整数通过 `rssn_value_new_int_decimal` 与 `rssn_value_to_text` 读写。

## 在线算子契约与宿主回调

- `rssn_op_info_json` 返回规范名称、重载签名、量纲元数据与 16 字节稳定签名摘要键；`rssn_value_call_key` 按选定键分派重载（Python/Julia/R 暴露 `op_info` 与 `call_key`，C++ 暴露 `Session::op_info_json` / `Session::op_key` / `call_key`）。
- 宿主回调通过 `rssn_callback_new_for`、`rssn_callback_new` 与 `rssn_callback_signatures` 注册到 `RssnSession::callbacks`；回调参数在调用期间借用，返回值须为通过 `*out` 交出的拥有所有权细胞（可用 `rssn_value_clone` 复制借用参数）。回调在进入计算的同一线程执行。

## 封装层职责边界

封装层不添加任何引擎之外的语义：算子门面、单位/常数与量纲文本均来自清单；仿真场景目录扫描自 `src/sim/scenario*.rs`；JSON 编解码使用各宿主语言标准实现。
