# RSSN Python 绑定

[English](README.md) | [简体中文](README.zh-CN.md)

基于 `ctypes` 的 rssn C ABI 绑定，为每个注册算子提供对应 Python 函数。无强制第三方依赖；安装 NumPy 时自动用于张量返回。

```python
import rssn

s = rssn.Session()
x = s.sym("x")

s.compute(rssn.diff(rssn.sin(x) * rssn.cos(x), x)).text      # 'cos(2*x)'
s.compute(rssn.diff(rssn.sin(x) * rssn.cos(x), x),
          rssn.Config().numeric(1e-12).bind("x", 0.0)).value  # 1.0
```

## 安装

```sh
pip install -e bindings/python          # 或：pip install bindings/python
export RSSN_LIB=/path/to/library-directory
```

未设置 `RSSN_LIB` 时依次探测 `rssn/_lib/`、分发根目录、Bazel FFI 输出目录、Cargo 输出目录及系统加载路径。可通过 `rssn.library_path()` 查看已加载路径，或在创建会话前调用 `rssn.load(path)`。

## API 概览

| 对象 | 职责 |
|---|---|
| `rssn.Session()` | 拥有所有项；单线程；支持 `with` 上下文或 `.close()` |
| `rssn.Term` | 句柄；`str()` 为中缀文本，另有 `.latex`/`.typst`/`.pretty`、`.value`、`.tensor()` |
| `rssn.Config()` | 目标阶段、绑定、假设、预算与量纲模式 |
| `rssn.Answer` | `.term`、`.text`、`.value`、`.error`、`.reduced`、`.iterations`、`.dimension`、`.diagnostic_codes` |
| `rssn.diff(f, x)` | 生成的算子函数，返回**未求值**项 |
| `rssn.OPS` / `NAMES` / `BY_RULE_SET` / `BY_CATEGORY` | 算子元数、标识符映射与分类索引 |
| `rssn.units` / `rssn.constants` | 清单单位与常数表（精确有理数） |
| `rssn.Value` | 动态值：`kind()`、`float()`、`text()`、`json()`、`to_term()`、`tensor()`、`len()`、`[i]` |
| `rssn.sim` | `run(name, params)`、`names()`、`params(name)`、`defaults(name)` |
| `rssn.manifest()` / `rssn.manifest_hash()` | 内核清单 |

## 常见用法与注意事项

- **命名约定**：推荐使用 `import rssn`（避免 `from rssn import *` 遮蔽 `abs`、`sum`、`pow` 等内置函数）；Python 关键字加下划线后缀（`and_`、`not_`、`or_`）。算术运算符走引擎 `sub`/`neg`/`divide`。
- **数值、张量与宽整数**：`s.compute_all(...)` 并行计算多条请求；`s.encode_npy`/`s.decode_npy` 读写 `.npy` 字节；超出有符号 64 位的 Python `int` 自动走精确十进制通道。
- **动态值、回调与后端**：原生对象（AD 磁带、求解器句柄）通过 `s.call_value(...)` 调用；高阶算子直接传入任意 Python 可调用对象（或长生命周期 `rssn.Callback`）；`checked=True` 直接按签名检查调用动态算子；`with_coherent_dimension` 附带 `DimExpr` 字典；`s.available_backends()`/`s.current_backend()`/`s.set_default_backend()` 管理后端；`rssn.op_help()`/`rssn.op_search()` 检索算子。

## 测试

```sh
export RSSN_LIB=/path/to/rssn.dll
python -m unittest discover -s bindings/python/tests -v
python bindings/tests/run_conformance.py
python bindings/tests/run_readme_examples.py python
```
