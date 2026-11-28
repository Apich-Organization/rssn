# RSSN Julia 绑定

[English](README.md) | [简体中文](README.zh-CN.md)

基于 `ccall` 的 rssn C ABI 绑定，为每个注册算子提供对应 Julia 函数，无额外包依赖（内置仿真所需的 JSON 编解码器）。

```julia
using RSSN

session = Session()
x = sym(session, "x")

compute(diff(sin(x) * cos(x), x)).text                 # "cos(2*x)"
cfg = bind!(numeric!(Config(), 1e-12), "x", 0.0)
compute(diff(sin(x) * cos(x), x), cfg).value           # 1.0
```

## 安装

```sh
julia --project=bindings/julia -e 'using Pkg; Pkg.instantiate()'
export RSSN_LIB=/path/to/library-directory
```

未设置 `RSSN_LIB` 时，包会依次查找 `lib/`、分发根目录、Bazel FFI 输出目录（`bazel-bin/ffi_distribution/ffi`、`bazel-bin/ffi_distribution_portable/ffi`、`bazel-bin/ffi_dev/ffi`）、Cargo 输出目录与系统加载路径。可通过 `RSSN.libpath()` 查看当前库路径，或在创建会话前调用 `RSSN.set_library!(path)`。

## API 概览

| 对象 | 职责 |
|---|---|
| `Session()` | 拥有所有项；单线程；`close(session)` 或终结器释放 |
| `Term` | 句柄；`text(t)`、`latex(t)`、`typst(t)`、`pretty(t)`、`value(t)`、`tensor(t)` |
| `Config()` | 目标阶段、绑定、假设、预算与量纲模式 |
| `Answer` | `term`、`text`、`value`、`error`、`reduced`、`iterations`、`dimension`、`diagnostic_codes` 及 `.latex`/`.typst`/`.pretty`/`.tensor` |
| `diff(f, x)` | 生成的算子函数，返回**未求值**项 |
| `OPS` / `NAMES` / `BY_RULE_SET` | 算子元数、绑定名称与按规则集索引 |
| `UNITS` / `CONSTANTS` | 清单单位与常数表（精确有理数以字符串表示） |
| `Value` / `call_value` | 动态值 ABI 与原生对象调用 |
| `run_scenario(name, params)` | 运行仿真（JSON 输入，`Dict` 输出） |
| `manifest_json_string()` / `manifest_hash()` | 内核清单 |

## 常见用法与命名规则

- **符号、数值与张量**：`Term` 重载算术运算符与 46 个 `Base` 数学函数（`sin`、`diff`、`sum` 等）；与 `Base` 冲突的非数学名称加下划线后缀（如 `RSSN.collect_`、`RSSN.pi_`）。`tensor(term)` 将行主序 C 数据重排为 Julia 列主序 `Array`；`encode_npy`/`decode_npy` 读写 NumPy 字节。
- **动态值与回调**：`value(session, big(...))` 支持任意精度 `BigInt`；高阶算子支持 Julia `do` 语法块或显式 `Callback`；`checked=true` 直接以原生输入调用动态算子；`with_coherent_dimension` 附带 `DimExpr` 字典。
- **后端与算子检索**：`available_backends`、`current_backend`、`set_default_backend!`、以及 `op_help`、`op_search`、`op_names`。`default_session()` 按任务隔离，`with_session(session) do ... end` 可切换作用域。

## 测试

```sh
export RSSN_LIB=/path/to/rssn.dll
julia --project=bindings/julia bindings/julia/test/runtests.jl
python bindings/tests/run_readme_examples.py julia
```
