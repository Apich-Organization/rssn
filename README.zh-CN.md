# RSSN 科学计算库

[English](README.md)

RSSN 为 Rust 程序提供符号变换、数值内核和物理仿真。符号请求以项表示，在哈希驻留的等价图中执行重写，再按目标提取答案。RSSN 也是 Aevia 的计算内核。

## 体验 3.0 beta API

在现有 Cargo 项目中，从 `v0.3` 分支的本地检出添加 RSSN。如果已检出该分支，使用已有检出的路径即可。

```console
git clone --branch v0.3 https://github.com/Apich-Organization/rssn.git ../rssn
cargo add --path ../rssn
```

将以下程序保存为 `src/main.rs`，然后运行 `cargo run`。

```rust
use rssn::api::*;

fn main() -> Result<(), ComputeError> {
    let session = Session::new();
    let x = session.sym("x");
    let request = diff(sin(x) * cos(x), x);
    let answer = session.compute(request, &Config::new())?;
    println!("{}", answer.term);

    Ok(())
}
```

`Config::new()` 采用符号目标。要在指定值上求数值，可用
`Config::new().numeric(1e-12).bind("x", 0.0)`。物理量纲支持建议检查和严格检查两种模式，
详情见[量纲约定](docs/DIMENSIONS.zh.md)。

| 资料 | 用途 |
|---|---|
| [示例](examples/README.md) | 比较符号结果、数值内核和仿真 |
| [系统架构](https://aevia.apich.org/architecture/03-rssn-engine?lang=zh) | 理解等价图、重写调度和执行后端 |
| [内核嵌入](docs/KERNEL.zh.md) | 查询内核清单、算子、回调和线格式 |
| [算法覆盖](docs/MATH_COVERAGE.zh.md) | 选择算法并了解边界行为 |
| [构建](docs/BUILDING.zh.md) | 构建 Rust 库、C 库和可选后端 |
| [文档索引](docs/README.zh.md) | 查找后端、回调和嵌入契约 |
| [变更日志](CHANGELOG.md) | 了解版本变化 |

## 从源码构建

使用 Rust 1.100 或更新版本。默认构建产生 Rust 库。`jit` 启用 Cranelift，`npy` 启用数组文件读写，
`gpu-wgpu` 启用 GPU 设备接口。要生成 C 库及头文件，请使用[构建指南](docs/BUILDING.zh.md)中的脚本。

```console
cargo build --locked
cargo test --locked --test dimensions
cargo run --locked --example quickstart
```

## 许可

[Apache License 2.0](LICENSE)。
