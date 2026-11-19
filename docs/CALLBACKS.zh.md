# 受检查的宿主回调

[English](CALLBACKS.md) · [数值回调算子](NUMERICAL_CALLBACKS.zh.md)

RSSN 动态 ABI 支持串行同步宿主回调。将 `HostCallback` 注册到 `CallCtx`，再把返回的句柄传给 `Session::call_ho`。这些句柄按所给顺序附在普通参数之后，对应算子签名末尾的函数参数。

```rust
use std::sync::Arc;
use rssn::api::Session;
use rssn::kernel::{CallCtx, OpError, Value};
use rssn::kernel::plugin::NativeCallback;
use rssn::kernel::signature::{Signature, TypeExpr};

let session = Session::new();
let mut context = CallCtx::new();
let square = context.add_callback(Arc::new(NativeCallback::new(
    |args: &[Value]| {
        let x = args.first().and_then(Value::as_f64)
            .ok_or_else(|| OpError::Host("expected a floating-point input".into()))?;
        Ok(Value::Float(x * x))
    },
    Signature::new(vec![TypeExpr::Float], TypeExpr::Float),
)))?;
let result = session.call_ho(
    "num_integrate",
    &[
        Value::Float(0.0), Value::Float(1.0),
        Value::Float(1e-10), Value::Float(1e-10), Value::int(128),
    ],
    &[square], &context,
)?;
# Ok::<(), OpError>(())
```

结果记录包含 `value`、`error`、`evaluations` 和 `converged`。`error` 是自适应 Gauss–Kronrod 的误差估计。子区间预算耗尽时 `converged` 为假。相等的上下限直接得到零积分。交换上下限会改变积分符号。首个宿主错误保留错误码和消息并向上传递，之后求值器不再调用宿主。

回调注册库在注册时保存并检查签名。签名必须具体、无类型约束且有固定参数个数。调用前检查参数，返回数值内核前检查结果。张量、记录和变体递归检查。兼容函数类型要求参数逆变、结果协变，参数个数相同。句柄类型取自注册库。

动态求值开始前会释放会话图的借用，因此回调可以在同一 `Session` 内嵌套调用。宿主提供签名或执行代码时，不持有图借用或回调库锁。使用句柄期间，其 `CallCtx` 必须保留回调对象。`NativeCallback` 要求 `Send + Sync`。清单声明使用回调的算子具有会话亲和性和宿主副作用。除非描述符提供对应实现，否则由解释器执行这些算子。

## 物理量纲

`num_integrate` 的量纲契约如下。

```
(a: Float[D], b: Float[D], abs_tol: Float[D*E], rel_tol: Float,
 max_panels: Int, f: (Float[D]) -> Float[E])
  -> {value: Float[D*E], error: Float[D*E], evaluations: Int, converged: Bool}
```

应给出具体回调量纲。适配器向宿主传入相干量值，并验证返回量纲。上下限须同量纲，相对容差无量纲，非零绝对容差和积分结果同量纲。绝对容差恰为零时可匹配任意量纲。两种容差均非负且不能同时为零。未知量纲不足以建立确定的回调契约。

语言桥接保留浮点表示、闭包生存期、线程约束以及具体类型和量纲依据。每次调用可创建新的作用域内句柄。源语言中的持久回调对象应与仅在单个 RSSN 调用上下文内有效的原始句柄分开。

## 清单兼容性

当前清单 API 版本为 2。算子元数据可含来源 `rust_path`，回调求值器标记为 `HostCallback`。解码前协商 API 版本，二进制记录增加字段会改变布局。读取较旧 JSON 时，可用默认值处理缺少的来源字段。来源信息描述适配器与 Rust 定义的关系。宿主应使用清单的签名、能力和后端字段选择实际可用路径。
