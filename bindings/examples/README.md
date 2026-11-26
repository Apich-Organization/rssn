# Examples: host callbacks

One program per language, all doing the same thing: hand a function of that
language to a higher-order operator (`num_brent`, `num_bfgs`) and read the
answer back.

| File | Run |
|---|---|
| [`python/callbacks.py`](python/callbacks.py) | `RSSN_LIB=… python bindings/examples/python/callbacks.py` |
| [`julia/callbacks.jl`](julia/callbacks.jl) | `RSSN_LIB=… julia --project=bindings/julia bindings/examples/julia/callbacks.jl` |
| [`cpp/callbacks.cpp`](cpp/callbacks.cpp) | `cmake --build build/cpp --target rssn_callbacks && ./build/cpp/rssn_callbacks` |

Build the shared library first:

```sh
RSSN_FFI_PROFILE=dev bash scripts/build_ffi.sh
export RSSN_LIB=$PWD/target-bindings/debug/rssn.dll   # librssn.so / librssn.dylib
```

What each program shows:

1. **`num_brent`** — `sqrt(2)` as the root of `x**2 - 2`, with the residual
   written as a lambda. The engine calls it once per probe.
2. **The generated facade** — `rssn.num_brent(...)` / `num_brent(...)` routes a
   higher-order operator through the value ABI, so the callable argument works
   there too.
3. **`num_bfgs`** — a two-dimensional objective *and* its gradient: two
   callbacks in one call, at argument positions 4 and 5.
4. **An explicit callback** — registered once, with its signature spelled out,
   for when the same function goes to several calls.
5. **A failure** — an exception in the callback aborts the operator and its
   message crosses as the `RssnStatus` error text.

The operator's declaration is the source of truth for each callback's type;
`rssn_callback_signatures` (Python `rssn.callback_signatures`, Julia
`callback_signatures`) prints what an operator accepts at a position.
