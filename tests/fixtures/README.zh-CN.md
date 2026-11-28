# 测试固件

[English](README.md) | [简体中文](README.zh-CN.md)

供 `io::hdf5`、`io::parquet` 与 `io::netcdf` 测试读取的小型二进制文件，已直接签入仓库以便 `cargo test` 无需额外工具即可运行。

所有固件均由 [`generate.py`](generate.py) 使用 `h5py` 与 `pyarrow` 生成；在仓库根目录运行以下命令可重新生成：

```sh
python tests/fixtures/generate.py
```

请保持文件体积精简（最大文件为 `hdf5/dense_nested.h5`，约 1.2 MB，用于触发 40,000 个链接的二级间接块分形堆）。

## `hdf5/`

| 文件 | 测试内容 |
|---|---|
| `dense_20.h5` | 单分形堆块内的稠密组 |
| `dense_100.h5` | 多个直接堆块 |
| `dense_1000.h5` | 深度为 2 的 v2 B 树 |
| `dense_order.h5` | `track_order=True`：链接信息消息中的创建顺序标志与双字节创建索引 |
| `dense_nested.h5` | 指向 100 个对象的 40,000 个链接（嵌套间接块） |
| `types.h5` | `enum` 数据集与变长字符串数据集（全局堆） |
| `compound.h5` | 复合数据集与引用数据集（读取器须跳过并正常读取旁边的普通数据集） |
| `dense_attrs.h5` | 14 个属性触发的稠密属性布局（分形堆 + v2 B 树） |
| `chunked_filters.h5` | 带 shuffle、deflate 与 fletcher32 过滤器的分块数据集 |
| `v4_chunked.h5` | v4 分块布局：single-chunk、fixed-array、extensible-array 与 v2 B 树 |
| `vlen_numeric.h5` | 变长 `int32` 与 `float64` 数据集（全局堆） |
| `soft_links.h5` | 绝对与相对 HDF5 软链接 |
| `filters_plugins.h5` | 插件过滤器（`LZ4` 32004、`Zstd` 32015、`Blosc` 32001） |

## `parquet/`

| 文件 | 测试内容 |
|---|---|
| `lists_v1.parquet` | v1 数据页下的 `list<double>`（含空列表/null）、`list<string>` 与 `map<string,int64>` |
| `lists_v2.parquet` | v2 数据页下的相同数据 |
| `codec_zstd.parquet`、`codec_lz4.parquet`、`codec_brotli.parquet` | Zstandard、原始 LZ4 与 Brotli 压缩的 2000 行字符串/浮点表 |
| `page_index.parquet` | 带 `ColumnIndex` 与 `OffsetIndex` 的 2000 行表 |
| `decimal.parquet` | `INT32` `DECIMAL(9,2)` 列 |
| `nested_lists.parquet`、`nested_lists_v2.parquet` | 嵌套 `list<list<int32>>`（含空列表、空内层列表与 null 元素） |
| `encodings_delta.parquet` | `DELTA_BYTE_ARRAY` 字符串列 |
| `encodings_bss.parquet` | `BYTE_STREAM_SPLIT` 双精度列 |

## `netcdf/`

| 文件 | 测试内容 |
|---|---|
| `cdf2.nc` | 带 64 位偏移的经典格式文件（版本字节 2） |
| `cdf5.nc` | 64 位数据格式文件（版本字节 5） |
