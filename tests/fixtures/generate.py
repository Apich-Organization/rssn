"""Regenerate the HDF5 and Parquet fixtures under tests/fixtures.

Run from the repository root:

    python tests/fixtures/generate.py

Requires h5py (with the HDF5 library) and pyarrow.
"""

import os

import h5py
import hdf5plugin
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

ROOT = os.path.dirname(os.path.abspath(__file__))
HDF5 = os.path.join(ROOT, "hdf5")
PARQUET = os.path.join(ROOT, "parquet")
NETCDF = os.path.join(ROOT, "netcdf")


def write_hdf5(name, body):
    path = os.path.join(HDF5, name)
    with h5py.File(path, "w", libver="latest") as handle:
        body(handle)
    print(f"{name}: {os.path.getsize(path)} bytes")


def dense(count):
    def body(handle):
        group = handle.create_group("many")
        for index in range(count):
            group.create_dataset(f"d{index}", data=np.arange(2, dtype=np.int32))

    return body


def main():
    os.makedirs(HDF5, exist_ok=True)
    os.makedirs(PARQUET, exist_ok=True)

    # A dense group is used from eight links upwards; `libver="latest"` keeps
    # the writer on the dense (fractal heap + version-2 B-tree) layout.
    write_hdf5("dense_20.h5", dense(20))
    write_hdf5("dense_100.h5", dense(100))
    write_hdf5("dense_1000.h5", dense(1000))

    def ordered(handle):
        group = handle.create_group("many", track_order=True)
        for index in range(300):
            group.create_dataset(f"d{index}", data=np.arange(3, dtype=np.int32))

    write_hdf5("dense_order.h5", ordered)

    def nested(handle):
        group = handle.create_group("many")
        targets = [
            group.create_dataset(f"o{index}", data=np.arange(2, dtype=np.int32))
            for index in range(100)
        ]
        for index in range(40000):
            group[f"l{index}"] = targets[index % len(targets)]

    write_hdf5("dense_nested.h5", nested)

    def types(handle):
        handle.create_dataset("vlen_str", data=np.array(["alpha", "beta", "gamma"], dtype=object))
        enum = h5py.enum_dtype({"RED": 0, "GREEN": 1, "BLUE": 2}, basetype="i4")
        handle.create_dataset("enum", data=np.array([0, 2, 1], dtype=enum))

    write_hdf5("types.h5", types)

    def compound(handle):
        record = np.dtype([("a", np.int32), ("b", np.float64)])
        handle.create_dataset(
            "pair", data=np.array([(1, 1.5), (2, 2.5), (3, 3.5)], dtype=record)
        )
        mixed = np.dtype([("label", "S4"), ("n", np.int32)])
        handle.create_dataset(
            "mixed", data=np.array([(b"aa", 7), (b"bb", 8)], dtype=mixed)
        )
        handle.create_dataset("plain", data=np.array([1, 2, 3], dtype=np.int32))
        handle.create_dataset("ref", data=np.array([handle["plain"].ref]))

    write_hdf5("compound.h5", compound)

    def dense_attributes(handle):
        data = handle.create_dataset("data", data=np.arange(4, dtype=np.float64))
        # More than eight attributes force the dense layout.
        for index in range(14):
            data.attrs[f"a{index:02d}"] = float(index) * 1.5
        handle.attrs["title"] = "dense"

    write_hdf5("dense_attrs.h5", dense_attributes)

    def chunked_filters(handle):
        data = np.arange(120, dtype=np.float64).reshape(10, 12)
        handle.create_dataset(
            "deflate_shuffle_f32",
            data=data,
            chunks=(5, 6),
            compression="gzip",
            compression_opts=6,
            shuffle=True,
            fletcher32=True,
        )

    write_hdf5("chunked_filters.h5", chunked_filters)

    def v4_chunked(handle):
        # Single-chunk, fixed-array, extensible-array, and v2 B-tree chunk indexes
        a = np.arange(24, dtype=np.float64).reshape(4, 6)
        handle.create_dataset("single_chunk", data=a, chunks=(4, 6), compression="gzip")
        handle.create_dataset("fixed_array", data=a, chunks=(2, 3), compression="gzip")
        handle.create_dataset(
            "ext_array",
            data=a,
            chunks=(2, 3),
            maxshape=(None, 6),
            compression="gzip",
        )
        handle.create_dataset(
            "btree_v2",
            data=a,
            chunks=(2, 3),
            maxshape=(None, None),
            compression="gzip",
        )

    write_hdf5("v4_chunked.h5", v4_chunked)

    def vlen_numeric(handle):
        dt_i32 = h5py.vlen_dtype(np.dtype("int32"))
        ds_i = handle.create_dataset("vlen_i32", (3,), dtype=dt_i32)
        ds_i[0] = np.array([1, 2, 3], dtype=np.int32)
        ds_i[1] = np.array([], dtype=np.int32)
        ds_i[2] = np.array([10, 20], dtype=np.int32)

        dt_f64 = h5py.vlen_dtype(np.dtype("float64"))
        ds_f = handle.create_dataset("vlen_f64", (2,), dtype=dt_f64)
        ds_f[0] = np.array([1.5, 2.5], dtype=np.float64)
        ds_f[1] = np.array([3.25], dtype=np.float64)

    write_hdf5("vlen_numeric.h5", vlen_numeric)

    def soft_links(handle):
        grp = handle.create_group("grp")
        grp.create_dataset("target", data=np.array([10.0, 20.0, 30.0], dtype=np.float64))
        handle["alias"] = h5py.SoftLink("/grp/target")
        grp["rel_alias"] = h5py.SoftLink("target")

    write_hdf5("soft_links.h5", soft_links)

    def filters_plugins(handle):
        arr = np.linspace(0.0, 99.0, 100, dtype=np.float64)
        handle.create_dataset("lz4", data=arr, chunks=(50,), **hdf5plugin.LZ4())
        handle.create_dataset("zstd", data=arr, chunks=(50,), **hdf5plugin.Zstd(clevel=3))
        handle.create_dataset(
            "blosc_lz4",
            data=arr,
            chunks=(50,),
            **hdf5plugin.Blosc(cname="lz4", clevel=5, shuffle=hdf5plugin.Blosc.SHUFFLE),
        )
        handle.create_dataset(
            "blosc_zstd",
            data=arr,
            chunks=(50,),
            **hdf5plugin.Blosc(cname="zstd", clevel=3, shuffle=hdf5plugin.Blosc.BITSHUFFLE),
        )

    write_hdf5("filters_plugins.h5", filters_plugins)

    scores = pa.array([[1.5, None], None, [], [2.5, 3.5]], type=pa.list_(pa.float64()))
    tags = pa.array([["a", "b"], ["c"], [], None], type=pa.list_(pa.string()))
    mapping = pa.array(
        [{"k1": 1, "k2": 2}, {"k3": 3}, None, {}], type=pa.map_(pa.string(), pa.int64())
    )
    table = pa.table({"scores": scores, "tags": tags, "m": mapping})
    pq.write_table(table, os.path.join(PARQUET, "lists_v1.parquet"))
    pq.write_table(table, os.path.join(PARQUET, "lists_v2.parquet"), data_page_version="2.0")
    print("lists_v1.parquet, lists_v2.parquet written")

    rows = 2000
    codecs = pa.table(
        {
            "label": pa.array([f"kind-{index % 6}" for index in range(rows)]),
            "value": pa.array([index * 0.5 for index in range(rows)]),
        }
    )
    pq.write_table(codecs, os.path.join(PARQUET, "codec_zstd.parquet"), compression="zstd")
    pq.write_table(codecs, os.path.join(PARQUET, "codec_lz4.parquet"), compression="lz4")
    pq.write_table(codecs, os.path.join(PARQUET, "codec_brotli.parquet"), compression="brotli")
    pq.write_table(
        codecs,
        os.path.join(PARQUET, "page_index.parquet"),
        write_page_index=True,
    )
    from decimal import Decimal
    decimal = pa.table(
        {
            "amount": pa.array(
                [Decimal("123.45"), Decimal("-0.50"), Decimal("0.00"), Decimal("9.99")],
                type=pa.decimal128(9, 2),
            ),
        }
    )
    pq.write_table(decimal, os.path.join(PARQUET, "decimal.parquet"))
    print("codec and decimal fixtures written")

    nested = pa.table(
        {
            "xs": pa.array(
                [[[1, 2], [3]], None, [], [[], [4, None]], [[5]]],
                type=pa.list_(pa.list_(pa.int32())),
            )
        }
    )
    pq.write_table(nested, os.path.join(PARQUET, "nested_lists.parquet"))
    pq.write_table(
        nested, os.path.join(PARQUET, "nested_lists_v2.parquet"), data_page_version="2.0"
    )
    print("nested list fixtures written")

    words = ["alpha", "alphabet", "alphanumeric", "beta", "betamax", "gamma", "gammadion", "delta"]
    encoded = pa.table(
        {
            "s": pa.array(words),
            "v": pa.array([float(index) for index in range(len(words))]),
        }
    )
    pq.write_table(
        encoded,
        os.path.join(PARQUET, "encodings_delta.parquet"),
        use_dictionary=False,
        column_encoding={"s": "DELTA_BYTE_ARRAY"},
    )
    floats = pa.table({"v": pa.array([float(index) * 0.5 for index in range(8)])})
    pq.write_table(
        floats,
        os.path.join(PARQUET, "encodings_bss.parquet"),
        use_dictionary=False,
        use_byte_stream_split=True,
    )
    print("encoding fixtures written")

    os.makedirs(NETCDF, exist_ok=True)
    try:
        from netCDF4 import Dataset
    except ImportError:
        print("netCDF4 is not installed; keeping the committed netcdf fixtures")
    else:
        for name, fmt in [("cdf2.nc", "NETCDF3_64BIT_OFFSET"), ("cdf5.nc", "NETCDF3_64BIT_DATA")]:
            with Dataset(os.path.join(NETCDF, name), "w", format=fmt) as handle:
                handle.createDimension("x", 3)
                variable = handle.createVariable("v", "f8", ("x",))
                variable[:] = [1.0, 2.0, 3.0]
                variable.units = "m"
            print(f"{name} written")

    # Media golden fixtures (PNG, GIF, WAV)
    MEDIA = os.path.join(ROOT, "media")
    os.makedirs(MEDIA, exist_ok=True)
    try:
        from PIL import Image
        import struct
        import wave
        from scipy.io import wavfile

        rgb = np.zeros((4, 6, 3), dtype=np.uint8)
        rgb[:, :, 0] = np.arange(6, dtype=np.uint8) * 40
        rgb[:, :, 1] = np.arange(4, dtype=np.uint8)[:, None] * 60
        rgb[:, :, 2] = 128
        Image.fromarray(rgb, mode="RGB").save(os.path.join(MEDIA, "golden_rgb8.png"))

        rgba = np.dstack([rgb, np.full((4, 6), 200, dtype=np.uint8)])
        Image.fromarray(rgba, mode="RGBA").save(os.path.join(MEDIA, "golden_rgba8.png"))

        gray16 = (np.arange(24, dtype="<u2").reshape(4, 6) * 2000)
        Image.frombytes("I;16", (6, 4), gray16.tobytes()).save(os.path.join(MEDIA, "golden_gray16.png"))

        idx_img = Image.fromarray(rgb, mode="RGB").convert("P")
        idx_img.save(os.path.join(MEDIA, "golden_indexed.png"))

        # Animated GIF
        f1 = Image.new("RGBA", (8, 8), (255, 0, 0, 255))
        f2 = Image.new("RGBA", (8, 8), (0, 255, 0, 255))
        f3 = Image.new("RGBA", (8, 8), (0, 0, 255, 255))
        f1.save(
            os.path.join(MEDIA, "golden_anim.gif"),
            save_all=True,
            append_images=[f2, f3],
            duration=[100, 150, 200],
            loop=2,
            disposal=2,
        )

        # WAV golden fixtures (8-bit, 16-bit stereo, 24-bit, 32-bit float, 64-bit float)
        with wave.open(os.path.join(MEDIA, "golden_pcm8.wav"), "wb") as w:
            w.setnchannels(1)
            w.setsampwidth(1)
            w.setframerate(8000)
            w.writeframes(bytes([128, 192, 64, 128]))

        with wave.open(os.path.join(MEDIA, "golden_pcm16_stereo.wav"), "wb") as w:
            w.setnchannels(2)
            w.setsampwidth(2)
            w.setframerate(16000)
            w.writeframes(struct.pack("<4h", 0, 16384, -16384, 8192))

        with wave.open(os.path.join(MEDIA, "golden_pcm24.wav"), "wb") as w:
            w.setnchannels(1)
            w.setsampwidth(3)
            w.setframerate(22050)
            # 3 samples: 0, +4194304, -4194304
            w.writeframes(b"\x00\x00\x00\x00\x00\x40\x00\x00\xc0")

        wavfile.write(
            os.path.join(MEDIA, "golden_f32.wav"),
            16000,
            np.array([0.0, 0.5, -0.5, 0.25], dtype=np.float32),
        )
        wavfile.write(
            os.path.join(MEDIA, "golden_f64.wav"),
            16000,
            np.array([0.0, 0.5, -0.5, 0.25], dtype=np.float64),
        )
        print("media golden fixtures written")
    except ImportError as exc:
        print(f"skipping media golden fixtures: {exc}")

    # VTK golden fixtures (ASCII & Binary)
    VTK = os.path.join(ROOT, "vtk")
    os.makedirs(VTK, exist_ok=True)
    with open(os.path.join(VTK, "golden_structured_ascii.vtk"), "w", encoding="ascii") as vf:
        vf.write("# vtk DataFile Version 3.0\nGolden ASCII\nASCII\nDATASET STRUCTURED_POINTS\n")
        vf.write("DIMENSIONS 3 2 1\nORIGIN 0 0 0\nSPACING 1 1 1\nPOINT_DATA 6\n")
        vf.write("SCALARS u double 1\nLOOKUP_TABLE default\n1\n2\n3\n4\n5\n6\n")
    with open(os.path.join(VTK, "golden_structured_binary.vtk"), "wb") as vf:
        vf.write(b"# vtk DataFile Version 3.0\nGolden Binary\nBINARY\nDATASET STRUCTURED_POINTS\n")
        vf.write(b"DIMENSIONS 3 2 1\nORIGIN 0 0 0\nSPACING 1 1 1\nPOINT_DATA 6\n")
        vf.write(b"SCALARS u double 1\nLOOKUP_TABLE default\n")
        vf.write(np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0], dtype=">f8").tobytes() + b"\n")
    print("vtk golden fixtures written")


if __name__ == "__main__":
    main()
