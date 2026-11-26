"""Cross-verify RSSN-written HDF5, Parquet, PNG, GIF, and WAV files using Python
reference libraries (h5py, hdf5plugin, pyarrow, Pillow, wave).

Usage:
    uv run --with h5py --with hdf5plugin --with pyarrow --with pillow --with numpy python tests/fixtures/verify_outputs.py <dir>
"""

import os
import sys
import wave
import h5py
import hdf5plugin  # noqa: F401
import numpy as np
import pyarrow.parquet as pq
from PIL import Image


def verify_hdf5(d):
    # 1. Chunked + shuffle + deflate + fletcher32
    p = os.path.join(d, "rssn_chunked.h5")
    with h5py.File(p, "r") as f:
        arr = f["/grid"][:]
        assert arr.shape == (6, 8), f"unexpected shape {arr.shape}"
        expected = np.arange(48, dtype=np.float64).reshape(6, 8)
        np.testing.assert_allclose(arr, expected)
        assert f["/grid"].chunks == (3, 4)

    # 2. Plugin filters (LZ4, Zstd, Blosc)
    p = os.path.join(d, "rssn_plugins.h5")
    with h5py.File(p, "r") as f:
        expected = np.arange(64, dtype=np.float64)
        for name in ["lz4", "zstd", "blosc_lz4"]:
            np.testing.assert_allclose(f[f"/{name}"][:], expected)

    # 3. Dense group (v2 superblock, OHDR v2, FRHP/FHDB, BTHD/BTLF)
    p = os.path.join(d, "rssn_dense.h5")
    with h5py.File(p, "r") as f:
        assert "/grp/d0" in f
        assert "/grp/d19" in f
        np.testing.assert_allclose(f["/grp/d5"][:], [5.0, 10.0])

    # 4. Vlen numeric
    p = os.path.join(d, "rssn_vlen.h5")
    with h5py.File(p, "r") as f:
        ds = f["/vlen_i32"][:]
        assert len(ds) == 3
        np.testing.assert_array_equal(ds[0], [1, 2, 3])
        assert len(ds[1]) == 0
        np.testing.assert_array_equal(ds[2], [10, 20])

    # 5. Soft links
    p = os.path.join(d, "rssn_soft.h5")
    with h5py.File(p, "r") as f:
        link = f.get("/alias", getlink=True)
        assert isinstance(link, h5py.SoftLink), f"expected SoftLink, got {type(link)}"
        assert link.path == "/data"
        np.testing.assert_allclose(f["/alias"][:], [1.0, 2.0, 3.0])

    print("HDF5 cross-verification PASSED")


def verify_parquet(d):
    for codec in ["brotli", "zstd", "lz4", "snappy", "gzip"]:
        p = os.path.join(d, f"rssn_{codec}.parquet")
        tbl = pq.read_table(p)
        assert tbl.num_rows == 100, f"{codec} rows={tbl.num_rows}"
        assert tbl.column("id").to_pylist() == list(range(100))

    p = os.path.join(d, "rssn_indexed.parquet")
    pf = pq.ParquetFile(p)
    tbl = pf.read()
    assert tbl.num_rows == 20
    # Verify ColumnIndex / OffsetIndex are readable by pyarrow
    col_meta = pf.metadata.row_group(0).column(0)
    assert col_meta.has_column_index, "missing column_index"
    assert col_meta.has_offset_index, "missing offset_index"
    print("Parquet cross-verification PASSED")


def verify_media(d):
    # PNG
    png_path = os.path.join(d, "rssn_field.png")
    im = Image.open(png_path)
    assert im.size == (16, 12), f"PNG size {im.size}"
    assert im.mode == "RGB"

    # Animated GIF
    gif_path = os.path.join(d, "rssn_anim.gif")
    gim = Image.open(gif_path)
    assert getattr(gim, "is_animated", False) and gim.n_frames == 4, f"GIF frames={getattr(gim, 'n_frames', 1)}"

    # WAV
    wav_path = os.path.join(d, "rssn_tone.wav")
    with wave.open(wav_path, "rb") as w:
        assert w.getnchannels() == 1
        assert w.getframerate() == 16000
        assert w.getnframes() == 1600

    # VTK Binary
    vtk_path = os.path.join(d, "rssn_bin.vtk")
    if os.path.exists(vtk_path):
        with open(vtk_path, "rb") as vf:
            hdr = [vf.readline().decode("ascii").strip() for _ in range(10)]
            assert hdr[0].startswith("# vtk DataFile Version")
            assert hdr[2] == "BINARY"
            assert hdr[3] == "DATASET STRUCTURED_POINTS"
            raw = vf.read(6 * 8)
            vals = np.frombuffer(raw, dtype=">f8")
            np.testing.assert_allclose(vals, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0])

    print("Media (PNG, GIF, WAV) and VTK cross-verification PASSED")


def main():
    if len(sys.argv) < 2:
        print("Usage: verify_outputs.py <dir>")
        sys.exit(2)
    d = sys.argv[1]
    verify_hdf5(d)
    verify_parquet(d)
    verify_media(d)
    print("ALL PYTHON CROSS-VERIFICATIONS PASSED")


if __name__ == "__main__":
    main()
