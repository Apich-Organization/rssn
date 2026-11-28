"""Exercise in-place Python SDK discovery from a stock release archive."""

import math
import os
from pathlib import Path
import sys


def main():
    package, library = (Path(arg).resolve(strict=True) for arg in sys.argv[1:])
    os.environ.pop("RSSN_LIB", None)
    sys.path.insert(0, str(package / "bindings/python"))
    import rssn

    assert Path(rssn.library_path()).resolve(strict=True) == library
    features = rssn.manifest()["features"]
    assert features["jit"] and features["npy"] and "gpu-wgpu" in features["devices"], features

    with rssn.Session() as session:
        automatic = session.current_backend()
        assert automatic in {"interpreter", "tiered"}, automatic
        assert {"auto", "interpreter"} <= set(session.available_backends())
        assert session.set_default_backend("interpreter") == "interpreter"
        assert session.current_backend() == "interpreter"
        assert session.set_default_backend("auto") == automatic
        assert session.current_backend() == automatic

        payload = session.encode_npy((2, 3), [-0.0, 1.25, 2.0, 3.0, 4.0, 5.0])
        assert payload.startswith(b"\x93NUMPY")
        with session.decode_npy(payload) as array:
            shape, samples = array._tensor_parts()
        assert shape == [2, 3] and samples == [0.0, 1.25, 2.0, 3.0, 4.0, 5.0]
        assert math.copysign(1, samples[0]) == -1
        print("Packaged Python discovered the stock library and verified backend selection and NumPy I/O")


if __name__ == "__main__":
    main()
