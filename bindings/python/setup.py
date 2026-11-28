"""Give wheels containing a bundled RSSN library an explicit platform tag."""
import os
import re
from pathlib import Path
from tempfile import TemporaryDirectory

from setuptools import setup
from setuptools.command.bdist_wheel import bdist_wheel
from setuptools.errors import OptionError


class RSSNWheel(bdist_wheel):
    def finalize_options(self):
        directory = Path(__file__).parent / "rssn" / "_lib"
        self.bundled_rssn = any(
            path.is_file() and (
                path.suffix in {".dll", ".so", ".dylib"} or ".so." in path.name
            )
            for path in directory.glob("*")
        )
        if self.bundled_rssn:
            name = os.environ.get("RSSN_WHEEL_PLATFORM", self.plat_name or "")
            platform = name.lower().replace("-", "_").replace(".", "_")
            if platform == "any" or not re.fullmatch(r"[a-z0-9_]+", platform):
                raise OptionError(
                    "set RSSN_WHEEL_PLATFORM to the platform of the bundled RSSN library"
                )
            self.plat_name = platform
        super().finalize_options()
        if self.bundled_rssn:
            self.root_is_pure = False

    def run(self):
        # A reused build/lib can still contain a native library removed
        # from the sources, which would make the next portable wheel invalid.
        with TemporaryDirectory(prefix="rssn-wheel-") as directory:
            workspace = Path(directory)
            self.bdist_dir = str(workspace / "wheel")
            build = self.get_finalized_command("build")
            build.build_lib = str(workspace / "lib")
            self.reinitialize_command("build_py")
            super().run()

    def get_tag(self):
        if not self.bundled_rssn:
            return super().get_tag()
        # RSSN uses ctypes with a C library and needs no Python extension ABI.
        return "py3", "none", self.plat_name


setup(cmdclass={"bdist_wheel": RSSNWheel})
