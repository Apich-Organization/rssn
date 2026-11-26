#!/usr/bin/env bash
# Run every suite: the unit tests of each binding, the shared conformance table
# in all four languages, and the C++ tour.
#
#   RSSN_LIB=/path/to/rssn.dll bash bindings/tests/run_all.sh
#
# Build the library first with `bash scripts/build_ffi.sh`.
set -euo pipefail

root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$root"

if [[ -z "${RSSN_LIB:-}" ]]; then
    for candidate in target-bindings/debug target/debug target/ffi target/release; do
        for name in rssn.dll librssn.so librssn.dylib; do
            if [[ -f "$candidate/$name" ]]; then
                export RSSN_LIB="$root/$candidate/$name"
                break 2
            fi
        done
    done
fi
if [[ -z "${RSSN_LIB:-}" ]]; then
    echo "error: no rssn library found; set RSSN_LIB" >&2
    exit 2
fi
echo "library: $RSSN_LIB"

echo "== python unit tests"
python -m unittest discover -s bindings/python/tests

echo "== conformance (python)"
python bindings/tests/run_conformance.py

echo "== conformance (julia)"
julia --project=bindings/julia bindings/julia/test/runtests.jl

if command -v Rscript >/dev/null 2>&1; then
    if Rscript -e 'quit(status = if (requireNamespace("rssn", quietly = TRUE)) 0L else 1L)'; then
        echo "== conformance (r)"
        Rscript bindings/tests/run_conformance.R
        echo "== r unit tests"
        Rscript bindings/r/tests/test-rssn.R
        Rscript bindings/r/tests/test-json.R
    else
        echo "== r (skipped: install the package first, see bindings/r/README.md)"
    fi
else
    echo "== r (skipped: Rscript not found)"
fi

echo "== c++ (cmake)"
build="$root/build-bindings/cpp-cmake"
libdir="$(dirname "$RSSN_LIB")"
import_library=""
for candidate in rssn.dll.lib librssn.dll.a rssn.lib librssn.a librssn.so librssn.dylib; do
    if [[ -f "$libdir/$candidate" ]]; then
        import_library="$libdir/$candidate"
        break
    fi
done
if [[ -z "$import_library" ]]; then
    echo "error: no import library next to $RSSN_LIB" >&2
    exit 2
fi
echo "import library: $import_library"
cmake -S bindings/cpp -B "$build" -DRSSN_ROOT="$root" \
      -DRSSN_LIBRARY="$import_library" >/dev/null
cmake --build "$build" >/dev/null
ctest --test-dir "$build" --output-on-failure

echo "== readme examples"
python bindings/tests/run_readme_examples.py

echo "all suites passed"
