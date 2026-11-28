#!/usr/bin/env bash
# Compatibility entry point for checkout-based binding conformance.
set -euo pipefail

root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"
bazel_executable=$(command -v bazel || command -v bazelisk) || {
    printf 'Install Bazel or Bazelisk to run the binding suites.\n' >&2
    exit 2
}
cd "$root"
exec "$bazel_executable" run //:bindings_legacy -- "$@"
