#!/usr/bin/env bash
# Explicit distribution build; keeps ordinary Cargo dependency builds lean.
set -euo pipefail
rssn_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
exec cargo rustc --manifest-path "$rssn_root/Cargo.toml" --lib \
    --profile "${RSSN_FFI_PROFILE:-release}" --crate-type staticlib,cdylib "$@"
