#!/usr/bin/env bash
# Explicit distribution build for rssn C/C++ FFI libraries; keeps ordinary Cargo
# dependency builds lean.
set -euo pipefail

rssn_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
profile="${RSSN_FFI_PROFILE:-release}"
target=""
package=false
dist_name=""
dist_dir="$rssn_root/dist"
version=""
cargo_args=()

while (($# > 0)); do
    case "$1" in
        --target)
            target="${2:?--target requires a value}"
            shift 2
            ;;
        --profile)
            profile="${2:?--profile requires a value}"
            shift 2
            ;;
        --package)
            package=true
            shift
            ;;
        --dist-name)
            dist_name="${2:?--dist-name requires a value}"
            package=true
            shift 2
            ;;
        --dist-dir)
            dist_dir="${2:?--dist-dir requires a value}"
            shift 2
            ;;
        --version)
            version="${2:?--version requires a value}"
            shift 2
            ;;
        *)
            cargo_args+=("$1")
            shift
            ;;
    esac
done

cmd=(
    cargo rustc
    --manifest-path "$rssn_root/Cargo.toml"
    --lib
    --profile "$profile"
    --crate-type staticlib,cdylib
)
if [[ -n "$target" ]]; then
    cmd+=(--target "$target")
fi
if ((${#cargo_args[@]} > 0)); then
    cmd+=("${cargo_args[@]}")
fi

"${cmd[@]}"

if [[ "$package" != true ]]; then
    exit 0
fi

if [[ -n "$target" ]]; then
    artifact_dir="$rssn_root/target/$target/$profile"
else
    artifact_dir="$rssn_root/target/$profile"
fi

if [[ -z "$dist_name" ]]; then
    if [[ -z "$version" ]]; then
        version="$(sed -n '/^\[package\]/,/^\[/{s/^version = "\(.*\)"$/\1/p}' "$rssn_root/Cargo.toml" | head -n1)"
    fi
    target_label="${target:-host}"
    dist_name="rssn-${version:-0.0.0}-${target_label}"
fi

stage_dir="$dist_dir/$dist_name"
rm -rf "$stage_dir"
mkdir -p "$stage_dir"

copied=0
for lib in librssn.a librssn.so librssn.dylib rssn.dll rssn.dll.lib rssn.lib; do
    if [[ -f "$artifact_dir/$lib" ]]; then
        cp "$artifact_dir/$lib" "$stage_dir/"
        copied=$((copied + 1))
    fi
done

if ((copied == 0)); then
    echo "error: no FFI library artifacts found in $artifact_dir" >&2
    exit 1
fi

for doc in rssn.h rssn.hpp LICENSE README.md; do
    if [[ -f "$rssn_root/$doc" ]]; then
        cp "$rssn_root/$doc" "$stage_dir/"
    fi
done

(
    cd "$dist_dir"
    tar -czf "${dist_name}.tar.gz" "$dist_name"
    if command -v sha256sum >/dev/null 2>&1; then
        sha256sum "${dist_name}.tar.gz" > "${dist_name}.tar.gz.sha256"
    else
        shasum -a 256 "${dist_name}.tar.gz" > "${dist_name}.tar.gz.sha256"
    fi
)
