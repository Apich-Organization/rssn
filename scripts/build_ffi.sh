#!/usr/bin/env bash
# Compatibility entry point for the RSSN Bazel FFI targets.
set -euo pipefail

rssn_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
profile="${RSSN_FFI_PROFILE:-ffi}"
target=""
package=false
dist_name=""
dist_dir="$rssn_root/dist"
version=""
llvm_overlay=false
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
        --llvm-overlay)
            llvm_overlay=true
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

bazel_executable=$(command -v bazel || command -v bazelisk) || {
    printf 'Install Bazel or Bazelisk to build FFI libraries.\n' >&2
    exit 2
}
bazel_flags=()
for variable in RSSN_TARGET_ENV RUSTC CARGO RSSN_LLVM_CONFIG RSSN_LLVM_IMPORT_LIBRARY; do
    if [[ -n "${!variable:-}" ]]; then bazel_flags+=("--repo_env=$variable=${!variable}"); fi
done
if [[ -n "$target" ]]; then bazel_flags+=("--define=rust_target=$target"); fi
if [[ "$dist_dir" != /* ]]; then dist_dir="$PWD/$dist_dir"; fi
cd "$rssn_root"

dist_flags=(--dist-dir "$dist_dir")
if [[ -n "$dist_name" ]]; then dist_flags+=(--dist-name "$dist_name"); fi
if [[ -n "$version" ]]; then dist_flags+=(--version "$version"); fi
case "$profile" in
    ffi|release)
        for arg in "${cargo_args[@]}"; do
            case "$arg" in
                --locked|--all-features) ;;
                *) printf 'Unsupported release Cargo argument %s.\n' "$arg" >&2; exit 2 ;;
            esac
        done
        if [[ "$package" == true ]]; then
            archive_target=//:package_ffi
            if [[ "$llvm_overlay" == true ]]; then archive_target=//:package_ffi_llvm; fi
            exec "$bazel_executable" run "${bazel_flags[@]}" "$archive_target" -- "${dist_flags[@]}"
        fi
        exec "$bazel_executable" build "${bazel_flags[@]}" //:ffi_distribution
        ;;
    dev|test)
        if [[ "$llvm_overlay" == true ]]; then
            printf 'The LLVM option requires profile ffi or release.\n' >&2
            exit 2
        fi
        if [[ "$package" == true ]]; then dist_flags+=(--package); fi
        exec "$bazel_executable" run "${bazel_flags[@]}" //:ffi_legacy_dev -- --profile "$profile" "${dist_flags[@]}" "${cargo_args[@]}"
        ;;
    *) printf 'Use profile ffi, release, dev, or test.\n' >&2; exit 2 ;;
esac
