#!/bin/bash

# --- Configuration Variables ---

# Rust Toolchains to iterate over
TOOLCHAINS=("nightly" "beta" "stable")

# CPU targets and their corresponding RUSTFLAGS value
# Note: 'x86v4' is a generic feature set for high-end micro-architectures (Zen 3/4, Ice Lake and newer)
#       The other flags target specific hardware vendors/generations.
declare -A CPU_TARGETS=(
    ["x86v3"]="x86-64-v3"
    ["x86v4"]="x86-64-v4"
    ["zen3"]="znver3"
    ["zen4"]="znver4"
    ["skylake"]="skylake"
    ["cascadelake"]="cascadelake"
    ["icelake"]="icelake"
    ["sapphirerapids"]="sapphirerapids"
    ["alderlake"]="alderlake"
    ["raptorlake"]="raptorlake"
    ["meteorlake"]="meteorlake"
    ["skylake-avx512"]="skylake-avx512"
    ["icelake-client"]="icelake-client"
    ["cooperlake"]="cooperlake"
    # Note: 'skylake' is a common micro-architecture base for x86-64-v3 optimization.
)

# The base target triplet (since we are only compiling for the current host's architecture family)
TARGET_TRIPLET="x86_64-unknown-linux-gnu"

# Output directory for all compiled artifacts
OUTPUT_DIR="build_artifacts_x86"

# IMPORTANT: Replace 'my_library_name' with your actual crate name!
CRATE_NAME="rssn"

# --- Setup ---
mkdir -p "$OUTPUT_DIR"
echo "Starting optimized x86_64 build process..."
echo "Artifacts will be collected in: $OUTPUT_DIR"

for TOOLCHAIN in "${TOOLCHAINS[@]}"; do
    # Ensure the base target is available for the current toolchain (usually automatic)
    # rustup target add "$TARGET_TRIPLET" --toolchain "$TOOLCHAIN"
    
    for BUILD_NAME_SUFFIX in "${!CPU_TARGETS[@]}"; do
        
        # Get the actual RUSTFLAGS value from the associative array
        RUST_TARGET_CPU_VALUE="${CPU_TARGETS[$BUILD_NAME_SUFFIX]}"
        
        # Set the full RUSTFLAGS command
        RUSTFLAGS_VALUE="-C target-cpu=${RUST_TARGET_CPU_VALUE}"

        echo "--- Building $TOOLCHAIN / $BUILD_NAME_SUFFIX (Flag: $RUSTFLAGS_VALUE) ---"
        
        # --- A. The Build Command ---
        
        # Run cargo build with the specific toolchain, target, and flags
        RUSTFLAGS="$RUSTFLAGS_VALUE" \
        cargo +$TOOLCHAIN build \
            --target $TARGET_TRIPLET \
            --all-features \
            --release || {
                echo "ERROR: Build failed for $TOOLCHAIN/$BUILD_NAME_SUFFIX. Skipping..."
                continue
            }

        # --- B. Collect and Rename Artifacts ---
        
        # The base path where Cargo places the build artifacts
        BUILD_PATH="target/$TARGET_TRIPLET/release"
        
        # Function to copy and rename an artifact
        copy_artifact() {
            local ext="$1"
            local default_name="lib${CRATE_NAME}.${ext}"
            
            if [[ -f "$BUILD_PATH/$default_name" ]]; then
                # New filename format: lib<CRATE>.<TOOLCHAIN>.<FEATURE_SUFFIX>.<EXT>
                NEW_NAME="lib${CRATE_NAME}.${TOOLCHAIN}.${BUILD_NAME_SUFFIX}.${ext}"
                cp "$BUILD_PATH/$default_name" "$OUTPUT_DIR/$NEW_NAME"
                echo "Copied: $NEW_NAME"
            fi
        }

        # Check and copy the common library types
        copy_artifact "rlib" # Rust Library
        copy_artifact "a"    # Static Library (e.g., on Linux/macOS)
        copy_artifact "so"   # Dynamic Library (e.g., on Linux)
        # If targeting Windows, you'd add:
        # copy_artifact "lib" # Static Library (Windows)
        # copy_artifact "dll" # Dynamic Library (Windows)

    done
done

echo "--- Build process finished. Optimized artifacts are in $OUTPUT_DIR/ ---"
