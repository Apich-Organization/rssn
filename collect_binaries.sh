#!/bin/bash
set -e

# --- Configuration ---

# 1. Output directory name
OUTPUT_DIR="dist"

# 2. Project Name / Executable Name
# The script attempts to read the package name from Cargo.toml. If unsuccessful, 
# it uses the current directory name.
if [ -f Cargo.toml ]; then
    EXECUTABLE_NAME=$(grep -E '^\s*name\s*=' Cargo.toml | head -1 | cut -d '=' -f 2 | tr -d ' '\"\' | tr -d '\r')
fi
if [ -z "$EXECUTABLE_NAME" ]; then
    EXECUTABLE_NAME=$(basename "$PWD")
    echo "⚠️ Warning: Could not extract project name from Cargo.toml. Using directory name: $EXECUTABLE_NAME"
fi

# 3. Target List (Must match the list used in the compilation script)
ALL_TARGETS=(
  x86_64-apple-darwin aarch64-apple-darwin
  x86_64-pc-windows-msvc aarch64-pc-windows-msvc i686-pc-windows-msvc
  x86_64-unknown-linux-gnu aarch64-unknown-linux-gnu
  x86_64-unknown-linux-musl armv7-unknown-linux-gnueabihf i686-unknown-linux-gnu
  x86_64-unknown-freebsd x86_64-unknown-netbsd
  aarch64-linux-android x86_64-linux-android
  wasm32-unknown-unknown
)

# --- Initialization ---
echo "--- 📦 Collecting Release Artifacts into $OUTPUT_DIR/ directory ---"
mkdir -p "$OUTPUT_DIR"
echo "Project Name: $EXECUTABLE_NAME"

# --- Loop through targets ---
for TARGET in "${ALL_TARGETS[@]}"; do
    SOURCE_DIR="target/$TARGET/release"
    
    # Determine file extension
    EXT=""
    if [[ "$TARGET" == *windows* ]]; then
        EXT=".exe"
    elif [[ "$TARGET" == *wasm* ]]; then
        EXT=".wasm"
    fi

    SOURCE_FILE="$SOURCE_DIR/$EXECUTABLE_NAME$EXT"
    DEST_FILE="$OUTPUT_DIR/$EXECUTABLE_NAME-$TARGET$EXT"

    if [ -f "$SOURCE_FILE" ]; then
        cp "$SOURCE_FILE" "$DEST_FILE"
        echo "✅ Copied successfully: $TARGET -> $DEST_FILE"
    else
        echo "⚠️ Warning: File not found at $SOURCE_FILE (Not compiled or failed compilation)."
    fi
done

echo -e "\n=== 🎉 All release artifacts collected ($OUTPUT_DIR) ==="