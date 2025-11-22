#! /bin/bash

# --- Configuration ---
# List of all targets to be compiled
ALL_TARGETS=(
  x86_64-apple-darwin aarch64-apple-darwin
  x86_64-pc-windows-msvc aarch64-pc-windows-msvc i686-pc-windows-msvc
  x86_64-unknown-linux-gnu aarch64-unknown-linux-gnu
  x86_64-unknown-linux-musl armv7-unknown-linux-gnueabihf i686-unknown-linux-gnu
  x86_64-unknown-freebsd x86_64-unknown-netbsd
  aarch64-linux-android x86_64-linux-android
  wasm32-unknown-unknown
)

# List of targets that cannot execute 'cargo test' or 'cargo bench'
NO_TEST_OR_BENCH=(
  wasm32-unknown-unknown
  aarch64-linux-android
  x86_64-linux-android
)

# --- User Input ---
echo "--- Rust Cross-Platform Operations Script ---"
echo "Please select the operation type:"
echo "1) Debug Build (cargo build)"
echo "2) Release Build (cargo build --release)"
echo "3) Test (cargo test)"
echo "4) Bench (cargo bench)"
read -p "Enter option (1-4): " CHOICE

case "$CHOICE" in
    1) CARGO_COMMAND="build --features full"; MODE_DESCRIPTION="Debug Build"; ;;
    2) CARGO_COMMAND="build --release --features full"; MODE_DESCRIPTION="Release Build"; ;;
    3) CARGO_COMMAND="test --features full"; MODE_DESCRIPTION="Test"; ;;
    4) CARGO_COMMAND="bench --features full"; MODE_DESCRIPTION="Bench"; ;;
    *) echo "Invalid option. Exiting."; exit 1 ;;
esac

echo -e "\nYou selected: ${MODE_DESCRIPTION}\n"

# --- Execution Loop ---
for TARGET in "${ALL_TARGETS[@]}"; do
    if [[ "$CARGO_COMMAND" == "test" || "$CARGO_COMMAND" == "bench" ]]; then
        # Check if the target should be skipped
        SKIP=false
        for NO_TARGET in "${NO_TEST_OR_BENCH[@]}"; do
            if [[ "$TARGET" == "$NO_TARGET" ]]; then
                echo "--- ⚠️ Skipping $TARGET: '$CARGO_COMMAND' is not supported for this target. ---"
                SKIP=true
                break
            fi
        done
        if $SKIP; then continue; fi
    fi

    # Execute Cargo Command
    echo "================================================================"
    echo "🚀 Running 'cargo $CARGO_COMMAND' for Target: $TARGET"
    echo "================================================================"

    # Run the command
    cargo $CARGO_COMMAND --target "$TARGET"
    
    # Check the exit status of the last command
    if [ $? -ne 0 ]; then
        echo "!!! ❌ $TARGET Build/Test failed. Continuing to the next target..."
    else
        echo "✅ $TARGET execution successful."
    fi
done

echo -e "\n=== All targets operations completed ($MODE_DESCRIPTION) ==="
