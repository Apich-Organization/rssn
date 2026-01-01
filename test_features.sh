#!/bin/bash

# Configuration: Adjust these variables as needed
CRATE_NAME="your-crate-name" # Replace with the actual name of your crate
CARGO_COMMAND="cargo build --features"
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
NC='\033[0m' # No Color

# --- Define the combinations to test ---
# 1. Individual features (excluding `default` and meta features)
INDIVIDUAL_FEATURES=(
    "compute"
    "output"
    "physics"
    "plugins"
    "jit"
    "ffi_api"
    "ffi_blinding"
)

# 2. Key combination features
COMBINATION_FEATURES=(
    "plugins,ffi_api"    # A logical combination for external/API usage
    "compute,jit"        # Combination of computation and JIT
    "output,ffi_blinding" # Combination of output and a specific FFI feature
)

# 3. Meta features
META_FEATURES=(
    "experimental"
    "nightly"
    "full"
)

# --- Helper function to run and check compilation ---
run_test() {
    FEATURES="$1"
    echo -e "${YELLOW}--- Testing: ${CRATE_NAME} with features: [${FEATURES}] ---${NC}"

    # Use "--no-default-features" only if the feature list doesn't include "default" 
    # and if we are testing a non-full/non-experimental set.
    # We use the full command here for simplicity and explicit feature setting.
    
    if $CARGO_COMMAND "$FEATURES" ; then
        echo -e "${GREEN}✅ SUCCESS: [${FEATURES}] compiled successfully.${NC}"
        return 0
    else
        echo -e "${RED}❌ FAILURE: [${FEATURES}] FAILED to compile!${NC}"
        return 1
    fi
    echo
}

# --- Main Test Execution ---

echo "Starting Feature Gate Compilation Tests for: ${CRATE_NAME}"
echo "--------------------------------------------------------"

# 1. Test the base compilation (default features)
run_test "default"

# 2. Test individual features
echo -e "\n## Testing Individual Features ##"
for feature in "${INDIVIDUAL_FEATURES[@]}"; do
    run_test "$feature"
done

# 3. Test key combinations
echo -e "\n## Testing Feature Combinations ##"
for features in "${COMBINATION_FEATURES[@]}"; do
    run_test "$features"
done

# 4. Test meta features (full, experimental, nightly)
echo -e "\n## Testing Meta Features ##"
for feature in "${META_FEATURES[@]}"; do
    run_test "$feature"
done

# 5. Test with no default features (using a complex set)
# This is a critical test to ensure features work when `default = []` is truly respected.
NO_DEFAULT_SET="output,plugins,compute,ffi_api"
echo -e "\n## Testing No-Default Features: ${NO_DEFAULT_SET} ##"
if cargo build --no-default-features --features "$NO_DEFAULT_SET" ; then
    echo -e "${GREEN}✅ SUCCESS: [--no-default-features + ${NO_DEFAULT_SET}] compiled successfully.${NC}"
else
    echo -e "${RED}❌ FAILURE: [--no-default-features + ${NO_DEFAULT_SET}] FAILED to compile!${NC}"
fi

echo -e "\n--------------------------------------------------------"
echo "Feature gate testing complete."
