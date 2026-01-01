#!/bin/bash

# --- Configuration ---
# 1. Directory to scan for item definitions (your main source code)
SOURCE_DIR="src" 

# 2. The path to your prelude file (the module you want to check against)
#    (e.g., src/lib.rs, src/prelude/mod.rs, or src/utils.rs)
PRELUDE_FILE="src/prelude.rs" 

# 3. Output file
OUTPUT_FILE="public_items_not_in_prelude.txt"

# --- Main Logic ---

# Check if the prelude file exists
if [ ! -f "$PRELUDE_FILE" ]; then
    echo "🛑 Error: Prelude file '$PRELUDE_FILE' not found. Please update PRELUDE_FILE variable."
    exit 1
fi

echo "🔎 Analyzing public items in '$SOURCE_DIR' not found in '$PRELUDE_FILE'..."

# Step 1: Find all public item names in the source directory
# - Uses grep to find lines starting with 'pub ' and containing item keywords.
# - Uses sed to extract ONLY the item name (e.g., 'MyStruct' from 'pub struct MyStruct {')
# - The names are written to a temporary file.
grep -r -h --include=\*.rs "^pub \(struct\|fn\|enum\|trait\|mod\)" "$SOURCE_DIR" | \
sed -E 's/^pub (struct|fn|enum|trait|mod) ([A-Za-z0-9_]+).*/\2/g' | \
sort -u > /tmp/public_items_list.txt

# Step 2: Identify which of these public items are NOT re-exported in the prelude file
# We iterate through the list of public items and check for their presence in the prelude.
MISSING_COUNT=0
echo "" > "$OUTPUT_FILE"

while IFS= read -r ITEM_NAME; do
    # Check if the item name is present in the prelude file.
    # We look for 'pub use .*::ITEM_NAME' or 'use .*::ITEM_NAME'
    # NOTE: This check is an approximation, as a re-export might be aliased (e.g., 'as MyNewName')
    if ! grep -q -E "(pub )?use .*::$ITEM_NAME(|;|,)" "$PRELUDE_FILE"; then
        
        # If the item is not found, find its full path for the output file
        # We search again for the definition line to get the file path.
        FULL_PATH=$(grep -r -l -h --include=\*.rs "^pub .* $ITEM_NAME" "$SOURCE_DIR" | head -n 1)
        
        if [ -n "$FULL_PATH" ]; then
            # Format: file_path::item_name
            echo "$FULL_PATH::$ITEM_NAME" >> "$OUTPUT_FILE"
            MISSING_COUNT=$((MISSING_COUNT + 1))
        fi
    fi
done < /tmp/public_items_list.txt

# 3. Report and Clean up
if [ "$MISSING_COUNT" -gt 0 ]; then
    echo "✅ Scan complete. Found $MISSING_COUNT public items not re-exported in '$PRELUDE_FILE'."
    echo "📝 Full paths written to '$OUTPUT_FILE'"
else
    echo "🎉 Scan complete. All public items found were re-exported in '$PRELUDE_FILE' (based on name match)."
fi

rm /tmp/public_items_list.txt
