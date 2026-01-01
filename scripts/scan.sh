#!/bin/bash

# --- Configuration ---
SOURCE_DIR="src" 
PRELUDE_FILE="src/prelude.rs" 
OUTPUT_FILE="public_items_not_in_prelude_fast.txt"

# --- Main Logic ---

echo "🔎 Analyzing public items in '$SOURCE_DIR' not found in '$PRELUDE_FILE'..."

if [ ! -f "$PRELUDE_FILE" ]; then
    echo "🛑 Error: Prelude file '$PRELUDE_FILE' not found. Please update PRELUDE_FILE variable."
    exit 1
fi

# Step 1: EFFICIENTLY Create a list of all public items with their full paths.
# This uses find and awk for a single, fast pass over the source files.
# The output format is: item_name::file_path
find "$SOURCE_DIR" -name "*.rs" -print0 | xargs -0 awk '
    # Set the field separator to newline to handle multi-line inputs cleanly
    BEGIN { RS="" }
    
    # Process code blocks starting with "pub"
    /pub (struct|fn|enum|trait|mod) / {
        # Iterate through lines in the block
        for (i=1; i<=NF; i++) {
            # Check for pub item definitions and extract the name
            if ($i ~ /^pub (struct|fn|enum|trait|mod)/) {
                # Extract the public item name (e.g., "MyStruct" from "pub struct MyStruct {")
                match($i, /pub (struct|fn|enum|trait|mod) ([A-Za-z0-9_]+)/, arr);
                ITEM_NAME = arr[2];
                
                # Print the result in the format: ITEM_NAME::FILE_PATH
                # FILENAME is a variable set by awk for the current file being processed.
                print ITEM_NAME "::" FILENAME;
            }
        }
    }
' > /tmp/public_items_with_paths.txt

# Step 2: Read all unique public item names from the generated list
cut -d':' -f1 /tmp/public_items_with_paths.txt | sort -u > /tmp/public_item_names_only.txt


# Step 3: Check which items are NOT in the prelude (FAST lookup)
MISSING_COUNT=0
echo "" > "$OUTPUT_FILE"

while IFS= read -r ITEM_NAME; do
    # Check for item name in the prelude file (checking for pub use or use)
    if ! grep -q -E "(pub )?use .*::$ITEM_NAME(|;|,)" "$PRELUDE_FILE"; then
        
        # Instead of running 'grep -r', we look up the full path instantly 
        # from the file generated in Step 1.
        # Format: ITEM_NAME::FILE_PATH
        FULL_DEFINITION=$(grep "^$ITEM_NAME::" /tmp/public_items_with_paths.txt | head -n 1)
        
        if [ -n "$FULL_DEFINITION" ]; then
            # Reformat to: file_path::item_name
            # Cut the part before the first '::' (the item name)
            FILE_PATH=$(echo "$FULL_DEFINITION" | cut -d':' -f3-)
            
            echo "$FILE_PATH::$ITEM_NAME" >> "$OUTPUT_FILE"
            MISSING_COUNT=$((MISSING_COUNT + 1))
        fi
    fi
done < /tmp/public_item_names_only.txt

# 4. Report and Clean up
if [ "$MISSING_COUNT" -gt 0 ]; then
    echo "✅ Scan complete. Found $MISSING_COUNT public items not re-exported in '$PRELUDE_FILE'."
    echo "📝 Results written to '$OUTPUT_FILE'"
else
    echo "🎉 Scan complete. All public items found were re-exported."
fi

rm /tmp/public_items_with_paths.txt /tmp/public_item_names_only.txt
