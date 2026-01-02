import os
import re


def promote_private_to_pub_crate(src_directory="src"):
    """
    Walks through all .rs files (excluding vendor/deps) and promotes private 'fn'
    to 'pub(crate) fn', but ONLY if the function starts at the beginning of a line (no indentation).
    This safely avoids indented functions inside 'impl' or other blocks.
    """

    # Directories to explicitly exclude from processing
    EXCLUDE_DIRS = ["vendor", "deps", "target"]

    # 1. Define the Safest Regular Expression

    # This pattern enforces:
    # 1. ^: Must be the absolute start of the line (no tabs/spaces before it).
    # 2. Group 1: Optional linkage (e.g., 'extern "C" ')
    # 3. Negative Lookahead: Must NOT be followed by 'pub' (to exclude existing pub functions).
    # 4. Group 2: The 'fn' keyword and the rest of the signature.

    VISIBILITY_REGEX = re.compile(
        # Must start the line (Multi-line flag makes ^ match start of line)
        r"^"
        # Group 1: Optional linkage (must immediately follow ^ if present)
        r'((?:extern\s+"[^"]+"\s*)?)'
        # Negative Lookahead: Ensure 'pub' is NOT present
        r"(?!pub\b)"
        # Group 2: The function definition body
        r"(fn\s+.*)",
        re.MULTILINE | re.DOTALL,  # MULTILINE is critical for '^' to work line-by-line
    )

    # 2. Replacement String
    # Group 1 (Linkage) + 'pub(crate) ' + Group 2 (fn ...)
    REPLACEMENT_STRING = (
        r"pub(crate) \1\2"  # Note: Group 1 is now linkage, Group 2 is fn...
    )

    # Let's adjust the grouping and replacement string for clarity:
    # Group 1: Linkage
    # Group 2: fn...
    # REPLACEMENT_STRING = r'pub(crate) \1\2'
    # Let's re-test the grouping on the new pattern:

    # Pattern: ^((?:extern\s+"[^"]+"\s*)?)(?!pub\b)(fn\s+.*)
    # Group 1: Linkage
    # Group 2: fn...

    # Replacement should be: pub(crate) + Linkage + fn...
    REPLACEMENT_STRING = r"pub(crate) \1\2"

    print(f"Starting highly-restricted file walk in directory: {src_directory}")
    print(
        "Promoting only functions that start at the beginning of a line (no indentation)."
    )

    # 3. Walk the Directory with Exclusion Logic
    for root, dirs, files in os.walk(src_directory):
        # Exclude directories
        dirs[:] = [d for d in dirs if d not in EXCLUDE_DIRS]

        for file_name in files:
            if file_name.endswith(".rs"):
                file_path = os.path.join(root, file_name)

                try:
                    with open(file_path, "r", encoding="utf-8") as f:
                        original_content = f.read()

                    new_content, count = VISIBILITY_REGEX.subn(
                        REPLACEMENT_STRING, original_content
                    )

                    if count > 0:
                        print(
                            f"  ✅ Promoted {count} private function(s) to pub(crate) in: {file_path}"
                        )

                        # Write the modified content back to the file
                        with open(file_path, "w", encoding="utf-8") as f:
                            f.write(new_content)
                    # else:
                    #     print(f"  - No changes needed in: {file_path}")

                except Exception as e:
                    print(f"  ❌ Error processing file {file_path}: {e}")

    print("\nModification complete.")


if __name__ == "__main__":
    # You must run 'git checkout .' or similar command to revert all previous changes
    # before running this new, safer script.
    promote_private_to_pub_crate()
