import re
import sys
from collections import OrderedDict

def remove_duplicate_pub_use(file_path):
    """
    Reads a file, removes duplicate 'pub use ...;' statements while
    preserving order, and writes the changes back to the file.
    """
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"Error: File not found at '{file_path}'")
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred while reading the file: {e}")
        sys.exit(1)

    # 1. Identify and group all 'pub use ...;' lines
    # The regex captures the entire statement, including leading/trailing whitespace,
    # and the important part (the module path) in a group.
    PUB_USE_PATTERN = re.compile(r'^\s*pub\s+use\s+(.+?);(?:\s*//.*)?\s*$', re.IGNORECASE)

    # OrderedDict to store unique pub use statements and their original line numbers
    # We use OrderedDict to maintain the original order of appearance.
    unique_pub_uses = OrderedDict()
    
    # List to hold lines that are *not* 'pub use' statements
    non_pub_use_lines = []
    
    pub_use_indices = []

    for i, line in enumerate(lines):
        match = PUB_USE_PATTERN.match(line)
        if match:
            # The key is the standardized import path (the content between 'pub use' and ';')
            # Stripping whitespace makes "  mod::item  " and "mod::item" the same.
            key = match.group(1).strip()
            
            if key not in unique_pub_uses:
                # Store the original line content and its index for the first occurrence
                unique_pub_uses[key] = line
                pub_use_indices.append(i)
        else:
            # Keep track of all other lines in their original order
            non_pub_use_lines.append(line)

    # 2. Reconstruct the file content
    
    # Get the unique 'pub use' lines in their order of first appearance
    unique_pub_use_lines = list(unique_pub_uses.values())
    
    # 3. Write the cleaned content back to the file
    
    # Simple strategy: Write back all the original lines, but skip any line 
    # that was identified as a duplicate 'pub use'.
    
    # We need to know which lines to keep. This set contains the indices of 
    # the *first* occurrence of each unique 'pub use' statement.
    indices_to_keep = set(pub_use_indices)
    
    new_lines = []
    
    # Re-iterate over the original lines to decide what to keep
    for i, line in enumerate(lines):
        if PUB_USE_PATTERN.match(line):
            # It's a 'pub use' line. Check if it's the first (and therefore unique) occurrence.
            if i in indices_to_keep:
                new_lines.append(line)
            # Else: it's a duplicate, so we skip it (don't append it to new_lines)
        else:
            # Not a 'pub use' line, so we keep it.
            new_lines.append(line)
            
    # Write the modified content back
    try:
        with open(file_path, 'w') as f:
            f.writelines(new_lines)
        
        original_count = len(lines)
        new_count = len(new_lines)
        removed_count = original_count - new_count
        
        print(f"Successfully processed: '{file_path}'")
        print(f"Total lines removed: {removed_count}")
        print(f"File updated in place.")
        
    except Exception as e:
        print(f"An error occurred while writing to the file: {e}")
        sys.exit(1)

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python script_name.py <path/to/prelude.rs>")
        sys.exit(1)
        
    file_path = sys.argv[1]
    remove_duplicate_pub_use(file_path)
