import re
import sys
from collections import defaultdict

def auto_rename_colliding_aliases(file_path):
    """
    Reads a file, identifies 'pub use' statements where the final exposed name
    (either the item name or its 'as' alias) leads to a collision.
    
    It then renames the colliding aliases by prepending the highest-level 
    module name (e.g., 'symbolic' or 'numerical') found in the path.
    """
    
    # 1. Read the file content
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"Error: File not found at '{file_path}'")
        sys.exit(1)

    # Regex to parse 'pub use' statements, including existing 'as' aliases.
    # Group 1: The full path including the item (e.g., crate::symbolic::...::euler_lagrange)
    # Group 2: The final 'as' alias name (optional)
    PUB_USE_PATTERN = re.compile(
        r'^\s*pub\s+use\s+(.+?)(?:\s+as\s+([a-zA-Z0-9_]+))?\s*;(?:\s*//.*)?\s*$',
        re.IGNORECASE
    )

    # Map to track the first line index where a name was exposed
    exposed_names = defaultdict(list)
    
    # List to store the new, modified lines
    new_lines = []
    
    # List to store information about lines that need renaming later (the colliders)
    colliding_lines = [] 
    
    print(f"--- Stage 1: Identifying Collisions in '{file_path}' ---")

    # --- STAGE 1: Identify all collisions ---
    for i, line in enumerate(lines):
        match = PUB_USE_PATTERN.match(line)
        
        if match:
            full_path_with_item = match.group(1)
            existing_alias = match.group(2)
            
            # 1a. Determine the name exposed to the module
            if existing_alias:
                exposed_name = existing_alias
            else:
                # If no alias, the exposed name is the final item name (e.g., euler_lagrange)
                exposed_name = full_path_with_item.split('::')[-1]

            if exposed_name in exposed_names:
                # Collision detected!
                colliding_lines.append({
                    'index': i,
                    'line': line,
                    'full_path': full_path_with_item,
                    'exposed_name': exposed_name
                })
            else:
                # First time this name is exposed.
                exposed_names[exposed_name].append(i)
                new_lines.append(line)
        else:
            # Not a 'pub use' line, just add it to the final list for now
            new_lines.append(line)

    # --- STAGE 2: Fix Collisions by Renaming ---
    print(f"--- Stage 2: Fixing {len(colliding_lines)} Collisions ---")
    
    modified_count = 0
    
    for collision in colliding_lines:
        i = collision['index']
        original_line = collision['line']
        full_path = collision['full_path']
        exposed_name = collision['exposed_name']
        
        # Heuristic: Find the highest-level namespace keyword for prefixing.
        # Check for 'symbolic', 'numerical', 'physics', etc.
        # Split the path and look for the second-to-last segment (after 'crate::').
        path_segments = full_path.split('::')
        
        # Find the index of the first segment that is *not* a standard keyword like 'crate' or 'super'
        # In 'crate::symbolic::...', this gets 'symbolic'.
        try:
            # Look for the first part after `crate::` or `super::`
            prefix_index = next(
                j for j, segment in enumerate(path_segments) 
                if segment not in ['crate', 'super', '']
            )
            # Use that segment as the prefix
            prefix = path_segments[prefix_index]
        except StopIteration:
            print(f"Warning: Could not determine prefix for line {i+1}. Skipping rename.")
            new_lines.insert(i, original_line)
            continue
            
        # Create the new, unique alias
        new_alias = f"{prefix}_{exposed_name}"
        
        # Check if the exposed name matches the new alias (e.g., if we were to change 'a_b' to 'a_a_b')
        if exposed_name == new_alias:
             print(f"Warning: Calculated new alias '{new_alias}' is the same as the old alias for line {i+1}. Skipping rename.")
             new_lines.insert(i, original_line)
             continue
             
        # Reconstruct the line with the new alias, preserving original whitespace/comments.
        # Find the original alias part (the " as ...;") and replace it.
        # First, find the path component (everything before the semicolon and comment).
        path_part_match = re.search(r'^\s*pub\s+use\s+(.+?);', original_line, re.IGNORECASE)
        if not path_part_match:
            # This shouldn't happen, but as a fallback
            new_lines.insert(i, original_line)
            continue

        base_import = path_part_match.group(1).strip()
        
        # Remove any existing 'as' clause from the base import for a clean overwrite
        base_import_clean = re.sub(r'\s+as\s+[a-zA-Z0-9_]+$', '', base_import)
        
        # Find the trailing part of the line (semicolon, comment, newline)
        trailing_content = original_line.split(';')[-1] if ';' in original_line else ''
        
        # Construct the final modified line
        modified_line = f"pub use {base_import_clean} as {new_alias};{trailing_content}"
        
        # Ensure it keeps the original indentation/leading space
        leading_space = re.match(r'^(\s*)', original_line).group(1)
        modified_line = leading_space + modified_line.lstrip()

        print(f"Collision at line {i+1}: '{exposed_name}'")
        print(f"  Old: {original_line.strip()}")
        print(f"  New: {modified_line.strip()}")
        
        new_lines.insert(i, modified_line)
        modified_count += 1
        
    # 3. Write the modified content back to the file
    try:
        # We need to re-generate the final lines because insertions were done based on original index.
        # A simpler way is to build the list in the second stage.
        
        # Since the lines were built sequentially, we use the list of lines 
        # from the first stage (`new_lines`), and update the colliding lines 
        # using their original index. This can get complicated.
        
        # Let's rebuild the `lines` array completely, replacing originals with modified ones.
        final_lines = list(lines)
        
        for collision in colliding_lines:
             i = collision['index']
             
             # Re-run the logic to get the final line string (easier than storing the string globally)
             
             # Heuristic: Find the highest-level namespace keyword for prefixing.
             path_segments = collision['full_path'].split('::')
             
             try:
                 prefix_index = next(
                     j for j, segment in enumerate(path_segments) 
                     if segment not in ['crate', 'super', '']
                 )
                 prefix = path_segments[prefix_index]
             except StopIteration:
                 # Should have been handled in the loop, but as a safety break:
                 continue

             new_alias = f"{prefix}_{collision['exposed_name']}"

             path_part_match = re.search(r'^\s*pub\s+use\s+(.+?);', collision['line'], re.IGNORECASE)
             if not path_part_match: continue
             base_import = path_part_match.group(1).strip()
             base_import_clean = re.sub(r'\s+as\s+[a-zA-Z0-9_]+$', '', base_import)
             trailing_content = collision['line'].split(';')[-1] if ';' in collision['line'] else ''
             
             modified_line = f"pub use {base_import_clean} as {new_alias};{trailing_content}"
             leading_space = re.match(r'^(\s*)', collision['line']).group(1)
             final_lines[i] = leading_space + modified_line.lstrip()
        
        with open(file_path, 'w') as f:
            f.writelines(final_lines)
        
        print(f"--- Processing Complete ---")
        print(f"Successfully processed: '{file_path}'")
        print(f"Total lines modified: {modified_count}")
        print(f"File updated in place.")
        
    except Exception as e:
        print(f"An error occurred while writing to the file: {e}")
        sys.exit(1)

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python script_name.py <path/to/prelude.rs>")
        sys.exit(1)
        
    file_path = sys.argv[1]
    auto_rename_colliding_aliases(file_path)
