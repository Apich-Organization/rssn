import re
import sys
from collections import defaultdict

def auto_rename_colliding_pub_use(file_path):
    """
    Reads a file, identifies 'pub use' statements that lead to name collisions,
    and automatically renames the imported item using the module name as a prefix 
    (e.g., `as physics_sm_solve_advection_diffusion_1d`).
    
    The first line to introduce a name is kept as-is; subsequent colliders are renamed.
    """
    
    # 1. Read the file content
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"Error: File not found at '{file_path}'")
        sys.exit(1)

    # Regex to parse 'pub use' statements
    # Group 1: The full path (e.g., crate::physics::physics_sm::solve_advection_diffusion_1d)
    # Group 2: The final item name (e.g., solve_advection_diffusion_1d)
    # Group 3: Optional ' as name' part
    PUB_USE_PATTERN = re.compile(
        r'^\s*pub\s+use\s+(.+?)::([a-zA-Z0-9_]+)(?:\s+as\s+([a-zA-Z0-9_]+))?\s*;(?:\s*//.*)?\s*$',
        re.IGNORECASE
    )

    # Map to track the first line index where a name was imported
    imported_names = defaultdict(list)
    
    # List to store the new, modified lines
    new_lines = []
    
    print(f"--- Analyzing '{file_path}' for name collisions ---")

    for i, line in enumerate(lines):
        match = PUB_USE_PATTERN.match(line)
        
        if match:
            full_path = match.group(1) + '::' + match.group(2)
            item_name = match.group(2)
            existing_alias = match.group(3)
            
            # Use the item_name unless an 'as' alias is already present
            exposed_name = existing_alias if existing_alias else item_name

            if exposed_name in imported_names:
                # --- COLLISION DETECTED! ---
                
                # We need to determine the module name for the prefix.
                # Example path: crate::physics::physics_sm::solve_advection_diffusion_1d
                # We want to extract 'physics_sm'
                path_parts = full_path.split('::')
                
                # Heuristic: The part before the item name is the module path.
                # Take the last segment of the path *before* the item name.
                module_name = path_parts[-2].replace('-', '_') 
                
                # Define the new alias based on the module and item name
                new_alias = f"{module_name}_{item_name}"
                
                # Construct the new line: "pub use full::path::item_name as new_alias;"
                # We ensure any existing 'as' alias is removed and replaced.
                modified_line = re.sub(
                    r'(\s*;\s*//.*)?\s*$',  # Find the semicolon and potential comment/end-of-line
                    f" as {new_alias};" + (match.group(0).split(';')[-1] if match.group(0).count(';') > 0 else ''),
                    line.split(';')[0].strip() + ";", # Only work with the path part of the line
                    1, # Replace only once
                    re.IGNORECASE
                )
                
                # Re-add trailing whitespace/newlines from original line
                trailing_content = lines[i].replace(line.split(';')[0].strip() + ";", "")
                
                new_line = modified_line.strip() + trailing_content
                
                print(f"Collision on '{exposed_name}' at line {i+1}.")
                print(f"  Old: {line.strip()}")
                print(f"  New: {new_line.strip()}")
                
                new_lines.append(new_line)
                
            else:
                # No collision, keep the line as is.
                imported_names[exposed_name].append(i)
                new_lines.append(line)
        else:
            # Not a 'pub use' line, keep it as is.
            new_lines.append(line)

    # 2. Write the cleaned content back to the file
    try:
        with open(file_path, 'w') as f:
            f.writelines(new_lines)
        
        # Calculate the number of lines that were modified
        modifications = sum(1 for original, new in zip(lines, new_lines) if original != new)
        
        print(f"--- Processing Complete ---")
        print(f"Successfully processed: '{file_path}'")
        print(f"Total lines modified: {modifications}")
        print(f"File updated in place.")
        
    except Exception as e:
        print(f"An error occurred while writing to the file: {e}")
        sys.exit(1)

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python script_name.py <path/to/prelude.rs>")
        sys.exit(1)
        
    file_path = sys.argv[1]
    auto_rename_colliding_pub_use(file_path)
