import os
from datetime import datetime

# --- Configuration ---
# The target directory to start the search from
TARGET_DIR = 'src'
# The output file name
OUTPUT_FILE = f'rust_file_list_{datetime.now().strftime("%Y%m%d_%H%M%S")}.txt'
# --- End Configuration ---

def list_rust_files_formatted(target_dir, output_file):
    """
    Finds all .rs files in the target_dir and writes the file path
    to the output_file in the format [ ] [ ] [ ] /path/to/file.rs
    """
    if not os.path.isdir(target_dir):
        print(f"❌ Error: Directory '{target_dir}' not found. Please ensure it exists.")
        return

    all_rs_files = []
    print(f"🔎 Scanning '{target_dir}' to find all .rs files...")
    
    # os.walk generates the file names in a directory tree
    for root, _, files in os.walk(target_dir):
        for file in files:
            if file.endswith('.rs'):
                # Store the full, absolute path of the file
                full_path = os.path.join(root, file)
                all_rs_files.append(full_path)

    total_files = len(all_rs_files)
    print(f"✅ Found {total_files} .rs files.")

    if total_files == 0:
        print("Nothing to process. Exiting.")
        return

    # Process and write to file
    
    # Open the output file for writing
    with open(output_file, 'w') as f:
        print(f"\n📝 Writing output to '{output_file}'...")
        
        # Define the static prefix for all lines
        STATIC_PREFIX = '[ ] [ ] [ ] '
        
        for file_path in all_rs_files:
            
            # To get the desired relative path starting from the TARGET_DIR:
            # - Use os.path.relpath(path, start=os.path.dirname(target_dir))
            # We assume 'src' is inside the current working directory.
            relative_path = os.path.relpath(file_path, start=os.path.dirname(target_dir))
            
            # Construct the line: [ ] [ ] [ ] /src/path/to/file.rs
            output_line = STATIC_PREFIX + '/' + relative_path
            
            # Write the formatted line to the file
            f.write(output_line + '\n')

    print(f"\n🎉 Successfully processed {total_files} files.")
    print(f"The results are saved in: **{output_file}**")

if __name__ == "__main__":
    list_rust_files_formatted(TARGET_DIR, OUTPUT_FILE)
