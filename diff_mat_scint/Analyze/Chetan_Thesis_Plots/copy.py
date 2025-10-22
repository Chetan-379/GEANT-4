import os
import shutil

def copy_index_php_to_subdirectories(target_directory='.', skip_keyword='res'):
    # Get the absolute path of the target directory
    target_directory = os.path.abspath(target_directory)

    # Get the absolute path of the source file (index.php)
    source_file_path = os.path.join(os.path.dirname(__file__), 'index.php')

    # Ensure the source file exists
    if not os.path.exists(source_file_path):
        print("Error: Source file 'index.php' does not exist.")
        return

    # Iterate through subdirectories and copy index.php
    for subdir, _, _ in os.walk(target_directory):
        # Skip directories that contain the specified keyword
        if skip_keyword in subdir:
            print(f"Skipping '{subdir}'")
            continue

        # Construct the target path for the copy
        target_path = os.path.join(subdir, 'index.php')

        # Check if the target path is a directory
        if os.path.isdir(subdir):
            # Skip if source and target paths are the same
            if os.path.exists(target_path) and os.path.samefile(source_file_path, target_path):
                continue

            # Create the directory if it does not exist
            os.makedirs(subdir, exist_ok=True)

            # Copy the file
            shutil.copy2(source_file_path, target_path)
            print(f"Copied 'index.php' to '{target_path}'.")
        else:
            print(f"Skipping '{subdir}' as it is not a directory.")

# Example usage:
copy_index_php_to_subdirectories('.', skip_keyword='res')
