#!/usr/bin/env python3

import os
import re
import sys

def replace_functions_in_file(file_path):
    with open(file_path, 'r') as file:
        content = file.read()

    # Replace F77_ prefixed functions with corresponding mpi_ functions
    updated_content = re.sub(r'F77_MPI_([A-Z_]+)', lambda m: f"mpi_{m.group(1).lower()}_", content)

    with open(file_path, 'w') as file:
        file.write(updated_content)

# Replace files in the current directory (or specify your directory)
# directory = './'
# for root, _, files in os.walk(directory):
#     for file in files:
#         if file.endswith('.c') or file.endswith('.h'):  # Adjust file types as needed
#             replace_functions_in_file(os.path.join(root, file))


def main():
    if len(sys.argv) != 2:
        print("Usage: refactorer.py <file>")
        sys.exit(1)
    replace_functions_in_file(sys.argv[1])

if __name__ == '__main__':
    main()
