#!/bin/bash

# Check if the user provided a directory path as an argument
if [ $# -lt 1 ]; then
    echo "Usage: $0 <output_file> [directory1] [directory2] [directory3] ..."
    exit 1
fi

output_file="$1"
shift   # Remove the first argument (output_file) from the list of arguments

# Loop through the remaining arguments, which are directory names
for directory in "$@"; do
    # Use the 'find' command to list all files in the specified directory and append to the output file
    find "$directory" -type f >> "$output_file"
done