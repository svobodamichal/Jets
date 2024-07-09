#!/bin/bash

# Check if the folder is passed as an argument
if [ -z "$1" ]; then
    echo "Usage: $0 <folder>"
    exit 1
fi

# Iterate over all .list files in the specified folder
for filepath in "$1"/*.list; do
    # Extract the filename without extension
    filename=$(basename "$filepath" .list)

    # Read all lines into an array
    IFS=$'\n' read -d '' -r -a lines < "$filepath"

    # Get the total number of lines
    total_lines=${#lines[@]}

    # Calculate the midpoint
    midpoint=$(( (total_lines + 1) / 2 ))

    # Write the first half (or one more if odd) to _1.list
    printf "%s\n" "${lines[@]:0:$midpoint}" > "$1/${filename}_1.list"

    # Write the second half to _2.list
    printf "%s\n" "${lines[@]:$midpoint}" > "$1/${filename}_2.list"
done
