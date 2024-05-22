#!/bin/bash

# Check if filename is provided
if [ $# -eq 0 ]
then
    echo "No filename provided. Usage: ./script.sh filename"
    exit 1
fi

# Filename
filename=$1

# Create a temporary file
temp_file=$(mktemp)

# Process the file
while IFS= read -r line
do
    # If line contains ==, keep it as is
    if [[ $line == *"=="* ]]; then
        echo "$line" >> "$temp_file"
    else
        # Remove everything after the second = if it is separated from the first =
        new_line=$(echo "$line" | sed 's/^\([^=]*=[^= ]*\)=.*/\1/')
        echo "$new_line" >> "$temp_file"
    fi
done < "$filename"

# Replace the original file with the processed file
mv "$temp_file" "$filename"