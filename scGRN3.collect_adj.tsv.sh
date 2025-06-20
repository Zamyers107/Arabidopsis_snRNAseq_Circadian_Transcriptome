#!/bin/bash

# Define target directory
TARGET_DIR="all_adj"
mkdir -p "$TARGET_DIR"

# Find all files named adj.tsv
find . -type f -name "adj.tsv" | while read -r filepath; do
    # Extract the cluster from the top-level folder name.
    # Assumes paths like ./GRNs.c0.ClusterOnlyCycling.c0OnlyCells/...
    cluster=$(echo "$filepath" | sed -r 's|^\./GRNs\.(c[^.]+).*|\1|')
    
    # Extract the identity from the immediate subdirectory.
    # If the file is in "exprMat_scenic", it's unshuffled.
    # If it's in a folder like "exprMat.shuffle3_scenic", extract the shuffle number.
    subdir=$(basename "$(dirname "$filepath")")
    if [ "$subdir" = "exprMat_scenic" ]; then
       identity="unshuffled"
    elif [[ "$subdir" =~ exprMat\.shuffle([0-9]+)_scenic ]]; then
       identity="shuffle${BASH_REMATCH[1]}"
    else
       identity="$subdir"
    fi

    # Construct the new filename
    new_filename="${cluster}.${identity}.adj.tsv"

    # Copy the file to the target directory with the new name
    cp "$filepath" "${TARGET_DIR}/${new_filename}"
    echo "Copied $filepath as ${TARGET_DIR}/${new_filename}"
done

