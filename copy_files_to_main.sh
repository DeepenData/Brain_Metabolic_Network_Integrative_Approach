#!/bin/bash

# Specify the branch names and parse the list of files from command-line arguments
source_branch="alejandro"
target_branch="main"
files=("$@")

# Switch to the source branch
git checkout "$source_branch"

# Copy the files from the source branch to the target branch
for file in "${files[@]}"; do
    git checkout "$target_branch" -- "$file"
    echo "Copied $file from the '$source_branch' branch to the '$target_branch' branch."
done

# Add a commit with the copied files
git add "${files[@]}"
git commit -m "Copied files from the '$source_branch' branch to the '$target_branch' branch"

# Pull the latest changes from the target branch
git pull origin "$target_branch"

# Push the changes to the target branch
git push origin "$target_branch"

# Switch back to the source branch
git checkout "$source_branch"


