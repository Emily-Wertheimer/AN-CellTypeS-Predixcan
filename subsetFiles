#!/bin/bash

#ensure you are in home directory (in this case: /home/ekw28/AN-CellTypeS-Predixcan/freeze3_cellTypes_repo)

# Set the directory path for models
directory='/gpfs/gibbs/pi/huckins/girgenti_celltype_models'

# Set the desired filename pattern
fileNameWant='Excitatory_Neurons'

# Define the output directory
outpath="/home/ekw28/AN-CellTypeS-Predixcan/freeze3_cellTypes_repo/${fileNameWant}"

# Create the output directory (uncomment if it doesn't exist)
# mkdir -p "$outpath"

# Loop through files in the directory
for file in "$directory"/*; do
    # Check if the filename contains the desired element
    if [[ "$file" == *"${fileNameWant}"* ]]; then
        # Copy the file to the output directory
        cp "$file" "$outpath"
    fi
done
