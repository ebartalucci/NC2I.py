#!/bin/bash

# BASH script for running Turbomole in automatic mode
# Ettore Bartalucci, Aachen 26.10.23


# Check if the user is in the correct directory
if [ ! -d "your_parent_folder" ]; then
    echo "Please navigate to the correct parent folder."
    exit 1
fi

# Loop through all subfolders
for folder in */; do
    if [ -d "$folder" ]; then
        cd "$folder" || exit 1

        # Search for an XYZ file in the current directory
        xyz_file=$(find . -maxdepth 1 -type f -name "*.xyz" | head -1)

        if [ -n "$xyz_file" ]; then
            # Execute the commands in the specified order
            vi INPUT
            ~/gcoord/gcoord.x < "$xyz_file"
            mv GCOORD coord

            # introduce a 1 minute delay
            sleep 60

            # run geometry optimization
            sbatch geoopt.in
        else
            echo "No XYZ file found in $folder"
        fi

        cd ..
    fi
done



# Usage for Linux HPC

# make it executable
#chmod +x turbomole_auto_run.sh

# run from cmd
#./turbomole_auto_run.sh