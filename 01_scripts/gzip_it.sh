# Change directory to the h5ad_convert directory
cd /data/williamsdrw/h5ad_convert/

# Loop through all subdirectories
for subdirectory in */; do

    # Change directory to the current subdirectory
    cd $subdirectory

    # Get a list of all files in the current directory
    files=$(ls)

    # Loop through all files in the current directory
    for file in $files; do

        # Skip the metadata.csv file
        if [[ $file == "metadata.csv" ]]; then
            continue
        fi

        # Gzip the file
        gzip $file
    done

    # Change directory back to the parent directory
    cd ..
done

# Change directory back to the current directory
cd $current_directory

# Print a message
echo "All files have been gzipped."

