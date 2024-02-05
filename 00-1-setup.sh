#!/usr/bin/env bash

#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=1G
#SBATCH --time=02:00:00
#SBATCH --job-name=setup
#SBATCH --output=/PATH/TO/OUTPUT/DIR/setup.o # !! CHANGE !!
#SBATCH --error=/PATH/TO/OUTPUT/DIR/setup.e # !! CHANGE !!
#SBATCH --partition=pall

#----------------------------------------------------
# !!! BITS TO CHANGE !!!
#----------------------------------------------------
# path to your raw paired end reads
    # (in the format {NAME]_R1.fastq.gz)
rawdata_path="/full/path/to/raw/data"

# path to your assemblies
    # format {NAME}.fa.gz
assemblies_path="/full/path/to/assemblies"

# Path for the working directory
    # i.e., the folder you want to store all the dependencies of this project in
workdir="full/path/to/work/dir"

#-----------------------------------------------------

# create date structure
mkdir ${workdir}

mkdir -p ${workdir}/data/assemblies
mkdir -p ${workdir}/data/raw
mkdir ${workdir}/errorOut
mkdir ${workdir}/results
mkdir ${workdir}/scripts
mkdir ${workdir}/singularity
mkdir ${workdir}/tempdir

#--- copy in the data ---
# (if space is an issue, linking or moving the data also works, but this option has more redundancy)
cp ${rawdata_path}/* ${workdir}/data/raw
cp ${assemblies_path}/* ${workdir}/data/assemblies

#--- extract the names of the datasets ---

# Array to store file names without extensions
file_array=()

for file in ${datadir}/*R1.fastq.gz; do
    file_name=$(basename "$file")
    file_array+=("$file_name")
done

datasets_array=()
for element in "${file_array[@]}"; do
    # Use IFS (Internal Field Separator) to split the element based on underscore
    IFS="_" read -ra parts <<< "$element"
    # Add the first part (before underscore) to the new array
    datasets_array+=("${parts[0]}")
done

printf "%s\n" "${datasets_array[@]}" > ${workdir}/datasets.txt
