#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --job-name=DAS_Prep
#SBATCH --output=/!!WORKDIR!!/errorOut/dasprep_out.o  # < CHANGE
#SBATCH --error=/!!WORKDIR!!/errorOut/dasprep_error.e # < CHANGE
#-------------------------------------------------------------
# STEP 5.1 -- DASTool Prep
# tools: dastool helper script 
# source: https://github.com/cmks/DAS_Tool/tree/master/src
#-------------------------------------------------------------
## !! BIT TO CHANGE !! 
workdir=/path/to/workdir                              # < CHANGE
##-------------------------------------------------------


#--- 1 --- preparation
# das tool requires a tsv input which we must generate from the binning directories 
# the helper script for this can be downloaded here 
wget -O ${workdir}/scripts/das_tool_helper_script.sh https://raw.githubusercontent.com/cmks/DAS_Tool/master/src/Fasta_to_Contig2Bin.sh 

    # i MetaBAT2
for dataset in "${datasets_array[@]}"; do
bins_dir="${workdir}/results/02-MetaBAT2/${dataset}_bins"
bash das_tool_helper_script.sh -e fa -i ${bins_dir} > ${bins_dir}/${dataset}.tsv
done 

    # ii Maxbin2
for dataset in "${datasets_array[@]}"; do
bins_dir="${workdir}/results/03-MaxBin2/${dataset}_bins"
bash das_tool_helper_script.sh -e fasta -i ${bins_dir} > ${bins_dir}/${dataset}.tsv
done 

    # iii concoct
for dataset in "${datasets_array[@]}"; do
bins_dir="${workdir}/results/04-concoct/${dataset}_bins/concoct_output/fasta_bins"
bash das_tool_helper_script.sh -e fa -i ${bins_dir} > ${bins_dir}/${dataset}.tsv
done 

    # this can be continued for any further binning tools included (see README)
