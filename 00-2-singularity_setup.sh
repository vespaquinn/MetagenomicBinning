#!/usr/bin/env bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --job-name=singularity_setup
#SBATCH --output=/PATH/TO/OUTPUT/FILE
#SBATCH --error=/PATH/TO/ERROR/FILE
#--------------------------------------------------------
# Singularity setup of software 
#--------------------------------------------------------
## !! BIT TO CHANGE !!
workdir=/PATH/TO/WORKING/DIRECTORY 
##-------------------------------------------------------
#--- 0 - Setup
singularity_dir=${workdir}/singularity
tempdir=${workdir}/tempdir
cd ${singularity_dir}
cd ${workdir}/singularity

#--- 1 - Getting .def files 
wget https://raw.githubusercontent.com/vespaquinn/MetagenomicBinning/main/helper_files/das_tool.def

#---  2 - building containers 
# - i samtools
singularity build --tmpdir=${tempdir} samtools.sif docker://pegi3s/samtools_bcftools

#- ii  Bowtie2 
singularity build --tmpdir=${tempdir} bowtie2.sif docker://biocontainers/bowtie2

#- iii Metabat2
singularity build --tmpdir=${tempdir} metabat2.sif docker://nanozoo/metabat2

#- iv Maxbin2 
singularity build --tmpdir=${tempdir} maxbin2.sif docker://nanozoo/maxbin2

#- v CONCOCT 
singularity build --tmpdir=${tempdir} concoct.sif docker://nanozoo/concoct 

#- vi DAS_Tool 

singularity build --tmpdir=${tempdir} das_tool.sif das_tool.def

#- vii checkm 
singularity build --tmpdir=${tempdir} checkm.sif docker://nanozoo/checkm
