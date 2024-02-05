#!/usr/bin/env bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --job-name=singularity_setup
#SBATCH --output=/PATH/TO/OUTPUT/FILE
#SBATCH --error=/PATH/TO/ERROR/FILE
#--------------------------------------------------------
# Singularity setup of software 
# NB. replace paths with the relevant paths for your setup
#--------------------------------------------------------
singularity_dir=PATH/TO/SINGULARITY/DIR

cd ${singularity_dir}

#1 samtools
singularity build samtools.sif docker://pegi3s/samtools_bcftools

#2 Bowtie2 
singularity build bowtie2.sif docker://biocontainers/bowtie2

#3 Metabat2
singularity build metabat2.sif docker://nanozoo/metabat2

#4 Maxbin2 
singularity build maxbin2.sif docker://nanozoo/maxbin2

#5 CONCOCT 
singularity build concoct.sif docker://nanozoo/concoct 

#6 DAS_Tool 

singularity build das_tool.sif docker://shengwei/das_tool 

#7 checkm 
singularity build checkm.sif docker://nanozoo/checkm
