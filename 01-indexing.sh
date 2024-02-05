#!/usr/bin/env bash

#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=1G
#SBATCH --time=1-00:00:00
#SBATCH --job-name=indexing
#SBATCH --output=/!!workdir!!/errorOut/indexing.o
#SBATCH --error=/!!workdir!!/indexing.e

#-------------------------------------------------------------
# SCRIPT 1 - alignment 
# tools: Bowtie2-2.5.1
#        Samtools 1.10 
#-------------------------------------------------------------
## !!! CHANGE THIS BIT !!!
## AND THE --output AND --error VARIABLES 
workdir=/PATH/TO/WORK/DIR
#-------------------------------------------------------------

# 0 --- setup 
# set variables

raw_dir="${workdir}/data/raw"
assembly_dir="${workdir}/data/assemblies"
outdir="${workdir}/results/01-indexing"
mkdir -p ${outdir}
cd ${workdir}
datasets_array=($(<datasets.txt))

cd ${outdir}

# Loop through datasets 
for dataset in "${datasets_array[@]}"; do
    contigs="${assemblies_dir}/${dataset}.fa.gz"
    reads_R1="${raw_dir}/${dataset}_R1.fastq.gz"
    reads_R2="${raw_dir}/${dataset}_R2.fastq.gz"

# Create output folder 
    mkdir ${outdir}/${dataset}_indices
    cd ${outdir}/${dataset}_indices

    # 1 --- indexing 
    singularity exec \
    --bind ${outdir} \
    --bind ${contigs} \
    --bind ${workdir} \
    ${singularity_dir}/${bowtie2.sif} bowtie2-build  ${contigs} ${dataset}_index 

    # 2 --- aligning 

    singularity exec \
    --bind ${outdir} \
    --bind ${contigs} \
    --bind ${workdir} \
    ${singularity_dir}/${bowtie2.sif} bowtie2 -p 16 -x ${dataset}_index -1 <(zcat ${reads_R1}) -2 <(zcat ${reads_R2}) -S ${outdir}/${dataset}_aligned.sam
        # -p 'performance' specifies number of threads 
        # -x identifies the index
        # -1 -2 specify the forward and backwards mates
        # -S specifies .sam output format 

    # 3 --- SAM to BAM
    singularity exec \
    --bind ${outdir} \
    --bind ${workdir} \
    ${singularity_dir}/${samtools.sif} samtools sort -@ 16 -o ${outdir}/${dataset}_aligned.sorted.bam ${outdir}/${dataset}_aligned.sam
        # -o specifies output file/format
        # -@ specifies number of threads (16)

    # 4 --- indexing .bam file (for Metabat2 input)
    singularity exec \
    --bind ${outdir} \
    --bind ${workdir} \
    ${singularity_dir}/${samtools.sif}samtools index ${outdir}/${dataset}_aligned.sorted.bam 
done
