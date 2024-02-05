#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --job-name=DASTool
#SBATCH --output=/!!WORKDIR!!/errorOut/dastool_out.o  # < CHANGE
#SBATCH --error=/!!WORKDIR!!/errorOut/dastool_error.e # < CHANGE
#-------------------------------------------------------------
# STEP 5.1 -- DASTool Run
# tools: DASTool v1.1.7
# image source: https://github.com/cmks/DAS_Tool/blob/master/DAS_Tool (built from .def file)
# software source: https://github.com/cmks/DAS_Tool
#-------------------------------------------------------------
## !! BIT TO CHANGE !! 
workdir=/path/to/workdir                              # < CHANGE
##-------------------------------------------------------
# 0 --- setup 
tmpdir="${workdir}/tempdir"
dastool_sif="${workdir}/singularity/das_tool.sif"
datasets_array=($(<datasets.txt))


# loop through datasets
for dataset in "${datasets_array[@]}"; do
    tmp_contigs=${tmpdir}/temp_${dataset}_contigs.fa
    outdir=${workdir}/results/05-DASTool/${dataset}
    contigs=${workdir}/data/assemblies/${dataset}.fa.gz
    bins_metabat="${workdir}/results/02-MetaBAT2/${dataset}_bins/${dataset}.tsv"
    bins_maxbin="${workdir}/results/03-MaxBin2/${dataset}_bins/${dataset}.tsv"
    bins_conoct="${workdir}/results/04-CONCOCT/${dataset}_bins/concoct_output/fasta_bins/${dataset}.tsv"

    #--- make a temporary contigs file to deal with the compression issue
    zcat ${contigs} > ${tmp_contigs}

    mkdir -p ${outdir}
    cd ${outdir}
    # 2 --- Run DAS_Tool 
        #NB. it's very sensitive to spaces!
    singularity exec \
    --bind ${workdir} \
    --bind ${outdir} \
    ${dastool_sif} DAS_Tool -i ${bins_metabat},${bins_maxbin},${bins_conoct} \
    -l metabat2,maxbin2,concoct \
    -c ${tmp_contigs} \
    -o ${dataset}_DASTool \
    --threads 16 

    # remove the temporart contigs file
    rm ${tmp_contigs}
done
