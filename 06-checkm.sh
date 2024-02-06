#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=06:00:00
#SBATCH --job-name=t-checkM
#SBATCH --output=/!!WORKDIR!!/errorOut/checkM_out.o  # < CHANGE
#SBATCH --error=/!!WORKDIR!!/errorOut/checkM_error.e # < CHANGE
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=user@unibe.ch                    # < CHANGE
#-------------------------------------------------------------
#-------------------------------------------------------------
# STEP 6 -  - CheckM
# tools: CheckM 1.0.13
# image source: https://hub.docker.com/r/nanozoo/checkm
# software source: https://ecogenomics.github.io/CheckM/
#-------------------------------------------------------------
## !! BIT TO CHANGE !! 
workdir=/path/to/workdir                              # < CHANGE
##-------------------------------------------------------
# 0 --- setup 
checkm_sif="${workdir}/singularity/checkm.sif"
tmpdir="${workdir}/tempdir"
datasets_array=($(<datasets.txt))

# 1 --- MetaBat2

# loop through datasets:
for dataset in "${datasets_array[@]}"; do

# set variables
bins_dir="${workdir}/results/02-MetaBAT2/${dataset}_bins" # Change for other binning outputs 
outdir="${workdir}/results/06-checkM/${dataset}/metabat2" 


mkdir -p ${outdir}
cd ${outdir}

singularity exec \
--bind ${outdir} \
--bind ${bins_dir} \
--bind ${workdir} \
${checkm_sif} checkm lineage_wf -x fasta ${bins_dir} ${outdir} --reduced_tree -t 16 --tmpdir ${tmpdir} 
# IMPORTANT: make sure you specify the tmpdir to a directory you have permission to write to and
# that has a large amount of memory, usually in /data for IBU and /storage for UBELIX 
# the default will be /home/tmp which is NOT ENOUGH 
# generate plots
plots_dir="${outdir}/plots"
mkdir ${plots_dir}

singularity exec \
--bind ${outdir} \
--bind ${bins_dir} \
--bind ${workdir} \
--bind ${plots_dir} \
${checkm_sif} checkm bin_qa_plot --image_type pdf -x fasta ${outdir} ${bins_dir} ${plots_dir}
done

# 2 --- MaxBin2

# loop through datasets:
for dataset in "${datasets_array[@]}"; do

# set variables
bins_dir="${workdir}/results/03-MaxBin2/${dataset}_bins" # Change for other binning outputs 
outdir="${workdir}/results/06-checkM/${dataset}/maxbin2" 


mkdir -p ${outdir}
cd ${outdir}

singularity exec \
--bind ${outdir} \
--bind ${bins_dir} \
--bind ${workdir} \
${checkm_sif} checkm lineage_wf -x fasta ${bins_dir} ${outdir} --reduced_tree -t 16 --tmpdir ${tmpdir} 
# IMPORTANT: make sure you specify the tmpdir to a directory you have permission to write to and
# that has a large amount of memory, usually in /data for IBU and /storage for UBELIX 
# the default will be /home/tmp which is NOT ENOUGH 
# generate plots
plots_dir="${outdir}/plots"
mkdir ${plots_dir}

singularity exec \
--bind ${outdir} \
--bind ${bins_dir} \
--bind ${workdir} \
--bind ${plots_dir} \
${checkm_sif} checkm bin_qa_plot --image_type pdf -x fasta ${outdir} ${bins_dir} ${plots_dir}
done

# 3 --- CONCOCT

# loop through datasets:
for dataset in "${datasets_array[@]}"; do

# set variables
bins_dir="${workdir}/results/04-CONCOCT/${dataset}_bins" # Change for other binning outputs 
outdir="${workdir}/results/06-checkM/${dataset}/concoct" 


mkdir -p ${outdir}
cd ${outdir}

singularity exec \
--bind ${outdir} \
--bind ${bins_dir} \
--bind ${workdir} \
${checkm_sif} checkm lineage_wf -x fasta ${bins_dir} ${outdir} --reduced_tree -t 16 --tmpdir ${tmpdir} 
# IMPORTANT: make sure you specify the tmpdir to a directory you have permission to write to and
# that has a large amount of memory, usually in /data for IBU and /storage for UBELIX 
# the default will be /home/tmp which is NOT ENOUGH 
# generate plots
plots_dir="${outdir}/plots"
mkdir ${plots_dir}

singularity exec \
--bind ${outdir} \
--bind ${bins_dir} \
--bind ${workdir} \
--bind ${plots_dir} \
${checkm_sif} checkm bin_qa_plot --image_type pdf -x fasta ${outdir} ${bins_dir} ${plots_dir}
done

# 4 --- DASTool

# loop through datasets:
for dataset in "${datasets_array[@]}"; do

# set variables
bins_dir="${workdir}/results/05-DASTool/${dataset}/${dataset}_DASTool_DASTool_bins" # Change for other binning outputs
outdir="${workdir}/results/06-checkM/${dataset}/das_tool" # Change for other binning outputs

mkdir -p ${outdir}
cd ${outdir}

singularity exec \
--bind ${outdir} \
--bind ${bins_dir} \
--bind ${workdir} \
--bind ${plots_dir} \
${checkm_sif} checkm bin_qa_plot --image_type pdf -x fa ${outdir} ${bins_dir} ${plots_dir}
done
