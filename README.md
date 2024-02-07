# Metagenomic Binning with DASTool

DAS Tool is an automated method that integrates the results of a flexible number of binning algorithms to calculate an optimized, non-redundant set of bins from a single assembly. More information on DASTool can be found here https://github.com/cmks/DAS_Tool.

Binning relies on a variety of tools, which are accessed through Singularity containers for ease, modularity, and to ensure future inter-compatibility. These tools can be installed from stable dockerhub images (recomended) or from the provided .def files. However, when building from .def files on the IBU and UBELIX servers, disk quota issues can be an issue for larger programmes. 

The binning process starts with metagenomic assemblies and paired-end raw reads, and produces fasta files of consensus bins, as well as a series of graphics for interpretting the binning quality (see Results). 

## Input 
The scripts are made to work specifically with the data available from the RAS sequencing. That is compressed assemblies in the format `{accession number}.fa.gz.` and quality-controlled paired-end reads in the format `{accession_number}_R1.fastq.gz`.

## Structure 
The scripts follow a specific structure setup in step1. All scripts are designed to need minimal user input, in each script at the beginning there is a `## !! BIT TO CHANGE !!##` section. Provided the file structure follows that setup in stage 0, changing just the variables covered here will be enough. For example 
```
#--------------------------------------------------------
## !! BIT TO CHANGE !!
workdir=/PATH/TO/WORKING/DIRECTORY 
##-------------------------------------------------------
```
## 0. Setup 
### File structure
Firstly set up the file structure you want to use for the project, this can be done manually but it is recomended to use the 00-1-setup.sh script. It's important to follow this structure as many of the scripts rely on it, however, symbolic links can be used in most places (for example - to the data), if space saving is a priority. 

This script also extracts the names of the accessions you're working with and stores them in the text file `datasets.txt`. If you want to work with just specific datasets, this file can be edited manually. 

### 0.2 Software installation 
The 7 tools used in the scripts are containerised to prevent combatability and versioning issues. The `00-2-singularity_setup.sh` script downloads stable versions of these containers from docker, except for DAS_Tool which is built from scratch using the das_tool.def file (itself downloaded during the running of the script). More tools can be added from docker using the syntax:
```
{
    singularity build --tmpdir=/workdir/tempdir {TOOL_NAME}.sif docker://PATH/TO/TOOL
}
```
For metagenomics binning tools, [dockerhub nanozoo](https://hub.docker.com/u/nanozoo) is a verified repository containing many relevant tools. 

## 1. Indexing and Aligning
This script prepares the reads and assemblies for binning in subsequent steps, and takes place in 3 stages. The outputs for each dataset will be produced in the `{workdir}/results/01-indexing/{dataset}_indices` directory.

### Stage 1: Indexing with bowtie2 
`bowtie2-build` builds a Bowtie index from a set of DNA sequences. `bowtie2-build `outputs a set of 6 files with suffixes .1.bt2, .2.bt2, .3.bt2, .4.bt2, .rev.1.bt2, and .rev.2.bt2. In the case of a large index these suffixes will have a bt2l termination. These files together constitute the index: they are all that is needed to align reads to that reference. The original sequence FASTA files are no longer used by Bowtie 2 once the index is built.

The basic syntax is `bowtie2-build [options]* <reference_in> <bt2_base>`. 
Where `<reference_in>` is the reference we are using (in this case, the assembly), and `<bt2_base>` is the basename for the index, which we use <datasetnumber>_index.

### Stage 2: Aligning with bowtie2
`bowtie2` takes a Bowtie 2 index and a set of sequencing read files and outputs a set of alignments in SAM format. It is used here as so: 
`bowtie2 -p 16 -x ${dataset}_index -1 <(zcat ${reads_R1}) -2 <(zcat ${reads_R2}) -S ${outdir}/${dataset}_aligned.sam`
With the following parameters:
-p 
Launch NTHREADS parallel search threads (default: 1). Threads will run on separate processors/cores and synchronize when parsing reads and outputting alignments. Searching for alignments is highly parallel, and speedup is close to linear. Increasing -p increases Bowtie 2's memory footprint. E.g. when aligning to a human genome index, increasing -p from 1 to 8 increases the memory footprint by a few hundred megabytes. 
- -x The basename of the index for the reference genome. The basename is the name of any of the index files up to but not including the final .1.bt2 / .rev.1.bt2 / etc. bowtie2 looks for the specified index first in the current directory, then in the directory specified in the BOWTIE2_INDEXES environment variable.#
- -1 <m1>
Comma-separated list of files containing mate 1s (filename usually includes _1), e.g. -1 flyA_1.fq,flyB_1.fq. Sequences specified with this option must correspond file-for-file and read-for-read with those specified in <m2>. Reads may be a mix of different lengths. If - is specified, bowtie2 will read the mate 1s from the "standard in" or "stdin" filehandle.

- -2 <m2>
Comma-separated list of files containing mate 2s (filename usually includes _2), e.g. -2 flyA_2.fq,flyB_2.fq. Sequences specified with this option must correspond file-for-file and read-for-read with those specified in <m1>. Reads may be a mix of different lengths. If - is specified, bowtie2 will read the mate 2s from the "standard in" or "stdin" filehandle.
- -S <sam> 
File to write SAM alignments to. By default, alignments are written to the "standard out" or "stdout" filehandle (i.e. the console).

### Stage 3 Converting SAM to BAM
Our alignment is in .sam (Sequence Alignment Map) format, but most of the downstream tools take input in .bam (Binary Alignment Map) format. This is bridged with Samtools-sort, which sorts alignments by leftmost coordinates.. An appropriate @HD SO sort order header tag will be added or an existing one updated if necessary, along with the @HD SS sub-sort header tag where appropriate.

The sorted output is written to standard output by default, or to the specified file (out.bam) when -o is used. This command will also create temporary files tmpprefix.%d.bam as needed when the entire alignment data cannot fit into memory (as controlled via the -m option).
Used here as: `samtools sort -@ 16 -o ${outdir}/${dataset}_aligned.sorted.bam ${outdir}/${dataset}_aligned.sam`
Where:
- -@
 number of threads to use, by default 1.
- -o
 the file to which the output is written (default, standard output)

### Stage 4 Indexing BAM files (needed for MetaBAT2)
For MetaBAT2 (and some other tools), an index file must be created alongside the BAM file. This will take the format {dataset}_aligned.sorted.bam.bai. Downstream tools will recognise this index as long as it is in the same directory as the sorted .bam file. 

## 2. MetaBAT2 
MetaBAT2 is run in two stages. 
1. Generating a depth file
2. Running MetaBAT2

The depth file is generated with the command `jgi_summarize_bam_contig_depths <bamfile>`. By default this writes the depth file to the standard output, however we specify it to output to a .txt file in the outdir.
`jgi_summarize_bam_contig_depths --outputDepth ${dataset}_depth.txt ${bamfile}` 
- --outputDepth 
 the file to which the depth file should be written. 

Once the depth file is generated, MetaBAT2 is run with the following command:
`metabat2 -t 16 -i ${contigs} -a ${dataset}_depth.txt -o ${outdir}/${dataset}_bin`
- -t  the number of threads to use 
- -i  the input file (contigs in FASTA format)
- -a  the depth file, i.e., a file having mean and variance of base coverage depth (tab delimited; 
 the first column should be contig names, and the first row will be considered as the header and be skipped).
- -o  the output file 

### Output
Each bin discovered by MetaBAT2 will be saved to the output directory `<workdir>/results/02-MetaBAT2/<dataset>_bins` in FASTA format. 

### More info
Full documentation for MetaBAT2 can be found [here](https://gensoft.pasteur.fr/docs/MetaBAT/2.15/)

## 3. MaxBin2 
MaxBin2 relies on the depthfile generated by MetaBAT2, and as such must be run after MetaBAT2. 
It is run in a single stage like so:
`run_MaxBin.pl -contig ${contigs} -abund ${depthfile} -out ${outdir}/${dataset}_maxbin2 -thread 16`
- -contig  the assembled contigs to be binned.
- -abund  the depthfile 
- -out  the output directory 
- -thread specfies the number of threads 

### Output
Each bin discovered by MaxBin2 will be saved to the output directory `<workdir>/results/03-MaxBin2/<dataset>_bins` in FASTA format, much like MetaBAt2.

### More info
Full documentation for MaxBin2 can be found [here](https://sourceforge.net/projects/maxbin2/files/)

## 4. CONCOCT 
CONCOCT “bins” metagenomic contigs. Metagenomic binning is the process of clustering sequences into clusters corresponding to operational taxonomic units of some level. It's implimentation is significantly more convoluted than that of the previous two binning tools and takes place in 5 stages.
1. Cutting up contigs into smaller parts 
` cut_up_fasta.py <( zcat ${contigs}) -c 10000 -o 0 --merge_last -b contigs_10k.bed > contigs_10k.fa`
- -c   chunk size to cut the contigs into
- -o   overlap size for the contigs 
- -b  he bed file to be created. 

2. Generating Coverage Table
Generate table with coverage depth information per sample and subcontig. This step assumes the directory ‘mapping’ contains sorted and indexed bam files where each sample has been mapped against the original contigs
`concoct_coverage_table.py contigs_10k.bed ${bamfile} > coverage_table.tsv`

3. Running CONCOCT
` concoct -t 8 --composition_file contigs_10k.fa --coverage_file coverage_table.tsv -b concoct_output/`
- -t number of threads
- -b Specify the basename for files or directory where outputwill be placed. Path to existing directory or basenamewith a trailing '/' will be interpreted as a directory.If not provided, current directory will be used.

4. Merge subcontig clustering into original contig clustering.
` merge_cutup_clustering.py concoct_output/clustering_gt1000.csv > concoct_output/clustering_merged.csv`

5. Extracting the bins into individual FASTA
`extract_fasta_bins.py <(zcat ${contigs}) concoct_output/clustering_merged.csv --output_path concoct_output/fasta_bins`

### Output
Unlike the previous two programmes, each bin discovered by CONCOCT will be saved to the output directory `<workdir>/results/04-CONCOCT/547_bins/concoct_output/fasta_bins` in .fa format.

### More info
Full documentation for CONCOCT can be found [here](https://concoct.readthedocs.io/en/latest/)

## 5 DASTool
DASTool relies on a TSV file with the first column corresponding the contig IDs and the second collumn corresponding to the bin in which that contig is placed. This can be generated from the FASTA bins with a simple bash script `das_tool_helper_script.sh `
### 1 DASPrep
The 05-1-DASPrep.sh downloads the helper script from github, and uses it to generate tsv files for each of the binning tool results. The syntax to run it is as follows:
`bash das_tool_helper_script.sh -e fa -i ${bins_dir} > ${bins_dir}/${dataset}.tsv`  
- -e  extension of the fasta files (fa for CONCOCT and MetaBAT2, fasta for MaxBin2)
- -i the input directory where the fasta bins are stored 
This can be executed for the results of any binning tool, so long as the bins are in fasta format in a single directory, and the extensions and paths are specified properly. 

### 2 Running DASTool
DAS Tool is an automated method that integrates the results of a flexible number of binning algorithms to calculate an optimized, non-redundant set of bins from a single assembly.
`DAS_Tool -i ${bins_metabat},${bins_maxbin},${bins_conoct} -l metabat2,maxbin2,concoct -c ${tmp_contigs} -o ${dataset}_DASTool --write_bins --threads 16`
- -i --bins=<contig2bin>                   Comma separated list of tab separated contigs to bin tables.
- -l --labels=<labels>                     Comma separated list of binning prediction names.
- -c --contigs=<contigs>                   Contigs in fasta format.
- -o --outputbasename=<outputbasename>     Basename of output files.
- -t --threads=<threads>                   Number of threads to use [default: 1].
- --write_bins                             Export bins as fasta files.

DASTool does not take the contigs in compressed (.gz) format, and due to issues with the singularity container using zcat within the container causes issues. For this reason, in the script, the contigs are temporarily unzipped into the temp directory, and then removed after each iteration. This prevents memory from becoming an issue. 

### Output 
The consensus bins from DASTool are output to `<workdir>/results/05-DASTool/<dataset>/<dataset>_DASTool_DASTool_bins.
The output directory also contains information about the consensus bins including:
- Summary of output bins including quality and completeness estimates (*_DASTool_summary.tsv).
- Contigs to bin file of output bins (*_DASTool_contigs2bin.tsv).

### More info
Full documentation for DASTool can be found [here] (https://github.com/cmks/DAS_Tool)

# Section 6: CheckM

## Overview
CheckM is utilized in this section to assess the quality and completeness of the bins generated by various binning tools, including MetaBAT2, MaxBin2, CONCOCT, and DASTool. It evaluates the bins based on lineage-specific marker genes, providing insights into their taxonomic composition and potential completeness.

### Tools and Resources
- **CheckM Version:** 1.0.13
- **Singularity Image Source:** [Nanozoo Docker Hub](https://hub.docker.com/r/nanozoo/checkm)
- **Software Source:** [CheckM Documentation](https://ecogenomics.github.io/CheckM/)

### 1. MetaBAT2
- Calculates quality and completeness metrics using CheckM lineage_wf.
- Generates quality assessment plots in PDF format.

### 2. MaxBin2
- Performs quality assessment and completeness estimation using CheckM lineage_wf.
- Produces quality assessment plots for visualization.

### 3. CONCOCT
- Evaluates bin quality and completeness with CheckM lineage_wf.
- Generates plots to visualize the assessment results.

### 4. DASTool
- Assesses bin quality and completeness using CheckM lineage_wf.
- Produces quality assessment plots for comprehensive evaluation.

## Output
For each dataset, the results are organized into dedicated directories within the CheckM output directory (`${workdir}/results/06-checkM/${dataset}`). The directories contain quality assessment plots and summary files providing detailed information on bin quality and completeness.

### Important Note
Ensure to specify the temporary directory (`--tmpdir`) with sufficient space and memory allocation for CheckM operations, as default locations may not provide adequate resources.

