#!/bin/bash

# 1 = sample names
# 2 = input directory
# 3 = output directory
# 4 = threads to use
# 5 = boolean filter a genome
# 6 = reference to omit 
# 7 = path to gtdbtk reference data

if [[ $1 == "-h" ]]
then
    echo " "
    echo "meta_compare: Hi-C deconvolution pipeline"
    echo " "
    echo "Note: This is an experimental pipeline, meant for exploratory use only."
    echo "      If you are working on a project of any amount of importance,"
    echo "      it is highly recommended that you use the ProxiMeta platform"
    echo "      available through Phase Genomics. https://phasegenomics.com/products/proximeta/"
    echo "  "
    echo "usage: ./meta_compare.sh -s <sample1,sample2,sampleN> -i <input_directory> -o <output_directory> -t <threads> -g <gtdbtk_reference_data> [-O True/False] [-R path/to/reference_genome]"
    echo " "
    echo "Options:"
    echo "-s    --sample_names      <sample1,sample2,sampleN>   comma separated base sample names."
    echo "-i    --input_directory   <directory name>            path to input directory"
    echo "-o    --output_directory  <directory name>            path to output directory"
    echo "-t    --threads           <thread number>             number of cpu cores to use"
    echo "-g    --gtdbtk_ref        <directory name>            path to gtdbtk reference data"
    echo "-O    --omit_genome       <True/False>                True or False option where True indicates that a dominating genome in the sample will be omitted and False specifies use of the entire provided sample in the analysis. a reference genome must be provided using the -R option if True option is chosen"
    echo "-R    --reference_genome  <directory name>            path to reference genome, which will be omitted from all provided samples. Useful if there is a known organism that is dominating a metagenomic sample."           
else
    SAMPLES=$1 INPUT=$2 OUTPUT=$3 THREADS=$4 OMIT=$5 REF=$6 GTDBTK_REF=$7 snakemake --cores $4 --use-conda
fi