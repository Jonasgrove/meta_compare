#!/bin/bash

Help()
{
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
}

while getopts s:i:o:t:g:O:R:h flag; do
	case "${flag}" in
		h) Help
			exit;; 
		s) sample_names=${OPTARG};;
		i) input_directory=${OPTARG};;
		o) output_directory=${OPTARG};;
		t) threads=${OPTARG};;
		g) gtdbtk_ref=${OPTARG};;
		O) omit_genome=${OPTARG};;
		R) reference_genome=${OPTARG};;
	esac
done


SAMPLES=$sample_names INPUT=$input_directory OUTPUT=$output_directory THREADS=$threads OMIT=$omit_genome REF=$reference_genome GTDBTK_REF=$gtdbtk_ref snakemake # --cores $4 --use-conda

