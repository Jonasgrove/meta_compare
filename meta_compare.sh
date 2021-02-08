#!/bin/bash

# 1 = sample names
# 2 = input directory
# 3 = output directory
# 4 = threads to use
# 5 = boolean filter a genome
# 6 = reference to omit 
# 7 = path to gtdbtk reference data


SAMPLES=$1 INPUT=$2 OUTPUT=$3 THREADS=$4 OMIT=$5 REF=$6 GTDBTK_REF=$7 snakemake --cores $4 --use-conda
