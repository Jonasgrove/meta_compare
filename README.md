# Meta Compare

The **Meta Compare** pipeline is based in Snakemake and analyzes a multi-sample metagenomic dataset of both Hi-C and shotgun reads, ultimately calculating the percentage abundances of species present in each sample.

## More information

To see a completed example of this pipeline, go to our main repository at https://github.com/Sam-Koehler/HiC-Metagenomics.

## Usage guide

```
./meta_compare.sh -h 

meta_compare: Hi-C deconvolution pipeline
 
Note: This is an experimental pipeline, meant for exploratory use only.
      If you are working on a project of any amount of importance,
      it is highly recommended that you use the ProxiMeta platform
      available through Phase Genomics. https://phasegenomics.com/products/proximeta/
  
usage: ./meta_compare.sh -s <sample1,sample2,sampleN> -i <input_directory> -o <output_directory> -t <threads> -g <gtdbtk_reference_data> [-O True/False] [-R path/to/reference_genome]
 
Options:
-s    --sample_names      <sample1,sample2,sampleN>   comma separated base sample names.
-i    --input_directory   <directory name>            path to input directory
-o    --output_directory  <directory name>            path to output directory
-t    --threads           <thread number>             number of cpu cores to use
-g    --gtdbtk_ref        <directory name>            path to gtdbtk reference data
-O    --omit_genome       <True/False>                True or False option where True indicates that a dominating genome in the sample will be omitted and False specifies use of the entire provided sample in the analysis. a reference genome must be provided using the -R option if True option is chosen
-R    --reference_genome  <directory name>            path to reference genome, which will be omitted from all provided samples. Useful if there is a known organism that is dominating a metagenomic sample.
```
