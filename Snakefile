import os

# retrieve environment variables
SAMPLE     = os.environ["SAMPLES"].split(",")
input_dir  = os.environ["INPUT"]
output_dir = os.environ["OUTPUT"]
threads    = int(os.environ["THREADS"])
omit       = os.environ["OMIT"]
reference  = os.environ["REF"]
gtdbtk_data= os.environ["GTDBTK_REF"]


# final target rule
rule all:
    input:
        expand(output_dir + "/quast_out/{samp}_sho", samp=SAMPLE),   ## Quast
        expand(output_dir + "/gtdbtk_out/{samp}/", samp=SAMPLE),     ## gtdbtk
        expand(output_dir + "/checkm_out/{samp}/", samp=SAMPLE),     ## checkm

#quality trim reads using fastp
#fastp version 0.20.1

rule fastp:
    input:
        # sho
        sho_r1 = expand(input_dir + "/{samp}_sho_R1.fastq", samp=SAMPLE),
        sho_r2 = expand(input_dir + "/{samp}_sho_R2.fastq", samp=SAMPLE),

        # hic
        hic_r1 = expand(input_dir + "/{samp}_hic_R1.fastq", samp=SAMPLE),
        hic_r2 = expand(input_dir + "/{samp}_hic_R2.fastq", samp=SAMPLE)

    output:
        # sho
        sho_r1 = output_dir + "/fastp/{samp}/sho/{samp}_sho_R1.fastq",       ## R1 Output
        sho_r2 = output_dir + "/fastp/{samp}/sho/{samp}_sho_R2.fastq",       ## R2 Output
        html_sho = output_dir + "/fastp/{samp}/sho/{samp}_summary.html",     ## HTML Summary
        json_sho = output_dir + "/fastp/{samp}/sho/{samp}_summary.json",     ## JSON Summary

        # sho
        hic_r1 = output_dir + "/fastp/{samp}/hic/{samp}_hic_R1.fastq",       ## R1 Output
        hic_r2 = output_dir + "/fastp/{samp}/hic/{samp}_hic_R2.fastq",       ## R2 Output
        html_hic = output_dir + "/fastp/{samp}/hic/{samp}_summary.html",     ## HTML Summary
        json_hic = output_dir + "/fastp/{samp}/hic/{samp}_summary.json",     ## JSON Summary

    conda:
        "envs/env_1.yaml"
    shell:
        """
        /usr/bin/time -v fastp -i {input.sho_r1} -I {input.sho_r2} -o {output.sho_r1} -O {output.sho_r2} -h {output.html_sho} -j {output.json_sho}
        /usr/bin/time -v fastp -i {input.hic_r1} -I {input.hic_r2} -o {output.hic_r1} -O {output.hic_r2} -h {output.html_hic} -j {output.json_hic}
        """

# GENOME FILTERING
if omit == "True":
    '''
    If there is one known species in the metagenomic samples that is in high abundance
    it can be filtered out using BWA alignment and filtering
    '''
    ## 01 - Read Processing, bwa mem & samtools
    rule bwa_samtools_1:
        input:
            ## Reference Genomes
            ref = reference + "/GCA_004519485.1_ASM451948v1_genomic.fna",
            ## Shotgun
            sho_r1 = expand(output_dir + "/fastp/{samp}/sho/{samp}_sho_R1.fastq", samp=SAMPLE),
            sho_r2 = expand(output_dir + "/fastp/{samp}/sho/{samp}_sho_R2.fastq", samp=SAMPLE),
            ## Hi-C
            hic_r1 = expand(output_dir + "/fastp/{samp}/hic/{samp}_hic_R1.fastq", samp=SAMPLE),
            hic_r2 = expand(output_dir + "/fastp/{samp}/hic/{samp}_hic_R2.fastq", samp=SAMPLE),
        output:
            ## Shotgun
            sho_r1 = output_dir + "/samtools_1/{samp}/sho/unmapped_{samp}_sho_R1.fastq",
            sho_r2 = output_dir + "/samtools_1/{samp}/sho/unmapped_{samp}_sho_R2.fastq",
            ## Hi-C
            hic_r1 = output_dir + "/samtools_1/{samp}/hic/unmapped_{samp}_hic_R1.fastq",
            hic_r2 = output_dir + "/samtools_1/{samp}/hic/unmapped_{samp}_hic_R2.fastq"
        conda:
            "envs/env_1.yaml"
        threads:
            threads
        shell:
            """
            ## Generate Index files for reference sequence genome
            bwa index {input.ref}

            ## Align Data: Shotgun
            /usr/bin/time -v bwa mem -t {threads} {input.ref} {input.sho_r1} {input.sho_r2} | \
            /usr/bin/time -v samtools view -@ {threads} -S -b | \
            /usr/bin/time -v samtools view -@ {threads} -b -f 4 | \
            /usr/bin/time -v samtools sort -@ {threads} -n | \
            /usr/bin/time -v samtools fastq -@ {threads} -1 {output.sho_r1} -2 {output.sho_r2} -0 /dev/null -s /dev/null -n
            
            ## Align Data: Hi-C
            /usr/bin/time -v bwa mem -5SP -t {threads} {input.ref} {input.hic_r1} {input.hic_r2} | \
            /usr/bin/time -v samtools view -@ {threads} -S -b  \
            /usr/bin/time -v samtools view -@ {threads} -b -f 4 | \
            /usr/bin/time -v samtools sort -@ {threads} -n | \
            /usr/bin/time -v samtools fastq -@ {threads} -1 {output.hic_r1} -2 {output.hic_r2} -0 /dev/null -s /dev/null -n
            """
        

else:

    rule make_dirs:
        output:
            output_dir + "/burn.txt"
        params:
            # main dir
            main = output_dir + "/samtools_1",

            # sho
            samp = expand(output_dir + "/samtools_1/{samp}", samp = SAMPLE),

            # sho
            samp_sho = expand(output_dir + "/samtools_1/{samp}/sho", samp = SAMPLE),

            # hic
            samp_hic = expand(output_dir + "/samtools_1/{samp}/hic", samp = SAMPLE)

        shell:
            '''
            touch {output}
            mkdir {params.main}
            mkdir {params.samp}
            mkdir {params.samp_sho}
            mkdir {params.samp_hic}
            '''
        
    rule movefiles:
        input:
            burn = output_dir + "/burn.txt",
            ## Shotgun
            sho_r1 = output_dir + "/fastp/{samp}/sho/{samp}_sho_R1.fastq",
            sho_r2 = output_dir + "/fastp/{samp}/sho/{samp}_sho_R2.fastq",
            ## Hi-C
            hic_r1 = output_dir + "/fastp/{samp}/hic/{samp}_hic_R1.fastq",
            hic_r2 = output_dir + "/fastp/{samp}/hic/{samp}_hic_R2.fastq"
        output:
            # sho 
            sho_r1 = output_dir + "/samtools_1/{samp}/sho/unmapped_{samp}_sho_R1.fastq",
            sho_r2 = output_dir + "/samtools_1/{samp}/sho/unmapped_{samp}_sho_R2.fastq",
            ## Hi-C
            hic_r1 = output_dir + "/samtools_1/{samp}/hic/unmapped_{samp}_hic_R1.fastq",
            hic_r2 = output_dir + "/samtools_1/{samp}/hic/unmapped_{samp}_hic_R2.fastq"
        threads:
            threads
        shell:
            '''
            mv -T {input.sho_r1} {output.sho_r1}
            mv -T {input.sho_r2} {output.sho_r2}
            mv -T {input.hic_r1} {output.hic_r1}
            mv -T {input.hic_r2} {output.hic_r2}
            '''



# 02 - Assembly, metaspades
rule metaspades:
    """
    Function:   Produce metagenome assemblies of all shotgun populations
    Parameters:
        > --meta        Indicates that spades should be ran as metaspades
        > -k            Specify kmer lengths to perform assembly with
        > -1            Input Read 1
        > -2            Input Read 2
        > --threads     Number of threads to use
        > -o            Output Directory
    Input:      Nannochloropsis-free reads (.fastq)
    Output:     Shotgun assembly (format = scaffolds.fasta)
    """
    input:
        r1 = expand(output_dir + "/samtools_1/{samp}/sho/unmapped_{samp}_sho_R1.fastq", samp=SAMPLE),
        r2 = expand(output_dir + "/samtools_1/{samp}/sho/unmapped_{samp}_sho_R2.fastq", samp=SAMPLE),
    output:
        vDir  = directory(output_dir + "/metaspades_out/{samp}_sho"),
        scaf=output_dir + "/metaspades_out/{samp}_sho/scaffolds.fasta",
        r1=output_dir + "/metaspades_out/{samp}_sho/corrected/R1.fastq", 
        r2=output_dir + "/metaspades_out/{samp}_sho/corrected/R2.fastq"

    conda:
        "envs/env_1.yaml"
    threads:
        threads
    shell:
        """
        /usr/bin/time -v spades.py --meta \
        -k 21,33,55 \
        --threads {threads} \
        -1 {input.r1} \
        -2 {input.r2} \
        -o {output.vDir}
        """

## 02 - Assembly, QUAST
rule quast:
    input:
        expand(output_dir + "/metaspades_out/{samp}_sho/scaffolds.fasta", samp=SAMPLE),
    output:
        directory(output_dir + "/quast_out/{samp}_sho")
    conda:
        "envs/env_2.yaml"
    threads:
        threads
    shell:
        """
        /usr/bin/time -v metaquast \
        --threads {threads} \
        -o {output} \
        {input}
        """

## 03 - Binning, MaxBin2
rule maxbin2:
    """
    Function:   Take MAG and group genomes using _.
    Resources:
        > https://sourceforge.net/projects/maxbin2/files/
    """
    input:
        contigs = expand(output_dir + "/metaspades_out/{samp}_sho/scaffolds.fasta", samp=SAMPLE),
        r1 = expand(output_dir + "/metaspades_out/{samp}_sho/corrected/R1.fastq", samp=SAMPLE),
        r2 = expand(output_dir + "/metaspades_out/{samp}_sho/corrected/R2.fastq", samp=SAMPLE),
    output:
        vOut = directory(output_dir + "/maxbin2_out/{samp}_sho/"), ## Defined only
        scaf = output_dir + "/maxbin2_out/{samp}_sho/scaffolds.fasta"
    conda:
        "envs/env_1.yaml"
    threads:
        threads
    shell:
        """
        /usr/bin/time -v run_MaxBin.pl \
        -contig {input.contigs} \
        -out {output.vOut}bin \
        -reads {input.r1} \
        -reads2 {input.r2} \
        -thread {threads}
        """

# 03 - Binning, Generate alignment (used by MetaBat2 and CONCOCT)
rule align_assembly_sho:
    """
    Function:
        1. Copy metagenome assemblies into a new directory, where they are indexed.
        2. Generate sho alignment BAM files for use with MetaBat2 and CONCOCT
        3. The MetaBat2 "depth" file is created
    Resources:
        > https://onestopdataanalysis.com/binning/
    """
    input:
        vRef     = expand(output_dir + "/metaspades_out/{samp}_sho/scaffolds.fasta", samp=SAMPLE),
        vR1_sho = expand(output_dir + "/metaspades_out/{samp}_sho/corrected/R1.fastq", samp=SAMPLE),
        vR2_sho = expand(output_dir + "/metaspades_out/{samp}_sho/corrected/R2.fastq", samp=SAMPLE),
    output:
        vDir_index     = directory(output_dir + "/mapping_out/index/{samp}_sho/"),
        vIndex         =           output_dir + "/mapping_out/index/{samp}_sho/scaffolds.fasta",
        vAlignment_sho =           output_dir + "/mapping_out/alignment/{samp}_sho.bam",
        vDepth         =           output_dir + "/mapping_out/depth/{samp}_depth.txt"
    conda:
        "envs/env_1.yaml"
    threads:
        threads
    shell:
        """
        ## 1.
        mkdir -p {output.vDir_index}
        cp {input.vRef} {output.vIndex}
        /usr/bin/time -v bwa index {output.vIndex}

        ## 2.
        /usr/bin/time -v bwa mem -t {threads} {output.vIndex} {input.vR1_sho} {input.vR2_sho} \
            | /usr/bin/time -v samtools sort -@ {threads} -o {output.vAlignment_sho}
        /usr/bin/time -v samtools index -@ {threads} {output.vAlignment_sho}

        ## 3.
        /usr/bin/time -v jgi_summarize_bam_contig_depths --outputDepth {output.vDepth} {output.vAlignment_sho}
        """

# 03 - Binning, Generate alignment (used by Bin3c)
rule align_assembly_hic:
    """
    Function: Generate hic alignment BAM files for use with Bin3c
    """
    input:
        vRef     = expand(output_dir + "/metaspades_out/{samp}_sho/scaffolds.fasta", samp=SAMPLE),
        vR1_hic  = expand(output_dir      + "/samtools_1/{samp}/hic/unmapped_{samp}_hic_R1.fastq", samp=SAMPLE),
        vR2_hic  = expand(output_dir      + "/samtools_1/{samp}/hic/unmapped_{samp}_hic_R2.fastq", samp=SAMPLE),
        vIndex   = expand(output_dir    + "/mapping_out/index/{samp}_sho/scaffolds.fasta", samp=SAMPLE),
    output:
        vAlignment_hic = output_dir + "/mapping_out/alignment/{samp}_hic.bam",
    conda:
        "envs/env_1.yaml"
    threads:
        threads
    shell:
        """
        /usr/bin/time -v bwa mem -t {threads} -5SP {input.vIndex} {input.vR1_hic} {input.vR2_hic} \
            | /usr/bin/time -v samtools view -@ {threads} -Sb - \
            | /usr/bin/time -v samtools sort -@ {threads} -o {output.vAlignment_hic} -n -
        """
            
## 03 - Binning, MetaBat2
rule metabat2:
    """
    Function: 
    Notes:
        > Double brackets ({{}}) used to escape bracket character in snakemake.
    """
    input:
        vAssembly  = expand(output_dir + "/metaspades_out/{samp}_sho/scaffolds.fasta", samp=SAMPLE),
        vDepth     = expand(output_dir + "/mapping_out/depth/{samp}_depth.txt", samp=SAMPLE),
    output:
        directory(output_dir + "/metabat2_out/{samp}_sho/"),
    conda:
        "envs/env_1.yaml"
    threads:
        threads
    shell:
        """
        /usr/bin/time -v metabat2 -t {threads} -i {input.vAssembly} -a {input.vDepth} -o {output}

        find {output} -type f -name '.*' -execdir sh -c 'mv -i "$0" "./${{0#./.}}"' {{}} \;
        """

## 03 - Binning, CONCOCT
rule concoct:
    """
    Resources:
        > https://onestopdataanalysis.com/binning/
    """
    input:
        vScaffold  = expand(output_dir + "/metaspades_out/{samp}_sho/scaffolds.fasta", samp=SAMPLE),
        vAlignment = expand(output_dir + "/mapping_out/alignment/{samp}_sho.bam", samp=SAMPLE),
    output:
        vDir  = directory(output_dir + "/concoct_out/{samp}_sho/"),
        vBins = directory(output_dir + "/concoct_out/{samp}_sho/fasta_bins"), ## Defined only
    conda:
        "envs/env_3_concoct.yaml"
    threads:
        threads
    shell:
        """
        mkdir {output.vDir}/fasta_bins

        /usr/bin/time -v cut_up_fasta.py {input.vScaffold} -c 10000 -o 0 --merge_last -b {output.vDir}/contigs_10K.bed \
            > {output.vDir}/contigs_10K.fa

        /usr/bin/time -v concoct_coverage_table.py {output.vDir}/contigs_10K.bed {input.vAlignment} \
            > {output.vDir}/coverage_table.tsv

        /usr/bin/time -v concoct --threads {threads} --composition_file {output.vDir}/contigs_10K.fa --coverage_file {output.vDir}/coverage_table.tsv -b {output.vDir}

        /usr/bin/time -v merge_cutup_clustering.py {output.vDir}/clustering_gt1000.csv \
            > {output.vDir}/clustering_merged.csv

        /usr/bin/time -v extract_fasta_bins.py {input.vScaffold} {output.vDir}/clustering_merged.csv --output_path {output.vDir}/fasta_bins
        """

## 03 - Binning, Bin3c
rule bin3c:
    """
    Notes
        > Walk-through: https://github.com/cerebis/bin3C
        > Environment Info: https://github.com/cerebis/proxigenomics_toolkit/blob/master/Pipfile
        > Environment Info: See github bin3C pipfile
    """
    input:
        vScaffold      = expand(output_dir + "/metaspades_out/{samp}_sho/scaffolds.fasta", samp=SAMPLE),
        vAlignment_hic = expand(output_dir + "/mapping_out/alignment/{samp}_hic.bam", samp=SAMPLE),
    output:
        vDir = directory(output_dir + "/bin3c_out/{samp}"),
        vDef = directory(output_dir + "/bin3c_out/{samp}/fasta/"), # Define only
    conda:
        "envs/env_4_bin3c.yaml"
    threads:
        threads
    shell:
        """
        ## Install bin3c
        git clone --recursive https://github.com/cerebis/bin3C

        ## Run Bin3c, 
        /usr/bin/time -v python bin3C/bin3C.py mkmap -e MluCI -e Sau3AI -v {input.vScaffold} {input.vAlignment_hic} map
        
        rm -r {output.vDir}
        /usr/bin/time -v python bin3C/bin3C.py cluster -v map/contact_map.p.gz {output.vDir}

        rm -r bin3C map
        """

rule dastool_tables:
    """
    Function: For each binning output, create a .tsv file associating bins with scaffolds
    """
    input:
        vIn_Maxbin2  = expand(output_dir  + "/maxbin2_out/{samp}_sho/", samp=SAMPLE),
        vIn_Metabat2 = expand(output_dir + "/metabat2_out/{samp}_sho/", samp=SAMPLE),
        vIn_Concoct  = expand(output_dir  + "/concoct_out/{samp}_sho/fasta_bins", samp=SAMPLE),
        vIn_Bin3c    = expand(output_dir    + "/bin3c_out/{samp}/fasta/", samp=SAMPLE),
    output:
        vOut_Maxbin2  = output_dir + "/dastool_out/table/{samp}/maxbin.scaffolds2bin.tsv",
        vOut_Metabat2 = output_dir + "/dastool_out/table/{samp}/metabat.scaffolds2bin.tsv",
        vOut_Concoct  = output_dir + "/dastool_out/table/{samp}/concoct.scaffolds2bin.tsv",
        vOut_Bin3c    = output_dir + "/dastool_out/table/{samp}/bin3c.scaffolds2bin.tsv",
    conda:
        "envs/env_5_dastool.yaml"
    threads:
        threads
    shell:
        """
        /usr/bin/time -v Fasta_to_Scaffolds2Bin.sh -i {input.vIn_Maxbin2}  -e fasta > {output.vOut_Maxbin2}
        /usr/bin/time -v Fasta_to_Scaffolds2Bin.sh -i {input.vIn_Metabat2} -e fa    > {output.vOut_Metabat2}
        /usr/bin/time -v Fasta_to_Scaffolds2Bin.sh -i {input.vIn_Concoct}  -e fa    > {output.vOut_Concoct}
        /usr/bin/time -v Fasta_to_Scaffolds2Bin.sh -i {input.vIn_Bin3c}    -e fna | ./scripts/sed_command.sh > {output.vOut_Bin3c} 
        """



rule dastool_predict:
    """
    Notes
        > Search engine "diamond" used for its easy conda installation
    """
    input:
        vScaffold    = expand(output_dir + "/metaspades_out/{samp}_sho/scaffolds.fasta", samp=SAMPLE),
        vIn_Maxbin2  = expand(output_dir + "/dastool_out/table/{samp}/maxbin.scaffolds2bin.tsv", samp=SAMPLE),
        vIn_Metabat2 = expand(output_dir + "/dastool_out/table/{samp}/metabat.scaffolds2bin.tsv", samp=SAMPLE),
        vIn_Concoct  = expand(output_dir + "/dastool_out/table/{samp}/concoct.scaffolds2bin.tsv", samp=SAMPLE),
        vIn_Bin3c    = expand(output_dir + "/dastool_out/table/{samp}/bin3c.scaffolds2bin.tsv", samp=SAMPLE),
    output:
        vOut_sho_hic = directory(output_dir + "/dastool_out/prediction_sho_hic/{samp}/"),
        vOut         = directory(output_dir + "/dastool_out/prediction_sho_hic/{samp}/_DASTool_bins/"), ## Define only
    conda:
        "envs/env_5_dastool.yaml"
    threads:
        threads
    shell:
        """
        ## Run WITH Bin3c data
        /usr/bin/time -v DAS_Tool \
            -i {input.vIn_Maxbin2},{input.vIn_Metabat2},{input.vIn_Concoct},{input.vIn_Bin3c} \
            -l maxbin2,metabat2,concoct,bin3c \
            -c {input.vScaffold} \
            --search_engine diamond \
            --write_bins 1 \
            --threads {threads} \
            -o {output.vOut_sho_hic}
        """

rule gtdbtk:
    """
    Notes:
        > GTDB-Tk requires ~27G of external data that needs to be downloaded and unarchived (performed in resources directory ahead of time):
            $ wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/auxillary_files/gtdbtk_r95_data.tar.gz
            $ tar xvzf gtdbtk_r95_data.tar.gz
    """
    input:
        expand(output_dir + "/dastool_out/prediction_sho_hic/{samp}/_DASTool_bins/", samp=SAMPLE),
    output:
        directory(output_dir + "/gtdbtk_out/{samp}/"),
    params:
        gtdbtk_data
    conda:
        "envs/env_6_gtdbtk.yaml"
    threads:
        threads
    shell:
        """
        GTDBTK_DATA_PATH={params} /usr/bin/time -v gtdbtk classify_wf \
            --cpus {threads} \
            --extension fa \
            --genome_dir {input} \
            --out_dir {output}
        """

rule checkm:
    input:
        vBins      = expand(output_dir + "/dastool_out/prediction_sho_hic/{samp}/_DASTool_bins/", samp=SAMPLE),
        vAlignment = expand(output_dir + "/mapping_out/alignment/{samp}_sho.bam", samp=SAMPLE),
    output:
        vDir = directory(output_dir + "/checkm_out/{samp}/"),
        vCov = output_dir + "/checkm_out/{samp}/coverage.tsv",
        vPro = output_dir + "/checkm_out/{samp}/profile.tsv",
        vSum = output_dir + "/checkm_out/{samp}/summary.tsv",
    conda:
        "envs/env_7_checkm.yaml"
    threads:
        threads
    shell:
        """
        ## Generate Coverage Table
        /usr/bin/time -v checkm coverage \
            --threads {threads} \
            --extension fa \
            {input.vBins} \
            {output.vCov} \
            {input.vAlignment}
        
        ## Generate Mapping % Table
        /usr/bin/time -v checkm profile \
            --tab_table \
            --file {output.vPro} \
            {output.vCov}

        ## Perform Taxonomic Classification
        /usr/bin/time -v checkm lineage_wf \
            --threads {threads} \
            --extension fa \
            --file {output.vSum} \
            {input.vBins} \
            {output.vDir}
        """

