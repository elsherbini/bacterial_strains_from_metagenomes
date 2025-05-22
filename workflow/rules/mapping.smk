# Snakemake rules for mapping metagenomic reads to the genome database

import os
import pandas as pd

# --- Rule to Map Reads to Database and Generate Sorted BAM files ---
rule map_reads_to_db:
    input:
        r1 = lambda wildcards: samples.loc[wildcards.sample, "fastq_1"],
        r2 = lambda wildcards: samples.loc[wildcards.sample, "fastq_2"],
        bowtie2_index_flag = os.path.join(BOWTIE2_INDEX_DIR, "concatenated_genomes.index_SUCCESS")
    output:
        sam = os.path.join(MAPPED_READS_DIR, "{sample}", "{sample}.sam"),
    log:
        os.path.join(RESULTS_DIR, "logs", "mapping", "{sample}.log")
    params:
        bowtie2_index = os.path.join(BOWTIE2_INDEX_DIR, "concatenated_genomes"),
        # Mapping quality filter for Samtools, if needed
    threads: 
        8
    conda:
        "../envs/bowtie2.yaml"  # Assuming bowtie2.yaml includes samtools, otherwise use a combined env
    shell:
        """
        bowtie2 \
            -x {params.bowtie2_index} \
            -1 {input.r1} \
            -2 {input.r2} \
            --threads {threads} \
            2> {log} > {output.sam}
        """
