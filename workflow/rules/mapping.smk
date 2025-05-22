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
        bam = os.path.join(MAPPING_RESULTS_DIR, "{sample}", "{sample}.sorted.bam"),
        bai = os.path.join(MAPPING_RESULTS_DIR, "{sample}", "{sample}.sorted.bam.bai")
    log:
        os.path.join(RESULTS_DIR, "logs", "mapping", "{sample}.log")
    params:
        bowtie2_index = os.path.join(BOWTIE2_INDEX_DIR, "concatenated_genomes"),
        bowtie2_options = config["params"]["bowtie2"]["options"],
        # Mapping quality filter for Samtools, if needed
        mapq_filter = config["params"]["samtools"]["mapq_filter"]
    threads: 
        config["params"]["bowtie2"]["threads"]
    conda:
        "../envs/bowtie2.yaml"  # Assuming bowtie2.yaml includes samtools, otherwise use a combined env
    shell:
        """
        # Map reads using Bowtie2 and pipe to samtools for BAM conversion and sorting
        bowtie2 \
            -x {params.bowtie2_index} \
            -1 {input.r1} \
            -2 {input.r2} \
            {params.bowtie2_options} \
            --threads {threads} \
            2> {log} | \
        samtools view -b -q {params.mapq_filter} - | \
        samtools sort -o {output.bam} -@ {threads}

        # Index the BAM file
        samtools index -@ {threads} {output.bam}
        """