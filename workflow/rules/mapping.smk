# Snakemake rules for mapping metagenomic reads to the genome database

import os
import pandas as pd

# --- Rule to Map Reads to Database and Generate Sorted BAM files ---
rule map_reads_to_db:
    input:
        r1 = lambda wildcards: SAMPLES_DF.loc[SAMPLES_DF['sample_id'] == wildcards.sample, "fastq_1"].iloc[0],
        r2 = lambda wildcards: SAMPLES_DF.loc[SAMPLES_DF['sample_id'] == wildcards.sample, "fastq_2"].iloc[0],
        bowtie2_index_flag = os.path.join(BOWTIE2_INDEX_DIR, "concatenated_genomes.index_SUCCESS")
    output:
        sam = os.path.join(MAPPED_READS_DIR, "{sample}.sam"),
    log:
        os.path.join(RESULTS_DIR, "logs", "mapping", "{sample}.log")
    params:
        bowtie2_index = os.path.join(BOWTIE2_INDEX_DIR, "concatenated_genomes"),
        # Mapping quality filter for Samtools, if needed
    threads: 
        8
    resources:
        cpus_per_task = 8, 
        runtime = "4h",
        mem_mb = 16000,
        partition = "short"
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
