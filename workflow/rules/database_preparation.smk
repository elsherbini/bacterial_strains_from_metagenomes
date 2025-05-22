# Snakemake rules for genome database preparation

import os
import pandas as pd

# --- Rule to Create Scaffold-to-Bin File ---
rule create_scaffold_to_bin:
    input:
        genome_paths = GENOME_PATHS, # List of paths to individual genome FASTA files
        taxonomy_file = TAXONOMY_FILE # This can be an empty file or a file with actual taxonomy info
    output:
        stb_file = os.path.join(GENOME_DB_RESULTS_DIR, "scaffold_to_bin.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs", "create_scaffold_to_bin.log")
    threads: 1
    conda:
        "../envs/instrain.yaml" # General python environment
    script:
        "../scripts/generate_stb.py"

# --- Rule to Run Prodigal on Each Genome and Concatenate Genes ---
# This rule will run prodigal on each genome individually.
rule run_prodigal_on_genome:
    input:
        genome_fasta = lambda wildcards: GENOME_PATHS[GENOME_BASENAMES.index(wildcards.genome_basename)]
    output:
        genes_fna = os.path.join(PRODIGAL_RESULTS_DIR, "{genome_basename}.genes.fna"),
        genes_faa = os.path.join(PRODIGAL_RESULTS_DIR, "{genome_basename}.genes.faa"),
        # proteins_gff = os.path.join(PRODIGAL_RESULTS_DIR, "{genome_basename}.genes.gff") # Optional output
    log:
        os.path.join(RESULTS_DIR, "logs", "prodigal", "{genome_basename}.log")
    params:
        mode = config["params"]["prodigal"]["mode"] # "single" or "meta"
    threads: 1 # Prodigal is usually single-threaded for one genome
    conda:
        "../envs/prodigal.yaml"
    shell:
        """
        prodigal -i {input.genome_fasta} \
                 -o {output.genes_fna}.tmp_gff_output_for_prodigal_change_this_to_a_real_gff_output_file \
                 -a {output.genes_faa} \
                 -d {output.genes_fna} \
                 -p {params.mode} \
                 -q  # Quiet mode
        # Prodigal's -o option is for GFF, -d for nuc, -a for prot.
        # The example in DESIGN.md suggests concatenating .fna and .faa. We do this in a later rule.
        # This rule creates per-genome gene files.
        # rm {output.genes_fna}.tmp_gff_output_for_prodigal_change_this_to_a_real_gff_output_file # Clean up temp GFF if not needed
        """

rule aggregate_prodigal_outputs:
    input:
        fnas = expand(os.path.join(PRODIGAL_RESULTS_DIR, "{gb}.genes.fna"), gb=GENOME_BASENAMES),
        faas = expand(os.path.join(PRODIGAL_RESULTS_DIR, "{gb}.genes.faa"), gb=GENOME_BASENAMES)
    output:
        all_genes_fna = os.path.join(GENOME_DB_RESULTS_DIR, "database_genes.fna"),
        all_genes_faa = os.path.join(GENOME_DB_RESULTS_DIR, "database_genes.faa")
    log:
        os.path.join(RESULTS_DIR, "logs", "aggregate_prodigal.log")
    threads: 1
    shell:
        """
        cat {input.fnas} > {output.all_genes_fna}
        cat {input.faas} > {output.all_genes_faa}
        """

# --- Rule to Concatenate Genomes and Build Bowtie2 Index ---
rule concatenate_genomes:
    input:
        GENOME_PATHS # List of paths to individual genome FASTA files
    output:
        concatenated_fasta = os.path.join(GENOME_DB_RESULTS_DIR, "concatenated_genomes.fasta")
    log:
        os.path.join(RESULTS_DIR, "logs", "concatenate_genomes.log")
    threads: 1
    shell:
        "cat {input} > {output} 2> {log}"

rule build_bowtie2_index:
    input:
        fasta = rules.concatenate_genomes.output.concatenated_fasta
    output:
        # Snakemake tracks index by prefix, but good to have a flag file or one specific output
        # BOWTIE2_INDEX_DIR is the directory, concatenated_genomes is the prefix
        #bowtie2_index_files = multiext(os.path.join(BOWTIE2_INDEX_DIR, "concatenated_genomes"), ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"),
        # Or for large indexes: multiext(os.path.join(BOWTIE2_INDEX_DIR, "concatenated_genomes"), ".1.bt2l", ".2.bt2l", ".3.bt2l", ".4.bt2l", ".rev.1.bt2l", ".rev.2.bt2l")
        # Using a touch file to mark completion, as Snakemake sometimes struggles with multi-file outputs like indexes
        completion_flag = touch(os.path.join(BOWTIE2_INDEX_DIR, "concatenated_genomes.index_SUCCESS"))
    log:
        os.path.join(RESULTS_DIR, "logs", "bowtie2_build.log")
    params:
        prefix = os.path.join(BOWTIE2_INDEX_DIR, "concatenated_genomes"),
        threads = config["params"]["bowtie2_build"]["threads"],
        # Consider adding --large-index if many/large genomes
        # This could be a function: lambda wildcards, input: "--large-index" if os.path.getsize(input.fasta) > 4*1024*1024*1024 else ""
        large_index_opt = "--large-index" # Default to large-index for robustness, can be configured
    threads:
        config["params"]["bowtie2_build"]["threads"]
    conda:
        "../envs/bowtie2.yaml"
    shell:
        """
        bowtie2-build --threads {threads} {params.large_index_opt} {input.fasta} {params.prefix} > {log} 2>&1
        """ 