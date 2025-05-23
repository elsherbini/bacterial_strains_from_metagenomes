"""
# Rules for inStrain analysis of mapped metagenomic samples
"""
rule instrain_profile:
    input:
        sam = "results/mapped_reads/{sample}.sam",
        fasta = "results/genome_database/concatenated_genomes.fasta",
        genes = "results/genome_database/database_genes.fna",
        stb = "results/genome_database/scaffold_to_bin.tsv"
    output:
        directory("results/instrain/profiles/{sample}")
    params:
        output_prefix = "results/instrain/profiles/{sample}/{sample}",
        min_read_ani = config.get("instrain", {}).get("min_read_ani", 0.95),
        min_mapq = config.get("instrain", {}).get("min_mapq", 2),
        extra = config.get("instrain", {}).get("profile_extra", "--database_mode")
    log:
        "logs/instrain/profile/{sample}.log"
    conda:
        "../envs/instrain.yaml"
    threads: 
        config.get("instrain", {}).get("threads", 8)
    resources:
        cpus_per_task = lambda wildcards, threads: threads,
        runtime = "8h",
        mem_mb = config.get("instrain", {}).get("memory_mb", 32000),
        partition = "short"
    shell:
        """
        inStrain profile {input.sam} {input.fasta} \
            -o {params.output_prefix} \
            -p {threads} \
            -g {input.genes} \
            -s {input.stb} \
            --min_read_ani {params.min_read_ani} \
            --min_mapq {params.min_mapq} \
            {params.extra} > {log} 2>&1
        """

rule instrain_compare:
    input:
        profiles = expand("results/instrain/profiles/{sample}", sample=SAMPLES)
    output:
        directory("results/instrain/comparison")
    params:
        output_prefix = "results/instrain/comparison/comparison",
        profiles_txt = "results/instrain/profile_list.txt",
        min_cov = config.get("instrain", {}).get("min_coverage", 5),
        extra = config.get("instrain", {}).get("compare_extra", "")
    log:
        "logs/instrain/compare.log"
    conda:
        "../envs/instrain.yaml"
    threads: 
        config.get("instrain", {}).get("threads", 8)
    resources:
        cpus_per_task = lambda wildcards, threads: threads,
        runtime = "6h",
        mem_mb = config.get("instrain", {}).get("memory_mb", 32000),
        partition = "short"
    shell:
        """
        # Create a list of profile paths
        echo {input.profiles} | tr ' ' '\\n' > {params.profiles_txt}
        
        # Run inStrain compare
        inStrain compare \
            -i {params.profiles_txt} \
            -o {params.output_prefix} \
            -p {threads} \
            --min_cov {params.min_cov} \
            {params.extra} > {log} 2>&1
        """
