# Bacterial Strains from Metagenomes Pipeline

A Snakemake workflow for analyzing bacterial strains from shotgun metagenomic data.

## Overview

This pipeline processes shotgun metagenomic reads to identify and analyze bacterial strains using inStrain. The workflow follows Snakemake best practices and is organized into modular components.

## Workflow Steps

1. **Genome Database Preparation**
   - Creates a scaffold-to-bin mapping file
   - Predicts genes using Prodigal
   - Builds a Bowtie2 index for read mapping

2. **Sample Processing and Read Mapping**
   - Maps preprocessed reads to the genome database
   - Produces sorted and indexed BAM files

3. **inStrain Analysis**
   - Profiles each sample's BAM file using `inStrain profile`
   - Compares all profiles using `inStrain compare`

## Input Requirements

### Genome Database
- Location: `bacterial_strains_from_metagenomes/genome_database/`
- Required files:
  - Individual genome FASTA files
  - `taxonomy_annotation.tsv`: Tab-separated file with columns:
    - tax_id: NCBI taxonomy ID
    - taxonomy: Full taxonomic classification
    - path: Path to genome FASTA file

### Sample Data
- Input: Samplesheet (CSV format) with columns:
  - sample_id: Unique identifier for each sample
  - fastq_1: Path to forward reads (R1)
  - fastq_2: Path to reverse reads (R2)
- Requirements:
  - Reads must be preprocessed and host-filtered
  - Paired-end reads in FASTQ format

## Usage

1. **Set up the environment**:
   ```
   conda env create -f workflow/envs/instrain.yaml
   conda activate instrain
   ```

2. **Configure the pipeline**:
   - Edit `config/config.yaml` to specify paths and parameters
   - Edit `config/samples.csv` to list your input samples

3. **Run the pipeline**:
   ```
   snakemake --cores <number_of_cores>
   ```

4. **Run with a cluster** (optional):
   ```
   snakemake --profile <your_cluster_profile>
   ```

## Output Structure

```
results/
├── genome_database/
│   ├── scaffold_to_bin.tsv
│   ├── database_genes.fna
│   ├── database_genes.faa
│   └── concatenated_genomes.fasta
├── mapped_reads/
│   ├── {sample_id}.bam
│   └── {sample_id}.bam.bai
└── instrain/
    ├── profiles/
    │   └── {sample_id}/
    └── comparison/
        └── comparison_*
```

## Configuration

The pipeline is configured through `config/config.yaml`, which includes:
- Paths to input files
- Tool parameters
- Resource requirements
- Output directories

## Dependencies

- Bowtie2
- Samtools
- Prodigal
- inStrain (v1.5.4 or higher)
- Python 3.9+ 