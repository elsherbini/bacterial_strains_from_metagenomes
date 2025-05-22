#!/usr/bin/env python3

import argparse
import os
import sys

def extract_scaffold_and_genome_id(fasta_file_path):
    """
    Extracts scaffold headers from a FASTA file and associates them with a genome ID.
    The genome ID is derived from the basename of the FASTA file.
    """
    scaffolds = []
    try:
        genome_id = os.path.splitext(os.path.basename(fasta_file_path))[0]
        with open(fasta_file_path, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    scaffold_name = line[1:].strip().split()[0] # Get the first part of the header
                    scaffolds.append((scaffold_name, genome_id))
    except FileNotFoundError:
        print(f"Error: File not found {fasta_file_path}", file=sys.stderr)
        return []
    except Exception as e:
        print(f"Error processing file {fasta_file_path}: {e}", file=sys.stderr)
        return []
    return scaffolds

def main():
    parser = argparse.ArgumentParser(description="Generate a scaffold-to-bin.tsv file from genome FASTA files.")
    parser.add_argument(
        "--genome-paths",
        nargs='+',
        required=False, # Changed to False to allow snakemake to provide it
        help="List of paths to input genome FASTA files. (Provided by Snakemake if run in a rule)"
    )
    parser.add_argument(
        "--output-file",
        required=False, # Changed to False to allow snakemake to provide it
        help="Path to the output scaffold_to_bin.tsv file. (Provided by Snakemake if run in a rule)"
    )
    parser.add_argument(
        "--taxonomy-file",
        help="Path to a taxonomy file (optional, currently not used for genome ID assignment but could be in future). Example format: genome_filename_stem<tab>taxonomy_string"
    )

    args = parser.parse_args()

    # Snakemake integration:
    if "snakemake" in globals():
        genome_files = snakemake.input.get("genome_paths", [])
        if not genome_files:
            # Compatibility with older Snakemake versions or direct input list
            genome_files = snakemake.input if isinstance(snakemake.input, list) else []
        output_file = str(snakemake.output.stb_file)
        # taxonomy_file = snakemake.input.get("taxonomy_file") # If needed
    else:
        # Command-line execution
        if not args.genome_paths:
            parser.error("argument --genome-paths is required when not run via Snakemake.")
        if not args.output_file:
            parser.error("argument --output-file is required when not run via Snakemake.")
        genome_files = args.genome_paths
        output_file = args.output_file
        # taxonomy_file = args.taxonomy_file

    if not genome_files:
        print("Error: No genome paths provided.", file=sys.stderr)
        sys.exit(1)

    all_mappings = []
    for genome_file in genome_files:
        mappings = extract_scaffold_and_genome_id(genome_file)
        all_mappings.extend(mappings)

    try:
        with open(output_file, 'w') as outfile:
            outfile.write("scaffold\tgenome_id\n") # Header
            for scaffold, genome_id in all_mappings:
                outfile.write(f"{scaffold}\t{genome_id}\n")
        print(f"Successfully wrote scaffold-to-bin mapping to {output_file}")
    except Exception as e:
        print(f"Error writing output file {output_file}: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    # Check if running under Snakemake and snakemake object is available
    # This helps avoid argparse errors when Snakemake provides parameters directly.
    if "snakemake" not in globals():
        # If not run by Snakemake, parse arguments. Otherwise, Snakemake handles inputs/outputs.
        pass # Argparse will be handled in main() if not snakemake
    main()
