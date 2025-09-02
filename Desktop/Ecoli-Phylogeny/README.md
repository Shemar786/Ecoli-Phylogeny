# E. coli Phylogeny with Parsnp + FastTree

This repository contains a pipeline and results for building a core-genome phylogenetic tree of *E. coli* genomes using Parsnp and FastTree.

## Files
- `build_tree_parsnp.py` – Python script to run Parsnp on a folder of genomes and generate SNP alignments.
- `results/fasttree_output.tree` – Phylogenetic tree (Newick format) of 94 genomes.
- `results/log/` – Log files from Parsnp/RAxML runs.

## Usage
1. Install Docker, Parsnp, and FastTree.
2. Place genome `.fasta` files in a folder.
3. Run:
   ```bash
   python3 build_tree_parsnp.py

