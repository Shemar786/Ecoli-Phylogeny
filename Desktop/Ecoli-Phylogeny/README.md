# E. coli Phylogeny Project

This repository contains my work on building a **phylogenetic tree of *Escherichia coli* genomes** as part of my research project. The goal was to analyze hundreds of *E. coli* genomes, align them, and reconstruct a tree showing their evolutionary relationships.

---

## What I Did

1. **Collected genomes**  
   - Started with ~625 *E. coli* genomes in FASTA format.  
   - Cleaned and organized them into a consistent format for analysis.  

2. **Built a pipeline**  
   - Wrote `build_tree_parsnp.py`, a Python script that:  
     - Filters genomes against a reference.  
     - Runs [Parsnp](https://github.com/marbl/parsnp) inside Docker to generate a core-genome alignment.  
     - Uses [FastTree](http://www.microbesonline.org/fasttree/) to reconstruct the phylogenetic tree from SNPs.  

3. **Generated results**  
   - After filtering, **94 genomes** were included in the final analysis.  
   - Parsnp outputs SNP alignments (`parsnp.snps.mblocks`).  
   - FastTree generated a Newick tree (`fasttree_output.tree`).  

4. **Visualized the tree**  
   - The `.tree` file can be visualized in tools like [iTOL](https://itol.embl.de/) or FigTree.  
   - Produced a circular phylogenetic tree showing how the 94 genomes cluster.

---

## Repository Structure

