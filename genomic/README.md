# Genomic

This file contains each scripts used for the genomic exctraction and analysis of the data.

The "Commandlinegrepgenomic" script is used to transform each breseq output into 3 lists of genes and locus tags. Each instance of either gene or locus tag represents a mutation on the gene.
Each list correspond to a type of mutation :
- Nonsynonymous
- Insertion/deletion
- Nonsense


The "transformgrepedtotablecounts" converts the obtained lists in table of counts similar to the transcriptomics table.

The "sumfusecounts" fuse nonsynonymous and Insertion/deletion mutation into one table of count.
