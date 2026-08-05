# SyToGen
#
# This is the main program for these tools. This will pull in motifs from parse_motif,
# and DNA annotated with motif_finder. 
#
# Sytogen will find the optimal mutation of a motif so that the motif is destroyed.
# Preferring wobble bases of codons (using the codon_bias of the DNA contruct as a 
# guide). In regions that don't fall inside a protein (as guided by genbank) the base
# change should prefer a GC/AT-preserving switch. All potential changes should be scored
# and the best scoring mutation should be chosen.
# 
# Output includes a decision matrix and the genbank, gff3/fasta
# 