# ncldv_markersearch
Script for annotating and merging phylogenetic marker proteins in nucleocytoplasmic large DNA virus genomes. 

This script will search a file of proteins against a set of 10 curated Hidden Markov Models for protein families prevalent in NCLDV. Proteins that have hits to the same HMM and form non-overlapping alignments with that HMM will be merged together and output as a single aa sequence assuming they are encoded close together on the genome. 

The program assumes protein IDs are provided in the Prodigal format (i.e., contigname_1, contigname_2, etc). 
The program requires HMMER3 and the SeqIO package of Biopython. 


### USAGE: python ncldv_markersearch.py <directory of protein .faa files> <proximity of genes to merge>

The proximity of genes to merge is the number of genes up- and down-stream of the initial hit that the program will search for additional hits to merge. Hits outside this range will be considered independent hits. 

### Output
ncldv_markersearch.py provides several output files:

full_output.txt         This is the main tab-delimited output file that provides the annotation results. 

ncldv_markersearch.faa  This is the protein file with all merged and unmerged proteins with best hits to the HMMs. Proteins are re-named to accommodate potential merged proteins. 

raw_output.txt          This is the parsed raw HMMER3 output. It can be used as a reference for debugging. 

cogs.txt                This is a cogs-formatted file that can be used as input for an ETE3 species tree workflow 
(http://etetoolkit.org/documentation/ete-build/).




