# ncldv_markersearch
Script for annotating and merging phylogenetic marker proteins in nucleocytoplasmic large DNA virus genomes. 
This script will search a file of proteins against a set of 10 curated Hidden Markov Models for protein families prevalent in NCLDV. 


# USAGE: python ncldv_markersearch.py <directory of protein .faa files> <proximity of genes to merge>

The proximity of genes to merge is the number of genes up- and down-stream of the initial hit that the program will search for additional hits to merge. 
