# ncldv_markersearch
NCLDV_Markersearch: A script for identifying phylogenetic marker genes in NCLDV and generating concatenated alignments.

This script will search a file of proteins against a set of 10 curated Hidden Markov Models for protein families prevalent in Nucleo-Cytoplasmic Large DNA Viruses (NCLDV). Proteins that have hits to the same HMM, form non-overlapping alignments with that HMM, and are located within a certain defined proximity on the same contig will be joined together and output as a single amino sequence. This facilitates the phylogenomic analysis of NCLDV that have split genes.  

The program assumes protein IDs are provided in the Prodigal format (i.e., contigname_1, contigname_2, etc). This is needed for sorting the output and ensuring that proteins are joined properly. The simple prodigal launcher script "prodigal_launcher.py" is provided to assist with this. 

The program requires Python >=3.5, HMMER3, and the SeqIO package of Biopython. Additinally, if you wish to produce concatenated alignment Clustal Omega is also required (in your PATH). 

For questions or comments please contact Frank Aylward at faylward at vt.edu

### MINIMAL USAGE: python ncldv_markersearch.py -i <directory of protein .faa files> -n project_name

### Options

**-p, --proximity**
The proximity of genes to merge is the number of genes up- and down-stream of the initial hit that the program will search for additional hits to merge. Hits outside this range will be considered independent hits. Default is 3. 

**-t, --cpus**
How many CPUs/threads to use for the hmmsearch and clustal steps

**-m, --markerset**
Markers to use. Must be comma-separated list of the following: A32,D5,SFII,mcp,mRNAc,PolB,RNAPL,RNAPS,RNR,VLTF3. The default is A32,mcp,SFII,PolB,VLTF3.

**-r, --redo**
If you have already run script and you want to re-run it with different parameters, you can use the -r flag to avoid re-running HMMER (this saves a bit of time if you're running multiple times)

**-c, --concat**
If this option is specified the script will also output a concatenated alignment of the chose marker genes (FASTA format, one entry per taxa, each marker gene aligned separately with Clustal Omega)

**-a, --allhits**
If this option is specified then all hits marker genes (that are above the predefined bit score thresholds) will be output. This option is not compatible with the -c option. This can be useful if you want to see if certain marker genes are present in multiple copies. 



### Output files
ncldv_markersearch.py provides several output files, all with the prefix designated with the -n option:

*full_output.txt         This is the main tab-delimited output file that provides the annotation results. 

*.faa  This is the protein file with all merged and unmerged proteins with best hits to the HMMs. Proteins are re-named to accommodate potential merged proteins. 

*raw_output.txt          This is the parsed raw HMMER3 output (no marker gene joining). It can be used as a reference for debugging, but you can usually ignore it. 

*.cogs.txt                This is a cogs-formatted file, in the same general format used by the ETE3 toolkit, and it's used as a reference for producing the concatenated alignment. You can also use it with the *.faa output file if you want to make a tree with the ETE3 toolkit instead (http://etetoolkit.org/documentation/ete-build/).. 

*.table.tsv              This is an occurrence table for the number of marker genes that were identified for each file in the input folder. Note that if several query proteins with hits to the same marker gene were joined they only counted once. Also, this table will only show multiple hits if the -a option is used, otherwise only best hits are recorded. 

*.concat.aln           If the -c option is chosen then a concatenated alignment is produced. This is FASTA formatted, each marker gene is aligned separately, and strings of X are used to fill spaces left by missing marker genes. This alignment is not trimmed in any way, so you may wish to process it further before phylogenetic analysis. 

log_file.txt          This is just a log file of some of the outputs produced by hmmsearch and clustal that is used for debugging purposes. 

In addition, for each .faa file in the input folder a .domout and .domout.parsed file is created. These are used if you re-run this tool with the -r flag. 


### Examples

To get the best hits to the markers A32,mcp,SFII,PolB,VLTF3, using 4 threads. 
>python ncldv_markersearch.py -i test_input -n test_run -t 4

To get all hits, including "secondary hits", or second-best hits:
>python ncldv_markersearch.py -i test_input -n test_run -t 4 -a

To get best hits and also generate a concatenated alignment: 
>python ncldv_markersearch.py -i test_input -n test_run -t 4 -c

To get best hits and generate an alignment for the markers A32,PolB,VLTF3,RNAPL, and RNAPS:
>python ncldv_markersearch.py -i test_input -n test_run -t 4 -c -m A32,PolB,VLTF3,RNAPL,RNAPS

#### Citation
M Moniruzzaman, CA Martinez-Gutierrez, AR Weinheimer, FO Aylward. Dynamic genome evolution and complex virocell metabolism of globally-distributed giant viruses. Nature Communications, 2020, 11(1): 1-11.


