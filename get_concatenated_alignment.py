import os
import sys
import subprocess
import re
import shlex
import operator
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

sequences = open(sys.argv[1], "r")
cogs = open(sys.argv[2], "r")
outfolder = sys.argv[3]
#trimvalue = sys.argv[4]
output = sys.argv[4]

record_dict = SeqIO.to_dict(SeqIO.parse(sequences, "fasta"))

taxon_list = []
for i in cogs.readlines():
	line = i.rstrip()
	tabs = line.split("\t")

	for j in tabs:	
		underscore = j.split("_")
		cog = underscore[1]
		taxon = underscore[0]
		taxon_list.append(taxon)
		#print underscore
		
	filename = os.path.join(outfolder, cog+".faa")
	handle = open(filename, "w")

	records = [record_dict[j] for j in tabs]
	SeqIO.write(records, handle, "fasta")
	handle.close()

taxon_set = set(taxon_list)
#print taxon_set

align_dict = defaultdict(list)
full_dict = {}
for i in os.listdir(outfolder):
	if i.endswith(".faa"):
		filename = os.path.join(outfolder, i)
		alignment = re.sub(".faa", ".aln", filename)
		trimmed = re.sub(".faa", ".trimal.aln", filename)
		cog = re.sub(".faa", "", i)
		#print "Aligning and trimming "+ cog +" and adding it to the concatenated alignment"
		cmd = "clustalo --force -i "+ filename +" -o "+ alignment 
		print cmd
		cmd2 = shlex.split(cmd)
		subprocess.call(cmd2, stdout=open("output.txt", "w"), stderr=open("err.txt", "w"))

	#	cmd = "trimal -gt "+ trimvalue +" -in "+ alignment +" -out "+ trimmed
	#	#print cmd
	#	cmd2 = shlex.split(cmd)
	#	subprocess.call(cmd2, stdout=open("output.txt", "w"), stderr=open("err.txt", "w"))

		seq_dict = SeqIO.to_dict(SeqIO.parse(alignment, "fasta"))
		values = seq_dict.values()
		first = values[0]
		length = len(first.seq)
		
		for taxon in taxon_set:
			entry = taxon +"_"+ cog
			if entry in seq_dict:
				align_dict[taxon].append(str(seq_dict[entry].seq))
			else:
				placeholder = "X" * length
				align_dict[taxon].append(placeholder)

outlist = []
for i in align_dict:
	record = SeqRecord(Seq("".join(align_dict[i]), IUPAC.protein), id=i)
	outlist.append(record)

SeqIO.write(outlist, output, "fasta")
print "Done"



















	


