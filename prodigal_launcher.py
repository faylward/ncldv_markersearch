import os
import sys
import subprocess
import re
import shlex

input_folder = sys.argv[1]
outfolder = sys.argv[2]

for i in os.listdir(input_folder):
	if i.endswith(".fna") or i.endswith(".fa"):
		name = re.sub("_genomic.fna", "", i)
		fasta = os.path.join(input_folder, i)
		faa = os.path.join(outfolder, name+".faa")
		gff = os.path.join(outfolder, name+".gff")
		genes = os.path.join(outfolder, name+".genes.fna")

		cmd = "prodigal -i "+ fasta +" -f gff -o "+ gff +" -a "+ faa +" -d "+ genes
		print i
		cmd2 = shlex.split(cmd)
		subprocess.call(cmd2, stdout=open("prodigal.out", "w"), stderr=open("prodigal.err", "w"))


