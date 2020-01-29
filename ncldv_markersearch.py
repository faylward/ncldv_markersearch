import os
import sys
import subprocess
import re
import shlex
import pandas
import glob
import operator
import numpy as np
from natsort import natsorted, ns
from collections import defaultdict
from operator import itemgetter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

working_dir = sys.argv[1]
prox = int(sys.argv[2])

speci_db = "hmm/NCLDV.hmm" # Please ensure the PATH to the HMM dataset is correct

cog_set = ["A32", "D5", "SFII", "mcp", "mRNAc", "PolB", "RNAPL", "RNAPS", "RNR", "VLTF3"]
#cog_set = ["A32", "SFII", "mcp", "PolB", "VLTF3"]

combined_output = open("raw_annotatino_output.txt", "w")

merged_proteins = open("ncldv_markersearch.faa", "w")
final_proteins = []

marker_tally = defaultdict(int)
exceptions = open("exceptions.txt", "w")
tally = 0

score_dict = {"A32":float(80), "D5":float(80), "SFII":float(100), "mcp":float(80), "mRNAc":float(80), "PolB":float(150), "RNAPL":float(200), "RNAPS":float(200), "RNR":float(80), "VLTF3":float(80)}

#################################################################
############# define hmm launcher function ######################
#################################################################
def hmm_launcher(folder):
	for files in os.listdir(folder):
		if files.endswith(".faa"):
			#print files
			input_file = os.path.join(folder, files)	
			dom_output = re.sub(".faa", ".domout", files)
			speci_dom_output = os.path.join(folder, dom_output)

			# run against the RNAP models
			cmd = "hmmsearch --cpu 16 -E 1e-5 --domtblout "+ speci_dom_output +" "+ speci_db + " " + input_file
			#print cmd
			cmd2 = shlex.split(cmd)
			subprocess.call(cmd2, stdout=open("hmm.out", 'w'), stderr=open("error_file.txt", 'a'))

# end
hmm_launcher(working_dir)


################################################################
###### Loop through and parse the checkm HMM output ############
################################################################

def hmm_parser(folder, suffix, output):
	
	record_list = []
	score_list = {}
	prot_list = []
	#combined_output = open(output, "w")
	combined_output.write("protein\tacc\thit\tstart\tend\taln_length\tscore\tcategory\n")
	hits = []
	bit_dict = {}

	for filenames in os.listdir(folder):
		if filenames.endswith(suffix):

			acc = re.sub(suffix, "", filenames)

			f = open(folder+"/"+filenames, 'r')
			o = open(folder+"/"+filenames+".parsed", 'w')

			faa_file = re.sub(".domout", ".faa", filenames)
			protein_dict = SeqIO.to_dict(SeqIO.parse(os.path.join(folder, faa_file), "fasta"))

			hit_dict = {}
			start_dict = {}
			end_dict = {}
			bit_dict = defaultdict(int)
			hit_type = {}
			marker_dict = {}
			position_dict = defaultdict(list)

			for line in f.readlines():
				if line.startswith("#"):
					pass
				else:
					newline = re.sub( '\s+', '\t', line)
					list1 = newline.split('\t')
					ids = list1[0]
					hit = list1[3]

					ids_hit = ids +"."+ hit

					start = int(list1[15])
					end   = int(list1[16])
					position_dict[ids_hit].append(start)
					position_dict[ids_hit].append(end)

					#print start, end

					score = float(list1[7])
					domain_evalue = float(list1[11])

					if score > bit_dict[ids] and domain_evalue < 1e-5:
						hit_dict[ids] = hit
						start_dict[ids] = start
						end_dict[ids] = end
						bit_dict[ids] = score

			bit_sorted = sorted(bit_dict.items(), key=operator.itemgetter(1), reverse=True)
			output_list = []
			for item in bit_sorted:
				entry = item[0]
				score = item[1]
				if score > 0:
					#print entry, item, filenames
					ids_hit = entry +"."+ hit_dict[entry]
					output_list.append(entry +"\t"+ str(hit_dict[entry]) +"\t"+ str(min(position_dict[ids_hit])) +"\t"+ str(max(position_dict[ids_hit])) +"\t"+ str(bit_dict[entry]) )

			hit_profile = defaultdict(int)
			done = []
			for line in output_list:
				line1 = line.rstrip()
				tabs = line1.split("\t")
				ids = tabs[0]
				
				record = protein_dict[ids]
				record_list.append(record)

				hits.append(ids)
				cog = tabs[1]
				start = tabs[2]
				end = tabs[3]
				aln_length = str(abs(float(end) - float(start)))

				score = tabs[4]
				nr = acc +"_"+ cog

				if nr in done:
					combined_output.write(ids +"\t"+ acc +"\t"+ cog +"\t"+ start +"\t"+ end +"\t"+ aln_length +"\t"+ score +"\tNH\n")
					o.write(ids +"\t"+ acc +"\t"+ cog +"\t"+ start +"\t"+ end +"\t"+ score +"\tNH\n")
				else:
					combined_output.write(ids +"\t"+ acc +"\t"+ cog +"\t"+ start +"\t"+ end +"\t"+ aln_length +"\t"+ score +"\tBH\n")
					o.write(ids +"\t"+ acc +"\t"+ cog +"\t"+ start +"\t"+ end +"\t"+ score +"\tBH\n")
					done.append(nr)
			o.close()
#	output = open("output.faa", "w")
#	SeqIO.write(record_list, output, "fasta")

#parse speci outputs
speci_df = hmm_parser(working_dir, ".domout", "all_hmmout.txt")


################################################################
########## Define function for parsing HMMER3 output ###########
################################################################
def parse_domout(path_to_parsed_hmmfile, acc, protein_dict, cog_name):
	parsed = open(path_to_parsed_hmmfile, "r")
	done = {}
	protein2coords = defaultdict(list)
	protein2align_length = {}

	main_hit = "NAN"
	rnap_hits = []
	protein2cog = defaultdict(lambda:"NA")
	protein2acc = {}
	protein2score = {}
	protein2category = {}
	protein2length = {}

	for n in parsed.readlines():
		line = n.rstrip()
		tabs = line.split("\t")
		protein = tabs[0]
		annot = tabs[2]
		if annot == cog_name:
			rnap_hits.append(protein)
			id_hit = protein +"|"+ annot
			hmm_score = float(tabs[5])
			category = tabs[6]

			if hmm_score > 1:

				start = int(tabs[3])
				end =   int(tabs[4])

				record = protein_dict[protein]
				prot_length = len(record.seq)

				nr = acc +"_"+ annot

				protein2cog[protein]    = annot
				protein2acc[protein]    = acc
				protein2score[protein]  = hmm_score
				protein2length[protein] = prot_length

				protein2coords[id_hit].append(start)
				protein2coords[id_hit].append(end)

				align_length = abs(end - start)
				protein2align_length[id_hit] = align_length

				if category == "BH" and annot == cog_name:
					main_hit = protein
					protein2dups[id_hit]

				protein2category[protein] = category

	parsed.close()
	return main_hit, protein2cog, protein2acc, protein2score, protein2length, protein2category, protein2coords, protein2align_length

def get_proteinsonreplicon(proteinid, seqdict):
	contig_name = re.sub("_\d*$", "", proteinid)
	final_list = []
	index = []
	indexzero=0
	ind = int(0)
	record_list = natsorted(seqdict.keys())
	#print record_list
	for record in record_list:
		#print record
		if contig_name in record:
			final_list.append(record)
			index.append(ind)
			if proteinid == record:
				indexzero = ind
			ind +=1

	#prox = int(5) # number of genes to look in front and in back of gene
	start = indexzero - prox
	end = indexzero + prox + 1

	if start < 0:
		start = 0
	if end > len(final_list):
		end = len(final_list)

	protein_list = final_list[start:end]

	#print proteinid, indexzero, protein_list	
	return(protein_list)

merged = open("full_output.txt", "w")
merged.write("New_protein_name\tseed_hit\tgenome\thit\tprotein_length\tbit_score\talignment_length\tnum_proteins_merged\tprotein_ids\thit_coords\n")
cog_out = open("cogs.txt", "w")
for i in os.listdir(working_dir):
	if i.endswith(".faa"):

		protein_file = os.path.join(working_dir, i)
		gff_file = re.sub(".faa", ".gff", protein_file)
		domout = re.sub(".faa", ".domout", protein_file)
		parsed = re.sub(".faa", ".domout.parsed", protein_file)
		acc = re.sub(".faa", "", i)		

		# get a dictionary of protein sequences
		seq_handle = open(protein_file, "r")
		seq_dict = SeqIO.to_dict(SeqIO.parse(seq_handle, "fasta"))
		#orf_set = [record.id for record in seq_dict.values()]

		# parse domout file and get protein hits and coordinates
		for cog in cog_set:
			protein2dups = defaultdict(lambda:"single_besthit")

			rnap, protein2cog, protein2acc, protein2score, protein2length, protein2category, protein2coords, protein2align_length = parse_domout(parsed, acc, seq_dict, cog)
			orf_set = get_proteinsonreplicon(rnap, seq_dict)

			if rnap == "NAN":
				#print acc, cog, rnap
				pass
			else:
				already_done = []
				num_proteins = defaultdict(lambda:int(1))
				prot2protlist = defaultdict(list)
				prot2loc = defaultdict(list)

				prot2protlist[rnap].append(rnap)

				id_hit1 = rnap +"|"+ cog
				range1 = protein2coords[id_hit1]
				r1 = range(range1[0], range1[1])
				meanloc1 = np.mean(range1)
				#meanloc1 = str(range1[0]) +"-"+ str(range1[1])
				prot2loc[rnap].append(meanloc1)

				orf_set.remove(rnap)

				for d in orf_set:
					if protein2cog[d] == cog:
						id_hit2 = d +"|"+ cog
						range2 = protein2coords[id_hit2]
						r2 = range(range2[0], range2[1])

						meanloc2 = np.mean(range2)
						#meanloc2 = str(range2[0]) +"-"+ str(range2[1])
						#prot2loc[rnap].append(meanloc2)
								
						set1 = set(r1)
						inter = set1.intersection(r2)

						if int(len(inter)) > 50:
							protein2dups[id_hit2] = "false_hit"
						else:
							protein2dups[id_hit1] = "main_hit"
							protein2dups[id_hit2] = "secondary_hit"

							minrange = min(range1 + range2)
							maxrange = max(range1 + range2)
							protein2coords[id_hit1] = [minrange, maxrange]
									
							protein2align_length[id_hit1] = abs(maxrange - minrange)
							protein2length[rnap] = int(protein2length[rnap]) + int(protein2length[d])
							protein2score[rnap] = float(protein2score[rnap]) + float(protein2score[d])
							prot2protlist[rnap].append(d)
							prot2loc[rnap].append(meanloc2)
							num_proteins[id_hit1] +=1



			for item in protein2dups:
				#print item
				if protein2dups[item] == "main_hit" or protein2dups[item] == "single_besthit": #or protein2dups[item] == "NEXT" or protein2dups[item] == "SECO":
					items = item.split("|")
					protein = items[0]
					hit = items[1]

					protlist = prot2protlist[protein]

					loc_list = [float(loc) for loc in prot2loc[protein]]
					index_list = [i[0] for i in sorted(enumerate(loc_list), key=lambda x:x[1])]
					sorted_loc_list = [i[1] for i in sorted(enumerate(loc_list), key=lambda x:x[1])]

					sorted_prot_list = [protlist[index] for index in index_list]
					prot_str = ";".join(sorted_prot_list)

					#print hit, sorted_prot_list, sorted_loc_list, sorted(enumerate(loc_list), key=lambda x:x[1])

					loc_str = ";".join([str(n) for n in sorted_loc_list])
					acc = protein2acc[protein]
					final_name = re.sub("_", ".", acc) +"_"+ hit

					merged.write(final_name +"\t"+ protein +"\t"+ acc +"\t"+ hit +"\t"+ str(protein2length[protein]) +"\t"+ str(protein2score[protein]) +"\t"+ str(protein2align_length[item]) +"\t"+ str(num_proteins[item]) +"\t"+ prot_str +"\t"+ loc_str +"\n")

					if protein2score[protein] > score_dict[hit]:

						if len(sorted_prot_list) > 1:
							tally = tally + len(sorted_prot_list)

							newrecord = SeqRecord(Seq("", IUPAC.protein), id=final_name, name=protein+" JOINED", description=protein2acc[protein])
							for fragment in sorted_prot_list:
								subrecord = seq_dict[fragment]
								subseq = subrecord.seq
								subseq = re.sub("\*", "", str(subseq))
								newrecord.seq = newrecord.seq +""+ subseq

							final_proteins.append(newrecord)

						else:
							tally +=1
							record = seq_dict[protein]
							record.id = final_name
							final_proteins.append(record)

names = [i.id for i in final_proteins]
for cog in cog_set:
	name_set = [i for i in names if cog in i]
	name_str = "\t".join(name_set)
	cog_out.write(name_str +"\n")

SeqIO.write(final_proteins, merged_proteins, "fasta")
#print tally





































