import os, sys, subprocess, re, shlex, pandas, glob, operator, argparse
import numpy as np
from natsort import natsorted, ns
from collections import defaultdict
from operator import itemgetter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
#from Bio.Alphabet import IUPAC

speci_db = "hmm/gvogs.common.hmm" # Please ensure the PATH to the HMM dataset is correct
#cog_set = ["A32", "D5", "SFII", "mcp", "mRNAc", "PolB", "RNAPL", "RNAPS", "RNR", "VLTF3"]
#cog_set = ["A32", "SFII","mcp", "PolB", "VLTF3"]
#gvogset = ["GVOGm0054", "GVOGm0020", "GVOGm0003", "GVOGm0013", "GVOGm0022", "GVOGm0760", "GVOGm0023", "GVOGm0239", "GVOGm0890", "GVOGm0461", "GVOGm0095", "GVOGm0172", "GVOGm0088", "GVOGm0056", "GVOGm0018", "GVOGm0214", "GVOGm1574", "GVOGm0048", "GVOGm0129", "GVOGm0787", "GVOGm0036", "GVOGm0004", "GVOGm0152", "GVOGm0028", "GVOGm0160", "GVOGm0263", "GVOGm0115", "GVOGm0031", "GVOGm0189", "GVOGm0083", "GVOGm0694", "GVOGm0001", "GVOGm0041", "GVOGm0798", "GVOGm0027"]
#gvogset = ["GVOGm0890", "GVOGm0760", "GVOGm0461", "GVOGm0172", "GVOGm0054", "GVOGm0023", "GVOGm0013"]
score_dict = {"GVOGm0890":float(50), "GVOGm0760":float(80), "GVOGm0461":float(500), "GVOGm0172":float(80), "GVOGm0054":float(300), "GVOGm0023":float(550), "GVOGm0022":float(550), "GVOGm0013":float(160), "GVOGm0003":float(130)}


#score_dict = {"A32":float(80), "D5":float(80), "SFII":float(100), "mcp":float(80), "mRNAc":float(80), "PolB":float(150), "RNAPL":float(200), "RNAPS":float(200), "RNR":float(80), "VLTF3":float(80)}
#score_dict = defaultdict(lambda:float(20))
#score_dict = {}
#gvog_file = open("gvog_scores.list", "r")
#for i in gvog_file:
#	line = i.rstrip()
#	tabs = line.split("\t")
#	score_dict[tabs[0]] = float(tabs[1])

#################################################################
############# define hmm launcher function ######################
#################################################################
def hmm_launcher(folder, redo):
	print("Performing HMM search...")
	for files in os.listdir(folder):
		if files.endswith(".faa"):
			#print files
			input_file = os.path.join(folder, files)	
			dom_output = re.sub(".faa", ".domout", files)
			speci_dom_output = os.path.join(folder, dom_output)

			# run against the RNAP models
			cmd = "hmmsearch --cpu 16 -E 1e-3 --domtblout "+ speci_dom_output +" "+ speci_db + " " + input_file
			#print(cmd)
			cmd2 = shlex.split(cmd)
			if redo:
				pass
			else:
				subprocess.call(cmd2, stdout=open("log_file.txt", 'w'), stderr=open("log_file.txt", 'a'))

# end

def getprot(item):
	items = item.split("~")
	protein = items[0]
	cog = items[1]
	return(protein)

################################################################
###### Loop through and parse the checkm HMM output ############
################################################################

def hmm_parser(folder, suffix, combined_output):
	
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
			o.write("protein_id\taccession\tbest_hit\taln_start\taln_end\tscore\ttype\n")

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
					hit = re.sub("_", "", list1[3])
					#print(hit)

					score = float(list1[7])
					domain_evalue = float(list1[11])

					if score > bit_dict[ids] and domain_evalue < 1e-3:
						ids_hit = ids +"."+ hit
						start = int(list1[15])
						end   = int(list1[16])
						position_dict[ids_hit].append(start)
						position_dict[ids_hit].append(end)
						#print(ids_hit, score, domain_evalue, start, end)
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

################################################################
########## Define function for parsing HMMER3 output ###########
################################################################
def parse_domout(path_to_parsed_hmmfile, acc, protein_dict, cog_name, protein2dups):
	parsed = open(path_to_parsed_hmmfile, "r")
	done = {}
	protein2coords = defaultdict(list)
	protein2align_length = {}

	main_hit = "NAN"
	rnap_hits = []
	main_hits = []
	protein2cog = defaultdict(lambda:"NA")
	protein2acc = {}
	protein2score = {}
	protein2category = {}
	protein2length = {}

	for n in parsed.readlines():
		line = n.rstrip()
		tabs = line.split("\t")
		protein = tabs[0]
		annot = tabs[2].rstrip(".trim")
		if annot == cog_name:
			rnap_hits.append(protein)
			id_hit = protein +"~"+ annot
			hmm_score = float(tabs[5])
			category = tabs[6]
			start = int(tabs[3])
			end =   int(tabs[4])
			align_length = abs(end - start)

			if hmm_score > 20 and align_length > 20:

				record = protein_dict[protein]
				prot_length = len(record.seq)

				nr = acc +"_"+ annot

				protein2cog[protein]    = annot
				protein2acc[protein]    = acc
				protein2score[protein]  = hmm_score
				protein2length[protein] = prot_length

				protein2coords[id_hit].append(start)
				protein2coords[id_hit].append(end)
				#if annot == "PolB":
				#	print(start, end, id_hit)

				protein2align_length[id_hit] = align_length

				if category == "BH" and annot == cog_name:
				#if annot == cog_name:
					#main_hit = protein
					protein2dups[id_hit] = "single_besthit"
					main_hits.append(protein)

				else:
					#main_hit = protein
					protein2dups[id_hit] = "secondary_hit"
					main_hits.append(protein)

				protein2category[protein] = category

	parsed.close()
	#main_hits = set(main_hits)
	#main_hits2 = []
	#main_hits2 = [main_hits2.append(x) for x in main_hits if x not in main_hits2]
	return main_hits, protein2cog, protein2acc, protein2score, protein2length, protein2category, protein2coords, protein2align_length

def get_proteinsonreplicon(proteinid, record_list, prox):
	contig_name = re.sub("_\d*$", "", proteinid)
	final_list = []
	index = []
	indexzero=0
	ind = int(0)
	#record_list = natsorted(seqdict.keys())
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


# main function that runs the program
def run_program(inputdir, project, prox, cpus, redo, allhits, markerset, concat):

	merged = open(project+".full_output.txt", "w")
	merged.write("New_protein_name\tgenome\thit\tprotein_length\tbit_score\tnum_proteins_merged\thit_type\tprotein_ids\thmm_aln_coords\n")
	cog_out = open(project+".cogs.txt", "w")
	merged_proteins = open(project+".faa", "w")
	raw_output = open(project+".rawout.txt", "w")

	cog_set = markerset.split(",")
	#cog_set = gvogset
	
	if allhits:
		hitset = ['single_besthit', 'main_hit', 'secondary_hit']
	else:
		hitset = ['single_besthit', 'main_hit']
	
	hmm_launcher(inputdir, redo)
	hmm_parser(inputdir, ".domout", raw_output)
	print("Compiling results...")

	final_proteins = []
	marker_tally = defaultdict(int)
	#exceptions = open("exceptions.txt", "w")
	tally = 0
	protein_tally = []

	df = pandas.DataFrame()
	
	for i in os.listdir(inputdir):
		if i.endswith(".faa"):

			markercount = defaultdict(float)
			protein_file = os.path.join(inputdir, i)
			gff_file = re.sub(".faa$", ".gff", protein_file)
			domout = re.sub(".faa$", ".domout", protein_file)
			parsed = re.sub(".faa$", ".domout.parsed", protein_file)
			acc = re.sub(".faa$", "", i)		
			merged_protein_list = []

			# get a dictionary of protein sequences
			seq_handle = open(protein_file, "r")
			seq_dict = SeqIO.to_dict(SeqIO.parse(seq_handle, "fasta"))
			record_list = natsorted(seq_dict.keys())
			#orf_set = [record.id for record in seq_dict.values()]
			prot2protlist = defaultdict(list)
			num_proteins = defaultdict(lambda:int(1))
			prot2loc = defaultdict(list)
			prot2locrange = defaultdict(list)
			# parse domout file and get protein hits and coordinates
			for cog in cog_set:
				#print(cog)
				protein2dups = defaultdict(lambda:"hits")

				main_hits, protein2cog, protein2acc, protein2score, protein2length, protein2category, protein2coords, protein2align_length = parse_domout(parsed, acc, seq_dict, cog, protein2dups)
				#print(cog, main_hits)
				
				for main_hit in main_hits:
					
					if main_hit in protein_tally:
						pass
					else:
					
						protein_tally.append(main_hit)
						#print("before")
						
						orf_set = get_proteinsonreplicon(main_hit, record_list, prox)
						#print("after")
						#print(main_hit, cog, len(main_hits), len(orf_set))
						if main_hit == "NAN":
							#print(acc, cog, main_hit)
							pass

						else:
							already_done = []

							#print(main_hit)
							prot2protlist[main_hit].append(main_hit)
							#print(len(prot2protlist[main_hit]))
							id_hit1 = main_hit +"~"+ cog
							range1 = protein2coords[id_hit1]

							#print(id_hit1, orf_set)
							
							r1 = range(range1[0], range1[1])
							meanloc1 = np.mean(range1)
							locrange1 = str(range1[0]) +"-"+ str(range1[1])
							prot2locrange[main_hit].append(locrange1)
							prot2loc[main_hit].append(meanloc1)

							orf_set.remove(main_hit)

							for m in protein_tally:
								if m in orf_set:
									orf_set.remove(m)


							orf_set = set(orf_set)
							for d in orf_set:

								if protein2cog[d] == cog:

									id_hit2 = d +"~"+ cog
									range2 = protein2coords[id_hit2]
									r2 = range(range2[0], range2[1])

									meanloc2 = np.mean(range2)
									locrange2 = str(range2[0]) +"-"+ str(range2[1])
									
									#prot2loc[rnap].append(meanloc2)
										
									set1 = set(r1)
									inter = set1.intersection(r2)

									#print(d, inter)
									#print(id_hit1, id_hit2, r1, r2, range1, range2, len(inter))

									if int(len(inter)) > 10:
										protein2dups[id_hit2] = "secondary_hit"
										#print(id_hit1, id_hit2, r1, r2, range1, range2, len(inter))
									else:
										protein_tally.append(d)
										
										protein2dups[id_hit1] = "main_hit"
										protein2dups[id_hit2] = "secondary_hit"
										#print(main_hit, id_hit1, id_hit2, d, protein2dups[id_hit1], protein2dups[id_hit2])
										protein_tally.append(id_hit1)
										protein_tally.append(id_hit2)

										minrange = min(range1 + range2)
										maxrange = max(range1 + range2)
										protein2coords[id_hit1] = [minrange, maxrange]
										#print(minrange, maxrange)
										prot2locrange[main_hit].append(locrange2)
									
										#protein2align_length[id_hit1] = abs(maxrange - minrange)
										#protein2align_length[main_hit] = int(protein2align_length[main_hit]) + int(protein2align_length[d])
										protein2length[main_hit] = int(protein2length[main_hit]) + int(protein2length[d])
										protein2score[main_hit] = float(protein2score[main_hit]) + float(protein2score[d])
										prot2protlist[main_hit].append(d)
										prot2loc[main_hit].append(meanloc2)
										num_proteins[id_hit1] +=1
										#print(id_hit1, id_hit2, prot2protlist, prot2protlist[main_hit])
										merged_protein_list.append(id_hit2)

										#if cog == "PolB":
										#	print(id_hit1, id_hit2, num_proteins[id_hit1])


				all_besthits = [p for p in protein2dups.keys() if protein2dups[p] in ["single_besthit", "main_hit"]] 
				all_scores = [protein2score[getprot(p)] for p in all_besthits]
				if len(all_scores) > 0:
					max_index, max_value = max(enumerate(all_scores), key=operator.itemgetter(1))
					best_hit = all_besthits[max_index]
					other_hits = [j for j in all_besthits if j != best_hit]
					for o in other_hits:
						protein2dups[o] = "secondary_hit"
					#all_scores   = [protein2score[p] for p in protein2dups.keys() if protein2dups[p] in ["single_besthit", "main_hit"]]
				
			#	print(all_besthits)
			#	print(all_scores)
			#	print(max_value, max_index)
							
				hit_tally = defaultdict(int)		
				best_hit_bit = defaultdict(float)
				for item in protein2dups:
					#print item
					if item in merged_protein_list:
						#print(item)
						pass
					
					elif protein2dups[item] in hitset:
						items = item.split("~")

						protein = items[0]
						hit = items[1]
						
			#			bit = protein2score[protein]
						protlist = prot2protlist[protein]
						
			#			if bit > best_hit_bit[hit]:
			#				protein2dups[protein] = 
						
						loc_list = [float(loc) for loc in prot2loc[protein]]
						index_list = [i[0] for i in sorted(enumerate(loc_list), key=lambda x:x[1])]
						sorted_loc_list = [i[1] for i in sorted(enumerate(loc_list), key=lambda x:x[1])]

						sorted_prot_list = [protlist[index] for index in index_list]
						prot_str = ";".join(sorted_prot_list)

						range_list = prot2locrange[protein]
						sorted_range_list = [range_list[index] for index in index_list]
						range_str = ";".join(sorted_range_list) 
						
						#print(items, protein, protein2dups[item], protein2score[protein], sorted_loc_list, prot2loc[protein], range_str)
						#print hit, sorted_prot_list, sorted_loc_list, sorted(enumerate(loc_list), key=lambda x:x[1])

						loc_str = ";".join([str(n) for n in sorted_loc_list])
						acc = protein2acc[protein]
						final_name_str = re.sub("_", ".", acc) +"_"+ hit
						hit_tally[final_name_str] +=1
						ptally = str(hit_tally[final_name_str])
						final_name = final_name_str +".copy"+ptally
						#print(final_name)
						if protein2score[protein] > score_dict[hit]:
							#merged.write(final_name +"\t"+ acc +"\t"+ hit +"\t"+ str(protein2length[protein]) +"\t"+ str(protein2score[protein]) +"\t"+ str(protein2align_length[item]) +"\t"+ str(num_proteins[item]) +"\t"+ protein2dups[item] +"\t"+ prot_str +"\t"+ range_str +"\n")
							merged.write(final_name +"\t"+ acc +"\t"+ hit +"\t"+ str(protein2length[protein]) +"\t"+ str(round(protein2score[protein], 1)) +"\t"+ str(num_proteins[item]) +"\t"+ protein2dups[item] +"\t"+ prot_str +"\t"+ range_str +"\n")


							markercount[hit] +=1
							if len(sorted_prot_list) > 1:
								#print(hit, protein, protlist)
								tally = tally + len(sorted_prot_list)

								newrecord = SeqRecord(Seq(""), id=final_name, name=protein+" JOINED", description=protein2acc[protein] +" JOINED_PROTEIN")
								for fragment in sorted_prot_list:
									subrecord = seq_dict[fragment]
									subseq = subrecord.seq
									subseq = re.sub("\*", "", str(subseq))
									newrecord.seq = newrecord.seq +""+ subseq

								final_proteins.append(newrecord)

							else:
								tally +=1
								record = seq_dict[protein]
								#print(protein, item, num_proteins[item])
								#record.description = protein2acc[protein] +" JOINED_PROTEIN" 
								record.id = final_name
								final_proteins.append(record)

			s1 = pandas.DataFrame(pandas.Series(markercount, name = acc))
			df = pandas.concat([df, s1], axis=1, sort=True)
								
	names = [i.id for i in final_proteins]
	for cog in cog_set:
		name_set = [i for i in names if cog in i]
		name_str = "\t".join(name_set)
		cog_out.write(name_str +"\n")

	final_records = []
	for seqrecord in final_proteins:
		seq = seqrecord.seq
		newseq = Seq("".join([n for n in seq if n != "*"]))
		#print(newseq)
		newrecord = SeqRecord(newseq, id=seqrecord.id, name=seqrecord.name, description=seqrecord.description)
		final_records.append(newrecord)

	SeqIO.write(final_records, merged_proteins, "fasta")

	merged.close()
	cog_out.close()
	merged_proteins.close()
	#print tally

	df2 = df.transpose()
	df2.fillna(0, inplace=True)
	df2.to_csv(project+".table.tsv", sep="\t", index_label="genome")

	if concat and allhits:
		print("Cannot choose -a and -c flags together! No concatenated alignment produced.")
	elif concat and not allhits:

		if os.path.isdir(project+"_alignments"):
			pass
		else:
			os.mkdir(project+"_alignments")
		
		record_dict = SeqIO.to_dict(SeqIO.parse(project+".faa", "fasta"))
		taxon_list = []
		cogs = open(project+".cogs.txt", "r")

		for i in cogs.readlines():
			line = i.rstrip()
			tabs = line.split("\t")

			for j in tabs:
				#print(j)
				underscore = j.split("_")
				cog = underscore[1]
				taxon = underscore[0]
				taxon_list.append(taxon)
				#print underscore
		
			filename = os.path.join(project+"_alignments", cog+".faa")
			handle = open(filename, "w")

			records = [record_dict[j] for j in tabs]
			SeqIO.write(records, handle, "fasta")
			handle.close()

		taxon_set = set(taxon_list)
		#print taxon_set

		print("Generating alignments with Clustal Omega...")
		align_dict = defaultdict(list)
		full_dict = {}
		for i in os.listdir(project+"_alignments"):
			if i.endswith(".faa"):
				filename = os.path.join(project+"_alignments", i)
				alignment = re.sub(".faa", ".aln", filename)
				cog = re.sub(".faa", "", i)
				#print "Aligning and trimming "+ cog +" and adding it to the concatenated alignment"
				cmd = "clustalo --threads "+ cpus +" --force -i "+ filename +" -o "+ alignment 
				#print(cmd)
				cmd2 = shlex.split(cmd)
				subprocess.call(cmd2, stdout=open("log_file.txt", "a"), stderr=open("log_file.txt", "a"))

				seq_dict = SeqIO.to_dict(SeqIO.parse(alignment, "fasta"))
				values = list(seq_dict.values())
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
			record = SeqRecord(Seq("".join(align_dict[i])), id=i)
			record.description = "concatenated alignment of "+ markerset
			outlist.append(record)

		SeqIO.write(outlist, project+".concat.aln", "fasta")




########################################################################
##### use argparse to run through the command line options given #######
########################################################################
def main(argv=None):

	args_parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="NCLDV_Markersearch: A script for identifying phylogenetic marker genes in NCLDV and generating concatenated alignments \nFrank O. Aylward, Virginia Tech Department of Biological Sciences <faylward at vt dot edu>", epilog='*******************************************************************\n\n*******************************************************************')
	args_parser.add_argument('-i', '--input', required=True, help='Input folder of FASTA file (ending in .fna, .fa, or .fasta)')
	args_parser.add_argument('-n', '--name', required=True, help='project name prefix for output files')
	args_parser.add_argument('-p', '--proximity', required=False, default=int(5), help='number of genes to look up- and downstream of hits to join genes (default=5)')
	args_parser.add_argument('-t', '--cpus', required=False, default=str(1), help='number of cpus to use for the HMMER3 search')
	args_parser.add_argument('-m', '--markerset', required=False, default=str('GVOGm0890,GVOGm0760,GVOGm0461,GVOGm0172,GVOGm0054,GVOGm0023,GVOGm0013'), help='Markers to use. Must be comma-separated list of the following: GVOGm0890,GVOGm0760,GVOGm0461,GVOGm0172,GVOGm0054,GVOGm0023,GVOGm0013,GVOGm0022,GVOGm0003. Default is GVOGm0890,GVOGm0760,GVOGm0461,GVOGm0172,GVOGm0054,GVOGm0023,GVOGm0013')
	args_parser.add_argument('-r', '--redo', type=bool, default=False, const=True, nargs='?', help='run without re-launching prodigal and HMMER3 (for quickly re-calculating outputs with different parameters if you have already run once)')
	args_parser.add_argument('-c', '--concat', type=bool, default=False, const=True, nargs='?', help='In addition to finding marker genes, generated a concatenated alignment of the best hits (not compatible with the -a option)')
	args_parser.add_argument('-a', '--allhits', type=bool, default=False, const=True, nargs='?', help='Provide all hits (default is to provide only best hits to each marker gene)')
	args_parser.add_argument('-v', '--version', action='version', version='ncldv_markersearch v. 1.1')
	args_parser = args_parser.parse_args()

	# set up object names for input/output/database folders
	inputdir = args_parser.input
	project = args_parser.name
	prox = int(args_parser.proximity)
	cpus = args_parser.cpus
	#contiglevel = args_parser.contiglevel
	redo = args_parser.redo
	allhits = args_parser.allhits
	markerset = args_parser.markerset
	concat = args_parser.concat
	
	run_program(inputdir, project, prox, cpus, redo, allhits, markerset, concat)
	return 0

if __name__ == '__main__':
	status = main()
	sys.exit(status)

# end








