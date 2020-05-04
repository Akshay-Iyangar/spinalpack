import sys
import networkx as nx

class GoConsistency:
	def __init__(self, match_filename, gene_to_go_filename1, gene_to_go_filename2, out_filename):

		matchfile = open(match_filename)
		gene_to_go_file1 = open(gene_to_go_filename1)
		gene_to_go_file2 = open(gene_to_go_filename2)
		outfile = open(out_filename, "w")

		#assumes gene_to_go file format is: X|X1|X2... GO1|GO2...
		#where X is the name matching the protein name in the matchfile
		#X1, X2... are assigned names to X, if any, and GO1, GO2... are GO ids
		gene_to_go1 = {}
		gene_naming1 = {}
		for line in gene_to_go_file1:
			cols = line.split()
			#store synonym info
			prot_names = cols[0].split("|")
			gene_naming1[prot_names[0]] = []
			for prot in prot_names:
				gene_naming1[prot_names[0]].append(prot)
			#map each gene to a set of GO ids
			gene_to_go1[prot_names[0]] = set()
			if len(cols) > 1:
				go_ids = cols[1].split("|")
				for go_id in go_ids:
					gene_to_go1[prot_names[0]].add(go_id)
		gene_to_go2 = {}
		gene_naming2 = {}
		for line in gene_to_go_file2:
			cols = line.split()
			#store synonym info
			prot_names = cols[0].split("|")
			gene_naming2[prot_names[0]] = []
			for prot in prot_names:
				gene_naming2[prot_names[0]].append(prot)
			#map each gene to a set of GO ids
			gene_to_go2[prot_names[0]] = set()
			if len(cols) > 1:
				go_ids = cols[1].split("|")
				for go_id in go_ids:
					gene_to_go2[prot_names[0]].add(go_id)

		#assumes matchfile format is: X Y
		#at each line where X is a protein from ppi1 Y is from ppi2	
		outlines = []
		overlap_sizes = {}
		for line in matchfile:
			if line[0] == "!":
				continue
			cols = line.split()	
			common_set = gene_to_go1[cols[0]].intersection(gene_to_go2[cols[1]])
			if len(common_set) in overlap_sizes:
				overlap_sizes[len(common_set)] += 1
			else:
				overlap_sizes[len(common_set)] = 1
			if len(common_set) >= 1:
				outlines.append('|'.join(gene_naming1[cols[0]]) + " " + \
				'|'.join(gene_naming2[cols[1]]) + " " + '|'.join(list(common_set)) + "\n")	

		outfile.write("!This file contains all matches that share at least one GO category.\n")
		overlap_score = 0
		for size in overlap_sizes.iterkeys():
			outfile.write("! %d matches with %d GO category overlaps.\n" % (overlap_sizes[size], size)) 
			overlap_score += size * overlap_sizes[size]
		outfile.write("!Total GO overlap score: %d\n" % overlap_score)
				
		for line in outlines:
			outfile.write(line)
		outfile.close()
					
if __name__ == "__main__":
	#execute go_association.py before this script if necessary
	#sample run:
	#python pyscripts/go_consistency.py data/hs-sc/hsapi-scere-spinalmatch-original7 data/hsapi/hsapi_gene_to_go.txt 
	#data/scere/scere_gene_to_go.txt data/hs-sc/hsapi-scere-spinalmatch-gooverlap7				
	inputdata = GoConsistency(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])		
																					
