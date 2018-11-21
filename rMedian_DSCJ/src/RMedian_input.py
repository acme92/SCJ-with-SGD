from gurobipy import *

def read_file(filename):
	string = open(filename, "r").read()
	string_list = string.split("\n")
	string_list = [line for line in string_list if line and line[0] != '#'] #Read line only if it is nonempty and not a comment.
	return string_list

#---------------------------------------
#SPECIES TREE
#---------------------------------------
def sptree(tree_file, tree_dict, rev_tree_dict):
	string_list = read_file(tree_file)
	for line in string_list:
		line = line.split("\t")
		ch, par = line[0], line[1]
		if par not in tree_dict:
			tree_dict[par] = [ch]
		else:
			tree_dict[par].append(ch)
		rev_tree_dict[ch] = par
	for species in rev_tree_dict:
		if species in tree_dict:
			median = species

	return tree_dict, rev_tree_dict, median

#---------------------------------------
#GENE CONTENT
#---------------------------------------
#Takes a line from input file and stores in two nested dictionaries
#1. cont_dict: Key = species, Value = dictionary of genes -> Key = gene name, Value = location in genome as tuple (scaffold index, gene index)
#Data types: species and gene names are str, scaffold and gene index are int
#
#Second dict required to store gene orders in extant genomes
#2. cont_dict_by_idx: Key = species, Value = dictionary of genes -> Key = scaffold index, Value = list of genes in the scaffold in order
##Data types: species and gene names are str, scaffold index is int
def update_gene_content(line, cont_dict, cont_dict_by_idx, prev_scidx, prev_species, gidx):
	species, scidx, gene_name, gene_orientation = line[0], int(line[1]), line[2], line[3]
	if species not in cont_dict:
		cont_dict[species] = {}
		cont_dict_by_idx[species] = {}

	if scidx != prev_scidx or species != prev_species:		#This check is only required for extant genomes since
		prev_scidx = scidx								
		#We are also storing the gene order.
		prev_species = species
		gidx = 0
		
		cont_dict_by_idx[species][scidx] = {}
		cont_dict_by_idx[species][scidx]['gene_list'] = []
		#cont_dict_by_idx[species][scidx]['type'] = 'L' 		#Redo this line

	if gene_name in cont_dict[species]:
		cont_dict_by_idx[species][scidx]['type'] = 'C'
	else:
		cont_dict_by_idx[species][scidx]['type'] = 'L'	
		cont_dict[species][gene_name] = (scidx, gidx)
	cont_dict_by_idx[species][scidx]['gene_list'].append(gene_orientation + gene_name)	
	gidx += 1

	return cont_dict, cont_dict_by_idx, prev_scidx, prev_species, gidx

def gene_cont(cont_file, cont_dict, cont_dict_by_idx, relevant_species):
	string_list = read_file(cont_file)
	prev_scidx = None
	prev_species = None
	gidx = 0
	for line in string_list:
		line = line.split('\t')
		if line[0] in relevant_species:
			cont_dict, cont_dict_by_idx, prev_scidx, prev_species, gidx \
				= update_gene_content(line, cont_dict, cont_dict_by_idx, prev_scidx, prev_species, gidx)

	return cont_dict, cont_dict_by_idx

#---------------------------------------
#ORTHOLOGY RELATIONS
#---------------------------------------
#Takes a line from input file and stores in two dictionaries
#1. tree_dict: Key = ancestor species (parent), Value = List of descendant species (children)
#2. rev_tree_dict: Key = descendant species (child), Value = ancestor species (parent)
def update_tree(line, tree_dict, rev_tree_dict):
	parent, child, gpar, gch, fam = line[0], line[1], line[2], line[3], line[4]
	if parent not in tree_dict:
		tree_dict[parent] = [child]
	else:	
		if child not in set(tree_dict[parent]):
			tree_dict[parent].append(child)
	if child not in rev_tree_dict:
		rev_tree_dict[child] = parent
	return tree_dict, rev_tree_dict

#Takes a line from input file and stores in two dictionaries
#1. ortho_reln: Key = gene name (gpar) in ancestor (parent), Value = List of gene names (gch) in descendant (child)
#2. rev_ortho_reln: Key = gene names (gch) in descendant (child), Value = gene name (gpar) in ancestor (parent)
def update_ortho_reln(line, ortho_reln, rev_ortho_reln, cont_dict):
	parent, child, gpar, gch, fam = line[0], line[1], line[2], line[3], line[4]
	if ((parent, child)) not in ortho_reln:
		ortho_reln[(parent, child)] = {}
		rev_ortho_reln[(child, parent)] = {}
	if gpar not in ortho_reln[(parent, child)]:
		ortho_reln[(parent, child)][gpar] = []
	
	p_flag, c_flag = 0, 0	#Existence of parent gene and child gene in the gene content of respective species
	if gpar in cont_dict[parent]:
		p_flag = 1
	if gch in cont_dict[child]:
		c_flag = 1
	if gpar not in ortho_reln[(parent, child)]:
		ortho_reln[(parent, child)][gpar] = []
	ortho_reln[(parent, child)][gpar].append(gch)
	rev_ortho_reln[(child, parent)][gch] = gpar

	return ortho_reln, rev_ortho_reln

def ortho(ortho_file, tree_dict, rev_tree_dict, ortho_reln, rev_ortho_reln, cont_dict):
	string_list = read_file(ortho_file)
	for line in string_list:
		line = line.split("\t")
		parent, child = line[0], line[1]
		if parent in tree_dict and child in rev_tree_dict:
			ortho_reln, rev_ortho_reln = update_ortho_reln(line, ortho_reln, rev_ortho_reln, cont_dict)

	return ortho_reln, rev_ortho_reln

#---------------------------------------
#ADJACENCY SETS
#---------------------------------------
#Takes a genome and output an adjacency set (for extant genomes)
#All adjacencies are a tuple:
#	( (gene1_name, location_1, orientation_1), (gene2_name, location_2, orientation_2) )
#	where location is an int tuple: (scaffold index, gene index) for the gene. 
def get_adj_set(species):
	adj_set = set()
	for cidx in species:
		gseq = species[cidx]['gene_list']
		for gidx in range(len(gseq)-1):
			if gseq[gidx][0] == '-':
				left = (gseq[gidx][1:], (cidx,gidx), 't')
			else:
				left = (gseq[gidx][1:], (cidx,gidx), 'h')	
			if gseq[gidx+1][0] == '-':
				right = (gseq[gidx+1][1:], (cidx,gidx+1), 'h')
			else:
				right = (gseq[gidx+1][1:], (cidx,gidx+1), 't')		
			adj_set.add((left, right))

	return adj_set