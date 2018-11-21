from gurobipy import *

def read_file(filename):
	string = open(filename, "r").read()
	string_list = string.split("\n")
	string_list = [line for line in string_list if line and line[0] != '#'] #Read line only if it is nonempty and not a comment.
	return string_list

#Takes a line from input file and stores in two nested dictionaries
#1. cont_dict: Key = species, Value = dictionary of genes -> Key = gene name, Value = location in genome as tuple (scaffold index, gene index)
#Data types: species and gene names are str, scaffold and gene index are int
#
#Second dict required to store gene orders in extant genomes
#2. cont_dict_by_idx: Key = species, Value = dictionary of genes -> Key = scaffold index, Value = list of genes in the scaffold in order
##Data types: species and gene names are str, scaffold index is int
#
def update_gene_content(line, cont_dict, cont_dict_by_idx, prev_scidx, prev_species, gidx):
	species, scidx, gene_name, gene_orientation = line[0], int(line[1]), line[2], line[3]
	if species not in cont_dict:
		cont_dict[species] = {}
		cont_dict_by_idx[species] = {}

	if scidx != prev_scidx or species != prev_species:	#Change of scaffold or change of species
		prev_scidx = scidx									
		prev_species = species
		gidx = 0
		
		cont_dict_by_idx[species][scidx] = {}
		cont_dict_by_idx[species][scidx]['gene_list'] = []
		cont_dict_by_idx[species][scidx]['type'] = 'L' 		#Since we are currently dealing with linear chr

	cont_dict[species][gene_name] = (scidx, gidx)
	cont_dict_by_idx[species][scidx]['gene_list'].append(gene_orientation + gene_name)	
	gidx += 1

	return cont_dict, cont_dict_by_idx, prev_scidx, prev_species, gidx

def gene_cont(cont_file, cont_dict, cont_dict_by_idx):
	string_list = read_file(cont_file)
	prev_scidx = None
	prev_species = None
	gidx = 0
	for line in string_list:
		line = line.split('\t')
		#line = line.split(' ')
		cont_dict, cont_dict_by_idx, prev_scidx, prev_species, gidx \
			= update_gene_content(line, cont_dict, cont_dict_by_idx, prev_scidx, prev_species, gidx)

	return cont_dict, cont_dict_by_idx

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
def update_ortho_reln(line, ortho_reln, rev_ortho_reln, cont_dict, logfile):
	parent, child, gpar, gch, fam = line[0], line[1], line[2], line[3], line[4]
	if ((parent, child)) not in ortho_reln:
		ortho_reln[(parent, child)] = {}
		rev_ortho_reln[(child, parent)] = {}
	if gpar not in ortho_reln[(parent, child)]:
		ortho_reln[(parent, child)][gpar] = []
	
	p_flag, c_flag = 0, 0	#Existence of parent gene and child gene in the gene content of respective species
	if parent in cont_dict:
		if gpar in cont_dict[parent]:
			p_flag = 1
	if child in cont_dict:	
		if gch in cont_dict[child]:
			c_flag = 1

	if p_flag != 1 or c_flag != 1:
		if p_flag != 1 and c_flag != 1:
			logfile.write(gpar + " not in " + parent + " & " + gch + " not in " + child + "\n")
		elif p_flag != 1:
			logfile.write(gpar + " not in " + parent + "\n")
		else:
			logfile.write(gch + " not in " + child + "\n")

	#if p_flag == 1 and gpar not in ortho_reln[(parent, child)]:
	if gpar not in ortho_reln[(parent, child)]:
		ortho_reln[(parent, child)][gpar] = []
	#if p_flag == 1 and c_flag == 1:
	ortho_reln[(parent, child)][gpar].append(gch)
	rev_ortho_reln[(child, parent)][gch] = gpar
	return ortho_reln, rev_ortho_reln

def ortho(ortho_file, tree_dict, rev_tree_dict, ortho_reln, rev_ortho_reln, cont_dict, logfile_1):
	string_list = read_file(ortho_file)
	for line in string_list:
		line = line.split("\t")
		#line = line.split(" ")
		tree_dict, rev_tree_dict = update_tree(line, tree_dict, rev_tree_dict)
		ortho_reln, rev_ortho_reln = update_ortho_reln(line, ortho_reln, rev_ortho_reln, cont_dict, logfile_1)
	return tree_dict, rev_tree_dict, ortho_reln, rev_ortho_reln



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






def get_gene_losses(gene_losses, nd_dict, ortho_reln, cont_dict):
	for branch in ortho_reln:
		gene_losses[branch] = set()
		nd_dict[branch] = 0
		parent, child = branch[0], branch[1]
		if parent in cont_dict:
			for gpar in cont_dict[parent]:
				lost = 1
				copy_no = 0
				if gpar in ortho_reln[branch]:
					for gch in ortho_reln[branch][gpar]:
						if child in cont_dict:
							if gch in cont_dict[child]:
								lost = 0
								copy_no += 1
				if lost == 1:
					gene_losses[branch].add(gpar)
				if copy_no >= 1:
					nd_dict[branch] += copy_no - 1
	return gene_losses, nd_dict

def get_gene_gains(gene_gains, ng_dict, rev_ortho_reln, cont_dict):
	for branch in rev_ortho_reln:
		gene_gains[branch] = set()
		ng_dict[branch] = 0
		child, parent = branch[0], branch[1]
		if child in cont_dict:
			for gch in cont_dict[child]:
				gained = 1
				if gch in rev_ortho_reln[branch]:
					gpar = rev_ortho_reln[branch][gch]
					if parent in cont_dict:
						if gpar in cont_dict[parent]:
							gained = 0
				if gained == 1:
					gene_gains[branch].add(gch)
					ng_dict[branch] += 1
	return gene_gains, ng_dict



#---------------------------------------------------------------------------------------------------------
def create_dicts(species, parent, pva_dict, wva_dict, C_dict):
	pva_dict[species] = {}
	wva_dict[species] = {}
	if parent != None:
		C_dict[(parent, species)] = {}
		C_dict[(species, parent)] = {}
	return pva_dict, wva_dict, C_dict	



#pva variables stored in a nested dictionary
#pva_dict: Key = Species, Val = dictionary of adjacencies -> Key = adjacency variable, Val = binary (0/1) (value will be assigned by optimizer)
#
#weights for every adjacency stored in a nested dictionary
#wva_dict: Key = species, Val = dictionary of adjacencies -> Key = adjacency variable, Val = weights (int), provided as input
def create_pva(m, pva_dict, wva_dict, adj, species, gene1, gene2, ext1, ext2, wt):
	pva_dict[species][adj] = m.addVar(vtype=GRB.BINARY, name='p_'+str(species)+'_'+str(gene1)+str(ext1)+str(gene2)+str(ext2))
	wva_dict[species][adj] = wt
	return m, pva_dict, wva_dict


#C_dict -> Key: (Species, Parent) and (Parent, Species) (For every species, two such keys are generated except for the root)
#		 Value:	Dictionary of adjacencies -> Key: Adjacency variable, Value: Binary (0/1) (value will be assigned by optimizer)	
#
#For every adjacency 'adj' in a species 'v' we create two types of C variables.
#For 'w' child species of 'v', create variable C_vw_adj and C_wv_adj
#For 'u' parent species of 'v', create variable C_uv_padj and C_vu_padj where 'padj' is the parent adjacency for 'adj'
def create_C_vw(m, C_dict, adj, species, gene1, gene2, ext1, ext2, tree_dict, ext_adj_set):

	if species not in ext_adj_set:
		for child in tree_dict[species]:
			if (species, child) in C_dict:
				if adj not in C_dict[(species, child)]:
					C_dict[(species, child)][adj] = m.addVar(vtype=GRB.BINARY, name='C_'+str(species)+str(child)+'_'+str(gene1)+str(ext1)+str(gene2)+str(ext2))
			if (child, species) in C_dict:
				if adj not in C_dict[(child, species)]:
					C_dict[(child, species)][adj] = m.addVar(vtype=GRB.INTEGER, name='C_'+str(child)+str(species)+'_'+str(gene1)+str(ext1)+str(gene2)+str(ext2))
	return m, C_dict

def create_C_uv(m, C_dict, species, gene1, gene2, ext1, ext2, logfile, cont_dict, rev_tree_dict, rev_ortho_reln):
	parent = rev_tree_dict[species] if species in rev_tree_dict else None
	if parent:
		gpar1, gpar2 = None, None
		if gene1 in rev_ortho_reln[(species, parent)]:
			gpar1 = rev_ortho_reln[(species, parent)][gene1]
		if gene2 in rev_ortho_reln[(species, parent)]:
			gpar2 = rev_ortho_reln[(species, parent)][gene2]
		pidx1, pidx2 = None, None
		#if gpar1 != None:
		if gpar1 in cont_dict[parent]:
			pidx1 = cont_dict[parent][gpar1]
		else:
			logfile.write("No parent for "+ gene1 + " in branch " + parent + " to " + species + "\n")
		#if gpar2 != None:
		if gpar2 in cont_dict[parent]:
			pidx2 = cont_dict[parent][gpar2]
		else:
			logfile.write("No parent for "+ gene2 + " in branch " + parent + " to " + species + "\n") 
		#if gpar1 != None and gpar2 != None:
		if pidx1 != None and pidx2 != None:
			padj = ((gpar1, pidx1, ext1),(gpar2, pidx2, ext2))
			if (species, parent) in C_dict:
				if padj not in C_dict[(species, parent)]:
					C_dict[(species, parent)][padj] = m.addVar(vtype=GRB.INTEGER, name='C_'+str(species)+str(parent)+'_'+str(gpar1)+str(ext1)+str(gpar2)+str(ext2))
			if (parent, species) in C_dict:
				if padj not in C_dict[(parent, species)]:
					C_dict[(parent, species)][padj] = m.addVar(vtype=GRB.BINARY, name='C_'+str(parent)+str(species)+'_'+str(gpar1)+str(ext1)+str(gpar2)+str(ext2))
	return m, C_dict




def create_vars(adjwt_file, logfile_2, m, pva_dict, C_dict, wva_dict, n_observed, alpha_dict, lambda_dict, \
	cont_dict, ext_adj_set, tree_dict, rev_tree_dict, rev_ortho_reln):
	string_list = read_file(adjwt_file)
	for line in string_list:
		line = line.split("\t")
		species, gene1, gene2, o1, o2, wt = line[0], line[1], line[2], line[3], line[4], float(line[5])

		idx1, idx2 = None, None
		if species in cont_dict:
			if gene1 in cont_dict[species]:
				idx1 = cont_dict[species][gene1]
			if gene2 in cont_dict[species]:
				idx2 = cont_dict[species][gene2]
		ext1 = 'h' if o1 == '+' else 't'
		ext2 = 't' if o2 == '+' else 'h'

		if idx1 != None and idx2 != None:
			adj = ((gene1, idx1, ext1),(gene2, idx2, ext2))
			m, pva_dict, wva_dict = create_pva(m, pva_dict, wva_dict, adj, species, gene1, gene2, ext1, ext2, wt)
			m, C_dict = create_C_vw(m, C_dict, adj, species, gene1, gene2, ext1, ext2, tree_dict, ext_adj_set)
		m, C_dict = create_C_uv(m, C_dict, species, gene1, gene2, ext1, ext2, logfile_2, cont_dict, rev_tree_dict, rev_ortho_reln)		
	i = 1
	for species in ext_adj_set:
		i += 1
		for adj in ext_adj_set[species]:
			gene1, gene2 = adj[0][0], adj[1][0]
			ext1, ext2 = adj[0][2], adj[1][2]
			if adj not in pva_dict[species] and adj[::-1] not in pva_dict[species]:
				wt = 1	#Assigned randomly, does not affect outcome.
				m, pva_dict, wva_dict = create_pva(m, pva_dict, wva_dict, adj, species, gene1, gene2, ext1, ext2, wt)
				m, C_dict = create_C_uv(m, C_dict, species, gene1, gene2, ext1, ext2, logfile_2, cont_dict, rev_tree_dict, rev_ortho_reln)

	for branch in rev_ortho_reln:
		branch = branch[::-1]
		lambda_dict[branch] = {}
		alpha_dict[branch] = {}
		n_observed[branch] = {}
		parent = branch[0]
		if parent in cont_dict:
			for gene in cont_dict[parent]:
				n_observed[branch][gene] = m.addVar(vtype= GRB.INTEGER, name='n_observed_dups'+str(parent)+str(species)+'_'+str(gene))
				alpha_dict[branch][gene] = m.addVar(vtype= GRB.BINARY, name='alpha'+str(parent)+str(species)+'_'+str(gene))
				lambda_dict[branch][gene] = {}
				lambda_dict[branch][gene]['ceil'] = m.addVar(vtype= GRB.BINARY, name='lambda'+'_'+str(gene)+str(parent)+str(species)+'_ceil')
				lambda_dict[branch][gene]['floor'] = m.addVar(vtype= GRB.BINARY, name='lambda'+'_'+str(gene)+str(parent)+str(species)+'_floor')
		n_observed[branch]['gains'] = m.addVar(vtype= GRB.INTEGER, name='n_observed_gains'+str(parent)+str(species))		
		
	return m, pva_dict, C_dict, wva_dict, n_observed, alpha_dict, lambda_dict