__author__ = 'amane'

from gurobipy import *
from sys import argv
import random
from matplotlib import pyplot as plt
from collections import Counter

import RMedian_input
import RMedian_geneorder


#---------------------------------------
#A: OBTAINING AND STORING INPUT
#---------------------------------------
input_folder = argv[1]

#---------------------------------------
#1. Getting species tree
tree_file = os.path.join(input_folder,argv[2])

tree_dict = {}
rev_tree_dict = {}
tree_dict, rev_tree_dict, median = RMedian_input.sptree(tree_file, tree_dict, rev_tree_dict)
relevant_species = set()
for species in tree_dict:
	relevant_species.add(species)
for species in rev_tree_dict:
	relevant_species.add(species)
#print(relevant_species)

#---------------------------------------
#2. Getting gene content
cont_file = os.path.join(input_folder,argv[3])

cont_dict = {}
cont_dict_by_idx = {}
cont_dict, cont_dict_by_idx = RMedian_input.gene_cont(cont_file, cont_dict, cont_dict_by_idx, relevant_species)

#---------------------------------------
#3. Getting orthology relations
ortho_file = os.path.join(input_folder,argv[4])

ortho_reln = {}
rev_ortho_reln = {}

ortho_reln, rev_ortho_reln \
	= RMedian_input.ortho(ortho_file, tree_dict, rev_tree_dict, ortho_reln, rev_ortho_reln, cont_dict)

#---------------------------------------
#4. Constructing adjacency sets
#-----------------------------------------------
print("Constructing adjacency set for genomes...")
adj_set = {}
for species in relevant_species:
	adj_set[species] = RMedian_input.get_adj_set(cont_dict_by_idx[species])

root = rev_tree_dict[median]

output_folder = argv[5]

logfile1 = open(os.path.join(output_folder,'all_adjacencies.log'), "w")
logfile2 = open(os.path.join(output_folder,'observed_dups.log'), "w")


for species in adj_set:
	if species != root:						#To find instances of observed duplication
		for adj in adj_set[species]:
			logfile1.write(str(species) + "\t"+str(adj)+"\n")
			g1, g2 = adj[0][0], adj[1][0]
			idx1, idx2 = adj[0][1], adj[1][1]
			ext1, ext2 = adj[0][2], adj[1][2]
			#print(ext1, ext2)

			parent = rev_tree_dict[species]
			if ext1 != ext2:
				gpar1 = rev_ortho_reln[(species,parent)][g1]
				gpar2 = rev_ortho_reln[(species,parent)][g2]
				logfile2.write("\n")
				logfile2.write(str(species)+"\t"+str(g1)+"\t"+str(g2)+"\n")
				logfile2.write(str(parent)+"\t"+str(gpar1)+"\t"+str(gpar2)+"\n")




#-----------------------------------------------
#B: MODELLING THE ILP
#-----------------------------------------------
print("ILP formulation starts.")
m = Model("RMedian")

def par_adj_exists(adj, species, rev_ortho_reln, cont_dict, adj_set):
	g1, g2 = adj[0][0], adj[1][0]
	ext1, ext2 = adj[0][2], adj[1][2]
	
	gpar1 = rev_ortho_reln[(median, species)][g1]
	gpar2 = rev_ortho_reln[(median, species)][g2]

	idx1 = cont_dict[species][gpar1]
	idx2 = cont_dict[species][gpar2]

	par_adj = ((gpar1, idx1, ext1),(gpar2, idx2, ext2))

	if par_adj in adj_set[species] or par_adj[::-1] in adj_set[species]:
		return 1
	else:
		return 0

def ch_adj_exists(adj, tree_dict, ortho_reln, cont_dict, adj_set):
	g1, g2 = adj[0][0], adj[1][0]
	ext1, ext2 = adj[0][2], adj[1][2]

	count = {}

	for species in tree_dict[median]:
		count[species] = 0
		gch_list1 = ortho_reln[(median, species)][g1]
		gch_list2 = ortho_reln[(median, species)][g2]

		for gch1 in gch_list1:
			if gch1 in cont_dict[species]:
				idx1 = cont_dict[species][gch1]
				for gch2 in gch_list2:
					if gch2 in cont_dict[species]:
						idx2 = cont_dict[species][gch2]

						ch_adj = ((gch1, idx1, ext1),(gch2, idx2, ext2))

						if ch_adj in adj_set[species] or ch_adj[::-1] in adj_set[species]:
							count[species] = 1

	c_count = 0				
	for species in count:				
		c_count += count[species]

	return c_count		





def check_if_candidate(adj, adj_set, tree_dict, ortho_reln, rev_ortho_reln, cont_dict):
	p_count = par_adj_exists(adj, root, rev_ortho_reln, cont_dict, adj_set)
	c_count = ch_adj_exists(adj, tree_dict, ortho_reln, cont_dict, adj_set)

	common = p_count + c_count

	#if common >= 2:
	#	print(adj, common)
	return common

def assign_weights(adj, median, adj_set, tree_dict, ortho_reln, cont_dict):
	g1, g2 = adj[0][0], adj[1][0]
	ext1, ext2 = adj[0][2], adj[1][2]
	gamma = {}
	wt = 0
	#print(adj)
	for child in tree_dict[median]:
		#print(child)
		gamma[child] = 0

		l1, l2 = None, None 
		if g1 in ortho_reln[(median, child)]:
			l1 = ortho_reln[(median, child)][g1]
		if g2 in ortho_reln[(median, child)]:
			l2 = ortho_reln[(median, child)][g2]
		fam = []
		if l1 != None and l2 != None:
			for c1 in l1:
				idx1 = None
				if c1 in cont_dict[child]:
					idx1 = cont_dict[child][c1]
				for c2 in l2:
					idx2 = None
					if c2 in cont_dict[child]:
						idx2 = cont_dict[child][c2]

					if idx1 != None and idx2 != None:	
						gene1 = (c1, idx1, ext1)
						gene2 = (c2, idx2, ext2)
						new_adj = (gene1, gene2)

						fam.append(new_adj)
		for member in fam:
			if member in adj_set[child] or member[::-1] in adj_set[child]:
				gamma[child] = 1

		wt += 2*gamma[child]
	k = len(gamma)
	#print(k)
	wt -= (k+1)

	return wt









#-----------------------------------------------
#1. Creating variables
median_genes = []
for gene in cont_dict[median]:
	median_genes.append(gene)
random.shuffle(median_genes)

x = {}
w = {}
occurence = {}
count = 0
candidates = 0
total_twos, total_threes = 0, 0
for g1 in median_genes:
	idx1 = cont_dict[median][g1]
	for g2 in median_genes:
		idx2 = cont_dict[median][g2]

		#USE CANDIDATE ADJACENCY FUNCTION HERE
		adj1 = ((g1, idx1, 'h'), (g2, idx2, 't'))
		common = check_if_candidate(adj1, adj_set, tree_dict, ortho_reln, rev_ortho_reln, cont_dict)

		#occurence = check_if_candidate(adj1, adj_set, tree_dict, ortho_reln, rev_ortho_reln, cont_dict)
		if common >= 2:
			if adj1 not in x and adj1[::-1] not in x:
				x[adj1] = m.addVar(vtype=GRB.BINARY, name='x_'+str(g1)+'_h'+str(g2)+'_h')
				w[adj1] = assign_weights(adj1, median, adj_set, tree_dict, ortho_reln, cont_dict)
				occurence[adj1] = common
				if common == 2:
					total_twos += 1
				else:
					total_threes += 1
				candidates += 1
		count += 1		
		
		adj2 = ((g1, idx1, 't'), (g2, idx2, 'h'))
		common = check_if_candidate(adj2, adj_set, tree_dict, ortho_reln, rev_ortho_reln, cont_dict)
		#occurence = check_if_candidate(adj2, adj_set, tree_dict, ortho_reln, rev_ortho_reln, cont_dict)
		if common >= 2:		
			if adj2 not in x and adj2[::-1] not in x:	
				x[adj2] = m.addVar(vtype=GRB.BINARY, name='x_'+str(g1)+'_t'+str(g2)+'_t')
				w[adj2] = assign_weights(adj2, median, adj_set, tree_dict, ortho_reln, cont_dict)
				occurence[adj2] = common
				if common == 2:
					total_twos += 1
				else:
					total_threes += 1				
				candidates += 1
		count += 1	
	
		if g1 != g2:
			adj3 = ((g1, idx1, 'h'), (g2, idx2, 'h'))			
			common = check_if_candidate(adj3, adj_set, tree_dict, ortho_reln, rev_ortho_reln, cont_dict)
			#occurence = check_if_candidate(adj3, adj_set, tree_dict, ortho_reln, rev_ortho_reln, cont_dict)
			if common >= 2:
				if adj3 not in x and adj3[::-1] not in x:
					x[adj3] = m.addVar(vtype=GRB.BINARY, name='x_'+str(g1)+'_h'+str(g2)+'_h')
					w[adj3] = assign_weights(adj3, median, adj_set, tree_dict, ortho_reln, cont_dict)
					occurence[adj3] = common
					if common == 2:
						total_twos += 1
					else:
						total_threes += 1					
					candidates += 1
			count += 1		

			adj4 = ((g1, idx1, 't'), (g2, idx2, 't'))
			common = check_if_candidate(adj4, adj_set, tree_dict, ortho_reln, rev_ortho_reln, cont_dict)
			#occurence = check_if_candidate(adj4, adj_set, tree_dict, ortho_reln, rev_ortho_reln, cont_dict)
			if common >= 2:
				if adj4 not in x and adj4[::-1] not in x:	
					x[adj4] = m.addVar(vtype=GRB.BINARY, name='x_'+str(g1)+'_t'+str(g2)+'_t')
					w[adj4] = assign_weights(adj4, median, adj_set, tree_dict, ortho_reln, cont_dict)
					occurence[adj4] = common
					if common == 2:
						total_twos += 1
					else:
						total_threes += 1					
					candidates += 1
			count += 1
#print(candidates)
#print(count)
root = rev_tree_dict[median]
c = {}
for adj in adj_set[root]:
	g1, ext1 = adj[0][0], adj[0][2]
	g2, ext2 = adj[1][0], adj[1][2]
	c[adj] = m.addVar(vtype=GRB.BINARY, name='c_'+str(g1)+'_'+str(ext1)+str(g2)+'_'+str(ext2))

lambda_dict = {}
alpha_dict = {}
n_observed = {}
for branch in rev_ortho_reln:
	branch = branch[::-1]
	lambda_dict[branch] = {}
	alpha_dict[branch] = {}
	n_observed[branch] = {}
	parent = branch[0]
	for gene in cont_dict[parent]:
		n_observed[branch][gene] = m.addVar(vtype= GRB.INTEGER, name='n_observed_dups'+str(parent)+str(species)+'_'+str(gene))
		alpha_dict[branch][gene] = m.addVar(vtype= GRB.BINARY, name='alpha'+str(parent)+str(species)+'_'+str(gene))
		lambda_dict[branch][gene] = {}
		lambda_dict[branch][gene]['ceil'] = m.addVar(vtype= GRB.BINARY, name='lambda'+'_'+str(gene)+str(parent)+str(species)+'_ceil')
		lambda_dict[branch][gene]['floor'] = m.addVar(vtype= GRB.BINARY, name='lambda'+'_'+str(gene)+str(parent)+str(species)+'_floor')
				




#-----------------------------------------------
#2. Defining objective function
expr = LinExpr()

for adj in x:
	expr.addTerms(w[adj], x[adj])

for adj in c:
	expr.addTerms(2, c[adj])

for branch in alpha_dict:	
	parent, species = branch[0], branch[1]
	for gene in alpha_dict[(parent, species)]:
		expr.addTerms(-2, alpha_dict[(parent, species)][gene])	

	for gene in n_observed[(parent, species)]:
		expr.addTerms(2, n_observed[(parent, species)][gene])
m.setObjective(expr, GRB.MAXIMIZE)

#-----------------------------------------------
#3. Defining constraints
#-----------------------------------------------
#For genome consistency constraints, for every extremity, we maintain a dictionary: extr_dict
#Key: extremity_name, Value: Adj_list, which contains a list of adjacencies of which the extremity is a part
#							 Species, to which the gene containing the extremity belongs
def update_extr_dict(adj, extr_dict):
	extr1, extr2 = adj[0], adj[1]
	if extr1 not in extr_dict:
		extr_dict[extr1] = {}
		extr_dict[extr1]['Adj_List'] = [adj]
	else:
		extr_dict[extr1]['Adj_List'].append(adj)
	if extr2 not in extr_dict:
		extr_dict[extr2] = {}
		extr_dict[extr2]['Adj_List'] = [adj]
	else:
		extr_dict[extr2]['Adj_List'].append(adj)
	return extr_dict

extr_dict = {}

for adj in x:
	extr_dict = update_extr_dict(adj, extr_dict)

for extr in extr_dict:		#For every extremity, the expression ext_expr
	extr_expr = LinExpr()	#Sums the pva variables for every adjacency containing a particular extremity
	for adj in extr_dict[extr]['Adj_List']:
		extr_expr.addTerms(1, x[adj])

	m.addConstr(extr_expr <= 1, "consistency")


def adj_fam(adj, ortho_reln, cont_dict):
	g1, g2 = adj[0][0], adj[1][0]
	ext1, ext2 = adj[0][2], adj[1][2]
	l1, l2 = None, None 
	if g1 in ortho_reln[(root,median)]:		
		l1 = ortho_reln[(root,median)][g1]
	if g2 in ortho_reln[(root,median)]:	
		l2 = ortho_reln[(root,median)][g2]
	fam = []
	if l1 != None and l2 != None:
		for c1 in l1:
			idx1 = None
			if c1 in cont_dict[median]:
				idx1 = cont_dict[median][c1]
			for c2 in l2:
				idx2 = None
				if c2 in cont_dict[median]:
					idx2 = cont_dict[median][c2]
				if idx1 != None and idx2 != None:		
					gene1 = (c1, idx1, ext1)
					gene2 = (c2, idx2, ext2)
					new_adj = (gene1, gene2)
					fam.append(new_adj)	
	return fam

for adj in adj_set[root]:
	#print(adj)
	c_expr = LinExpr()
	fam = adj_fam(adj, ortho_reln, cont_dict)
	
	for member in fam:
		#print(member)
		if member in x:
			c_expr.addTerms(1, x[member])
		elif member[::-1] in x:
			c_expr.addTerms(1, x[member[::-1]])	

	if len(fam) > 0:		
		c_expr = c_expr/len(fam)
	else:
		c_expr = 0

	m.addConstr(c_expr <= c[adj], "lbd_color")
	m.addConstr(c[adj] <= c_expr + 0.99, "ubd_color")

'''
'''
observed_dups = {}
for species in rev_tree_dict:
	parent = rev_tree_dict[species]
	if parent:
		branch = (parent, species)
		observed_dups[branch] = {}
		if species == median:
			for adj in x:
				gene1, gene2 = adj[0][0], adj[1][0]	#Getting adj info
				ext1, ext2 = adj[0][2], adj[1][2]
				if ext1 != ext2:					#One ext is head, other tail	
					par1, par2 = None, None
					if gene1 in rev_ortho_reln[(species, parent)]:		#Getting parents
						par1 = rev_ortho_reln[(species, parent)][gene1]
					if gene2 in rev_ortho_reln[(species, parent)]:
						par2 = rev_ortho_reln[(species, parent)][gene2]
					if par1 == par2 and par1 != None:					#Checking if both genes have same parents
						if par1 in cont_dict[parent]:
							if par1 not in observed_dups[branch]:
								observed_dups[branch][par1] = [adj]
							else:
								observed_dups[branch][par1].append(adj)
		else:
			for adj in adj_set[species]:
				gene1, gene2 = adj[0][0], adj[1][0]	#Getting adj info
				ext1, ext2 = adj[0][2], adj[1][2]
				if ext1 != ext2:					#One ext is head, other tail	
					par1, par2 = None, None
					if gene1 in rev_ortho_reln[(species, parent)]:		#Getting parents
						par1 = rev_ortho_reln[(species, parent)][gene1]
					if gene2 in rev_ortho_reln[(species, parent)]:
						par2 = rev_ortho_reln[(species, parent)][gene2]
					if par1 == par2 and par1 != None:					#Checking if both genes have same parents
						if par1 in cont_dict[parent]:
							if par1 not in observed_dups[branch]:
								observed_dups[branch][par1] = [adj]
							else:
								observed_dups[branch][par1].append(adj)				


for branch in n_observed:
	parent, child = branch[0], branch[1]
	for gene in cont_dict[parent]:
		expr = LinExpr()
		copy_no = len(ortho_reln[branch][gene])
		expr.addConstant(0)
		if gene in observed_dups[branch]:
			for adj in observed_dups[branch][gene]:
				if child == median:
					expr.addTerms(1,x[adj])
				else:
					expr.addConstant(1)						

		if copy_no > 0:	
			expr_ubd = LinExpr()
			expr1 = expr/copy_no
			m.addConstr(lambda_dict[branch][gene]['floor'] >= expr1 - 0.99, "floor_lower_bound")
			m.addConstr(lambda_dict[branch][gene]['floor'] <= expr1, "floor_upper_bound")
			m.addConstr(lambda_dict[branch][gene]['ceil'] >= expr1, "ceil_upper_bound")
			m.addConstr(lambda_dict[branch][gene]['ceil'] <= expr1 + 0.99, "ceil_upper_bound")
		else:
			m.addConstr(lambda_dict[branch][gene]['floor'] == 0, "no_g_chr" )
			m.addConstr(lambda_dict[branch][gene]['ceil'] == 0, "no_g_chr" )

		p_idx = cont_dict[parent][gene]
		p_adj = ((gene, p_idx, 'h'),(gene, p_idx, 't'))

		if parent == root:
			if p_adj in adj_set[root]:
				m.addConstr(alpha_dict[branch][gene] <= 1, "alpha"+str(gene))
			elif p_adj[::-1] in adj_set[root]:
				m.addConstr(alpha_dict[branch][gene] <= 1, "alpha"+str(gene))
			else:
				m.addConstr(alpha_dict[branch][gene] <= 0, "alpha"+str(gene))
		else:			
			if p_adj in x:	
				m.addConstr(alpha_dict[branch][gene] <= x[p_adj], "alpha"+str(gene))
			elif p_adj[::-1] in x:
				m.addConstr(alpha_dict[branch][gene] <= x[p_adj[::-1]], "alpha"+str(gene))
			else:
				m.addConstr(alpha_dict[branch][gene] <= 0, "alpha"+str(gene))
		m.addConstr(alpha_dict[branch][gene] <= lambda_dict[branch][gene]['ceil'] - lambda_dict[branch][gene]['floor'], "alpha"+str(gene))

		expr.addTerms(-1, lambda_dict[branch][gene]['floor'])
		m.addConstr(n_observed[branch][gene] == expr, "observed_dups")		
'''

'''

# Limit how many solutions to collect
m.setParam(GRB.Param.PoolSolutions, 100000)
# Limit the search space by setting a gap for the worst possible solution that will be accepted
m.setParam(GRB.Param.PoolGap, 0)
# do a systematic search for the k- best solutions
m.setParam(GRB.Param.PoolSearchMode, 2)



#-----------------------------------------------
#4. Optimizing
m.optimize()

#-----------------------------------------------
#4. Output and Evaluation
output_file = open(os.path.join(output_folder,'final_solution.out'), "w")
genome_file = open(os.path.join(output_folder,'all_solutions.log'), "w")
stats_file = open(os.path.join(output_folder,'stats_only.txt'), "w") 



soln_dict = {}
nSolutions = m.SolCount
print('Number of solutions found : ', nSolutions)
stats_file.write("Multiplicity\t"+str(nSolutions)+"\n")
for i in range(nSolutions):
	soln_dict[i] = {}
	m.setParam(GRB.Param.SolutionNumber, i)
	#print ('Selected elements in fourth best solution :')
	#print ('\t')
	for e in x:
		soln_dict[i][e] = x[e].Xn
		#print (x[e].Xn)

#Max diff between two solutions (Diameter)
max_diff = [None, None, 0]
for i in soln_dict:
	for j in soln_dict:
		diff = 0
		for e in x:
			if soln_dict[i][e] > soln_dict[j][e]:
				diff += (soln_dict[i][e] - soln_dict[j][e])
			else:
				diff += (soln_dict[j][e] - soln_dict[i][e])
		#print(i, j, diff)
		if diff >= max_diff[2]:
			max_diff[0] = i
			max_diff[1] = j
			max_diff[2] = diff
#print(' ')			
print(max_diff)

solution_set = set()
soln_ext_dict = {}
all_ext_dict = {}
unchosen = set()
for adj in x:
	if x[adj].x > 0:
		solution_set.add(adj)
		soln_ext_dict[adj[0]] = adj[1]
		soln_ext_dict[adj[1]] = adj[0]
	else:
		unchosen.add(adj)
	if adj[0] not in all_ext_dict:	
		all_ext_dict[adj[0]] = [adj[1]]
	else:
		all_ext_dict[adj[0]].append(adj[1])
	if adj[1] not in all_ext_dict:	
		all_ext_dict[adj[1]] = [adj[0]]
	else:
		all_ext_dict[adj[1]].append(adj[0])

output_file.write("Unchosen adjacencies in first solution: \n")
for adj in unchosen:
	#print("\n")
	#print(adj)	

	output_file.write(str(adj)+ "\n")

	l, r = adj[0], adj[1]
	for ext in all_ext_dict[l]:
		l_adj = (l, ext)
		if l_adj != adj and l_adj[::-1] != adj:
			if l_adj in occurence:
				print(l_adj, occurence[l_adj])
			elif l_adj[::-1] in occurence:
				print(l_adj[::-1], occurence[l_adj[::-1]])
	for ext in all_ext_dict[r]:
		r_adj = (r, ext)
		if r_adj != adj and r_adj[::-1] != adj:			
			if r_adj in occurence:
				print(r_adj, occurence[r_adj])
			elif r_adj[::-1] in occurence:
				print(r_adj[::-1], occurence[r_adj[::-1]])

print(len(unchosen), " candidate adjacencies not chosen.")			
soln_chr_dict = {}
soln_chr_dict = RMedian_geneorder.get_gene_order(soln_chr_dict, soln_ext_dict, solution_set, x)

for chromosome in soln_chr_dict:
	#genome_file.write(soln_chr_dict[chromosome]['Chr'])
	for gene in soln_chr_dict[chromosome]['Chr']:
		genome_file.write(str(gene)+"\n")
	#if soln_chr_dict[chromosome]['Type'] == 'C':
	#	genome_file.write(soln_chr_dict[chromosome]['Chr'][0])
	genome_file.write("\n")




#Diff against actual median
output_file.write("\nComparison with actual median:\n")
output_file.write("Adjacencies in actual median: " + str(len(cont_dict[median])) + "\n")
stats_file.write("Actual_adjs\t"+ str(len(cont_dict[median])) + "\n")
output_file.write("Number of candidate adjacencies: " + str(candidates) + "\n\n")
stats_file.write("Candidate_adjs\t" + str(candidates) + "\n")	
output_file.write("Cuts: Adjacencies in ILP median but not in actual median\n")
output_file.write("Joins: Adjacencies in actual median but not in ILP median\n")
output_file.write("\nIndex\tCuts\tJoins\tTotal Distance\n")
cuts, joins, distance = {}, {}, {}
dist_array = []

total_cuts = 0

for i in soln_dict:
	cuts[i], joins[i], distance[i] = 0,0,0
	for adj in soln_dict[i]:
		if adj not in adj_set[median] and adj[::-1] not in adj_set[median]:
			cuts[i] += 1
			distance[i] += 1

	for adj in adj_set[median]:
		if adj not in soln_dict[i] and adj[::-1] not in soln_dict[i]:
			joins[i] += 1
			distance[i] += 1
		elif adj in soln_dict[i]:
			if soln_dict[i][adj] == 0:
				#print(adj)
				joins[i] += 1
				distance[i] += 1
		elif adj[::-1] in soln_dict[i]:
			if soln_dict[i][adj[::-1]] == 0:
				joins[i] += 1
				distance[i] += 1
	output_file.write(str(i) + "\t" + str(cuts[i]) + "\t" + str(joins[i]) + "\t"+ str(distance[i]) + "\n")
	total_cuts += cuts[i]
	dist_array.append(distance[i])
	#print('Cuts for solution number ', str(i), ': ', cuts[i])
	#print('Joins for solution number ', str(i), ': ', joins[i])
	#print('Total distance for solution number ', str(i), ': ', distance[i])
avg_cuts = total_cuts/nSolutions


closest = min(dist_array)
farthest = max(dist_array)
freq = Counter(dist_array)
#for i in xrange(closest,farthest+1):
#	print(freq[i])

plt.hist(dist_array)
plt.title("Freq dist")
plt.xlabel("Distance")
plt.ylabel("Frequency")
#plt.show()





#for branch in n_observed:
#	for gene in n_observed[branch]:
#		if n_observed[branch][gene].x > 0:
#			print(branch, gene)

'''
# Print number of solutions stored
nSolutions = m.SolCount
print('Number of solutions found : ', nSolutions)

# Print objective values of solutions
for e in range(nSolutions):
	m.setParam(GRB.Param.SolutionNumber, e)
	print(m.PoolObjVal)
	if e % 15 == 14:
		print (' ')
print(' ')
'''



#output_file.write("\nCandidate adjacencies not chosen: ")
unchosen_adjs = {}
twos = {}
threes = {}

#print("Overall candidate adjacencies common to 2 genomes: ", total_twos)
#print("Overall candidate adjacencies common to 3 genomes: ", total_threes)
output_file.write("\nOverall candidate adjacencies common to 2 genomes: " + str(total_twos))
stats_file.write("Cand_twos\t"+str(total_twos)+"\n")
output_file.write("\nOverall candidate adjacencies common to 3 genomes: " + str(total_threes) + "\n")
stats_file.write("Cand_threes\t"+str(total_threes)+"\n")


total_adjs = 0
for i in soln_dict:
	output_file.write("\nCandidate adjacencies not chosen in solution " + str(i)+ ": \n")
	unchosen_adjs[i] = []
	twos[i], threes[i] = 0, 0 
	for adj in soln_dict[i]:
		if soln_dict[i][adj] == 1:
			if occurence[adj] == 3:
				threes[i] += 1
			elif occurence[adj] == 2:
				twos[i] += 1
			else:
				print("Unexpected edge: ", adj, occurence[adj])
		else:
			unchosen_adjs[i].append(adj)
			output_file.write("\n"+str(adj))		

	output_file.write("\n\nAdjacencies common to 2 genomes: " + str(twos[i]))
	output_file.write("\nAdjacencies common to 3 genomes: " + str(threes[i]))
	total_adjs += twos[i]
	total_adjs += threes[i]

avg_adjs = total_adjs/nSolutions
stats_file.write("ILP_adjs\t"+str(avg_adjs)+"\n")
true_adjs = avg_adjs - avg_cuts
precision = float(true_adjs)/float(avg_adjs)
actual_adjs = len(cont_dict[median])
recall = float(true_adjs)/float(actual_adjs)
print("recall: ", recall)
print("true_adjs: ", true_adjs)
print("actual_adjs: ", actual_adjs)

stats_file.write("Precision\t"+str(precision)+"\n")
stats_file.write("Recall\t"+str(recall)+"\n")

#print("avg_adjs: ", avg_adjs)
#print("cuts: ", avg_cuts)
#print("ILP_adjs: ", len(cont_dict[median]))




