from __future__ import division
__author__ = 'amane'

from gurobipy import *
from sys import argv
import collections
from random import *

import time

import SPP_preprocessing
import SPP_rmcircular
import SPP_eval

full_start = time.time()

home_folder = "../../../"
exp_folder = "exp/2019-11-13__SPP_DeClone/"

mode = argv[1]
run = argv[2]
temp = argv[3]
rep = argv[4]
alpha = argv[5]

#Main program
output_folder = os.path.join(home_folder, exp_folder, "output",mode,"Run_"+run+"_"+temp+rep+"_"+alpha)

#os.makedirs(output_folder, exists_ok=True)
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

alpha = float(alpha)
#-----------------------------------------------
#Storing gene content
cont_file = os.path.join(home_folder, exp_folder,"formatted_input",mode,"Run_"+run+"_T0_"+temp+rep,"GeneOrders.txt")
cont_dict = {}
cont_dict_by_idx = {}
cont_dict, cont_dict_by_idx = SPP_preprocessing.gene_cont(cont_file, cont_dict, cont_dict_by_idx)

#print(cont_dict)
#Storing adjacency set of all species using gene orders (for precision recall statistics)
#actual_adj_set = {}
#actual_adj_set = SPP.get_actual_adj(cont_file, actual_adj_set)

#-----------------------------------------------
#Storing orthology relations
ortho_file = os.path.join(home_folder, exp_folder,"formatted_input",mode,"Run_"+run+"_T0_"+temp+rep,"Orthology.txt")

tree_dict = {}
rev_tree_dict = {}
ortho_reln = {}
rev_ortho_reln = {}

logfile_1 = open(os.path.join(output_folder,'OrthoRelns_Extra.log'),"w")	#Orthology relations where either the anc or desc gene is missing in gene content
tree_dict, rev_tree_dict, ortho_reln, rev_ortho_reln \
	= SPP_preprocessing.ortho(ortho_file, tree_dict, rev_tree_dict, ortho_reln, rev_ortho_reln, cont_dict, logfile_1)

all_parents = set()
for branch in ortho_reln:
	parent = branch[0]
	all_parents.add(parent)

#-----------------------------------------------
#Segregating gene losses and de novo gene creations. 
#Also maintaining a dict of no. of duplicates for every branch: nd_dict
#and a dict of no. of gene gains: ng_dict
#Key = (ancestor, descendant), Value = no. of duplicate genes in child
gene_losses = {}
gene_gains = {}
nd_dict = {}
ng_dict = {}
gene_losses, nd_dict = SPP_preprocessing.get_gene_losses(gene_losses, nd_dict, ortho_reln, cont_dict)
gene_gains, ng_dict = SPP_preprocessing.get_gene_gains(gene_gains, ng_dict, rev_ortho_reln, cont_dict)

total_gene_loss = 0
for x in gene_losses:
	total_gene_loss += len(gene_losses[x])
	#print(len(gene_losses[x]))
	#print(gene_losses[x])
#print(total_gene_loss)



#-----------------------------------------------
print("\nConstructing adjacency set for extant genomes...")
ext_adj_set = {}
initial_adjs = {}
initial_adjs_total = 0
for species in rev_tree_dict:
	initial_adjs[species] = 0
	if species not in tree_dict:
		if species in cont_dict_by_idx:
			ext_adj_set[species] = list(SPP_preprocessing.get_ext_adj_set(cont_dict_by_idx[species]))
			initial_adjs[species] += len(ext_adj_set[species])
		#print(ext_adj_set[species])
	initial_adjs_total += initial_adjs[species]

#for x in ext_adj_set:
#extant_count = 0
#for x in cont_dict_by_idx:
#	extant_count +=1
#	print(x)	
#print(extant_count,"\n")
print("Storing adjacency set for ancestral genomes...")
anc_adj_set = {}
actual_anc_adjs = {}
actual_anc_adjs_total = 0 
for species in tree_dict:
	actual_anc_adjs[species] = 0
	if species in cont_dict_by_idx:
		#print(species, tree_dict[species], len(cont_dict_by_idx[species][1]['gene_list']))
		anc_adj_set[species] = list(SPP_preprocessing.get_anc_adj_set(cont_dict_by_idx[species]))
		actual_anc_adjs[species] += len(anc_adj_set[species])
		#print(anc_adj_set[species])
	actual_anc_adjs_total += actual_anc_adjs[species]


	



initial_scaffolds = {}
for species in rev_tree_dict:
	if species not in tree_dict: 
		if species in cont_dict_by_idx:
			initial_scaffolds[species] = len(cont_dict_by_idx[species])	#No. of scaffolds initial




#-----------------------------------------------
#Modelling the ILP
print("ILP formulation starts...")
m = Model("SPP")

#-----------------------------------------------
#Initializing nested dictionaries for variables p_v(adj) and C_uv(adj)/C_vu(adj). 
#Every adjacency is of the type: ((gene1, idx1, 'h/t'), (gene2, idx2, 'h/t')) where idx is an int tuple (scaffold idx, gene index)
pva_dict = {}
C_dict = {}

#Every adjacency has a weight associated with it, to be stored in wva_dict
wva_dict = {}

for species in rev_tree_dict:
	parent = rev_tree_dict[species]
	pva_dict, wva_dict, C_dict = SPP_preprocessing.create_dicts(species, parent, pva_dict, wva_dict, C_dict)
for species in tree_dict:
	if species not in rev_tree_dict:
		parent = None
		pva_dict, wva_dict, C_dict = SPP_preprocessing.create_dicts(species, parent, pva_dict, wva_dict, C_dict)


print("Assigning adjacency weights according to input...")
logfile_2 = open(os.path.join(output_folder,'OrthoRelns_Missing.log'),"w")


adjwt_file = os.path.join(home_folder, exp_folder,"formatted_input",mode,"Run_"+run+"_T0_"+temp+rep,"Adjacencies.txt")
#adjwt_file = os.path.join(home_folder,"exp","2019-09-03__DeCoSTAR_test3","output",mode,"Run_"+run,"DeClone","DeClone_"+temp+rep,"formatted_adjacencies")
lambda_dict = {}
alpha_dict = {}
n_observed = {}

m, pva_dict, C_dict, wva_dict, n_observed, alpha_dict, lambda_dict \
	= SPP_preprocessing.create_vars(adjwt_file, logfile_2, m, pva_dict, C_dict, wva_dict, n_observed, alpha_dict, lambda_dict, \
		cont_dict, ext_adj_set, tree_dict, rev_tree_dict, rev_ortho_reln)
print("Reaching past create_vars")

#-----------------------------------------------
#Setting up the expression for the objective function
print("Setting up expression for objective function...")
#alpha = 0
expr = LinExpr()

for species in pva_dict:
	for adj in pva_dict[species]:
		expr.addConstant(alpha*wva_dict[species][adj])
		expr.addTerms(-alpha*wva_dict[species][adj], pva_dict[species][adj])
		
	parent = rev_tree_dict[species] if species in rev_tree_dict else None
	if parent:
		for adj in C_dict[(parent, species)]:
			expr.addTerms(1-alpha, C_dict[(parent, species)][adj])
		for adj in C_dict[(species, parent)]:
			expr.addTerms(1-alpha, C_dict[(species, parent)][adj])

		expr.addConstant(2*(1-alpha)*nd_dict[(parent, species)])

		expr.addConstant((1-alpha)*len(gene_losses[(parent, species)]))		

		expr.addConstant(2*(1-alpha)*ng_dict[(species, parent)])

		for gene in n_observed[(parent, species)]:
			expr.addTerms(-2*(1-alpha), n_observed[(parent, species)][gene])
		
		for gene in alpha_dict[(parent, species)]:
			expr.addTerms(2*(1-alpha), alpha_dict[(parent, species)][gene])



m.setObjective(expr, GRB.MINIMIZE)

#-----------------------------------------------
#Setting up constraints

#-----------------------------------------------
#c_uv_a and c_vu_a constraints 
def adj_fam(parent, species, adj, ortho_reln, cont_dict):
	l1, l2 = None, None 
	if adj[0][0] in ortho_reln[(parent,species)]:		
		l1 = ortho_reln[(parent, species)][adj[0][0]]
	if adj[1][0] in ortho_reln[(parent,species)]:	
		l2 = ortho_reln[(parent, species)][adj[1][0]]
	fam = []
	if l1 != None and l2 != None:
		for g1 in l1:
			idx1 = None
			if species in cont_dict:
				if g1 in cont_dict[species]:
					idx1 = cont_dict[species][g1]
			for g2 in l2:
				idx2 = None
				if species in cont_dict:
					if g2 in cont_dict[species]:
						idx2 = cont_dict[species][g2]
				if idx1 != None and idx2 != None:		
					gene1 = (g1, idx1, adj[0][2])
					gene2 = (g2, idx2, adj[1][2])
					new_adj = (gene1, gene2)
					fam.append(new_adj)	
	return(fam)

def assign_premul(adj, parent, species, gene_losses):
	branch = (parent, species)
	ext1, ext2 = adj[0], adj[1]
	gene1, gene2 = ext1[0], ext2[0]
	if gene1 in gene_losses[branch] and gene2 in gene_losses[branch]:
		premul = 1
	elif gene1 in gene_losses[branch] or gene2 in gene_losses[branch]:
		premul = 1
	else:
		premul = 1
	return premul


def C_expr(adj, species, parent, ortho_reln, cont_dict, i, premul):
	fam = adj_fam(parent, species, adj, ortho_reln, cont_dict)
	expr = LinExpr()
	for member in fam:
		if member in pva_dict[species]:
			expr.addTerms(-i, pva_dict[species][member])
		elif member[::-1] in pva_dict[species]:
			expr.addTerms(-i, pva_dict[species][member[::-1]])
	if adj in pva_dict[parent]:	
		expr.addTerms(i*premul, pva_dict[parent][adj])
	return expr

#-----------------------------------------------
#For genome consistency constraints, for every extremity, we maintain a dictionary: extr_dict
#Key: extremity_name, Value: Adj_list, which contains a list of adjacencies of which the extremity is a part
#							 Species, to which the gene containing the extremity belongs
def update_extr_dict(adj, species, extr_dict):
	extr1, extr2 = adj[0], adj[1]
	if extr1 not in extr_dict:
		extr_dict[extr1] = {}
		extr_dict[extr1]['Adj_List'] = [adj]
		extr_dict[extr1]['Species'] = species
	else:
		extr_dict[extr1]['Adj_List'].append(adj)
	if extr2 not in extr_dict:
		extr_dict[extr2] = {}
		extr_dict[extr2]['Adj_List'] = [adj]
		extr_dict[extr2]['Species'] = species
	else:
		extr_dict[extr2]['Adj_List'].append(adj)
	return extr_dict		
#-----------------------------------------------
print("Setting up constraints...")

nconstr = 0
#p_v_a = 1 for all known extant adjacencies 
for species in ext_adj_set:	
	for adj in ext_adj_set[species]:
		if adj in pva_dict[species]:
			m.addConstr(pva_dict[species][adj] == 1, "extant")
			nconstr += 1
		elif adj[::-1]:
			m.addConstr(pva_dict[species][adj[::-1]] == 1, "extant")
			nconstr += 1
print("\nExtant adj constraints done.")
print("No. of constraints = ", nconstr)

#Constraints accounting for the SCJ distance
for species in rev_tree_dict:																			
	parent = rev_tree_dict[species]																			
	if parent:																							
		for adj in C_dict[(parent,species)]:															
			m.addConstr(C_dict[(parent,species)][adj] >= 0, "non-negativity"+parent+"to"+species)
			premul = assign_premul(adj, parent, species, gene_losses)
			expr = C_expr(adj, species, parent, ortho_reln, cont_dict, 1, premul)
			m.addConstr(C_dict[(parent,species)][adj] >= expr, parent+"to"+species)
			nconstr += 1

		for adj in C_dict[(species,parent)]:
			m.addConstr(C_dict[(species,parent)][adj] >= 0, "non-negativity"+species+"to"+parent)
			premul = assign_premul(adj, parent, species, gene_losses)
			expr = C_expr(adj, species, parent, ortho_reln, cont_dict, -1, premul)
			m.addConstr(C_dict[(species,parent)][adj] >= expr, species+"to"+parent)
			nconstr += 1

print("\nChange constraints done.")
print("No. of constraints = ", nconstr)


#Constraints confirming the consistency of genomes
extr_dict = {}
for species in pva_dict:
	for adj in pva_dict[species]:
		extr_dict = update_extr_dict(adj, species, extr_dict)

for ext in extr_dict:		#For every extremity, the expression ext_expr
	ext_expr = LinExpr()	#Sums the pva variables for every adjacency containing a particular extremity
	species = extr_dict[ext]['Species']
	for adj in extr_dict[ext]['Adj_List']:
		ext_expr.addTerms(1,pva_dict[species][adj])

	m.addConstr(ext_expr <= 1, "consistency")
	nconstr += 1		
print("\nConsistency constraints done.")
print("No. of constraints = ", nconstr)

#Constraints to count the number of observed duplications.
#These constraints have been implemented according to the new SCJ-TD-FD formula in the median paper.
observed_dups = {}
observed_gains = {}
for species in rev_tree_dict:
	parent = rev_tree_dict[species]
	if parent:
		branch = (parent, species)
		observed_dups[branch] = {}
		observed_gains[branch] = {}
		for adj in pva_dict[species]:
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
				elif gene1 == gene2:
					observed_gains[branch] = [adj]

for branch in n_observed:
	parent, child = branch[0], branch[1]
	if parent in cont_dict:
		for gene in cont_dict[parent]:
			expr = LinExpr()
			if gene in gene_losses[branch]:
				expr.addConstant(0)
				copy_no = 0
			else:	
				copy_no = len(ortho_reln[branch][gene])
				expr.addConstant(0)
				if branch in observed_dups:
					if gene in observed_dups[branch]:
						for adj in observed_dups[branch][gene]:
							expr.addTerms(1,pva_dict[child][adj])

		
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

			p_idx = cont_dict[parent][gene]			#Finding if pva_adj = 1 in ancestor genome
			p_adj = ((gene, p_idx, 'h'),(gene, p_idx, 't'))
			if p_adj in pva_dict[parent]:	
				m.addConstr(alpha_dict[branch][gene] <= pva_dict[parent][p_adj], "alpha"+str(gene))
			elif p_adj[::-1] in pva_dict[parent]:
				m.addConstr(alpha_dict[branch][gene] <= pva_dict[parent][p_adj[::-1]], "alpha"+str(gene))
			else:
				m.addConstr(alpha_dict[branch][gene] <= 0, "alpha"+str(gene))
			m.addConstr(alpha_dict[branch][gene] <= lambda_dict[branch][gene]['ceil'] - lambda_dict[branch][gene]['floor'], "alpha"+str(gene)) 

			expr.addTerms(-1, lambda_dict[branch][gene]['floor'])
			m.addConstr(n_observed[branch][gene] == expr, "observed_dups")
	gains_expr = LinExpr()
	gains_expr.addConstant(0)
	if branch in observed_gains:
		for adj in observed_gains[branch]:
			gains_expr.addTerms(1,pva_dict[child][adj])
		m.addConstr(n_observed[branch]['gains'] == gains_expr, "observed_gains")


#m.optimize()

#m.computeIIS()
#m.write("m.ilp")


#-----------------------------------------------
#Iteratively remove circular chromosomes if any and re-optimize
logfile_3 = open(os.path.join(output_folder,'Circ_chromosomes.log'),"w")
time_file = open(os.path.join(output_folder,'Time_per_iteration.log'),"w")
circular = 1
i = 1
while circular == 1:
	start = time.time()
	m.optimize()
	stop = time.time()
	duration = stop - start
	time_file.write("Iteration "+str(i)+":\t"+str(duration)+"\n")


	soln_ext_dict = {}
	soln_adj_dict = {}
	soln_ext_dict, soln_adj_dict = SPP_rmcircular.get_adj_soln(soln_ext_dict, soln_adj_dict, pva_dict)		

	reached = {}
	soln_chr_dict = {}
	circ_chr = {}
	c = 1 
	logfile_3.write("Iteration "+ str(i)+"\n")

	#logfile_4 = open(os.path.join(output_folder,'Final_soln.log'),"w")
	circular = 0	

	circular, circ_chr, soln_chr_dict = SPP_rmcircular.check_if_circular(reached, soln_chr_dict, circ_chr, c, i, soln_ext_dict, soln_adj_dict, pva_dict, logfile_3)
	if circular == 1:	#If chr is circular, adds corresponding constraint
		if i <= 50:
			x = randint(1, len(circ_chr))
			species, chosen_chr = circ_chr[x]['Species'], circ_chr[x]['Chr']
			logfile_3.write("\nChosen seed at iteration "+str(i)+": "+str(x))
			logfile_3.write("\nChosen chr at iteration "+str(i)+"\n")
			expr = LinExpr()
			for adj in chosen_chr:
				if adj in pva_dict[species]:
					logfile_3.write(str(species)+": "+ str(adj)+"\n")
					expr.addTerms(1, pva_dict[species][adj])
				elif adj[::-1] in pva_dict[species]:
					logfile_3.write(str(species)+": "+ str(adj[::-1])+"\n")
					expr.addTerms(1, pva_dict[species][adj[::-1]])
			logfile_3.write("\n\n")		
			m.addConstr(expr <= len(chosen_chr) - 1, "circular")
					
			print("\n", i, "th circular constraint satisfied\n")
			i += 1
		else:
			break
				




		
#-----------------------------------------------
#POST-PROCESSING

logfile_4 = open(os.path.join(output_folder,'Final_soln.log'),"w")
#logfile_4.write("#Species_name\n")
#logfile_4.write("#Chr_1\tSpName@GeneFam_GeneId+ SpName@GeneFam_GeneId- SpName@GeneFam_GeneId+\n")
for species in soln_chr_dict:
	logfile_4.write("Species:\t"+species+"\n")
	gene_set = set()
	for chrom in soln_chr_dict[species]:
		logfile_4.write(str(chrom)+":\t")
		#print(chrom)
		chr_type = soln_chr_dict[species][chrom]['Type']
		for gene in soln_chr_dict[species][chrom]['Chr']:
			orient = gene[1]
			gname = gene[0][0]
			gene_set.add(gname)
			logfile_4.write(gname+orient+" ")
		logfile_4.write("\n")	
		print(species, chrom, soln_chr_dict[species][chrom])
	for gene in cont_dict[species]:
		if gene not in gene_set:
			logfile_4.write("single-gene:\t"+gene+"+ \n")	
	logfile_4.write("\n")

final_log = open(os.path.join(output_folder,'stats_'+str(alpha)+'.out'),"w")

final_log.write("Total number of circular chromosomes removed: "+str(i-1)+"\n\n")

#Improvements in extant genomes: Gained adjacencies and Final no. of scaffolds
gained_adjs = {}
gained_adjs_total = 0
final_adjs = {}
final_adjs, gained_adjs, gained_adjs_total = SPP_eval.adj_gains(ext_adj_set, initial_adjs, final_adjs, gained_adjs, gained_adjs_total, pva_dict, final_log)
#print(species, initial_adjs[species], gained_adjs[species])


final_scaffolds = {}
anc_scaffolds = {}
final_scaffolds, anc_scaffolds = SPP_eval.count_scaffolds(ext_adj_set, tree_dict, initial_scaffolds, final_scaffolds, anc_scaffolds, soln_chr_dict, final_log)

SPP_eval.selected_adjacencies(pva_dict, wva_dict, final_log)


SPP_eval.compute_stats(pva_dict, anc_adj_set, ext_adj_set, actual_anc_adjs_total, final_log)


cuts_and_joins = open(os.path.join(output_folder,'Cuts_and_joins.log'),"w")
	
SPP_eval.get_distance_details(pva_dict, rev_tree_dict, C_dict, ng_dict, nd_dict, n_observed, gene_gains, gene_losses, final_log, cuts_and_joins)




full_stop = time.time()
total_duration = full_stop - full_start
final_log.write("\n\nTime taken = "+str(total_duration))


