from __future__ import division

def adj_gains(ext_adj_set, initial_adjs, final_adjs, gained_adjs, gained_adjs_total, pva_dict, final_log):
	initial_adjs_total = 0
	for species in ext_adj_set:
		final_adjs[species] = 0
		for adj in pva_dict[species]:
			if pva_dict[species][adj].x > 0:
				final_adjs[species] += 1
		gained_adjs[species] = final_adjs[species] - initial_adjs[species]
		initial_adjs_total += initial_adjs[species]
		gained_adjs_total += gained_adjs[species]

	final_log.write("Initial adjacencies in extant genomes: "+ str(initial_adjs_total)+"\n")
	final_log.write("Gained adjacencies in extant genomes: "+ str(gained_adjs_total)+"\n")		
	return final_adjs, gained_adjs, gained_adjs_total	

def count_scaffolds(ext_adj_set, tree_dict, initial_scaffolds,  final_scaffolds, anc_scaffolds, soln_chr_dict, final_log):
	for species in ext_adj_set:
		final_scaffolds[species] = len(soln_chr_dict[species])
		#print(species, initial_scaffolds[species], gained_adjs[species], final_scaffolds[species])
		#print("This is the soln chr dict:", species, soln_chr_dict[species])	
	for species in tree_dict:
		if species not in ext_adj_set:
			anc_scaffolds[species] = len(soln_chr_dict[species])
			#print("This is the soln chr dict:", species, soln_chr_dict[species])	
	initial_scaffolds_total = 0
	final_scaffolds_total = 0
	for species in ext_adj_set:
		initial_scaffolds_total += initial_scaffolds[species]
		final_scaffolds_total += final_scaffolds[species]
	#print("Total:", initial_scaffolds_total, gained_adjs_total, final_scaffolds_total)

	final_log.write("Initial scaffolds in extant genomes: "+ str(initial_scaffolds_total)+"\n")
	final_log.write("Final scaffolds in extant genomes: "+ str(final_scaffolds_total)+"\n")

	final_log.write("\nSpecies\tNo. of scaffolds\n")
	anc_scaffolds = {}
	for species in tree_dict:
		if species not in ext_adj_set:
			anc_scaffolds[species] = len(soln_chr_dict[species])
			final_log.write(str(species)+"\t"+str(len(soln_chr_dict[species]))+"\n")

	for species in ext_adj_set:
		final_log.write(str(species)+"\t"+str(len(soln_chr_dict[species]))+"\n")			
	return final_scaffolds, anc_scaffolds

def selected_adjacencies(pva_dict, wva_dict, final_log):
	#Adjacencies selected from input
	total_adjs = {}
	total_adjwts = {}
	selected_adjs = {}
	selected_adjwts = {}
	all_adjs_total, all_adjs_selected = 0, 0
	all_adjwts_total, all_adjwts_selected = 0, 0
	for species in pva_dict:
		total_adjs[species] = 0
		total_adjwts[species] = 0
		selected_adjs[species] = 0
		selected_adjwts[species] = 0
		for adj in pva_dict[species]:
			total_adjs[species] += 1
			total_adjwts[species] += wva_dict[species][adj]
			if pva_dict[species][adj].x > 0:
				selected_adjs[species] += 1
				selected_adjwts[species] += wva_dict[species][adj]
		all_adjs_total += total_adjs[species]
		all_adjwts_total += total_adjwts[species]
		all_adjs_selected += selected_adjs[species]
		all_adjwts_selected += selected_adjwts[species]		

	final_log.write("\n\n\nTotal Adj Wts: "+str(all_adjwts_total)+"\n")
	final_log.write("Selected Adj Wts: "+str(all_adjwts_selected)+"\n")
	final_log.write("Total Adjs: "+str(all_adjs_total)+"\n")
	final_log.write("Selected Adjs: "+str(all_adjs_selected)+"\n")	

def compute_stats(pva_dict, anc_adj_set, ext_adj_set, actual_anc_adjs_total, final_log):
	#Listing down selected adjacencies
	selected_adj_set = {}
	TP_total = set()
	FP_total = set()
	TP = {}
	FP = {}
	for species in pva_dict:
		selected_adj_set[species] = set()
		TP[species] = set()
		FP[species] = set()
		for adj in pva_dict[species]:
			if pva_dict[species][adj].x > 0:
				#c1+=1
				selected_adj_set[species].add(adj)
				if species in anc_adj_set:
					
					if adj in set(anc_adj_set[species]) or adj[::-1] in set(anc_adj_set[species]):
						TP_total.add(adj)
						TP[species].add(adj)
					else:	
						FP_total.add(adj)
						FP[species].add(adj)
				else:
					
					if adj in set(ext_adj_set[species]) or adj[::-1] in set(ext_adj_set[species]):
				#		TP_total.add(adj)
						TP[species].add(adj)
					else:	
				#		FP_total.add(adj)
						FP[species].add(adj)
	anc_precision = len(TP_total)/(len(FP_total)+len(TP_total))

	anc_recall = len(TP_total)/actual_anc_adjs_total

	final_log.write("\n\n\nPrecision: "+str(anc_precision)+"\n")
	final_log.write("Recall: "+str(anc_recall)+"\n")
	F1_score = 2*anc_precision*anc_recall/(anc_precision+anc_recall)
	final_log.write("F1_score: "+str(F1_score)+"\n")	


def list_cuts_and_joins(parent, species, distance, C_dict, gene_gains, gene_losses, cuts_and_joins):
	cuts_and_joins.write("\nFrom "+parent+ " to "+ species+"\n")
	cuts_and_joins.write("Cuts:\n")
	for adj in C_dict[(parent, species)]:
		distance[(parent, species)]['Cuts'] += C_dict[(parent, species)][adj].x
		distance[(parent, species)]['SCJTDFD'] += C_dict[(parent, species)][adj].x

				
		if C_dict[(parent, species)][adj].x == 1:
			cuts_and_joins.write(str(adj)+"\n")

		gene1, gene2 = adj[0][0], adj[1][0]
		if gene1 in gene_losses[(parent, species)] or gene2 in gene_losses[(parent, species)]:
			distance[(parent, species)]['Cuts_by_Losses'] += C_dict[(parent, species)][adj].x

	cuts_and_joins.write("\nJoins:\n")		
	for adj in C_dict[(species, parent)]:
		distance[(parent, species)]['Joins'] += C_dict[(species, parent)][adj].x
		distance[(parent, species)]['SCJTDFD'] += C_dict[(species, parent)][adj].x

				
		if C_dict[(species, parent)][adj].x == 1:
			cuts_and_joins.write(str(adj)+"\n")			

		gene1, gene2 = adj[0][0], adj[1][0]
		if gene1 in gene_gains[(species, parent)] or gene2 in gene_gains[(species, parent)]:
			distance[(parent, species)]['Joins_by_Gains'] += C_dict[(species, parent)][adj].x
	cuts_and_joins.write("\n\n")

	return distance	

def get_distance_details(pva_dict, rev_tree_dict, C_dict, ng_dict, nd_dict, n_observed, gene_gains, gene_losses, final_log, cuts_and_joins):	
	#Distance between branches
	distance = {}
	for species in pva_dict:
		parent = rev_tree_dict[species] if species in rev_tree_dict else None
		distance[(parent, species)] = {}
		distance[(parent, species)]['Cuts'] = 0
		distance[(parent, species)]['Joins'] = 0
		distance[(parent, species)]['Cuts_by_Losses'] = 0
		distance[(parent, species)]['Dups'] = 0
		distance[(parent, species)]['Gains'] = 0
		distance[(parent, species)]['Joins_by_Gains'] = 0
		distance[(parent, species)]['Observed_Dups'] = 0
		distance[(parent, species)]['SCJTDFD'] = 0


		if parent:
			distance = list_cuts_and_joins(parent, species, distance, C_dict, gene_gains, gene_losses, cuts_and_joins)

			distance[(parent, species)]['Dups'] += nd_dict[(parent, species)]	
			distance[(parent, species)]['SCJTDFD'] += 2*nd_dict[(parent, species)]

			distance[(parent, species)]['Gains'] += ng_dict[(species, parent)]
			distance[(parent, species)]['SCJTDFD'] += 2*ng_dict[(species, parent)]

			if (parent, species) in n_observed:
				for gene in n_observed[(parent, species)]:
					distance[(parent, species)]['Observed_Dups'] += n_observed[(parent, species)][gene].x
					distance[(parent, species)]['SCJTDFD'] += -2*n_observed[(parent, species)][gene].x

	final_log.write("\n\n\nBranch\tSCJTDFD\tCuts\tJoins\tDuplications\tObserved Dups\n")
	total_dist, total_cuts, total_joins, cuts_by_losses, joins_by_gains = 0, 0, 0, 0, 0
	#x = branch in this case
	for x in distance:
		final_log.write("("+str(x[0])+","+str(x[1])+")"+"\t"+str(distance[x]['SCJTDFD'])+"\t"+str(distance[x]['Cuts'])+"\t"+str(distance[x]['Joins'])+"\t"+str(distance[x]['Dups'])+"\t"+str(distance[x]['Observed_Dups'])+"\n")
		total_dist += distance[x]['SCJTDFD']
		total_cuts += distance[x]['Cuts']
		total_joins += distance[x]['Joins']
		cuts_by_losses += distance[x]['Cuts_by_Losses']
		joins_by_gains += distance[x]['Joins_by_Gains']
	final_log.write("Overall\t"+str(total_dist)+"\t"+str(total_cuts)+"\t"+str(total_joins)+"\n")	
	final_log.write("Cuts involving lost genes\t"+"\t"+str(cuts_by_losses)+"\n")	
	final_log.write("Joins involving gained genes\t"+"\t"+str(joins_by_gains))
