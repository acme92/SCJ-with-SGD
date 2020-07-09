from gurobipy import *
import collections

def get_adj_soln(soln_ext_dict, soln_adj_dict, pva_dict):
	for species in pva_dict:	#Separating those with pva = 1
		soln_ext_dict[species] = {}
		soln_adj_dict[species] = []
		for adj in pva_dict[species]:
			var = pva_dict[species][adj]
			if var.x > 0:
				soln_adj_dict[species].append(adj) 
				soln_ext_dict[species][adj[0]] = adj[1]
				soln_ext_dict[species][adj[1]] = adj[0]
	return soln_ext_dict, soln_adj_dict

def other(ext):
	c = 't' if ext == 'h' else 'h'
	return c

def get_next_adj(gene, ext, ext_dict):
	next_gene, next_ext, next_adj = None, None, None
	other_ext = other(ext)[0]
	l = (gene[0], gene[1], other_ext)
	if l in ext_dict:
		next_gene = ext_dict[l][0:2]
		next_ext = ext_dict[l][2]
	if next_gene != None:
		r = (next_gene[0], next_gene[1], next_ext)
		next_adj = (l,r)
	return next_adj, next_gene, next_ext

def check_if_circular(reached, soln_chr_dict, circ_chr, c, i, soln_ext_dict, soln_adj_dict, pva_dict, logfile_3):
	circular = 0
	for species in soln_adj_dict:
		reached[species] = {}
		soln_chr_dict[species] = {}
		chr_no = 0

		for adj in soln_adj_dict[species]:
			#circular = 0
			l, r = adj[0], adj[1]
			g1, ext1, g2, ext2 = l[0:2], l[2], r[0:2], r[2]
			if g1 not in reached[species] and g2 not in reached[species]:		#Every gene would be a part of exactly one chr
				chr_no += 1
				soln_chr_dict[species][chr_no] = {}
				d = collections.deque()			#List (queue) of (gene, orientation) tuples
				curr_chr = collections.deque()	#List (queue) of adjacecencies

				reached[species][g1] = 1
				reached[species][g2] = 1
		
				d.append((g1, '+')) if ext1 == 'h' else d.append((g1, '-'))
				d.append((g2, '+')) if ext2 == 't' else d.append((g2, '-'))
					
				curr_chr.append(adj)

				telomere = 0

				#Follow g2 dir
				gene, ext = g2, ext2
				while telomere == 0:
					next_adj = None
					next_adj, next_gene, next_ext = get_next_adj(gene, ext, soln_ext_dict[species])

					if next_gene != None and next_gene in reached[species]:
						circular = 1
						chr_type = 'C'
						curr_chr.append(next_adj)

						circ_chr[c] = {}
						circ_chr[c]['Species'] = species
						circ_chr[c]['Chr'] = curr_chr
						logfile_3.write(str(i)+"\t"+str(curr_chr)+"\n")
						c += 1						

						break	#As circular chrom reached. Breaks from while.

					elif next_gene != None:
						if next_ext == 't':
							d.append((next_gene, '+'))
						else:
							d.append((next_gene, '-'))
						curr_chr.append(next_adj)

						reached[species][next_gene] = 1
						gene = next_gene
						ext = next_ext	

					else:
						telomere = 1	#As telomere reached

				#Follow g1 dirn if not circular					
				if telomere == 1:
					chr_type = 'L'
					#telomere == 0	
					gene, ext = g1, ext1
					while telomere == 1:
						next_adj = None
						next_adj, next_gene, next_ext = get_next_adj(gene, ext, soln_ext_dict[species])

						if next_gene != None:
							if next_ext == 'h':
								d.appendleft((next_gene, '+'))
							else:
								d.appendleft((next_gene, '-'))
							reached[species][next_gene] = 1
							gene = next_gene
							ext = next_ext

						else:
							telomere = 2

				soln_chr_dict[species][chr_no]['Chr'] = d
				soln_chr_dict[species][chr_no]['Type'] = chr_type
				#for adj in d:
				#	if adj in pva_dict[species]:
				#		logfile_4.write(str(species)+": "+ str(adj)+"\n")
				#	elif adj[::-1] in pva_dict[species]:
				#		logfile_4.write(str(species)+": "+ str(adj[::-1])+"\n")	

	#logfile_4.close()
	return circular, circ_chr, soln_chr_dict		