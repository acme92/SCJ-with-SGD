from gurobipy import *
import collections

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


def get_gene_order(soln_chr_dict, soln_ext_dict, solution_set, x):
	reached = {}
	soln_chr_dict = {}
	chr_no = 0

	for adj in solution_set:
		l, r = adj[0], adj[1]
		g1, ext1, g2, ext2 = l[0:2], l[2], r[0:2], r[2]

		if g1 not in reached and g2 not in reached:								#Every gene would be a part of exactly one chr
			chr_no += 1
			soln_chr_dict[chr_no] = {}
			d = collections.deque()
			curr_chr = collections.deque()

			reached[g1] = 1
			reached[g2] = 1

			d.append((g1, '+')) if ext1 == 'h' else d.append((g1, '-'))
			d.append((g2, '+')) if ext2 == 't' else d.append((g2, '-'))

			curr_chr.append(adj)

			telomere = 0

			#Follow g2 dir
			gene, ext = g2, ext2
			while telomere == 0:
				next_adj = None
				next_adj, next_gene, next_ext = get_next_adj(gene, ext, soln_ext_dict)

				if next_gene != None and next_gene in reached:
					chr_type = 'C'
					curr_chr.append(next_adj)

				elif next_gene != None:
					if next_ext == 't':
						d.append((next_gene, '+'))
					else:
						d.append((next_gene, '-'))
					curr_chr.append(next_adj)

					reached[next_gene] = 1
					gene = next_gene
					ext = next_ext

				else:
					telomere = 1	#As telomere reached	

			#Follow g1 dirn if not circular					
			if telomere == 1:
				chr_type = 'L'
				gene, ext = g1, ext1
				while telomere == 1:
					next_adj = None
					next_adj, next_gene, next_ext = get_next_adj(gene, ext, soln_ext_dict)

					if next_gene != None:
						if next_ext == 'h':
							d.appendleft((next_gene, '+'))
						else:
							d.appendleft((next_gene, '-'))

						reached[next_gene] = 1
						gene = next_gene
						ext = next_ext

					else:
						telomere = 2

			soln_chr_dict[chr_no]['Chr'] = d
			soln_chr_dict[chr_no]['Type'] = chr_type
			#for adj in d:
			#	if adj in x:
			#		logfile.write(str(adj)+"\n")
			#	elif adj[::-1] in x:
			#		logfile.write(str(adj[::-1])+"\n")										


	#logfile.close()
	return soln_chr_dict	