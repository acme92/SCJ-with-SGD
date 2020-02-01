## rMedian_DSCJ: the rooted median problem under the Single Cut or Join model with Tandem Duplications and Floating Duplications 

Aniket Mane, Manuel Lafond, Pedro Feijao, Cedric Chauve  
Department of Mathematics, Simon Fraser University
Department of Computer Science, Univeristy  of Sherbrooke
School of Computing Science, Simon Fraser University 
Contact: amane@sfu.ca

### Overview
RMedian_DSCJ is a program that computes the rooted median of a given instance under the SCJTDFD model. It uses an Integer Linear Program (ILP) framework, described in (link to paper).
An instance of the rooted median problem involves an ancestor A and k descendant genomes D1, ..., Dk. Given these genomes, it is desired to compute a median genome M of the (k+1) genomes.

### Requirements
The code is composed of a set of Python scripts. It requires the following to be available:
* python (3)
* Gurobi optimization solver (version 7.5.2 or higher), which is available at: http://www.gurobi.com/downloads/download-center

### Usage
Use the following command to run the script:
python RMedian_main.py tree geneorders orthology outputfolder
  
The input files should have the following format.
1. Any line that is a comment should start with a '#'.
2. The tree file should contain the details of the species tree. Each line should have the species and the parents of the species, in that order, separated by a tab ("\t").
3. The geneorders file should contain the gene orders of the known genomes (an ancestor A and k descendant genomes D1, ..., Dk). Each line  should be tab-separated and should contain the details of one gene in the genome namely the species name, scaffold/chromosome index, gene name and orientation (+/-).
4. All genes from the same (species, scaffold) combination should appear in subsequent lines, maintaining the order from the scaffold in which they feature. 
5. The orthology file should contain the details of the orthology relations. Each line should have exactly one relation, presenting the ancestor species name, descendant species name, ancestor gene name, descendant gene name and gene family name all separated by tabs.


Example:
1. GeneOrders:<br/>
#Species Scaffold Gene Orientation<br/>
A	1	a1	+<br/>
A	1	b1	+<br/>
A	1	c1	+<br/>
D2	1	b7	-<br/>
D2	1	a4	+<br/>
D2	1	c5	+<br/>
D2	2	b8	+<br/>
D2	2	a5	+<br/>
D2	2	b9	-<br/>

2. Orthology:
#Ancestor_species_id Descendant_species_id Ancestor_gene_name Descendant_gene_name Gene_tree<br/>
A	M	b1	b3	b<br/>
A	M	c1	c2	c<br/>
M	D1	a2	a3	a<br/>
M	D2	a2	a4	a<br/>
