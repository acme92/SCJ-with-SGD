# SCJ-with-SGD
The directed distance, unrooted median and the rooted median problems under the SCJTDFD model

## dist_DSCJ: the directed Single Cut or Join distance with Tandem Duplications and Floating Duplications
Pedro Feijao, Aniket Mane, Cedric Chauve  
School of Computing Science and Department of Mathematics, Simon Fraser University  
Contact: cedric.chauve@sfu.ca

### Overview

SCJTDFD is a program that implements genome rearrangement algorithms in the directed Single Cut or Join with Tandem Duplications and Floating Duplications (SCJTDFD) model. This model and the related algorithms are described in the paper "A tractable variant of the Single Cut or Join distance with duplicated genes" published in the proceedings of [RECOMB-CG 2017](http://www.crg.eu/en/event/15th-recomb-comparative-genomics-satellite-workshop), Lecture Notes in Boinformatics (LNBI), volume 10562, Springer (DOI: 10.1007/978-3-319-67979-2).

It is implements three algorithms:
* Computing the directed SCJTDFD distance between an ancestral genome A and a descendant genome D.
* Computing a parsimonious directed SCJTDFD scenario between an ancestral genome A and a descendant genome D.
* Computing a parsimonious directed median genome A of k descendant genomes D1, ..., Dk.

### Requirements

The code is composed of a set of Python scripts. It requires the following to be available:

* python (3)
* networkx library, which is available at: https://networkx.github.io/

### Usage

Use the following command to run the script: <br>

./DSCJ.sh dist -s/-m/-d inputfile outputfile <br>
OR <br>
bash ./DSCJ.sh dist -s/-m/-d inputfile outputfile <br>

Use -d for computing the dSCJTDFD distance between the two input genomes <br>
    -s for computing a parsimonious dSCJTDFD scenario between the two input genomes <br>
    -m for computing a parsimonious dSCJTDFD median of k genomes
    
The input file should have the following format.
1. Each genome should start with the genome name followed by the genome itself in the next non-empty and non-commented line.
2. Each line that is read contains exactly one chromosome.
3. Linear chromosomes should end with '|' while circular chromosomes with ')'. 
4. Any gene in reverse orientation should be preceded by a '-' sign.
5. Lines containing the genome name should not end with '|' or ')'.
6. Any line that is a comment should start with a '#'. 
7. For -d and -s, the input file should contain exactly two genomes, both having the same set of gene families.
8. For -m, the input file should have at least two genomes and all genomes should have the same set of gene families.

Example: <br>
Genome 1: <br>
Chr1_name a -b c | <br>
Chr2_name -d e f g ) <br>

Genome 2: <br>
Chr1_name a -b c d ) <br>
Chr2_name e -f g | <br>

###Output file format
Every implementation (-d/-s/-m) of DSCJ.py yields two files - a log file and an output file.
The output files contain the following information:

* Distance (-d):
- 	The SCJTDFD distance and a breakdown of events (SCJ, TD, FD, tandem arrays and single gene circular chromosomes)
- 	Two lists of adjacencies, Cuts and Joins. Every element of the list is of the format: 
	[('gene1_name','extremity'),('gene2_name','extremity')], where extremity = 'h' or 't'.

* Scenario (-s):
-	The SCJTDFD distance and a breakdown of events (SCJ, TD, FD, tandem arrays and single gene circular chromosomes)
-	Two lists of genes, Floating Duplications and Tandem Duplications. Every gene in the list has a copy number associated with it.
	For instance, the 3rd copy of gene 'a' seen in D will be named as 'acopy3'. This is done to differentiate genes from the same family.
- 	Two lists of adjacencies, Cuts and Joins. Every element of the list is of the format: 
	[('gene1_name','extremity'),('gene2_name','extremity')], where extremity = 'h' or 't'.

* Median (-m):
-	List of adjacencies retained in the median genome.
	(('gene1_name','extremity'),('gene2_name','extremity')), where extremity = 'h' or 't'.
-	Score of the median genome.	

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
./DSCJ.sh rMedian inputfolder tree geneorders orthology outputfolder <br>
OR <br>
bash ./DSCJ.sh rMedian inputfolder tree geneorders orthology outputfolder 
Note: It is assumed that the 3 input files (species tree, orthology and gene orders) are all in the same input folder.  
The input files should have the following format.
1. Any line that is a comment should start with a '#'.
2. The tree file should contain the details of the species tree. Each line should have the species and the parents of the species, in that order, separated by a tab ("\t").
3. The geneorders file should contain the gene orders of the known genomes (an ancestor A and k descendant genomes D1, ..., Dk). Each line  should be tab-separated and should contain the details of one gene in the genome namely the species name, scaffold/chromosome index, gene name and orientation (+/-).
4. All genes from the same (species, scaffold) combination should appear in subsequent lines, maintaining the order from the scaffold in which they feature. 
5. The orthology file should contain the details of the orthology relations. Each line should have exactly one relation, presenting the ancestor species name, descendant species name, ancestor gene name, descendant gene name and gene family name all separated by tabs.


Example: <br>
1. GeneOrders: <br>
#Species Scaffold Gene Orientation <br> 
A	1	a1	+ <br>
A	1	b1	+ <br>
A	1	c1	+ <br>
D2	1	b7	- <br>
D2	1	a4	+ <br>
D2	1	c5	+ <br>
D2	2	b8	+ <br>
D2	2	a5	+ <br>
D2	2	b9	- <br>

2. Orthology: <br>
#Ancestor_species_id Descendant_species_id Ancestor_gene_name Descendant_gene_name Gene_tree <br>
A	M	b1	b3	b <br>
A	M	c1	c2	c <br>
M	D1	a2	a3	a <br>
M	D2	a2	a4	a <br>
