## SPP_DSCJ: The Small Parsimony Problem (SPP) under the Single Cut or Join model with Tandem Duplications and Floating Duplications 

Aniket Mane, Pedro Feijao, Cedric Chauve  
Department of Mathematics, Simon Fraser University
School of Computing Science, Simon Fraser University 
Contact: amane@sfu.ca

### Overview
SPP_DSCJ is a program that solves the Small Parsimony Problem under the SCJTDFD model, using an Integer Linear Program (ILP) framework.
In an instance of the SPP, we are provided with a phylogeny and the gene orders of extant genomes. Given these genomes, it is desired to compute the gene orders at the internal nodes of the tree.

### Requirements
The code is composed of a set of Python scripts. It requires the following to be available:
* python (3)
* Gurobi optimization solver (version 7.5.2 or higher), which is available at: http://www.gurobi.com/downloads/download-center

### Usage
Use the following command to run the script:
python SPP_main.py geneorders orthology adjacencies alpha outputfolder
  
The input files should have the following format.
1. Any line that is a comment should start with a '#'. 

Example:
