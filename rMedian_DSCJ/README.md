## RMedian_DSCJ: the rooted median problem under the Single Cut or Join model with Tandem Duplications and Floating Duplications 

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

Example:
