#!/bin/sh

if [ "$1" == "rMedian" ];
then
	echo "Rooted median code to be executed"
	if [ "$#" -eq 6 ];
	then
		python2.7 rMedian_DSCJ/src/RMedian_main.py $2 $3 $4 $5 $6
		#RMedian_main.py inputfolder tree geneorders orthology outputfolder
		#Path to inputfolder: $2
		#Path to species tree: $3
		#Path to gene orders: $4
		#Path to orthology file: $5
		#Path to output folder: $6
	else
		echo "Error: 6 arguments needed for the rooted median code"	
	fi
elif [ "$1" == "dist" ];
then	
	echo "Distance code to be executed"
	if [ "$#" -eq 4 ];
	then
		python dist_DSCJ/src/DSCJ.py $2 $3 $4
		#Mode: $2 (-s/-d/-m)
		#Path to input file: $3
		#Path to output file: $4
	else
		echo "Error: 4 arguments needed for the distance/scenario/unrooted median code"		
	fi		
fi
