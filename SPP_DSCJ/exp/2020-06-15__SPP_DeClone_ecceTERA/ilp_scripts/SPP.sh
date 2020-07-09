#!/bin/bash
for mode in "without_LT_high_rearr"
do
	for run in {11..20}
	do
		for temp in 0.1 1
		do
			for rep in "a"
			do
				for alpha in 0 0.25 0.5 0.75 1
				do	
					echo $mode $run $temp $rep $alpha
                	time python2.7 SPP_main.py $mode $run $temp $rep $alpha
				done
			done
		done
	done
done
