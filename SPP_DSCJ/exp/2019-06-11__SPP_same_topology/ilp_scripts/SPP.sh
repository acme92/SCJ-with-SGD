#!/bin/bash
for mode in "with_L" "without_LT" "without_LT_high_rearr"
do
	for run in {1..20}
	do
		for temp in 01 1
		do
			for rep in "a" "b"
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
