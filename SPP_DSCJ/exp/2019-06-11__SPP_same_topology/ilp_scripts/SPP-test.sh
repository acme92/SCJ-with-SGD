#!/bin/bash
for count in 1 2 3
do
        for alpha in 0 0.25 0.5 0.75 1
        do
                echo $count $alpha
                time python2.7 SPP_test.py $count $alpha
        done        
done