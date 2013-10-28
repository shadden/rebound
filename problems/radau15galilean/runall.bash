#!/bin/bash
function runtest {
	points=100
	min=-2
	max=3 
	for i in $(seq 0 $points)
	do 
		exp=$(echo "scale=10; ($max-($min))/$points*$i+($min) " |bc)
		v=$(echo "scale=10; e($exp*l(10))"  | bc -l )
		./nbody --e=0.01 --v=$v --epsilon=1e-2
	#	echo $v
	done
}

rm energy_radau15.txt
make radau15
runtest 
