#!/bin/bash
function runtest {
	points=100
	min=2
	max=-6
	for i in $(seq $points)
	do 
		exp=$(echo "scale=10; ($max-($min))/$points*$i+($min) " |bc)
		epsilon=$(echo "scale=10; e($exp*l(10))"  | bc -l )
		./nbody --epsilon=$epsilon
		#echo $epsilon
	done
}

rm energy_*.txt
rm energy.txt
rm energy.plot

make radau15
runtest
