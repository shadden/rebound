#!/bin/bash
function runtest {
	points=100
	max=10 
	for i in $(seq $points)
	do 
		exp=$(echo "scale=10; $max/$points*$i " |bc)
		epsilon=$(echo "scale=10; e(-$exp*l(10))"  | bc -l )
		./nbody --epsilon=$epsilon
		#echo $timesteps
	done
}

rm energy_*.txt

make radau15
runtest
