#!/bin/bash
function runtest {
	points=1000
	max=5 
	for i in $(seq $points)
	do 
		exp=$(echo "scale=10; $max/$points*$i " |bc)
		timesteps=$(echo "scale=10; e($exp*l(10))"  | bc -l )
		./nbody --timesteps=$timesteps
		#echo $timesteps
	done
}

rm energy_*.txt

make radau15
runtest
