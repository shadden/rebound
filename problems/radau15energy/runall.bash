#!/bin/bash
function runtest {
	points=100
	min=0
	max=3 
	for i in $(seq $points)
	do 
		exp=$(echo "scale=10; ($max-$min)/$points*$i+$min " |bc)
		timesteps=$(echo "scale=10; e($exp*l(10))"  | bc -l )
		./nbody --timesteps=$timesteps
		#echo $timesteps
	done
}

rm energy_*.txt
rm energy.txt
rm energy.plot

make radau15
runtest
