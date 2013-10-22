#!/bin/bash
function runtest {
	points=512
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
rm lastfrac.txt 

make radau15
runtest
