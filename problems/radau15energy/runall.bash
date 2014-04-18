#!/bin/bash
function runtest {
	points=20
	min=0
	max=-5 
	for i in $(seq 0 $points)
	do 
		exp=$(echo "scale=10; ($max-$min)/$points*$i+$min " |bc)
		e=$(echo "scale=10; 1.-e($exp*l(10))"  | bc -l )
		./nbody --epsilon=$1 --e=$e
		#echo $1
	done
}

rm energy_radau15.txt
make radau15
runtest "1.0e-0"
runtest "1.0e-2"
runtest "1.0e-4"
