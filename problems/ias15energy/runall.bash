#!/bin/bash
function runepsilon {
	points=100
	min=-5
	max=1 
	for i in $(seq 0 $points)
	do 
		exp=$(echo "scale=10; ($max-($min))/$points*$i+($min) " |bc)
		e=$(echo "scale=10; e($exp*l(10))"  | bc -l )
		./nbody --epsilon=$e --e=0.1
	done
}

function rundt {
	points=200
	min=4
	max=7
	for i in $(seq 0 $points)
	do 
		exp=$(echo "scale=10; ($max-($min))/$points*$i+($min) " |bc)
		e=$(echo "scale=10; e($exp*l(10))"  | bc -l )
		./nbody --timesteps=$e --e=0.1
	done
}

rm energy_*.txt

make ias15
runepsilon
mv energy_ias15.txt energy_ias15_variable.txt
rundt

make wh
rundt
