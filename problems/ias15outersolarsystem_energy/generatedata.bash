#!/bin/bash
function rundt {
	points=100
	min=2
	max=4
	for i in $(seq 0 $points)
	do 
		exp=$(echo "scale=10; ($max-($min))/$points*$i+($min) " |bc)
		e=$(echo "scale=10; e($exp*l(10))"  | bc -l )
		./nbody --dt=$e 
	done
}

rm energy_*.txt
make wh
rundt
mv energy.txt energy_wh.txt

make ias15
rundt
mv energy.txt energy_ias15.txt

 
