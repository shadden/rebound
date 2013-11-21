#!/bin/bash
function runtest {
	points=50
	min=-2
	max=3 
	for i in $(seq 0 $points)
	do 
		exp=$(echo "scale=10; ($max-($min))/$points*$i+($min) " |bc)
		timesteps=$(echo "scale=10; e($exp*l(10))"  | bc -l )
		./nbody --timesteps=$timesteps
		#echo $timesteps
	done
}

rm energy_*.txt
make leapfrog
runtest
make radau15



points=50
min=-5
max=-1
for i in $(seq 0 $points)
do 
	exp=$(echo "scale=10; ($max-($min))/$points*$i+($min) " |bc)
	timesteps=$(echo "scale=10; e($exp*l(10))"  | bc -l )
	./nbody --epsilon=$timesteps
	#echo $timesteps
done
mv energy_radau15.txt energy_radau15_adaptive.txt


runtest




make wh
runtest
