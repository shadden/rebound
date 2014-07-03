#!/bin/bash
function runepsilon {
	echo "Running mercury epsilon $1"
	rm -f ../energy_$1.txt
	points=20
	min=-10
	max=-2
	for i in $(seq 0 $points)
	do 
		exp=$(echo "scale=10; ($max-($min))/$points*$i+($min) " |bc)
		e=$(echo "scale=10; e($exp*l(10))"  | bc -l )
		rm -f *.tmp
		rm -f *.dmp
		rm -f *.out
		sed "s/EPSILON/$e/g" param.template > param.template2
		sed "s/SCHEME/$1/g" param.template2 > param.template3
		sed "s/STEPSIZE/1/g" param.template3 > param.in
		utime="$( TIMEFORMAT='%R';time ( ./mercury6 ) 2>&1 1>/dev/null )"
		energy="$(../mercury_read/mercury_energy big.in big.dmp)"
		echo "$utime $energy $e" >> ../energy_$1.txt 
		echo "$utime $energy $e"  

	done
}
function rundt {
	echo "Running mercury dt $1"
	rm -f ../energy_$1.txt
	points=4
	min=1
	max=4
	for i in $(seq 0 $points)
	do 
		exp=$(echo "scale=10; ($max-($min))/$points*$i+($min) " |bc)
		e=$(echo "scale=10; e($exp*l(10))"  | bc -l )
		rm -f *.tmp
		rm -f *.dmp
		rm -f *.out
		sed "s/EPSILON/1/g" param.template > param.template2
		sed "s/SCHEME/$1/g" param.template2 > param.template3
		sed "s/STEPSIZE/$e/g" param.template3 > param.in
		utime="$( TIMEFORMAT='%R';time ( ./mercury6 ) 2>&1 1>/dev/null )"
		energy="$(../mercury_read/mercury_energy big.in big.dmp)"
		echo "$utime $energy $e" >> ../energy_$1.txt 
		echo "$utime $energy $e"  

	done
}
function runepsilonnbody {
	echo "Running REBOUND epsilon"
	rm -f ../energy_$1.txt
	points=20
	min=-3
	max=3
	for i in $(seq 0 $points)
	do 
		exp=$(echo "scale=10; ($max-($min))/$points*$i+($min) " |bc)
		e=$(echo "scale=10; e($exp*l(10))"  | bc -l )
		utime="$( TIMEFORMAT='%R';time ( ./nbody --integrator_epsilon=$e 2>&1 ) 2>&1 1>/dev/null )"
		energy="$(cat energy.txt)"
		echo "$utime $energy $e" >> energy_ias15.txt 
		echo "$utime $energy $e"  

	done
}

rm -f energy_*.txt
make ias15
runepsilonnbody


pushd mercury
rm -f *.tmp
rm -f *.dmp
rm -f *.out
rm -f output.txt
runepsilon bs 
runepsilon bs2
runepsilon radau 
runepsilon hybrid 
rundt mvs 
popd

