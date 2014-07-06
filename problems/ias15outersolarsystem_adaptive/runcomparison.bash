#!/bin/bash
function runepsilon {
	echo "Running mercury epsilon $1"
	rm -f ../energy_$1.txt
	points=20
	min=-10
	max=-5
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
		utime="$( TIMEFORMAT='%R';time ( doalarm 10  ./mercury6 ) 2>&1 1>/dev/null )"
		energy="$(../mercury_read/mercury_energy big.in big.dmp)"
		if [[ $utime == *Alarm* ]]; then
			echo "Did not finish in time."
		else
			echo "$utime $energy $e" >> ../energy_$1.txt 
			echo "$utime $energy $e"  
		fi
		rm -f param.template?

	done
}
function rundt {
	echo "Running mercury dt $1"
	rm -f ../energy_$1.txt
	points=7
	min=0
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
		utime="$( TIMEFORMAT='%R';time ( doalarm 10 ./mercury6 ) 2>&1 1>/dev/null )"
		energy="$(../mercury_read/mercury_energy big.in big.dmp)"
		if [[ $utime == *Alarm* ]]; then
			echo "Did not finish in time."
		else
			echo "$utime $energy $e" >> ../energy_$1.txt 
			echo "$utime $energy $e"  
		fi
		rm -f param.template?

	done
}
function runepsilonnbody {
	echo "Running REBOUND epsilon $1"
	rm -f energy_$1.txt
	points=10
	min=$2
	max=$3
	for i in $(seq 0 $points)
	do 
		exp=$(echo "scale=10; ($max-($min))/$points*$i+($min) " |bc)
		e=$(echo "scale=10; e($exp*l(10))"  | bc -l )
		utime="$( TIMEFORMAT='%R';time ( doalarm 10 ./nbody --integrator_epsilon=$e 2>&1 ) 2>&1 1>/dev/null )"
		energy="$(cat energy.txt)"
		if [[ $utime == *Alarm* ]]; then
			echo "Did not finish in time."
		else
			echo "$utime $energy $e" >> energy_$1.txt 
			echo "$utime $energy $e"  
		fi
		rm -f energy.txt

	done
}

make problemgenerator



for i in $(seq 0 3)
do
	echo "Running test case $i"

	./problemgenerator --testcase=$i

	make ra15
	runepsilonnbody ra15 -10 -8

	make ias15
	runepsilonnbody ias15 -3 0

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

	rm -rf testcase_$1
	mkdir testcase_$1
	mv energy*.txt testcase_$1/

done
exit

