#!/bin/gnuplot
set output "plot.pdf" 
set terminal pdf color enhanced size 6in,5in
set xlabel "time to complete run [s]"
set ylabel "relative energy after 10000 orbits"
set logscale xy
set autoscale fix
set yrange [1e-14:1]
set xrange [:10]
set multiplot layout 3,2
unset key
set st d p

do for [i=0:5]{
	if (i==2) {
		set key bottom left  box; 
	}else{
		unset key;
	}

	set title "Testcase ".i
	plot \
	"testcase_".i."/energy_ias15.txt" t "REBOUND IAS15", \
	"testcase_".i."/energy_ra15.txt" t "REBOUND RA15", \
	"testcase_".i."/energy_wh.txt" t "      REBOUND WH",  \
	"testcase_".i."/energy_bs2.txt" t "MERCURY BS2",  \
	"testcase_".i."/energy_radau.txt" t "MERCURY RADAU",  \
	"testcase_".i."/energy_mvs.txt" t "      MERCURY MVS",  \
}
