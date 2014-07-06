#!/bin/gnuplot
set output "plot.pdf" 
set terminal pdf color enhanced size 6in,4in
set xlabel "time to complete run [s]"
set ylabel "relative energy after 10000 yrs"
set logscale xy
set autoscale fix
set key below
set yrange [:1]
set multiplot layout 2,2

set st d p
do for [i=0:3]{

	set title "Testcase ".i
	plot \
	"testcase_".i."/energy_ias15.txt" t "REBOUND IAS15", \
	"testcase_".i."/energy_ra15.txt" t "REBOUND RA15", \
	"testcase_".i."/energy_bs.txt" t "MERCURY BS",  \
	"testcase_".i."/energy_bs2.txt" t "MERCURY BS2",  \
	"testcase_".i."/energy_radau.txt" t "MERCURY RADAU",  \
	"testcase_".i."/energy_hybrid.txt" t "MERCURY HYBRID",  \
	"testcase_".i."/energy_mvs.txt" t "MERCURY MVS",  \

}
