#!/bin/gnuplot
set output "plot.pdf" 
set terminal pdf color enhanced size 6in,5in
set xlabel "time to complete run [s]"
set logscale xyx2y2
set autoscale fix
set yrange [1e-14:1]
set xrange [0.05:10]
set multiplot layout 3,2
unset key
set st d p


do for [i=0:5]{
	if (i==2) {
		set key bottom left  box; 
	}else{
		unset key;
	}


	if (i%2==0){
		set rmargin 0
		set lmargin 12
		set ytics 
		set ylabel "relative energy after 10000 orbits"
	}else{
		set lmargin 0
		set rmargin 12
		unset ylabel
		unset ytics 
	}

	titfile=system("sed '".(i+1)."q;d' titles.txt");
	set label 1 titfile at graph 0.99,0.9 right


	plot \
	"testcase_".i."/energy_ias15.txt" t "REBOUND IAS15", \
	"testcase_".i."/energy_ra15.txt" t "REBOUND RA15", \
	"testcase_".i."/energy_wh.txt" t "      REBOUND WH",  \
	"testcase_".i."/energy_bs2.txt" t "MERCURY BS2",  \
	"testcase_".i."/energy_radau.txt" t "MERCURY RADAU",  \
	"testcase_".i."/energy_mvs.txt" t "      MERCURY MVS" lt 7,  \
	"testcase_".i."/energy_ias15_canonical.txt" notit ps 5 lt 6, \
}
