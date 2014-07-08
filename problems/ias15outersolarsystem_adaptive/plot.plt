#!/bin/gnuplot
set output "plot.pdf" 
set terminal pdf color enhanced size 6in,5in
set logscale xyx2y2
set autoscale fix
set yrange [1e-14:1]
set xrange [0.05:10]
set multiplot layout 3,2
unset key
set st d p

bottommargin = 0.08
keymargin = 0.1
topmargin = 0.05
ny = 3
px = 0

do for [i=0:5]{
	if (i==2) {
		set key at screen 0.5,0.05 center horizontal 
	}else{
		unset key;
	}


	if (i%2==0){
		set rmargin 0
		set lmargin 12
		set format y "%g" 
		if (i/2==(ny-1)/2){
			set ylabel "relative energy"
		}else{
			unset ylabel
		}
		py = 0
	}else{
		set lmargin 0
		set rmargin 12
		unset ylabel
		set format y  ""
		py = 0.5
	}
	if (i/2==ny-1){
		px=keymargin
	}else{
		px=(1.-bottommargin-topmargin-keymargin)/ny*(ny-i/2-1)+bottommargin+keymargin

	}

	set origin py,px

	if (i/2==ny-1){
		set size 0.5,(1.-bottommargin-topmargin-keymargin)/ny+bottommargin
		set xlabel "time to complete run [s]"
		set bmargin 5
		set format x "%g" 
	}else{
		set size 0.5,(1.-bottommargin-topmargin-keymargin)/ny
		set format x ""
		set bmargin 0
		set tmargin 0
	}

	titfile=system("sed '".(i+1)."q;d' titles.txt");
	set label 1 titfile at graph 0.99,0.1 right


	plot \
	"testcase_".i."/energy_ias15.txt" t "REBOUND IAS15", \
	"testcase_".i."/energy_ias15_canonical.txt" u (0.0001):(1.) t "IAS15 DEFAULT" ps 2 lt 6, \
	"testcase_".i."/energy_ra15.txt" t "REBOUND RA15", \
	"testcase_".i."/energy_wh.txt" t "      REBOUND WH",  \
	"testcase_".i."/energy_bs2.txt" t "MERCURY BS2",  \
	"testcase_".i."/energy_radau.txt" t "MERCURY RADAU",  \
	"testcase_".i."/energy_mvs.txt" t "      MERCURY MVS" lt 7,  \
	"testcase_".i."/energy_ias15_canonical.txt" notit ps 4 lt 6, \
}
