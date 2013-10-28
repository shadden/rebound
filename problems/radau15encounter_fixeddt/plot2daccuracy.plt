#!/bin/gnuplot
set term pdf enhanced color
set output "plot2daccuracy.pdf"
load "default.plt"
set logscale cb
set autoscale fix
set format y "10^{%.0f}"
set ytics 1
set ylabel "timestep [orbits]"
set xlabel "time [orbits]"
set cblabel "\nmagnitude of last term in series {/Symbol e}_{UB}"
set format cb "10^{%T}"; 
set cbrange [1e-13:1]
set yrange [-3.5:-0.5]

set key below
set pal gray
plot "lastfrac.txt" w image notit, \
"../radau15encounter_adaptivedt/dt_1e-2.txt" u ($1/2./pi):(log10($2/2./pi)) w l lt 3 t "{/Symbol e}=10^{-2}",\
"../radau15encounter_adaptivedt/dt_1e-4.txt" u ($1/2./pi):(log10($2/2./pi)) w l lt 1 t "{/Symbol e}=10^{-4}",\
"../radau15encounter_adaptivedt/dt_1e-6.txt" u ($1/2./pi):(log10($2/2./pi)) w l lt 2 t "{/Symbol e}=10^{-6}",\
"../radau15encounter_adaptivedt/dt_1e-8.txt" u ($1/2./pi):(log10($2/2./pi)) w l lt 9 t "{/Symbol e}=10^{-8}",\



