#!/bin/gnuplot
set output "speed.pdf"
set term pdf color enhanced size 3in,2in
set xlabel "time to calculate 10^{4} years [s]"
set ylabel "relative energy error"
set logscale xy
set autoscale fix
set format y "10^{%T}"
set xrange[1e-3:10]
set st d lp
plot \
"energy_leapfrog.txt" u 3:2 t "leap frog", \
"energy_radau15.txt" u 3:2 t "radau15", \
"energy_wh.txt" u 3:2 t "wh", \

