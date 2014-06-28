#!/bin/gnuplot
set output "plot.pdf"
set terminal pdf color enhanced
set xlabel "timestep [days]"
set ylabel "relative energy after 1Myr"
set logscale xy
set autoscale fix

set st d lp

plot \
"energy_wh.txt" t "Wisdom Holman",  \
"energy_ias.txt" t "IAS15, \
