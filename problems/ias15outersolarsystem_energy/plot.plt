#!/bin/gnuplot
set output "plot.pdf" 
set terminal pdf color enhanced size 4in,3in
set xlabel "timestep [days]"
set ylabel "relative energy after 1000 yrs"
set logscale xy
set autoscale fix
set yrange [1e-17:0.1]
set key right center

set st d lp

plot \
"energy_wh.txt" t "Wisdom Holman",  \
"energy_ias15.txt" t "IAS15", \
1e-3*(x/1000.)**2 t "dt^{2}", \
1e-12*(x/1000.)**15 t "dt^{15}", \
1e-16 t "machine precision"
