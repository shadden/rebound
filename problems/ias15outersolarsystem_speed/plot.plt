#!/bin/gnuplot
set output "plot.pdf" 
set terminal pdf color enhanced size 3in,2in
set xlabel "time to complete run [s]"
set ylabel "relative energy after 1000 yrs"
set logscale xy
set autoscale fix
set key center left
set yrange [1e-17:1]

set st d lp

plot \
"energy_wh.txt" t "Wisdom Holman",  \
"energy_ias15.txt" t "IAS15", \
1e-16 t "Machine precision", \
