#!/bin/gnuplot
set output "plot.pdf"
set terminal pdf enhanced dashed size 4 in,3 in
set st d lp
set key top left
set logscale xy
set xlabel "timestepis/orbits"
set ylabel "relative error"
set yrange [:1]

plot \
"energy_ias15.txt" u (1./$1):(abs($5)) t "ias15 energy error", \
'' u (1./$1):(abs($6)) t "ias15 phase error", \
"energy_wh.txt" u (1./$1):(abs($5)) t "wh energy error", \
'' u (1./$1):(abs($6)) t "wh phase error"
