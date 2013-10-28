#!/bin/gnuplot
set output "convergence.pdf" 
set term pdf enhanced color size 2.5in,2in
set logscale xy
set xlabel "timestep, dt"
set ylabel "relative energy error, E"
set autoscale fix
set yrange [1e-16:100]
set xrange [0.02:3]
set ytics 1000
set key top left
set format y "10^{%T}"

set style arrow 7 nohead ls 1
set palette gray negative
set cblabel "\nfraction of timesteps rejected"

plot "energy_radau15.txt" u 1:(1e-16):(0):(100):(1*$3/(2.*pi/$1)) with vectors nohead ls 1  lc palette notit, \
"energy_radau15.txt" w p pt 7 ps 0.4 lt -1 notit, x**15*1e5 t "dt^{15}" lt 1, \
