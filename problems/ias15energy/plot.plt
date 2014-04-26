#!/bin/gnuplot
set output "plot.pdf"
set terminal pdf enhanced dashed size 4 in,3 in
set xlabel "cputime [s]"
set ylabel "maximum relative energy error"
set st d lp
set key top right
set logscale x
set yrange [-16:0]
set autoscale fix
set format y "10^{%.0f}" 
set ytics 4


set arrow from 2, -15 to 4, -13 
set label 'machine precission ' at 2, -15 right
plot "energy_ias15_variable.txt" u 5:(log10($4)) t "IAS15, variable timestep" , "energy_ias15.txt" u 5:(log10($4)) t "IAS15, fixed timestep" , "energy_wh.txt" u 5:(log10($4)) t "WH, fixed timestep, symplectic, mixed variable integrator", log10(0.000001*x**(-2)) t "2nd order",  log10(0.001*x**(-15)) t "15th order"
