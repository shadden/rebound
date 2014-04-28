#!/bin/gnuplot
set output "plot.pdf"
set terminal pdf enhanced size 4 in,3 in
set xlabel "simulation time [orbits]"
set ylabel "relative energy error"
set st d l
set key bottom right
set yrange [-16:-2]
set autoscale fix
set format y "10^{%.0f}" 
set ytics 4

set arrow from 4000, -12 to 5000, -13 
set label 'machine precission, error grows like t^{1/2} ' at 4000, -12 right
set arrow from 3000, -9 to 4000, -8 
set arrow from 3000, -9 to 4000, -9.2 
set label 'errors grow linearly ' at 3000, -9 right
set arrow from 3000, -3 to 4000, -4 
set label 'bound error ' at 3000, -3 right

set arrow from 6500, -12 to 6500, -13 
set label '{/Symbol e} = 10^{-4} ' at 6500, -12 right

set arrow from 6500, -9 to 6500, -8 
set label '{/Symbol e} = 1 ' at 6500, -9 right

plot "energy_ias15_variable.txt" every 20 u ($1/2./pi):(log10($5)) t "IAS15, variable timestep" , "energy_ias15.txt" u ($1/2./pi):(log10($5)) every 20  t "IAS15, fixed timestep" , "energy_wh.txt" u ($1/2./pi):(log10($5)) every 20 t "WH, fixed timestep, symplectic, mixed variable integrator"
