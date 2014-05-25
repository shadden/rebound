#!/bin/gnuplot
set output "plot.pdf"
set terminal pdf enhanced dashed size 4 in,4 in
set st d lp
set key top right
set logscale xy
set xlabel "timesteps/orbits"
set yrange [:1]

set ylabel "relative error after 10 orbits"
plot \
"energy_ias15.txt" u (1./$1):(abs($5)) t "ias15(fixed) offtrack error", \
'' u (1./$1):(abs($6)) t "ias15(fixed) ontrack error", \
'' u (1./$1):(abs($8)) t "ias15(fixed) energy error", \
"energy_ias15variable.txt" u (1./$1):(abs($5)) t "ias15(variable) offtrack error", \
'' u (1./$1):(abs($6)) t "ias15(variable) ontrack error", \
'' u (1./$1):(abs($8)) t "ias15(variable) energy error", \
"energy_wh.txt" u (1./$1):(abs($5)) t "wh offtrack error", \
'' u (1./$1):(abs($6)) t "wh ontrack error"

set ylabel "relative error per timestep"
plot \
"energy_ias15.txt" u (1./$1):(abs($5/$7)) t "ias15(fixed) offtrack error", \
'' u (1./$1):(abs($6/$7)) t "ias15(fixed) ontrack error", \
'' u (1./$1):(abs($8/$7)) t "ias15(fixed) energy error", \
"energy_ias15variable.txt" u (1./$1):(abs($5)/$7) t "ias15(variable) offtrack error", \
'' u (1./$1):(abs($6/$7)) t "ias15(variable) ontrack error", \
'' u (1./$1):(abs($8/$7)) t "ias15(variable) energy error", \
"energy_wh.txt" u (1./$1):(abs($5/$7)) t "wh offtrack error", \
'' u (1./$1):(abs($6/$7)) t "wh ontrack error"

set xlabel "epsilon"

set ylabel "relative error after 10 orbits"
plot \
"energy_ias15variable.txt" u ($3):(abs($5)) t "ias15(variable) offtrack error", \
'' u ($3):(abs($6)) t "ias15(variable) ontrack error", \
'' u ($3):(abs($8)) t "ias15(variable) energy error", \
1e-15*x**(15./7.) t "equation 20"

set ylabel "relative error per timestep"
plot \
"energy_ias15variable.txt" u ($3):(abs($5)/$7) t "ias15(variable) offtrack error", \
'' u ($3):(abs($6/$7)) t "ias15(variable) ontrack error", \
'' u ($3):(abs($8/$7)) t "ias15(variable) energy error", \
1e-15*x**(15./7.) t "equation 20"
