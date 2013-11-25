#!/bin/gnuplot
set output "prdrag.pdf"
set key top left
set terminal pdf monochrome dashed enhanced size 3in,3in
set xlabel "time [years]"
set ylabel "semimajor axis [AU]"
set multiplot layout 2,1
beta = 0.01
set lmargin 12
set ytics 0.1
k = 2.497557889905430e-03*beta
a(t) = sqrt(1.-k*t)
set st d  l
set autoscale xfix 
set xtics 10000

plot "r.txt" u ($1/2./pi):($2) notit

set ytics 1e-4
set ylabel "semimajor axis error [AU]"
plot "r.txt" u ($1/2./pi):($2-a($1/2./pi)) notit


