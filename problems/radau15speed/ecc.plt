#!/bin/gnuplot
set output "ecc.pdf"
set term pdf color enhanced size 3in,2in
set xlabel "1-e"
set ylabel "relative energy error after 10^{4} orbits "

set logscale xy
set autoscale xfix
set format y "10^{%T}"

set st d lp 

plot \
"< awk '{if ($2==\"1.000000e+00\") print $0}' energy_radau15.txt" u (1.-$3):4 t "{/Symbol e}=10^{0}", \
"< awk '{if ($2==\"1.000000e-02\") print $0}' energy_radau15.txt" u (1.-$3):4 t "{/Symbol e}=10^{-2}", \
"< awk '{if ($2==\"1.000000e-04\") print $0}' energy_radau15.txt" u (1.-$3):4 t "{/Symbol e}=10^{-4}"
