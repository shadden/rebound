#!/bin/gnuplot
set output "plot.pdf" 
set terminal pdf color enhanced size 3in,2in
set xlabel "time to complete run [s]"
set ylabel "relative energy after 10000 yrs"
set logscale xy
set autoscale fix
set key below
set yrange [:1]

set st d p

plot \
"energy_ias15.txt" t "REBOUND IAS15", \
"energy_ra15.txt" t "REBOUND RA15", \
"energy_bs.txt" t "MERCURY BS",  \
"energy_bs2.txt" t "MERCURY BS2",  \
"energy_radau.txt" t "MERCURY RADAU",  \
"energy_hybrid.txt" t "MERCURY HYBRID",  \
"energy_mvs.txt" t "MERCURY MVS",  \
