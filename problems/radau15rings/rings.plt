#!/bin/gnuplot
set output "rings.pdf" 
set term pdf enhanced color size 2.5in,2.5in
set size ratio -1
set xtics 0.0005
set ytics 0.0005
set xlabel "x [AU]"
set ylabel "y [AU]"

set object 1 circle at  0,0 size first 0.00038925688 fc rgb "black" fs solid 

plot "out__/ascii.txt" u 2:3 notit lt 7 lc -1 ps 0.4
