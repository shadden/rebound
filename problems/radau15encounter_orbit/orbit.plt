#!/bin/gnuplot
set terminal pdf enhanced dashed size  1.8in,1.8in
set size square 
set output "orbit.pdf"
set xlabel "x"
set ylabel "y"

set xrange [-0.20:0.20]
set yrange [-0.20:0.20]
set label "CM" at 0,0  center
set xtics 0.1
set ytics 0.1



plot "position1.txt" u 1:2 w l lt 2 notit, '' u 3:4 w l lt 3 notit

