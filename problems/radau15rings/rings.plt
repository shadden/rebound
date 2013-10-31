#!/bin/gnuplot
set output "rings.pdf" 
set term pdf enhanced color size 10in,2.5in
set multiplot layout 1,4
set size ratio -1
set xtics 0.0005
set ytics 0.0005
set xlabel "x [AU]"
set ylabel "y [AU]"
set xrange [-0.001:0.001]
set yrange [-0.001:0.001]

set object 1 circle at  0,0 size first 0.00038925688 fc rgb "black" fs solid 

set title "t = 0 orbits"
plot "out__betaparticles_1.000e-02/ascii_0000.txt" u 2:3 notit lt 7 lc -1 ps 0.4
set title "t = 1200 orbits"
plot "out__betaparticles_1.000e-02/ascii_0013.txt" u 2:3 notit lt 7 lc -1 ps 0.4
set title "t = 2400 orbits"
plot "out__betaparticles_1.000e-02/ascii_0024.txt" u 2:3 notit lt 7 lc -1 ps 0.4
set title "t = 3600 orbits"
plot "out__betaparticles_1.000e-02/ascii_0036.txt" u 2:3 notit lt 7 lc -1 ps 0.4
