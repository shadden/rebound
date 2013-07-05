#!/usr/bin/gnuplot
set output "radii2d.png"
set key top left
set terminal png enhanced 
set xlabel "time [t_0]"
set ylabel "% mass"
set logscale cb
set cbrange [0.5:4]
set autoscale fix
plot "radii.txt" w image 

