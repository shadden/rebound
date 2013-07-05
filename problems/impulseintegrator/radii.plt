#!/usr/bin/gnuplot
set output "radii.png"
set key top left
set terminal png enhanced 
set xlabel "time [t_0]"
set ylabel "r [r_0]"
set logscale y
plot "radii.txt" u 1:2 t "10%",  "radii.txt" u 1:3 t "50%",  "radii.txt" u 1:4 t "90%" 

