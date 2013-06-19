#!/usr/bin/gnuplot
set output "compare.pdf"
set terminal pdf enhanced
set logscale xy
plot "error_1storder.txt" u 1:2 w lp t "1st order", "error.txt" u 1:2 w lp t "2nd order", x, x**2

