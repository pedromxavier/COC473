 #plot/ex1.pdf
set size 1,1
set terminal pdf
set output "plot/ex1.pdf"
set title "y'(t) = -2 t (y)²; y(0) = 1"
set grid
set xlabel "t"
set ylabel "y(t)"
xydata = "plot/ex1.dat"

plot xydata u 1:2 t "Euler" with lines,\
 "" u 1:3 t "Runge-Kutta II" with lines,\
 "" u 1:4 t "Runge-Kutta IV" with lines,\
  "" u 1:5 t "y(t)" dashtype '.-_' with lines

plot x t "f(x) = x" with lines
