 #plot/plotfile.pdf
set terminal pdf
set output 'plot/plotfile.pdf'
set title ''
set nokey
set grid
set xlabel 'x'
set ylabel 'y'
xydata = 'plot/plotfile.dat'
plot xydata using 1:2 with linespoints
