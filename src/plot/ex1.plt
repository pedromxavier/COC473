 #plot/ex1.pdf
set terminal pdf
set output "plot/ex1.pdf"

set grid
set xlabel "x"
set ylabel "y"

h = 1.0;
n = 3.0;

set multiplot;                          # get into multiplot mode

set size 1, (h / n);  
##
set origin 0.0,(0.0 * h);
set title "y = sin(x)";
plot sin(x); 
##
set origin 0.0,((1.0/n) * h);
set title "y = cos(x)";
plot cos(x);
##
set origin 0.0,((2.0/n) * h);
set title "y = tan(x)";
plot tan(x);
unset multiplot                         # exit multiplot mode
