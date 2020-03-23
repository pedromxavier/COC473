# COC473 - Álgebra Linear Computacional @ ECI/UFRJ
# Macros
f95c = gfortran
args = -O2 -ffpe-summary=none
# Makefile
matrix: matrixlib.o main.o
	$(f95c) -o matrix $(args) matrixlib.o main.o

%.mod: %.o %.f95
	$(f95c) -c $(args) $<

%.o: %.f95
	$(f95c) -c $(args) $<

clean:
	rm *.o *.mod matrix
# End Makefile
