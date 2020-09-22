# COC473 - Álgebra Linear Computacional @ ECI/UFRJ
# Macros
f95c = gfortran
args = -O2 -ffpe-summary=none
# Makefile
matrix1: matrixlib.o main1.o
	$(f95c) -o matrix1 $(args) matrixlib.o main1.o

matrix2: matrixlib.o main2.o
	$(f95c) -o matrix2 $(args) matrixlib.o main2.o

%.mod: %.o %.f95
	$(f95c) -c $(args) $<

%.o: %.f95
	$(f95c) -c $(args) $<

clean:
	rm *.o *.mod matrix
# End Makefile