# COC473 - Álgebra Linear Computacional @ ECI/UFRJ
# Macros
f95c = gfortran
args = -O2 -ffpe-summary=none
blas = -llapack -lblas
# Makefile

all: matrix1 matrix2 matrix3 calc1

matrix1: matrixlib.o main1.o
	$(f95c) -o matrix1 $(args) matrixlib.o main1.o

matrix2: matrixlib.o main2.o
	$(f95c) -o matrix2 $(args) matrixlib.o main2.o

matrix3: matrixlib.o main3.o
	$(f95c) -o matrix3 $(args) matrixlib.o main3.o

calc1: calclib.o funcmod.o main4.o
	$(f95c) -o calc1 $(args) calclib.o funcmod.o main4.o $(blas)
%.mod: %.o %.f95
	$(f95c) -c $(args) $<

%.o: %.f95
	$(f95c) -c $(args) $<

clean:
	rm -f *.o *.mod matrix[1-9] calc[1-9]
# End Makefile