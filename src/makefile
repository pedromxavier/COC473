# COC473 - Álgebra Linear Computacional @ ECI/UFRJ
# Macros
f95c = gfortran
blas = -llapack -lblas
args = -O2 -ffpe-summary=none -fno-range-check $(blas)

# Begin Makefile
all: matrix1 matrix2 matrix3 calc1 calc2

matrix1: utillib.o matrixlib.o main1.o
	$(f95c) -o matrix1 $(args) utillib.o matrixlib.o main1.o

matrix2: utillib.o matrixlib.o main2.o
	$(f95c) -o matrix2 $(args) utillib.o matrixlib.o main2.o

matrix3: utillib.o matrixlib.o main3.o
	$(f95c) -o matrix3 $(args) utillib.o matrixlib.o main3.o

calc1: utillib.o matrixlib.o calclib.o funcmod.o main4.o
	$(f95c) -o calc1 $(args) utillib.o matrixlib.o calclib.o funcmod.o main4.o

calc2: utillib.o matrixlib.o calclib.o funcmod.o main5.o
	$(f95c) -o calc2 $(args) utillib.o matrixlib.o calclib.o funcmod.o main5.o

calc3: utillib.o matrixlib.o calclib.o plotlib.o funcmod.o main6.o
	$(f95c) -o calc3 $(args) utillib.o matrixlib.o calclib.o plotlib.o funcmod.o main6.o

%.mod: %.o %.f95
	$(f95c) -c $(args) $<

%.o: %.f95
	$(f95c) -c $(args) $<

clean:
	rm -f *.o *.mod matrix[1-9] calc[1-9]
# End Makefile