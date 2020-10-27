# COC473 - Álgebra Linear Computacional @ ECI/UFRJ
# Macros
f95c = gfortran
blas = -llapack -lblas
args = -O2 -ffpe-summary=none -fno-range-check $(blas)

lib = utillib.o matrixlib.o calclib.o plotlib.o funcmod.o

# Begin Makefile
all: matrix1 matrix2 matrix3 calc1 calc2

matrix1: $(lib) main1.o
	$(f95c) -o matrix1 $(args) $(lib) main1.o

matrix2: $(lib) main2.o
	$(f95c) -o matrix2 $(args) $(lib) main2.o

matrix3: $(lib) main3.o
	$(f95c) -o matrix3 $(args) $(lib) main3.o

calc1: $(lib) main4.o
	$(f95c) -o calc1 $(args) $(lib) main4.o

calc2: $(lib) main5.o
	$(f95c) -o calc2 $(args) $(lib) main5.o

calc3: $(lib) main6.o
	$(f95c) -o calc3 $(args) $(lib) main6.o

plottest: utillib.o plotlib.o plottest.o
	$(f95c) -o plottest $(args) utillib.o plotlib.o plottest.o

lib: $(lib)

%.mod: %.o %.f95
	$(f95c) -c $(args) $<

%.o: %.f95
	$(f95c) -c $(args) $<

clean:
	rm -f *.o *.mod plottest matrix[1-9] calc[1-9] plot/*.dat plot/*.plt
# End Makefile