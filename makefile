# COC473 - Álgebra Linear Computacional @ ECI/UFRJ
#Makefile
matrix: matrixlib.o main.o
	gfortran -o matrix matrixlib.o main.o

matrixlib.mod: matrixlib.o matrixlib.f95
	gfortran -c matrixlib.f95

matrixlib.o: matrixlib.f95
	gfortran -c matrixlib.f95

main.o: matrixlib.mod main.f95
	gfortran -c main.f95

clean:
	rm matrixlib.o main.o
#End Makefile
