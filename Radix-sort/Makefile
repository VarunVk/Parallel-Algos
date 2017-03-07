.PHONY: all pthread openmp clean

all: mpi_clean  openmp_clean mpi openmp 

mpi: rs_mpi.c
	mpicc -o rs_mpi rs_mpi.c -lm  -Wfatal-errors -Werror  -march=native -ffast-math

openmp: rs_openmp.c
	gcc -o rs_openmp rs_openmp.c -lm  -g -Wfatal-errors -Werror -march=native -ffast-math -fopenmp 

mpi_clean: 
	rm -f rs_mpi

openmp_clean:
	rm -f rs_openmp
