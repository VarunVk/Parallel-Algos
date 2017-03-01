.PHONY: all pthread openmp clean

all: 

o: o_clean openmp

m: m_clean mpi

s: clean rs_openmp.c
	gcc -o rs_openmp rs_openmp.c -g -lm  -Wfatal-errors -Werror -march=native -ffast-math -O3
	
# PUt -O3
mpi: rs_mpi.c
	mpicc -o rs_mpi rs_mpi.c -lm  -Wfatal-errors -Werror  -march=native -ffast-math

openmp: rs_openmp.c
	gcc -o rs_openmp rs_openmp.c -lm  -g -Wfatal-errors -Werror  -O3 -march=native -ffast-math -fopenmp

m_clean: 
	rm -f rs_mpi

o_clean:
	rm -f rs_openmp
