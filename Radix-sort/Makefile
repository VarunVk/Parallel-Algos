.PHONY: all pthread openmp clean

all: clean pthread openmp 

s: clean rs_openmp.c
	gcc -o rs_openmp rs_openmp.c -g -lm  -Wfatal-errors -Werror -march=native -ffast-math -O3
	
mpi: rs_mpi.c
	gcc -o rs_mpi rs_mpi.c -lm  -Wfatal-errors -Werror  -O3 -march=native -ffast-math -pthread -DUSE_PTHREADS

openmp: rs_openmp.c
	gcc -o rs_openmp rs_openmp.c -lm  -Wfatal-errors -Werror  -O3 -march=native -ffast-math -fopenmp -DUSE_OPENMP

clean: 
	rm -f rs_openmp rs_mpi
