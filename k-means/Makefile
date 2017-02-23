.PHONY: all pthread openmp clean

all: clean pthread openmp 

serial: km_threads.c
	gcc -o km_serial km_pthreads.c -lm  -Wfatal-errors -Werror -O3 -march=native -ffast-math
	
pthread: km_pthreads.c
	gcc -o km_pthreads km_pthreads.c -lm  -Wfatal-errors -Werror  -O3 -march=native -ffast-math -pthread -DUSE_PTHREADS

openmp: km_openmp.c
	gcc -o km_openmp km_openmp.c -lm  -Wfatal-errors -Werror  -O3 -march=native -ffast-math -fopenmp -DUSE_OPENMP

clean: 
	rm -f km_pthreads km_openmp
