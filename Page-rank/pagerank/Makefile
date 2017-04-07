pagerank: clean pagerank.c pr_graph.c
	mpicc -o pagerank *.c -lm  -march=native -ffast-math -O3 -std=c99

clean:
	rm -f pagerank
