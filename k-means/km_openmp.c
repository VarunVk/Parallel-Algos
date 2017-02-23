#include <stdio.h> 		// For printf();
#include <stdlib.h>     // For exit();
#include <pthread.h> 	// For pthread functions 
#include <errno.h> 		// For errno for error handling
#include <string.h> 	// To include memset
#include <math.h> 		// To use sqrt()
#include "tmr.h"        // To use the timer functions
#include <omp.h>

// Global declarations 
double *clusAvg;   // Holds the Cluster Averages - (Size = numClusters x Dim)
float *data; 	   // Holds the data - (Size= numDatapoints x Dim) 
int *cluster; 	   // Holds the index to which clsuter the data belongs to - (Size= numDataPoints x 1)
int numClusters, numThreads, dim,numDataPoints;
int *counter;
typedef enum {CHANGED=0, NO_CHANGES=1} stabilty;

#define CLUSTERS_FILE  "clusters.txt"
#define CENTROID_FILE  "centroids.txt"
#define MAX_ATTEMPTS   20

// Calculate the euclidean distance 
double dist(double *C, float *data)
{
	double sum=0.0;
	for(int i=0;i<dim;i++)
		sum+=(*(C+i)-*(data+i))*(*(C+i)-*(data+i));
	return sqrt(sum);
}

// Return the closet cluster id to data 
int closestClus(float *data, double *clusAvg)
{
	int idx=0;
	double closestDist=INFINITY ,newDist;
	for(int i=0;i<numClusters;i++)
	{
		newDist = dist(clusAvg+i*dim, data);
		if(newDist < closestDist){
			idx=i;
			closestDist = newDist;
		}
	}
	return idx;
}

int assignClus(int startIdx, int endIdx)
{
	int stable=NO_CHANGES, clusIdx=0;
	for(int i= startIdx; i<endIdx;i++)
	{
		clusIdx = closestClus(data+i*dim, clusAvg);
		if(*(cluster+i) != clusIdx) {
			// Mark as unstable
			stable = CHANGED;
			*(cluster+i) = clusIdx;
		}
	}
	return stable;
}

// Initialising the fisrt 'n'  data as cluster average, where n=numClusters
void init_clusAvg()
{
	for(int i=0; i<numClusters; i++)
		for(int j=0; j<dim; j++)
			*(clusAvg+i*dim+j)= (double)*(data+i*dim+j);
}

void updateClusAvg(double *clusAvg,float* data)
{
	int idx=0;
	memset(clusAvg, 0, sizeof(double)*numClusters*dim);

	for(int i=0; i<numDataPoints; i++)
	{
		idx= *(cluster+i);
		*(counter+idx)=*(counter+idx)+1;
		for(int j=0; j<dim ;j++)
		{
			*(clusAvg+idx*dim+j) += (double)*(data+i*dim+j);
		}
	}
	for(int i=0; i<numClusters; i++)
	{
		if(*(counter+i))
			for(int j=0; j<dim; j++)
			{
				*(clusAvg+i*dim+j) /= *(counter+i); 
			}
	}
	free(counter);
}

void writeToResult(double time)
{
	FILE* fResult= fopen("result.txt" , "a");
	if(fResult !=NULL) {
#ifdef USE_OPENMP
		fprintf(fResult, "OpenMP  C=%d Nthreads=%d time: %0.04fs\n",numClusters, numThreads, time);
#else
		fprintf(fResult, "Serial  C=%d Nthreads=%d time: %0.04fs\n",numClusters, numThreads, time);
#endif
		fclose(fResult);
	}
}

void kmeansOpenMP()
{
	int *counter=malloc(sizeof(int)*numClusters);
	int clusIdx, stable=NO_CHANGES;
	double end, start=monotonic_seconds();
	for(int k=0; k< MAX_ATTEMPTS; k++)
	{
#pragma omp parallel for num_threads(numThreads) private(clusIdx) shared(numDataPoints, cluster, data, clusAvg) schedule(dynamic, 1) reduction(*:stable)
		for(int j= 0;j< numDataPoints;j++)
		{
			clusIdx =closestClus(data+j*dim, clusAvg);
			if(clusIdx != *(cluster+j)) {
				// Mark as unstable
				stable = CHANGED;
				*(cluster+j) = clusIdx;
			}
		}
		if(stable == NO_CHANGES)
			goto done;

		int idx=0;
		memset(clusAvg, 0, sizeof(double)*numClusters*dim);
		memset(counter, 0, sizeof(int)*numClusters);

#pragma omp parallel for num_threads(numThreads) shared(numDataPoints, cluster, data, clusAvg) schedule(dynamic, 1)
		for(int i=0; i<numDataPoints; i++)
		{
			idx= *(cluster+i);
			__sync_fetch_and_add(counter+idx,1);

			for(int j=0; j<dim ;j++)
			{
#pragma omp atomic
				*(clusAvg+idx*dim+j) += (double)*(data+i*dim+j);
			}
		}

#pragma omp parallel for num_threads(numThreads) shared(numDataPoints, cluster, data, clusAvg) schedule(dynamic, 1)
		for(int i=0; i<numClusters; i++)
		{
			if(*(counter+i))
				for(int j=0; j<dim; j++)
				{
					*(clusAvg+i*dim+j) /= *(counter+i); 
				}
		}
	}
done: 
	end=monotonic_seconds();
	free(counter);
	print_time(end-start);
	//writeToResult(end-start);
}

void main(int argc, char **argv)
{
	FILE *fDataPoints=NULL;

	if(argc-1 < 3 ){
		printf("Error : Please use ./km_pthreads <path to data file> <Number of clusters> <Number of threads>.\n");
		exit(0);
	}

	fDataPoints = fopen(argv[1],"r");
	numClusters=atoi(argv[2]);
	numThreads=atoi(argv[3]);

	if(fDataPoints == NULL) {
		printf("Error while opening the file %s.Errno = %d.\n", argv[1], errno);
		exit(1);
	} else {
		//printf("Data file       = \"%s\".\n"
	 	//  		"Num of clusters = %d.\n"
		//		"Num of threads  = %d.\n" , argv[1], numClusters, numThreads);
	}

	if(fscanf(fDataPoints,"%d %d", &numDataPoints, &dim) < 2) {
		printf("Error while getting the dimension and number of Data points. errno= %d\n", errno);
	}
	//printf("Number of data points = %d \n"
	//		"Dimension             = %d \n",numDataPoints, dim);

	// Declare and initialise the cluster list 
	cluster= malloc(sizeof(int)*numDataPoints);
	if(cluster==NULL){
		printf("Out of memory! While allocating memory for clusters.\n");
		exit(1);
	}
	memset(cluster, 0, sizeof(int)*numDataPoints);

	data=malloc(sizeof(float)*numDataPoints*dim);
	if(data==NULL){
		printf("Out of memory! While allocating mem for data. errno %d\n",errno);
		exit(1);
	}
	memset(data,0,sizeof(float)* numDataPoints* dim);

	// Initialize the clusters with the first numClusters 
	clusAvg= malloc(sizeof(double)*numClusters*dim);
	if(clusAvg==NULL){
		printf("Out of memory! While allocating mem for cluster Averages. Errno %d \n",errno);
		exit(1);
	}
	memset(clusAvg, 0, sizeof(double)*numClusters*dim);

	// Read the data 
	for(int i=0; i< numDataPoints; i++) {
		for(int j=0; j<dim; j++) {
			if(fscanf(fDataPoints,"%f ", data+i*dim+j) != 1) {
				printf("Error while getting the (%d)data . errno= %d\n", i, errno);
				exit(0);
			}
		}
	}
	fclose(fDataPoints);
	// Assign to the respective clusters
	init_clusAvg();
	
#ifdef USE_OPENMP
	kmeansOpenMP();
#else
	double start=monotonic_seconds();
	for(int i=0;i<20;i++)
	{
		if(1 == assignClus(0, numDataPoints))
			break;
		updateClusAvg(clusAvg, data);
	}
	double end=monotonic_seconds();
	print_time(end-start);
	//writeToResult(end-start);

#endif 

	// write results to output files 
	FILE *fCluster= fopen(CLUSTERS_FILE, "w+");
	FILE *fCentroid= fopen(CENTROID_FILE, "w+");

	fprintf(fCentroid,"%d %d\n", numDataPoints, dim);
	for(int i=0; i<numClusters; i++)
	{
		for(int j=0; j<dim; j++)
			fprintf(fCentroid,"%.10e ", *(clusAvg+i*dim+j) );
		fprintf(fCentroid,"\n");
	}

	for(int i=0; i<numDataPoints; i++)
		fprintf(fCluster,"%d\n", *(cluster+i));
		
	// Cleanup code - Free all allocated memory 
	if(fCluster)
		fclose(fCluster);
	if(fCentroid)
		fclose(fCentroid);
	if(clusAvg)
		free(clusAvg);
	if(cluster)
		free(cluster);
	if(data)
		free(data);
}
