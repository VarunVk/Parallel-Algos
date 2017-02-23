#include <stdio.h> 		// For printf();
#include <stdlib.h>     // For exit();
#include <pthread.h> 	// For pthread functions 
#include <errno.h> 		// For errno for error handling
#include <string.h> 	// To include memset
#include <math.h> 		// To use sqrt()
#include "tmr.h" 		// To use the timer functions 
#include <pthread.h>    // For pthread APIs

// Global declarations 
double *clusAvg;   // Holds the Cluster Averages - (Size = numClusters x Dim)
float *data; 	   // Holds the data - (Size= numDatapoints x Dim) 
int *cluster; 	   // Holds the index to which clsuter the data belongs to - (Size= numDataPoints x 1)
int numClusters, numThreads, dim,numDataPoints;
int *counter;
pthread_mutex_t *lockClusAvg;   // Mutex lock for clusAvg data structure 
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

// Return the cluster id of the closest cluster to data 
int closestClus(float *data, double *clusAvg)
{
	int idx=0;
	double closestDist=INFINITY, newDist;
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

// Get the closest cluster for a set of data 
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
				*(clusAvg+i*dim+j) /= *(counter+i); 
	}
	free(counter);
}

void *pThreadUpdateClusAvg(void* s)
{
	int i=*(int*)s;
	int idx,end;
	int blockSize = numDataPoints/numThreads;

	int start=i*blockSize;
	if(i == (numThreads-1))
		end = numDataPoints;
	else 
		end = (i+1)*blockSize;

	for(int i=start; i<end; i++)
	{
		idx= *(cluster+i);

		pthread_mutex_lock(lockClusAvg+idx);
		*(counter+idx)=*(counter+idx)+1;
		pthread_mutex_unlock(lockClusAvg+idx);

		pthread_mutex_lock(lockClusAvg+idx);
		for(int j=0; j<dim ;j++)
		{
			(*(clusAvg+idx*dim+j)) += (double)*(data+i*dim+j);
		}
		pthread_mutex_unlock(lockClusAvg+idx);
	}

}

// main Thread function to Re-assign the clusters for data points 
void *updateClusters(void *s)
{
	int idx= *(int *)s;
	int blockSize = numDataPoints/numThreads;
	if(idx!=(numThreads-1))
		*(int *)s= assignClus(idx*blockSize, (idx+1)*blockSize);
	else 
		*(int *)s= assignClus(idx*blockSize, numDataPoints);
	pthread_exit(s);
}

// Used for Unit Testing 
void writeToResult(double time)
{
	FILE* fResult= fopen("result.txt" , "a");
	if(fResult !=NULL) {
#ifdef USE_PTHREADS
		fprintf(fResult, "pThread C=%d Nthreads=%d time: %0.04fs\n",numClusters, numThreads, time);
#else
		fprintf(fResult, "Serial  C=%d Nthreads=%d time: %0.04fs\n",numClusters, numThreads, time);
#endif
		fclose(fResult);
	}
}

void kmeansPThreads()
{
	// pThread setup 
	pthread_t * tHandle= malloc(sizeof(pthread_t)*numThreads);
	memset(tHandle, 0, sizeof(pthread_t)*numThreads);

	int *ret_val, stable;

	int *idx= calloc(numThreads,sizeof(int));

	counter=malloc(sizeof(int)*numClusters);

	lockClusAvg =malloc(sizeof(pthread_mutex_t)*numClusters);
	memset(lockClusAvg, 0, sizeof(pthread_mutex_t)*numClusters);

	for(int i=0; i<numClusters;i++)
		pthread_mutex_init(lockClusAvg+i, NULL);

	double end, start=monotonic_seconds();
	for(int k=0; k< MAX_ATTEMPTS ; k++)
	{
		stable =NO_CHANGES;
		for(int i=0; i<numThreads; i++)
		{
			idx[i]=i;
			if(pthread_create(tHandle+i, NULL, updateClusters, (void *)(idx+i))==0) {
			}
		}
		for(int i=0; i<numThreads; i++)
			if(pthread_join(*(tHandle+i), (void **)&ret_val)==0) {
				if(*ret_val == CHANGED && stable==NO_CHANGES)
					stable=CHANGED;
			}
		if(stable == NO_CHANGES)
			goto done;

		memset(tHandle, 0, sizeof(pthread_t)*numThreads);
		memset(clusAvg, 0, sizeof(double)*numClusters*dim);
		memset(counter, 0, sizeof(int)*numClusters);

		// update the averages 
		for(int i=0; i<numThreads; i++)
		{
			idx[i]=i;
			if(pthread_create(tHandle+i, NULL, pThreadUpdateClusAvg, (void *)(idx+i))==0) {
			}
		}

		for(int i=0; i<numThreads; i++)
			if(pthread_join(*(tHandle+i), NULL)==0) {
			}

		for(int i=0; i<numClusters; i++)
		{
			if(*(counter+i))
				for(int j=0; j<dim; j++)
				{
					pthread_mutex_lock(lockClusAvg+i);
					*(clusAvg+i*dim+j) /= *(counter+i); 
					pthread_mutex_unlock(lockClusAvg+i);
				}
		}
	}
done:
	end=monotonic_seconds();
	print_time(end-start);
	//writeToResult(end-start);

	free(counter);
	free(lockClusAvg);

	if(tHandle)
		free(tHandle);
	if(idx)
		free(idx);
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
	
#ifdef USE_PTHREADS
	kmeansPThreads();
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

	// Include N and D in centroid.txt, this is only for plotting. 
	// fprintf(fCentroid,"%d %d\n", numDataPoints, dim);

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
