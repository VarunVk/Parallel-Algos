#include <stdio.h>      // For printf();
#include <stdlib.h>     // For exit();
#include <errno.h>      // For errno for error handling
#include <string.h> 	// To include memset
#include <math.h>       // To use sqrt()
#include "util.h"       // To use the timer functions
#include <omp.h>

/* GLOBAL DECLARATIONS */
#define TID omp_get_thread_num()
int numData, numThreads;

unsigned int * radixSort(unsigned int *A,unsigned int *tmp_data,  int numData);

void main(int argc, char **argv)
{
    if(argc-1 < 3 ){
        printf("Error : Please use ./rs_openmp <path to data file> <Number of threads> <path to o/p file>.\n");
        exit(1);
    }
    FILE *fData = fopen(argv[1],"r");
    numThreads=atoi(argv[2]);
    omp_set_num_threads(numThreads);

    if(fData==NULL) {
        printf("Error while opening the file %s.Errno = %d.\n", argv[1], errno);
        exit(1);
    } 
    
    if(fscanf(fData,"%d", &numData) < 1) {
        printf("Error while getting the number of Data . errno= %d\n", errno);
        exit(1);
    }
   
    // Holds the original data
    unsigned int *data=calloc(numData, sizeof(unsigned int));
    // Used to hold temporary sorted data 
    unsigned int *tmp_data=calloc(numData, sizeof(unsigned int));
    
    if(data==NULL){
        printf("Out of memory! while allocating mem for Data. errno %d \n",errno);
        exit(1);
    }

    // Read the data from input file 
    for(int i=0; i< numData; i++) {
        if(fscanf(fData,"%u", &data[i]) != 1) {
            printf("Error while getting the (%d)data . errno= %d\n", i, errno);
            exit(0);
        }
    }

    double start=monotonic_seconds();
    data=radixSort(data,tmp_data, numData);
    double end=monotonic_seconds();

    print_time(end-start);

    print_numbers(argv[3], data, numData);

    if(tmp_data)
        free(tmp_data);
    if(data)
        free(data);
    if(fData)
        fclose(fData);
}

unsigned int * radixSort(unsigned int *data,unsigned int *tmp_data,  int numData)
{
    int b=32;       // unsigned int=32 bits
    int r=8;        // Choosing r-bit blocks per loop
    int numOfValues=pow(2,r);         // Range of values for r-bit blocks
    unsigned int MASK=pow(2,r)-1;     // Mask to get the r-bit block from the number

    /* Local count for the bits per thread*/
    unsigned int count[numThreads][numOfValues];

    /* Sum of the r-bit counts across threads */
    unsigned int Pcount[numOfValues];

    for(int i=0; i<(b/r); i++)
    {
      /* All OpenMP segments are executed which static scheduling so that
	 threads access the same elements to increase locality */
#pragma omp parallel for schedule(static)
        for(int j=0; j<numThreads; j++)
            for(int k=0; k<numOfValues; k++)
                count[j][k]=0;
	
	/* Count the number of data which have the same r-bit value */
#pragma omp parallel for schedule(static)
        for(int k=0; k<numData; k++)
            count[TID][(data[k]>>(i*r))&MASK]++;

#pragma omp parallel for schedule(static)
        for(int k=0; k<numOfValues; k++)
            Pcount[k]=0;
	
	/*This loop is to evaluate Prefix Scan for each r-bit value across
	  all threads. Calculation of each bit reduction is independent of each 
	  other, thus can be done in parallel.	*/
#pragma omp parallel for schedule(static)
        for(int j=0; j<numOfValues; j++)
        {
            for(int k=0; k<numThreads; k++)
            {
                int tmp=count[k][j];
                count[k][j]=Pcount[j];
                Pcount[j] += tmp; 
            }
        }

	/* Serially converting the Reduction the r-bit sums to Pcount, 
       as numOfValues is a small number  */
        int prev=0;
        for(int j=0; j<numOfValues; j++){
            int tmp=Pcount[j];
            Pcount[j]=prev;
            prev +=tmp;
        }

	/* Put the data back into a temporary data array, using the 
	   reduction variables calcluated from above.	 */
#pragma omp parallel for schedule(static)
        for(int k=0; k<numData; k++)
        {
            int pos = count[TID][(data[k]>>(i*r))&MASK]++;
            int start = Pcount[(data[k]>>(i*r))&MASK];
            tmp_data[start+pos]=data[k];
        }

	/* The temporary sorted data will now become the actual data as we 
	   swap the addresses */
	unsigned int *tmp=data;
	data=tmp_data;
	tmp_data=tmp;
    }
    return data;
}
