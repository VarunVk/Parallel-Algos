#include <stdio.h> 		// For printf();
#include <stdlib.h>     // For exit();
#include <errno.h> 		// For errno for error handling
#include <string.h> 	// To include memset
#include <math.h> 		// To use sqrt()
#include "util.h"        // To use the timer functions
#include <omp.h>

/* GLOBAL DECLARATIONS */
#define NUM_DIGITS 16
#define TID omp_get_thread_num()
int numData, numThreads;

//typdef unsigned int Data;

void printData(char *name, unsigned int *data, int numData);

void radix_sort_A(unsigned int *A,unsigned int *tmp_data,  int numData)
{
    int b=32; 
    int r=8;
    int bits=pow(2,r);
    unsigned int MASK=pow(2,r)-1;
    unsigned int Pcount[bits];
    unsigned int count[numThreads][bits];

    for(int i=0; i<(b/r); i++)
    {
#pragma omp parallel for schedule(static)
        for(int j=0; j<numThreads; j++)
            for(int k=0; k<bits; k++)
                count[j][k]=0;

#pragma omp parallel for schedule(static)
        for(int k=0; k<numData; k++)
            count[TID][(A[k]>>(i*r))&MASK]++;

#pragma omp parallel for schedule(static)
        for(int k=0; k<bits; k++)
            Pcount[k]=0;

#pragma omp parallel for schedule(static)
        for(int j=0; j<bits; j++)
        {
            for(int k=0; k<numThreads; k++)
            {
                int tmp=count[k][j];
                count[k][j]=Pcount[j];
                Pcount[j] += tmp; 
            }
        }

        int prev=0;
        for(int j=0; j<bits; j++){
            int tmp=Pcount[j];
            Pcount[j]=prev;
            prev +=tmp;
        }

#pragma omp parallel for schedule(static)
        for(int k=0; k<numData; k++)
        {
            int pos = count[TID][(A[k]>>(i*r))&MASK]++;
            int start = Pcount[(A[k]>>(i*r))&MASK];
            tmp_data[start+pos]=A[k];
        }

#pragma omp parallel for schedule(static)
        for(int k=0; k<numData; k++)
            A[k]=tmp_data[k];
    }
}

void printData(char *name, unsigned int *data, int numData)
{
    printf("%10s(%10d): ", name, numData);
    for(int j=0; j<numData; j++)
        printf("%10d ", data[j]);
    printf("\n");
}

void main(int argc, char **argv)
{
    if(argc-1 < 2 ){
        printf("Error : Please use ./rs_openmp <path to data file> <Number of threads> <path to o/p file>.\n");
        exit(0);
    }
    FILE *fData = fopen(argv[1],"r");
    numThreads=atoi(argv[2]);
    omp_set_num_threads(numThreads);
    FILE *fOutput= fopen(argv[3],"w+");

    if(fData== NULL) {
        printf("Error while opening the file %s.Errno = %d.\n", argv[1], errno);
        exit(1);
    } else {
        printf("Data file       = \"%s\".\n"
                "Num of threads  = %d.\n" 
                "Output file     = \"%s\".\n", argv[1], numThreads, argv[3]);
    }
    if(fscanf(fData,"%d", &numData) < 1) {
        printf("Error while getting the number of Data . errno= %d\n", errno);
        exit(1);
    }
    printf("Number of data points = %d \n", numData);

    unsigned int *data ;

    if(numData*sizeof(unsigned int)> (100*1000*1000)) {
        printf("Input data size(%.4fMB) is greater than allowed(100MB).\n", ((double)numData*sizeof(unsigned int)/(1000*1000)));
        //exit(1);
    }

    data=calloc(numData, sizeof(unsigned int));

    if(data==NULL ){
        printf("Out of memory! while allocating mem for Data. errno %d \n",errno);
        exit(1);
    }

    // Read the data 
    for(int i=0; i< numData; i++) {
        if(fscanf(fData,"%u", &data[i]) != 1) {
            printf("Error while getting the (%d)data . errno= %d\n", i, errno);
            exit(0);
        }
    }
    //printf("\n*****************  Data elements *************************\n");
    //for(int i=0; i< numData; i++)
    //  printf("%u\t", data[i]);
    //printf("\n**********************************************************\n");

    unsigned int *tmp_data=calloc(numData, sizeof(unsigned int));
    double start=monotonic_seconds();
    radix_sort_A(data,tmp_data, numData);
    double end=monotonic_seconds();

    print_time(end-start);

    if(argv[3])
        print_numbers(argv[3], data, numData);

    if(tmp_data)
        free(tmp_data);
    if(data)
        free(data);
    if(fData)
        fclose(fData);
    if(fOutput)
        fclose(fOutput);
}
