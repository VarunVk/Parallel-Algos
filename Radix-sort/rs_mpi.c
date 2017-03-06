#include <stdio.h> 		// For printf();
#include <stdlib.h>     // For exit();
#include <errno.h> 		// For errno for error handling
#include <string.h> 	// To include memset
#include <math.h> 		// To use sqrt()
#include "util.h"        // To use the timer functions
#include <omp.h>
#include <mpi.h>

/* GLOBAL DECLARATIONS */
#define NUM_DIGITS 4
#define ROOT 0
int numData, numThreads;

void convertToExclusiveScan(unsigned int *a, int size);
void printData(char *name, unsigned int *data, int numData);
void exclusiveScan(unsigned int *new, unsigned int *old, int size);

void radix_sort_m(int *data, int numData)
{
    int b=32; 
    int r=4;
    int bits=pow(2,r);
    unsigned int MASK=pow(2,r)-1;
    int digit=0;

    unsigned int Lcount[bits];
    unsigned int Scan_Lcount[bits];
    unsigned int Scan_Gcount[bits];

    unsigned int recvcounts[numThreads];
    unsigned int displs[numThreads];

    unsigned int *localData=calloc(numData,sizeof(unsigned int));
    unsigned int *sortedLocalData =calloc(numData,sizeof(unsigned int));

    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    char str[256];

    for(int i=0; i<(b/r); i++)
    {
        /* Scatter the data to various threads */
        MPI_Scatter(data, numData, MPI_INT, 
                localData, numData, MPI_INT, 0, MPI_COMM_WORLD);

        /* Reset count and read if data falls under that digit */
        memset(Lcount, 0, sizeof(unsigned int)*bits);
        for(int k=0; k<numData; k++)
            Lcount[(localData[k]>>(i*r))&MASK]++;

        /* Generate Scan for the count */ 
        memset(Scan_Lcount, 0, sizeof(unsigned int)*bits);
        exclusiveScan(Scan_Lcount, Lcount, bits);

        memset(Lcount, 0, sizeof(unsigned int)*bits);
        for(int k=0; k<numData; k++)
        {
            digit = (localData[k]>>(i*r)) & MASK; 
            sortedLocalData[Scan_Lcount[digit]+Lcount[digit]++]= localData[k];
        }

        if(MPI_SUCCESS != MPI_Allreduce(Scan_Lcount, Scan_Gcount, bits, MPI_INT, 
                    MPI_SUM, MPI_COMM_WORLD))
            printf("***************** We fucked!!!!!!!*************\n");

        for(int k=0; k<bits; k++){
            int tmp=0;
            memset(recvcounts, 0, sizeof(unsigned int)*numThreads);
            memset(displs, 0, sizeof(unsigned int)*numThreads);

            MPI_Gather(&Lcount[k], 1, MPI_INT, recvcounts, 
                    1, MPI_INT, ROOT, MPI_COMM_WORLD); 
            MPI_Scan(&Lcount[k], &tmp, 1, MPI_INT, 
                    MPI_SUM, MPI_COMM_WORLD); 

            MPI_Gather(&tmp, 1, MPI_INT, 
                    displs, 1, MPI_INT, ROOT,MPI_COMM_WORLD);

            if(myrank==0)
                convertToExclusiveScan(displs, numThreads);

            MPI_Gatherv(sortedLocalData+Scan_Lcount[k], Lcount[k], MPI_INT,
                    data+Scan_Gcount[k], recvcounts, displs, 
                    MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }
}

void exclusiveScan(unsigned int *res, unsigned int *data, int size)
{
    res[0]=0;
    for(int i=1; i<size; i++)
        res[i]=res[i-1]+data[i-1];
}

void convertToExclusiveScan(unsigned int *a, int size)
{
    for(int i=size-1; i>0; i--)
        a[i]=a[i-1];
    a[0]=0;
}

void printData(char *name, unsigned int *data, int numData)
{
    printf("%s(%d): ", name, numData);
    for(int j=0; j<numData; j++)
        printf("%10d ", data[j]);
    printf("\n");
}

void main(int argc, char **argv)
{
    double start;
    FILE *fData;
    unsigned int *data ;
    int myrank;
    MPI_Init(NULL,NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &numThreads);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    if(myrank==0) {
        if(argc-1 < 1 ){
            printf("Error : Please use ./rs_mpi <path to data file> <path to o/p file>.\n");
            MPI_Finalize();
            exit(0);
        }
        fData = fopen(argv[1],"r");

        if(fData== NULL) {
            printf("Error while opening the file %s.Errno = %d.\n", argv[1], errno);
            exit(1);
        } else {
            printf("Data file       = \"%s\".\n"
                    "Num of threads  = %d.\n" 
                    "Output file     = \"%s\".\n", argv[1], numThreads, argv[2]);
        }
        if(fscanf(fData,"%d", &numData) < 1) {
            printf("Error while getting the number of Data . errno= %d\n", errno);
            exit(1);
        }
        printf("Number of data points = %d \n", numData);


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
        if(fData)
            fclose(fData);
        //printf("\n*****************  Data elements *************************\n");
        //for(int i=0; i< numData; i++)
        //  printf("%u\t", data[i]);
        //printf("\n**********************************************************\n");

        start=monotonic_seconds();
    }
    MPI_Bcast(&numData, sizeof(numData), MPI_INT, 0, MPI_COMM_WORLD);

    //printf("T%d. NumThreads %d.\n", myrank,numThreads);
    radix_sort_m(data, numData/numThreads);

    if(myrank==0) {
        double end=monotonic_seconds();
        print_time(end-start);
        if(argv[2])
            print_numbers(argv[2], data, numData);
    }

    /*
    if(myrank==0) {
        printf("(%d) Im the Freeing Master!!!!\n", myrank);
        if(data)
            free(data);
    }
    */
exit:
    MPI_Finalize();
}
