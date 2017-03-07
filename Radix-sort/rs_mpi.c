#include <stdio.h> 		// For printf();
#include <stdlib.h>     // For exit();
#include <errno.h> 		// For errno for error handling
#include <string.h> 	// To include memset
#include "util.h"       // To use the timer functions
#include <math.h>       // To use pow()
#include <mpi.h>        // For MPI functionality 

/* GLOBAL DECLARATIONS */
#define ROOT 0
int numData;
int numThreads;
int myrank;

/* This fucntion converts a Inclusive Scan to a Exclusive by adding shifting the elements */
void convertToExclusiveScan(unsigned int *a, int size);

/* This function is used to compute exclusive Scan for an unsigne int 
   Array of size 'size'. */
void exclusiveScan(unsigned int *new, unsigned int *old, int size);

/* Core fucntion to sort the numbers using Radix Sort technique
   Inputs - pointer to data and number of data elements.  */ 
void radixSort(unsigned int *data, int numData); 

void main(int argc, char **argv)
{
    double start, end;
    FILE *fData=NULL;
    unsigned int *data ;
    int error=0;

    /* MPI Initializations */
    MPI_Init(NULL,NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &numThreads);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    /* Input routine to be done by single process i.e rank=0*/
    if(myrank==ROOT) {
        /*XXX:*/
        if(argc-1 < 1 ){
            printf("Error : Please use ./rs_mpi <path to data file> <path to o/p file>.\n");
            error=1;
        } else if ((fData=fopen(argv[1],"r")) ==NULL) {
            printf("Error while opening the file %s.Errno = %d.\n", argv[1], errno);
            error=2;
        } else if(fscanf(fData,"%d", &numData) < 1) {
            printf("Error while getting the number of Data . errno= %d\n", errno);
            error=3;
            /*XXX:*/
        } else if((data=calloc(numData, sizeof(unsigned int))) == NULL) {
            printf("Out of memory! while allocating mem for Data. errno %d \n",errno);
            error=4;
        } else {
            for(int i=0; i< numData; i++) {
                if(fscanf(fData,"%u", &data[i]) != 1) {
                    printf("Error while getting the (%d)data . errno= %d\n", i, errno);
                    error=5;
                }
            }
        }
        if(fData)
            fclose(fData);
    }

    /* Check if there was an error while inputing data,
       if Yes, fload it to all processes to exit. */
    MPI_Bcast(&error, 1,MPI_INT, ROOT, MPI_COMM_WORLD);
    if(error!=0) {
        printf("Errorcode %d while inputing data. Exting process(%d)\n", error, myrank);
        MPI_Finalize();
        return;
    }

    /* Fload the Number of Data (i.e. variable 'numData') to all processes */
    MPI_Bcast(&numData, sizeof(numData), MPI_INT, ROOT, MPI_COMM_WORLD);

    /* Only the ROOT monitors the time */
    if(myrank==0)
        start=monotonic_seconds();

    /* Divide the Data by number of processes */ 
    radixSort(data, numData/numThreads);

    MPI_Barrier(MPI_COMM_WORLD);
    if(myrank==0) {
        end=monotonic_seconds();
        print_time(end-start);
        /* Ouput to file done by the ROOT process */ 
        /*XXX:*/ 
        if(argv[2])
            print_numbers(argv[2], data, numData);
    }
    MPI_Finalize();
}

void radixSort(unsigned int *data, int numData)
{
    int b=32;       // unsigned int=32 bits
    int r=8;        // r-bit blocks per loop
    int numOfValues=pow(2,r);   // Range of values for r-bit blocks
    unsigned int MASK=pow(2,r)-1;     // Mask to get the r-bit block from the number

    /* Local count of the r-bit value of the data */
    unsigned int Lcount[numOfValues];

    /* Local Scan count of the r-bit value of the data */
    unsigned int Scan_Lcount[numOfValues];

    /* Global Scan count of the r-bit value of the data. Evaluated only at the ROOT*/
    unsigned int Scan_Gcount[numOfValues];

    /*This 2 variables while sending data back to ROOT using MPI_Gatherv
      holding the counts of the array elements and the displacements */
    unsigned int recvcounts[numThreads];
    unsigned int displs[numThreads];

    /* This holds the unsorted data received from ROOT */
    unsigned int *localData=calloc(numData,sizeof(unsigned int));

    /* This holds the sorted data of the above array*/
    unsigned int *sortedLocalData =calloc(numData,sizeof(unsigned int));

    for(int i=0; i<(b/r); i++)
    {
        /* Scatter the data to various processes from ROOT */
        MPI_Scatter(data, numData, MPI_INT, 
                localData, numData, MPI_INT, 0, MPI_COMM_WORLD);

        /* Reset Local count and increment if data falls under that r-bit value */
        memset(Lcount, 0, sizeof(unsigned int)*numOfValues);
        for(int k=0; k<numData; k++)
            Lcount[(localData[k]>>(i*r))&MASK]++;

        /* Generate Exclusive Scan for the local count to put in the local sorted Array */ 
        memset(Scan_Lcount, 0, sizeof(unsigned int)*numOfValues);
        exclusiveScan(Scan_Lcount, Lcount, numOfValues);

        /* Put the unsorted data into the sorted list using Scan offsets */
        memset(Lcount, 0, sizeof(unsigned int)*numOfValues);
        for(int k=0; k<numData; k++)
        {
            int r_bitValue= (localData[k]>>(i*r)) & MASK; 
            sortedLocalData[Scan_Lcount[r_bitValue]+Lcount[r_bitValue]++]= localData[k];
        }

        /* Now to put back the values into the global data in ROOT 
           we have to compute the Reduction for the Local Scan counts */
        if(MPI_SUCCESS != MPI_Reduce(Scan_Lcount, Scan_Gcount, numOfValues, MPI_INT, 
                    MPI_SUM, ROOT, MPI_COMM_WORLD))
            printf("(P%d): Unable to reduce data.\n", myrank);

        /* This is the for loop responsible to transfer each r-bit value from local sorted array
           global data in ROOT using MPI_Gatherv*/
        for(int k=0; k<numOfValues; k++){
            int local_sum=0;
            memset(recvcounts, 0, sizeof(unsigned int)*numThreads);
            memset(displs, 0, sizeof(unsigned int)*numThreads);

            MPI_Gather(&Lcount[k], 1, MPI_INT, recvcounts, 
                    1, MPI_INT, ROOT, MPI_COMM_WORLD); 
            MPI_Scan(&Lcount[k], &local_sum, 1, MPI_INT, 
                    MPI_SUM, MPI_COMM_WORLD); 

            MPI_Gather(&local_sum, 1, MPI_INT, 
                    displs, 1, MPI_INT, ROOT,MPI_COMM_WORLD);

            /* MPI_Scan provides inclusive Scan convert that to Exclusive Scan */
            if(myrank==ROOT)
                convertToExclusiveScan(displs, numThreads);

            MPI_Gatherv(sortedLocalData+Scan_Lcount[k], Lcount[k], MPI_INT,
                    data+Scan_Gcount[k], recvcounts, displs, 
                    MPI_INT, ROOT, MPI_COMM_WORLD);
        }
    }
}

void exclusiveScan(unsigned int *res, unsigned int *data, int size)
{
    res[0]=0;
    for(int i=1; i<size; i++)
        res[i]=res[i-1]+data[i-1];
}

/* Conversion is easily done by shifting the data by one place */ 
void convertToExclusiveScan(unsigned int *a, int size)
{
    for(int i=size-1; i>0; i--)
        a[i]=a[i-1];
    a[0]=0;
}
