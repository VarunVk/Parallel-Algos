#include <stdio.h> 		// For printf();
#include <stdlib.h>     // For exit();
#include <errno.h> 		// For errno for error handling
#include <string.h> 	// To include memset
#include <math.h> 		// To use sqrt()
#include "util.h"        // To use the timer functions
#include <omp.h>

/* GLOBAL DECLARATIONS */
#define NUM_DIGITS 16
int numData, numThreads;

//typdef unsigned int Data;

void printData(char *name, unsigned int *data, int numData);

void radix_sort_A(unsigned int *A, int numData)
{
    int b=32; 
    int r=2;
    int r2=pow(2,r);
    unsigned int MASK=pow(2,r)-1;
    unsigned int count[r2];
    unsigned int Pcount[r2];
    unsigned int *tmp_data;
    tmp_data=calloc(numData, sizeof(unsigned int));

    for(int i=0; i<(b/r2); i++)
    {
        memset(count, 0, sizeof(unsigned int)*r2);
//#pragma omp parallel for num_threads(numThreads) shared(numData, A,count) schedule(dynamic,1)
        for(int k=0; k<numData; k++)
        {
//#pragma omp atomic
            __sync_fetch_and_add(&count[A[k]&(MASK<<(i*r2))],1);
        }

        memset(Pcount, 0, sizeof(unsigned int)*r2);
        for(int j=0; j<r2-1; j++)
        {
            Pcount[j+1]= Pcount[j]+count[j];
            printf("%d\t", Pcount[j]);
            fflush(stdout);
        }
        printf("\tTotal %d\n", Pcount[r2-1]+count[r2-1]);

        memset(count, 0, sizeof(unsigned int)*r2);
        memset(tmp_data, 0, numData*sizeof(unsigned int));
//#pragma omp parallel for num_threads(numThreads) shared(A, tmp_data, count, Pcount) schedule(dynamic,1)
        for(int j=0; j<numData; j++)
        {
//#pragma omp crtical 
            tmp_data[Pcount[A[j]&(MASK<<(i*r2))]+__sync_fetch_and_add(&count[A[j]&(MASK<<(i*r2))],1)]=A[j];
        }
        print_numbers("A", A, numData);
        print_numbers("tmp_data", tmp_data, numData);
        memcpy(A, tmp_data, sizeof(unsigned int)*numData);
    }
    if(tmp_data)
        free(tmp_data);
}
#if 0  
void radix_sort_s(unsigned int *data, int numData)
{
    int m=1;
    unsigned int count[NUM_DIGITS]={0};
    unsigned int Pcount[NUM_DIGITS]={0};
    int MSB=0;

    printData("Input data", data, numData);

    unsigned int largest=0;
#pragma omp parallel for num_threads(numThreads) shared(data) schedule(dynamic,1)  reduction(max:largest) 
    for(unsigned int j=0; j<numData; j++)
    {
        if(data[j]>largest)
            largest=data[j];
    }
    printf("Largest number is %d.\n", largest);

    for(int i=(sizeof(unsigned int)*8)-1; i>0; i--)
        if(largest&(0x01<<i)){
            MSB=i;
            break;
        }
    printf("MSB is %d (0x%X)\n",MSB, largest);

    for(int i=0; i<=MSB;i++)
    {
        memset(count, 0, NUM_DIGITS*sizeof(unsigned int));
#pragma omp parallel for num_threads(numThreads) shared(data, count, base) schedule(dynamic,1)   
        for(unsigned int j=0; j<numData; j++)
        {
            for(int i=0; i<NUM_DIGITS; i++)
                base[i][j]=0;
        }

#pragma omp parallel for num_threads(numThreads) shared(data, count, base) schedule(dynamic,1)   
        for(unsigned int j=0; j<numData; j++)
        {
            if(((data[j]/m)%NUM_DIGITS)==0) {
#pragma omp critical (b0)
                base[0][count[0]++]= data[j];
            } else if(((data[j]/m)%NUM_DIGITS)==1) {
#pragma omp critical (b1)
                base[1][count[1]++]= data[j];
            } else if(((data[j]/m)%NUM_DIGITS)==2) {
#pragma omp critical (b2)
                base[2][count[2]++]= data[j];
            } else if(((data[j]/m)%NUM_DIGITS)==3) {
#pragma omp critical (b3)
                base[3][count[3]++]= data[j];
            } else if(((data[j]/m)%NUM_DIGITS)==4) {
#pragma omp critical (b4)
                base[4][count[4]++]= data[j];
            } else if(((data[j]/m)%NUM_DIGITS)==5) {
#pragma omp critical (b5)
                base[5][count[5]++]= data[j];
            } else if(((data[j]/m)%NUM_DIGITS)==6) {
#pragma omp critical (b6)
                base[6][count[6]++]= data[j];
            } else if(((data[j]/m)%NUM_DIGITS)==7) {
#pragma omp critical (b7)
                base[7][count[7]++]= data[j];
            } else if(((data[j]/m)%NUM_DIGITS)==8) {
#pragma omp critical (b8)
                base[8][count[8]++]= data[j];
            } else if(((data[j]/m)%NUM_DIGITS)==9) {
#pragma omp critical (b9)
                base[9][count[9]++]= data[j];
            } else if(((data[j]/m)%NUM_DIGITS)==10) {
#pragma omp critical (b10)
                base[10][count[10]++]= data[j];
            } else if(((data[j]/m)%NUM_DIGITS)==11) {
#pragma omp critical (b11)
                base[11][count[11]++]= data[j];
            } else if(((data[j]/m)%NUM_DIGITS)==12) {
#pragma omp critical (b12)
                base[12][count[12]++]= data[j];
            } else if(((data[j]/m)%NUM_DIGITS)==13) {
#pragma omp critical (b13)
                base[13][count[13]++]= data[j];
            } else if(((data[j]/m)%NUM_DIGITS)==14) {
#pragma omp critical (b14)
                base[14][count[14]++]= data[j];
            } else if(((data[j]/m)%NUM_DIGITS)==15) {
#pragma omp critical (b15)
                base[15][count[15]++]= data[j];
            }
        }
        unsigned int tot_count=0;
        for(int i=1; i<NUM_DIGITS; i++)
            Pcount[i] = Pcount[i-1]+ count[i-1];

        tot_count = Pcount[NUM_DIGITS-1]+count[NUM_DIGITS-1];
        if(tot_count != numData) {
            printf("There is a BUG in segregation(%d!=%d).\n", tot_count, numData);
        }

#if 0 
        unsigned int j=0;
#pragma omp parallel for num_threads(numThreads) shared(data, count, base) schedule(dynamic,1)   
        for(j=0; j<count[0]; j++)
                data[j]=base[0][j];
#pragma omp parallel for num_threads(numThreads) shared(data, count, base) schedule(dynamic,1)   
        for(j=count[0]; j<count[1]; j++)
                data[j]=base[1][j-count[0]];
#pragma omp parallel for num_threads(numThreads) shared(data, count, base) schedule(dynamic,1)   
        for(j=count[1]; j<count[2]; j++)
                data[j]=base[2][j-count[1]];
#pragma omp parallel for num_threads(numThreads) shared(data, count, base) schedule(dynamic,1)   
        for(j=count[2]; j<count[3]; j++)
                data[j]=base[3][j-count[2]];
#pragma omp parallel for num_threads(numThreads) shared(data, count, base) schedule(dynamic,1)   
        for(j=count[3]; j<count[4]; j++)
                data[j]=base[4][j-count[3]];
#pragma omp parallel for num_threads(numThreads) shared(data, count, base) schedule(dynamic,1)   
        for(j=count[4]; j<count[5]; j++)
                data[j]=base[5][j-count[4]];
#pragma omp parallel for num_threads(numThreads) shared(data, count, base) schedule(dynamic,1)   
        for(j=count[5]; j<count[6]; j++)
                data[j]=base[6][j-count[5]];
#pragma omp parallel for num_threads(numThreads) shared(data, count, base) schedule(dynamic,1)   
        for(j=count[6]; j<count[7]; j++)
                data[j]=base[7][j-count[6]];
#pragma omp parallel for num_threads(numThreads) shared(data, count, base) schedule(dynamic,1)   
        for(j=count[7]; j<count[8]; j++)
                data[j]=base[7][j-count[7]];
#pragma omp parallel for num_threads(numThreads) shared(data, count, base) schedule(dynamic,1)   
        for(j=count[8]; j<count[9]; j++)
                data[j]=base[1][j-count[8]];
#pragma omp parallel for num_threads(numThreads) shared(data, count, base) schedule(dynamic,1)   
        for(j=count[9]; j<count[10]; j++)
                data[j]=base[10][j-count[9]];
#pragma omp parallel for num_threads(numThreads) shared(data, count, base) schedule(dynamic,1)   
        for(j=count[10]; j<count[11]; j++)
                data[j]=base[11][j-count[10]];
#pragma omp parallel for num_threads(numThreads) shared(data, count, base) schedule(dynamic,1)   
        for(j=count[11]; j<count[12]; j++)
                data[j]=base[12][j-count[11]];
#pragma omp parallel for num_threads(numThreads) shared(data, count, base) schedule(dynamic,1)   
        for(j=count[12]; j<count[13]; j++)
                data[j]=base[13][j-count[12]];
#pragma omp parallel for num_threads(numThreads) shared(data, count, base) schedule(dynamic,1)   
        for(j=count[13]; j<count[14]; j++)
                data[j]=base[14][j-count[13]];
#pragma omp parallel for num_threads(numThreads) shared(data, count, base) schedule(dynamic,1)   
        for(j=count[14]; j<count[15]; j++)
                data[j]=base[15][j-count[14]];
#endif 

#pragma omp parallel for num_threads(numThreads) shared(data, count, base, Pcount) schedule(dynamic,1)   
        for(int i=0; i<NUM_DIGITS; i++)
            for(int j=0; j<count[i]; j++)
                data[Pcount[i]+j]=base[i][j];

        m*=NUM_DIGITS;
    }
}
#endif

void printData(char *name, unsigned int *data, int numData)
{
    return;
    printf("%s(%d): ", name, numData);
    //for(int j=0; j<numData; j++)
    //    printf("%d \t", data[j]);
    printf("\n");
}

void main(int argc, char **argv)
{
    if(argc-1 < 3 ){
        printf("Error : Please use ./rs_openmp <path to data file> <Number of threads> <path to o/p file>.\n");
        exit(0);
    }
    FILE *fData = fopen(argv[1],"r");
    numThreads=atoi(argv[2]);
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

    double start=monotonic_seconds();
    radix_sort_A(data, numData);
    double end=monotonic_seconds();

    print_time(end-start);

    print_numbers(argv[3], data, numData);

    if(data)
        free(data);
    if(fData)
        fclose(fData);
    if(fOutput)
        fclose(fOutput);
}
