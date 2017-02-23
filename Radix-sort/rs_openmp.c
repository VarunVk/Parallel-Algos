#include <stdio.h> 		// For printf();
#include <stdlib.h>     // For exit();
#include <errno.h> 		// For errno for error handling
#include <string.h> 	// To include memset
#include <math.h> 		// To use sqrt()
#include "util.h"        // To use the timer functions
#include <omp.h>

/* GLOBAL DECLARATIONS */
int numData, numThreads;

//typdef unsigned int Data;

unsigned int *data;
#define NUM_DIGITS 10
typedef struct Node {
	unsigned int data;
	struct Node *pNext; 
}Nd;
Nd bin[NUM_DIGITS];

void addData(Nd *N, unsigned int data)
{
	Nd *tmp= malloc(sizeof(Nd));
  if(tmp==NULL) {
    printf("Unable to allocate memory.\n");
    exit(1);
  }
  tmp->data=data;
	tmp->pNext=N->pNext;
	N->pNext=tmp;
}

void printBin();
void printData(unsigned int *D);
unsigned int* loadBin(Nd *tmp, unsigned int *D);
void radix_sort_s(unsigned int *data, int numData);
void radix_sort_s(unsigned int *data, int numData)
{
  int m=1,n=1;
  for(int i=0; i<10;i++)
    bin[i].pNext=NULL;

  for(int i=0; i<NUM_DIGITS;i++)
  {
    for(unsigned int j=0; j<numData; j++){
      addData(&bin[((*(data+j))/m)%10], *(data+j));
      //printf("%d ",j);
    }
    m*=10;
    printBin();
    //printData(data);
  }
}

void printData(unsigned int *D)
{
  printf("Data is : ");
  for(int i=0;i<numData;i++)
  {
    printf("%u\t",*D);
    D++;
  }
  printf("\n");
}

void printBin()
{
  Nd *tmp;
  unsigned int *D=data;
  for(int i=0; i<NUM_DIGITS; i++)
  {
    tmp=bin[i].pNext;
    printf("Bin[%d] = ", i);
    fflush(stdout);
    D=loadBin(tmp, D);
    bin[i].pNext=NULL;
    printf("\n");
  }
  printf("\n");
}

unsigned int* loadBin(Nd *tmp, unsigned int *D)
{
    //static unsigned int a=1;
    //printf("%u %p %u\n",a++, tmp->pNext, tmp->data);
    //fflush(stdout);
    if(tmp != NULL)
    {
      D=loadBin(tmp->pNext, D);
      printf("%u\t", tmp->data);
      *D=tmp->data;
      D++;
      free(tmp);
    }
    return D;
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

  if(numData*sizeof(unsigned int)> (100*1000*1000)) {
    printf("Input data size(%.4fMB) is greater than allowed(100MB).\n", ((double)numData*sizeof(unsigned int)/(1000*1000)));
    exit(1);
  }

  data= calloc(numData,sizeof(unsigned int));
	if(data==NULL){
		printf("Out of memory! while allocating mem for Data. errno %d \n",errno);
		exit(1);
	}

	// Read the data 
	for(int i=0; i< numData; i++) {
			if(fscanf(fData,"%u", data+i) != 1) {
				printf("Error while getting the (%d)data . errno= %d\n", i, errno);
				exit(0);
		}
	}
  printf("\n*****************  Data elements *************************\n");
	//for(int i=0; i< numData; i++)
  //  printf("%u\t", *(data+i));
  printf("\n**********************************************************\n");

  memset(bin, 0, sizeof(Nd)*NUM_DIGITS);

  double start=monotonic_seconds();
  radix_sort_s(data, numData);
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
