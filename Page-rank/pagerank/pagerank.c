
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <string.h>

#include "pr_graph.h"


/* GLOBAL DECLARATIONS */
#define ROOT 0
int numThreads;
int myrank;
pr_graph * graph; 

/**
* @brief Compute the PageRank (PR) of a graph.
*
* @param graph The graph.
* @param damping Damping factor (or, 1-restart). 0.85 is typical.
* @param max_iterations The maximium number of iterations to perform.
*
* @return A vector of PR values.
*/
double * pagerank(
    pr_graph * graph,
    double const damping,
    int const max_iterations);


int main(
    int argc,
    char * * argv)
{
  if(argc == 1) {
    fprintf(stderr, "usage: %s <graph> [output file]\n", *argv);
    return EXIT_FAILURE;
  }

  /* MPI Initializations */
  MPI_Init(NULL,NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &numThreads);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);


  char * ifname = argv[1];

  graph = pr_graph_load(ifname);
  if(!graph) {
      return EXIT_FAILURE;
  }

  double * PR = pagerank(graph, 0.85, 100); 

#if 0
  if(myrank==ROOT) {
      printf("Writing to output file .\n");
      char * ofname = NULL;
      if(argc > 2) {
          ofname = argv[2];
      }
      /* write pagerank values */
      if(ofname) {
          FILE * fout = fopen(ofname, "w");
          if(!fout) {
              fprintf(stderr, "ERROR: could not open '%s' for writing.\n", ofname);
              return EXIT_FAILURE;
          }
          for(pr_int v=0; v < graph->nvtxs; ++v) {
              fprintf(fout, "%0.3e\n", PR[v]);
          }
          fclose(fout);
      }
      free(PR);
  }
#endif
  MPI_Finalize();
  return EXIT_SUCCESS;
}



double * pagerank(
    pr_graph *Rgraph,
    double const damping,
    int const max_iterations)
{
  pr_graph * Lgraph = malloc(sizeof(*graph));
  memset(Lgraph, 0, sizeof(*graph));
  pr_int tot_nvtxs;

  if(myrank==ROOT) 
      tot_nvtxs=Rgraph->nvtxs;
  MPI_Bcast(&tot_nvtxs, 1,MPI_UINT64_T, ROOT, MPI_COMM_WORLD);
  Lgraph->nvtxs=tot_nvtxs/numThreads;


  // Scatter the adj data 
  Lgraph->xadj = malloc((Lgraph->nvtxs+1)*sizeof(*Lgraph->xadj));
  if(Lgraph->xadj==NULL)
  {
      printf("Unable to allocate mem for xadj");
      return NULL;
  }
  memset(Lgraph->xadj, 0, (Lgraph->nvtxs+1)*sizeof(*Lgraph->xadj));

  // Scatter the xadj to all ranks 
  int sendcounts[numThreads];
  int displs[numThreads];
  for(int i=0; i<numThreads; i++)
  {
      sendcounts[i]=Lgraph->nvtxs+1;
      displs[i]=(i*Lgraph->nvtxs);
  }
  void *send;
  if(myrank==ROOT)
      send=Rgraph->xadj;
  else 
      send=NULL;

  if(MPI_SUCCESS != MPI_Scatterv(send, sendcounts, displs, MPI_UINT64_T, 
                                Lgraph->xadj, Lgraph->nvtxs+1, MPI_UINT64_T,
                                ROOT, MPI_COMM_WORLD)) {
      printf("Scattering the xadj failed.!!!\n");
      return NULL;
  }
  // Time to scatter nbrs to all ranks 
  int nnbrs=(Lgraph->xadj[Lgraph->nvtxs]-Lgraph->xadj[0])/*-1 Should we subtract one here*/;
  Lgraph->nbrs = malloc(nnbrs*sizeof(*Lgraph->nbrs));
  if(Lgraph->nbrs==NULL)
  {
      printf("Unable to allocate mem for nbrs");
      return NULL;
  }
  memset(Lgraph->nbrs, 0, nnbrs * sizeof(*Lgraph->nbrs));


  if(myrank==ROOT){
      for(int i=0; i<numThreads; i++)
      {
          sendcounts[i]=Rgraph->xadj[(i+1)*Lgraph->nvtxs]- Rgraph->xadj[i*Lgraph->nvtxs];
          displs[i]=Rgraph->xadj[i*Lgraph->nvtxs];
      }
  }

  if(myrank==ROOT)
      send=Rgraph->nbrs;
  else 
      send=NULL;
  if(MPI_SUCCESS != MPI_Scatterv(send, sendcounts, displs, MPI_UINT64_T, 
                                Lgraph->nbrs, nnbrs, MPI_UINT64_T,
                                ROOT, MPI_COMM_WORLD)) {
      printf("Scattering the xadj failed.!!!\n");
      return NULL;
  }

  /* Initialize pageranks to be a probability distribution. */
  double * PR = malloc(Lgraph->nvtxs * sizeof(*PR));
  for(pr_int v=0; v < Lgraph->nvtxs; ++v) {
    PR[v] = 1. / (double) tot_nvtxs;
  }

  /* Probability of restart */
  double const restart = (1 - damping) / (double) tot_nvtxs;


  /* Convergence tolerance. */
  double const tol = 1e-12;

  double * PR_accum = malloc(Lgraph->nvtxs * sizeof(*PR));

  /* Setup MPI win for Inter process communication*/
  MPI_Win win;

  /* int MPI_Win_create( void *base, MPI_Aint size, int disp_unit, MPI_Info info, MPI_Comm comm, MPI_Win *win); */
  MPI_Win_create(PR_accum, Lgraph->nvtxs * sizeof(*PR), sizeof(*PR), MPI_INFO_NULL, MPI_COMM_WORLD, &win);

  /* int MPI_Win_fence( int assert, MPI_Win win); */
  MPI_Win_fence(0, win);

  pr_int out_nnbrs=0;
  pr_int out_nbrs_loc[nnbrs];
  double out_nbrs_val[nnbrs];
  memset(out_nbrs_loc, 0, nnbrs*sizeof(pr_int));

  /* Main loop */
  for(int i=0; i < 1 /*max_iterations*/; ++i) {
      memset(out_nbrs_val, 0, nnbrs*sizeof(double));

    for(pr_int v=0; v < Lgraph->nvtxs; ++v) {
      PR_accum[v] = 0.;
    }

    /* Each vertex pushes PR contribution to all outgoing links */
    for(pr_int v=0; v < Lgraph->nvtxs; ++v) {
      double const num_links = (double)(Lgraph->xadj[v+1] - Lgraph->xadj[v]);
      double const pushing_val = PR[v] / num_links;

      for(pr_int e=Lgraph->xadj[v]; e < Lgraph->xadj[v+1]; ++e) {
          // Check if the nbr is in our rank 
          if((Lgraph->nbrs[e]> myrank*Lgraph->nvtxs) && (Lgraph->nbrs[e]<(myrank+1)*Lgraph->nvtxs))
              PR_accum[Lgraph->nbrs[e]] += pushing_val;
          else { 
              continue;
              int dup=0;
              for(pr_int x=0; x<out_nnbrs; x++)
              {
                  if(out_nbrs_loc[x]==Lgraph->nbrs[e]) {
                      out_nbrs_val[x] += pushing_val;
                      out_nnbrs++;
                      dup=1;
                      break;
                  }
              }
              if(dup==0) {
                  out_nbrs_loc[out_nnbrs] = Lgraph->nbrs[e];
                  out_nbrs_val[out_nnbrs] += pushing_val;
                  out_nnbrs++;
              }

          }
      }
    }
    /* Send out the Pushing_val to the outside nbrs */

    for(pr_int x=0; x<out_nnbrs; x++)
    {
        // MPI Unicast the value to the resp rank 
        int target_rank = (int )(out_nbrs_loc[x]/Lgraph->nvtxs); 
        int target_disp = (out_nbrs_loc[x]%Lgraph->nvtxs);
        MPI_Accumulate(&out_nbrs_val[x], 1, MPI_DOUBLE, target_rank, target_disp, 1, MPI_DOUBLE, MPI_SUM, win);
        MPI_Win_fence(0, win);
    }
    /* Finalize new PR values */
    double norm_changed = 0.;
    for(pr_int v=0; v < Lgraph->nvtxs; ++v) {
      double const old = PR[v];
      PR[v] = restart + (damping * PR_accum[v]);

      norm_changed += (PR[v] - old) * (PR[v] - old);
    }
    norm_changed = sqrt(norm_changed);

    // For MPI, lets get the max(norm_changed) from all ranks and check if it less than tol
    // then terminate the processing .
    MPI_Allreduce(MPI_IN_PLACE, &norm_changed, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    if(i > 1 && norm_changed < tol) {
      break;
    }
    printf("%d. %f\n", myrank, norm_changed);
    MPI_Barrier(MPI_COMM_WORLD);
  }

  //free(PR_accum);
  /* Only on ROOT gather all the PR values. */
  double *G_PR;
  if(myrank == ROOT) {
        G_PR=malloc(Rgraph->nvtxs*sizeof(*PR));
  }
  MPI_Gather(PR, Lgraph->nvtxs, MPI_DOUBLE, G_PR, Lgraph->nvtxs, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
  //free(PR);
  MPI_Win_free(&win);
  printf("%d I'm done. \n", myrank);
  MPI_Barrier(MPI_COMM_WORLD);
  return G_PR;
}
