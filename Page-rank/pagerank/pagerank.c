
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <string.h>

#include "pr_graph.h"

/* GLOBAL DECLARATIONS */
int numThreads;
int myrank;

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

void write_results(char *name, double *PR, pr_int );

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

  pr_graph * graph; 
  graph = pr_graph_load(ifname);
  if(!graph) {
      return EXIT_FAILURE;
  }

  double * PR = pagerank(graph, 0.85, 100);

  char * ofname = NULL;
  if(argc > 2) 
      ofname = argv[2];

  /* write pagerank values */
  if(ofname) 
      write_results(ofname, PR, graph->nvtxs);
  free(PR);
  pr_graph_free(graph);
  MPI_Finalize();
  return EXIT_SUCCESS;
}



double * pagerank(
    pr_graph *graph,
    double const damping,
    int const max_iterations)
{
    int tot_nvtxs = graph->nvtxs*numThreads;
    int nnbrs = graph->xadj[graph->nvtxs];

    /* Initialize pageranks to be a probability distribution. */
    double * PR = malloc(graph->nvtxs * sizeof(*PR));
    double * PR_accum = malloc(graph->nvtxs * sizeof(*PR));
    for(pr_int v=0; v < graph->nvtxs; ++v) {
        PR[v] = 1. / (double) tot_nvtxs;
    }

    /* Probability of restart */
    double const restart = (1 - damping)/(double)tot_nvtxs;

    /* Convergence tolerance. */
    double const tol = 1e-9;

    /* Setup MPI RMA win for Inter rank communication*/
    MPI_Win win;

    /* int MPI_Win_create( void *base, MPI_Aint size, int disp_unit, MPI_Info info, MPI_Comm comm, MPI_Win *win); */
    MPI_Win_create(PR_accum, graph->nvtxs*sizeof(*PR_accum), sizeof(*PR_accum), MPI_INFO_NULL, MPI_COMM_WORLD, &win);

    /* int MPI_Win_fence( int assert, MPI_Win win); */
    MPI_Win_fence(0, win);

    /* Data structure to maintain the list of data outside our rank */ 
    pr_int out_nnbrs=0;                                    // Total numbers of out nbrs
    pr_int *out_nbrs_loc= malloc(nnbrs*sizeof(pr_int));    // Which rank it belongs to 
    double *out_nbrs_val= malloc(nnbrs*sizeof(double));    // Offset in that remote rank 
    memset(out_nbrs_loc, 0, nnbrs*sizeof(pr_int));

    double start;
    if(myrank==ROOT) {
        start=MPI_Wtime();
    }
    int i;
    /* Main loop */
    for(i=0; i < max_iterations; ++i) {
        memset(out_nbrs_val, 0, nnbrs*sizeof(double));
        out_nnbrs=0;

        for(pr_int v=0; v < graph->nvtxs; ++v) {
            PR_accum[v] = 0.;
        }

        /* Each vertex pushes PR contribution to all outgoing links */
        for(pr_int v=0; v < graph->nvtxs; ++v) {
            double const num_links = (double)(graph->xadj[v+1] - graph->xadj[v]);
            double const pushing_val = (num_links)?(PR[v]/num_links):0.;

            for(pr_int e=graph->xadj[v]; e < graph->xadj[v+1]; ++e) {
                // Check if the nbr is in our rank 
                if((graph->nbrs[e] >= myrank*graph->nvtxs) && (graph->nbrs[e]<(myrank+1)*graph->nvtxs)) {
                    PR_accum[(graph->nbrs[e])-(myrank*graph->nvtxs)] += pushing_val;
                
                // Save the remote rank and offset     
                } else { 
                    int dup=0;
                    // Accumulate the value if already present. 
                    for(pr_int x=0; x<out_nnbrs; x++)
                    {
                        if(out_nbrs_loc[x]==graph->nbrs[e]) {
                            out_nbrs_val[x] += pushing_val;
                            dup=1;
                            break;
                        }
                    }
                    // add the new value to out_nbrs 
                    if(dup==0) {
                        out_nbrs_loc[out_nnbrs] = graph->nbrs[e];
                        out_nbrs_val[out_nnbrs] += pushing_val;
                        out_nnbrs++;
                    }
                }
            }
        }

        /* Send out the Pushing_val to the outside nbrs */
        for(pr_int x=0; x<out_nnbrs; x++)
        {
            // MPI Unicast the value to the resp rank using MPI_Accumulate() RMA operation 
            int target_rank = (int )(out_nbrs_loc[x]/graph->nvtxs); 
            int target_disp = (out_nbrs_loc[x]%graph->nvtxs);
            MPI_Request req;
            MPI_Accumulate(&out_nbrs_val[x], 1, MPI_DOUBLE, target_rank, target_disp, 1, MPI_DOUBLE, MPI_SUM, win);
        }
        // Wait till RMA is complete 
        MPI_Win_fence(0, win);

        /* Finalize new PR values */
        double norm_changed = 0.;
        for(pr_int v=0; v < graph->nvtxs; ++v) {
            double const old = PR[v];
            PR[v] = restart + (damping * PR_accum[v]);
            norm_changed += (PR[v] - old) * (PR[v] - old);
        }

        // For MPI, lets get the sum(norm_changed) from all ranks and check if it less than tol
        // then terminate the processing .
        MPI_Allreduce(MPI_IN_PLACE, &norm_changed, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        norm_changed = sqrt(norm_changed);

        if(i > 1 && norm_changed < tol) {
            break;
        }
    }
    if(myrank==ROOT) {
        double seconds=MPI_Wtime()-start;
        printf("Number of iterations: %d average time %0.03fs\n",i, seconds/i);
    }

    free(PR_accum);
    free(out_nbrs_val);
    free(out_nbrs_loc);
    MPI_Win_free(&win);
    return PR;
}

void write_results(char *name, double *PR, pr_int size)
{
    MPI_Request req;
    int MAX_CHUNK=(size/CHUNK_SIZE);
    for(int i=0; i<MAX_CHUNK; i++)
        MPI_Isend(&PR[i*CHUNK_SIZE], CHUNK_SIZE, MPI_DOUBLE, ROOT, myrank*(MAX_CHUNK+1)+i, MPI_COMM_WORLD, &req);
    if(size%CHUNK_SIZE)
        MPI_Isend(&PR[MAX_CHUNK*CHUNK_SIZE], size%CHUNK_SIZE, MPI_DOUBLE, ROOT, myrank*(MAX_CHUNK+1)+MAX_CHUNK, MPI_COMM_WORLD, &req);


    if(myrank==ROOT) {
        double *TMP_PR=malloc(size*sizeof(*TMP_PR));
        FILE * fout = fopen(name, "w");
        if(!fout) {
            fprintf(stderr, "ERROR: could not open '%s' for writing.\n", name);
            return;
        }
        for(int tid=0; tid<numThreads; ++tid)
        {
            MPI_Request recv_req[MAX_CHUNK+1];
            MPI_Status   status;
            int chunk=0;
            for(; chunk<(size/CHUNK_SIZE); chunk++)
                MPI_Irecv(&TMP_PR[chunk*CHUNK_SIZE], CHUNK_SIZE, MPI_DOUBLE, tid, tid*(MAX_CHUNK+1)+chunk, MPI_COMM_WORLD, &recv_req[chunk]); 
            // Check size is not divisible by CHUNK_SIZE
            if(size%CHUNK_SIZE)
                MPI_Irecv(&TMP_PR[chunk*CHUNK_SIZE], size%CHUNK_SIZE, MPI_DOUBLE, tid, tid*(MAX_CHUNK+1)+MAX_CHUNK, MPI_COMM_WORLD, &recv_req[chunk]); 
            for(int j=0; j<MAX_CHUNK; j++) 
            {
                MPI_Wait(&recv_req[j], &status);
                for(int k=0; k<CHUNK_SIZE; ++k)
                    fprintf(fout, "%0.3e\n", TMP_PR[((j*CHUNK_SIZE)+k)]);
            }
            if(size%CHUNK_SIZE)
            {
                MPI_Wait(&recv_req[MAX_CHUNK], &status);
                for(int k=0; k<size%CHUNK_SIZE; ++k) {
                    fprintf(fout, "%0.3e\n", TMP_PR[((MAX_CHUNK*CHUNK_SIZE)+k)]);
                }
            }
        }
        free(TMP_PR);
        fclose(fout);
    }
}
