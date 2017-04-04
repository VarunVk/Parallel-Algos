

/* ensure we have `getline()` */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>

#include "pr_graph.h"

#define ROOT 0
extern int myrank;
int numThreads;

pr_graph * pr_graph_load(
    char const * const ifname)
{
    pr_graph * graph;
    FILE * fin; 
    if(myrank==ROOT) {

        fin = fopen(ifname, "r");
        if(!fin) {
            fprintf(stderr, "ERROR: could not open '%s' for reading.\n", ifname);
            return NULL;
        }

        graph = malloc(sizeof(*graph));

        /* read nvtxs and nedges */
        fscanf(fin, "%lu", &(graph->nvtxs));
        fscanf(fin, "%lu", &(graph->nedges));
        fscanf(fin, "\n"); /* make sure we process the newline, too. */
    }

    pr_graph * Lgraph = malloc(sizeof(*Lgraph));
    memset(Lgraph, 0, sizeof(*Lgraph));
    pr_int tot_nvtxs;

    if(myrank==ROOT) 
        tot_nvtxs=graph->nvtxs;

    MPI_Bcast(&tot_nvtxs, 1,MPI_UINT64_T, ROOT, MPI_COMM_WORLD);

    Lgraph->nvtxs=tot_nvtxs/numThreads;
    printf("(%d): Everyone gets %lu.\n", myrank, Lgraph->nvtxs);

    if(myrank==ROOT) {
        graph->xadj = malloc((Lgraph->nvtxs + 1) * sizeof(*graph->xadj));
        graph->nbrs = malloc(graph->nedges * sizeof(*graph->nbrs));
        if(graph->nbrs==NULL)
            printf("Error while allocating mem for graph->nbrs on ROOT.\n");
    }
    Lgraph->xadj = malloc((Lgraph->nvtxs + 1) * sizeof(*Lgraph->xadj));

    char * line = malloc(1024 * 1024);
    size_t len = 0;

    if(myrank==ROOT)
    {
        for(int x=0; x<numThreads ;x++)
        {
#if 1
            pr_int v;
            /* How many edges we have read. */
            pr_int edge_ptr = 0;

            /* Read in graph one vertex at a time. */
            for(v=0; v < Lgraph->nvtxs; ++v) {
                ssize_t read = getline(&line, &len, fin);
                if(read == -1) {
                    fprintf(stderr, "ERROR: premature EOF at line %lu\n", v+1);
                    pr_graph_free(graph);
                    return NULL;
                }

                /* Store the beginning of the adjacency list. */
                graph->xadj[v] = edge_ptr;

                /* Check for sinks -- these make pagerank more difficult. */
                if(read == 1) {
                    fprintf(stderr, "WARNING: vertex '%lu' is a sink vertex.\n", v+1);
                    continue;
                }

                /* Foreach edge in line. */
                char * ptr = strtok(line, " ");
                while(ptr != NULL) {
                    char * end = NULL;
                    pr_int const e_id = strtoull(ptr, &end, 10);
                    /* end of line */
                    if(ptr == end) {
                        break;
                    }
                    assert(e_id > 0 && e_id <= graph->nvtxs);

                    graph->nbrs[edge_ptr++] = e_id - 1; /* 1 indexed */
                    ptr = strtok(NULL, " ");
                }
            }
            //assert(edge_ptr == graph->nedges);
            graph->xadj[v] = edge_ptr;
#endif 

            MPI_Request send_req;
            MPI_Status   status;
            MPI_Isend(graph->xadj, Lgraph->nvtxs+1, MPI_UINT64_T, x, x, MPI_COMM_WORLD, &send_req); 
            /* Send edge_ptr and then the actual data */
            MPI_Wait(&send_req, &status);
            printf("edge_ptr %lu == nvtxs %lu\n", edge_ptr, graph->xadj[Lgraph->nvtxs]);

            MPI_Isend(graph->nbrs, edge_ptr, MPI_UINT64_T, x, 2*x, MPI_COMM_WORLD, &send_req); 
            MPI_Wait(&send_req, &status);
        } 
        /* MPI_Wait should come here  */
    } 
    {
        MPI_Request recv_req;
        MPI_Status   status;
        /* int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source,
           int tag, MPI_Comm comm, MPI_Request *request) */ 
        MPI_Irecv(Lgraph->xadj, Lgraph->nvtxs+1,MPI_UINT64_T, ROOT, myrank, MPI_COMM_WORLD, &recv_req); 
        MPI_Wait(&recv_req, &status);
        printf("(%d): Recd %lu xadj values.\n", myrank, Lgraph->nvtxs);
        for(pr_int y=0; y<Lgraph->nvtxs+1; y++)
            printf("(%d): xadj[%lu]=%lu .\t", myrank,y, Lgraph->xadj[y]);
        printf("\n");

        Lgraph->nbrs=malloc(Lgraph->xadj[Lgraph->nvtxs]*sizeof(*Lgraph->nbrs));
        MPI_Irecv(Lgraph->nbrs, Lgraph->xadj[Lgraph->nvtxs]+1,MPI_UINT64_T, ROOT, 2*myrank, MPI_COMM_WORLD, &recv_req); 
        MPI_Wait(&recv_req, &status);

        printf("(%d): Recd %lu nbr values.\n", myrank, Lgraph->nvtxs);
        for(pr_int y=0; y<Lgraph->xadj[Lgraph->nvtxs]; y++)
            printf("(%d): nbrs[%lu]=%lu .\t", myrank,y, Lgraph->nbrs[y]);
        printf("\n");
    }
    pr_graph_free(graph);
    free(line);
    return Lgraph;
}


void pr_graph_free(
        pr_graph * const graph)
{
  free(graph->xadj);
  free(graph->nbrs);
  free(graph);
}


