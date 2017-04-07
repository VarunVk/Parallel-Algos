

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

extern int myrank;
extern int numThreads;

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

    // Lgraph - holds the graph Locally on all ranks. 
    // graph  - holds the partial graph on ROOT, which is sent as CHUNKS
    pr_graph * Lgraph = malloc(sizeof(*Lgraph));
    memset(Lgraph, 0, sizeof(*Lgraph));
    pr_int tot_nvtxs;

    if(myrank==ROOT) 
        tot_nvtxs=graph->nvtxs;

    // Broadcast the total number of vertices
    MPI_Bcast(&tot_nvtxs, 1,MPI_UINT64_T, ROOT, MPI_COMM_WORLD);
    Lgraph->nvtxs=tot_nvtxs/numThreads;

    if(myrank==ROOT) {
        // Only allocate data required for one rank
        graph->xadj = malloc((Lgraph->nvtxs + 1) * sizeof(*graph->xadj));
        graph->nbrs = malloc(graph->nedges * sizeof(*graph->nbrs));
    }

    char * line = malloc(1024 * 1024);
    size_t len = 0;
    int MAX_CHUNKS=Lgraph->nvtxs/CHUNK_SIZE;

    // On ROOT takes care of reading and sending in CHUNK_SIZE to all ranks in parallel
    if(myrank==ROOT)
    {
        MPI_Request send_req;
        MPI_Status   status;
        for(pr_int x=0; x<numThreads ;x++)
        {
            pr_int v=0;
            /* How many edges we have read. */
            pr_int edge_ptr = 0;
            pr_int index=0;   // xadj index 
            pr_int old_end=0; // nbr index  

            /* Read in graph one vertex at a time. */
            for(v=0; v < Lgraph->nvtxs; ++v) {
                ssize_t read = getline(&line, &len, fin);
                if(read == -1) {
                    fprintf(stderr, "ERROR: premature EOF at line %lu\n", v+1);
                    pr_graph_free(graph);
                    return NULL;
                }

                /* Store the beginning of the adjacency list. */
                graph->xadj[index++] = edge_ptr;

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
                // Send out the data if we have reached CHUNK_SIZE or end of data for that rank 
                if((v+1)%CHUNK_SIZE==0 || ((v+1)==Lgraph->nvtxs)) {

                    int chunk_num= (int)(v/CHUNK_SIZE);
                    graph->xadj[index] = edge_ptr;
                    if(x==ROOT) {
                        if(old_end==0)
                            Lgraph->xadj = malloc((Lgraph->nvtxs + 1) * sizeof(*Lgraph->xadj));
                        memcpy(&Lgraph->xadj[CHUNK_SIZE*chunk_num], &graph->xadj[CHUNK_SIZE*chunk_num], (((v%CHUNK_SIZE)+1)+1)*sizeof(*Lgraph->xadj));

                        Lgraph->nbrs=(pr_int *)realloc(Lgraph->nbrs, edge_ptr*sizeof(*Lgraph->nbrs));
                        memcpy(&Lgraph->nbrs[old_end], &graph->nbrs[old_end], (edge_ptr-old_end)*sizeof(*Lgraph->nbrs)); 

                    } else {
                        // Send the xadj chunk in non blocking mode
                        MPI_Isend(&graph->xadj[CHUNK_SIZE*chunk_num], ((v%CHUNK_SIZE)+1)+1, MPI_UINT64_T, x, MAX_CHUNKS+chunk_num, MPI_COMM_WORLD, &send_req); 
                        // Send the number of nbrs and then the nbrs list
                        MPI_Send(&edge_ptr, 1, MPI_UINT64_T, x, 2*MAX_CHUNKS+chunk_num, MPI_COMM_WORLD);
                        MPI_Isend(&graph->nbrs[old_end], edge_ptr-old_end, MPI_UINT64_T, x, 3*MAX_CHUNKS+chunk_num, MPI_COMM_WORLD, &send_req);
                    }
                    old_end=edge_ptr;
                }
            }
            /* Before going to next rank, make sure the transfer is complete, so that we can reuse the buffer.*/
            if(x!=0)
                MPI_Wait(&send_req, &status);
        } 
    } else {
        // Other than ROOT, everyone waits for data from ROOT
        MPI_Request recv_req[MAX_CHUNKS];
        MPI_Status   status;

        Lgraph->xadj = malloc((Lgraph->nvtxs + 1) * sizeof(*Lgraph->xadj));
        pr_int old_size=0;
        pr_int cur_size=0;
        pr_int x;
        /* MPI_Recv for nbrs list */
        for(x=0; x< (Lgraph->nvtxs/CHUNK_SIZE); x++)
        {
            MPI_Recv(&Lgraph->xadj[x*CHUNK_SIZE], CHUNK_SIZE+1, MPI_UINT64_T, ROOT, MAX_CHUNKS+x, MPI_COMM_WORLD, &status /*, &recv_req[x]*/); 
            /* First recv the size of the chunk to be recvd */
            MPI_Recv(&cur_size, 1, MPI_UINT64_T, ROOT, 2*MAX_CHUNKS+x, MPI_COMM_WORLD, &status);
            Lgraph->nbrs=(pr_int *)realloc(Lgraph->nbrs, cur_size*sizeof(*Lgraph->nbrs));
            MPI_Recv(&Lgraph->nbrs[old_size], cur_size-old_size, MPI_UINT64_T, ROOT, 3*MAX_CHUNKS+x, MPI_COMM_WORLD , &status); 
            old_size=cur_size;
        }
        // If the size is not divisible by CHUNK_SIZE
        if(Lgraph->nvtxs%CHUNK_SIZE){
            MPI_Recv(&Lgraph->xadj[x*CHUNK_SIZE], (Lgraph->nvtxs%CHUNK_SIZE)+1, MPI_UINT64_T, ROOT, MAX_CHUNKS+x, MPI_COMM_WORLD, &status /*, &recv_req[x]*/); 
            /* First recv the size of the chunk to be recvd */
            MPI_Recv(&cur_size, 1, MPI_UINT64_T, ROOT, 2*MAX_CHUNKS+x, MPI_COMM_WORLD, &status);
            Lgraph->nbrs=(pr_int *)realloc(Lgraph->nbrs, cur_size*sizeof(*Lgraph->nbrs));
            MPI_Recv(&Lgraph->nbrs[old_size], cur_size-old_size, MPI_UINT64_T, ROOT, 3*MAX_CHUNKS+x, MPI_COMM_WORLD , &status); 
        } 
    }
    free(line);
    if(myrank==ROOT)
        pr_graph_free(graph);
    return Lgraph;
}

void pr_graph_free(
        pr_graph * const graph)
{
  free(graph->xadj);
  free(graph->nbrs);
  free(graph);
}
