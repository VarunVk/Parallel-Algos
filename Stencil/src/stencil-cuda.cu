#include <stdlib.h>
#include <stdio.h>

extern "C"
{
#include "image.h"
#include "stencil.h"
}

#define TILE_SIZE 16
#define BLOCK_SIZE 16

__constant__ float  G_St[3][3];

__global__ void stencil_kernel (int *h, int *w, float *St,
        float *InR, float *InG, float *InB,
        float *OutR, float *OutG, float *OutB)
{
    int tx=threadIdx.x;
    int ty=threadIdx.y;
    int x=blockIdx.x*BLOCK_SIZE+tx;
    int y=blockIdx.y*BLOCK_SIZE+ty;

#if 0    // Simple implementation with Constant memory - Working like a Don!!
    if(x<*w && y<*h) {
        float R=0.,G=0.,B=0.;
        for(int i=0; i<3; ++i)
        {
            for(int j=0; j<3; ++j)
            {
                int const xpos = x+(i-1);
                int const ypos = y+(j-1);
                //float const x = *(St+3*i+j);
                float const x = G_St[i][j];
                R += *(InR+(ypos**w+xpos)) * x;
                G += *(InG+(ypos**w+xpos)) * x;
                B += *(InB+(ypos**w+xpos)) * x;
            }
        }
        *(OutR+(y**w+x)) = R;
        *(OutG+(y**w+x)) = G;
        *(OutB+(y**w+x)) = B;
    }
#elif 1
    __shared__ float sR[BLOCK_SIZE+2][BLOCK_SIZE+2];
    __shared__ float sG[BLOCK_SIZE+2][BLOCK_SIZE+2];
    __shared__ float sB[BLOCK_SIZE+2][BLOCK_SIZE+2];

    // Adjust the co-ordinates
    x -=1;    y -=1;
    // Put Input data in shared memory
    if(x>=0 && x<*w && y>=0 && y<*h) {
        sR[tx][ty]=*(InR+(y**w+x));
        sG[tx][ty]=*(InR+(y**w+x));
        sB[tx][ty]=*(InR+(y**w+x));
    } else {
        sR[tx][ty] = 0.;
        sG[tx][ty] = 0.;
        sB[tx][ty] = 0.;
    }
    __syncthreads();

    if(tx<BLOCK_SIZE && ty<BLOCK_SIZE) {
        float R=0.,G=0.,B=0.;
        for(int i=0; i<3; ++i)
            for(int j=0; j<3; ++j)
            {
                float const x = G_St[i][j];
                R += sR[tx+(i-1)][ty+(j-1)] * x;
                G += sG[tx+(i-1)][ty+(j-1)] * x;
                B += sB[tx+(i-1)][ty+(j-1)] * x;
            }
        __syncthreads();
        if(x<*w && y<*h) {
        *(OutR+(y**w+x)) = R;
        *(OutG+(y**w+x)) = G;
        *(OutB+(y**w+x)) = B;
      }
    }
#elif 0
    // Put Stencil in constant memory
        float R=0.,G=0.,B=0.;
        for(int i=0; i<3; ++i)
            for(int j=0; j<3; ++j)
            {
                int xpos = x+(i-1);
                int ypos = y+(j-1);
                R += *(InR+(ypos**w+xpos)) * G_St[i][j];
                G += *(InG+(ypos**w+xpos)) * G_St[i][j];
                B += *(InB+(ypos**w+xpos)) * G_St[i][j];
            }
        *(OutR+(y**w+x)) = R;
        *(OutG+(y**w+x)) = G;
        *(OutB+(y**w+x)) = B;
#elif 0
    // Process a bunch of data points in one thread
#endif
}

image_t * stencil_cuda(
    image_t const * const input,
    float stencil[3][3],
    int const num_times)
{
  float *Ired, *Iblue, *Igreen;
  float *Ored, *Oblue, *Ogreen;
  int *d_h, *d_w;
  float *d_stencil;

  const int h=input->height;
  const int w=input->width;
  cudaError_t rv;

  // Allocate memory on the device
  if(cudaSuccess != (rv = cudaMalloc((void **)&d_h, sizeof(int))))              printf("rv = %d = %s.\n", rv, cudaGetErrorString(rv));
  if(cudaSuccess != (rv = cudaMalloc((void **)&d_w, sizeof(int))))              printf("rv = %d = %s.\n", rv, cudaGetErrorString(rv));
  if(cudaSuccess != (rv = cudaMalloc((void **)&d_stencil, sizeof(float)*3*3)))  printf("rv = %d = %s.\n", rv, cudaGetErrorString(rv));

  if(cudaSuccess != (rv = cudaMalloc((void **)&Ired,   h*w*sizeof(float))))      printf("rv = %d = %s.\n", rv, cudaGetErrorString(rv));
  if(cudaSuccess != (rv = cudaMalloc((void **)&Igreen, h*w*sizeof(float))))      printf("rv = %d = %s.\n", rv, cudaGetErrorString(rv));
  if(cudaSuccess != (rv = cudaMalloc((void **)&Iblue,  h*w*sizeof(float))))      printf("rv = %d = %s.\n", rv, cudaGetErrorString(rv));

  if(cudaSuccess != (rv = cudaMalloc((void **)&Ored,   h*w*sizeof(float))))      printf("rv = %d = %s.\n", rv, cudaGetErrorString(rv));
  if(cudaSuccess != (rv = cudaMalloc((void **)&Ogreen, h*w*sizeof(float))))      printf("rv = %d = %s.\n", rv, cudaGetErrorString(rv));
  if(cudaSuccess != (rv = cudaMalloc((void **)&Oblue,  h*w*sizeof(float))))      printf("rv = %d = %s.\n", rv, cudaGetErrorString(rv));

  // Copy data to GPU memory
  // Just copy the width and height
  cudaMemcpy(d_h, &input->height, sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_w, &input->width,  sizeof(int), cudaMemcpyHostToDevice);

  cudaMemcpy(Ired,   input->red,   h*w*sizeof(*input->red),   cudaMemcpyHostToDevice);
  cudaMemcpy(Igreen, input->green, h*w*sizeof(*input->green), cudaMemcpyHostToDevice);
  cudaMemcpy(Iblue,  input->blue,  h*w*sizeof(*input->blue),  cudaMemcpyHostToDevice);

  if(cudaSuccess != (rv=cudaMemcpy(d_stencil, &stencil,  sizeof(float)*3*3,  cudaMemcpyHostToDevice)))     printf("rv = %d = %s.\n", rv, cudaGetErrorString(rv));

  if(cudaSuccess != (rv=cudaMemcpyToSymbol(G_St, stencil, 3*3*sizeof(float))))
      printf("Const mem error! sizeof(%lu): Err msg %s  \n", 3*3*sizeof(float), cudaGetErrorString(rv));

  cudaDeviceSynchronize();

  dim3 dimBlock(BLOCK_SIZE+2, BLOCK_SIZE+2,1);
  dim3 dimGrid(ceil(float (w)/BLOCK_SIZE), ceil(float (h)/BLOCK_SIZE),1);

  printf("Image (w=%d x h=%d).\n", w, h);
  printf("Grids %dx%d Blocks %dx%d. St[1][1] %f\n", dimGrid.x, dimGrid.y, dimBlock.x, dimBlock.y, stencil[1][1]);

  for(int i=0; i<num_times ;++i)
  {
      printf("Calling kernel .... \n");
      stencil_kernel<<<dimGrid, dimBlock>>>(d_h, d_w, d_stencil,
              Ired, Igreen, Iblue,
              Ored, Ogreen, Oblue);
      if ( cudaSuccess != (rv=cudaGetLastError()) )
          printf( "**** Error! %d: %s ****\n", rv, cudaGetErrorString(rv) );
  }

  image_t * output = image_alloc(input->width, input->height);

  // Get the result from GPU
  cudaMemcpy(output->red,   Ored,   h*w*sizeof(*input->red),   cudaMemcpyDeviceToHost);
  cudaMemcpy(output->green, Ogreen, h*w*sizeof(*input->green), cudaMemcpyDeviceToHost);
  cudaMemcpy(output->blue,  Oblue,  h*w*sizeof(*input->blue),  cudaMemcpyDeviceToHost);

  cudaFree(d_w); cudaFree(d_h); cudaFree(d_stencil);
  cudaFree(Ired); cudaFree(Iblue); cudaFree(Igreen);
  cudaFree(Ored); cudaFree(Oblue); cudaFree(Ogreen);
  return output;
}
