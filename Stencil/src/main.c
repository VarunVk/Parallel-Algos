


/*
 * For high-resolution timers.
 */
#ifndef _POSIX_C_SOURCE
  #define _POSIX_C_SOURCE 199309L
#endif

#include <stdio.h>
#include <stdlib.h>

#include "image.h"
#include "stencil.h"
#include "timer.h"



/*
 * Stencils.
 */
float identity[3][3] = {
  {0.,  0.,  0.},
  {0.,  1.,  0.},
  {0.,  0.,  0.}
};

float blur[3][3] = {
  {0.111,  0.111,  0.111},
  {0.111,  0.111,  0.111},
  {0.111,  0.111,  0.111}
};


float sharp[3][3] = {
  { 0.0,  -1.0,   0.0},
  {-1.0,   5.0,  -1.0},
  { 0.0,  -1.0,   0.0},
};


float smooth[3][3] = {
  {1.0,  2.0,  1.0},
  {2.0,  4.0,  2.0},
  {1.0,  2.0,  1.0},
};

float emboss[3][3] = {
  {2.,  0.,  0.},
  {0., -1.,  0.},
  {0.,  0., -1.}
};


float dither[3][3] = {
  {6.,  8.,  4.},
  {1.,  0.,  3.},
  {5.,  2.,  7.}
};



int main(
    int argc,
    char * * argv)
{
  if(argc < 3) {
    fprintf(stderr, "usage: %s <input-image.{jpg,bmp,png}> <# repeats>  [output-image.bmp]\n", *argv);
    return EXIT_FAILURE;
  }

  image_t * im = image_load(argv[1]);
  int const num_times = atoi(argv[2]);

  printf("#REPEATS=%d\n", num_times);

  size_t const flops = image_stencil_flops(im, 3) * num_times;

  /* Apply the stencil. */
  double start_time = monotonic_seconds();
  image_t * output_omp = stencil_omp(im, emboss, num_times);
  double elapsed_time = monotonic_seconds() - start_time;
  double gflops = 1e-9 * (double) flops / elapsed_time;
  printf("CPU time: %0.3fs  GFLOPS: %0.2f\n", elapsed_time, gflops);


  start_time = monotonic_seconds();
  image_t * output_cuda = stencil_cuda(im, emboss, num_times);
  elapsed_time = monotonic_seconds() - start_time;
  gflops = 1e-9 * (double) flops / elapsed_time;
  printf("GPU time: %0.3fs  GFLOPS: %0.2f\n", elapsed_time, gflops);

  /* Write the output if requested. */
  if(argc == 4) {
    image_write_bmp(argv[3], output_cuda);
  }

  image_free(im);
  image_free(output_omp);
  image_free(output_cuda);

  return EXIT_SUCCESS;
}



