
#include <assert.h>

#include "image.h"

static void p_apply_stencil(
    float const * const restrict input,
    float       * const restrict output,
    float stencil[3][3],
    int const width,
    int const height)
{
  #pragma omp parallel
  {
    /* just copy boundaries */
    #pragma omp for schedule(static) nowait
    for(int j=0; j < width; ++j) {
      output[j + (0)] = input[j + (0)];
      output[j + ((height-1) * width)] = input[j + ((height-1) * width)];
    }
    #pragma omp for schedule(static) nowait
    for(int i=0; i < height; ++i) {
      output[0 + (i*width)]  = input[0 + (i*width)];
      output[(width-1) + (i*width)]  = input[(width-1) + (i*width)];
    }

    #pragma omp for schedule(static)
    for(int i=1; i < height-1; ++i) {
      for(int j=1; j < width-1; ++j) {

        float v = 0.;
        for(int is=0; is < 3; ++is) {
          for(int js=0; js < 3; ++js) {
            /* global image indices */
            int const g_i = i + (is-1);
            int const g_j = j + (js-1);

            v += stencil[is][js] * input[g_j + (g_i * width)];
          }
        }

        output[j + (i*width)] = v;
      }
    }
  } /* end omp parallel */
}



image_t * stencil_omp(
    image_t const * const input,
    float stencil[3][3],
    int const num_times)
{
  image_t * output = image_alloc(input->width, input->height);

  for(int i=0;i < num_times; ++i) {
    /* Apply stencil to each channel separately. */
    p_apply_stencil(input->red,   output->red,   stencil, input->width, input->height);
    p_apply_stencil(input->green, output->green, stencil, input->width, input->height);
    p_apply_stencil(input->blue,  output->blue,  stencil, input->width, input->height);
  }

  return output;
}

