
#include <stdio.h>

extern "C"
{
#include "image.h"
#include "stencil.h"
}



image_t * stencil_cuda(
    image_t const * const input,
    float stencil[3][3],
    int const num_times)
{
  image_t * output = image_alloc(input->width, input->height);
  return output;
}


