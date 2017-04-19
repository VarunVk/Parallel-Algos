#ifndef STENCIL_STENCIL_H
#define STENCIL_STENCIL_H


#include "image.h"


/**
* @brief Apply a 3x3 stencil to an image, returning the new image.
*
* @param input The input image.
* @param stencil The stencil.
* @param num_times How many times to apply the stencil.
*
* @return The new image.
*/
image_t * stencil_omp(
    image_t const * const input,
    float stencil[3][3],
    int const num_times);

/**
* @brief Apply a 3x3 stencil to an image using CUDA, returning the new image.
*
* @param input The input image.
* @param stencil The stencil.
* @param num_times How many times to apply the stencil.
*
* @return The new image.
*/
image_t * stencil_cuda(
    image_t const * const input,
    float stencil[3][3],
    int const num_times);


#endif
