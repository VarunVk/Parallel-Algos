#ifndef STENCIL_IMAGE_H
#define STENCIL_IMAGE_H

#include <stdlib.h>


/**
* @brief A structure representing an RGB image. Pixels are stored row-major.
*/
typedef struct
{
  int width;  /** The width of the image in pixels. */
  int height; /** The height of the image in pixels. */

  /*
   * Actual image data. These arrays are all of length (width x height). All
   * arrays are row-major!
   */
  float * red;   /** Red-channel values. */
  float * green; /** Green-channel values. */
  float * blue;  /** Blue-channel values. */
} image_t;



/**
* @brief Count the number of FLOPs performed during an arbitrary stencil
*        application.
*
* @param image       The image we are evaluating.
* @param stencil_dim The number of rows or columns in the stencil.
*
* @return The number of FLOPs.
*/
size_t image_stencil_flops(
    image_t const * const image,
    int const stencil_dim);


/**
* @brief Load an image from a file. The image can be JPEG, PNG, or BMP.
*
* @param filename The name of the file to load.
*
* @return  The image loaded into an `image_t`. Returns NULL on error.
*/
image_t * image_load(
    char const * const filename);


/**
* @brief Write an image to a bitmap file. Note that this will map the pixels
*        back to a valid range if necessary, so the image data *may* change.
*
* @param filename The name of the file to create.
* @param image The image structure to write.
*/
void image_write_bmp(
    char const * const filename,
    image_t * const image);


/**
* @brief Allocate an image.
*
* @param width The width, in pixels.
* @param height The height, in pixels.
*/
image_t * image_alloc(
    int width,
    int height);


/**
* @brief Free the memory allocated by `image_load()`.
*
* @param im The image to free.
*/
void image_free(
    image_t * im);


#endif
