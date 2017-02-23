/* Gives us high-resolution timers. */
#include <time.h>

/**
 * @brief Return the number of seconds since an unspecified time (e.g., Unix
 *        epoch). This is accomplished with a high-resolution monotonic timer,
 *        suitable for performance timing.
 *
 * @return The number of seconds.
 */
static inline double monotonic_seconds()
{
  /* Linux systems */
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/**
 * @brief Output the seconds elapsed while sorting. This excludes input and
 *        output time. This should be wallclock time, not CPU time.
 *
 * @param seconds Seconds spent sorting.
 */
static void print_time(
    double const seconds)
{
  printf("Sort Time: %0.04fs\n", seconds);
}

/**
 * @brief Write an array of integers to a file.
 *
 * @param filename The name of the file to write to.
 * @param numbers The array of numbers.
 * @param nnumbers How many numbers to write.
 */
static void print_numbers(
    char const * const filename,
    int const * const numbers,
    size_t const nnumbers)
{
  size_t i;
  FILE * fout;

  /* open file */
  if((fout = fopen(filename, "w")) == NULL) {
    fprintf(stderr, "error opening '%s'\n", filename);
    abort();
  }

  /* write numbers to fout */
  for(i = 0; i < nnumbers; ++i) {
    fprintf(fout, "%d\n", numbers[i]);
  }

  fclose(fout);
}
