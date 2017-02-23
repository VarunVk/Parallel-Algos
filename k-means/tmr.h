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
  * @brief Output the seconds elapsed while clustering.
  *
  * @param seconds Seconds spent on k-means clustering, excluding IO.
  */
static void print_time(double const seconds)
{
	  printf("k-means clustering time: %0.04fs\n", seconds);
}
