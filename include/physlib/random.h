/*
 * physlib/random.h — thin wrappers around the C standard PRNG.
 *
 * All physlib modules use these functions so that the seed is set once
 * and shared consistently across the library.
 */
#ifndef PHYSLIB_RANDOM_H
#define PHYSLIB_RANDOM_H

/* Seed the PRNG. Pass 0 to seed from time(NULL). */
void phys_srand(unsigned int seed);

/* Uniform double in [0, 1). */
double phys_rand(void);

/* Uniform integer in [0, n-1]. Returns 0 for n <= 0. */
long phys_rand_range(long n);

#endif /* PHYSLIB_RANDOM_H */
