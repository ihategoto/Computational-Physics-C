#include "physlib/random.h"
#include <stdlib.h>
#include <time.h>

void phys_srand(unsigned int seed)
{
    srand(seed ? seed : (unsigned int)time(NULL));
}

/* [0, 1) — RAND_MAX+1 avoids the edge case rand()==RAND_MAX mapping to 1. */
double phys_rand(void)
{
    return (double)rand() / ((double)RAND_MAX + 1.0);
}

long phys_rand_range(long n)
{
    if (n <= 0) return 0;
    return (long)(phys_rand() * (double)n);
}
