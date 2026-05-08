/*
 * 1D random walk — computes <x>, <x²>, <x⁴> over n_traj trajectories.
 *
 * Usage: random_walk <t_max> <n_traj> [d] [mu]
 *   d   diffusion coefficient  (default 1.0 — simple walk)
 *   mu  drift coefficient      (default 0.0 — unbiased)
 *
 * Output columns: t  <x>  <x²>  <x⁴>
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "physlib/walker.h"
#include "physlib/random.h"

int main(int argc, char *argv[])
{
    if (argc < 3 || argc > 5) {
        fprintf(stderr, "Usage: %s t_max n_traj [d] [mu]\n", argv[0]);
        return EXIT_FAILURE;
    }

    unsigned long long t_max  = (unsigned long long)atoll(argv[1]);
    unsigned long long n_traj = (unsigned long long)atoll(argv[2]);
    double d  = (argc >= 4) ? atof(argv[3]) : 1.0;
    double mu = (argc >= 5) ? atof(argv[4]) : 0.0;

    if (walker_check_params(d, mu) != 0) {
        fprintf(stderr, "Invalid parameters: need d+mu²∈[0,1] and |mu|≤1\n");
        return EXIT_FAILURE;
    }

    phys_srand(0);

    double *x_av  = calloc(t_max + 1, sizeof(double));
    double *x2_av = calloc(t_max + 1, sizeof(double));
    double *x4_av = calloc(t_max + 1, sizeof(double));
    if (!x_av || !x2_av || !x4_av) {
        perror("calloc");
        return EXIT_FAILURE;
    }

    printf("# t_max=%llu n_traj=%llu d=%g mu=%g\n", t_max, n_traj, d, mu);
    printf("# t  <x>  <x^2>  <x^4>\n");

    for (unsigned long long traj = 0; traj < n_traj; traj++) {
        long long pos = 0;
        x_av[0]  += 0.0;
        x2_av[0] += 0.0;
        x4_av[0] += 0.0;
        for (unsigned long long t = 1; t <= t_max; t++) {
            pos += walker_1d_step(d, mu);
            double p2 = (double)(pos * pos);
            x_av[t]  += (double)pos;
            x2_av[t] += p2;
            x4_av[t] += p2 * p2;
        }
    }

    double inv = 1.0 / (double)n_traj;
    for (unsigned long long t = 0; t <= t_max; t++)
        printf("%llu %g %g %g\n", t, x_av[t]*inv, x2_av[t]*inv, x4_av[t]*inv);

    free(x_av);
    free(x2_av);
    free(x4_av);
    return EXIT_SUCCESS;
}
