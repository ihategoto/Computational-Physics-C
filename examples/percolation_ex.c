/*
 * Site percolation on an L×L square lattice.
 *
 * Usage: percolation_ex <L> <n_configs> <rho_min> <rho_max> <d_rho>
 *
 * Output columns: L  n_sites  percolates  n_clusters  mean_cluster_size
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "physlib/percolation.h"
#include "physlib/random.h"

int main(int argc, char *argv[])
{
    if (argc != 6) {
        fprintf(stderr,
                "Usage: %s L n_configs rho_min rho_max d_rho\n", argv[0]);
        return EXIT_FAILURE;
    }

    size_t L        = (size_t)atol(argv[1]);
    size_t n_cfg    = (size_t)atol(argv[2]);
    double rho_min  = atof(argv[3]);
    double rho_max  = atof(argv[4]);
    double d_rho    = atof(argv[5]);

    phys_srand(0);

    Lattice lat;
    if (lattice_init(&lat, L) != 0) {
        fprintf(stderr, "Failed to allocate lattice\n");
        return EXIT_FAILURE;
    }

    printf("# L=%zu n_configs=%zu rho_min=%g rho_max=%g d_rho=%g\n",
           L, n_cfg, rho_min, rho_max, d_rho);
    printf("# L  n_sites  percolates  n_clusters  mean_cluster_size\n");

    for (size_t cfg = 0; cfg < n_cfg; cfg++) {
        for (double rho = rho_min; rho < rho_max + 1e-9; rho += d_rho) {
            lattice_clear(&lat);
            lattice_fill_rho(&lat, rho);
            bool perc  = lattice_percolates(&lat);
            size_t nc  = lattice_n_clusters(&lat);
            double mcs = lattice_mean_cluster_size(&lat);
            printf("%zu %zu %d %zu %g\n",
                   L, lat.n_filled, (int)perc, nc, mcs);
        }
    }

    lattice_free(&lat);
    return EXIT_SUCCESS;
}
