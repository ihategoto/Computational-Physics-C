/*
 * physlib/percolation.h — site percolation on a 2D square lattice.
 *
 * Uses a Union-Find structure with path compression.  Sites are identified
 * by their flat index i = row * l + col (row 0 = top, row l-1 = bottom).
 *
 * Typical usage for a fixed density scan:
 *
 *   Lattice lat;
 *   lattice_init(&lat, L);
 *   for each density rho {
 *       lattice_clear(&lat);
 *       lattice_fill_rho(&lat, rho);
 *       bool perc = lattice_percolates(&lat);
 *   }
 *   lattice_free(&lat);
 *
 * For incremental filling (adding one site at a time) use lattice_add_site()
 * or lattice_fill_n().
 */
#ifndef PHYSLIB_PERCOLATION_H
#define PHYSLIB_PERCOLATION_H

#include <stddef.h>
#include <stdbool.h>

typedef struct {
    size_t  l;        /* side length */
    size_t *sites;    /* sites[i]: 0 = empty, root+1 = occupied */
    size_t *pos_av;   /* scratch pool for Fisher-Yates sampling */
    size_t  n_filled; /* number of occupied sites */
} Lattice;

/* Allocate a lattice of side l. Returns 0 on success, -1 on failure. */
int    lattice_init(Lattice *lat, size_t l);

/* Free all resources. */
void   lattice_free(Lattice *lat);

/* Reset to empty (resets pos_av pool too). */
void   lattice_clear(Lattice *lat);

/* Find the root of site idx (path-compressed). Returns (size_t)-1 if empty. */
size_t lattice_find(Lattice *lat, size_t idx);

/* Merge the clusters containing sites a and b. */
void   lattice_union(Lattice *lat, size_t a, size_t b);

/* Occupy site idx and merge with any adjacent occupied sites. */
void   lattice_add_site(Lattice *lat, size_t idx);

/*
 * Randomly fill exactly n sites using sampling without replacement
 * (incremental Fisher-Yates via pos_av).  Call lattice_clear() first.
 */
void   lattice_fill_n(Lattice *lat, size_t n);

/* Fill to density rho = n_sites / l². */
void   lattice_fill_rho(Lattice *lat, double rho);

/* Return true if the top row is connected to the bottom row. */
bool   lattice_percolates(Lattice *lat);

/* Number of distinct clusters (occupied sites only). */
size_t lattice_n_clusters(Lattice *lat);

/* Mean cluster size S = (Σ s²) / (Σ s), excluding the percolating cluster.
   Returns 0 if no occupied sites. */
double lattice_mean_cluster_size(Lattice *lat);

#endif /* PHYSLIB_PERCOLATION_H */
