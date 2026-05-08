#include "physlib/percolation.h"
#include "physlib/random.h"
#include <stdlib.h>
#include <string.h>

int lattice_init(Lattice *lat, size_t l)
{
    lat->l        = l;
    lat->n_filled = 0;
    lat->sites  = calloc(l * l, sizeof(size_t));
    lat->pos_av = malloc(l * l * sizeof(size_t));
    if (!lat->sites || !lat->pos_av) {
        free(lat->sites);
        free(lat->pos_av);
        return -1;
    }
    for (size_t i = 0; i < l * l; i++)
        lat->pos_av[i] = i;
    return 0;
}

void lattice_free(Lattice *lat)
{
    free(lat->sites);
    free(lat->pos_av);
    lat->sites = lat->pos_av = NULL;
}

void lattice_clear(Lattice *lat)
{
    size_t l2 = lat->l * lat->l;
    memset(lat->sites, 0, l2 * sizeof(size_t));
    for (size_t i = 0; i < l2; i++)
        lat->pos_av[i] = i;
    lat->n_filled = 0;
}

/* Path-compressed find. */
size_t lattice_find(Lattice *lat, size_t idx)
{
    if (lat->sites[idx] == 0) return (size_t)-1;

    /* Find root. */
    size_t root = idx;
    while (lat->sites[root] != root + 1)
        root = lat->sites[root] - 1;

    /* Path compression. */
    while (lat->sites[idx] != root + 1) {
        size_t next = lat->sites[idx] - 1;
        lat->sites[idx] = root + 1;
        idx = next;
    }
    return root;
}

/* Union by larger index (preserves the original scheme from the course). */
void lattice_union(Lattice *lat, size_t a, size_t b)
{
    size_t ra = lattice_find(lat, a);
    size_t rb = lattice_find(lat, b);
    if (ra == (size_t)-1 || rb == (size_t)-1 || ra == rb) return;

    size_t big   = ra > rb ? ra : rb;
    size_t small = ra > rb ? rb : ra;
    lat->sites[small] = big + 1;
    lat->sites[big]   = big + 1;
}

void lattice_add_site(Lattice *lat, size_t idx)
{
    size_t l  = lat->l;
    size_t l2 = l * l;

    lat->sites[idx] = idx + 1; /* self-root */
    lat->n_filled++;

    if (idx >= l    && lat->sites[idx - l] != 0) lattice_union(lat, idx, idx - l);
    if (idx + l < l2 && lat->sites[idx + l] != 0) lattice_union(lat, idx, idx + l);
    if (idx % l != 0 && lat->sites[idx - 1] != 0) lattice_union(lat, idx, idx - 1);
    if ((idx + 1) % l != 0 && idx + 1 < l2 && lat->sites[idx + 1] != 0)
        lattice_union(lat, idx, idx + 1);
}

void lattice_fill_n(Lattice *lat, size_t n)
{
    size_t l2 = lat->l * lat->l;
    if (n > l2) n = l2;

    for (size_t i = lat->n_filled; i < n; i++) {
        long slot = phys_rand_range((long)(l2 - i));
        size_t pos = lat->pos_av[slot];
        lat->pos_av[slot] = lat->pos_av[l2 - 1 - i];
        lattice_add_site(lat, pos);
    }
}

void lattice_fill_rho(Lattice *lat, double rho)
{
    size_t n = (size_t)(rho * (double)(lat->l * lat->l));
    lattice_fill_n(lat, n);
}

bool lattice_percolates(Lattice *lat)
{
    size_t l = lat->l;
    for (size_t col = 0; col < l; col++) {
        if (lat->sites[col] == 0) continue;
        size_t top_root = lattice_find(lat, col);
        for (size_t bcol = 0; bcol < l; bcol++) {
            size_t bot = (l - 1) * l + bcol;
            if (lat->sites[bot] == 0) continue;
            if (lattice_find(lat, bot) == top_root) return true;
        }
    }
    return false;
}

size_t lattice_n_clusters(Lattice *lat)
{
    size_t l2 = lat->l * lat->l;
    size_t count = 0;
    for (size_t i = 0; i < l2; i++) {
        if (lat->sites[i] != 0 && lattice_find(lat, i) == i)
            count++;
    }
    return count;
}

double lattice_mean_cluster_size(Lattice *lat)
{
    size_t l2 = lat->l * lat->l;
    if (lat->n_filled == 0) return 0.0;

    /* Count cluster sizes using a temporary array indexed by root. */
    size_t *sz = calloc(l2, sizeof(size_t));
    if (!sz) return 0.0;

    for (size_t i = 0; i < l2; i++) {
        if (lat->sites[i] != 0)
            sz[lattice_find(lat, i)]++;
    }

    double sum_s2 = 0.0;
    for (size_t i = 0; i < l2; i++) {
        if (sz[i] > 0)
            sum_s2 += (double)sz[i] * (double)sz[i];
    }
    free(sz);
    return sum_s2 / (double)lat->n_filled;
}
