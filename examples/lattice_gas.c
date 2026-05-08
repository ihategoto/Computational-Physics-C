/*
 * 2D lattice gas: N hard-core particles on an L×L toroidal grid.
 * At each step every particle attempts a move to a random neighbour;
 * the move is rejected if the target cell is occupied (excluded volume).
 *
 * Usage: lattice_gas <L> <N> <steps>
 *
 * Output columns: step  <r²>/step   (mean squared displacement per step)
 */
#include <stdio.h>
#include <stdlib.h>
#include "physlib/random.h"

typedef struct {
    int x, y;   /* current position on the torus */
    int xt, yt; /* total displacement (unwrapped) */
} Particle;

static int L, N;

static void fill_random(int **grid, Particle *p)
{
    /* Fisher-Yates: pick N distinct cells. */
    int *pool = malloc((size_t)(L * L) * sizeof(int));
    if (!pool) { perror("malloc"); exit(EXIT_FAILURE); }
    for (int i = 0; i < L * L; i++) pool[i] = i;

    for (int i = 0; i < N; i++) {
        long j = phys_rand_range((long)(L * L - i));
        int pos = pool[j];
        pool[j] = pool[L * L - 1 - i];
        int px = pos % L, py = pos / L;
        grid[py][px] = i;
        p[i].x = px;  p[i].y = py;
        p[i].xt = 0;  p[i].yt = 0;
    }
    free(pool);
}

static void step(int **grid, Particle *p)
{
    /* Randomise update order. */
    int *order = malloc((size_t)N * sizeof(int));
    if (!order) { perror("malloc"); exit(EXIT_FAILURE); }
    for (int i = 0; i < N; i++) order[i] = i;
    for (int i = 0; i < N; i++) {
        long j = phys_rand_range((long)(N - i));
        int tmp = order[j]; order[j] = order[N - 1 - i]; order[N - 1 - i] = tmp;
    }

    for (int k = 0; k < N; k++) {
        int id = order[k];
        int cx = p[id].x, cy = p[id].y;
        int dir = (int)phys_rand_range(4);
        int nx = cx, ny = cy;
        int dtx = 0, dty = 0;
        if      (dir == 0) { ny = (cy + 1) % L; dty =  1; }
        else if (dir == 1) { ny = (cy - 1 + L) % L; dty = -1; }
        else if (dir == 2) { nx = (cx + 1) % L; dtx =  1; }
        else               { nx = (cx - 1 + L) % L; dtx = -1; }

        if (grid[ny][nx] != -1) continue; /* occupied */

        grid[cy][cx]  = -1;
        grid[ny][nx]  = id;
        p[id].x  = nx;  p[id].y  = ny;
        p[id].xt += dtx; p[id].yt += dty;
    }
    free(order);
}

static double mean_sq_displacement(const Particle *p)
{
    double sum = 0.0;
    for (int i = 0; i < N; i++)
        sum += (double)(p[i].xt * p[i].xt + p[i].yt * p[i].yt);
    return sum / N;
}

int main(int argc, char *argv[])
{
    if (argc != 4) {
        fprintf(stderr, "Usage: %s L N steps\n", argv[0]);
        return EXIT_FAILURE;
    }
    L = atoi(argv[1]);
    N = atoi(argv[2]);
    int steps = atoi(argv[3]);

    if (N >= L * L) {
        fprintf(stderr, "N must be < L*L\n");
        return EXIT_FAILURE;
    }

    phys_srand(0);

    int **grid = malloc((size_t)L * sizeof(int *));
    for (int i = 0; i < L; i++) {
        grid[i] = malloc((size_t)L * sizeof(int));
        for (int j = 0; j < L; j++) grid[i][j] = -1;
    }
    Particle *p = malloc((size_t)N * sizeof(Particle));
    if (!grid || !p) { perror("malloc"); return EXIT_FAILURE; }

    fill_random(grid, p);

    printf("# L=%d N=%d steps=%d rho=%.4f\n", L, N, steps, (double)N/(L*L));
    printf("# step  <r^2>/step\n");

    for (int s = 1; s <= steps; s++) {
        step(grid, p);
        printf("%d %g\n", s, mean_sq_displacement(p) / s);
    }

    for (int i = 0; i < L; i++) free(grid[i]);
    free(grid);
    free(p);
    return EXIT_SUCCESS;
}
