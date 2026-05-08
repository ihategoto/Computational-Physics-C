/*
 * Rössler attractor: dx/dt = -y-z,  dy/dt = x+ay,  dz/dt = b+(x-c)z
 *
 * Usage: rossler <a> <b> <c> <x0> <y0> <z0> <dt> <t_max> <method>
 *   method: 0 = RK2, 1 = RK4
 *
 * Output columns: t  x  y  z
 */
#include <stdio.h>
#include <stdlib.h>
#include "physlib/ode3d.h"

typedef struct { double a, b, c; } RosslerParams;

static void rossler_deriv(const State3D *s, State3D *ds, void *ctx)
{
    RosslerParams *p = ctx;
    ds->x = -s->y - s->z;
    ds->y =  s->x + p->a * s->y;
    ds->z =  p->b + (s->x - p->c) * s->z;
}

int main(int argc, char *argv[])
{
    if (argc != 10) {
        fprintf(stderr,
                "Usage: %s a b c x0 y0 z0 dt t_max method\n"
                "  method: 0=RK2  1=RK4\n", argv[0]);
        return EXIT_FAILURE;
    }

    RosslerParams p = { atof(argv[1]), atof(argv[2]), atof(argv[3]) };
    State3D s = { atof(argv[4]), atof(argv[5]), atof(argv[6]) };
    double dt    = atof(argv[7]);
    double t_max = atof(argv[8]);
    int    method = atoi(argv[9]);

    if (method != 0 && method != 1) {
        fprintf(stderr, "method must be 0 (RK2) or 1 (RK4)\n");
        return EXIT_FAILURE;
    }

    printf("# a=%g b=%g c=%g x0=%g y0=%g z0=%g dt=%g method=%s\n",
           p.a, p.b, p.c, s.x, s.y, s.z, dt, method == 0 ? "RK2" : "RK4");
    printf("# t  x  y  z\n");

    int steps = (int)(t_max / dt + 0.5);
    for (int i = 0; i <= steps; i++) {
        printf("%.10g %.10g %.10g %.10g\n", i * dt, s.x, s.y, s.z);
        if (method == 0) rk2_3d_step(&s, dt, rossler_deriv, &p);
        else             rk4_3d_step(&s, dt, rossler_deriv, &p);
    }

    return EXIT_SUCCESS;
}
