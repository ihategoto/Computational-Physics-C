/*
 * Harmonic oscillator: ẍ = -(k/m) x
 *
 * Usage: harmonic_oscillator <t_max> <dt> <x0> <v0> <k> <m> <method>
 *
 * method: euler | euler_cromer | leapfrog | midpoint | verlet_v | pc | rk2 | rk4
 *
 * Output columns: t  x  v  (E - E0) / E0
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "physlib/ode.h"

typedef struct { double k, m; } SHOParams;

static double sho_accel(double x, double v, void *ctx)
{
    (void)v;
    SHOParams *p = ctx;
    return -(p->k / p->m) * x;
}

static double energy(PhysState s, double k, double m)
{
    return 0.5 * m * s.v * s.v + 0.5 * k * s.x * s.x;
}

int main(int argc, char *argv[])
{
    if (argc != 8) {
        fprintf(stderr, "Usage: %s t_max dt x0 v0 k m method\n"
                        "  method: euler | euler_cromer | leapfrog | midpoint"
                        " | verlet_v | pc | rk2 | rk4\n", argv[0]);
        return EXIT_FAILURE;
    }

    double t_max = atof(argv[1]);
    double dt    = atof(argv[2]);
    SHOParams p  = { .k = atof(argv[5]), .m = atof(argv[6]) };

    PhysState s = { .x = atof(argv[3]), .v = atof(argv[4]) };

    typedef void (*StepFn)(PhysState *, double, AccelFn, void *);
    StepFn step = NULL;
    int is_leapfrog = 0;

    const char *method = argv[7];
    if      (!strcmp(method, "euler"))        step = euler_step;
    else if (!strcmp(method, "euler_cromer")) step = euler_cromer_step;
    else if (!strcmp(method, "leapfrog"))   { step = leapfrog_step; is_leapfrog = 1; }
    else if (!strcmp(method, "midpoint"))     step = midpoint_step;
    else if (!strcmp(method, "verlet_v"))     step = verlet_v_step;
    else if (!strcmp(method, "pc"))           step = predictor_corrector_step;
    else if (!strcmp(method, "rk2"))          step = rk2_step;
    else if (!strcmp(method, "rk4"))          step = rk4_step;
    else {
        fprintf(stderr, "Unknown method: %s\n", method);
        return EXIT_FAILURE;
    }

    double e0 = energy(s, p.k, p.m);
    printf("# method=%s dt=%g k=%g m=%g x0=%g v0=%g\n",
           method, dt, p.k, p.m, s.x, s.v);
    printf("# t  x  v  (E-E0)/E0\n");

    if (is_leapfrog)
        leapfrog_init(&s, dt, sho_accel, &p);

    int steps = (int)(t_max / dt + 0.5);
    for (int i = 0; i <= steps; i++) {
        double v_out = is_leapfrog
                       ? leapfrog_velocity(&s, dt, sho_accel, &p)
                       : s.v;
        PhysState s_out = { s.x, v_out };
        double e_err = (energy(s_out, p.k, p.m) - e0) / e0;
        printf("%g %g %g %g\n", i * dt, s.x, v_out, e_err);
        step(&s, dt, sho_accel, &p);
    }

    return EXIT_SUCCESS;
}
