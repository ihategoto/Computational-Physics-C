/*
 * physlib/ode3d.h — generic 3D ODE integrators for autonomous systems
 *
 *   d/dt (x, y, z) = f(x, y, z, ctx)
 *
 * Designed for flows like the Rössler attractor.
 */
#ifndef PHYSLIB_ODE3D_H
#define PHYSLIB_ODE3D_H

/* Three-dimensional state vector. */
typedef struct {
    double x, y, z;
} State3D;

/*
 * Derivative function pointer.
 * Given the current state s, writes the time derivatives into ds.
 * ctx carries any system parameters (e.g. Rössler a, b, c).
 */
typedef void (*Deriv3DFn)(const State3D *s, State3D *ds, void *ctx);

/* Runge-Kutta 2nd order (midpoint) — O(dt²). */
void rk2_3d_step(State3D *s, double dt, Deriv3DFn f, void *ctx);

/* Runge-Kutta 4th order — O(dt⁴). */
void rk4_3d_step(State3D *s, double dt, Deriv3DFn f, void *ctx);

#endif /* PHYSLIB_ODE3D_H */
