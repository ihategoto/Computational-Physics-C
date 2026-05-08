/*
 * physlib/ode.h — generic 1D ODE integrators for systems of the form
 *
 *   dx/dt = v
 *   dv/dt = a(x, v, ctx)
 *
 * All integrators advance the state by one time step dt in-place.
 */
#ifndef PHYSLIB_ODE_H
#define PHYSLIB_ODE_H

/* Phase-space state: position x and velocity v. */
typedef struct {
    double x;
    double v;
} PhysState;

/* Acceleration function pointer: returns a(x, v) given a user-supplied context. */
typedef double (*AccelFn)(double x, double v, void *ctx);

/* Euler method — O(dt). */
void euler_step(PhysState *s, double dt, AccelFn accel, void *ctx);

/* Euler-Cromer (symplectic Euler) — O(dt), conserves a shadow Hamiltonian. */
void euler_cromer_step(PhysState *s, double dt, AccelFn accel, void *ctx);

/*
 * Leapfrog (drift-kick form) — O(dt²), symplectic.
 *
 * The algorithm is NOT self-starting: before the integration loop call
 * leapfrog_init() once to advance v to the half-step v_{1/2}.
 * Inside the loop use leapfrog_step(). To recover v at integer steps
 * (e.g. for energy evaluation or output) call leapfrog_velocity().
 */
void leapfrog_init(PhysState *s, double dt, AccelFn accel, void *ctx);
void leapfrog_step(PhysState *s, double dt, AccelFn accel, void *ctx);
double leapfrog_velocity(const PhysState *s, double dt, AccelFn accel, void *ctx);

/* Midpoint method — O(dt²). */
void midpoint_step(PhysState *s, double dt, AccelFn accel, void *ctx);

/* Velocity Verlet — O(dt²), symplectic. */
void verlet_v_step(PhysState *s, double dt, AccelFn accel, void *ctx);

/*
 * Predictor-corrector (Euler predictor + trapezoidal corrector, iterated to
 * convergence with tolerance 1e-10) — effectively O(dt²).
 */
void predictor_corrector_step(PhysState *s, double dt, AccelFn accel, void *ctx);

/* Runge-Kutta 2nd order (midpoint method) — O(dt²). */
void rk2_step(PhysState *s, double dt, AccelFn accel, void *ctx);

/* Runge-Kutta 4th order — O(dt⁴). */
void rk4_step(PhysState *s, double dt, AccelFn accel, void *ctx);

#endif /* PHYSLIB_ODE_H */
