#include "physlib/ode.h"
#include <math.h>

void euler_step(PhysState *s, double dt, AccelFn accel, void *ctx)
{
    double a = accel(s->x, s->v, ctx);
    s->x += s->v * dt;
    s->v += a * dt;
}

void euler_cromer_step(PhysState *s, double dt, AccelFn accel, void *ctx)
{
    s->v += accel(s->x, s->v, ctx) * dt;
    s->x += s->v * dt;
}

/*
 * Leapfrog (drift-kick form):
 *   init:  v_{1/2} = v_0 + a(x_0) * dt/2
 *   step:  x_{n+1} = x_n + v_{n+1/2} * dt
 *          v_{n+3/2} = v_{n+1/2} + a(x_{n+1}) * dt
 */
void leapfrog_init(PhysState *s, double dt, AccelFn accel, void *ctx)
{
    s->v += accel(s->x, s->v, ctx) * dt * 0.5;
}

void leapfrog_step(PhysState *s, double dt, AccelFn accel, void *ctx)
{
    s->x += s->v * dt;
    s->v += accel(s->x, s->v, ctx) * dt;
}

/* Reconstruct v at the integer step from the current half-step state.
   After leapfrog_step, s->x already holds x_{n+1}, so this returns v_{n+1}. */
double leapfrog_velocity(const PhysState *s, double dt, AccelFn accel, void *ctx)
{
    return s->v - accel(s->x, s->v, ctx) * dt * 0.5;
}

void midpoint_step(PhysState *s, double dt, AccelFn accel, void *ctx)
{
    double v_old = s->v;
    s->v += accel(s->x, s->v, ctx) * dt;
    s->x += (v_old + s->v) * 0.5 * dt;
}

void verlet_v_step(PhysState *s, double dt, AccelFn accel, void *ctx)
{
    double a0 = accel(s->x, s->v, ctx);
    s->x += s->v * dt + 0.5 * a0 * dt * dt;
    double a1 = accel(s->x, s->v, ctx);
    s->v += 0.5 * (a0 + a1) * dt;
}

/*
 * Predictor-corrector: Euler predictor, trapezoidal corrector iterated to
 * convergence (tolerance 1e-10 on both x and v).
 */
void predictor_corrector_step(PhysState *s, double dt, AccelFn accel, void *ctx)
{
    double a0    = accel(s->x, s->v, ctx);
    double x_p   = s->x + s->v * dt;   /* Euler prediction */
    double v_p   = s->v + a0 * dt;

    double diff;
    do {
        double x_prev = x_p, v_prev = v_p;
        double a1 = accel(x_p, v_p, ctx);
        x_p = s->x + (s->v + v_p) * 0.5 * dt;
        v_p = s->v + (a0 + a1)    * 0.5 * dt;
        diff = fabs(v_p - v_prev) + fabs(x_p - x_prev);
    } while (diff > 1e-10);

    s->x = x_p;
    s->v = v_p;
}

void rk2_step(PhysState *s, double dt, AccelFn accel, void *ctx)
{
    /* Midpoint (Heun) RK2:
       k1 = f(y_n)
       k2 = f(y_n + dt/2 * k1)
       y_{n+1} = y_n + dt * k2                */
    double k1_v = accel(s->x, s->v, ctx);
    double mx   = s->x + s->v  * dt * 0.5;
    double mv   = s->v + k1_v  * dt * 0.5;
    double k2_v = accel(mx, mv, ctx);
    s->x += mv   * dt;
    s->v += k2_v * dt;
}

void rk4_step(PhysState *s, double dt, AccelFn accel, void *ctx)
{
    double k1x = s->v,                              k1v = accel(s->x,                s->v,                ctx);
    double k2x = s->v + k1v * dt * 0.5,            k2v = accel(s->x + k1x*dt*0.5,  s->v + k1v*dt*0.5,  ctx);
    double k3x = s->v + k2v * dt * 0.5,            k3v = accel(s->x + k2x*dt*0.5,  s->v + k2v*dt*0.5,  ctx);
    double k4x = s->v + k3v * dt,                  k4v = accel(s->x + k3x*dt,      s->v + k3v*dt,      ctx);

    s->x += dt * (k1x + 2.0*k2x + 2.0*k3x + k4x) / 6.0;
    s->v += dt * (k1v + 2.0*k2v + 2.0*k3v + k4v) / 6.0;
}
