#include "physlib/ode3d.h"

void rk2_3d_step(State3D *s, double dt, Deriv3DFn f, void *ctx)
{
    State3D k1, mid;
    f(s, &k1, ctx);

    mid.x = s->x + k1.x * dt * 0.5;
    mid.y = s->y + k1.y * dt * 0.5;
    mid.z = s->z + k1.z * dt * 0.5;

    State3D k2;
    f(&mid, &k2, ctx);

    s->x += k2.x * dt;
    s->y += k2.y * dt;
    s->z += k2.z * dt;
}

void rk4_3d_step(State3D *s, double dt, Deriv3DFn f, void *ctx)
{
    State3D k1, k2, k3, k4, tmp;

    f(s, &k1, ctx);

    tmp.x = s->x + k1.x * dt * 0.5;
    tmp.y = s->y + k1.y * dt * 0.5;
    tmp.z = s->z + k1.z * dt * 0.5;
    f(&tmp, &k2, ctx);

    tmp.x = s->x + k2.x * dt * 0.5;
    tmp.y = s->y + k2.y * dt * 0.5;
    tmp.z = s->z + k2.z * dt * 0.5;
    f(&tmp, &k3, ctx);

    tmp.x = s->x + k3.x * dt;
    tmp.y = s->y + k3.y * dt;
    tmp.z = s->z + k3.z * dt;
    f(&tmp, &k4, ctx);

    s->x += dt * (k1.x + 2.0*k2.x + 2.0*k3.x + k4.x) / 6.0;
    s->y += dt * (k1.y + 2.0*k2.y + 2.0*k3.y + k4.y) / 6.0;
    s->z += dt * (k1.z + 2.0*k2.z + 2.0*k3.z + k4.z) / 6.0;
}
