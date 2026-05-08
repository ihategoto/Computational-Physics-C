# physlib

A mini C library of numerical methods written for a **Computational Physics** course during my Physics undergraduate degree. The code in this repository is a refactored and consolidated version of the exercises and exam problems I wrote throughout the course.

## Contents

| Module | Header | What it provides |
|--------|--------|-----------------|
| ODE integrators (1D) | `physlib/ode.h` | Euler, Euler-Cromer, Leapfrog, Midpoint, Velocity Verlet, Predictor-Corrector, RK2, RK4 |
| ODE integrators (3D) | `physlib/ode3d.h` | RK2 and RK4 for autonomous 3D flows (e.g. Rössler attractor) |
| Random numbers | `physlib/random.h` | Seeding, uniform double, uniform integer range |
| Random walk | `physlib/walker.h` | 1D lattice walk with diffusion coefficient and drift |
| Percolation | `physlib/percolation.h` | Site percolation on a square lattice (Union-Find with path compression) |

Example programs using the library are in `examples/`.

## Requirements

- C99-compatible compiler (tested with GCC and Clang)
- GNU Make
- Standard C math library (`-lm`)

## Building

```bash
# Build the static library (lib/libphyslib.a)
make

# Build all example programs
make examples

# Remove compiled artifacts
make clean
```

## Usage

Link your program against `libphyslib.a` and include the headers:

```bash
gcc -Ipath/to/physics_C/include my_program.c \
    -Lpath/to/physics_C/lib -lphyslib -lm -o my_program
```

### Quick example — harmonic oscillator with RK4

```c
#include <stdio.h>
#include "physlib/ode.h"

typedef struct { double k, m; } Params;

static double accel(double x, double v, void *ctx) {
    (void)v;
    Params *p = ctx;
    return -(p->k / p->m) * x;
}

int main(void) {
    Params p = { .k = 1.0, .m = 1.0 };
    PhysState s = { .x = 1.0, .v = 0.0 };
    double dt = 0.01;

    for (int i = 0; i < 1000; i++) {
        printf("%g %g %g\n", i * dt, s.x, s.v);
        rk4_step(&s, dt, accel, &p);
    }
}
```

### Leapfrog note

Leapfrog is not self-starting — call `leapfrog_init()` once before the loop, then use `leapfrog_velocity()` inside the loop to reconstruct v at integer time steps:

```c
leapfrog_init(&s, dt, accel, &p);
for (int i = 0; i <= steps; i++) {
    double v = leapfrog_velocity(&s, dt, accel, &p);
    printf("%g %g %g\n", i * dt, s.x, v);
    leapfrog_step(&s, dt, accel, &p);
}
```

### 3D integrators

The `ode3d` module takes a derivative function `Deriv3DFn` that fills a `State3D` struct:

```c
#include "physlib/ode3d.h"

typedef struct { double a, b, c; } RosslerParams;

static void rossler(const State3D *s, State3D *ds, void *ctx) {
    RosslerParams *p = ctx;
    ds->x = -s->y - s->z;
    ds->y =  s->x + p->a * s->y;
    ds->z =  p->b + (s->x - p->c) * s->z;
}
```

### Percolation

```c
#include "physlib/percolation.h"
#include "physlib/random.h"

phys_srand(0);
Lattice lat;
lattice_init(&lat, 64);

for (double rho = 0.0; rho <= 1.0; rho += 0.05) {
    lattice_clear(&lat);
    lattice_fill_rho(&lat, rho);
    printf("rho=%.2f percolates=%d\n", rho, lattice_percolates(&lat));
}
lattice_free(&lat);
```

## Example programs

| Binary | Source | Description |
|--------|--------|-------------|
| `harmonic_oscillator` | `examples/harmonic_oscillator.c` | Compare all integrators on the simple harmonic oscillator; prints `t x v (E-E0)/E0` |
| `rossler` | `examples/rossler.c` | Rössler attractor trajectory with RK2 or RK4 |
| `random_walk` | `examples/random_walk.c` | 1D random walk moments `<x>`, `<x²>`, `<x⁴>` over many trajectories |
| `percolation_ex` | `examples/percolation_ex.c` | Site percolation probability scan |
| `lattice_gas` | `examples/lattice_gas.c` | 2D hard-core lattice gas, mean squared displacement |

Run any example without arguments to see its usage string.

## Project structure

```
physics_C/
├── include/physlib/   # Public headers
├── src/               # Library source files
├── examples/          # Standalone example programs
├── lib/               # Static library output (after make)
└── Makefile
```

## Background

The algorithms implemented here were studied in the course in roughly this order:

1. **ODE integration** — Euler and symplectic variants (Euler-Cromer, leapfrog, Verlet) for the harmonic oscillator and nonlinear pendulum; predictor-corrector and Runge-Kutta methods; application to chaotic systems (Rössler attractor).
2. **Random numbers and stochastic processes** — pseudo-random number generation, 1D random walks with drift and diffusion, 2D lattice gas with excluded volume.
3. **Percolation** — site percolation on a square lattice, cluster labelling, critical threshold, mean cluster size.

## License

Personal academic project — no formal license. Feel free to use the code for learning purposes.
