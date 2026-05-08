/*
 * physlib/walker.h — 1D lattice random walk with diffusion and drift.
 *
 * At each step the walker moves +1 with probability p+ = (d + μ² + μ)/2,
 * moves -1 with probability p- = (d + μ² - μ)/2, and stays with the
 * remaining probability 1 - d - μ².  The constraint d + μ² ∈ [0,1] and
 * |μ| ≤ 1 must hold; call walker_check_params() to verify.
 *
 * For a simple unbiased walk use d = 1, mu = 0.
 */
#ifndef PHYSLIB_WALKER_H
#define PHYSLIB_WALKER_H

/*
 * Take one step. Returns +1, -1, or 0.
 * Behaviour is undefined if parameters violate the constraint.
 */
int walker_1d_step(double d, double mu);

/*
 * Validate parameters. Returns 0 if valid, -1 otherwise.
 */
int walker_check_params(double d, double mu);

#endif /* PHYSLIB_WALKER_H */
