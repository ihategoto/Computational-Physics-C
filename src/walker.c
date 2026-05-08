#include "physlib/walker.h"
#include "physlib/random.h"
#include <math.h>

int walker_1d_step(double d, double mu)
{
    double r      = phys_rand();
    double p_minus = 0.5 * (d + mu*mu - mu);
    double p_plus  = 0.5 * (d + mu*mu + mu);
    if (r < p_minus)              return -1;
    if (r < p_minus + p_plus)     return +1;
    return 0;
}

int walker_check_params(double d, double mu)
{
    double total = d + mu * mu;
    if (total < 0.0 || total > 1.0 || fabs(mu) > 1.0) return -1;
    return 0;
}
