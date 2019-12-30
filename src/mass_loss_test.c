#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

void rebx_mass_loss(struct reb_simulation* const sim, struct rebx_operator* const operator, const double dt){