#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>
#include "reboundx.h"


int main(int argc, char* argv[]) {
        struct reb_simulation* r = reb_create_simulation();
        r->dt = 0.1;
        r->integrator = REB_INTEGRATOR_WHFAST;

        struct reb_particle p1 = {0};
        p1.m = 1.;
        reb_add(r, p1);

        struct reb_particle p2 = {0};
        p2.x = 1;
        p2.vy = 1;
        p2.m = 0.;
        reb_add(r, p2);
        
        struct rebx_extras* rebx = rebx_attach(r);  // first initialize rebx
        struct rebx_operator* mass_loss_test = rebx_load_operator(rebx, "mass_loss_test"); // add our new force
        rebx_add_operator(rebx, mass_loss_test);
        reb_move_to_com(r);
        reb_integrate(r,10000.);
        
        printf("%f\n",p1.m);
        printf("%f\n", r->t);
        printf("%f", p2.x);
        
}

