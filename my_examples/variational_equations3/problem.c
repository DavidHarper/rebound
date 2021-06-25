/**
 * Variational Equations
 *
 * This example shows how to use first and second
 * order variational equations.
 * See also https://github.com/hannorein/rebound/blob/master/ipython_examples/VariationalEquations.ipynb and Rein and Tamayo (2016).
 */
#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>

// This function creates a simulation with one star, one planet and one test particle.
struct reb_simulation* create_sim(int nTestParticles){
    struct reb_simulation* r = reb_create_simulation();
    // r->integrator = REB_INTEGRATOR_WHFAST;  Only first order variational equations supported in WHFast.
    struct reb_particle star = {0.};
    star.m = 1;
    reb_add(r, star);
    struct reb_particle planet = reb_tools_orbit_to_particle(1.,star,1e-3,1.,0.,0.,0.,0.,0.);
    reb_add(r, planet);

    r->N_active = r->N;

    for (int n = 0; n < nTestParticles; n++) {
        struct reb_particle testparticle = reb_tools_orbit_to_particle(1.,star,0.,1.7,0.1,0.2,0.3,0.4,0.5);
        reb_add(r, testparticle);
    }

    return r;
}

int main(int argc, char* argv[]) {
    struct reb_simulation* r;

    double DeltaX = argc > 1 ? atof(argv[1]) : 0.001;
    printf("\nShifting testparticle's x coordinate by %f.\n", DeltaX);

    int doDV = argc > 2;

    double DeltaVX = doDV ? atof(argv[2]) : 0.0;

    if (doDV)
        printf("\nShifting testparticle's x velocity by %f.\n", DeltaVX);

    r = create_sim(doDV ? 4 : 2);

    printf("The simulation has %d particles before adding the variational equations.\n", r->N);

    r->particles[3].x += DeltaX;

    if (doDV)
        r->particles[5].vx += DeltaVX;

    int var_i_dx = reb_add_var_1st_order(r, 2); // The 2 corresponds to the index of the testparticle that we vary.
    int var_ii_dx = reb_add_var_2nd_order(r, 2, var_i_dx, var_i_dx);

    r->particles[var_i_dx].x = 1.;

    int var_i_dv = -1, var_ii_dv = -1;

    if (doDV) {
        var_i_dv = reb_add_var_1st_order(r, 4); // The 2 corresponds to the index of the testparticle that we vary.

        var_ii_dv = reb_add_var_2nd_order(r, 4, var_i_dv, var_i_dv);

        r->particles[var_i_dv].vx = 1.;
    }

    printf("The simulation has %d particles.\n", r->N);
    printf("There are %d variational equations.\n", r->N_var);
    printf("There are %d variational particle configurations.\n", r->var_config_N);
    printf("\n");

    reb_integrate(r,100.);

    printf("Position of testparticle at t=100 in shifted simulation:       %.8f %.8f\n",
        r->particles[3].x,
        r->particles[3].y);

    printf("Position of testparticle at t=100 using 1st order var. eqs.:   %.8f %.8f\n",
        r->particles[2].x+DeltaX*r->particles[var_i_dx].x,
        r->particles[2].y+DeltaX*r->particles[var_i_dx].y);

    printf("Position of testparticle at t=100 using 2nd order var. eqs.:   %.8f %.8f\n",
        r->particles[2].x+DeltaX*r->particles[var_i_dx].x+DeltaX*DeltaX/2.*r->particles[var_ii_dx].x,
        r->particles[2].y+DeltaX*r->particles[var_i_dx].y+DeltaX*DeltaX/2.*r->particles[var_ii_dx].y);

    if (doDV) {
        printf("Position of testparticle at t=100 in shifted simulation:       %.8f %.8f\n",
            r->particles[5].x,
            r->particles[5].y);
        printf("Position of testparticle at t=100 using 1st order var. eqs.:   %.8f %.8f\n",
            r->particles[4].x+DeltaVX*r->particles[var_i_dv].x,
            r->particles[4].y+DeltaVX*r->particles[var_i_dv].y);

        printf("Position of testparticle at t=100 using 2nd order var. eqs.:   %.8f %.8f\n",
            r->particles[4].x+DeltaVX*r->particles[var_i_dv].x+DeltaVX*DeltaVX/2.*r->particles[var_ii_dv].x,
            r->particles[4].y+DeltaVX*r->particles[var_i_dv].y+DeltaVX*DeltaVX/2.*r->particles[var_ii_dv].y);
    }

    reb_free_simulation(r);
}
