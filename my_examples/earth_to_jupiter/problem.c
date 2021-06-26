/**
 * Earth to Jupiter Hohman transfer orbit
 *
 * David Harper <david@obliquity.net>
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "tools.h"

void heartbeat(struct reb_simulation* const r);

int main(int argc, char* argv[]){
    if (argc < 5) {
        fprintf(stderr, "Usage: %s tEnd L[Jupiter] L[Saturn] L[Uranus] L[Neptune]\n", argv[0]);
        return 1;
    }

    double tmax = atof(argv[1]);

    struct reb_simulation* r = reb_create_simulation();

    const double k = 0.01720209895; // Gaussian constant
    r->G = k * k;            // These are the same units as used by the mercury6 code.

    r->dt = 1.0;
    r->integrator = REB_INTEGRATOR_IAS15;
    r->heartbeat        = heartbeat;
    r->exact_finish_time = 1; // Finish exactly at tmax in reb_integrate(). Default is already 1.

    // Sun + inner planets
    struct reb_particle pSun = {0};
    pSun.m      = 1.00000597682;
    reb_add(r, pSun);

    // Adding a particle with orbital elements requires the following 7 things
    // There's more flexibility in python for passing different orbital elements
    // for edge cases.  The user has to calculate these manually in C and pass
    // the elements below

    struct reb_particle primary = r->particles[0];

    double masses[4] = {
      1. / 1047.355,
      1. / 3501.6,
      1. / 22869.,
      1. / 19314.
    };

    double sma[4] = {  5.2, 9.6, 19.2, 30.1 };

    double e = 0.0;
    double inc = 0.0;
    double Omega = 0.0;
    double omega = 0.0;

    for (int i = 0; i < 4; i++) {
      double f = atof(argv[i+2]) * M_PI/180.0;

      struct reb_particle pPlanet = reb_tools_orbit_to_particle(r->G, primary, masses[i], sma[i], e, inc, Omega, omega, f);

      reb_add(r, pPlanet);
    }

    r->N_active = r->N;

    // Spacecraft
    double a = (1.0 + 5.2)/2.0;
    e = 1.0 - 1.0/a;

    struct reb_particle pSpacecraft = reb_tools_orbit_to_particle(r->G, primary, 0.0, a, e, inc, Omega, omega, 0.0);
    reb_add(r, pSpacecraft);

    int kSpacecraft = r->N - 1;

    int var_i_dx = reb_add_var_1st_order(r, kSpacecraft);
    r->particles[var_i_dx].vx = 1.;

    int var_i_dy = reb_add_var_1st_order(r, kSpacecraft);
    r->particles[var_i_dy].vy = 1.;

    //reb_integrate(r, tmax);

    while (r->t < tmax) {
      reb_step(r);
      heartbeat(r);
    }
}

void heartbeat(struct reb_simulation* const r) {
  printf("%14.6f %14.6f", r->t, r->dt);

  //for (int i = 0; i < 5; i++)
  //  printf("  %10.6f %10.6f", r->particles[i].x, r->particles[i].y);

  double dx = r->particles[5].x - r->particles[1].x;
  double dy = r->particles[5].y - r->particles[1].y;

  double dvx = r->particles[5].vx - r->particles[1].vx;
  double dvy = r->particles[5].vy - r->particles[1].vy;

  double rr = sqrt(dx*dx+dy*dy);
  double vdot = (dx*dvx + dy*dvy)/rr;

  printf("  %15.8f %15.8f", rr, vdot);

  printf("  %15.4f %15.4f", r->particles[6].x, r->particles[6].y);
  printf("  %15.4f %15.4f", r->particles[7].x, r->particles[7].y);

  printf("\n");
}
