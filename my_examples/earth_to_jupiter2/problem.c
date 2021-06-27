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

void runSimulation(double tmax, double lJupiter0);

int main(int argc, char* argv[]){
    if (argc < 4) {
        fprintf(stderr, "Usage: %s tEnd LJ[start] LJ[end] LJ[step]\n", argv[0]);
        return 1;
    }

    double tmax = atof(argv[1]);
    double lJupiterStart = atof(argv[2]);
    double lJupiterEnd = atof(argv[3]);
    double lJupiterStep = atof(argv[4]);

    for (double lJupiter0 = lJupiterStart; lJupiter0 <= lJupiterEnd; lJupiter0 += lJupiterStep)
      runSimulation(tmax, lJupiter0);
}

void runSimulation(double tmax, double lJupiter0) {
    struct reb_simulation* r = reb_create_simulation();

    const double k = 0.01720209895; // Gaussian constant
    r->G = k * k;            // These are the same units as used by the mercury6 code.

    r->dt = 1.0;
    r->integrator = REB_INTEGRATOR_IAS15;
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

    double mJupiter = 1. / 1047.355;

    double sma = 5.2;

    double e = 0.0;
    double inc = 0.0;
    double Omega = 0.0;
    double omega = 0.0;

    double f = lJupiter0 * M_PI/180.0;

    struct reb_particle pJupiter = reb_tools_orbit_to_particle(r->G, primary, mJupiter, sma, e, inc, Omega, omega, f);

    reb_add(r, pJupiter);

    r->N_active = r->N;

    int kJupiter = r->N - 1;

    // Spacecraft
    double a = (1.0 + 5.2)/2.0;
    e = 1.0 - 1.0/a;

    struct reb_particle pSpacecraft = reb_tools_orbit_to_particle(r->G, primary, 0.0, a, e, inc, Omega, omega, 0.0);
    reb_add(r, pSpacecraft);

    int kSpacecraft = r->N - 1;

    double trr = -1.0, rr, vdot;

    while (r->t < tmax && vdot <= 0.0) {
      reb_step(r);

      double dx = r->particles[kSpacecraft].x - r->particles[kJupiter].x;
      double dy = r->particles[kSpacecraft].y - r->particles[kJupiter].y;

      double dvx = r->particles[kSpacecraft].vx - r->particles[kJupiter].vx;
      double dvy = r->particles[kSpacecraft].vy - r->particles[kJupiter].vy;

      trr = r->t;
      rr = sqrt(dx*dx+dy*dy);
      vdot = (dx*dvx + dy*dvy)/rr;
    }

    printf("%8.3f", lJupiter0);

    double aj = 71492.0/149597870.700;

    if (trr > 0.0) {
      double dx = r->particles[kSpacecraft].x - r->particles[kJupiter].x;
      double dy = r->particles[kSpacecraft].y - r->particles[kJupiter].y;

      double dvx = r->particles[kSpacecraft].vx - r->particles[kJupiter].vx;
      double dvy = r->particles[kSpacecraft].vy - r->particles[kJupiter].vy;

      double xj = r->particles[1].x - r->particles[0].x;
      double yj = r->particles[1].y - r->particles[0].y;

      double rj = sqrt(xj * xj + yj * yj);

      xj /= rj;
      yj /= rj;

      double yr = xj * dx + yj * dy;
      double xr = yj * dx - xj * dy;

      double yvr = xj * dvx + yj * dvy;
      double xvr = yj * dvx - xj * dvy;

      printf(" %14.6f %15.8f %15.3f %15.8f %15.8f %15.8f %15.8f\n", trr, rr, rr/aj,
        xr, yr, xvr, yvr);
    } else
      printf(" NO CLOSEST APPROACH\n");
}
