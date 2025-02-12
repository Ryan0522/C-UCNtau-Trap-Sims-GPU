#ifndef TRACKUCN_HPP
#define TRACKUCN_HPP

#include "../inc/particle_types.h"
#include "../inc/fields_nate.h"
#include <limits>

#ifdef __CUDACC__
#define HD __host__ __device__
#else
#define HD
#endif

typedef struct fixedResult {
    double energy;
    double theta;
    double t;
    double settlingT;
    double ePerp;
    double x;
    double y;
    double z;
    double zOff;
    int nHit;
    int nHitHouseLow;
    int nHitHouseHigh;
    double eStart;
    double deathTime;
} fixedResult;

HD fixedResult fixedEffDaggerHitTime(double state[6], double dt, const trace &tr);

#endif // TRACKUCN_HPP