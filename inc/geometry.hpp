#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <cmath>

// Cross product for 3D vectors (all arrays of length 3)
inline void cross(const double a[3], const double b[3], double result[3]) {
    result[0] = a[1]*b[2] - a[2]*b[1];
    result[1] = a[2]*b[0] - a[0]*b[2];
    result[2] = a[0]*b[1] - a[1]*b[0];
}

// Normalize a 3D vector (in-place)
inline void normalize(double a[3]) {
    double len = std::sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
    a[0] /= len;
    a[1] /= len;
    a[2] /= len;
}

// Calculate the dip offset from time t.
double zOffDipCalc(double t);

// Reflect the momentum part (elements 3-5 of state) given a normal and tangent (both 3D arrays)
void reflect(double state[6], const double norm[3], const double tang[3]);

// Check whether a dagger hit occurs.
bool checkDagHit(double x, double y, double z, double zOff);

// Check cleaning conditions (returns 1 or 2 if certain conditions are met, 0 otherwise).
int checkClean(const double state[6], const double prevState[6], double cleanHeight);

// Calculate the “dagger zeta” value.
double calcDagZeta(double x, double y, double z, double zOff);

// Check whether a low house hit occurs.
bool checkHouseHitLow(double x, double y, double z, double zOff);

// Check whether a high house hit occurs.
bool checkHouseHitHigh(double x, double y, double z, double zOff);

// Initialize a new Lyapunov perturbation state based on ref state (both are arrays of length 6).
// The result is written to pair.
void initializeLyapState(const double ref[6], double pair[6]);

// Reset the perturbed state pair so that the separation from ref is exactly EPSILON.
// Both ref and pair are assumed to be arrays of length 6.
void resetStates(const double ref[6], double pair[6]);

// Compute the “distance” between two 6-element states.
double distance(const double ref[6], const double pair[6]);

#endif /* GEOMETRY_HPP */
