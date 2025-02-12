#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdlib>
#include <assert.h>

#include "../setup.h"
#include "../inc/constants.h"
#include "../inc/geometry.hpp"

extern "C" {
    #include "../inc/xorshift.h"
}

// --- zOffDipCalc ---
// This function computes a dip offset given the current time t.
// (The details are taken from your original code.)
double zOffDipCalc(double t) {
    double acc = 6.0;
    double vel = 1.4; // was 1.6
    double spr = 51200 / 2.0;
    
    double holdT = HOLDTIME;  // assumed defined in constants.h
    // These arrays are assumed to be defined (for example, via macros or constants)
    double dipHeights[NDIPS] = HEIGHTS;   // for example, { ... }
    double dipEnds[NDIPS] = ENDTIMES;       // for example, { ... }
    
    if (t > dipEnds[NDIPS - 1]) {
        return 0.01;
    }
    
    int i;
    for (i = 0; i < NDIPS; i++) {
        if (dipEnds[i] > t) break;
    }
    
    if(i == 0) {
        return dipHeights[0];
    }
    
    double target = dipHeights[i];
    double start = dipHeights[i - 1];
    double distance = std::fabs(target - start) * 1000000;
    double time = (distance - vel * (vel / acc) * spr) / (vel * spr);
    
    double sign = (dipHeights[i] - dipHeights[i - 1] > 0) ? 1 : -1;
    if (dipHeights[i] == dipHeights[i - 1]) sign = 0;
    
    if (time > 0) {
        if ((t - dipEnds[i - 1]) < (vel / acc)) {
            return dipHeights[i - 1] + sign * (acc * spr * (t - dipEnds[i - 1]) * (t - dipEnds[i - 1]) / 2.0) / 1000000;
        }
        else if ((t - dipEnds[i - 1]) < ((vel / acc) + time)) {
            return dipHeights[i - 1] + sign * ((vel / acc) * vel * spr / 2.0 + (t - dipEnds[i - 1] - vel / acc) * vel * spr) / 1000000;
        }
        else if ((t - dipEnds[i - 1]) < (2 * (vel / acc) + time)) {
            return dipHeights[i - 1] + sign * ((vel / acc) * vel * spr / 2.0 + time * vel * spr +
                vel * spr * (t - dipEnds[i - 1] - time - vel / acc) - acc * spr * (t - dipEnds[i - 1] - time - vel / acc) * (t - dipEnds[i - 1] - time - vel / acc) / 2.0) / 1000000;
        }
        else {
            return dipHeights[i];
        }
    }
    else {
        double halfTime = std::sqrt(distance / (acc * spr));
        if (t - dipEnds[i - 1] < halfTime) {
            return dipHeights[i - 1] + sign * (acc * spr * (t - dipEnds[i - 1]) * (t - dipEnds[i - 1]) / 2.0) / 1000000;
        }
        else if (t - dipEnds[i - 1] < 2 * halfTime) {
            return dipHeights[i - 1] + (dipHeights[i] - dipHeights[i - 1]) / 2 +
                sign * ((acc * spr * halfTime * (t - dipEnds[i - 1] - halfTime)) - acc * spr * (t - dipEnds[i - 1] - halfTime) * (t - dipEnds[i - 1] - halfTime) / 2.0) / 1000000;
        }
        else {
            return dipHeights[i];
        }
    }
}

// --- reflect ---
// Reflect the momentum portion of the state using given normal and tangent vectors.
// The state is assumed to be an array of 6 elements (first 3: position, next 3: momentum).
// This function uses nextU01() to generate randomness (assumed available).
void reflect(double state[6], const double norm[3], const double tang[3]) {
    double pTarget = std::sqrt(state[3]*state[3] + state[4]*state[4] + state[5]*state[5]);
    double tangPrime[3];
    cross(norm, tang, tangPrime);
    
    double u1 = nextU01();
    double u2 = nextU01();
    
    double theta = std::asin(std::sqrt(u1));
    double phi = 2 * M_PI * u2;
    
    double pN = std::cos(theta);
    double pT = std::sin(theta) * std::cos(phi);
    double pTprime = std::sin(theta) * std::sin(phi);
    
    double newPdir[3];
    newPdir[0] = pN * norm[0] + pT * tang[0] + pTprime * tangPrime[0];
    newPdir[1] = pN * norm[1] + pT * tang[1] + pTprime * tangPrime[1];
    newPdir[2] = pN * norm[2] + pT * tang[2] + pTprime * tangPrime[2];
    
    double pLen = std::sqrt(newPdir[0]*newPdir[0] + newPdir[1]*newPdir[1] + newPdir[2]*newPdir[2]);
    
    state[3] = newPdir[0] * pTarget / pLen;
    state[4] = newPdir[1] * pTarget / pLen;
    state[5] = newPdir[2] * pTarget / pLen;
}

// --- checkDagHit ---
// Returns true if a dagger hit condition is met.
bool checkDagHit(double x, double y, double z, double zOff) {
    double zeta;
    if (x > 0) {
        zeta = 0.5 - std::sqrt(x*x + std::pow(std::fabs(z - zOff) - 1.0, 2));
    } else {
        zeta = 1.0 - std::sqrt(x*x + std::pow(std::fabs(z - zOff) - 0.5, 2));
    }
    if ((x > -0.3524) && (x < 0.0476) && (zeta > 0.0) && (z < (-1.5 + zOff + 0.2))) {
        return true;
    }
    return false;
}

// --- checkClean ---
// Returns 1 or 2 if a cleaning condition is met, 0 otherwise.
int checkClean(const double state[6], const double prevState[6], double cleanHeight) {
    if (((prevState[2] < -1.5 + cleanHeight && state[2] > -1.5 + cleanHeight) ||
         (prevState[2] > -1.5 + cleanHeight && state[2] < -1.5 + cleanHeight))) {
        if (state[1] > 0) {
            return 1;
        }
        if (state[1] > -(0.218041/2 + 0.335121 + 0.3556) &&
            state[1] < -(0.218041/2 + 0.335121) &&
            state[0] > (0.115529/2 + 0.212841 - 0.6604) &&
            state[0] < (0.115529/2 + 0.212841)) {
            return 2;
        }
    }
    return 0;
}

// --- calcDagZeta ---
double calcDagZeta(double x, double y, double z, double zOff) {
    double zeta;
    if (x > 0) {
        zeta = 0.5 - std::sqrt(x*x + std::pow(std::fabs(z - zOff) - 1.0, 2));
    } else {
        zeta = 1.0 - std::sqrt(x*x + std::pow(std::fabs(z - zOff) - 0.5, 2));
    }
    return zeta;
}

// --- checkHouseHitLow ---
bool checkHouseHitLow(double x, double y, double z, double zOff) {
    if (z >= (-1.5 + zOff + 0.2) &&
        z < (-1.5 + zOff + 0.2 + 0.14478) &&
        std::fabs(x + 0.1524) < (0.40 + 2.0179*(z + 1.5 - zOff - 0.2)) / 2.0) {
        return true;
    }
    return false;
}

// --- checkHouseHitHigh ---
bool checkHouseHitHigh(double x, double y, double z, double zOff) {
    if (z >= (-1.5 + zOff + 0.2 + 0.14478) &&
        z < (-1.5 + zOff + 0.2 + 0.2667) &&
        std::fabs(x + 0.1524) < 0.69215/2.0) {
        return true;
    }
    return false;
}

// --- initializeLyapState ---
// Given a reference state (double ref[6]), produce an initial perturbed state in pair[6].
void initializeLyapState(const double ref[6], double pair[6]) {
    double lenP = std::sqrt(ref[3]*ref[3] + ref[4]*ref[4] + ref[5]*ref[5]);
    
    double randomP[3];
    double paraP[3];
    randomP[0] = nextU01();
    randomP[1] = nextU01();
    randomP[2] = nextU01();
    
    paraP[0] = ref[3];
    paraP[1] = ref[4];
    paraP[2] = ref[5];
    
    normalize(randomP);
    normalize(paraP);
    
    double perpP[3];
    cross(paraP, randomP, perpP);
    normalize(perpP);
    
    double alpha = EPSILON*EPSILON*PSCALE*PSCALE / (lenP * 2);
    double beta = std::sqrt(EPSILON*EPSILON*PSCALE*PSCALE -
                   (EPSILON*EPSILON*PSCALE*PSCALE/(2*lenP))*(EPSILON*EPSILON*PSCALE*PSCALE/(2*lenP)));
    
    // Copy positions
    pair[0] = ref[0];
    pair[1] = ref[1];
    pair[2] = ref[2];
    // Update momenta with perturbations
    pair[3] = ref[3] - alpha * paraP[0] + beta * perpP[0];
    pair[4] = ref[4] - alpha * paraP[1] + beta * perpP[1];
    pair[5] = ref[5] - alpha * paraP[2] + beta * perpP[2];
}

// --- resetStates ---
// Adjust pair so that its distance from ref is exactly EPSILON.
void resetStates(const double ref[6], double pair[6]) {
    double dist[6];
    for (int i = 0; i < 6; i++) {
        dist[i] = pair[i] - ref[i];
    }
    double scaling = EPSILON / distance(ref, pair);
    for (int i = 0; i < 6; i++) {
        pair[i] = ref[i] + scaling * dist[i];
    }
}

// --- distance ---
// Compute a weighted Euclidean distance between two 6-element states.
double distance(const double ref[6], const double pair[6]) {
    double dist = 0.0;
    dist += (pair[0]-ref[0])*(pair[0]-ref[0]) / (XSCALE*XSCALE);
    dist += (pair[1]-ref[1])*(pair[1]-ref[1]) / (YSCALE*YSCALE);
    dist += (pair[2]-ref[2])*(pair[2]-ref[2]) / (ZSCALE*ZSCALE);
    dist += (pair[3]-ref[3])*(pair[3]-ref[3]) / (PSCALE*PSCALE);
    dist += (pair[4]-ref[4])*(pair[4]-ref[4]) / (PSCALE*PSCALE);
    dist += (pair[5]-ref[5])*(pair[5]-ref[5]) / (PSCALE*PSCALE);
    return std::sqrt(dist);
}
