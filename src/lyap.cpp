#include "../inc/lyap.hpp"
#include "../setup.h"
#include "../inc/constants.h"
#include "../inc/symplectic.hpp"
#include "../inc/geometry.hpp"
#include <cmath>

lyapResult calcLyap(std::vector<double> ref, double dt, trace tr, double tStart) {
    lyapResult res;
    double t = 0;
    double eStart;
    potential(&ref[0], &ref[1], &ref[2], &(eStart), &t, &tr);
    eStart = eStart - MINU + (ref[3]*ref[3] + ref[4]*ref[4] + ref[5]*ref[5])/(2*MASS_N);
    res.eStart = eStart;
    res.theta = acos(ref[5]/sqrt(ref[3]*ref[3] + ref[4]*ref[4] + ref[5]*ref[5]));
    
    std::vector<double> pair = initializeLyapState(ref);
    double x = 0;
    double y = 0;
    double z = 0;
    
    //Create a trace that will only give the shift at time t.
    shift(&x, &y, &z, tStart, &tr);
    trace trStop;
    trStop.num = 1;
    trStop.x = &x;
    trStop.y = &y;
    trStop.z = &z;
    
    double sumDist = 0;
    
    int numSteps = LYAPTIME/dt;

    double energy;

    for(int i = 0; i < NUMSEP; i++) {
        for(int j = 0; j < numSteps; j++) {
            symplecticStep(pair, dt, energy, tStart, trStop);
            symplecticStep(ref, dt, energy, tStart, trStop);
        }
        double dist = distance(ref, pair);
        sumDist += log(dist/EPSILON);
        resetStates(ref, pair);
    }
    
    double lce = (1.0/(NUMSEP*LYAPTIME))*sumDist;
    res.lce = lce;
    res.eEnd = energy;
    return res;
}
