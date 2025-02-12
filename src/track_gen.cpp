#include "../inc/track_gen.hpp"
#include "../inc/constants.h"
#include "../setup.h"
#include <sstream>
#include <cstdlib>
#include <cmath>

#include "../inc/fields_nate.h"

extern "C" {
    #include "../inc/xorshift.h"
}

std::vector<double> readPENTrack(const std::string &s, char delimiter) {
    std::vector<double> state(6);
    std::string value;
    std::istringstream valueStream(s);
    unsigned int index = 0;
    while (std::getline(valueStream, value, delimiter)) {
        if (index < 3) {
            state.at(index) = std::stod(value);
        } else {
            state.at(index) = std::stod(value) * MASS_N;
        }
        index++;
    }
    return state;
}

std::vector<double> randomPointTrapOptimum(trace tr) {
    std::vector<double> state(6);
    double maxEnergy = GRAV*MASS_N*0.3444;
    double maxP = sqrt(2*MASS_N*maxEnergy);
    
    double t = 0.0;    
    
    double energy;
    while(true) {
        energy = maxEnergy * nextU01();
        if(energy < ECUT/JTONEV) {
            continue;
        }
        if(nextU01() < pow(energy/maxEnergy, EPOW)) {
            break;
        }
    }
    
    double totalU;
    do {
        state[2] = -1.464413669130002;
        state[0] = nextU01()*0.15 - 0.075;
        state[1] = nextU01()*0.15 - 0.075;
        potential(&state[0], &state[1], &state[2], &totalU, &t, &tr);
        totalU = totalU - MINU;
    } while(totalU >= energy);
    
    double targetP = sqrt(2.0*MASS_N*(energy - totalU));
    
    double theta;
    while(true) {
        double u1 = nextU01();
        theta = asin(sqrt(u1));
        if(nextU01() < pow(cos(theta), THETAPOW)) {
            break;
        }
    }
    
    double u2 = nextU01();
    double phi = 2 * M_PI * u2;
    
    state[3] = sin(theta)*cos(phi);
    state[4] = sin(theta)*sin(phi);
    state[5] = cos(theta);
    
    double pLen = sqrt(state[3]*state[3] + state[4]*state[4] + state[5]*state[5]);
    
    state[3] = (targetP/pLen)*state[3];
    state[4] = (targetP/pLen)*state[4];
    state[5] = (targetP/pLen)*state[5];
    
    return state;
}

std::vector<double> TRACKGENERATOR(const trace &tr) {
    // Generate a default initial state.
    std::vector<double> state(6);
    state[0] = 0.0;  // x
    state[1] = 0.0;  // y
    state[2] = -2.0; // z (for example)
    state[3] = 1.0;  // momentum in x
    state[4] = 0.0;  // momentum in y
    state[5] = 0.0;  // momentum in z
    return state;
}
