#include "../inc/symplectic.hpp"
#include "../inc/constants.h"

#include "../setup.h"

#include <vector>
#include <math.h>

#include "../inc/fields_nate.h"

extern "C" {
    #include "../inc/xorshift.h"
}

HD void symplecticStep(double state[6], double deltaT, double &energy, double &t, const trace &tr) {
    double fx, fy, fz, totalU;
    int n = 0;
    
    const double a[4] = {
        0.5153528374311229364,
        -0.085782019412973646,
        0.4415830236164665242,
        0.1288461583653841854,
    };
    const double b[4] = {
        0.1344961992774310892,
        -0.2248198030794208058,
        0.7563200005156682911,
        0.3340036032863214255,
    };

    force(&state[0], &state[1], &state[2], &fx, &fy, &fz, &totalU, &t, (trace *)&tr);
    energy = totalU - MINU + (state[3]*state[3] + state[4]*state[4] + state[5]*state[5])/(2.0*MASS_N);

    bool can_have_defect = ((totalU - GRAV*MASS_N*state[2]) * 10000 / MU_N) >= 155.340528314;

    state[3] = state[3] + b[n]*fx*deltaT;
    state[4] = state[4] + b[n]*fy*deltaT;
    state[5] = state[5] + b[n]*fz*deltaT;
    state[0] = state[0] + a[n]*state[3]*deltaT/MASS_N;
    state[1] = state[1] + a[n]*state[4]*deltaT/MASS_N;
    state[2] = state[2] + a[n]*state[5]*deltaT/MASS_N;
    t = t + a[n]*deltaT;

    for(n = 1; n < 4; n++) {
        force(&state[0], &state[1], &state[2], &fx, &fy, &fz, &totalU, &t, (trace *)&tr);
        state[3] = state[3] + b[n]*fx*deltaT;
        state[4] = state[4] + b[n]*fy*deltaT;
        state[5] = state[5] + b[n]*fz*deltaT;
        state[0] = state[0] + a[n]*state[3]*deltaT/MASS_N;
        state[1] = state[1] + a[n]*state[4]*deltaT/MASS_N;
        state[2] = state[2] + a[n]*state[5]*deltaT/MASS_N;
        t = t + a[n]*deltaT;
    }
}

HD bool symplecticStep_Defect(double state[6], double deltaT, double &energy, double &t, const trace &tr) {
    double fx, fy, fz, totalU;
    int n = 0;
    
    const double a[4] = {
        0.5153528374311229364,
        -0.085782019412973646,
        0.4415830236164665242,
        0.1288461583653841854,
    };
    const double b[4] = {
        0.1344961992774310892,
        -0.2248198030794208058,
        0.7563200005156682911,
        0.3340036032863214255,
    };

    force(&state[0], &state[1], &state[2], &fx, &fy, &fz, &totalU, &t, (trace *)&tr);
    energy = totalU - MINU + (state[3]*state[3] + state[4]*state[4] + state[5]*state[5])/(2.0*MASS_N);

    bool can_have_defect = ((totalU - GRAV*MASS_N*state[2]) * 10000 / MU_N) >= 155.340528314;

    state[3] = state[3] + b[n]*fx*deltaT;
    state[4] = state[4] + b[n]*fy*deltaT;
    state[5] = state[5] + b[n]*fz*deltaT;
    state[0] = state[0] + a[n]*state[3]*deltaT/MASS_N;
    state[1] = state[1] + a[n]*state[4]*deltaT/MASS_N;
    state[2] = state[2] + a[n]*state[5]*deltaT/MASS_N;
    t = t + a[n]*deltaT;

    for(n = 1; n < 4; n++) {
        force(&state[0], &state[1], &state[2], &fx, &fy, &fz, &totalU, &t, (trace *)&tr);
        state[3] = state[3] + b[n]*fx*deltaT;
        state[4] = state[4] + b[n]*fy*deltaT;
        state[5] = state[5] + b[n]*fz*deltaT;
        state[0] = state[0] + a[n]*state[3]*deltaT/MASS_N;
        state[1] = state[1] + a[n]*state[4]*deltaT/MASS_N;
        state[2] = state[2] + a[n]*state[5]*deltaT/MASS_N;
        t = t + a[n]*deltaT;
    }

    if (can_have_defect) {
        if (nextU01() <= DEFECT) {
            double magnitude = sqrt(state[3]*state[3] + state[4]*state[4] + state[5]*state[5]);
            double phi = nextU01() * M_PI;
            double theta = nextU01() * 2 * M_PI;
            state[3] = magnitude * sin(phi) * cos(theta);
            state[4] = magnitude * sin(phi) * sin(theta);
            state[5] = magnitude * cos(phi);
            return true;
        }
    }
    return false;
}