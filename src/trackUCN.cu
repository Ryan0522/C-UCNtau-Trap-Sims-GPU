#include "../inc/trackUCN.hpp"
#include "../inc/constants.h"
#include "../inc/symplectic.hpp"
#include "../inc/geometry.hpp"
#include "../inc/quant_refl.hpp"
#define _USE_MATH_DEFINES
#include <cmath>
#include <assert.h>

#include "../setup.h"

#include "../inc/fields_nate.h"

extern "C" {
    #include "../inc/xorshift.h"
}

#ifdef __CUDACC__
#define HD __host__ __device__
#else
#define HD
#endif

HD fixedResult fixedEffDaggerHitTime(double state[6], double dt, const trace &tr) {
    double tang[3]    = {0.0, 0.0, 1.0};
    double normPlus[3]= {0.0, 1.0, 0.0};
    double normMinus[3]={0.0, -1.0, 0.0};
    fixedResult res;
    double pMag = sqrt(state[3]*state[3] + state[4]*state[4] + state[5]*state[5]);
    res.theta = acos(state[5] / pMag);
    double t = 0;
    double pot;
    potential(&state[0], &state[1], &state[2], &pot, &t, (trace *)&tr);
    res.eStart = pot - MINU + (state[3]*state[3] + state[4]*state[4] + state[5]*state[5])/(2*MASS_N);
    
    double deathTime = 9999;
    double settlingTime = CLEANINGTIME;
    res.settlingT = settlingTime;
    
    double prevState[6];
    int numSteps = (int)(settlingTime/dt);
    bool defect = false;
    double energy;
    for (int i = 0; i < numSteps; i++) {
        for (int j = 0; j < 6; j++) prevState[j] = state[j];
        defect = symplecticStep_Defect(state, dt, energy, t, tr);
        if (((prevState[2] < -1.5 + CLEANINGHEIGHT && state[2] > -1.5 + CLEANINGHEIGHT) ||
             (prevState[2] > -1.5 + CLEANINGHEIGHT && state[2] < -1.5 + CLEANINGHEIGHT)) &&
            (state[1] > 0)) {
            res.energy = energy;
            res.t = t - settlingTime;
            res.ePerp = state[5]*state[5] / (2*MASS_N);
            res.x = state[0];
            res.y = state[1];
            res.z = state[2];
            res.zOff = -2;
            res.nHit = 0;
            res.nHitHouseLow = 0;
            res.nHitHouseHigh = 0;
            res.deathTime = deathTime;
            return res;
        }
        t += dt;
    }
    
    int nHit = 0;
    int nHitHouseLow = 0;
    int nHitHouseHigh = 0;
    while (true) {
        for (int j = 0; j < 6; j++) prevState[j] = state[j];
        defect = symplecticStep_Defect(state, dt, energy, t, tr);
        if (defect) {
            nHitHouseHigh += 1;
        }
        t += dt;
        if (t - settlingTime > deathTime) {
            res.energy = energy;
            res.t = t - settlingTime;
            res.ePerp = state[4]*state[4] / (2*MASS_N);
            res.x = state[0];
            res.y = state[1];
            res.z = state[2];
            res.zOff = -1;
            res.nHit = nHit;
            res.nHitHouseLow = nHitHouseLow;
            res.nHitHouseHigh = nHitHouseHigh;
            res.deathTime = deathTime;
            return res;
        }
        if (((prevState[2] < -1.5 + RAISEDCLEANINGHEIGHT && state[2] > -1.5 + RAISEDCLEANINGHEIGHT) ||
             (prevState[2] > -1.5 + RAISEDCLEANINGHEIGHT && state[2] < -1.5 + RAISEDCLEANINGHEIGHT)) &&
            (state[1] > 0)) {
            res.energy = energy;
            res.t = t - settlingTime;
            res.ePerp = state[5]*state[5] / (2*MASS_N);
            res.x = state[0];
            res.y = state[1];
            res.z = state[2];
            res.zOff = -3;
            res.nHit = nHit;
            res.nHitHouseLow = nHitHouseLow;
            res.nHitHouseHigh = nHitHouseHigh;
            res.deathTime = deathTime;
            return res;
        }
        if (isnan(energy)) {
            res.energy = energy;
            res.t = t - settlingTime;
            res.ePerp = 0.0;
            res.x = state[0];
            res.y = state[1];
            res.z = state[2];
            res.zOff = -4;
            res.nHit = nHit;
            res.nHitHouseLow = nHitHouseLow;
            res.nHitHouseHigh = nHitHouseHigh;
            res.deathTime = deathTime;
            return res;
        }
        if ((prevState[1] < 0 && state[1] > 0) || (prevState[1] > 0 && state[1] < 0)) {
            double fracTravel = fabs(prevState[1]) / (fabs(state[1]) + fabs(prevState[1]));
            double predX = prevState[0] + fracTravel * (state[0] - prevState[0]);
            double predZ = prevState[2] + fracTravel * (state[2] - prevState[2]);
            
            double zOff = zOffDipCalc(t - settlingTime);
            
            if (checkDagHit(predX, 0.0, predZ, zOff)) {
                nHit += 1;
                if (((predX > -0.3302 && predX < -0.2794) ||
                     (predX > -0.2286 && predX < -0.1778) ||
                     (predX > -0.127 && predX < -0.0762) ||
                     (predX > -0.0254 && predX < 0.0254)) &&
                    absorbMultilayer((state[3]*state[3]+state[4]*state[4]+state[5]*state[5]) *
                                     cos(nextU01() * M_PI / 2) / (2*MASS_N),
                                     BTHICK, predX, 0.0, predZ, zOff)) {
                    res.energy = energy;
                    res.t = t - settlingTime;
                    res.ePerp = state[4]*state[4] / (2*MASS_N);
                    res.x = predX;
                    res.y = 0.0;
                    res.z = predZ;
                    res.zOff = zOff;
                    res.nHit = nHit;
                    res.nHitHouseLow = nHitHouseLow;
                    res.nHitHouseHigh = nHitHouseHigh;
                    res.deathTime = deathTime;
                    return res;
                } else {
                    if (nextU01() <= 0.0003) {
                        res.energy = energy;
                        res.t = t - settlingTime;
                        res.ePerp = state[4]*state[4] / (2*MASS_N);
                        res.x = predX;
                        res.y = 0.0;
                        res.z = predZ;
                        res.zOff = -5;
                        res.nHit = nHit;
                        res.nHitHouseLow = nHitHouseLow;
                        res.nHitHouseHigh = nHitHouseHigh;
                        res.deathTime = deathTime;
                        return res;
                    }
                }
                if (prevState[1] > 0 && prevState[4] < 0)
                    reflect(prevState, normPlus, tang);
                else
                    reflect(prevState, normMinus, tang);
                for (int j = 0; j < 6; j++) state[j] = prevState[j];
            }
            else if (checkHouseHitLow(predX, 0.0, predZ, zOff)) {
                nHitHouseLow += 1;
                if (nextU01() <= 0.0003) {
                    res.energy = energy;
                    res.t = t - settlingTime;
                    res.ePerp = state[4]*state[4] / (2*MASS_N);
                    res.x = predX;
                    res.y = 0.0;
                    res.z = predZ;
                    res.zOff = -6;
                    res.nHit = nHit;
                    res.nHitHouseLow = nHitHouseLow;
                    res.nHitHouseHigh = nHitHouseHigh;
                    res.deathTime = deathTime;
                    return res;
                }
                if (prevState[1] > 0 && prevState[4] < 0)
                    reflect(prevState, normPlus, tang);
                else
                    reflect(prevState, normMinus, tang);
                for (int j = 0; j < 6; j++) state[j] = prevState[j];
            }
            else if (checkHouseHitHigh(predX, 0.0, predZ, zOff)) {
                if (nextU01() <= 0.0003) {
                    res.energy = energy;
                    res.t = t - settlingTime;
                    res.ePerp = state[4]*state[4] / (2*MASS_N);
                    res.x = predX;
                    res.y = 0.0;
                    res.z = predZ;
                    res.zOff = -7;
                    res.nHit = nHit;
                    res.nHitHouseLow = nHitHouseLow;
                    res.nHitHouseHigh = nHitHouseHigh;
                    res.deathTime = deathTime;
                    return res;
                }
                if (prevState[1] > 0 && prevState[4] < 0)
                    reflect(prevState, normPlus, tang);
                else
                    reflect(prevState, normMinus, tang);
                for (int j = 0; j < 6; j++) state[j] = prevState[j];
            }
        }
    }
}