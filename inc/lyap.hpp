#ifndef LYAP_H
#define LYAP_H

#include <vector>

typedef struct lyapResult {
    float eStart;
    float eEnd;
    float theta;
    float lce;
} lyapResult;

#include "../inc/fields_nate.h"

lyapResult calcLyap(std::vector<double> ref, double dt, trace tr, double tStart = 0);

#endif /* LYAP_H */