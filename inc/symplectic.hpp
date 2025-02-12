#ifndef SYMPLECTIC_HPP
#define SYMPLECTIC_HPP

#include <vector>
#include <cmath>

#include "../inc/fields_nate.h"

#ifdef __CUDACC__
#define HD __host__ __device__
#else
#define HD
#endif

HD void symplecticStep(double state[6], double deltaT, double &energy, double &t, const trace &tr);

HD bool symplecticStep_Defect(double state[6], double deltaT, double &energy, double &t, const trace &tr);

#endif // SYMPLECTIC_HPP