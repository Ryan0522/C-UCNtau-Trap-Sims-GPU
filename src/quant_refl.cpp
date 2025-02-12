#include "../inc/quant_refl.hpp"
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>
#include <complex>
#include <assert.h>
#include "../setup.h"

extern "C" {
    #include "../inc/xorshift.h"
}

std::vector<std::complex<double>> matmul(std::vector<std::complex<double>> a, std::vector<std::complex<double>> b) {
    std::vector<std::complex<double>> res(4);
    assert(a.size() == 4 && b.size() == 4);
    res[0] = a[0]*b[0] + a[1]*b[2];
    res[1] = a[0]*b[1] + a[1]*b[3];
    res[2] = a[2]*b[0] + a[3]*b[2];
    res[3] = a[2]*b[1] + a[3]*b[3];
    return res;
}

std::complex<double> k(double ePerp, std::complex<double> u) {
    return std::sqrt((2*MASS_N/(HBAR*HBAR))*(ePerp - u));
}

std::complex<double> gamma(std::complex<double> kn, std::complex<double> knm1) {
    return knm1/kn;
}

std::vector<std::complex<double>> m(std::complex<double> kn, std::complex<double> knm1, double z) {
    std::vector<std::complex<double>> res(4);
    res[0] = (1.0/2.0)*(1.0 + gamma(kn,knm1))*std::exp(std::complex<double>(0,1)*(knm1-kn)*z);
    res[1] = (1.0/2.0)*(1.0 - gamma(kn,knm1))*std::exp(-std::complex<double>(0,1)*(knm1+kn)*z);
    res[2] = (1.0/2.0)*(1.0 - gamma(kn,knm1))*std::exp(std::complex<double>(0,1)*(knm1+kn)*z);
    res[3] = (1.0/2.0)*(1.0 + gamma(kn,knm1))*std::exp(-std::complex<double>(0,1)*(knm1-kn)*z);
    return res;
}

double absorbProbQuantOxide(double ePerp, double thickBoron) {
//    const double voxide = (2*M_PI*(HBAR*HBAR)/MASS_N)*ABORON*NBORONB2O3 + (2*M_PI*(HBAR*HBAR)/MASS_N)*AOXYGEN*NOXYGENB2O3;
//    const double woxide = (HBAR/2)*NBORONB2O3*SIGMABORON + (HBAR/2)*NOXYGENB2O3*SIGMAOXYGEN;
    const double vboron = (2*M_PI*(HBAR*HBAR)/MASS_N)*ABORON*NBORON;
    const double wboron = (HBAR/2)*NBORON*SIGMABORON;
    const double vzns = (2*M_PI*(HBAR*HBAR)/MASS_N)*AZINC*NZINC + (2*M_PI*(HBAR*HBAR)/MASS_N)*ASULFUR*NSULFUR;
    const double wzns = (HBAR/2)*NZINC*SIGMAZINC + (HBAR/2)*NSULFUR*SIGMASULFUR;
    
    std::vector<std::complex<double>> pots = {std::complex<double>(0, 0),
//                                              std::complex<double>(vcarbon, -wcarbon),
//                                              std::complex<double>(voxide, -woxide),
                                              std::complex<double>(vboron, -wboron),
                                              std::complex<double>(vzns, -wzns)};
    std::vector<std::complex<double>> mbar = {std::complex<double>(1,0), std::complex<double>(0,0), std::complex<double>(0,0), std::complex<double>(1,0)};
    std::vector<double> zs = {0.0, thickBoron*1e-9, 10000e-9};
    
    for(int i = pots.size()-1; i > 0; i--) {
        mbar = matmul(mbar, m(k(ePerp, pots[i]), k(ePerp, pots[i-1]), zs[i-1]));
    }
    
    return 1.0 - (std::conj(-mbar[2]/mbar[3])*-mbar[2]/mbar[3]).real();
}

bool absorbMultilayer(double ePerp, double thickBoron, double x, double y, double z, double zOff) {
    double zeta;
    if(x > 0) {
        zeta = 0.5 - sqrt(x*x + pow(fabs(z - zOff) - 1.0, 2));
    }
    else {
        zeta = 1.0 - sqrt(x*x + pow(fabs(z - zOff) - 0.5, 2));
    }
    double boronProb = zeta > ZETACUT ? 1 : zeta/ZETACUT;
    double u = nextU01();
    if(u < boronProb*absorbProbQuantOxide(ePerp, thickBoron)) {
        return true;
    }
    return false; 
}