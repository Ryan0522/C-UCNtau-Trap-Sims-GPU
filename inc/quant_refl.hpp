#ifndef QUANT_REFL_H
#define QUANT_REFL_H

#include <vector>
#include <complex>

#define MASS_N 1.674927471e-27
#define GRAV 9.80665e0
#define JTONEV 6.2415091e27
#define HBAR 1.054571800e-34

#define NBORONB2O3 8.824e27
#define NBORON 1.37e29
#define ABORON -0.1e-15
#define SIGMABORON 2200*3.835e-25

#define NOXYGENB2O3 1.32e28
#define AOXYGEN 5.803e-15
#define SIGMAOXYGEN 4.232e-28
//#define SIGMAOXYGEN 0.0

#define NCARBON 1.133e29
#define ACARBON 6.6460e-15
#define SIGMACARBON 5.551e-28
//#define SIGMACARBON 0.0

#define NZINC 2.527e28
#define AZINC 5.68e-15
#define SIGMAZINC 5.241e-28
//#define SIGMAZINC 0.0

#define NSULFUR 2.527e28
#define ASULFUR 2.847e-15
#define SIGMASULFUR 1.556e-28
//#define SIGMASULFUR 0.0

std::vector<std::complex<double>> matmul(std::vector<std::complex<double>> a, std::vector<std::complex<double>> b);

std::complex<double> k(double ePerp, std::complex<double> u);

std::complex<double> gamma(std::complex<double> kn, std::complex<double> knm1);

std::vector<std::complex<double>> m(std::complex<double> kn, std::complex<double> knm1, double z);

double absorbProbQuantOxide(double ePerp, double thickBoron);

bool absorbMultilayer(double ePerp, double thickBoron, double x, double y, double z, double zOff);

#endif /* QUANT_REFL_H */