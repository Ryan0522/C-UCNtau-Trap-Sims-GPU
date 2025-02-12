#include "../inc/fields_nate.h"
#include "../setup.h"
#include <math.h>
#include <stdio.h>

#define B_HOLD 0.005           // holding field strength
#define B_REM 1.35             // remnant magnet field strength
#define MAG_THICK 0.0254       // thickness of layer of PM array
#define MAG_SPACE 0.05114      // characteristic spacing of magnets
#define MU_N -9.6623647e-27
#define MASS_N 1.674927471e-27
#define GRAV 9.80665e0
#define N_TERMS 3              // how far out to go in field ripple expansion
#define KAPPA 1000


HD void shift(double *x, double *y, double *z, double t, trace* tr) {
    int iLow = (int)(t / SAMPDT);
    double frac = (t - iLow * SAMPDT) / (SAMPDT);
    int iHi = iLow + 1;
    iLow = (iLow % tr->num + tr->num) % tr->num;
    iHi = (iHi % tr->num + tr->num) % tr->num;

    *x = *x + HEATMULT * (tr->x[iLow] + frac * (tr->x[iHi] - tr->x[iLow]));
    *y = *y + HEATMULT * (tr->y[iLow] + frac * (tr->y[iHi] - tr->y[iLow]));
    *z = *z + HEATMULT * (tr->z[iLow] + frac * (tr->z[iHi] - tr->z[iLow]));
}

HD void force(double *x_in, double *y_in, double *z_in, double *fx, double *fy, double *fz, double *totalU, double* t, trace* tr) {
    double A = 4 * B_REM / (M_PI * sqrt(2.0));

    double x = *x_in;
    double y = *y_in;
    double z = *z_in;
    // Apply the shift to x,y,z
    shift(&x, &y, &z, *t, tr);
    double z_grav = *z_in;

    double gx = 0.0, gy = 0.0, gz = 0.0, R, r, Rprime, rprime;
    
    R = 0.5 + 0.5 / (1 + exp(-KAPPA * x));
    r = 1.0 - 0.5 / (1 + exp(-KAPPA * x));
    Rprime = 0.5 * KAPPA * (1.0 - 1.0 / (1 + exp(-KAPPA * x))) / (1 + exp(-KAPPA * x));
    rprime = -Rprime;

    double rho = sqrt(y * y + z * z);
    double r_zeta = sqrt((rho - R) * (rho - R) + x * x);

    if (z < -1.0 && r_zeta < r) {
        double eta = r * atan(x / (sqrt(y * y + z * z) - R));
        double zeta = r - sqrt(x * x + (sqrt(y * y + z * z) - R) * (sqrt(y * y + z * z) - R));
        double sum_cos = 0.0, sum_sin = 0.0, sum_k_cos = 0.0, sum_k_sin = 0.0;
        double cos_term = 0.0, sin_term = 0.0;
        double k_n;

        for (int n = 1; n <= N_TERMS; n++) {
            k_n = 2 * M_PI * (4.0 * n - 3.0) / MAG_SPACE;
            
            double sign = (n % 2 == 0) ? 1.0 : -1.0;
            cos_term = sign / (4.0 * n - 3.0) * (1 - exp(-k_n * MAG_THICK)) * exp(-k_n * zeta) * cos(k_n * eta);
            sin_term = sign / (4.0 * n - 3.0) * (1 - exp(-k_n * MAG_THICK)) * exp(-k_n * zeta) * sin(k_n * eta);
            
            sum_cos += cos_term;
            sum_k_cos += k_n * cos_term;
            sum_sin += sin_term;
            sum_k_sin += k_n * sin_term;
        }
        
        double b_zeta = A * sum_cos;
        double b_eta  = A * sum_sin;
        double b_hold = B_HOLD * (r + R) / sqrt(y * y + z * z);
        double b_tot = sqrt(b_zeta * b_zeta + b_eta * b_eta + b_hold * b_hold);
    
        double d_BZeta[2] = { -A * sum_k_cos, -A * sum_k_sin };
        double d_BEta[2] = { -A * sum_k_sin,  A * sum_k_cos };

        double d_Bh[3] = { 0.0,
                           -B_HOLD * (r + R) * y / pow(y * y + z * z, 3.0 / 2.0),
                           -B_HOLD * (r + R) * z / pow(y * y + z * z, 3.0 / 2.0) };
        double d_Zeta[3] = { rprime + (Rprime * (sqrt(y*y + z*z) - R) - x) / 
                              sqrt(x*x + (R - sqrt(y*y+z*z)) * (R - sqrt(y*y+z*z))),
                             -y * (sqrt(y*y+z*z) - R) / (sqrt(y*y+z*z) *
                              sqrt(x*x + (R - sqrt(y*y+z*z))*(R - sqrt(y*y+z*z)))),
                             -z * (sqrt(y*y+z*z) - R) / (sqrt(y*y+z*z) *
                              sqrt(x*x + (R - sqrt(y*y+z*z))*(R - sqrt(y*y+z*z))) ) };
        double d_Eta[3] = { -rprime * atan(x / (R - sqrt(y*y+z*z))) - (r * (R - sqrt(y*y+z*z) - x*Rprime)) / 
                             (x*x + y*y + z*z + R*R - 2*R*sqrt(y*y+z*z)),
                             -r * x * y / (sqrt(y*y+z*z) * pow((R - sqrt(y*y+z*z)), 2) *
                             (1.0 + x*x / pow((R - sqrt(y*y+z*z)), 2))),
                             -r * x * z / (sqrt(y*y+z*z) * pow((R - sqrt(y*y+z*z)), 2) *
                             (1.0 + x*x / pow((R - sqrt(y*y+z*z)), 2))) };
        
        gx = MU_N / b_tot * (b_zeta * (d_BZeta[0] * d_Zeta[0] + d_BZeta[1] * d_Eta[0])
                             + b_eta  * (d_BEta[0] * d_Zeta[0] + d_BEta[1] * d_Eta[0]));
        
        gy = MU_N / b_tot * (b_zeta * (d_BZeta[0] * d_Zeta[1] + d_BZeta[1] * d_Eta[1])
                             + b_eta  * (d_BEta[0] * d_Zeta[1] + d_BEta[1] * d_Eta[1])
                             + b_hold * d_Bh[1]);
        
        gz = MU_N / b_tot * (b_zeta * (d_BZeta[0] * d_Zeta[2] + d_BZeta[1] * d_Eta[2])
                             + b_eta  * (d_BEta[0] * d_Zeta[2] + d_BEta[1] * d_Eta[2])
                             + b_hold * d_Bh[2]);
        
        *totalU = -MU_N * b_tot + GRAV * MASS_N * z_grav;
        gz -= GRAV * MASS_N;
    } else {
        gx = NAN;
        gy = NAN;
        gz = NAN;
        *totalU = NAN;
    }
    
    *fx = gx;
    *fy = gy;
    *fz = gz;
}

HD void fieldstrength(double *x_in, double *y_in, double *z_in, double *totalB, double* t, trace* tr) {
    double A = 4 * B_REM / (M_PI * sqrt(2.0));
    double x = *x_in;
    double y = *y_in;
    double z = *z_in;
    shift(&x, &y, &z, *t, tr);
    double z_grav = *z_in;
    
    double R, r, Rprime, rprime;
    R = 0.5 + 0.5 / (1 + exp(-KAPPA * x));
    r = 1.0 - 0.5 / (1 + exp(-KAPPA * x));
    Rprime = 0.5 * KAPPA * (1.0 - 1.0 / (1 + exp(-KAPPA * x))) / (1 + exp(-KAPPA * x));
    rprime = -Rprime;
    
    double rho = sqrt(y*y + z*z);
    double r_zeta = sqrt((rho-R)*(rho-R) + x*x);
    
    if (z < -1.0 && r_zeta < r) {
        double eta = r * atan(x / (sqrt(y*y + z*z) - R));
        double zeta = r - sqrt(x*x + (sqrt(y*y + z*z) - R)*(sqrt(y*y + z*z) - R));
        double sum_cos = 0.0, sum_sin = 0.0;
        double cos_term = 0.0, sin_term = 0.0;
        double k_n;
        
        for (int n = 1; n <= N_TERMS; n++) {
            k_n = 2 * M_PI * (4.0*n - 3.0) / MAG_SPACE;
            cos_term = (n % 2 == 0 ? 1.0 : -1.0) / (4.0*n - 3.0) *
                       (1 - exp(-k_n * MAG_THICK)) * exp(-k_n * zeta) * cos(k_n * eta);
            sin_term = (n % 2 == 0 ? 1.0 : -1.0) / (4.0*n - 3.0) *
                       (1 - exp(-k_n * MAG_THICK)) * exp(-k_n * zeta) * sin(k_n * eta);
            sum_cos += cos_term;
            sum_sin += sin_term;
        }
        double b_zeta = A * sum_cos;
        double b_eta = A * sum_sin;
        double b_hold = B_HOLD * (r + R) / sqrt(y*y + z*z);
        double b_tot = sqrt(b_zeta*b_zeta + b_eta*b_eta + b_hold*b_hold);
        *totalB = b_tot;
    } else {
        *totalB = NAN;
    }
}

HD void potential(double *x_in, double *y_in, double *z_in, double *totalU, double* t, trace* tr) {
    double A = 4 * B_REM / (M_PI * sqrt(2.0));
    double x = *x_in;
    double y = *y_in;
    double z = *z_in;
    shift(&x, &y, &z, *t, tr);
    double z_grav = *z_in;
    
    double R, r, Rprime, rprime;
    R = 0.5 + 0.5 / (1 + exp(-KAPPA * x));
    r = 1.0 - 0.5 / (1 + exp(-KAPPA * x));
    Rprime = 0.5 * KAPPA * (1.0 - 1.0 / (1 + exp(-KAPPA * x))) / (1 + exp(-KAPPA * x));
    rprime = -Rprime;
    
    double rho = sqrt(y*y + z*z);
    double r_zeta = sqrt((rho-R)*(rho-R) + x*x);
    
    if (z < -1.0 && r_zeta < r) {
        double eta = r * atan(x / (sqrt(y*y + z*z) - R));
        double zeta = r - sqrt(x*x + (sqrt(y*y+z*z) - R) * (sqrt(y*y+z*z) - R));
        double sum_cos = 0.0, sum_sin = 0.0;
        double cos_term = 0.0, sin_term = 0.0;
        double k_n;
        
        for (int n = 1; n <= N_TERMS; n++) {
            k_n = 2 * M_PI * (4.0*n - 3.0) / MAG_SPACE;
            cos_term = (n % 2 == 0 ? 1.0 : -1.0) / (4.0*n - 3.0) *
                       (1 - exp(-k_n * MAG_THICK)) * exp(-k_n * zeta) * cos(k_n * eta);
            sin_term = (n % 2 == 0 ? 1.0 : -1.0) / (4.0*n - 3.0) *
                       (1 - exp(-k_n * MAG_THICK)) * exp(-k_n * zeta) * sin(k_n * eta);
            sum_cos += cos_term;
            sum_sin += sin_term;
        }
        double b_zeta = A * sum_cos;
        double b_eta = A * sum_sin;
        double b_hold = B_HOLD * (r + R) / sqrt(y*y + z*z);
        double b_tot = sqrt(b_zeta*b_zeta + b_eta*b_eta + b_hold*b_hold);
        
        *totalU = -MU_N * b_tot + GRAV * MASS_N * z_grav;
    } else {
        *totalU = NAN;
    }
}
