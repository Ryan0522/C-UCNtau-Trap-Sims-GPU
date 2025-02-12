#ifndef UCNTFIELDS_F_H 
#define UCNTFIELDS_F_H

#define SAMPDT 0.0004

typedef struct trace {
    double* x;
    double* y;
    double* z;
    int num;
} trace;

#ifdef __CUDACC__
#define HD __host__ __device__
#else
#define HD
#endif

HD void shift(double *x, double *y, double *z, double t, trace* tr);
HD void force(double *x_in, double *y_in, double *z_in, double *fx, double *fy, double *fz, double *totalU, double* t, trace* tr);
HD void fieldstrength(double *x_in, double *y_in, double *z_in, double *totalB, double* t, trace* tr);
HD void potential(double *x_in, double *y_in, double *z_in, double *totalU, double* t, trace* tr);

#endif /* UCNTFIELDS_F_H */
