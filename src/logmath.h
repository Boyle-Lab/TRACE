#ifndef _LOGMATH_H_
#define _LOGMATH_H_
#include <math.h>
#include <stdio.h>

double log_2(const double x);

double _logadd(const double p, const double q);

double logadd(const double p, const double q); 

double logCheckAdd(const double p, const double q);
double NormDist(double *mean, int j, double *var, double data);
double BiVarNormDist(double **mean, int j, double **var, double **corr, double *data);
double MultiVarNormDist(double **mean, int j, double **covar, 
                        int size, double *data);
double MultiVarNormDist_2(double **mean, int j, double **inv, double det, 
                        int size, double *data);

double determinant(double **matrix, int size);
double determinant_2(double **matrix, int size);

void cofactor(double **matrix, int size, double **fac);

void transpose(double **matrix, int size1, int size2, double **trans);

void inverse(double **matrix, int size, double **inv);
double inverse_det(double **matrix, int size, double **inv);

void matrixMultip(int size1, int size2, int size3, double **matrixL, 
                  double **matrixR, double **product);

double matrixMultip_1(int size, double **matrixL, double **matrixR);

void ludcmp(double **a, int n, int *indx, double *d);  
void lubksb(double **a, int n, int *indx, double b[]);
void inverse_det_lu(double  **a, int n, double **inv, double *det);
#endif
