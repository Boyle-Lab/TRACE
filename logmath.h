#ifndef _LOGMATH_H_
#define _LOGMATH_H_
#include <math.h>
#include <stdio.h>

double _logadd(const double p, const double q);

double logadd(const double p, const double q); 

double logCheckAdd(const double p, const double q);

double MultiVarNormDist(double **mean, int j, double **covar, 
                        int size, double *data);

double determinant(double **matrix, int size);

void cofactor(double **matrix, int size, double **fac);

void transpose(double **matrix, int size1, int size2, double **trans);

void inverse(double **matrix, int size, double **inv);

void matrixMultip(int size1, int size2, int size3, double **matrixL, 
                  double **matrixR, double **product);

double matrixMultip_1(int size, double **matrixL, double **matrixR);

#endif