#include <math.h>
#include <stdio.h>
#include "logmath.h"
#include "nrutil.h"
#include "hmm_edit.h"
#include "const.h"


double _logadd(const double p, const double q) {
  return p + log1p(exp(q - p));
}

double logadd(const double p, const double q) {
  return (p > q) ? _logadd(p, q) : _logadd(q, p);
}

double logCheckAdd(const double p, const double q) {
  if (p == log(0.0)){
    return q;
  }
  else{
    return logadd(p,q);
  }
}


double MultiVarNormDist(double **mean, int j, double **covar, 
                        int size, double *data) 
{
                        
  double pNorm;
  int i;
  double **diff = dmatrix(1, size, 1, 1); /*vector of X - mu*/
  double det = determinant(covar, size);
  double **inv = dmatrix(1, size, 1, size);
  inverse(covar, size, inv);
  for (i = 1; i <= size; i++) {
    diff[i][1] = data[i] - mean[i][j];
  }
  double **diffTrans = dmatrix(1, 1, 1, size);
  transpose(diff, size, 1, diffTrans);
  double **temp = dmatrix(1, 1, 1, size);
  matrixMultip(1, size, size, diffTrans, inv, temp);
  /*multivariance normal distribution density function*/
  pNorm = exp(-0.5 * matrixMultip_1(size, temp, diff))/(pow(SQRT_TWO_PI, size) * pow(det, 0.5));
  /*
                        for (i = 1; i <= size; i++) {
                          fprintf(stdout,"%lf\t", data[i]);
                        }
                        for (i = 1; i <= size; i++) {
                          fprintf(stdout,"%lf\t", mean[i][j]);
                        }
                        printMatrix(covar, size, size);*/
  return pNorm;
}

  
double determinant(double **matrix, int size)
{
  int pr,j,p,q,t;
  double ** b = dmatrix(1, size, 1, size);
  double * c = dvector(1, size);
  double d = 0.0;
  if (size == 1){
    d = matrix[1][1];
    return(d);
  }
  else if (size==2){
    d=(matrix[1][1]*matrix[2][2])-(matrix[1][2]*matrix[2][1]);
    return(d);
  }
  else {
    for(j = 1; j <= size; j++) {       
      int r = 1, s = 1;
      for(p = 1; p <= size; p++) {
        for(q = 1; q <= size; q++) {
          if(p != 1 && q != j) {
            b[r][s] = matrix[p][q];
            s++;
            if(s > size -1) {
              r++;
              s=1;
            }
          }
        }
      }
      for(t = 1, pr = 1; t <= (1+j); t++) {
        pr = (-1) * pr;
      }
      c[j] = pr * determinant(b,size-1);
    }
    for(j = 1, d = 0.0; j <= size; j++) {
      d=d+(matrix[1][j]*c[j]);
    }
    return(d);
  }
}

void cofactor(double **matrix, int size, double **fac)
{
  double ** b = dmatrix(1, size-1, 1, size-1);
  //double ** fac = dmatrix(1, size, 1, size);
 
  int p,q,m,n,i,j;
  for (q = 1; q <= size; q++) {
    for (p = 1; p <= size; p++) {
      m=1;
      
      for (i = 1; i <= size; i++) {
        if (i != q){
          n=1;
          for (j = 1; j <= size; j++) {
            if (j != p) {
              b[m][n] = matrix[i][j];
              n++;
            }
          }     
          m++;
        }
        
      }
      fac[q][p]=pow(-1,q + p) * determinant(b,size-1);
    }
  }
      
    
}


/*Finding transpose of matrix of (size1 x size2)*/
void transpose(double **matrix, int size1, int size2, double **trans)
{
  int i,j;
  
  for (i = 1;i <= size2; i++)
  {
    for (j = 1; j<= size1; j++)
    {
      trans[i][j]= matrix[j][i];
    }
  }
}

/*Finding the inverse matrix*/
void inverse(double **matrix, int size, double **inv)
{
  double ** cof = dmatrix(1, size, 1, size);
  double ** trans = dmatrix(1, size, 1, size);
  double det = determinant(matrix, size);
  cofactor(matrix, size, cof);
  transpose(cof, size, size, trans);
  int i,j;
  for (i = 1;i <= size; i++)
  {
    for (j = 1; j<= size; j++)
    {
      inv[i][j]= trans[i][j]/det;
    }
  }
}
  
//void matrixMulti(double **matrix, int size, double **prod)  
  
/*calculate the dot product of two matrix, but not check the size of the two input matrix*/  
/*matrixL: size1 x size2*/
/*matrixR: size2 x size3*/
void matrixMultip(int size1, int size2, int size3, double **matrixL, 
                  double **matrixR, double **product)
{
  int i,j,k;
  for (k = 1; k <= size3; k++){
    for (i = 1; i <= size1; i++)
    {  
      product[i][k] = 0.0;
      for (j = 1; j<= size2; j++)
      {
        product[i][k] += matrixL[i][j] * matrixR[j][k];
      }
    }
  }
}
  
/*calculate the dot product of two matrix (1 x size) and (size x 1)*/   
double matrixMultip_1(int size, double **matrixL, double **matrixR)
{
  int i;
  double product = 0.0;
  for (i = 1; i <= size; i++)
  {  
    product += matrixL[1][i] * matrixR[i][1]; 
  }
  return product;
}