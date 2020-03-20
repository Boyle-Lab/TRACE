/*
 *  File: logmath.c
 *
 *  math functions and matrix processing functions
 *
 */

#include <math.h>
#include <stdio.h>
#include "logmath.h"
#include "nrutil.h"
#include "hmm.h"
#include <omp.h>

double log_2(const double x){
  return log(x) * D_LOG2E;
}

double _logadd(const double p, const double q) {
  return p + log1p(exp(q - p));
}

double logadd(const double p, const double q) {
  return (p > q) ? _logadd(p, q) : _logadd(q, p);
}

double logCheckAdd(const double p, const double q) {
  if (p != -INFINITY && q != -INFINITY){
    return logadd(p,q);
  }
  else if (p == -INFINITY && q == -INFINITY){
    return -INFINITY;
  }
  else if (p == -INFINITY){
    return q;
  }
  else if (q == -INFINITY){
    return p;
  }
  
}

/*calculate probability from normal distribution*/
double NormDist(double *mean, int j, double *var, double data)
{
  double pNorm;
  pNorm=(1.0/(SQRT_TWO_PI*var[j]))* exp((-0.5/(var[j]*var[j]))*
        ((data-mean[j])*(data-mean[j])));
  return pNorm;
}


/*calculate probability from bivariance normal distribution*/
double BiVarNormDist(double **mean, int j, double **var, 
                     double **corr, double *data)
{
  double pNorm;
  pNorm=(1.0/(SQRT_TWO_PI*SQRT_TWO_PI*var[0][j]*var[1][j]*
        sqrt(1.0-corr[0][j]*corr[0][j])))*
        exp((-0.5/(1-corr[0][j]*corr[0][j]))*(
        (pow((data[0]-mean[0][j]),2.0)/(var[0][j]*var[0][j]))
        +(pow((data[1]-mean[1][j]),2.0)/(var[1][j]*var[1][j]))
        -(2.0*corr[0][j]*(data[0]-mean[0][j])*(data[1]-mean[1][j])
        /(var[0][j]*var[1][j]))));
  if (!(pNorm==pNorm)){
    fprintf(stdout, "BiVarNormDist NAN %d %lf\t",j,corr[0][j]);
  }
  return pNorm;
}

/*calculate probability from multivariance normal distribution*/
double MultiVarNormDist(double **mean, int j, double **covar,
                        int size, double *data) 
{
                        
  double pNorm;
  
  double **diff = dmatrix(size, 1); /*vector of X - mu*/
  double det = determinant(covar, size);
  double **inv = dmatrix(size, size);
  inverse(covar, size, inv);
  int i;

  for (i = 0; i < size; i++) {
    diff[i][0] = data[i] - mean[i][j];
  }
  double **diffTrans = dmatrix(1, size);
  transpose(diff, size, 1, diffTrans);
  double **temp = dmatrix(1, size);
  matrixMultip(1, size, size, diffTrans, inv, temp);
  /*multivariance normal distribution density function*/
  pNorm = exp(-0.5 * matrixMultip_1(size, temp, diff))/
          (pow(SQRT_TWO_PI, size) * sqrt(det));
  
  if (!(pNorm==pNorm)){
    fprintf(stdout, "MultiVarNormDistNAN %lf\t",det);
    printMatrix(stdout,covar, size, size);  
  }
  
  free_dmatrix(diff,size, 1);
  free_dmatrix(inv, size,  size);
  free_dmatrix(temp, 1, size);
  free_dmatrix(diffTrans, 1, size);
  return pNorm;
}

double MultiVarNormDist_2(double **mean, int j, double **inv, double det,
                        int size, double *data) 
{               
  double pNorm;
  double **diff = dmatrix(size, 1); /*vector of X - mu*/
  double **diffTrans = dmatrix(1, size);
  int i;
 
  for (i = 0; i < size; i++) {
    diffTrans[0][i] = diff[i][0] = data[i] - mean[i][j];
  }
  
  double **temp = dmatrix(1, size);
  matrixMultip(1, size, size, diffTrans, inv, temp);
  /*multivariance normal distribution density function*/
  //pNorm = (exp(-0.5 * matrixMultip_1(size, temp, diff))/
  //        (pow(SQRT_TWO_PI, size) * sqrt(det))) * 100000000 + 0.000000000000001; 
                         //pNorm can be extremely small that system cannot compare,
                         //so *100000000 to make sure they won't be the same tiny number
  pNorm = (-0.5 * matrixMultip_1(size, temp, diff))-
          log(pow(SQRT_TWO_PI, size) * sqrt(det));
  if (!(pNorm==pNorm)){
    fprintf(stdout, "MultiVarNormDistNAN %lf\t",det);
    printMatrix(stdout,inv, size, size);  
    exit(0);
  }    
  
  free_dmatrix(diff,size, 1);
  free_dmatrix(temp, 1, size);
  free_dmatrix(diffTrans, 1, size);
 
  return pNorm;
}

void ludcmp(double **a, int n, int *indx, double *d)  
{
  int i,imax,j,k;
  double big,dum,sum,temp;
  double *vv; //vv stores the implicit scaling of each row.
  vv=dvector(n);
  *d=1.0; //No row interchanges yet.
  for (i=0;i<n;i++) { //Loop over rows to get the implicit scaling information.
    big=0.0;
    for (j=0;j<n;j++)
      if ((temp=fabs(a[i][j])) > big) big=temp;
    if (big == 0.0) nrerror("Singular matrix in routine ludcmp"); //No nonzero largest element.
    vv[i]=1.0/big; //Save the scaling.
  }
  for (j=0;j<n;j++) { //This is the loop over columns of Crout�s method.
    for (i=0;i<j;i++) { //except for i = j.
      sum=a[i][j];
      for (k=0;k<i;k++) sum -= a[i][k]*a[k][j];
        a[i][j]=sum;
    }
    big=0.0; //Initialize for the search for largest pivot element.  
    for (i=j;i<n;i++) { //This is i = j of equation (2.3.12) and i = j +1 ...N of equation (2.3.13).
      sum=a[i][j]; 
      for (k=0;k<j;k++)
        sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
      if ( (dum=vv[i]*fabs(sum)) >= big) {//Is the figure of merit for the pivot better than the best so far?
        big=dum;
        imax=i;
      }
    }
    if (j != imax) { //Do we need to interchange rows?
      for (k=0;k<n;k++) { 
        dum=a[imax][k];
        a[imax][k]=a[j][k];
        a[j][k]=dum;
      }
      *d = -(*d); //and change the parity of d.
      vv[imax]=vv[j]; //Also interchange the scale factor.
    }
    indx[j]=imax;
    if (a[j][j] == 0.0) a[j][j]=TINY;
//If the pivot element is zero the matrix is singular (at least to the precision of the
//algorithm). For some applications on singular matrices, it is desirable to substitute
//TINY for zero.
    if (j != n-1) { //Now, finally, divide by the pivot element.
      dum=1.0/(a[j][j]);
      for (i=j+1;i<n;i++) a[i][j] *= dum;
    }
  } //Go back for the next column in the reduction.
  free_dvector(vv,n);
}

void lubksb(double **a, int n, int *indx, double b[])
{
  int i,ii=-1,ip,j;
  double sum;
  for (i=0;i<n;i++) { //When ii is set to a positive value, it will become the
                       //index of the first nonvanishing element of b. We now
                       //do the forward substitution, equation (2.3.6). The
                       //only new wrinkle is to unscramble the permutation as we go.
    ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    if (ii>=0)
      for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
    else if (sum) ii=i; //A nonzero element was encountered, so from now on we
    b[i]=sum;           //will have to do the sums in the loop above.
  }
  for (i=n-1;i>=0;i--) { //Now we do the backsubstitution, equation (2.3.7).
    sum=b[i];
    for (j=i+1;j<n;j++) sum -= a[i][j]*b[j];
    b[i]=sum/a[i][i]; //Store a component of the solution vector X.
  } //All done!
}

void inverse_det_lu(double  **a, int n, double **inv, double *det)
{
  double **y,d,*col;
  int i,j,*indx;
  indx=ivector(n);
  ludcmp(a,n,indx,&d); //This returns d as �1
  col=dvector(n);
  for(j=0;j<n;j++) { //Find inverse by columns.
    d *= a[j][j];
    for(i=0;i<n;i++) col[i]=0.0;
    col[j]=1.0;
    lubksb(a,n,indx,col);
    for(i=0;i<n;i++) inv[i][j]=col[i];
  }
  *det = d;
  
}

double determinant(double **matrix, int size)
{
  
  int i,j,j1,j2;
  double d = 0;
  double **m = NULL;
  double x;
  if (size == 1){
    d = matrix[0][0];
    return(d);
  }
  else if (size==2){
    d=(matrix[0][0]*matrix[1][1])-(matrix[0][1]*matrix[1][0]);
    return(d);
  }
  
  else if (size <= 20){
    if (size > 10){
    printf("det: %d ", size);
    fflush(stdout);
    }
    return(determinant_2(matrix,size));
    
  }
  printf("det: %d ", size);
  fflush(stdout);
  d = 0;
  #pragma omp task shared(d) private(j1, i, j2, m)
         {
    
      for (j1=0;j1<size;j1++) {
         m = malloc((size-1)*sizeof(double *));
         for (i=0;i<size-1;i++)
            m[i] = malloc((size-1)*sizeof(double));
         for (i=1;i<size;i++) {
            j2 = 0;
            for (j=0;j<size;j++) {
               if (j == j1)
                  continue;
               m[i-1][j2] = matrix[i][j];
               j2++;
            }
         }

         d += pow(-1.0,1.0+j1+1.0) * matrix[0][j1] * determinant(m,size-1);
         
         
         for (i=0;i<size-1;i++)
            free(m[i]);
         free(m);
      }
      }
      #pragma omp taskwait
    return(d);
}

double determinant_2(double **matrix, int size)
{
  int i,j,j1,j2;
  double d = 0;
  double **m = NULL;
  
  if (size == 1){
    d = matrix[0][0];
    return(d);
  }
  else if (size==2){
    d=(matrix[0][0]*matrix[1][1])-(matrix[0][1]*matrix[1][0]);
    return(d);
  }
  else {
    d = 0;
      for (j1=0;j1<size;j1++) {
         m = malloc((size-1)*sizeof(double *));
         for (i=0;i<size-1;i++)
            m[i] = malloc((size-1)*sizeof(double));
         for (i=1;i<size;i++) {
            j2 = 0;
            for (j=0;j<size;j++) {
               if (j == j1)
                  continue;
               m[i-1][j2] = matrix[i][j];
               j2++;
            }
         }
         d += pow(-1.0,1.0+j1+1.0) * matrix[0][j1] * determinant(m,size-1);
         for (i=0;i<size-1;i++)
            free(m[i]);
         free(m);
      }
    return(d);
  }
}


void cofactor(double **matrix, int size, double **fac)
{
  double ** b = dmatrix(size-1, size-1);
  double x;
  int p,q,m,n,i,j;
  for (q = 0; q < size; q++) {
   printf("cov q: %d ", q);
   fflush(stdout);
    for (p = 0; p < size; p++) {
    printf(" p: %d ", p);
    fflush(stdout);
      m=0;
      
      for (i = 0; i < size; i++) {
        if (i != q){
          n=0;
          for (j = 0; j < size; j++) {
            if (j != p) {
              b[m][n] = matrix[i][j];
              n++;
            }
          }     
          m++;
        }
        
      }
  #pragma omp parallel
  {
  #pragma omp single
  {
      x = determinant(b,size-1);
  }
  }
  #pragma omp taskwait
      fac[q][p]=pow(-1,q + p) * x;
    }
  }
  free_dmatrix(b, size-1, size-1);    
    
}

/*Finding transpose of matrix of (size1 x size2)*/
void transpose(double **matrix, int size1, int size2, double **trans)
{
  int i,j;
  
  for (i = 0;i < size2; i++)
  {
    for (j = 0; j< size1; j++)
    {
      trans[i][j]= matrix[j][i];
    }
  }
}

/*Finding the inverse matrix*/
void inverse(double **matrix, int size, double **inv)
{
  double ** cof = dmatrix(size, size);
  double ** trans = dmatrix(size, size);
  double det = determinant(matrix, size);
  cofactor(matrix, size, cof);
  transpose(cof, size, size, trans);
  int i,j;
  for (i = 0;i < size; i++)
  {
    for (j = 0; j< size; j++)
    {
      inv[i][j]= trans[i][j]/det;
    }
  }
  free_dmatrix(cof, size, size);    
  free_dmatrix(trans, size, size);   
}
  
double inverse_det(double **matrix, int size, double **inv)
{
  double ** cof = dmatrix(size, size);
  double ** trans = dmatrix(size, size);
  double det;
  #pragma omp parallel
  {
  #pragma omp single
  {
  det = determinant(matrix, size);
  }
  }
  #pragma omp taskwait
  printf("det ");
  fflush(stdout); 
  
  cofactor(matrix, size, cof);
  printf("cof ");
  fflush(stdout); 
  transpose(cof, size, size, trans);
  printf("trans ");
  fflush(stdout); 
  int i,j;
  for (i = 0;i < size; i++)
  {
    for (j = 0; j< size; j++)
    {
      inv[i][j]= trans[i][j]/det;
    }
  }
  free_dmatrix(cof, size, size);    
  free_dmatrix(trans, size, size);   
  printf("inverse_det\t");
  fflush(stdout); 
  return det;
}
  
/*calculate the dot product of two matrix, but not check the size of the two input matrix*/  
/*matrixL: size1 x size2*/
/*matrixR: size2 x size3*/
void matrixMultip(int size1, int size2, int size3, double **matrixL, 
                  double **matrixR, double **product)
{
  int i,j,k;
  for (i = 0; i < size1; i++){
    for (k = 0; k < size3; k++){
      product[i][k] = 0.0;
      for (j = 0; j< size2; j++){
        product[i][k] += (matrixL[i][j] * matrixR[j][k]);
      }
    }
  }
}
  
/*calculate the dot product of two matrix (1 x size) and (size x 1)*/   
double matrixMultip_1(int size, double **matrixL, double **matrixR)
{
  int i;
  double product = 0.0;
  for (i = 0; i < size; i++)
  {  
    product += (matrixL[0][i] * matrixR[i][0]); 
  }
  return product;
}
