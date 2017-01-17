#include <stdio.h>
#include "hmm_edit.h"
#include "const.h"

void BackwardWithScale(HMM *phmm, int T, int *O1, double *O2, double **beta, double *scale, double *pprob)
{
  int     i, j;   /* state indices */
  int     t;      /* time index */
  double sum;
  double pNorm;
 
  /* 1. Initialization */
  for (i = 1; i <= phmm->N; i++){
    //beta[T][i] = 1.0 / scale[T]; //why using the same scale factor as forward???
    beta[T][i] = 1.0 / phmm->N;
  }            
  /* silent state 
  for (j = 2; j <= phmm->N; i++) {
    if (j <= (phmm->D + 1)){
      pNorm = (1/(SQRT_TWO_PI * phmm->mu[2]))*exp(-0.5 * (1/pow(phmm->mu[2],2.0))*pow((O2[T] - phmm->mu[2]), 2.0)); 
      beta[T][1] += beta[T][j] * (phmm->A[1][j]) * pNorm * phmm->pwm[j][O1[T]];
    }
    if (j > (phmm->D + 1)){
      pNorm = (1/(SQRT_TWO_PI * phmm->mu[j+1-phmm->D]))*exp(-0.5 * (1/pow(phmm->mu[j+1-phmm->D],2.0))*pow((O2[T] - phmm->mu[j+1-phmm->D]), 2.0)); 
      beta[T][1] += beta[T][j] * (phmm->A[1][j]) * pNorm * phmm->CG[O1[T]];
    }
  }
    beta[T][1] /= scale[T]; 
  */
  
  /* 2. Induction */
 
  for (t = (T - 1); t >= 1; t--) {
    scale[t] = 0.0;
    for (i = 1; i <= phmm->N; i++) {
      beta[t][i] = 0.0;
      for (j = 1; j <= phmm->N; j++) {
        beta[t][i] += beta[t+1][j] * (phmm->A[i][j]) * ComputeEmission(phmm, j, (t + 1), O1, O2);
      }
      //beta[t][i] /= scale[t];
      scale[t] += beta[t][i];
      //fprintf(stdout, "beta: %lf \t", beta[t][i]);
    }
    
    for (i = 1; i <= phmm->N; i++) {
      beta[t][i] /= scale[t];
      //fprintf(stdout, "beta: %lf \t", beta[t][i]);
    } 
    
    
    /* silent state 
    beta[t][1] = 0.0;
    for (j = 2; j <= phmm->N; i++) {
      if (j <= (phmm->D + 1)){
      pNorm = (1/(SQRT_TWO_PI * phmm->mu[2]))*exp(-0.5 * (1/pow(phmm->mu[2],2.0))*pow((O2[t] - phmm->mu[2]), 2.0)); 
      beta[t][1] += beta[t][j] * (phmm->A[j][1]) * pNorm * phmm->pwm[j][O1[t]];
      }
      if (j > (phmm->D + 1)){
        pNorm = (1/(SQRT_TWO_PI * phmm->mu[j+1-phmm->D]))*exp(-0.5 * (1/pow(phmm->mu[j+1-phmm->D],2.0))*pow((O2[t] - phmm->mu[j+1-phmm->D]), 2.0)); 
        beta[t][1] += beta[t][j] * (phmm->A[j][1]) * pNorm * phmm->CG[O1[t]];
      }
      //beta[t][1] += beta[t][j]* (phmm->A[j][1]) * (phmm->B[j][O[t]]); //why same t???
    }
    
    //scale[t] += beta[t][1]
    beta[t][1] /= scale[t];
    */
  }
 
  /* 3. Termination */
  //*pprob = 0.0;
  //for (i = 1; i <= phmm->N; i++) {
    //*pprob += beta[1][i];
  //}
}