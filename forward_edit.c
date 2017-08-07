
#include <stdio.h>
#include <math.h>

#include "hmm_edit.h"
#include "const.h" 
#include "logmath.h"


void ForwardMultiSeq(HMM *phmm, int T, int *O1, double *O2, double *S, double **alpha, double *pprob, int P, int *peakPos)
{
  int i, j, k;   /* state indices */
  int t;      /* time index */
  double sum;     /* partial sum */
  double LOGZERO = log(0.0);
  
/* 1. Initialization */
  for (k = 1; k <= P; k++){
    for (i = 1; i <= phmm->N; i++){
      if (phmm->pi[i] == 0.0){
        alpha[peakPos[k]][i] = LOGZERO;
      } 
      else {
        alpha[peakPos[k]][i] = log(phmm->pi[i] * ComputeEmission(phmm, i, peakPos[k], S, O2));
      }
   }
  
/* 2. Induction */
 
    for (t = peakPos[k]; t < peakPos[k+1]-1; t++) {
      for (j = 1; j <= phmm->N; j++) {
        //fprintf(stdout, "%d,%d\t", t,j);
        sum = LOGZERO;
        for (i = 1; i <= phmm->N; i++) {
          if (alpha[t][i] != LOGZERO && log(phmm->A[i][j])!= LOGZERO) {
            if (sum != LOGZERO) {
              sum = logadd(sum, alpha[t][i] + log(phmm->A[i][j]));  
            }
            else{
              sum =  alpha[t][i] + log(phmm->A[i][j]);
            }
          }
        }
        alpha[t+1][j] = sum + log(ComputeEmission(phmm, j, (t+1), S, O2));
      }
    }
/* 3. Termination */
    pprob[k] = LOGZERO;
    for (i = 1; i <= phmm->N; i++) {
      if (alpha[peakPos[k+1]-1][i] != LOGZERO && alpha[peakPos[k+1]-1][i] == alpha[peakPos[k+1]-1][i]) {
        if (pprob[k] != LOGZERO) {
          pprob[k] = logadd(pprob[k], alpha[peakPos[k+1]-1][i]);
        }
        else {
          pprob[k] = alpha[peakPos[k+1]-1][i];
        }
      }
      
    }
  }
}

double ComputeEmission(HMM *phmm, int j, int t, double *S, double *O2) {
  double pNorm;    /* multivariance normal distribution */ 

  pNorm = (1/((pow(SQRT_TWO_PI, 2.0) * phmm->sigma[j][1] * phmm->sigma[j][2]) * pow((1 - pow(phmm->rho[j], 2.0)), 0.5))) * 
	  exp( (-0.5) / (1 - pow(phmm->rho[j], 2.0))*(pow((S[t] - phmm->mu[j][1]), 2.0)/pow(phmm->sigma[j][1],2.0)+
	  pow((O2[t] - phmm->mu[j][2]), 2.0)/pow(phmm->sigma[j][2],2.0)-2*phmm->rho[j]*(S[t] - phmm->mu[j][1])*(O2[t] - 
	  phmm->mu[j][2])/(phmm->sigma[j][1]*phmm->sigma[j][2]))); 
  //fprintf(stdout, "pNorm: %lf \t", pNorm);  
  return pNorm;
}

/*
double ComputeEmission(HMM *phmm, int j, int t, int *O1, double *O2) {
  double pNorm;    // P(O2) ~ N(mu, gamma) 
  if (j <= phmm->D){ //first D states are TF motif
    pNorm = (1/(SQRT_TWO_PI * phmm->sigma[1])) * exp( (-0.5) * (1/pow(phmm->sigma[1],2.0))*pow((O2[t] - phmm->mu[1]), 2.0)); 
    return (pNorm * (phmm->pwm[j][O1[t]]));
  } 
  
  if (j > (phmm->D)){
    pNorm = (1/(SQRT_TWO_PI * (phmm->sigma[j+1-(phmm->D)]))) * exp((-0.5) * (1/pow(phmm->sigma[j+1-(phmm->D)],2.0))*pow((O2[t] - phmm->mu[j+1-(phmm->D)]), 2.0)); 
    return (pNorm * (phmm->CG[O1[t]]));
  }
  
}
*/
