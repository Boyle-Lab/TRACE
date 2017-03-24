
#include <stdio.h>
#include <math.h>

#include "hmm_edit.h"
#include "const.h" 
#include "logmath.h"

//XXX: use log add or scaling factors???
//didn't include silence state


void ForwardMultiSeq(HMM *phmm, int T, int *O1, double *O2, double **alpha, double *pprob, int P, int *peakPos)
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
        alpha[peakPos[k]][i] = log(phmm->pi[i] * ComputeEmission(phmm, i, peakPos[k], O1, O2));
      }
    }
  
/* 2. Induction */
 
    for (t = peakPos[k]; t < peakPos[k+1]-1; t++) {
      for (j = 1; j <= phmm->N; j++) {
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
      
        alpha[t+1][j] = sum + log(ComputeEmission(phmm, j, (t+1), O1, O2));
      }
    }

/* 3. Termination */
    pprob[k] = LOGZERO;
    for (i = 1; i <= phmm->N; i++) {
      if (alpha[peakPos[k+1]-1][i] != LOGZERO) {
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
