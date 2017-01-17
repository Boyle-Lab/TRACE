
#include <stdio.h>
#include <math.h>

#include "hmm_edit.h"
#include "const.h" 


//XXX: use log add or scaling factors???
//didn't include silence state

void ForwardWithScale(HMM *phmm, int T, int *O1, double *O2, double **alpha, double *scale, double *pprob)
{
  int i, j;   /* state indices */
  int t;      /* time index */
  double sum;     /* partial sum */
  //double pNorm;    /* P(O2) ~ N(mu, gamma) */
  
/* 1. Initialization */
  scale[1] = 0.0;
  for (i = 1; i <= phmm->N; i++){
    alpha[1][i] = phmm->pi[i] * ComputeEmission(phmm, i, 1, O1, O2);
    scale[1] += alpha[1][i];
    
  }
  
  for (i = 1; i <= phmm->N; i++){
    alpha[1][i] /= scale[1];
  }
  
  
/* 2. Induction */
 
  for (t = 1; t < T; t++) {
    scale[t+1] = 0.0;
    for (j = 1; j <= phmm->N; j++) {
      sum = 0.0;
      for (i = 1; i <= phmm->N; i++) {
        sum += alpha[t][i] * (phmm->A[i][j]);   
      }
      
      alpha[t+1][j] = sum * ComputeEmission(phmm, j, (t+1), O1, O2);
      scale[t+1] += alpha[t+1][j];
    }
    /* silent state 
    alpha[t+1][1] = 0.0;
    for (i = 2; i <= phmm->N; i++) {
      alpha[t+1][1] += alpha[t+1][i]* (phmm->A[i][1]); //why use t+1 instead of t as others???
    }
    scale[t+1] += alpha[t+1][1];
    */ 
    
    for (j = 1; j <= phmm->N; j++) {
      //fprintf(stdout, "%lf \t", alpha[t+1][j]);
			alpha[t+1][j] /= scale[t+1];
      //fprintf(stdout, "%lf \t", alpha[t+1][j]);
      
    }
    
    //fprintf(stdout, "%lf \t", alpha[t+1][phmm->N]);
  }
 
/* 3. Termination */
  *pprob = 0.0;
  for (t = 1; t <= T; t++){
		*pprob += log(scale[t]);
  }
}


double ComputeEmission(HMM *phmm, int j, int t, int *O1, double *O2) {
  double pNorm;    // P(O2) ~ N(mu, gamma) 
  if (j <= phmm->D){ //first D states are TF motif
    pNorm = (1/(SQRT_TWO_PI * phmm->sigma[1])) * exp( (-0.5) * (1/pow(phmm->sigma[1],2.0))*pow((O2[t] - phmm->mu[1]), 2.0)); 
    //fprintf(stdout, "pNorm: %lf \t", pNorm);
    //return (1.0 - (1.0-pNorm) * (1.0-phmm->pwm[j][O1[t]]));
    return (pNorm * (phmm->pwm[j][O1[t]]));
  }
  if (j > phmm->D){
    pNorm = (1/(SQRT_TWO_PI * (phmm->sigma[j+1-(phmm->D)]))) * exp( (-0.5) * (1/pow(phmm->sigma[j+1-(phmm->D)],2.0))*pow((O2[t] - phmm->mu[j+1-(phmm->D)]), 2.0)); 
    //fprintf(stdout, "pNorm: %lf \t", pNorm);
    //return (1.0 - (1.0-pNorm) * (1.0-phmm->CG[O1[t]]));
    return (pNorm * (phmm->CG[O1[t]]));
  }
  
}
