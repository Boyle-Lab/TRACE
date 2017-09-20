
#include <stdio.h>
#include <math.h>
#include "nrutil.h"
#include "hmm_edit.h"
#include "const.h" 
#include "logmath.h"


void ForwardMultiSeq(HMM *phmm, int T, double **S, double *O2, double **alpha, 
                     double *pprob, int P, int *peakPos)
{
  int i, j, k, n;   /* state indices */
  int t;      /* time index */
  double sum;     /* partial sum */
  double LOGZERO = log(0.0);
  double *Obs = dvector(1, phmm->K);
/* 1. Initialization */
  for (k = 1; k <= P; k++){
    for (i = 1; i <= phmm->N; i++){
      for (n = 1; n <= phmm->M; n++){
        Obs[n] = S[n][t + 1];
      }
      Obs[phmm->K] = O2[t + 1];
      if (phmm->pi[i] == 0.0){
        alpha[peakPos[k]][i] = LOGZERO;
      } 
      else {
        alpha[peakPos[k]][i] = log(phmm->pi[i] * ComputeEmission(phmm, i, Obs));
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
        for (n = 1; n <= phmm->M; n++){
        Obs[n] = S[n][t + 1];
        }
        Obs[phmm->K] = O2[t + 1];
        alpha[t+1][j] = sum + log(ComputeEmission(phmm, j, Obs));
      }
    }
/* 3. Termination */
    pprob[k] = LOGZERO;
    for (i = 1; i <= phmm->N; i++) {
      if (alpha[peakPos[k+1]-1][i] != LOGZERO && alpha[peakPos[k+1]-1][i] 
          == alpha[peakPos[k+1]-1][i]) {
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

