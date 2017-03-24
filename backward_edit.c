#include <stdio.h>
#include "hmm_edit.h"
#include "const.h"
#include "logmath.h"

void BackwardMultiSeq(HMM *phmm, int T, int *O1, double *O2, double **beta, double *pprob, int P, int *peakPos)
{
  int     i, j, k;   /* state indices */
  int     t;      /* time index */
  double sum;
  double LOGZERO = log(0.0);
  
  /* 1. Initialization */
  for (k = 1; k <= P; k++){
    for (i = 1; i <= phmm->N; i++){
      beta[peakPos[k+1]-1][i] = LOGZERO;
    }            
    beta[peakPos[k+1]-1][phmm->N-1] = 0.0;
  
  /* 2. Induction */
    for (t = peakPos[k+1]-2; t >= peakPos[k]; t--) {
      for (i = 1; i <= phmm->N; i++) {
        beta[t][i] = LOGZERO;
        for (j = 1; j <= phmm->N; j++) {
          if (beta[t+1][j] != LOGZERO && phmm->A[i][j] != LOGZERO) {
            if (beta[t][i] != LOGZERO) {
              beta[t][i] = logadd(beta[t][i], beta[t+1][j] + log(phmm->A[i][j]) + log(ComputeEmission(phmm, j, (t + 1), O1, O2)));  
            }
            else{
              beta[t][i] =  beta[t+1][j] + log(phmm->A[i][j]) + log(ComputeEmission(phmm, j, (t + 1), O1, O2));
            }
          }
        }
      //fprintf(stdout, "beta: %lf \t", beta[t][i]);
      }
    } 
  }
 
  /* 3. Termination */
  //*pprob = 0.0;
  //for (i = 1; i <= phmm->N; i++) {
    //*pprob += beta[1][i];
  //}
}
