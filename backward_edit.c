#include <stdio.h>
#include "hmm_edit.h"
#include "const.h"
#include "logmath.h"
#include "nrutil.h"

void BackwardMultiSeq(HMM *phmm, int T, double **S, double *O2, double **beta, 
                      double *pprob, int P, int *peakPos)
{
  int i, j, k, n;   /* state indices */
  int t;      /* time index */
  double sum;
  double LOGZERO = log(0.0);
  double *Obs = dvector(1, phmm->K);
  /* 1. Initialization */
  for (k = 1; k <= P; k++){
    for (i = 1; i <= phmm->N; i++){
      beta[peakPos[k+1]-1][i] = LOGZERO;
    }            
    beta[peakPos[k+1]-1][phmm->N-1] = 0.0;
    //fprintf(stdout, "beta: %d  ", k);
  /* 2. Induction */
    for (t = peakPos[k+1]-2; t >= peakPos[k]; t--) {
      for (n = 1; n <= phmm->M; n++){
        Obs[n] = S[n][t + 1];
      }
      Obs[phmm->K] = O2[t + 1];
      //fprintf(stdout, "beta: %lf %lf  ", Obs[1], Obs[2]);
      for (i = 1; i <= phmm->N; i++) {
        beta[t][i] = LOGZERO;
        for (j = 1; j <= phmm->N; j++) {
          if (beta[t+1][j] != LOGZERO && phmm->A[i][j] != LOGZERO) {
            
            if (beta[t][i] != LOGZERO) {
              beta[t][i] = logadd(beta[t][i], beta[t+1][j] + log(phmm->A[i][j]) + log(ComputeEmission(phmm, j, Obs)));  
            }
            else{
              beta[t][i] = beta[t+1][j] + log(phmm->A[i][j]) + log(ComputeEmission(phmm, j, Obs));
            }
          }
        }
      }
    } 
  }
 
  /* 3. Termination */
  //*pprob = 0.0;
  //for (i = 1; i <= phmm->N; i++) {
    //*pprob += beta[1][i];
  //}
}
