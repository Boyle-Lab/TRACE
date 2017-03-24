#include <math.h>
#include "hmm_edit.h"
#include "nrutil.h"
static char rcsid[] = "$Id: viterbi.c,v 1.1 1999/05/06 05:25:37 kanungo Exp kanungo $";

#define VITHUGE  100000000000.0

void ViterbiMultiSeq(HMM *phmm, int T, int *O1, double *O2, double **delta, int **psi,
        int *q, double *vprob, double *pprob, int P, int *peakPos)
{
        int     i, j, k;   /* state indices */
        int     t;      /* time index */
 
        int     maxvalind;
        double  maxval, val;
	double  **biot;

	/* 0. Preprocessing */

	biot = dmatrix(1, phmm->N, 1, T);
	for (i = 1; i <= phmm->N; i++) {
		for (t = 1; t <= T; t++) {
			biot[i][t] = log(ComputeEmission(phmm, i, t, O1, O2));
		}
  } 
  for (k = 1; k <= P; k++){
    /* 1. Initialization  */     
    for (i = 1; i <= phmm->N; i++) {
      delta[peakPos[k]][i] = log(phmm->pi[i]) + biot[i][peakPos[k]];
      psi[peakPos[k]][i] = 0;
    }
 
    /* 2. Recursion */
    for (t = peakPos[k]+1; t < peakPos[k+1]; t++) {
      for (j = 1; j <= phmm->N; j++) {
        maxval = -VITHUGE;
        maxvalind = 1;
        for (i = 1; i <= phmm->N; i++) {
          val = delta[t-1][i] + log(phmm->A[i][j]);
          if (val > maxval) {
            maxval = val;
            maxvalind = i;
          }
        }
        delta[t][j] = maxval + biot[j][t]; 
        psi[t][j] = maxvalind;
 
      }
    }
 
    /* 3. Termination */
    *pprob = -VITHUGE;
    q[peakPos[k+1]-1] = 1;
    for (i = 1; i <= phmm->N; i++) {
      if (delta[peakPos[k+1]-1][i] > *pprob) {
        *pprob = delta[peakPos[k+1]-1][i];
        q[peakPos[k+1]-1] = i;
      }
    }
    vprob[peakPos[k+1]-1] = *pprob; 
	  /* 4. Path (state sequence) backtracking */

  	for (t = peakPos[k+1] - 2; t >= peakPos[k]; t--){
	  	q[t] = psi[t+1][q[t+1]];
      vprob[t] = delta[t+1][q[t+1]];
    }
  }
}
 


  