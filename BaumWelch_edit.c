#include <stdio.h> 
#include "nrutil.h"
#include "hmm_edit.h"
#include <math.h>
#include "const.h"

static char rcsid[] = "$Id: baumwelch.c,v 1.6 1999/04/24 15:58:43 kanungo Exp kanungo $";

#define DELTA 0.001 
void BaumWelch(HMM *phmm, int T, int *O1, double *O2, double **alpha, 
  double **beta, double **gamma, int *pniter, 
  double *plogprobinit, double *plogprobfinal)
{
  //fprintf(stdout, "pass1: %d \n", T);
  int	i, j, k;
  int	t, l = 0;
  double	logprobf, logprobb,  threshold;
  double	numeratorA, denominatorA;
  double	mean, variance, pNorm;
  double sumGamma;
  double ***xi, *scale;
  double delta, deltaprev, logprobprev;

  deltaprev = 10e-70;
  //fprintf(stdout, "pass1.0: %d \n", T);
  xi = AllocXi(T, phmm->N);
  scale = dvector(1, T);
  //fprintf(stdout, "pass1.1: %d \n", T);
  
  for (t = 1; t <= (T - 1); t++){
    for (i = 1; i <= phmm->N; i++){
      for (j = 1; j <= phmm->N; j++){
        //fprintf(stdout, "A: %lf \n", phmm->A[i][j]);
        xi[t][i][j] = phmm->A[i][j];
      }
    }
  }
  //fprintf(stdout, "pass2: %d \n", T);
  
  ForwardWithScale(phmm, T, O1, O2, alpha, scale, &logprobf);
  *plogprobinit = logprobf; /* log P(O |intial model) */
  //fprintf(stdout, "pass3: %lf %lf %lf \n", alpha[1][1],alpha[1][2],alpha[1][3]);
  
  BackwardWithScale(phmm, T, O1, O2, beta, scale, &logprobb);
  //fprintf(stdout, "pass4: %lf %lf %lf  \n", beta[1][1],beta[1][2],beta[1][3]);
  
  ComputeGamma(phmm, T, alpha, beta, gamma);
  //fprintf(stdout, "pass5: %d \n", T);
  ComputeXi(phmm, T, O1, O2, alpha, beta, xi);
  //fprintf(stdout, "pass6: %d \n", T);
  logprobprev = logprobf;
 
  /*
  for (t = 1; t <= T; t++){
      fprintf(stdout, "\n t:");
      for (i = 1; i <= phmm->N; i++){
        fprintf(stdout, "alpha:%lf bata:%lf gamma: %lf \n", alpha[t][i], beta[t][i], gamma[t][i]);
      }
  }
  */
  
  do {	
    
		/* frequency of state i in time t=1 never change,
       has to start with silence state*/
		//for (i = 1; i <= phmm->N; i++) 
			//phmm->pi[i] = .001 + .999*gamma[1][i];
    
		/* reestimate transition matrix  and symbol prob in each state */
    mean = 0.0;
    variance = 0.0;
    sumGamma = 0.0;
    
    for (i = 1; i <= phmm->D; i++) { 
      for (t = 1; t <= T; t++) {
        sumGamma += gamma[t][i];
      }
    }
    for (i = 1; i <= phmm->D; i++) { 
			//denominatorA = 0.0;  
      for (t = 1; t <= T; t++) {
				//denominatorA += gamma[t][i];
        mean += (gamma[t][i] / sumGamma) * O2[t];
      }
    }
   /*
    for (i = 1; i <= phmm->D; i++) { 
      
      sumGamma = 0.0;
      for (t = 1; t <= T; t++) {
        sumGamma += gamma[t][i];
      }
      for (t = 1; t <= T; t++) {
				//denominatorA += gamma[t][i];
        mean += (gamma[t][i] / sumGamma) * O2[t];
      }
    } 
    */
    phmm->mu[1] = mean;
    
    for (i = 1; i <= phmm->D; i++) { 
      for (t = 1; t <= T; t++) {
        variance += (gamma[t][i] / sumGamma) * pow((O2[t] - mean), 2.0);
      }
    }
    phmm->sigma[1] = 0.00001 + 0.99999 * sqrt(variance);
    //fprintf(stdout, "%lf %lf \t", phmm->mu[1], phmm->sigma[1]);
    
    denominatorA = 0.0;
    numeratorA = 0.0;
    for (t = 1; t <= (T-1); t++) {
      denominatorA += gamma[t][phmm->D];
      numeratorA += xi[t][phmm->D][phmm->D + 1];
    }
    phmm->A[phmm->D][phmm->D + 1] = numeratorA/denominatorA;
    phmm->A[phmm->D][phmm->N - 1] = 1.0 - phmm->A[phmm->D][phmm->D + 1]; 
    //phmm->A[phmm->D][phmm->N] = 1.0 - phmm->A[phmm->D][phmm->D + 1];
    //fprintf(stdout, "%lf %lf \t", phmm->A[phmm->D][phmm->D + 1], phmm->A[phmm->D][phmm->N - 1]);
    
    for (i > phmm->D; i <= phmm->N; i++) { 
      denominatorA = 0.0;
      mean = 0.0;
      variance = 0.0;
      sumGamma = 0.0;
      for (t = 1; t <= T; t++) {
        sumGamma += gamma[t][i];
      }
      for (t = 1; t <= (T-1); t++) {
        denominatorA += gamma[t][i];
        mean += (gamma[t][i] / sumGamma) * O2[t]; 
      }
      mean += (gamma[T][i] / sumGamma) * O2[T];
      for (j = 1; j <= phmm->N; j++) {
        numeratorA = 0.0;
        for (t = 1; t <= (T - 1); t++) {
          numeratorA += xi[t][i][j];
        }
        phmm->A[i][j] = numeratorA/denominatorA;
      }
      phmm->mu[i- phmm->D + 1] = mean;
      for (t = 1; t <= T; t++) {
        variance += (gamma[t][i] / sumGamma) * pow((O2[t] - mean), 2.0);
      } 
      phmm->sigma[i- phmm->D + 1] = 0.00001 + 0.99999 * sqrt(variance);
      //fprintf(stdout, "%lf %lf \t", phmm->mu[i- phmm->D + 1], phmm->sigma[i- phmm->D + 1]);
    }
    ForwardWithScale(phmm, T, O1, O2, alpha, scale, &logprobf);
    BackwardWithScale(phmm, T, O1, O2, beta, scale, &logprobb);
    ComputeGamma(phmm, T, alpha, beta, gamma);
    ComputeXi(phmm, T, O1, O2, alpha, beta, xi);
   /*
    fprintf(stdout, "\n run:");
    for (t = 1; t <= T; t++){
      fprintf(stdout, "\n t:");
      for (i = 1; i <= phmm->N; i++){
        fprintf(stdout, "alpha:%lf bata:%lf gamma: %lf \t", alpha[t][i], beta[t][i], gamma[t][i]);
      }
    }
    */
		/* compute difference between log probability of 
		   two iterations */
    delta = logprobf - logprobprev; 
    logprobprev = logprobf;
    l++;	
	}
  while (delta > DELTA); /* if log probability does not 
                                  change much, exit */ 
 
  *pniter = l;
  *plogprobfinal = logprobf; /* log P(O|estimated model) */
  FreeXi(xi, T, phmm->N);
  free_dvector(scale, 1, T);
}

//XXX:may need to change to i and i-1, not sure

void ComputeGamma(HMM *phmm, int T, double **alpha, double **beta, 
	double **gamma)
{

  int i, j;
  int	t;
  double denominator;

  for (t = 1; t <= T; t++) {
    denominator = 0.0;
    for (j = 1; j <= phmm->N; j++) {
      gamma[t][j] = alpha[t][j]*beta[t][j];
      denominator += gamma[t][j];
    }
    //fprintf(stdout, "gamma: ");
    for (i = 1; i <= phmm->N; i++) {
      gamma[t][i] /= denominator;
      //fprintf(stdout, "ab: %lf %lf \t", alpha[t][i], beta[t][i]);
      //fprintf(stdout, "%lf \t",gamma[t][i]);
    }
  }
}

void ComputeXi(HMM* phmm, int T, int *O1, double *O2, double **alpha, double **beta, 
	double ***xi)
{
  int i, j;
  int t;
  double sum;
  double pNorm;

  for (t = 1; t <= (T - 1); t++) {
    sum = 0.0;	
    for (i = 1; i <= phmm->N; i++){ 
      for (j = 1; j <= phmm->N; j++) {
        if (xi[t][i][j] != 0.0){
          //fprintf(stdout, "t,i,j: %d %d  %d \n", t,i,j);
          xi[t][i][j] = alpha[t][i]*beta[t+1][j] *(phmm->A[i][j]) * ComputeEmission(phmm, j, (t+1), O1, O2);
          sum += xi[t][i][j];
        }
      }
    }

    for (i = 1; i <= phmm->N; i++) {
      for (j = 1; j <= phmm->N; j++){ 
        xi[t][i][j] /= sum;
        //fprintf(stdout, "xi: %lf \t", xi[t][i][j]);
      }
    }
  }
}

double *** AllocXi(int T, int N)
{
  int t;
  double ***xi;

  xi = (double ***) malloc(T*sizeof(double **));

  xi --;

  for (t = 1; t <= T; t++){
    xi[t] = dmatrix(1, N, 1, N);
  }
  return xi;
}

void FreeXi(double *** xi, int T, int N)
{
  int t;

  for (t = 1; t <= T; t++) {
    free_dmatrix(xi[t], 1, N, 1, N);
  }
  xi ++;
  free(xi);

}
