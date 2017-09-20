#include <stdio.h> 
#include "nrutil.h"
#include "hmm_edit.h"
#include <math.h>
#include "const.h"
#include "logmath.h"

static char rcsid[] = "$Id: baumwelch.c,v 1.6 1999/04/24 15:58:43 kanungo Exp kanungo $";

#define DELTA 0.001 
#define LOGZERO log(0.0)


void BaumWelchBVN(HMM *phmm, int T, int *O1, double *O2, double **alpha, 
  double **beta, double **gamma, int *pniter, int P, int *peakPos)
{
  int	i, j, k, n, m;
  int	t, l = 0;
  
  /*the following are temp values used when calcaulating variables in BW*/
  double	*logprobf, *logprobb, plogprobinit, plogprobfinal;
  double	numeratorA, denominatorA;
  double	*mean, variance;
  double	cov; /*covariance*/
  double sumGamma, tempSum, tempNumeratorA;
  double ***xi; /*the probability of being in state i at time t and in state j at time t + 1*/
  double delta, deltaprev, logprobprev, totalProb;
  int numNonZero;
  
  fprintf(stdout, "passa");
  
  deltaprev = 10e-70;
  xi = AllocXi(T, phmm->N);
  
  /*initialize xi values*/
  for (t = 1; t <= (T - 1); t++){
    for (i = 1; i <= phmm->N; i++){
      for (j = 1; j <= phmm->N; j++){
        xi[t][i][j] = phmm->A[i][j];
      }
    }
  }

  double **S = dmatrix(1, phmm->M, 1, T);/*containing PWM score for each TF*/
  CalMotifScore(phmm, S, O1, P, peakPos);
  
  logprobf = logprobb = dvector(1,P); /*vector containing log likelihood 
                                        for each peak*/
  //double **Obs = (double **)malloc(phmm->K*sizeof(double*));
  //double **Obs = dmatrix(1, phmm->K, 1, T);
  //for (i = 1; i <= phmm->M; i++){
    //Obs[i] = S[i];
  //}
  //Obs[phmm->K] = O2;/*only using 1+N known and unchangeable data vectors now,
                      //mignt need to change the structure in the future*/
  //free_dvector(S, 1, T);  
  
  fprintf(stdout, "passb");
  
  //double **ObsTrans = dmatrix(1, T, 1, phmm->K);
  //transpose(Obs, phmm->K, T, ObsTrans);
  //free_dmatrix(Obs, 1, phmm->K, 1, T);
  //printMatrix(ObsTrans, 10,phmm->K);
  ForwardMultiSeq(phmm, T, S, O2, alpha, logprobf, P, peakPos);
  //printMatrix(alpha, 10,phmm->N);
  BackwardMultiSeq(phmm, T, S, O2, beta, logprobb, P, peakPos);
  //printMatrix(alpha, 10,phmm->N);
  //printMatrix(beta, 10,phmm->N);
  ComputeGammaWithLog(phmm, T, alpha, beta, gamma);
  //printMatrix(gamma, 10,phmm->N);
  ComputeXiWithLog(phmm, T, S, O2, alpha, beta, xi);
  //printMatrix(xi[1], phmm->N,phmm->N);
  logprobprev = 0.0;
  for (i = 1; i <= P; i ++) {
    logprobprev += logprobf[i];
  }
  mean = dvector(1, phmm->K);
  //variance = dvector(1, phmm->K);

  do {
/* reestimate transition matrix  and symbol prob in each state */
    for (i = 1; i <= phmm->N; i++) { 
      denominatorA = LOGZERO; /*denominator when calculating aij, 
                               which is sum of gamma(i)*/
      for (n = 1; n <= phmm->K; n++){
        mean[n] = 0.0;
        //variance[n] = 0.0;
      }
      sumGamma = LOGZERO; /* sum of gamma used to calculate mean and variance 
                             but not aij*/
      numNonZero = 0;
      for (k = 1; k <= P; k++) {
        tempSum = LOGZERO;
        for (t = peakPos[k]; t < peakPos[k+1] - 1; t++) {
          if (gamma[t][i] != LOGZERO){
            tempSum = logCheckAdd(tempSum, gamma[t][i]);
          }
        }
        denominatorA = logCheckAdd(denominatorA, tempSum);
        /*the last base at each peak is used when calculating 
          mean and variance but now aij*/
        if (gamma[peakPos[k+1]-1][i] != LOGZERO){
          sumGamma = logCheckAdd(sumGamma, logCheckAdd(tempSum, 
                                 gamma[peakPos[k+1]-1][i]));
        }
        else {
          sumGamma = logCheckAdd(sumGamma, tempSum);
        }
      }
      /*calculate mean of each data for each hidden state*/
      for (k = 1; k <= P; k++) {
        for (t = peakPos[k]; t < peakPos[k+1]; t++) {
				//denominatorA += gamma[t][i];
          for (n = 1; n <= phmm->M; n++){
            mean[n] += exp(gamma[t][i] - sumGamma) * S[n][t];
          }
          mean[phmm->K] += exp(gamma[t][i] - sumGamma) * O2[t];
        }
      }
      
      for (j = 1; j <= phmm->N; j++) {
        if (phmm->A[i][j] != 0.0){
          numNonZero ++;
        }
      }
      for (j = 1; j <= phmm->N; j++) {
        numeratorA = LOGZERO; /*nominator when calculate aij,
                                which is the sum of xi(i,j)*/
        for (k = 1; k <= P; k++) {
          tempNumeratorA = LOGZERO;
          for (t = peakPos[k]; t < peakPos[k+1] - 1; t++) {
            if (xi[t][i][j] != LOGZERO){
              tempNumeratorA = logCheckAdd(tempNumeratorA, xi[t][i][j]);
            }
            
          }
          numeratorA = logCheckAdd(numeratorA, tempNumeratorA);
          
        }
        if (phmm->A[i][j] != 0.0){
          phmm->A[i][j] = (0.999999 * exp(numeratorA - denominatorA) + 0.000001) 
                           / (0.999999 + 0.000001 * numNonZero);  
        }
        
      }
      for (n = 1; n <= phmm->K; n++){
        phmm->mu[n][i] = mean[n];
      }
      //printMatrix(phmm->mu, phmm->K, phmm->N);
      /*calculate variance fisrt*/
      for (n = 1; n <= phmm->K; n++) {
        cov = 0.0;
          for (k = 1; k <= P; k++) {
            for (t = peakPos[k]; t < peakPos[k+1]; t++) {
              if (n <= phmm->M){
                cov += exp(gamma[t][i] - sumGamma) * pow((S[n][t] 
                       - mean[n]), 2.0);
              }
              else {
                cov += exp(gamma[t][i] - sumGamma) * pow((O2[t] 
                       - mean[n]), 2.0);
              }
            }
          }
          phmm->sigma[n][i] = 0.00001 + 0.99999 * sqrt(cov);   
      }
      /*then calculate covariance and correlation*/
      for (n = 1; n <= phmm->K - 1; n++){
        for (m = n + 1; m <= phmm->K; m++){
          cov = 0.0;   
          for (k = 1; k <= P; k++) {
            for (t = peakPos[k]; t < peakPos[k+1]; t++) {
              if (m <= phmm->M){
                cov += exp(gamma[t][i] - sumGamma) * (S[n][t] - mean[n]) 
                       * (S[m][t] - mean[m]);
              }
              else {
                cov += exp(gamma[t][i] - sumGamma) * (S[n][t] - mean[n]) 
                       * (O2[t] - mean[m]);
              }
            }
          }
          /*covariance matrix[n][m] should match to rho[phmm->K*(n-1)+m-n*(n+1)/2][n] */
          phmm->rho[phmm->K*(n-1)+m-n*(n+1)/2][i] = cov / (phmm->sigma[n][i] 
                                                    * phmm->sigma[m][i]); 
        }
      }
    }
    ForwardMultiSeq(phmm, T, S, O2, alpha, logprobf, P, peakPos);
    BackwardMultiSeq(phmm, T, S, O2, beta, logprobb, P, peakPos);
    ComputeGammaWithLog(phmm, T, alpha, beta, gamma);
    ComputeXiWithLog(phmm, T, S, O2, alpha, beta, xi);
    
    PrintHMM(stdout, phmm);
    /* compute difference between log probability of two iterations */

    totalProb = 0.0;
    
    for (i = 1; i <= P; i ++) {
      if (logprobf[i] != log(0.0)){
        totalProb += logprobf[i];
      }
    }
    
    fprintf(stdout, "\n %d %lf %lf", l, totalProb, logprobprev);
    delta = totalProb - logprobprev; 
    logprobprev = totalProb;
    l++;	
    fprintf(stdout, "\n %d %lf", l, delta);
    
    
	}
  while (fabs(delta) > DELTA); /* if log probability does not 
                                  change much, exit */ 
  
  *pniter = l;
  //*plogprobfinal = totalProb; /* log P(O|estimated model) */
  FreeXi(xi, T, phmm->N);
}
	
void ComputeGammaWithLog(HMM *phmm, int T, double **alpha, double **beta, 
	                       double **gamma)
{

  int i, j;
  int	t;
  double denominator;

  for (t = 1; t <= T; t++) {
    denominator = LOGZERO;
    for (j = 1; j <= phmm->N; j++) {
      gamma[t][j] = alpha[t][j]+beta[t][j];
      if (gamma[t][j] != LOGZERO){
        if (denominator != LOGZERO){
          denominator = logadd(denominator, gamma[t][j]);
        }
        else{
          denominator = gamma[t][j];
        }
      }
    }
    //fprintf(stdout, "gamma: ");
    for (i = 1; i <= phmm->N; i++) {
      gamma[t][i] -= denominator;
      //fprintf(stdout, "ab: %lf %lf \t", alpha[t][i], beta[t][i]);
      //fprintf(stdout, "%lf \t",gamma[t][i]);
    }
  }
}


void ComputeXiWithLog(HMM* phmm, int T, double **S, double *O2,
                      double **alpha, double **beta, double ***xi)
{
  int i, j, k;
  int t;
  double sum;
  double pNorm;
  double *Obs = dvector(1, phmm->K);
  for (t = 1; t <= (T - 1); t++) {
    sum = -INFINITY;
    for (i = 1; i <= phmm->N; i++){ 
      for (j = 1; j <= phmm->N; j++) {
        for (k = 1; k <= phmm->M; k++){
          //fprintf(stdout, "t,i,j: %d %d  %d \n", t,i,j);
          Obs[k] = S[k][t + 1];
        }
        Obs[phmm->K] = O2[t + 1];
        xi[t][i][j] = alpha[t][i] + beta[t+1][j] + log(phmm->A[i][j]) + 
                      log(ComputeEmission(phmm, j, Obs));
        if (xi[t][i][j] != LOGZERO){
          if (sum != LOGZERO){
            sum = logadd(sum, xi[t][i][j]);
          }
          else{
            sum = xi[t][i][j];
          }
        }
      }
    }
    for (i = 1; i <= phmm->N; i++) {
      for (j = 1; j <= phmm->N; j++){ 
        xi[t][i][j] -= sum;
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

/*
void BaumWelchMultiSeq(HMM *phmm, int T, int *O1, double *O2, double **alpha, 
  double **beta, double **gamma, int *pniter, int P, int *peakPos)
{
  //fprintf(stdout, "pass1: %d \n", T);
  int	i, j, k;
  int	t, l = 0;
  double	*logprobf, *logprobb, plogprobinit, plogprobfinal;
  double	numeratorA, denominatorA;
  double	mean, variance, pNorm;
  double sumGamma, tempSum, tempNumeratorA;
  double ***xi, *scale;
  double delta, deltaprev, logprobprev, totalProb;
  int numNonZero;
  
  deltaprev = 10e-70;
  xi = AllocXi(T, phmm->N);
  
  for (t = 1; t <= (T - 1); t++){
    for (i = 1; i <= phmm->N; i++){
      for (j = 1; j <= phmm->N; j++){
        xi[t][i][j] = phmm->A[i][j];
      }
    }
  }
  logprobf = logprobb = dvector(1,P);
  ForwardMultiSeq(phmm, T, O1, O2, alpha, logprobf, P, peakPos);
  BackwardMultiSeq(phmm, T, O1, O2, beta, logprobb, P, peakPos);
  ComputeGammaWithLog(phmm, T, alpha, beta, gamma);
  ComputeXiWithLog(phmm, T, O1, O2, alpha, beta, xi);
  logprobprev = 0.0;
  for (i = 1; i <= P; i ++) {
    logprobprev += logprobf[i];
  }
  
  do {  
    mean = 0.0;
    variance = 0.0;
    sumGamma = LOGZERO;
    
    for (i = 1; i <= (phmm->D); i++) { 
      for (k = 1; k <= P; k++) {
        tempSum = LOGZERO;
        for (t = peakPos[k]; t < peakPos[k+1]; t++) {
          if (gamma[t][i] != LOGZERO){
            tempSum = logCheckAdd(tempSum, gamma[t][i]);
          }
        }
        sumGamma = logCheckAdd(sumGamma, tempSum);
      }
    }
    
    for (i = 1; i <= (phmm->D); i++) { 
			//denominatorA = 0.0;  
      for (k = 1; k <= P; k++) {
        for (t = peakPos[k]; t < peakPos[k+1]; t++) {
				//denominatorA += gamma[t][i];
          mean += exp(gamma[t][i] - sumGamma) * O2[t];

        }
      }
    }

    phmm->mu[1] = mean;
    
    for (i = 1; i <= (phmm->D); i++) { 
      for (k = 1; k < P; k++) {
        for (t = peakPos[k]; t < peakPos[k+1]; t++) {
          variance += exp(gamma[t][i] - sumGamma) * pow((O2[t] - mean), 2.0);
          
        }
      }
    }

    phmm->sigma[1] = 0.00001 + 0.99999 * sqrt(variance);
    
    for (i = (phmm->D + 1); i <= phmm->N; i++) { 
      denominatorA = LOGZERO;
      mean = 0.0;
      variance = 0.0;
      sumGamma = LOGZERO;
      numNonZero = 0;
      for (k = 1; k <= P; k++) {
        tempSum = LOGZERO;
        for (t = peakPos[k]; t < peakPos[k+1] - 1; t++) {
          if (gamma[t][i] != LOGZERO){
            tempSum = logCheckAdd(tempSum, gamma[t][i]);
          }
        }
        denominatorA = logCheckAdd(denominatorA, tempSum);
        if (gamma[peakPos[k+1]-1][i] != LOGZERO){
          sumGamma = logCheckAdd(sumGamma, logCheckAdd(tempSum, gamma[peakPos[k+1]-1][i]));
        }
        else {
          sumGamma = logCheckAdd(sumGamma, tempSum);
        }
      }
      for (k = 1; k <= P; k++) {
        for (t = peakPos[k]; t < peakPos[k+1]; t++) {
				//denominatorA += gamma[t][i];
          mean += exp(gamma[t][i] - sumGamma) * O2[t];
        }
      }
      
      for (j = 1; j <= phmm->N; j++) {
        if (phmm->A[i][j] != 0.0){
          numNonZero ++;
        }
      }
      for (j = 1; j <= phmm->N; j++) {
        numeratorA = LOGZERO;
        for (k = 1; k <= P; k++) {
          tempNumeratorA = LOGZERO;
          for (t = peakPos[k]; t < peakPos[k+1] - 1; t++) {
            if (xi[t][i][j] != LOGZERO){
              tempNumeratorA = logCheckAdd(tempNumeratorA, xi[t][i][j]);
            }
            
          }
          numeratorA = logCheckAdd(numeratorA, tempNumeratorA);
          
        }
        if (phmm->A[i][j] != 0.0){
          phmm->A[i][j] = (0.999999 * exp(numeratorA - denominatorA) + 0.000001) / (0.999999 + 0.000001 * numNonZero);
          
        }
        
      }
        phmm->mu[i- (phmm->D) + 1] = mean;
        for (k = 1; k <= P; k++) {
          for (t = peakPos[k]; t < peakPos[k+1]; t++) {
            variance += exp(gamma[t][i] - sumGamma) * pow((O2[t] - mean), 2.0);
          }
        } 
        phmm->sigma[i- (phmm->D) + 1] = 0.00001 + 0.99999 * sqrt(variance);
    }
    ForwardMultiSeq(phmm, T, O1, O2, alpha, logprobf, P, peakPos);
    BackwardMultiSeq(phmm, T, O1, O2, beta, logprobb, P, peakPos);
    ComputeGammaWithLog(phmm, T, alpha, beta, gamma);
    ComputeXiWithLog(phmm, T, O1, O2, alpha, beta, xi);
    

    totalProb = 0.0;
    for (i = 1; i <= P; i ++) {
      totalProb += logprobf[i];
      fprintf(stdout, "%lf ", logprobf[i]);
    }
    fprintf(stdout, "\n %d %lf %lf", l, totalProb, logprobprev);
    delta = totalProb - logprobprev; 
    logprobprev = totalProb;
    l++;	
    fprintf(stdout, "\n %d %lf", l, delta);
    PrintHMM(stdout, phmm);
    
	}
  while (abs(delta) > DELTA); 
  
  *pniter = l;
  //*plogprobfinal = totalProb;
  FreeXi(xi, T, phmm->N);
}
*/

