/*
 *  File: fwd_bwd.c
 *
 *  forward-backward algorithm
 *
 *  The HMM structure and some codes are borrowed and modified from Kanungo's
 *  original HMM program.
 *  Tapas Kanungo, "UMDHMM: Hidden Markov Model Toolkit," in "Extended Finite State Models of Language," A. Kornai (editor), Cambridge University Press, 1999. http://www.kanungo.com/software/software.html.
 *
 */

#include <stdio.h>
#include <math.h>
#include "nrutil.h"
#include "hmm.h"
#include "logmath.h"
#include <omp.h>

void Forward_P(HMM *phmm, int T, double **alpha, double *pprob, int P,
               int *peakPos, gsl_matrix * emission_matrix)
{
  int thread_id, nloops;
  int i, j, k;   /* state indices */
  int t;      /* time index */
  int n;
  int start, end;
  int *TFlist, TF;
  double sum;     /* partial sum */
  TFlist = ivector(phmm->M * (phmm->inactive+1));
    TF = 0;
    for (j = 0; j < phmm->M; j++){
      TFlist[j * (phmm->inactive+1)] = TF;
      TF += phmm->D[j];
      if (phmm->inactive == 1){
        TFlist[j * (phmm->inactive+1) + 1] = TF;
        TF += phmm->D[j];
      }
    }
    TF -= 1;
/* 1. Initialization */
#pragma omp parallel num_threads(THREAD_NUM) \
  private(thread_id, nloops, start, end, sum, i, j, t, n)
  {
    nloops = 0;
#pragma omp for
  for (k = 0; k < P; k++){
    ++nloops;
    for (i = 0; i < phmm->N; i++){
      if (phmm->pi[i] == 0.0){
        alpha[i][peakPos[k]-1] = -INFINITY;
      } 
      else {
        alpha[i][peakPos[k]-1] = log(phmm->pi[i]) + gsl_matrix_get(emission_matrix, i, peakPos[k]-1);
      }
   }
   
/* 2. Induction */
    start = peakPos[k];
    end = peakPos[k+1] - 1;
    for (t = start-1; t < end-1; t++) {
        for (j = 1; j <= TF; j++) {
          alpha[j][t+1] = alpha[j-1][t] + gsl_matrix_get(phmm->log_A_matrix, j-1, j) + gsl_matrix_get(emission_matrix, j, t+1);
        }
        for (n = 0; n < phmm->M * (phmm->inactive+1); n++) {
          sum = -INFINITY;
          j = TFlist[n];
          for (i = TF + 1 ; i < phmm->N; i++) {
            if (alpha[i][t] != -INFINITY &&  gsl_matrix_get(phmm->log_A_matrix, i, j)!= -INFINITY) {
              if (sum != -INFINITY) {
                sum = logadd(sum, alpha[i][t] + gsl_matrix_get(phmm->log_A_matrix, i, j));  
              }
              else{
                sum = alpha[i][t] + gsl_matrix_get(phmm->log_A_matrix, i, j);
              }
            }
          }
          alpha[j][t+1] = sum + gsl_matrix_get(emission_matrix, j, t+1);
        }
        for (j = TF + 1; j < phmm->N; j++) {
          sum = -INFINITY;
          for (n = 1; n < phmm->M * (phmm->inactive+1); n++) {
            i = TFlist[n] - 1;
            if (alpha[i][t] != -INFINITY && gsl_matrix_get(phmm->log_A_matrix, i, j)!= -INFINITY) {
              if (sum != -INFINITY) {
                sum = logadd(sum, alpha[i][t] + gsl_matrix_get(phmm->log_A_matrix, i, j));  
              }
              else{
                sum = alpha[i][t] + gsl_matrix_get(phmm->log_A_matrix, i, j);
              }
            }
          }
          for (i = 0; i < phmm->N; i++) {
            if (alpha[i][t] != -INFINITY && gsl_matrix_get(phmm->log_A_matrix, i, j)!= -INFINITY) {
              if (sum != -INFINITY) {
                sum = logadd(sum, alpha[i][t] + gsl_matrix_get(phmm->log_A_matrix, i, j));  
              }
              else{
                sum = alpha[i][t] + gsl_matrix_get(phmm->log_A_matrix, i, j);
              }
            }
          }
          alpha[j][t+1] = sum + (gsl_matrix_get(emission_matrix, j, t+1));
        }
    }
    
/* 3. Termination */
    pprob[k] = -INFINITY;
    for (i = 0; i < phmm->N; i++) {
      if (alpha[i][end-1] != -INFINITY && alpha[i][end-1] 
          == alpha[i][end-1]) {
        if (pprob[k] != -INFINITY) {
          pprob[k] = logadd(pprob[k], alpha[i][end-1]);
        }
        else {
          pprob[k] = alpha[i][end-1];
        }
      }
      
    }
    
  }
  thread_id = omp_get_thread_num();
  }
  free_ivector(TFlist, MAX(phmm->M, 2));
}


void Backward_P(HMM *phmm, int T, double **beta, int P, int *peakPos, 
                gsl_matrix * emission_matrix)
{
  int thread_id, nloops;
  int i, j, k, n;   /* state indices */
  int t;  /* time index */
  int start, end;
  int *TFstartlist, *TFendlist;
  double sum, temp;
  int TF = 0;
  TFstartlist = ivector(phmm->M * (phmm->inactive+1));
  TFendlist = ivector(phmm->M * (phmm->inactive+1));
  for (j = 0; j < phmm->M; j++){
      TFstartlist[j * (phmm->inactive+1)] = TF;
      TF += phmm->D[j];
      TFendlist[j * (phmm->inactive+1)] = TF - 1;
      if (phmm->inactive == 1){
        TFstartlist[j * (phmm->inactive+1) + 1] = TF;
        TF += phmm->D[j];
        TFendlist[j * (phmm->inactive+1) + 1] = TF - 1;
      } 
    }
  TF -= 1;
#pragma omp parallel num_threads(THREAD_NUM) \
  private(thread_id, nloops, start, end, sum, i, j, t, n)
  {
    nloops = 0;
#pragma omp for
  /* 1. Initialization */
  for (k = 0; k < P; k++){
    ++nloops;
    for (i = 0; i < phmm->N; i++){
      beta[i][peakPos[k+1]-2] = -INFINITY;
    }            
    beta[phmm->N-2][peakPos[k+1]-2] = 0.0;
  /* 2. Induction */
    start = peakPos[k];
    end = peakPos[k+1] - 1;
    for (t = end - 2; t >= start-1; t--) {
    

      for (i = 0; i < TF; i++) {
        beta[i][t] = beta[i+1][t+1] + gsl_matrix_get(phmm->log_A_matrix, i, i+1) + gsl_matrix_get(emission_matrix, i+1, t+1);  
      }
      for (n = 0; n < phmm->M * (phmm->inactive+1); n++) {
        i = TFendlist[n];
        beta[i][t] = -INFINITY;
        
        for (j = TF + 1; j < phmm->N; j++) {
          beta[i][t] = logCheckAdd(beta[i][t], beta[j][t+1] + gsl_matrix_get(phmm->log_A_matrix, i, j) + gsl_matrix_get(emission_matrix, j, t+1));  
        }
      }
      for (i = TF + 1; i < phmm->N; i++) {
        beta[i][t] = -INFINITY;
        for (n = 0; n < phmm->M * (phmm->inactive+1); n++) {
          j = TFstartlist[n];
          beta[i][t] = logCheckAdd(beta[i][t], beta[j][t+1] + gsl_matrix_get(phmm->log_A_matrix, i, j) + gsl_matrix_get(emission_matrix, j, t+1));  
        }
        for (j = TF + 1; j < phmm->N; j++) {
          beta[i][t] = logCheckAdd(beta[i][t], beta[j][t+1] + gsl_matrix_get(phmm->log_A_matrix, i, j) + gsl_matrix_get(emission_matrix, j, t+1));  
        }
      }
    }
    
  }
  thread_id = omp_get_thread_num();
  }
  free_ivector(TFstartlist, phmm->M * (phmm->inactive+1));
  free_ivector(TFendlist, phmm->M * (phmm->inactive+1));
}
