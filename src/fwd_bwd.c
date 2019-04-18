
#include <stdio.h>
#include <math.h>
#include "nrutil.h"
#include "hmm.h"
#include "logmath.h"
#include <omp.h>


void Forward(HMM *phmm, int T, gsl_matrix * alpha_matrix, double *pprob, 
             int P, int *peakPos, gsl_matrix * emission_matrix)
{
  int thread_id, nloops;
  int i, j, k;   /* state indices */
  int t;      /* time index */
  int n;
  int start, end;
  int *TFlist, TF;
  double sum;     /* partial sum */
  double **alpha = dmatrix(phmm->N, T);
  TFlist = ivector(phmm->M);
  TF = 0;
  for (j = 0; j < phmm->M; j++){
    TFlist[j] = TF;
    //fprintf(stdout,"%d ", TFlist[j]);
    TF += phmm->D[j];
  }
  TF -= 1;
  //fprintf(stdout,"%d \n", TF);
/* 1. Initialization */
#pragma omp parallel num_threads(THREAD_NUM) \
  private(thread_id, nloops, start, end, sum, i, j, t, n) \
  //shared (P, peakPos, phmm, alpha,emission, pprob)
  {
    nloops = 0;
#pragma omp for
  for (k = 0; k < P; k++){
    //printf("1: %d ", i);
    //fflush(stdout);   
    ++nloops;
    for (i = 0; i < phmm->N; i++){
      
      if (phmm->pi[i] == 0.0){
        alpha[i][peakPos[k]-1] = -INFINITY;
      } 
      else {
        alpha[i][peakPos[k]-1] = log(phmm->pi[i]) + gsl_matrix_get (emission_matrix, i, peakPos[k]-1);
      }
   }
   
/* 2. Induction */
    start = peakPos[k];
    end = peakPos[k+1] - 1;
    for (t = start-1; t < end-1; t++) {
        for (j = 1; j <= TF; j++) {
          alpha[j][t+1] = alpha[j-1][t] + gsl_matrix_get(phmm->log_A_matrix, j-1, j) + gsl_matrix_get(emission_matrix, j, t+1);
        }
        for (n = 0; n < phmm->M; n++) {
          sum = -INFINITY;
          j = TFlist[n];
          for (i = TF + 1 ; i < phmm->N; i++) {
            if (alpha[i][t] != -INFINITY && gsl_matrix_get(phmm->log_A_matrix, i, j)!= -INFINITY) {
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
          for (n = 1; n < phmm->M; n++) {
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
          for (i = TF; i < phmm->N; i++) {
            if (alpha[i][t] != -INFINITY && gsl_matrix_get(phmm->log_A_matrix, i, j)!= -INFINITY) {
              if (sum != -INFINITY) {
                sum = logadd(sum, alpha[i][t] + gsl_matrix_get(phmm->log_A_matrix, i, j));  
              }
              else{
                sum = alpha[i][t] + gsl_matrix_get(phmm->log_A_matrix, i, j);
              }
            }
          }
          //alpha[j][t+1] = sum + (emission[j][t+1]);
          alpha[j][t+1] = sum + gsl_matrix_get(emission_matrix, j, t+1);
          //gsl_matrix_set (alpha_matrix, i, t, gsl_matrix_get (emission_matrix, j, t+1));
          //if (end==peakPos[2] - 1 && t >= end - 5) fprintf(stdout, "t,j: %d %d %lf %lf %lf\t", t, j, alpha[j][t+1], sum, (emission[j][t+1]));
        }
         
    }
   
/* 3. Termination */
    pprob[k] = -INFINITY;
    for (i = 0; i < phmm->N; i++) {
      if (alpha[i][end-1] != -INFINITY && alpha[i][end-1] 
          == alpha[i][end-1]) {
        if (pprob[k] != -INFINITY) {
          pprob[k] = logadd(pprob[k], alpha[i][end-1]);
          if (end==peakPos[2] - 1 && i > phmm->N - 6) fprintf(stdout, "a:%e %e\t", alpha[i][end-1], pprob[k]);
        }
        else {
          pprob[k] = alpha[i][end-1];
        }
      }
      
    }
    //fprintf(stdout, "%lf\t", pprob[k]);
    
  }
  thread_id = omp_get_thread_num();
  //printf("Thread %d performed %d iterations of the forward.\n",
    //thread_id, nloops );
  }
  for (i = 0; i < phmm->N; i ++){ 
    for (j = 0; j < T; j ++){
      gsl_matrix_set(alpha_matrix, i, j, alpha[i][j]);
    }
  }
  free_dmatrix(alpha, phmm->N, T);
  free_ivector(TFlist, phmm->M);
}


void Backward(HMM *phmm, int T, gsl_matrix * beta_matrix, int P, 
              int *peakPos, gsl_matrix * emission_matrix)
{
  int thread_id, nloops;
  int i, j, k, n;   /* state indices */
  int t;      /* time index */
  int start, end;
  int TF, *TFstartlist, *TFendlist;
  double sum, temp;
  double **beta = dmatrix(phmm->N, T);
  TFstartlist = ivector(phmm->M);
  TFendlist = ivector(phmm->M);
  TFstartlist[0] = 0;
  TFendlist[0] = phmm->D[0] - 1;
  for (j = 1; j < phmm->M; j++){
    TFstartlist[j] = TFstartlist[j - 1] + phmm->D[j - 1];
    TFendlist[j] = TFendlist[j - 1] + phmm->D[j];
  }
 
#pragma omp parallel num_threads(THREAD_NUM) \
  private(thread_id, nloops, start, end, sum, i, j, t, n) \
  //shared (P, peakPos, phmm, alpha,emission, pprob)
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
    //fprintf(stdout, "beta: %d  ", k);
  /* 2. Induction */
    start = peakPos[k];
    end = peakPos[k+1] - 1;
    for (t = end - 2; t >= start-1; t--) {

      for (i = 0; i < TFendlist[phmm->M - 1]; i++) {
        beta[i][t] = beta[i+1][t+1] + gsl_matrix_get(phmm->log_A_matrix, i, i+1) + 
                     gsl_matrix_get(emission_matrix, i+1, t+1);  
      }
      for (n = 0; n < phmm->M; n++) {
        i = TFendlist[n];
        beta[i][t] = -INFINITY;
        
        for (j = TFendlist[phmm->M - 1] + 1; j < phmm->N; j++) {
          beta[i][t] = logCheckAdd(beta[i][t], beta[j][t+1] + 
                       gsl_matrix_get(phmm->log_A_matrix, i, j)) + 
                       gsl_matrix_get (emission_matrix, j, t+1);  
        }
      }
      for (i = TFendlist[phmm->M - 1] + 1; i < phmm->N; i++) { 
        beta[i][t] = -INFINITY;
        for (n = 0; n < phmm->M; n++) {
          j = TFstartlist[n];
          beta[i][t] = logCheckAdd(beta[i][t], beta[j][t+1] + 
                       gsl_matrix_get(phmm->log_A_matrix, i, j)) + 
                       gsl_matrix_get (emission_matrix, j, t+1);  
        }
        for (j = TFendlist[phmm->M - 1] + 1; j < phmm->N; j++) { 
          beta[i][t] = logCheckAdd(beta[i][t], beta[j][t+1] + 
                       gsl_matrix_get(phmm->log_A_matrix, i, j)) + 
                       gsl_matrix_get (emission_matrix, j, t+1);  
        }
      }
    }
  }
  thread_id = omp_get_thread_num();
  
  }
  for (i = 0; i < phmm->N; i ++){ 
    for (j = 0; j < T; j ++){
      gsl_matrix_set (beta_matrix, i, j, beta[i][j]);
    }
  }
  free_dmatrix(beta, phmm->N, T);
  free_ivector(TFstartlist, phmm->M);
  free_ivector(TFendlist, phmm->M);
}


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
  /*
  if (phmm->M > 1){
    TFlist = ivector(phmm->M);
    TF = 0;
    for (j = 0; j < phmm->M; j++){
      TFlist[j] = TF;
      TF += phmm->D[j];
    }
    TF -= 1;
  }
  else{
    TFlist = ivector(2);
    TFlist[0] = 0;  
    TFlist[1] = phmm->D[0];
    TF = phmm->D[0] * 2 - 1;
  }
  */
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
  //fprintf(stdout,"%d \n", TF);
/* 1. Initialization */
#pragma omp parallel num_threads(THREAD_NUM) \
  private(thread_id, nloops, start, end, sum, i, j, t, n) \
  //shared (P, peakPos, phmm, alpha,emission, pprob)
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
          //if (alpha[j][t+1] != alpha[j][t+1]) fprintf(stdout, "t,j: %d %d %lf %lf %lf\t", t, j, alpha[j][t+1], sum,  gsl_matrix_get(emission_matrix, j, t+1));
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
          //if (end==peakPos[2] - 1 && (t <= start + 5 || t >= end - 5)) fprintf(stdout, "t,j: %d %d %lf %lf %lf\t", t, j, alpha[j][t+1], sum,  gsl_matrix_get(emission_matrix, j, t+1));
          
        }
        
    }
    
/* 3. Termination */
    pprob[k] = -INFINITY;
    for (i = 0; i < phmm->N; i++) {
      if (alpha[i][end-1] != -INFINITY && alpha[i][end-1] 
          == alpha[i][end-1]) {
        if (pprob[k] != -INFINITY) {
          pprob[k] = logadd(pprob[k], alpha[i][end-1]);
          if (end==peakPos[2] - 1 && i > phmm->N - 6) fprintf(stdout, "a:%e %e\t", alpha[i][end-1], pprob[k]);
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
  int t;      /* time index */
  int start, end;
  int *TFstartlist, *TFendlist;
  double sum, temp;
  /*
  if (phmm->M > 1){
    TFstartlist = ivector(phmm->M);
    TFendlist = ivector(phmm->M);
    TFstartlist[0] = 0;
    TFendlist[0] = phmm->D[0] - 1;
    for (j = 1; j < phmm->M; j++){
      TFstartlist[j] = TFstartlist[j - 1] + phmm->D[j - 1];
      TFendlist[j] = TFendlist[j - 1] + phmm->D[j]; 
    }
  }
  else{
    TFstartlist = ivector(2);
    TFendlist = ivector(2);
    TFstartlist[0] = 0;
    TFstartlist[1] = phmm->D[0];
    TFendlist[0] = phmm->D[0] - 1;
    TFendlist[1] = phmm->D[0] * 2 - 1;
  }
 */
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
  private(thread_id, nloops, start, end, sum, i, j, t, n) \
  //shared (P, peakPos, phmm, alpha,emission, pprob)
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