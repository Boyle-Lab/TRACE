/*
 * File: BaumWelch.c
 *
 * Baum-Welch algorithm for HMM training.
 *
 * The HMM structure and some codes are borrowed and modified from Kanungo's
 * original HMM program.
 * Tapas Kanungo, "UMDHMM: Hidden Markov Model Toolkit," in "Extended Finite State Models of Language," A. Kornai (editor), Cambridge University Press, 1999. http://www.kanungo.com/software/software.html.
 *
 */

#include <stdio.h>
#include "nrutil.h"
#include "hmm.h"
#include <math.h>
#include "logmath.h"
#include <time.h>
#include <omp.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#define DELTA 0.001

void BaumWelch(HMM *phmm, int T, gsl_matrix * obs_matrix, int *pniter, int P,
               int *peakPos, double *logprobf, double **alpha, double **beta, 
               double **gamma, gsl_matrix * emission_matrix)
{
  int i, j, k, n, m, x, TF;
  int start, end;
  int t, l = 0;
  int thread_id, nloops;
  int *stateList; /* list of states that have none-one transition probability */
  int *motifList; /* list of index of motifs */
  //int blocklimit; /*for loop unrolling */
  
  /*the following are temp values used when calcaulating variables in BW*/
  int numNonZero;
  double *logprobb, *freq, *plogprobinit, plogprobfinal;
  double numeratorA, denominatorA;
  double *mean;
  double cov; /*covariance*/
  double delta, logprobprev, totalProb;
  double deltaprev = 10e-70;
  double sumGamma, tempSum, tempNumeratorA, xi;
  double **tmp_A;
  double tmp;

  stateList = ivector(phmm->M * (phmm->inactive+1) + phmm->extraState);
  TF = 0;
  /* stateList contains the first state of each motif and all additional states
   besides binding sites*/
  for (j = 0; j < phmm->M; j++){
    stateList[j * (phmm->inactive+1)] = TF;
    TF += phmm->D[j];
    if (phmm->inactive == 1){
      stateList[j * (phmm->inactive+1) + 1] = TF;
      TF += phmm->D[j];
    }
  }
  TF -= 1;
  for (j = phmm->M * (phmm->inactive+1); j < phmm->M * (phmm->inactive+1) + phmm->extraState; j++){
    stateList[j] = TF + j - phmm->M * (phmm->inactive+1) + 1;
  }
  TF = 0;
  motifList = ivector(phmm->N);
  /* motifList contians the index of motif that each state is part of */
  if (phmm->inactive == 1) {
    for (j = 0; j < phmm->M; j++) {
      for (i = TF; i < TF + phmm->D[j]; i++) {
        motifList[i] = 2*j;
      }
      TF += phmm->D[j];
      for (i = TF; i < TF + phmm->D[j]; i++) {
        motifList[i] = 2*j+1;
      }
      TF += phmm->D[j];
    }
  }
  /* for all other states, motifList will contain their actual state number */
  for (i = TF; i < TF + phmm->extraState; i++) {
    motifList[i] = i;
  }
  TF -= 1;
    
  /* choose from independent model, covariant model and covariant model with state dropping*/
  if (phmm->model == 0) EmissionMatrix(phmm, obs_matrix, P, peakPos, emission_matrix, T);
  if (phmm->model == 1) EmissionMatrix_mv(phmm, obs_matrix, P, peakPos, emission_matrix, T);
  if (phmm->model == 2) EmissionMatrix_mv_reduce(phmm, obs_matrix, P, peakPos, emission_matrix, T);
  /* forward-backward algorithm, calculate alpha, beta and gamma*/
  double *xi_sum = dvector(T);
  Forward_P(phmm, T, alpha, logprobf, P, peakPos, emission_matrix);
  Backward_P(phmm, T, beta, P, peakPos, emission_matrix);
  ComputeGamma_P(phmm, T, alpha, beta, gamma);
  ComputeXi_sum_P(phmm, alpha, beta, xi_sum, emission_matrix, T);
  logprobprev = 0.0;
  /* total log probability is the sum of log probability of each peak */
  for (i = 0; i < P; i ++) {
    if (logprobf[i] != -INFINITY){
      logprobprev += logprobf[i];
    }
  }
  /* some matrices and vectors that will be used in reestimating parameter */
  gsl_matrix * prob_matrix;
  gsl_matrix * prob_trans;
  gsl_matrix * post_obs;
  gsl_matrix * obs_obs;
  gsl_matrix * tmp_matrix, *tmp_matrix_2;
  gsl_vector * prob_sum, *post_sum, *prob_vector;
  gsl_vector * prob_log_sum;
  gsl_vector * mul_vector;
  gsl_vector * tmp_vector, *tmp_vector_2;
  gsl_matrix * obs_trans;
  mul_vector = gsl_vector_alloc(phmm->N);
  gsl_vector_set_all(mul_vector, -INFINITY);

  do {
  /* store transition matrix in a temporary matrix */
  tmp_A = dmatrix(phmm->N, phmm->N);
  for (i = 0; i < phmm->N; i ++){
    for (j = 0; j < phmm->N; j ++){
      tmp_A[i][j] = (gsl_matrix_get(phmm->log_A_matrix, i, j));
    }
  }
  prob_matrix = gsl_matrix_alloc(phmm->N, T);
/* reestimate transition matrix and symbol prob in each state */
#pragma omp parallel num_threads(THREAD_NUM)\
  private(thread_id, nloops, n, x, k, j, m, xi, numeratorA, denominatorA, \
  mean, cov, tempSum, sumGamma, t, start, end, \
  numNonZero,tempNumeratorA)
  {
  #pragma omp for
    for (i = 0; i < phmm->N; i++) {
      denominatorA = -INFINITY; /* denominator when calculating aij,
                                   which is sum of gamma(i) */
      sumGamma = -INFINITY; /* sum of gamma used to calculate mean and variance
                               but not aij*/
      numNonZero = 0;
      for (k = 0; k < P; k++) {
        tempSum = -INFINITY;
        start = peakPos[k];
        end = peakPos[k+1] - 1;
        for (t = start-1; t < end-1; t++) {
          tempSum = logCheckAdd(tempSum, gamma[i][t]);
        }
        denominatorA = logCheckAdd(denominatorA, (tempSum)); //XXX - logprobf[k]
        /* the last base at each peak is used when calculating
           mean and variance but not aij */
        sumGamma = logCheckAdd(sumGamma, (logCheckAdd(tempSum, gamma[i][end-1])));
      }
      for (k = 0; k < P; k++) {
        start = peakPos[k];
        end = peakPos[k+1] - 1;
        for (t = start-1; t < end; t++) {
          gsl_matrix_set(prob_matrix, i, t, exp(gamma[i][t] - sumGamma)); 
        }
      }
      /* only the aij in which state i outside binding sites need to be reestimated */
      if (i > TF){
        for (n = 0; n < phmm->M * (phmm->inactive+1) + phmm->extraState; n++) {
          j = stateList[n];
          numeratorA = -INFINITY; /* numerator when calculate aij,
                                     which is the sum of xi(i,j) */
          if (tmp_A[i][j] != -INFINITY){
            for (k = 0; k < P; k++) {
              tempNumeratorA = -INFINITY;
              start = peakPos[k];
              end = peakPos[k+1] - 1;
              for (t = start-1; t < end-1; t++) {
                xi = alpha[i][t] + beta[j][t+1] + (tmp_A[i][j]) +
                     gsl_matrix_get(emission_matrix, j, t+1) - xi_sum[t];
                tempNumeratorA = logCheckAdd(tempNumeratorA, xi);
              }
              numeratorA = logCheckAdd(numeratorA, (tempNumeratorA)); //XXX - logprobf[k]
            }
            gsl_matrix_set(phmm->log_A_matrix, i , j, numeratorA - denominatorA);
            if (numeratorA == -INFINITY && denominatorA == -INFINITY) {
              gsl_matrix_set(phmm->log_A_matrix, i, j, -INFINITY);
            }
          }
        }
      }
    }
    thread_id = omp_get_thread_num();
  }
  free_dmatrix(tmp_A, phmm->N, phmm->N);

#pragma omp parallel num_threads(THREAD_NUM)\
  private(tmp_vector, prob_vector, tmp)
  {
#pragma omp for
    for (i = 0; i < phmm->N; i++){               
      prob_vector = gsl_vector_alloc(T);
      tmp_vector = gsl_vector_alloc(phmm->K);    
      gsl_matrix_get_row(prob_vector, prob_matrix, i);
      gsl_blas_dgemv(CblasNoTrans, 1.0, obs_matrix, prob_vector, 0.0, tmp_vector);   
      gsl_matrix_set_col(phmm->mean_matrix, i, tmp_vector);
      gsl_vector_free(prob_vector);
      prob_vector = gsl_vector_alloc(phmm->K);
      gsl_vector_set_all(prob_vector, 1.0);
      gsl_blas_ddot(tmp_vector, prob_vector, &tmp);
      gsl_vector_free(tmp_vector);
      gsl_vector_free(prob_vector);
    }
  }

  /* reestimate the variance or covariance matrix */
  prob_sum = gsl_vector_alloc(phmm->N);
  if (phmm->model == 0) UpdateVariance_2(phmm, obs_matrix, prob_sum, prob_matrix, T, TF);
  if (phmm->model == 1 || phmm->model == 2) UpdateCovariance_2(phmm, obs_matrix, prob_sum, prob_matrix, T, TF);
  gsl_matrix_free(prob_matrix);
  gsl_vector_free(prob_sum);
  /* calculate new emission probs based in reestimated parameters */
  if (phmm->model == 0) EmissionMatrix(phmm, obs_matrix, P, peakPos, emission_matrix, T);
  if (phmm->model == 1) EmissionMatrix_mv(phmm, obs_matrix, P, peakPos, emission_matrix, T);
  if (phmm->model == 2) EmissionMatrix_mv_reduce(phmm, obs_matrix, P, peakPos, emission_matrix, T);
  /* calculate new alpha, beta, gamma and xi for next iteration of forward-backward algorithm*/
  Forward_P(phmm, T, alpha, logprobf, P, peakPos, emission_matrix);
  Backward_P(phmm, T, beta, P, peakPos, emission_matrix);
  ComputeGamma_P(phmm, T, alpha, beta, gamma);
  ComputeXi_sum_P(phmm, alpha, beta, xi_sum, emission_matrix, T);
    
  /* compute difference between log probability of two iterations */
  totalProb = 0.0;
  for (i = 0; i < P; i ++) {
    if (logprobf[i] != log(0.0)){
      totalProb += logprobf[i];
    }
  }
  delta = totalProb - logprobprev;
    //fprintf(stdout, "\n logprobf: %d %lf %lf", l, totalProb, logprobprev);
  logprobprev = totalProb;
  l++;
  /* print out the teh number of current iteration and prob change */
  if (l == 1||ceilf((l*10.0)/MAXITERATION)==(l*10.0)/MAXITERATION) {
    fprintf(stdout, "iteration: %d delta: %lf \n", l, delta);
  }
      
  }
  while ((fabs(delta) > DELTA) && l <= MAXITERATION); /* if log probability does not 
                                                         change much, exit */
                                                      /* defult max iterations is 200 */
  gsl_vector_free(mul_vector);
  *pniter = l;
}

/*
void UpdateVariance(HMM *phmm, gsl_matrix * obs_matrix, gsl_matrix * post_obs,
                    gsl_vector * prob_sum, gsl_matrix *prob_matrix, int T, int TF)
{
  int i, j;
  int thread_id, nloops;
  gsl_matrix * obs_trans, * tmp_matrix, * tmp_matrix_2, * obs_obs;
  obs_trans = gsl_matrix_alloc(T, phmm->K);
  gsl_matrix_transpose_memcpy(obs_trans, obs_matrix); 
  obs_obs = gsl_matrix_alloc(phmm->K, phmm->N);
  tmp_matrix_2 = gsl_matrix_alloc(phmm->N, phmm->K);
  tmp_matrix = gsl_matrix_alloc(T, phmm->K);
  gsl_matrix_memcpy(tmp_matrix, obs_trans);
  gsl_matrix_mul_elements(tmp_matrix, obs_trans);
  gsl_matrix_free(obs_trans);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
                 1.0, prob_matrix, tmp_matrix,
                 0.0, tmp_matrix_2);
  gsl_matrix_transpose_memcpy(obs_obs, tmp_matrix_2);
  gsl_matrix_free(tmp_matrix);
  gsl_matrix_free(tmp_matrix_2);
  
  tmp_matrix = gsl_matrix_alloc(phmm->K, phmm->N);
  gsl_matrix_memcpy(tmp_matrix, phmm->mean_matrix);
  gsl_matrix_scale(tmp_matrix, 2.0);
  tmp_matrix_2 = gsl_matrix_alloc(phmm->K, phmm->N);
  gsl_matrix_transpose_memcpy(tmp_matrix_2, post_obs);
  gsl_matrix_mul_elements(tmp_matrix, tmp_matrix_2);    
  gsl_matrix_sub(obs_obs, tmp_matrix);
  gsl_matrix_memcpy(tmp_matrix, phmm->mean_matrix);
  gsl_matrix_mul_elements(tmp_matrix, phmm->mean_matrix);  
  
  for (i = 0; i < phmm->K; i++){
     gsl_matrix_set_row(tmp_matrix_2, i, prob_sum);
  }
  gsl_matrix_mul_elements(tmp_matrix, tmp_matrix_2); 
  gsl_matrix_add(obs_obs, tmp_matrix);
  gsl_matrix_div_elements(obs_obs, tmp_matrix_2);
  for (i = 0; i < phmm->K; i++){
    for (j = 0; j < phmm->N; j++){
      gsl_matrix_set(phmm->var_matrix, i, j, sqrt(gsl_matrix_get(obs_obs, i, j)));
    }
  }
  gsl_matrix_free(tmp_matrix);
  gsl_matrix_free(tmp_matrix_2);
  gsl_matrix_free(obs_obs);
}


void UpdateCovariance(HMM *phmm, gsl_matrix * obs_matrix, gsl_matrix * post_obs,
                      gsl_vector * prob_sum, gsl_matrix *prob_matrix, int T, int TF)
{
  int i, j;
  int thread_id, nloops;
  gsl_matrix * obs_trans, * tmp_matrix, * obs_obs;
  gsl_vector * tmp_vector, * tmp_vector_2;
  obs_trans = gsl_matrix_alloc(T, phmm->K);
  gsl_matrix_transpose_memcpy(obs_trans, obs_matrix); 
  
#pragma omp parallel num_threads(THREAD_NUM) \
  private(thread_id, nloops, i, j, tmp_vector, tmp_vector_2, tmp_matrix, obs_obs) 
  {
    nloops = 0;
#pragma omp for   
    for (i = 0; i < phmm->N; i++){
      
      tmp_vector = gsl_vector_alloc(T);
      tmp_vector_2 = gsl_vector_alloc(phmm->K);
      tmp_matrix = gsl_matrix_alloc(T, phmm->K);
      obs_obs = gsl_matrix_alloc(phmm->K, phmm->K);
      gsl_matrix_get_row(tmp_vector, prob_matrix, i);
      for (j = 0; j < phmm->K; j ++){
        gsl_matrix_set_col(tmp_matrix, j, tmp_vector);
      }
      gsl_matrix_mul_elements(tmp_matrix, obs_trans);
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
                     1.0, obs_matrix, tmp_matrix,
                     0.0, obs_obs); //stats['obs*obs.T']
      gsl_vector_free(tmp_vector);
      gsl_matrix_free(tmp_matrix);
      tmp_vector = gsl_vector_alloc(phmm->K);
      tmp_matrix = gsl_matrix_alloc(phmm->K, phmm->K);
      gsl_matrix_get_row(tmp_vector_2, post_obs, i);
      for (j = 0; j < phmm->K; j ++){
        gsl_matrix_get_col(tmp_vector, phmm->mean_matrix, i);
        gsl_vector_scale(tmp_vector, gsl_vector_get(tmp_vector_2, j));
        gsl_matrix_set_row(tmp_matrix, j, tmp_vector); //obsmean 
      }
      gsl_matrix_sub(obs_obs, tmp_matrix);
      
      gsl_matrix_get_col(tmp_vector, phmm->mean_matrix, i);
      for (j = 0; j < phmm->K; j ++){
        gsl_matrix_get_row(tmp_vector_2, post_obs, i);
        gsl_vector_scale(tmp_vector_2, gsl_vector_get(tmp_vector, j));
        gsl_matrix_set_row(tmp_matrix, j, tmp_vector_2); //obsmean 
      }
      gsl_matrix_sub(obs_obs, tmp_matrix);
      
      gsl_matrix_get_col(tmp_vector_2, phmm->mean_matrix, i);
      for (j = 0; j < phmm->K; j ++){
        gsl_matrix_get_col(tmp_vector, phmm->mean_matrix, i);
        gsl_vector_scale(tmp_vector, gsl_vector_get(tmp_vector_2, j));
        gsl_matrix_set_row(tmp_matrix, j, tmp_vector); //obsmean 
      }
      gsl_matrix_scale(tmp_matrix, gsl_vector_get(prob_sum, i));
      gsl_matrix_add(obs_obs, tmp_matrix);
      gsl_matrix_scale(obs_obs, 1.0/gsl_vector_get(prob_sum, i));
      gsl_matrix_memcpy(phmm->cov_matrix[i], obs_obs);
      gsl_vector_free(tmp_vector);
      gsl_vector_free(tmp_vector_2);
      gsl_matrix_free(tmp_matrix);
      gsl_matrix_free(obs_obs);
    }
    thread_id = omp_get_thread_num();
  }
}
*/

void UpdateVariance_2(HMM *phmm, gsl_matrix * obs_matrix, 
                      gsl_vector * prob_sum, gsl_matrix *prob_matrix, 
                      int T, int TF)
{
  int i, j;
  int thread_id, nloops;
  double tmp;
  gsl_matrix * obs_trans, * tmp_matrix;
  gsl_vector * tmp_vector, * tmp_vector_2;
#pragma omp parallel num_threads(THREAD_NUM) \
  private(thread_id, nloops, i, j, tmp_vector, tmp_vector_2, tmp_matrix, obs_trans, tmp) 
  {
    nloops = 0;
#pragma omp for   
    for (i = 0; i < phmm->N; i++){
      tmp_vector = gsl_vector_alloc(T);
      tmp_vector_2 = gsl_vector_alloc(T);
      for (j = 0; j < phmm->K; j ++){
        gsl_matrix_get_row(tmp_vector_2, prob_matrix, i);
        gsl_matrix_get_row(tmp_vector, obs_matrix, j);
        tmp = gsl_matrix_get(phmm->mean_matrix, j, i);
        gsl_vector_add_constant(tmp_vector, tmp * -1.0);
        gsl_vector_mul(tmp_vector_2, tmp_vector);
        gsl_blas_ddot (tmp_vector_2, tmp_vector, &tmp);
        gsl_matrix_set(phmm->var_matrix, j, i, sqrt(tmp));
      }
      gsl_vector_free(tmp_vector);
      gsl_vector_free(tmp_vector_2);
    }
    thread_id = omp_get_thread_num();
  }
}


void UpdateCovariance_2(HMM *phmm, gsl_matrix * obs_matrix, 
                        gsl_vector * prob_sum, gsl_matrix *prob_matrix, 
                        int T, int TF)
{
  int i, j;
  int thread_id, nloops;
  double tmp;
  gsl_matrix * obs_trans, * tmp_matrix;
  gsl_vector * tmp_vector, * tmp_vector_2;
#pragma omp parallel num_threads(THREAD_NUM) \
  private(thread_id, nloops, i, j, tmp_vector, tmp_vector_2, tmp_matrix, obs_trans, tmp) 
  {
    nloops = 0;
#pragma omp for   
    for (i = 0; i < phmm->N; i++){
      tmp_vector = gsl_vector_alloc(T);
      tmp_vector_2 = gsl_vector_alloc(T);
      tmp_matrix = gsl_matrix_alloc(phmm->K, T);
      for (j = 0; j < phmm->K; j ++){
        gsl_matrix_get_row(tmp_vector, obs_matrix, j);
        tmp = gsl_matrix_get(phmm->mean_matrix, j, i);
        gsl_vector_add_constant(tmp_vector, tmp * -1.0);
        gsl_matrix_set_row(tmp_matrix, j, tmp_vector);
      }
      obs_trans = gsl_matrix_alloc(T, phmm->K);
      gsl_matrix_transpose_memcpy(obs_trans, tmp_matrix); 
      
      gsl_matrix_get_row(tmp_vector, prob_matrix, i);
      for (j = 0; j < phmm->K; j ++){
        gsl_matrix_get_row(tmp_vector_2, tmp_matrix, j);
        gsl_vector_mul(tmp_vector_2, tmp_vector);
        gsl_matrix_set_row(tmp_matrix, j, tmp_vector_2);
      }
      gsl_vector_free(tmp_vector);
      gsl_vector_free(tmp_vector_2);
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
                     1.0, tmp_matrix, obs_trans,
                     0.0, phmm->cov_matrix[i]);
      
      gsl_matrix_free(tmp_matrix);
      gsl_matrix_free(obs_trans);
    }
    thread_id = omp_get_thread_num();
  }
}


void ComputeGamma(HMM *phmm, int T, gsl_matrix * alpha_matrix, 
                  gsl_matrix * beta_matrix, gsl_matrix * gamma_matrix)
{
  int thread_id, nloops;
  int i, j;
  int	t;
  double denominator;
  double **gamma = dmatrix(phmm->N, T);
#pragma omp parallel num_threads(THREAD_NUM)\
  private(thread_id, nloops, j, denominator, i) 
  {
    nloops = 0;
#pragma omp for
  for (t = 0; t < T; t++) { 
    ++nloops;
    denominator = -INFINITY;
    for (j = 0; j < phmm->N; j++) {
      gamma[j][t] = gsl_matrix_get(alpha_matrix, j, t) + 
                    gsl_matrix_get(beta_matrix, j, t);
      denominator = logCheckAdd(denominator, gamma[j][t]);
    }
    if (denominator == -INFINITY){
      fprintf(stdout, "ERROR: invalid gamma value at position:%d ",t);
      //fprintf(stdout, "ab: %d   %lf %lf %lf\t", t,
              //gsl_matrix_get(alpha_matrix, j, t),
              //gsl_matrix_get(beta_matrix, j, t), denominator);
      exit(1);
    }
    for (i = 0; i < phmm->N; i++) {
      gamma[i][t] -= denominator;
    }
    
  }
  thread_id = omp_get_thread_num();
  }
  for (i = 0; i < phmm->N; i ++){ 
    for (j = 0; j < T; j ++){
      gsl_matrix_set(gamma_matrix, i, j, gamma[i][j]);
    }
  }
  free_dmatrix(gamma, phmm->N, T);
}

void ComputeXi_sum(HMM* phmm, gsl_matrix * alpha_matrix, 
                   gsl_matrix * beta_matrix, gsl_vector * xi_sum_vector, 
                   gsl_matrix * emission_matrix, int T)
{
  int thread_id, nloops;
  int i, j, k, n;
  int t, index;
  int TF, *TFstartlist, *TFendlist;
  double *xi_sum = dvector(T);
  TFstartlist = ivector(phmm->M);
  TFendlist = ivector(phmm->M);
  TFstartlist[0] = 0;
  TFendlist[0] = phmm->D[0] - 1;
  for (j = 1; j < phmm->M; j++){
    TFstartlist[j] = TFstartlist[j - 1] + phmm->D[j - 1];
    TFendlist[j] = TFendlist[j - 1] + phmm->D[j];
  }
  double sum;
  double xi; 
  double *sum_list = dvector(T);

#pragma omp parallel num_threads(THREAD_NUM) \
  private(thread_id, nloops, i, j, sum, index, xi, n) 
  {
    nloops = 0;
#pragma omp for
  for (t = 0; t < (T - 1); t++) {
    ++nloops;
    sum = -INFINITY;
    for (n = 0; n < phmm->M; n++) {
      for (i = TFstartlist[n]; i < TFendlist[n]; i++) {
        j = i + 1;
        xi = (gsl_matrix_get(alpha_matrix, i, t) + 
              gsl_matrix_get(beta_matrix, j, t+1)) + 
              gsl_matrix_get(phmm->log_A_matrix, i, j)+
              gsl_matrix_get(emission_matrix, j, t+1);
        sum = logCheckAdd(sum, xi);
      }
      i = TFendlist[n];
      j = TFendlist[phmm->M - 1] + 1;
      xi = (gsl_matrix_get(alpha_matrix, i, t) + 
              gsl_matrix_get(beta_matrix, j, t+1)) + 
              gsl_matrix_get(phmm->log_A_matrix, i, j)+
              gsl_matrix_get(emission_matrix, j, t+1);
      sum = logCheckAdd(sum, xi);
    }
    for (i = TFendlist[phmm->M - 1] + 1; i < phmm->N; i++) {
      for (n = 0; n < phmm->M; n++) {
        j = TFstartlist[n];
        xi = (gsl_matrix_get(alpha_matrix, i, t) + 
              gsl_matrix_get(beta_matrix, j, t+1)) + 
              gsl_matrix_get(phmm->log_A_matrix, i, j)+
              gsl_matrix_get(emission_matrix, j, t+1);
        sum = logCheckAdd(sum, xi);
      }
      for (j = TFendlist[phmm->M - 1] + 1; j < phmm->N; j++) {
        xi = (gsl_matrix_get(alpha_matrix, i, t) + 
              gsl_matrix_get(beta_matrix, j, t+1)) + 
              gsl_matrix_get(phmm->log_A_matrix, i, j)+
              gsl_matrix_get(emission_matrix, j, t+1);
        sum = logCheckAdd(sum, xi);
      }
    }
    sum_list[t] = sum;
    
  }
  thread_id = omp_get_thread_num();
  }
  for (i = 0; i < T; i ++){
    gsl_vector_set(xi_sum_vector, i, sum_list[i]);
  }
  free_dvector(sum_list, T);
}


void ComputeGamma_P(HMM *phmm, int T, double **alpha, double **beta, 
                    double **gamma)
{
  int thread_id, nloops;
  int i, j;
  int	t;
  double denominator;
#pragma omp parallel num_threads(THREAD_NUM)\
  private(thread_id, nloops, j, denominator, i)
  {
    nloops = 0;
#pragma omp for
  for (t = 0; t < T; t++) {
    ++nloops;
    denominator = -INFINITY;
    for (j = 0; j < phmm->N; j++) {
      gamma[j][t] = alpha[j][t]+beta[j][t];
          denominator = logCheckAdd(denominator, gamma[j][t]);
    }
    if (denominator == -INFINITY){
      fprintf(stdout, "ERROR: -inf gamma t:%d ",t);
      //fprintf(stdout, "ab: %d   %lf %lf %lf\t", t, alpha[1][t], beta[1][t],denominator);
      exit(1);
    }
    if (denominator != denominator){
      fprintf(stdout, "ERROR: na gamma t:%d ",t);
      //fprintf(stdout, "ab: %d   %lf %lf %lf\t", t, alpha[1][t], beta[1][t],denominator);
      exit(1);
    }
    for (i = 0; i < phmm->N; i++) {
      gamma[i][t] -= denominator;
      if (gamma[i][t] != gamma[i][t]){
        fprintf(stdout, "ERROR: gamma NA t:%d ",t);
        //fprintf(stdout, "ab: %d   %lf %lf %lf\t", t, alpha[i][t], beta[i][t], gamma[i][t]);
        exit(1);
      }
    }
  }
  thread_id = omp_get_thread_num();
  }
}

void ComputeXi_sum_P(HMM* phmm, double **alpha, double **beta, double *xi_sum, 
                     gsl_matrix * emission_matrix, int T)
{
  int thread_id, nloops;
  int i, j, k, n;
  int t, index;
  int TF, *TFstartlist, *TFendlist;
  TF = 0;
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
    
  double sum;
  double xi; 
#pragma omp parallel num_threads(THREAD_NUM) \
  private(thread_id, nloops, i, j, sum, index, xi, n)
  {
    nloops = 0;
#pragma omp for
  for (t = 0; t < (T - 1); t++) {
    ++nloops;
    sum = -INFINITY;
    for (n = 0; n < phmm->M * (phmm->inactive+1); n++) {
      for (i = TFstartlist[n]; i < TFendlist[n]; i++) {
        j = i + 1;
        xi = (alpha[i][t] + beta[j][t+1]) + 
             gsl_matrix_get(phmm->log_A_matrix, i, j) + 
             gsl_matrix_get(emission_matrix, j, t+1);
        sum = logCheckAdd(sum, xi);
      }
      i = TFendlist[n];
      j = TFendlist[phmm->M - 1] + 1;
      if (phmm->extraState > (phmm->M * (phmm->inactive+1)*phmm->lPeak)) { //TODO:change how to determine j
        j = TFendlist[phmm->M * (phmm->inactive+1) - 1] + 1 + phmm->lPeak * n;
      }
      else if (phmm->extraState >= (2*phmm->lPeak+4)) {
        if(n % 2 == 0) j = TFendlist[phmm->M * (phmm->inactive+1) - 1] + 1;
        else j = TFendlist[phmm->M * (phmm->inactive+1) - 1] + 1 + phmm->lPeak;
      }
      xi = (alpha[i][t] + beta[j][t+1]) + 
            gsl_matrix_get(phmm->log_A_matrix, i, j) + 
            gsl_matrix_get(emission_matrix, j, t+1);
      sum = logCheckAdd(sum, xi);
    }
    for (i = TFendlist[phmm->M * (phmm->inactive+1) - 1] + 1; i < phmm->N; i++) {
      for (n = 0; n < phmm->M * (phmm->inactive+1); n++) {
        j = TFstartlist[n];
        xi = (alpha[i][t] + beta[j][t+1]) + 
              gsl_matrix_get(phmm->log_A_matrix, i, j) + 
              gsl_matrix_get(emission_matrix, j, t+1);
        sum = logCheckAdd(sum, xi);
      }
      for (j = TFendlist[phmm->M * (phmm->inactive+1) - 1] + 1; j < phmm->N; j++) {
        xi = (alpha[i][t] + beta[j][t+1]) + 
              gsl_matrix_get(phmm->log_A_matrix, i, j) + 
              gsl_matrix_get(emission_matrix, j, t+1);
        sum = logCheckAdd(sum, xi);
      }
    }
    xi_sum[t] = sum;
    
  }
  thread_id = omp_get_thread_num();
  }
}

