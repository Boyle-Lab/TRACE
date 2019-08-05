/*
 *  File: hmmutils.c
 *
 *  Some functions to handel model file.
 *
 *  The HMM structure and some codes are borrowed and modified from Kanungo's
 *  original HMM program.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nrutil.h"
#include "hmm.h"

/* Read the M parameter in model file */
void ReadM(FILE *fp, HMM *phmm)
{
  fscanf(fp, "M= %d\n", &(phmm->M)); 
}

/* Read parameters in a initial model file */
void ReadInitHMM(FILE *fp, HMM *phmm)
{
  int i, j, k, n, totalStates = 0;
  double sum;
  double tmp;

  fscanf(fp, "N= %d\n", &(phmm->N));
  /* Read transition matrix */
  fscanf(fp, "A:\n");
  phmm->log_A_matrix = gsl_matrix_alloc(phmm->N, phmm->N);
  for (i = 0; i < phmm->N; i++) {
    for (j = 0; j < phmm->N; j++) {
      fscanf(fp, "%lf", &tmp);
      gsl_matrix_set(phmm->log_A_matrix, i, j, log(tmp));
    }
    fscanf(fp,"\n");
  }
  /* Read PWMs */
  phmm->D = (int *) ivector(phmm->M);
  phmm->pwm = (double ***)malloc(phmm->M*sizeof(double**));
  for (i = 0; i < phmm->M; i++) {
    fscanf(fp, "PWM: n=%d\n", &(phmm->D[i]));
    totalStates += phmm->D[i];
    phmm->pwm[i] = (double **) dmatrix(phmm->D[i], 4);
    for (j = 0; j < phmm->D[i]; j++) {
      for (k = 0; k < 4; k++) {
        fscanf(fp, "%lf", &(phmm->pwm[i][j][k]));
      }
      fscanf(fp,"\n");
    }
  }
  totalStates *= (phmm->inactive * 2);
  phmm->extraState = phmm->N - totalStates;
  phmm->mean_matrix = gsl_matrix_alloc(phmm->K, phmm->N);
  phmm->var_matrix = gsl_matrix_alloc(phmm->K, phmm->N);
  phmm->pi = (double *) dvector(phmm->N);
  phmm->bg = (double *) dvector(4);
  /* Read initial probability for each hidden state */
  fscanf(fp, "pi:\n");
  for (i = 0; i < phmm->N; i++){
    fscanf(fp, "%lf\t", &(phmm->pi[i]));
  }
  /* Read initial means for each hidden state and each feature */
  fscanf(fp, "mu:\n");
  for (i = 0; i < phmm->K; i++) {
    for (j = 0; j < phmm->N; j++) {
      fscanf(fp, "%lf", &tmp);
      gsl_matrix_set(phmm->mean_matrix, i, j, tmp/2.0);
    }
    fscanf(fp,"\n");
  }

  /* Set initial covariance matrix */
  for (i = 0; i < phmm->K; i++) {
    for (j = 0; j < phmm->N ; j++) {
      gsl_matrix_set(phmm->var_matrix, i, j, 5.0);
    }
  }
  for (i = phmm->M; i < phmm->K; i++) {
    for (j = 0; j < phmm->N ; j++) {
      gsl_matrix_set(phmm->var_matrix, i, j, 5.0);
    }
  }
  if (phmm->model == 1 || phmm->model == 2){
    phmm->rho = (double **) dmatrix(phmm->K*(phmm->K-1)/2, phmm->N);
    phmm->cov_matrix = (gsl_matrix **)malloc((phmm->N)*sizeof(gsl_matrix *));
    for (i = 0; i < phmm->K*(phmm->K-1)/2; i++) {
      for (j = 0; j < phmm->N ; j++) {
        phmm->rho[i][j] = 0.5;
      }
    }
    for (i = 0; i < phmm->N; i++){
      phmm->cov_matrix[i] = gsl_matrix_alloc(phmm->K, phmm->K);
      covarMatrix_GSL(phmm, i, phmm->cov_matrix[i]);
    }
  }
}

/* Read parameters in a trained model file */
void ReadHMM(FILE *fp, HMM *phmm)
{
  int i, j, k, n, totalStates = 0;
  double sum;
  double tmp;  
  
  fscanf(fp, "N= %d\n", &(phmm->N));
  /* Read transition matrix */
  fscanf(fp, "A:\n");
  phmm->log_A_matrix = gsl_matrix_alloc(phmm->N, phmm->N);
  for (i = 0; i < phmm->N; i++) {
    for (j = 0; j < phmm->N; j++) {
      fscanf(fp, "%lf", &tmp);
      gsl_matrix_set(phmm->log_A_matrix, i, j, log(tmp));
    }
    fscanf(fp,"\n");
  }
  /* Read PWMs */
  phmm->D = (int *) ivector(phmm->M);
  phmm->pwm = (double ***)malloc(phmm->M*sizeof(double**));
  for (i = 0; i < phmm->M; i++) {
    fscanf(fp, "PWM: n=%d\n", &(phmm->D[i]));
    totalStates += phmm->D[i];
    phmm->pwm[i] = (double **) dmatrix(phmm->D[i], 4);
    for (j = 0; j < phmm->D[i]; j++) {
      for (k = 0; k < 4; k++) {
        fscanf(fp, "%lf", &(phmm->pwm[i][j][k]));
      }
      fscanf(fp,"\n");
    }
  }
  totalStates *= (phmm->inactive * 2);
  phmm->extraState = phmm->N - totalStates;
  phmm->mean_matrix = gsl_matrix_alloc(phmm->K, phmm->N);
  phmm->var_matrix = gsl_matrix_alloc(phmm->K, phmm->N);
  phmm->pi = (double *) dvector(phmm->N);
  phmm->bg = (double *) dvector(4);
  /* Read initial probability for each hidden state */
  fscanf(fp, "pi:\n");
  for (i = 0; i < phmm->N; i++){
    fscanf(fp, "%lf\t", &(phmm->pi[i])); 
  }
  /* Read means for each hidden state and each feature */
  fscanf(fp, "mu:\n");
  for (i = 0; i < phmm->K; i++) {
    for (j = 0; j < phmm->N; j++) {
      fscanf(fp, "%lf", &tmp);
      gsl_matrix_set(phmm->mean_matrix, i, j, tmp);  
    }
    fscanf(fp,"\n");
  }
  /* Read variance, correlation and set covariance matrix */
  fscanf(fp, "sigma:\n");
  for (i = 0; i < phmm->K; i++) {
    for (j = 0; j < phmm->N ; j++) {
      fscanf(fp, "%lf", &tmp); 
      gsl_matrix_set(phmm->var_matrix, i, j, tmp);
    }
    fscanf(fp,"\n");
  }
  if (phmm->model == 1 || phmm->model == 2){
  phmm->rho = (double **) dmatrix(phmm->K*(phmm->K-1)/2, phmm->N);
  phmm->cov_matrix = (gsl_matrix **)malloc((phmm->N)*sizeof(gsl_matrix *));
  fscanf(fp, "rho:\n");
  for (i = 0; i < phmm->K*(phmm->K-1)/2; i++) {
    for (j = 0; j < phmm->N ; j++) {
      fscanf(fp, "%lf", &(phmm->rho[i][j])); 
    }
    fscanf(fp,"\n");
  }
  for (i = 0; i < phmm->N; i++){
    phmm->cov_matrix[i] = gsl_matrix_alloc(phmm->K, phmm->K);
    covarMatrix_GSL(phmm, i, phmm->cov_matrix[i]);
  }
  }
}

/* Free menmory */
void FreeHMM(HMM *phmm)
{
  int i;
  for (i = 0; i < phmm->M; i++) {
    free_dmatrix(phmm->pwm[i], phmm->D[i], 4);
  }
  free_dvector(phmm->pi, phmm->N);
  free_dvector(phmm->bg, 4);
  gsl_matrix_free(phmm->mean_matrix);
  gsl_matrix_free(phmm->log_A_matrix);
  if (phmm->model == 0){
    gsl_matrix_free(phmm->var_matrix);
  }
  else if (phmm->model == 1 || phmm->model == 2){
    free_dmatrix(phmm->rho, phmm->K*(phmm->K-1)/2, phmm->N);
    for (i = 0; i < phmm->N; i++){
      gsl_matrix_free(phmm->cov_matrix[i]);
    }
  }
}

/*print the entire model as the correct format that ReadOutHMM can read*/
void PrintHMM(FILE *fp, HMM *phmm)
{
  int i, j, k, n, m, x;
	fprintf(fp, "M= %d\n", phmm->M);
	fprintf(fp, "N= %d\n", phmm->N);
	fprintf(fp, "A:\n");
  for (i = 0; i < phmm->N; i++) {
    for (j = 0; j < phmm->N; j++) {
      fprintf(fp, "%e ", exp(gsl_matrix_get(phmm->log_A_matrix, i, j)));
    }
    fprintf(fp, "\n");
  }
  for (i = 0; i < phmm->M; i++) {
	  fprintf(fp, "PWM: n=%d\n", (phmm->D[i]));
    for (j = 0; j < phmm->D[i]; j++) {
      for (k = 0; k < 4; k++){
        fprintf(fp, "%f ", phmm->pwm[i][j][k]);
	    }
      fprintf(fp, "\n");
    }
  }
  fprintf(fp, "pi:\n");
  for (i = 0; i < phmm->N; i++) {
	  fprintf(fp, "%f ", phmm->pi[i]);
	}
  fprintf(fp, "\n");
  
  fprintf(fp, "mu:\n");
  for (i = 0; i < phmm->K; i++) {
    for (j = 0; j < phmm->N; j++) {
      fprintf(fp, "%e ", gsl_matrix_get(phmm->mean_matrix, i, j));
    }
    fprintf(fp, "\n");
  }
  if (phmm->model == 0){
  fprintf(fp, "sigma:\n");
  for (i = 0; i < phmm->K; i++) {
    for (j = 0; j < phmm->N ; j++) {
      fprintf(fp, "%e ", gsl_matrix_get(phmm->var_matrix, i, j));
    }
    fprintf(fp, "\n");
  }
  }
  else if (phmm->model == 1 || phmm->model == 2){
  for (i = 0; i < phmm->N; i++) {
    for (n = 0; n < phmm->K; n++) {
      gsl_matrix_set(phmm->var_matrix, n, i, sqrt(gsl_matrix_get(phmm->cov_matrix[i], n, n)));
    }
  }
  for (i = 0; i < phmm->N; i++) {
    x = 0;
    for (n = 0; n < phmm->K - 1; n++){
      for (m = n + 1; m < phmm->K; m++){
        phmm->rho[x][i] = (gsl_matrix_get(phmm->cov_matrix[i], n, m)) / 
                          (gsl_matrix_get(phmm->var_matrix, n, i) * 
                          gsl_matrix_get(phmm->var_matrix, m, i)); 
        x++;//phmm->K*(n-1)+m-n*(n+1)/2
      }
    }
  }
  fprintf(fp, "sigma:\n");
  for (i = 0; i < phmm->K; i++) {
    for (j = 0; j < phmm->N ; j++) {
      fprintf(fp, "%e ", gsl_matrix_get(phmm->var_matrix, i, j));
    }
    fprintf(fp, "\n");
  }
  fprintf(fp, "rho:\n");
  for (i = 0; i < phmm->K*(phmm->K-1)/2; i++) {
    for (j = 0; j < phmm->N ; j++) {
      fprintf(fp, "%e ", phmm->rho[i][j]);
    }
    fprintf(fp, "\n");
  }
  fprintf(fp, "\n");
  }
}

/* Calculate corralation */
void getRho(HMM *phmm){
  int i ,n, m, x;
  for (i = 0; i < phmm->N; i++) {
    for (n = 0; n < phmm->K; n++) {
      gsl_matrix_set(phmm->var_matrix, n, i, sqrt(gsl_matrix_get(phmm->cov_matrix[i], n, n)));
    }
  }
  
  for (i = 0; i < phmm->N; i++) {
    x = 0;
    for (n = 0; n < phmm->K - 1; n++){
      for (m = n + 1; m < phmm->K; m++){
        phmm->rho[x][i] = (gsl_matrix_get(phmm->cov_matrix[i], n, m)) / 
                          (gsl_matrix_get(phmm->var_matrix, n, i) * 
                          gsl_matrix_get(phmm->var_matrix, m, i)); 
        x++;//phmm->K*(n-1)+m-n*(n+1)/2
      }
    }
  }
}




  