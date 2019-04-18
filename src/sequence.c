#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_randist.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>

#include "nrutil.h"
#include "hmm.h"
#include "logmath.h"
#include <omp.h>


static char rcsid[] = "$Id: sequence.c,v 1.2 1998/02/23 06:19:41 kanungo Exp kanungo $";

int gsl_ran_multivariate_gaussian_log_pdf (const gsl_vector * x,
                                           const gsl_vector * mu,
                                           const gsl_matrix * L,
                                           double * result,
                                           gsl_vector * work);

int gsl_linalg_cholesky_decomp_check (gsl_matrix * A, int *error_row, char *tmp_str, gsl_matrix * covar, HMM* phmm, int state);
static inline 
double
quiet_sqrt (double x)  
     /* avoids runtime error, for checking matrix for positive definiteness */
{
  return (x >= 0) ? sqrt(x) : GSL_NAN;
}

void ReadSequence(FILE *fp, int *pT, double *GC, int **pO, int *pP, int **peakPos)
{
  int *O, *peaks;
  //double *X;
  int i;
  fscanf(fp, "T= %d\n", pT);
  //X = dvector(4);
  fscanf(fp, "GC: ");
  for (i = 0; i < 4; i++) {
	  fscanf(fp, "%lf\t", &GC[i]); 
	}
	fscanf(fp,"\n");
  //*GC = X;
  O = ivector(*pT);
  for (i=0; i < *pT; i++) {
    fscanf(fp,"%d", &O[i]);
  }
  fscanf(fp,"\n");
  *pO = O;
  fscanf(fp, "P= %d\n", pP);
  peaks = ivector(*pP + 1);
  for (i=0; i < *pP + 1; i++){
    fscanf(fp,"%d", &peaks[i]);
  }
  *peakPos = peaks;
  //fprintf(stdout, "c: %d %d ", *pP, *pT);
}


void ReadTagFile(FILE *fp, int T, gsl_vector * data_vector, double adjust)
{
  double tmp;
  int i;
 
  for (i=0; i < T; i++) {
    fscanf(fp,"%lf\t", &tmp);
    gsl_vector_set(data_vector, i, tmp*adjust);
  }
}


void CalMotifScore_P(HMM *phmm, gsl_matrix * S, int *O1, int P, int *peakPos)
{
  int thread_id, nloops;
  int k, t, j, i, n, m, D;
  int start, end;
  double tempF, tempR, tmpL, tmpR, tmp;
  gsl_vector * tempList;

  for (m = 0; m < phmm->M; m++) {
    D = phmm->D[m];
    
#pragma omp parallel num_threads(THREAD_NUM) \
  private(thread_id, nloops, start, end, tempF, tempR, tmpL, tmpR, tmp, \
  t, i, j, n, tempList) 
  {
    nloops = 0;
#pragma omp for
    for (k = 0; k < P; k++) {
      ++nloops;
      tempList = gsl_vector_alloc(D);
      start = peakPos[k];
      end = peakPos[k+1] - 1;

      tmpL = -INFINITY;
      t = start - 1;
      i = 1;
      tempF = tempR = 0.0;
      for (j = 1; j <= D; j++){
        tempF += log_2(phmm->pwm[m][j-1][O1[t-i+j]]);
        tempR += log_2(phmm->pwm[m][D-j][3-O1[t-i+j]]);
        tempF -= log_2(phmm->bg[O1[t-i+j]]);
        tempR -= log_2(phmm->bg[3-O1[t-i+j]]);
      }
      tmpL=MAX(tmpL, MAX(tempF, tempR));
      tmpR = -INFINITY;
      t = end - D;
      i = 1;
      tempF = tempR = 0.0;
      for (j = 1; j <= D; j++){
        tempF += log_2(phmm->pwm[m][j-1][O1[t-i+j]]);
        tempR += log_2(phmm->pwm[m][D-j][3-O1[t-i+j]]);
        tempF -= log_2(phmm->bg[O1[t-i+j]]);
        tempR -= log_2(phmm->bg[3-O1[t-i+j]]);
      }
      tmpR=MAX(tmpR, MAX(tempF, tempR));
        
      t = start + D - 2;
      for (i = 1; i <= D; i++){
        tempF = tempR = 0.0;
        for (j = 1; j <= D; j++){
          tempF += log_2(phmm->pwm[m][j-1][O1[t-i+j]]);
          tempR += log_2(phmm->pwm[m][D-j][3-O1[t-i+j]]);
          tempF -= log_2(phmm->bg[O1[t-i+j]]);
          tempR -= log_2(phmm->bg[3-O1[t-i+j]]);
        }  
        gsl_vector_set(tempList, i-1, MAX(tempF, tempR));
      } 
      for (t = start - 1; t < end; t++) {
        
        if (t < (start + D - 1)){
          gsl_matrix_set (S, m, t, tmpL);
        }
        else if (t > (end - D - 1)) {
          gsl_matrix_set (S, m, t, tmpR);
        }
        else{
        //
          tempF = tempR = 0.0;
          i = 1;
          for (j = 1; j <= D; j++){
            tempF += log_2(phmm->pwm[m][j-1][O1[t-i+j]]);
            tempR += log_2(phmm->pwm[m][D-j][3-O1[t-i+j]]);
            tempF -= log_2(phmm->bg[O1[t-i+j]]);
            tempR -= log_2(phmm->bg[3-O1[t-i+j]]);
          }
          for (n = D-1; n > 0; n--){
            gsl_vector_swap_elements(tempList, n-1, n);
          } 
          gsl_vector_set(tempList, 0, MAX(tempF, tempR));
          gsl_matrix_set (S, m, t, gsl_vector_max(tempList));
          if (!(gsl_matrix_get(S, m, t)==gsl_matrix_get(S, m, t))){
            fprintf(stdout, "pwmNAN\t %d %d %d %lf ",t,m,phmm->D[m],gsl_matrix_get(S, m, t));
          }
          //
        }
      }
      free(tempList);
    }
    thread_id = omp_get_thread_num();
  }
    
  }
  
}


void EmissionMatrix(HMM* phmm, gsl_matrix * obs_matrix, int P, int *peakPos, 
                    gsl_matrix * emission_matrix, int T)
{
  int i, j, nloops, thread_id;
  gsl_vector * tmp_vector;// = gsl_vector_alloc(phmm->K);
  gsl_vector * tmp_vector_2;// = gsl_vector_alloc(T);
  gsl_vector * tmp_vector_3;
  gsl_vector * mean_vector;// = gsl_vector_alloc(phmm->K);
  gsl_vector * var_vector;// = gsl_vector_alloc(phmm->K);
  gsl_matrix * tmp_matrix;// = gsl_matrix_alloc(phmm->K, T);
  gsl_matrix * tmp_matrix_2;// = gsl_matrix_alloc(phmm->K, T);
  
  //gsl_matrix * tmp_mean_matrix;// = gsl_matrix_alloc(phmm->K, T);
  //gsl_matrix * tmp_var_matrix;// = gsl_matrix_alloc(phmm->K, T);
#pragma omp parallel num_threads(THREAD_NUM) \
  private(thread_id, nloops, i, j, tmp_vector, mean_vector, var_vector, tmp_matrix, tmp_matrix_2, tmp_vector_2, tmp_vector_3)
  { 
#pragma omp for    
  for (i = 0; i < phmm->N; i++){
    //printf("0: %d ", i);
    //fflush(stdout); 
    mean_vector = gsl_vector_alloc(phmm->K);
    var_vector = gsl_vector_alloc(phmm->K);
    tmp_vector = gsl_vector_alloc(phmm->K);
    tmp_vector_2 = gsl_vector_alloc(phmm->K);
    tmp_vector_3 = gsl_vector_alloc(phmm->K);
    //tmp_mean_matrix = gsl_matrix_alloc(phmm->K, T);
    //tmp_var_matrix = gsl_matrix_alloc(phmm->K, T);
    
    gsl_matrix_get_col(mean_vector, phmm->mean_matrix, i);  
    gsl_matrix_get_col(var_vector, phmm->var_matrix, i); 
    for (j = 0; j < phmm->K; j++){
      gsl_vector_set(tmp_vector_3, j, log((1.0/gsl_vector_get(var_vector, j))/SQRT_TWO_PI));
    }
    //printf("1: %d ", i);
    //fflush(stdout); 
    
    //printf("2: %d ", i);
    //fflush(stdout); 
    tmp_matrix = gsl_matrix_alloc(phmm->K, T);
    for (j = 0; j < T; j++){
      //gsl_matrix_set_col(tmp_mean_matrix, j, mean_vector); 
      gsl_matrix_get_col(tmp_vector, obs_matrix, j);  
      gsl_vector_sub(tmp_vector, mean_vector);
      gsl_vector_memcpy(tmp_vector_2, tmp_vector);
      gsl_vector_mul(tmp_vector, tmp_vector_2);
      gsl_vector_scale(tmp_vector, -0.5);
      gsl_vector_div(tmp_vector, var_vector);
      gsl_vector_div(tmp_vector, var_vector);
      gsl_vector_add(tmp_vector, tmp_vector_3);
      gsl_matrix_set_col(tmp_matrix, j, tmp_vector);
      //gsl_matrix_set_col(tmp_var_matrix, j, var_vector); 
      //gsl_matrix_set_col(tmp_matrix_2, j, tmp_vector); 
    }
    gsl_vector_free(mean_vector);
    //gsl_vector_free(tmp_vector);
    gsl_vector_free(tmp_vector_2);
    gsl_vector_free(tmp_vector_3);
    //gsl_vector_free(var_vector);
    //printf("3: %d ", i);
    //fflush(stdout); 
    //tmp_matrix = gsl_matrix_alloc(phmm->K, T);
    
    //gsl_matrix_memcpy(tmp_matrix, obs_matrix);
    //gsl_matrix_sub(tmp_matrix, tmp_mean_matrix);
    //gsl_matrix_free(tmp_mean_matrix);
    
    //gsl_matrix_mul_elements(tmp_matrix, tmp_matrix);
    //printf("4: %d ", i);
    //fflush(stdout); 
   // for (j = 0; j < T; j++){
   //   gsl_matrix_set_col(tmp_var_matrix, j, var_vector); 
   // }
    
    //gsl_vector_free(var_vector);
    //gsl_matrix_div_elements(tmp_matrix, tmp_var_matrix);
    //gsl_matrix_div_elements(tmp_matrix, tmp_var_matrix);
    //gsl_matrix_scale(tmp_matrix, -0.5);
    
    //gsl_matrix_free(tmp_mean_matrix);
    //gsl_matrix_free(tmp_var_matrix);
    //printf("5: %d ", i);
    //fflush(stdout); 
    //tmp_matrix_2 = gsl_matrix_alloc(phmm->K, T);
    //for (j = 0; j < T; j++){
      //gsl_matrix_set_col(tmp_matrix_2, j, tmp_vector); 
    //}
    //gsl_matrix_add(tmp_matrix, tmp_matrix_2);
    
    //gsl_matrix_free(tmp_matrix_2);
    //tmp_vector_2 = gsl_vector_alloc(T);
    //tmp_vector = gsl_vector_alloc(T);
    gsl_vector_set_all(tmp_vector, 1.0);
    //printf("6: %d ", i);
    //fflush(stdout); 
    tmp_vector_2 = gsl_vector_alloc(T);
    gsl_blas_dgemv(CblasTrans, 1.0, tmp_matrix, tmp_vector, 0.0, tmp_vector_2);
    gsl_matrix_set_row(emission_matrix, i, tmp_vector_2); 
    
    gsl_vector_free(tmp_vector);
    gsl_vector_free(tmp_vector_2);
    gsl_matrix_free(tmp_matrix);
    
    
  }
  }
    /*
    for (j = 0; j < T; j++){    
    gsl_matrix_get_col(tmp_vector, obs_matrix, j);  
    gsl_vector_sub(tmp_vector, mean_vector);
    gsl_vector_memcpy(tmp_vector_2, tmp_vector);
    gsl_vector_mul(tmp_vector, tmp_vector_2);
    gsl_vector_scale(tmp_vector, -1.0);
    gsl_vector_memcpy(tmp_vector_2, var_vector);
    gsl_vector_mul(tmp_vector_2, var_vector);
    gsl_vector_scale(tmp_vector_2, 2.0);
    gsl_vector_div(tmp_vector, tmp_vector_2);
    
    gsl_matrix_set_col(tmp_matrix, j, tmp_vector);
    */
    //gsl_vector_free(tmp_vector);
    //gsl_vector_free(tmp_vector_2);
    //gsl_matrix_free(tmp_matrix);
    //gsl_matrix_free(tmp_matrix_2);
}
  
     
                    
void EmissionMatrix_GSL(HMM* phmm, gsl_matrix * obs_matrix, int P, int *peakPos, 
                    gsl_matrix * emission_matrix, int T)
{
  int i, j, t, nloops, thread_id;
  double mean, sd, data_mean, *emission;
  
#pragma omp parallel num_threads(THREAD_NUM) \
  private(thread_id, nloops, i, j, t, emission)
  {
    nloops = 0;
    
#pragma omp for        
    for (i = 0; i < phmm->N; i++){
      emission = dvector(T);
      for (j = 0; j < phmm->K; j++){  
        mean =  gsl_matrix_get(phmm->mean_matrix, j, i);
        sd = gsl_matrix_get(phmm->var_matrix, j, i);
        for (t = 0; t < T; t++){                
          data_mean = gsl_matrix_get(obs_matrix, j, t) - mean;    
          emission[t] += log(gsl_ran_gaussian_pdf(data_mean, sd));  
        } 
      }
      for (t = 0; t < T; t++){
        gsl_matrix_set(emission_matrix, i, t, emission[t]);  
      }
      free_dvector(emission, T);
    }
  }
}
                    
void EmissionMatrix_mv(HMM* phmm, gsl_matrix * obs_matrix, int P, int *peakPos, 
                       gsl_matrix * emission_matrix, int T)
{
  int thread_id, nloops;
  int k, t, i, j, start, end, m, n;
  //XXX:work on bivariance later
  /*
  if (phmm->K == 1){
    data = dvector(phmm->K);
    for (k = 0; k < P; k++){
      start = peakPos[k];
      end = peakPos[k+1] - 1;
      for (t = start-1; t < end; t++) {
        data[0] = Obs[0][t];
        for (i = 0; i < phmm->N; i++){
          //emission[t][i]=NormDist(phmm->mu[0], i, phmm->sigma[0], Obs[t][0])+0.000000000000001;
          emission_GSL[i][t]=NormDist(phmm->mu[0], i, phmm->sigma[0], data[0])+0.000000000000001;
          //ComputeEmission(phmm, i, Obs[t]);
        }
      }
    }
  }
  else if (phmm->K == 2){
    data = dvector(phmm->K);
    for (k = 0; k < P; k++){
      start = peakPos[k];
      end = peakPos[k+1] - 1;
      for (t = start-1; t < end; t++) {
        data[0] = Obs[0][t];
        data[1] = Obs[1][t];
        for (i = 0; i < phmm->N; i++){
          //emission[t][i]=BiVarNormDist(phmm->mu, i, phmm->sigma, phmm->rho, Obs[t])+0.000000000000001;
          emission_GSL[i][t]=BiVarNormDist(phmm->mu, i, phmm->sigma, phmm->rho, data)+0.000000000000001;
          //ComputeEmission(phmm, i, Obs[t]);
        }
      }
    }
  }
  */
  if (phmm->K > 2){
    gsl_matrix * cov_matrix_tmp[phmm->N];
    gsl_matrix * tmp_matrix;
    gsl_vector * mu_vector[phmm->N]; 
    int n = -1;
    int l;
    char tmp_str[1000]; 
    int error_row;
    int x, y;
    //for (i = 0; i < phmm->N; i++) error_row[i] = -1;
#pragma omp parallel num_threads(THREAD_NUM) \
  private(thread_id, nloops, j, n, l, x, y, tmp_str, tmp_matrix, error_row)
  {
    nloops = 0;
    
#pragma omp for   
    for (i = 0; i < phmm->N; i++){
      //printf("1: %d ", i);
      //fflush(stdout);   
      cov_matrix_tmp[i] = gsl_matrix_alloc(phmm->K, phmm->K);
      gsl_matrix_memcpy(cov_matrix_tmp[i], phmm->cov_matrix[i]);
      
      mu_vector[i] = gsl_vector_alloc(phmm->K);
      gsl_matrix_get_col(mu_vector[i], phmm->mean_matrix, i);
      
      sprintf(tmp_str, "%s%d", "./cov/",i);

      strcat(tmp_str,".txt");
      
      gsl_set_error_handler_off ();
      gsl_linalg_cholesky_decomp_check(cov_matrix_tmp[i], &error_row, tmp_str, phmm->cov_matrix[i], phmm, i);
      //printf("5: %d ", i);
      //fflush(stdout);   
      
      //gsl_linalg_cholesky_decomp(cov_matrix_tmp[i]);
    } 
  }
  
  gsl_vector * data_vector;  
  gsl_vector * workspace; 
  double tmp;
  
#pragma omp parallel num_threads(THREAD_NUM) \
  private(thread_id, nloops, start, end, i, t,workspace, tmp, data_vector) 
  {
    nloops = 0;
#pragma omp for   
    for (k = 0; k < P; k++){
      //printf("k: %d ", k);
      //fflush(stdout); 
      ++nloops;
      start = peakPos[k];
      end = peakPos[k+1] - 1;
      data_vector = gsl_vector_alloc(phmm->K);
      workspace = gsl_vector_alloc(phmm->K);
      for (t = start-1; t < end; t++) {
        gsl_matrix_get_col(data_vector, obs_matrix, t);
        for (i = 0; i < phmm->N; i++){
          gsl_ran_multivariate_gaussian_log_pdf(data_vector, mu_vector[i], cov_matrix_tmp[i], &tmp, workspace);
          gsl_matrix_set (emission_matrix, i, t, tmp);
        }
      }
      gsl_vector_free(workspace);
      gsl_vector_free(data_vector);  
    }
    thread_id = omp_get_thread_num();
  }
  
    for (i = 0; i < phmm->N; i++){
      gsl_vector_free(mu_vector[i]);
      gsl_matrix_free(cov_matrix_tmp[i]);
    }
    
  } 
}

void EmissionMatrix_mv_reduce(HMM* phmm, gsl_matrix * obs_matrix, int P, int *peakPos, 
                       gsl_matrix * emission_matrix, int T)
{
  int thread_id, nloops;
  int k, t, i, j, start, end, m, n;
  int TF = 0;
  int *stateList = ivector(phmm->N); //TODO: adapt to one states
  for (j = 0; j < phmm->M; j++){
    for (i = TF; i < TF + phmm->D[j]; i++) {
      stateList[i] = j;
    }
    TF += phmm->D[j];
    fprintf(stdout,"%d ", stateList[TF-1]);
    if (phmm->inactive == 1){
      for (i = TF; i < TF + phmm->D[j]; i++) {
        stateList[i] = j;
      }
      TF += phmm->D[j];
      fprintf(stdout,"%d ", stateList[TF-1]);
    }
  }
  TF -= 1;

  for (j = phmm->N - phmm->extraState; j < phmm->N; j++){
    stateList[j] = phmm->M;
    fprintf(stdout,"%d ", stateList[j]);
  }

  if (phmm->K > 1) {
    gsl_matrix *cov_matrix_tmp[phmm->N];
    gsl_matrix *tmp_matrix;
    gsl_vector *mu_vector[phmm->N];
    gsl_vector *tmp_vector;
    gsl_vector *tmp_vector2;
    gsl_vector *deleted_vector[phmm->N];
    const int FULL = -1;
    const int EMPTY = -2;
    int n;
    int l;
    char tmp_str[100];
    int error_row;
    int x, y;

#pragma omp parallel num_threads(THREAD_NUM) \
  private(thread_id, nloops, j, n, l, x, y, tmp_str, tmp_matrix, error_row, tmp_vector, tmp_vector2)
  {
    nloops = 0;
    
#pragma omp for   
    for (i = 0; i < phmm->N; i++){
      
      n = -1;
      
      cov_matrix_tmp[i] = gsl_matrix_alloc(phmm->K, phmm->K);
      gsl_matrix_memcpy(cov_matrix_tmp[i], phmm->cov_matrix[i]);
      
      mu_vector[i] = gsl_vector_alloc(phmm->K);
      gsl_matrix_get_col(mu_vector[i], phmm->mean_matrix, i);
      
      deleted_vector[i] = gsl_vector_alloc(1);
      gsl_vector_set(deleted_vector[i], 0, -1);
      
      sprintf(tmp_str, "%s%d", "./cov/",i);
      strcat(tmp_str,".txt");
      //if (phmm->K > 2) {
        do {
          n++;
          //gsl_matrix_memcpy(cov_matrix_tmp[i], tmp_matrix);

          if (n != 0) {
            if (n == 1) {
              gsl_vector_set(deleted_vector[i], 0, error_row);
            } else {
              tmp_vector2 = gsl_vector_alloc(n - 1);
              gsl_vector_memcpy(tmp_vector2, deleted_vector[i]);
              gsl_vector_free(deleted_vector[i]);
              deleted_vector[i] = gsl_vector_alloc(n);
              for (j = 0; j < n - 1; j++) gsl_vector_set(deleted_vector[i], j,
                                                         gsl_vector_get(
                                                                 tmp_vector2,
                                                                 j));
              gsl_vector_free(tmp_vector2);
              gsl_vector_set(deleted_vector[i], n - 1, error_row);
            }
            //printf("1: %d %d ", i, error_row);
            //fflush(stdout);
            gsl_matrix_free(cov_matrix_tmp[i]);
            cov_matrix_tmp[i] = gsl_matrix_alloc(phmm->K - n, phmm->K - n);
            tmp_vector = gsl_vector_alloc(phmm->K - n + 1);
            gsl_vector_memcpy(tmp_vector, mu_vector[i]);
            gsl_vector_free(mu_vector[i]);
            mu_vector[i] = gsl_vector_alloc(phmm->K - n);
            x = 0;
            for (j = 0; j < phmm->K - n + 1; j++) {
              if (j != error_row) {
                gsl_vector_set(mu_vector[i], x, gsl_vector_get(tmp_vector, j));
                y = 0;
                for (l = 0; l < phmm->K - n + 1; l++) {
                  if (l != error_row) {
                    gsl_matrix_set(cov_matrix_tmp[i], x, y,
                                   gsl_matrix_get(tmp_matrix, j, l));
                    y++;
                  }
                }
                x++;
              }
            }
            gsl_matrix_free(tmp_matrix);
            gsl_vector_free(tmp_vector);

          }

          tmp_matrix = gsl_matrix_alloc(phmm->K - n, phmm->K - n);
          gsl_matrix_memcpy(tmp_matrix, cov_matrix_tmp[i]);
          //gsl_matrix_free(tmp_matrix);


          error_row = FULL;
          gsl_set_error_handler_off();
          gsl_linalg_cholesky_decomp_check(cov_matrix_tmp[i], &error_row,
                                           tmp_str, phmm->cov_matrix[i], phmm,
                                           i);
          if (error_row == 0) {
            gsl_vector_set(deleted_vector[i], 0, EMPTY);
            error_row = FULL;
          }
        } while (error_row != FULL);
        gsl_matrix_free(tmp_matrix);
      //}
      //else{
        //gsl_linalg_cholesky_decomp_check(cov_matrix_tmp[i], &error_row,
          //                               tmp_str, phmm->cov_matrix[i], phmm,
            //                             i);
        //gsl_vector_set(deleted_vector[i], 0, FULL);
      //}
      //gsl_linalg_cholesky_decomp(cov_matrix_tmp[i]);
    } 
  }
  
  //gsl_vector_free(tmp_vector);
  gsl_vector * data_vector;  
  gsl_vector * deleted_data_vector;  
  //gsl_vector * tmp_vector;  
  gsl_vector * workspace; 
  double tmp;
  double mean, sd, data_mean, emission;

#pragma omp parallel num_threads(THREAD_NUM) \
  private(thread_id, nloops, start, end, i, j, t, m, x, workspace, tmp, data_vector, tmp_vector, deleted_data_vector,mean,sd,data_mean,emission) 
  {
    nloops = 0;
#pragma omp for   
    for (k = 0; k < P; k++){
      //printf("k: %d ", k);
      //fflush(stdout); 
      ++nloops;
      start = peakPos[k];
      end = peakPos[k+1] - 1;
      data_vector = gsl_vector_alloc(phmm->K);
      workspace = gsl_vector_alloc(phmm->K);
      for (t = start-1; t < end; t++) {
        gsl_matrix_get_col(data_vector, obs_matrix, t);
        for (i = 0; i < phmm->N; i++){
          if(phmm->thresholds[stateList[i]] != -INFINITY && i < (phmm->N - phmm->extraState) && gsl_vector_get(data_vector, stateList[i]) < phmm->thresholds[stateList[i]]){
            tmp = -INFINITY;
          }
          else if (gsl_vector_get(deleted_vector[i], 0) == FULL) {
            gsl_ran_multivariate_gaussian_log_pdf(data_vector, mu_vector[i], cov_matrix_tmp[i], &tmp, workspace);
          }
          else if (gsl_vector_get(deleted_vector[i], 0) == EMPTY) {
            tmp = -INFINITY;
          }
          else {
          /*mean =  gsl_matrix_get(phmm->mean_matrix, j, i);
        sd = gsl_matrix_get(phmm->var_matrix, j, i);
        for (t = 0; t < T; t++){                
          data_mean = gsl_matrix_get(obs_matrix, j, t) - mean;    
          emission[t] += log(gsl_ran_gaussian_pdf(data_mean, sd));  
          */
            emission = 0.0;
            deleted_data_vector = gsl_vector_alloc(phmm->K-deleted_vector[i]->size);
            //deleted_data_vector = gsl_vector_alloc(phmm->K);
            //gsl_vector_memcpy(deleted_data_vector, data_vector);
          
            //for (j = 0; j < deleted_vector[i]->size; j++){
              //tmp_vector = gsl_vector_alloc(phmm->K - j - 1);
              x = 0;
              j = 0;
              //for (m = 0; m < deleted_data_vector->size; m++) {
              for (m = 0; m < data_vector->size; m++) {
                //if (m != gsl_vector_get(deleted_vector[i], j)){
                if ((m-j) != gsl_vector_get(deleted_vector[i], j)){
                  //gsl_vector_set(tmp_vector, x, gsl_vector_get(deleted_data_vector, m));
                  gsl_vector_set(deleted_data_vector, x, gsl_vector_get(data_vector, m));
                  x++;
                }
                else {
                  mean =  gsl_matrix_get(phmm->mean_matrix, m-j, i);
                  sd = sqrt(gsl_matrix_get(phmm->cov_matrix[i], m-j, m-j));
                  if (sd == 0.0) {emission += -INFINITY;}
                  else {
                  data_mean = gsl_matrix_get(obs_matrix, m-j, t) - mean;    
                  emission += log(gsl_ran_gaussian_pdf(data_mean, sd));  
                  }
                  j++;
                   
                }
                
              }
              //gsl_vector_free(deleted_data_vector);
              //deleted_data_vector = gsl_vector_alloc(phmm->K - j - 1);
              //gsl_vector_memcpy(deleted_data_vector, tmp_vector);
              //gsl_vector_free(tmp_vector);
              
            gsl_ran_multivariate_gaussian_log_pdf(deleted_data_vector, mu_vector[i], cov_matrix_tmp[i], &tmp, workspace);  
            gsl_vector_free(deleted_data_vector);
            tmp += emission;
          }
          gsl_matrix_set (emission_matrix, i, t, tmp);

        }
        
      }
      gsl_vector_free(workspace);
      gsl_vector_free(data_vector);  
    }
    thread_id = omp_get_thread_num();
  }
  
    for (i = 0; i < phmm->N; i++){
      gsl_vector_free(mu_vector[i]);
      gsl_matrix_free(cov_matrix_tmp[i]);
      gsl_vector_free(deleted_vector[i]);
    }
    
  } 
}


void covarMatrix_GSL(HMM *phmm, int state, gsl_matrix * cov_matrix)
{
  int i, j, n, k;
  double corr;
  //for (n = 1; n <= phmm->N; n++){
  k = 0;
    for (i = 0; i < phmm->K; i++) {
      for (j = 0; j < phmm->K; j++) {
        if (i == j) {
          gsl_matrix_set (cov_matrix, i, j, 
                          (gsl_matrix_get(phmm->var_matrix, i, state) * 
                          gsl_matrix_get(phmm->var_matrix, i, state)));
          //matrix[i][j] = phmm->sigma[i][state] * phmm->sigma[j][state];
        }
        else if (i < j) {
        //phmm->K*(i-1)+j-i*(i+1)/2
          corr = phmm->rho[k][state];
          gsl_matrix_set (cov_matrix, i, j, 
                          (gsl_matrix_get(phmm->var_matrix, i, state) * 
                          gsl_matrix_get(phmm->var_matrix, j, state) * corr));
          //matrix[i][j] = phmm->sigma[i][state] * phmm->sigma[j][state] * corr;
          k++;
        }
        else {
          //corr = phmm->rho[phmm->K*(j-1)+i-j*(j+1)/2][state];
          gsl_matrix_set(cov_matrix, i, j, gsl_matrix_get(cov_matrix, j, i));
          //matrix[i][j] = matrix[j][i];//phmm->sigma[i][state] * phmm->sigma[j][state] * corr;
        }
      }
    }
  //}
}


int
gsl_ran_multivariate_gaussian_log_pdf (const gsl_vector * x,
                                       const gsl_vector * mu,
                                       const gsl_matrix * L,
                                       double * result,
                                       gsl_vector * work)
{
  const size_t M = L->size1;
  const size_t N = L->size2;

  if (M != N)
    {
      GSL_ERROR("requires square matrix", GSL_ENOTSQR);
    }
  else if (mu->size != M)
    {
      GSL_ERROR("incompatible dimension of mean vector with variance-covariance matrix", GSL_EBADLEN);
    }
  else if (x->size != M)
    {
      GSL_ERROR("incompatible dimension of quantile vector", GSL_EBADLEN);
    }
  else if (work->size != M)
    {
      GSL_ERROR("incompatible dimension of work vector", GSL_EBADLEN);
    }
  else
    {
      size_t i;
      double quadForm;        /* (x - mu)' Sigma^{-1} (x - mu) */
      double logSqrtDetSigma; /* log [ sqrt(|Sigma|) ] */

      /* compute: work = x - mu */
      for (i = 0; i < M; ++i)
        {
          double xi = gsl_vector_get(x, i);
          double mui = gsl_vector_get(mu, i);
          gsl_vector_set(work, i, xi - mui);
        }

      /* compute: work = L^{-1} * (x - mu) */
      gsl_blas_dtrsv(CblasLower, CblasNoTrans, CblasNonUnit, L, work);

      /* compute: quadForm = (x - mu)' Sigma^{-1} (x - mu) */
      gsl_blas_ddot(work, work, &quadForm);

      /* compute: log [ sqrt(|Sigma|) ] = sum_i log L_{ii} */
      logSqrtDetSigma = 0.0;
      for (i = 0; i < M; ++i)
        {
          double Lii = gsl_matrix_get(L, i, i);
          logSqrtDetSigma += log(Lii);
        }

      *result = -0.5*quadForm - logSqrtDetSigma - 0.5*M*log(2.0*M_PI);

      return GSL_SUCCESS;
    }
}

void PrintSequenceProb(FILE *fp, int T, int *O, double *vprob, double *g, 
                       double **posterior, int indexTF)
{
  int i;
  fprintf(fp,"%d\t%lf\t%lf\t%lf", O[0], vprob[0], g[0], posterior[0][indexTF]);
  for (i=1; i < T; i++) {
    fprintf(fp,"\n%d\t%lf\t%lf\t%lf", O[i], vprob[i], g[i], posterior[i][indexTF]);
  }
}


int
gsl_linalg_cholesky_decomp_check (gsl_matrix * A, int *error_row, char *tmp_str, gsl_matrix * covar, HMM* phmm, int state)
{
  const size_t M = A->size1;
  const size_t N = A->size2;
  int m, n;
  FILE	* tmp_fp;
  if (M != N)
    {
      GSL_ERROR("cholesky decomposition requires square matrix", GSL_ENOTSQR);
    }
  else
    {
      size_t i,j,k;
      int status = 0;

      /* Do the first 2 rows explicitly.  It is simple, and faster.  And
       * one can return if the matrix has only 1 or 2 rows.  
       */

      double A_00 = gsl_matrix_get (A, 0, 0);
      if (A_00 != A_00)
       {
          status = GSL_EDOM ;
          *error_row = 0;
          printf("nan %d\t", state);
          return GSL_SUCCESS;
       }
      double L_00 = quiet_sqrt(A_00);
      
      if (A_00 <= 0)
        {
          status = GSL_EDOM ;
          *error_row = 0;
          return GSL_SUCCESS;
        }

      gsl_matrix_set (A, 0, 0, L_00);
  
      if (M > 1)
        {
          double A_10 = gsl_matrix_get (A, 1, 0);
          double A_11 = gsl_matrix_get (A, 1, 1);
          
          double L_10 = A_10 / L_00;
          double diag = A_11 - L_10 * L_10;
          double L_11 = quiet_sqrt(diag);
          
          if (diag <= 0)
            {
              printf("%d %d\t", state,1);

              status = GSL_EDOM;
              *error_row = 1;
              return GSL_SUCCESS;
            }

          gsl_matrix_set (A, 1, 0, L_10);        
          gsl_matrix_set (A, 1, 1, L_11);
        }
      
      for (k = 2; k < M; k++)
        {
          double A_kk = gsl_matrix_get (A, k, k);
          
          for (i = 0; i < k; i++)
            {
              double sum = 0;

              double A_ki = gsl_matrix_get (A, k, i);
              double A_ii = gsl_matrix_get (A, i, i);

              gsl_vector_view ci = gsl_matrix_row (A, i);
              gsl_vector_view ck = gsl_matrix_row (A, k);

              if (i > 0) {
                gsl_vector_view di = gsl_vector_subvector(&ci.vector, 0, i);
                gsl_vector_view dk = gsl_vector_subvector(&ck.vector, 0, i);
                
                gsl_blas_ddot (&di.vector, &dk.vector, &sum);
              }

              A_ki = (A_ki - sum) / A_ii;
              gsl_matrix_set (A, k, i, A_ki);
            } 

          {
            gsl_vector_view ck = gsl_matrix_row (A, k);
            gsl_vector_view dk = gsl_vector_subvector (&ck.vector, 0, k);
            
            double sum = gsl_blas_dnrm2 (&dk.vector);
            double diag = A_kk - sum * sum;

            double L_kk = quiet_sqrt(diag);
            
            if (diag <= 0)
              {
                //printf("in1: %d ", state);
                //fflush(stdout);   
                //getRho(phmm);
                
                //printf("%d %d\t", state,k);
                
                //fflush(stdout);   
                //tmp_fp = fopen(tmp_str, "w");
                //fprintf(tmp_fp, "error_row %d \n", k);
                
                //printf("in2: %d ", state);
                //fflush(stdout);  
                /* 
                for (m = 0; m <= k; m++){
                  for (n = 0; n < M; n++){
                    //fprintf(tmp_fp, "cov %d %e\t", n, gsl_matrix_get(cov_matrix[i], m, n));
                    fprintf(tmp_fp, "cov %d %e\t", n, gsl_matrix_get(covar, m, n));
                  }
                  fprintf(tmp_fp, "\n");
                }
                //printf("in3: %d ", state);
                //fflush(stdout);   
                for (m = 0; m <= M; m++){
                  fprintf(tmp_fp, "cor %d %e\t", m, phmm->rho[phmm->K*k+m-(k+2)*(k+1)/2][state]);
                }
                */
                //fclose(tmp_fp);
                
                status = GSL_EDOM;
                *error_row = k;
                return GSL_SUCCESS;
              }
            
            gsl_matrix_set (A, k, k, L_kk);
          }
        }

      /* Now copy the transposed lower triangle to the upper triangle,
       * the diagonal is common.  
       */
      
      for (i = 1; i < M; i++)
        {
          for (j = 0; j < i; j++)
            {
              double A_ij = gsl_matrix_get (A, i, j);
              gsl_matrix_set (A, j, i, A_ij);
            }
        } 
      
      if (status == GSL_EDOM)
        {
          GSL_ERROR ("matrix must be positive definite", GSL_EDOM);
          
          
        }
      
      return GSL_SUCCESS;
    }
}
