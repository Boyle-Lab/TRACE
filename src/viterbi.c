#include <math.h>
#include "hmm.h"
#include "nrutil.h"
#include <omp.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

static char rcsid[] = "$Id: viterbi.c,v 1.1 1999/05/06 05:25:37 kanungo Exp kanungo $";

#define VITHUGE  100000000000.0

void Viterbi(HMM *phmm, int T, double *g, double  **alpha, double	**beta, 
             double	**gamma, double  *logprobf, double **delta, 
             int **psi, int *q, double *vprob, double *pprob, 
             double **posterior, int P, int *peakPos,
             gsl_matrix * emission_matrix)
{
  int thread_id, nloops;
  int     i, j, k, m, x, y, TF;   /* state indices */
  int     t;      /* time index */
 
  int     maxvalind;
  int	l = 0;
  int *TFlist, *lengthList;
  int nonInf;
  double  maxval, val;
	//double  **biot;
 
  
  double plogprobinit, plogprobfinal;
  t = 0;
  
  double temp;
  /*
  TFlist = ivector(phmm->M);
  TF = -1;
  for (j = 0; j < phmm->M; j++){
    TF += phmm->D[j];
    TFlist[j] = TF;
    fprintf(stdout,"%d ", TFlist[j]);
  }
  fprintf(stdout,"\n");  

  */
  TF = 0;
  TFlist = ivector(phmm->M * (phmm->inactive+1));
    for (j = 0; j < phmm->M; j++){
      TF += phmm->D[j];
      TFlist[j * (phmm->inactive+1)] = TF - 1;
      if (phmm->inactive == 1){
        TF += phmm->D[j];
        TFlist[j * (phmm->inactive+1) + 1] = TF - 1;
      } 
    }
  lengthList = ivector(phmm->M * (phmm->inactive+1));
  if (phmm->inactive == 1){
    for (j = 0; j < phmm->M; j++){
      lengthList[j*2] = phmm->D[j];
      lengthList[j*2+1] = phmm->D[j];
    }
  }
  else lengthList = phmm->D;
  //fflush(stdout);
/* 0. Preprocessing */
/*
  biot = dmatrix(phmm->N, T);
  for (i = 0; i < phmm->N; i++) {
    for (t = 0; t < T; t++) {
      biot[i][t] = log(ComputeEmission(phmm, i, Obs[t]));
    }
  } 
*/
#pragma omp parallel num_threads(THREAD_NUM) \
  private(thread_id, nloops, val, maxval, maxvalind, t, j, i, temp, x, y, nonInf) \
  //shared (P, peakPos, phmm, alpha,emission, pprob)
  {
    nloops = 0;
#pragma omp for
  for (k = 0; k < P; k++){
    ++nloops;
    /* 1. Initialization  */     
    for (i = 0; i < phmm->N; i++) {
      delta[peakPos[k]-1][i] = log(phmm->pi[i]) + 
                               gsl_matrix_get(emission_matrix, i, peakPos[k]-1);
      psi[peakPos[k]-1][i] = 0;
      posterior[peakPos[k]-1][i] = alpha[i][peakPos[k]-1] + 
                                   beta[i][peakPos[k]-1] - logprobf[k];
    }
    for (x = 0; x < phmm->M * (phmm->inactive+1); x++){
        posterior[peakPos[k]-1][TFlist[x]] = -100000000000.0;
    }  
    /* 2. Recursion */
    for (t = peakPos[k]; t < peakPos[k+1] - 1; t++) {
      for (j = 0; j < phmm->N; j++) {
        maxval = -VITHUGE;
        maxvalind = 1;
        for (i = 0; i < phmm->N; i++) {
          val = delta[t-1][i] + gsl_matrix_get(phmm->log_A_matrix, i, j);
          if (val > maxval) {
            maxval = val;
            maxvalind = i;
          }
        }
        delta[t][j] = maxval + gsl_matrix_get(emission_matrix, j, t); 
        psi[t][j] = maxvalind;
        
        posterior[t][j] = alpha[j][t] + beta[j][t] - logprobf[k];
        
        for (x = 0; x < phmm->M * (phmm->inactive+1); x++){
          if (j == TFlist[x]){
            if (t <= (peakPos[k] + lengthList[x])){
              posterior[t][j] = -50000000000.0;
            }
            else {
              temp=0.0;
              nonInf = 0;
              for (i = 0; i < lengthList[x]; i++) {
                if ((alpha[j-i][t-i] + beta[j-i][t-i] - logprobf[k]) != -INFINITY){
                  temp += (alpha[j-i][t-i] + beta[j-i][t-i] - logprobf[k]); 
                  nonInf += 1;  
                }
                else{
                 if (k == 10&& i == 0 && t > peakPos[k] + 10 && t < peakPos[k+1] - 10) printf("-INF:%d %d %d %f %f %f\t", t, x, i, alpha[j-i][t-i], beta[j-i][t-i], logprobf[k]);
                  //temp -= logprobf[k]; 
                }
              }
              //if (phmm->M > 1){
                //posterior[t][j] = temp/(phmm->D[x]);
                //posterior[t][j] = temp/nonInf;
                //posterior[t][j] = temp;
              for (y = t - lengthList[x] + 1; y <= t; y++) posterior[y][j] = temp;
              //}
              //else{
                //posterior[t][j] = temp;
              //}
              if (nonInf== 0){
                //posterior[t][j] = -5000000000000.0;
                for (y = t - lengthList[x] + 1; y <= t; y++) posterior[y][j] = -5000000000000.0;
              }
            }
            break;
          }
        }
        //if (posterior[t][j] == -INFINITY) printf("-INF:%f %f %f\t", alpha[j][t], beta[j][t], logprobf[k]); 
        /*
        if (phmm->M == 1 && j == TFlist[0] + phmm->D[0]){
          if (t <= (peakPos[k] + phmm->D[0])){
              posterior[t][j] = -1000000000.0;
          }
          else {
              temp=0.0;
              for (i = 0; i < phmm->D[0]; i++) {
                temp += (alpha[j-i][t-i] + beta[j-i][t-i] - logprobf[k]);   
              }
              posterior[t][j] = temp;
          }    
        }
        */
      }
    }
 
    /* 3. Termination */
    pprob[k] = -VITHUGE;
    q[peakPos[k+1]-2] = 1;
    for (i = 0; i < phmm->N; i++) {
      if (delta[peakPos[k+1]-2][i] > pprob[k]) {
        pprob[k] = delta[peakPos[k+1]-2][i];
        q[peakPos[k+1]-2] = i;
      }
    }
    g[peakPos[k+1]-2] = gamma[q[peakPos[k+1]-2]][peakPos[k+1]-2];
    vprob[peakPos[k+1]-2] = pprob[k]; 
	  /* 4. Path (state sequence) backtracking */

  	for (t = peakPos[k+1] - 3; t >= peakPos[k] - 1; t--){
	  	q[t] = psi[t+1][q[t+1]];
      vprob[t] = delta[t+1][q[t+1]];
      g[t] = gamma[q[t]][t];
    }
  }
  thread_id = omp_get_thread_num();
  }

}
  
void getPosterior_P(FILE *fpIn, FILE *fpOut,int T, int *peakPos, 
                    double **posterior, int indexTF, HMM *phmm)
{
  int *O, *peaks, start, end, TFstart, TFend, length, init, t;
  int old_start = -1;
  int i= -1;
  char chr[4];
  
  while (fscanf(fpIn, "%s\t%d\t%d\t%s\t%d\t%d\t%d", chr, &start, &end, chr, &TFstart, &TFend, &length) != EOF) {
    if (start != old_start){
      i++;
    }
    if (TFstart != -1){
      if (length + 1 > (TFend-TFstart)) {
        init = peakPos[i] - 1;
        t = init + TFstart - start + phmm->D[0] - 2;
        fprintf(fpOut,"%s\t%d\t%d\t%lf\n", chr, TFstart, TFend, posterior[t][indexTF]);
      }
    }
    old_start = start;
	}
	
}

int getPosterior_all_P(FILE *fpIn, FILE *fpOut1, FILE *fpOut2, int T, 
                        int *peakPos, double **posterior,
                        HMM *phmm, int *q, double *vprob)
{
  int *O, *peaks, start, end, TFstart, TFend, length, init, t, j, m, n;
  int old_start = -1;
  //double *X;
  int i= -1;
  int TF, maxTF, indexTF_end, state, pos;
  int half;
  int ifProb;
  double prob;
  char chr[8];
  char chr_[8];
  int *lengthList = ivector(phmm->M * (phmm->inactive+1));
  if (phmm->inactive == 1){
    for (j = 0; j < phmm->M; j++){
      lengthList[j*2] = phmm->D[j];
      lengthList[j*2+1] = phmm->D[j];
    }
  }
  else lengthList = phmm->D;
  fprintf(stdout,"scanning list file and getting posterior\n");
  /*
  if (phmm->M > 1){
    int * stateList = ivector(phmm->M + phmm->extraState);//TO DO: need to change the number 5
    TF = 0;
    for (j = 0; j < phmm->M; j++){
      stateList[j] = TF;
      fprintf(stdout,"%d ", stateList[j]);
      TF += phmm->D[j];
    }
    indexTF_end = stateList[indexTF+1] - 1;
  }
  else {
    int * stateList = ivector(2+8);
    stateList[0] = 0;
    stateList[1] = phmm->D[0];
    fprintf(stdout,"%d %d ", stateList[0], stateList[1]);
    indexTF_end = phmm->D[0] * (indexTF + 1) - 1;
  }
  */
  while (fscanf(fpIn, "%s\t%d\t%d\t%s\t%d\t%d\t%d", chr, &start, &end, chr_, &TFstart, &TFend, &length) != EOF) {
    if (start != old_start) {
      i++;
    }
    if (TFstart != -1){
      if (length + 1 > (TFend-TFstart) && TFend <= end) {
        init = peakPos[i] - 1;
        /*
        prob = -INFINITY;
        TF = -1;
        for (j = 0; j < phmm->M; j++){
          t = init + TFstart - start + phmm->D[j] - 2;
          TF += phmm->D[j];
          if (posterior[t][TF] > prob){
            prob = posterior[t][TF];
            maxTF = j+1;
          }
        }
        */
        state = 100000;
        for (m = init + TFstart - start - 1; m <= init + TFend - start - 1; m++) {
          if (state >= q[m]) {
            state = q[m];
            pos = m;
          }
        }
        ifProb = -1;
        TF = -1;
        ///t = init + TFstart - start - 1;
        half = length / 2;
        t = init + TFstart - start - 1 + half;
        for (j = 0; j < phmm->M * (phmm->inactive+1); j++){
          TF += lengthList[j];
          //if (q[t] <= TF){
          if (state <= TF){
            maxTF = j + 1;
            //prob = posterior[init + TFstart - start + phmm->D[j] - 2][TF];
            //prob = posterior[t + TF - q[t]][TF];
            prob = posterior[pos][state];
            ifProb = 1;
            break;
          }
        }
        if (ifProb == -1){
        /*
        //if (phmm->M == 1 && q[t] > phmm->D[0] - 1 && q[t] <= phmm->D[0] * 2 - 1){
          if (phmm->M == 1 && state > phmm->D[0] - 1 && state <= phmm->D[0] * 2 - 1){
            maxTF = 2;
            //prob = posterior[init + TFstart - start + phmm->D[0] - 2][phmm->D[0] * 2 - 1];
            prob = posterior[pos][state];
          }
          else if (q[t] > TF){
          */
            maxTF = q[t];
          //prob = posterior[init + TFstart - start - 1][maxTF];
            prob = -1000000000.0;
          //}
        }
        fprintf(fpOut2,"%s\t%d\t%d\t%d\t%e\n", chr, TFstart, TFend, maxTF, prob);
        //fprintf(fpOut,"%s\t%d\t%d\t%e\n", chr, TFstart, TFend, posterior[t][indexTF]);
        //fprintf(stdout,"%s\t%d\t%d\t%lf\n", chr, TFstart, TFend, posterior[t][indexTF]);
        //fprintf(fpOut1,"%s\t%d\t%d\t%e\n", chr, TFstart, TFend, posterior[(init + TFstart - start + phmm->D[indexTF] - 2)][indexTF_end]);
        fprintf(fpOut1,"%s\t%d\t%d", chr, TFstart, TFend);
        TF = -1;
        /*
        if (phmm->M == 1) {
          fprintf(fpOut1,"\t%e", posterior[(init + TFstart - start + MIN(phmm->D[0], length+1) - 2)][phmm->D[0] - 1]);
          fprintf(fpOut1,"\t%e", posterior[(init + TFstart - start + MIN(phmm->D[0], length+1) - 2)][phmm->D[0]*2 - 1]);
        
        else {
        */
          for (j = 0; j < phmm->M * (phmm->inactive+1); j++){
            TF += lengthList[j];
            prob = -INFINITY;
            for (m = 0; m <= length+MIN(lengthList[j]-1,end-TFend); m ++) prob = MAX(prob, posterior[(init + TFstart - start + m - 1)][TF]);
            //fprintf(fpOut1,"\t%e", posterior[(init + TFstart - start + MIN(phmm->D[j], length+1) - 2)][TF]);
            fprintf(fpOut1,"\t%e", prob);
          }
        //}
        prob = 0.0;
        for (m = 0; m <= length; m ++) prob += posterior[(init + TFstart - start + m - 1)][phmm->N-4];
        fprintf(fpOut1,"\t%e", prob/length);
        prob = 0.0;
        for (m = 0; m <= length; m ++) prob += posterior[(init + TFstart - start + m - 1)][phmm->N-3];
        fprintf(fpOut1,"\t%e", prob/length);
        fprintf(fpOut1,"\n");
        //fprintf(fpOut2,"%s\t%d\t%d\t%d\t%lf\n", chr, TFstart, TFend, maxTF, prob);
        
        
      }
    }
    //fprintf(stdout,"%d %s %d %d %s %d %d\t", i, chr, start, end, chr, TFstart, TFend);
    //fflush(stdout);
    old_start = start;
	}
	return indexTF_end;
}

int getPosterior_one_P(FILE *fpIn, FILE *fpOut1, FILE *fpOut2, int index, int T,
                       int *peakPos, double **posterior,
                       HMM *phmm, int *q, double *vprob)
{
  int *O, *peaks, start, end, TFstart, TFend, length, init, t, j, m, n;
  int old_start = -1;
  //double *X;
  int i= -1;
  int TF, maxTF, indexTF_end, state, pos;
  int half;
  int ifProb;
  double prob;
  char chr[8];
  char chr_[8];
  int *lengthList = ivector(phmm->M * (phmm->inactive+1));
  lengthList[0] = phmm->D[0];
  lengthList[1] = phmm->D[0]*2;

  for (j = 1; j < phmm->M; j++){
    if (phmm->inactive == 1){
      lengthList[j*2] = lengthList[j*2-1]+phmm->D[j];
      lengthList[j*2+1] = lengthList[j*2]+phmm->D[j];
    }
    else{
      lengthList[j] = lengthList[j-1]+phmm->D[j];
    }
  }

  fprintf(stdout,"scanning list file and getting posterior\n");

  while (fscanf(fpIn, "%s\t%d\t%d\t%s\t%d\t%d\t%d", chr, &start, &end, chr_, &TFstart, &TFend, &length) != EOF) {
    if (start != old_start){
      i++;
    }
    if (TFstart != -1){
      //fprintf(stdout,"%d %s %d %d %lf\n", i, chr, TFstart, TFend, posterior[t][indexTF]);

      if (length + 1 > (TFend-TFstart) && TFend <= end) {
        init = peakPos[i] - 1;
        state = 100000;
        for (m = init + TFstart - start - 1; m <= init + TFend - start - 1; m++) {
          if (state >= q[m]) {
            state = q[m];
            pos = m;
          }
        }
        ifProb = -1;
        TF = -1;
        half = length / 2;
        t = init + TFstart - start - 1 + half;
        for (j = 0; j < phmm->M * (phmm->inactive+1); j++){
          TF += lengthList[j];
          if (state <= TF){
            maxTF = j + 1;
            prob = posterior[pos][state];
            ifProb = 1;
            break;
          }
        }
        if (ifProb == -1){
          maxTF = q[t];
          prob = -1000000000.0;
        }
        fprintf(fpOut2,"%s\t%d\t%d\t%d\t%e\n", chr, TFstart, TFend, maxTF, prob);
        fprintf(fpOut1,"%s\t%d\t%d", chr, TFstart, TFend);
        TF = -1;

        TF = lengthList[(index-1) * (phmm->inactive+1)]-1;
        prob = -INFINITY;
        for (m = 0; m <= length+MIN(phmm->D[index-1]-1,end-TFend); m ++) prob = MAX(prob, posterior[(init + TFstart - start + m - 1)][TF]);
        fprintf(fpOut1,"\t%e", prob);

        TF = lengthList[(index-1) * (phmm->inactive+1)+1]-1;
        prob = -INFINITY;
        for (m = 0; m <= length+MIN(phmm->D[index-1]-1,end-TFend); m ++) prob = MAX(prob, posterior[(init + TFstart - start + m - 1)][TF]);
        fprintf(fpOut1,"\t%e", prob);

        prob = 0.0;
        for (m = 0; m <= length; m ++) prob += posterior[(init + TFstart - start + m - 1)][phmm->N-4];
        fprintf(fpOut1,"\t%e", prob/length);
        prob = 0.0;
        for (m = 0; m <= length; m ++) prob += posterior[(init + TFstart - start + m - 1)][phmm->N-3];
        fprintf(fpOut1,"\t%e", prob/length);
        fprintf(fpOut1,"\n");

      }
    }
    old_start = start;
  }
  return indexTF_end;
}


void getPosterior_labels(FILE *fpIn, FILE *fpOut, int T, int *q,
                         int *peakPos, double **posterior, HMM *phmm)
{
  int *O, *peaks, start, end, TFstart, TFend, length, init, t, i, j;
  int stateStart, stateEnd, dataStart, dataEnd, stateLength;
  int old_start = -1;
  //double *X;
  int TF, maxTF;
  double prob;
  char chr[20];

  fprintf(stdout,"scanning peak file and getting posterior for all positions\n");

  int * stateList = ivector(phmm->N);
  TF = 0;
  for (j = 0; j < phmm->M; j++){
    for (i = TF; i < TF + phmm->D[j]; i++) {
      stateList[i] = j * (phmm->inactive + 1);
    }
    TF += phmm->D[j];
    fprintf(stdout,"%d ", stateList[TF-1]);
    //fflush(stdout);
    if (phmm->inactive == 1){
      for (i = TF; i < TF + phmm->D[j]; i++) {
        stateList[i] = j * (phmm->inactive + 1) + 1;
      }
      TF += phmm->D[j];
      fprintf(stdout,"%d ", stateList[TF-1]);
      //fflush(stdout);
    }
  }
  TF -= 1;
  //fprintf(stdout,"%d ", TF);
  //fflush(stdout);
  for (j = phmm->N - phmm->extraState; j < phmm->N; j++){
    stateList[j] = j;
    fprintf(stdout,"%d ", stateList[j]);
    //fflush(stdout);
  }

  //fprintf(stdout,"check0 ");
  //fflush(stdout);
  i = -1;

  while (fscanf(fpIn, "%s\t%d\t%d", chr, &start, &end) != EOF) {
    //fprintf(stdout,"check1 ");
    //fflush(stdout);
    if (start != old_start){
      i++;
      //fprintf(stdout,"check1 ");
      //fflush(stdout);
    dataStart = peakPos[i] - 1;
    dataEnd = peakPos[i+1] - 2;
    stateStart = 0;
    stateEnd = stateStart;
      t = dataStart;
      if (posterior[t][q[t]] != -INFINITY) {
        prob = posterior[t][q[t]];
        stateLength = 1;
      }
      else {
      }
      //fprintf(stdout,"check2 ");
        prob = 0;
        stateLength = 0;
      //fflush(stdout);
    for (t = dataStart+1; t <= dataEnd; t++){

      if (stateList[q[t]] == stateList[q[t-1]]){
        stateEnd ++;
        if (posterior[t][q[t]] != -INFINITY){
          prob += posterior[t][q[t]];
          stateLength ++;
        }
        maxTF = stateList[q[t]] + 1;
        if (t == dataEnd) fprintf(fpOut,"%s\t%d\t%d\t%d\t%e\t%e\n", chr, start + stateStart, start + stateEnd, maxTF, prob/stateLength, prob/stateLength);
        //fprintf(stdout,"check %d ",t);
        //fflush(stdout);
      }
      else {
        maxTF = stateList[q[t-1]] + 1;
        //if (maxTF <= phmm->M * (phmm->inactive+1)) stateLength *= stateLength;
        //fprintf(fpOut,"%s\t%d\t%d\t%d\t%e\n", chr, start + stateStart, start + stateEnd, maxTF, prob/stateLength);
        if (maxTF <= phmm->M * (phmm->inactive+1)) {
            if (maxTF % 2 == 0) {
                //prob = posterior[t-1][q[t-1]] - posterior[t-1][q[t-1]-phmm->D[maxTF/2-1]];
                fprintf(fpOut, "%s\t%d\t%d\t%d\t%e\t%e\n", chr, start + stateStart,
                        start + stateEnd, maxTF, posterior[t - 1][q[t - 1]],
                        posterior[t - 1][q[t - 1] - phmm->D[maxTF / 2 - 1]]);
            }
            else {
                //prob = posterior[t-1][q[t-1]] - posterior[t-1][q[t-1]+phmm->D[(maxTF-1)/2]];
                fprintf(fpOut, "%s\t%d\t%d\t%d\t%e\t%e\n", chr, start + stateStart,
                        start + stateEnd, maxTF, posterior[t - 1][q[t - 1]],
                        posterior[t - 1][q[t - 1] + phmm->D[(maxTF - 1) / 2]]);

            }
        }
        else fprintf(fpOut,"%s\t%d\t%d\t%d\t%e\t%e\n", chr, start + stateStart, start + stateEnd, maxTF, prob/stateLength, prob/stateLength);
        stateEnd ++;
        stateStart = stateEnd;
        if (posterior[t][q[t]] != -INFINITY) {
          prob = posterior[t][q[t]];
          stateLength = 1;
        }
        else {
          prob = 0;
          stateLength = 0;
        }
        if (t == dataEnd) fprintf(fpOut,"%s\t%d\t%d\t%d\t%e\t%e\n", chr, start + stateStart, start + stateEnd, stateList[q[t]] + 1, prob/stateLength, prob/stateLength);

          //fprintf(stdout,"check %d ",t);
        //fflush(stdout);
      }
    }
    }
    old_start = start;
	}
	
}


