
#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "nrutil.h"
#include "hmm_edit.h"
#include "logmath.h"
#include "const.h"

static char rcsid[] = "$Id: sequence.c,v 1.2 1998/02/23 06:19:41 kanungo Exp kanungo $";

void ReadSequence(FILE *fp, int *pT, double **GC, int **pO, int *pP, int **peakPos)
{
        int *O, *peaks;
        double *X;
        int i;
 
        fscanf(fp, "T= %d\n", pT);
        X = dvector(1,4);
        fscanf(fp, "GC: ");
        for (i = 1; i <= 4; i++) {
	  fscanf(fp, "%lf\t", &(X[i])); 
	}
	fscanf(fp,"\n");
        *GC = X;
        O = ivector(1,*pT);
        for (i=1; i <= *pT; i++) {
          fscanf(fp,"%d", &O[i]);
        }
        fscanf(fp,"\n");
        *pO = O;
        fscanf(fp, "P= %d\n", pP);
        peaks = ivector(1, *pP + 1);
        for (i=1; i <= *pP + 1; i++){
          fscanf(fp,"%d", &peaks[i]);
        }
        *peakPos = peaks;
        
}

/*using the highest PWM score for sequence centere at each base in + or - strand, but not giving good performance*/
/* only for motif length of odd number right now, need to change later */
void CalMotifScore(HMM *phmm, double **S, int *O1, int P, int *peakPos)
{

  int k, t, j, i, n;
  double *X, tempF, tempR, bgF, bgR;
  //X = dvector(1,peakPos[P+1]-1);/*temp store score vector*/
  for (n = 1; n < phmm->K; n++) {
    for (k = 1; k <= P; k++) {
      for (t = peakPos[k]; t < peakPos[k+1]; t++) {
        if (t < (peakPos[k]+phmm->D[n]) || t >= (peakPos[k+1]-1-phmm->D[n])){
          S[n][t] = -100;
        }
        else{
          S[n][t] = log(0.0);
          i = (phmm->D[n] + 1) / 2;
          //for (i = 1; i <= phmm->D; i++){
          tempF = tempR = 0.0;
          bgF = bgR = 0.0;
          for (j = 1; j <= phmm->D[n]; j++){
            tempF += log(phmm->pwm[n][j][O1[t-i+j]]);
            tempR += log(phmm->pwm[n][phmm->D[n]-j+1][5-O1[t-i+j]]);
            tempF -= log(phmm->bg[O1[t-i+j]]);
            tempR -= log(phmm->bg[5-O1[t-i+j]]);
          //fprintf(stdout, "X: %lf ", X[t]);
          }
          S[n][t]=MAX(S[n][t], tempF);
          S[n][t]=MAX(S[n][t], tempR);
          //fprintf(stdout, "X: %lf %lf %lf ", tempF, tempR, X[t]);
        //}
        //X[t] = exp(X[t]) * pow(10,phmm->D);
        }
      }
     
    }
    //S[n] = X;
  }
  
}
      
/*calculate emission probability for data vector in state j*/
double ComputeEmission(HMM *phmm, int j, double *data) 
{
  
  double pNorm;    /* multivariance normal distribution */ 
  
  double **covar = dmatrix(1, phmm->K, 1, phmm->K);
  covarMatrix(phmm, j, covar);
  //printMatrix(covar, phmm->K, phmm->K);
  pNorm = MultiVarNormDist(phmm->mu, j, covar, phmm->K, data);

  /*this was used for test case of bivariance
  pNorm = (1/((pow(SQRT_TWO_PI, 2.0) * phmm->sigma[j][1] * phmm->sigma[j][2]) * pow((1 - pow(phmm->rho[j], 2.0)), 0.5))) * 
	  exp( (-0.5) / (1 - pow(phmm->rho[j], 2.0))*(pow((S[t] - phmm->mu[j][1]), 2.0)/pow(phmm->sigma[j][1],2.0)+
	  pow((O2[t] - phmm->mu[j][2]), 2.0)/pow(phmm->sigma[j][2],2.0)-2*phmm->rho[j]*(S[t] - phmm->mu[j][1])*(O2[t] - 
	  phmm->mu[j][2])/(phmm->sigma[j][1]*phmm->sigma[j][2]))); */
  //fprintf(stdout, "pNorm: %lf \t", pNorm);  
  return pNorm;
}

void covarMatrix(HMM *phmm, int state, double **matrix)
{
  int i, j, n;
  double corr;
  //for (n = 1; n <= phmm->N; n++){
    for (i = 1; i <= phmm->K; i++) {
      for (j = 1; j <= phmm->K; j++) {
        if (i == j) {
          matrix[i][j] = phmm->sigma[i][state] * phmm->sigma[j][state];
        }
        else if (i < j) {
          corr = phmm->rho[phmm->K*(i-1)+j-i*(i+1)/2][state];
          matrix[i][j] = phmm->sigma[i][state] * phmm->sigma[j][state] * corr;
        }
        else {
          //corr = phmm->rho[phmm->K*(j-1)+i-j*(j+1)/2][state];
          matrix[i][j] = matrix[j][i];//phmm->sigma[i][state] * phmm->sigma[j][state] * corr;
        }
      }
    }
  //}
}
          
          
/*
double ComputeEmission(HMM *phmm, int j, int t, int *O1, double *O2) {
  double pNorm;    // P(O2) ~ N(mu, gamma) 
  int i;
   
  
  
  if (j <= phmm->D){ //first D states are TF motif
    pNorm = (1/(SQRT_TWO_PI * phmm->sigma[1])) * exp( (-0.5) * (1/pow(phmm->sigma[1],2.0))*pow((O2[t] - phmm->mu[1]), 2.0)); 
    return (pNorm * (phmm->pwm[j][O1[t]]));
  } 
  
  if (j > (phmm->D)){
    pNorm = (1/(SQRT_TWO_PI * (phmm->sigma[j+1-(phmm->D)]))) * exp((-0.5) * (1/pow(phmm->sigma[j+1-(phmm->D)],2.0))*pow((O2[t] - phmm->mu[j+1-(phmm->D)]), 2.0)); 
    return (pNorm * (phmm->CG[O1[t]]));
  }
  
}
*/

    
void ReadSlop(FILE *fp, int *pT, double **pO)
{
  double *O;
  int i;
  fscanf(fp, "T= %d\n", pT);
  O = dvector(1,*pT);
  for (i=1; i <= *pT; i++) {
    fscanf(fp,"%lf\t", &O[i]);
  }
  *pO = O;
}

 
void PrintSequence(FILE *fp, int T, int *O)
{
        int i;
        fprintf(fp,"%d", O[1]);
        for (i=2; i <= T; i++) 
                fprintf(fp,"\t%d", O[i]);
	printf("\n");
 
}

void PrintSequenceProb(FILE *fp, int T, int *O, double *vprob, double *g)
{
  int i;
  fprintf(fp,"%d\t%lf\t%lf", O[1], vprob[1], g[1]);
  for (i=2; i <= T; i++) {
    fprintf(fp,"\n%d\t%lf\t%lf", O[i], vprob[i], g[i]);
  }
}

/*
void GenSequenceArray(HMM *phmm, int seed, int T, int *O, int *q)
{
        int     t = 1;
        int     q_t, o_t;

	hmmsetseed(seed); 
 
        q[1] = GenInitalState(phmm);
        O[1] = GenSymbol(phmm, q[1]);
 
        for (t = 2; t <= T; t++) {
                q[t] = GenNextState(phmm, q[t-1]);
                O[t] = GenSymbol(phmm, q[t]);
        }
}

int GenInitalState(HMM *phmm)
{
        double val, accum;
        int i, q_t;
 
        val = hmmgetrand();
        accum = 0.0;
        q_t = phmm->N;
        for (i = 1; i <= phmm->N; i++) {
                if (val < phmm->pi[i] + accum) {
                        q_t = i;
                        break;
                }
                else {
                                accum += phmm->pi[i];
                }
        }
 
        return q_t;
}

int GenNextState(HMM *phmm, int q_t)
{
        double val, accum;
        int j, q_next;
 
        val = hmmgetrand();
        accum = 0.0;
        q_next = phmm->N;
        for (j = 1; j <= phmm->N; j++) {
                if ( val < phmm->A[q_t][j] + accum ) {
                        q_next = j;
                        break;
                }
                else
                        accum += phmm->A[q_t][j];
        }
 
        return q_next;
}

//XXX:change it
int GenSymbol(HMM *phmm, int q_t)
{
        double val, accum;
        int j, o_t;
 
        val = hmmgetrand();
        accum = 0.0;
        o_t = phmm->M;
        for (j = 1; j <= phmm->M; j++) {
                if ( val < phmm->pwm[q_t][j] + accum ) {
                       o_t = j;
                       break;
                }
                else
                        accum += phmm->pwm[q_t][j];
        }
 
        return o_t;
}
 
*/