
#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "nrutil.h"
#include "hmm_edit.h"

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

//using the highest PWM score for each base in + or - strand, but not giving good performance
void CalMotifScore(HMM *phmm, double **S, int *O1, int P, int *peakPos)
{
/* only for odd number right now */
  int k, t, j, i;
  double *X, tempF, tempR, bgF, bgR;
  X = dvector(1,peakPos[P+1]-1);
  for (k = 1; k <= P; k++) {
    
    for (t = peakPos[k]; t < peakPos[k+1]; t++) {
      if (t < (peakPos[k]+phmm->D) || t >= (peakPos[k+1]-1-phmm->D)){
        X[t] = -100;
      }
      else{
        X[t] = log(0.0);
        i = (phmm->D + 1) / 2;
        //for (i = 1; i <= phmm->D; i++){
          tempF = tempR = 0.0;
          bgF = bgR = 0.0;
          for (j = 1; j <= phmm->D; j++){
            tempF += log(phmm->pwm[j][O1[t-i+j]]);
            tempR += log(phmm->pwm[phmm->D-j+1][5-O1[t-i+j]]);
            tempF -= log(phmm->CG[O1[t-i+j]]);
            tempR -= log(phmm->CG[5-O1[t-i+j]]);
          //fprintf(stdout, "X: %lf ", X[t]);
          }
          X[t]=MAX(X[t], tempF);
          X[t]=MAX(X[t], tempR);
          //fprintf(stdout, "X: %lf %lf %lf ", tempF, tempR, X[t]);
        //}
        //X[t] = exp(X[t]) * pow(10,phmm->D);
      }
     
    }
  }
  *S = X;
}
      
    
void ReadSlop(FILE *fp, int *pT, double **pO)
{
        double *O;
        int i;
 
        fscanf(fp, "T= %d\n", pT);
        O = dvector(1,*pT);
        for (i=1; i <= *pT; i++)
                fscanf(fp,"%lf\t", &O[i]);
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
 
