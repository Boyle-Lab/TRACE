
#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "nrutil.h"
#include "hmm_edit.h"

static char rcsid[] = "$Id: sequence.c,v 1.2 1998/02/23 06:19:41 kanungo Exp kanungo $";

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
        //fprintf(stdout, "pass1: %d %d %d %d %d\n", *pP, peaks[1],peaks[2],peaks[3],peaks[4]);
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

void PrintSequenceProb(FILE *fp, int T, int *O, double *vprob)
{
        int i;
        fprintf(fp,"%d\t%lf", O[1], vprob[1]);
        for (i=2; i <= T; i++) {
                fprintf(fp,"\n%d\t%lf", O[i], vprob[i]);
        }
}
/*
void ReadBedFile(FILE *fp, int range, int **pos)
{
  int i;
  for (i = 1, i <= range, i ++){
    fscanf(fp, "chr%d\t%d\t%d", pos[i][1], pos[i][2], pos[i][3]);
  
  }
}
*/