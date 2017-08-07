#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "nrutil.h"
#include "hmm_edit.h"
static char rcsid[] = "$Id: testvit.c,v 1.3 1998/02/23 07:39:07 kanungo Exp kanungo $";

int main (int argc, char **argv)
{
  int 	t, T, range; 
  HMM  	hmm;
  int	*O1;	/* observation sequence O[1..T] */
  double *O2;
  double *GC;
  int P;
  int *peakPos;
  int	*q;	/* state sequence q[1..T] */
  double **delta, *vprob;
  int	**psi;
  double 	proba, logproba; 
  FILE	*fp,*fp2;
  double *S, *g; 
  
  if (argc != 6) {
    printf("Usage error \n");
    printf("Usage: viterbi <model.hmm> <obs.seq> <obs.slop> \n");
    exit (1);
  }
	
  fp = fopen(argv[1], "r");
  if (fp == NULL) {
    fprintf(stderr, "Error: File %s not found\n", argv[1]);
    exit (1);
  }
  ReadOutHMM(fp, &hmm);
  fclose(fp);
  fp = fopen(argv[2], "r");
  if (fp == NULL) {
    fprintf(stderr, "Error: File %s not found\n", argv[2]);
    exit (1);
  }
  ReadSequence(fp, &T, &GC, &O1, &P, &peakPos); 
  fclose(fp);
  fprintf(stdout, "pass1");

  fp = fopen(argv[3], "r");
  if (fp == NULL) {
    fprintf(stderr, "Error: File %s not found\n", argv[3]);
    exit (1);
  }
  ReadSlop(fp, &T, &O2); 
  fclose(fp);
   
  hmm.CG[2] = hmm.CG[3] = GC[2];
  hmm.CG[1] = hmm.CG[4] = GC[1];
  
  fp = fopen(argv[4], "w");
  if (fp == NULL) {
    fprintf(stderr, "Error: File %s not valid \n", argv[4]);
    exit (1);
  }
  PrintHMM(fp, &hmm);
  
  fp2 = fopen(argv[5], "w");
  if (fp2 == NULL) {
    fprintf(stderr, "Error: File %s not valid \n", argv[5]);
    exit (1);
  }
 
  q = ivector(1,T);
  g = dvector(1,T);
  vprob = dvector(1,T);
  delta = dmatrix(1, T, 1, hmm.N);
  psi = imatrix(1, T, 1, hmm.N);

  CalMotifScore(&hmm, &S, O1, P, peakPos);
  
  ViterbiMultiSeq(&hmm, T, O1, O2, S, g, delta, psi, q, vprob, &logproba, P, peakPos); 
  fprintf(fp, "Viterbi  MLE log prob = %E\n", logproba);
  fprintf(fp, "Optimal state sequence:\n");
  fprintf(fp, "T= %d\n", T);
  PrintSequence(fp, T, q);
  fprintf(fp,"------------------------------------\n");
  fprintf(fp,"The two log probabilites and optimal state sequences\n");
  fprintf(fp,"should identical (within numerical precision). \n");
  PrintSequenceProb(fp2, T, q, vprob, g);
  free_ivector(q, 1, T);
  free_ivector(O1, 1, T);
  free_dvector(O2, 1, T);
  free_imatrix(psi, 1, T, 1, hmm.N);
  free_dmatrix(delta, 1, T, 1, hmm.N);
  FreeHMM(&hmm);
}

