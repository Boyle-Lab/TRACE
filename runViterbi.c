#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "nrutil.h"
#include "hmm_edit.h"
static char rcsid[] = "$Id: testvit.c,v 1.3 1998/02/23 07:39:07 kanungo Exp kanungo $";

int main (int argc, char **argv)
{
	int 	t, T; 
	HMM  	hmm;
	int	*O1;	/* observation sequence O[1..T] */
  double *O2;
  double *GC;
	int	*q;	/* state sequence q[1..T] */
	double **delta;
	int	**psi;
	double 	proba, logproba; 
	FILE	*fp;

	if (argc != 4) {
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
	ReadSequence(fp, &T, &GC, &O1); 
	fclose(fp);
  fprintf(stdout, "pass1");
  
  
  fp = fopen(argv[3], "r");
	if (fp == NULL) {
		fprintf(stderr, "Error: File %s not found\n", argv[3]);
		exit (1);
	}
  ReadSlop(fp, &T, &O2); 
  fclose(fp);
  
  fprintf(stdout, "pass2 %d %d %lf ", T, O1[T], O2[T]);
  
  hmm.CG[2] = hmm.CG[3] = GC[2];
  hmm.CG[1] = hmm.CG[4] = GC[1];
  PrintHMM(stdout, &hmm);
  
	q = ivector(1,T);

	delta = dmatrix(1, T, 1, hmm.N);
	psi = imatrix(1, T, 1, hmm.N);

/*
	printf("------------------------------------\n");
	printf("Viterbi using direct probabilities\n");
	Viterbi(&hmm, T, O1, O2, delta, psi, q, &proba); 
	fprintf(stdout, "Viterbi  MLE log prob = %E\n", log(proba));
	fprintf(stdout, "Optimal state sequence:\n");
	PrintSequence(stdout, T, q);

	printf("------------------------------------\n");
*/ 
	//printf("Viterbi using log probabilities\n");
	/* note: ViterbiLog() returns back with log(A[i][j]) instead
	** of leaving the A matrix alone. If you need the original A,
	** you can make a copy of hmm by calling CopyHMM */ 

	ViterbiLog(&hmm, T, O1, O2, delta, psi, q, &logproba); 

	fprintf(stdout, "Viterbi  MLE log prob = %E\n", logproba);
	fprintf(stdout, "Optimal state sequence:\n");
	PrintSequence(stdout, T, q);
	printf("------------------------------------\n");
	printf("The two log probabilites and optimal state sequences\n");
	printf("should identical (within numerical precision). \n");
	
	free_ivector(q, 1, T);
	free_ivector(O1, 1, T);
  free_dvector(O2, 1, T);
	free_imatrix(psi, 1, T, 1, hmm.N);
	free_dmatrix(delta, 1, T, 1, hmm.N);
	FreeHMM(&hmm);
}

