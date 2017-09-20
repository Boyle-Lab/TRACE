#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "nrutil.h"
#include "hmm_edit.h"
#include <sys/types.h>
#include <unistd.h>

static char rcsid[] = "$Id: esthmm.c,v 1.1 1998/02/23 07:49:45 kanungo Exp kanungo $";

void Usage(char *name);

int main (int argc, char *argv[])
{
  /*the following are read from input files*/
  FILE	*fp;
  int T; /*total length*/
  double *GC;/*GC content*/
  HMM  	hmm;
  double  *O2; /* slops */
  int *O1; /* sequence, represented by number, 1=A, 2=C, 3=G, 4=T */
           /*right now, there are only two observations: slop and sequence
           might need to change the structure if more data are used in thefuture*/
  int P; /* total number of peaks */
  int *peakPos; /*starting location of each peaks*/
                /*there two are used to seperate values from each peak in concatenated data*/
  int	seed; /* seed for random number generator */
  
  int	niter; /*numbers of iterations, might need to set a max later*/
  
  /*these are calculated variables in BW*/
  double  logprobinit, logprobfinal;
	double  **alpha; 
  double  **beta;
  double  **gamma; 
  
  if (argc != 6 && argc != 5) {
		Usage(argv[0]);
		exit (1);
  }
  /* read the observed sequence */
  fp = fopen(argv[1], "r");
  if (fp == NULL) {
    fprintf(stderr, "Error: File %s not found \n", argv[1]);
    exit (1);
  }
  ReadSequence(fp, &T, &GC, &O1, &P, &peakPos); 
  fclose(fp);
 
  fp = fopen(argv[2], "r");
  if (fp == NULL) {
    fprintf(stderr, "Error: File %s not found \n", argv[2]);
    exit (1);
  }
  ReadSlop(fp, &T, &O2); 
  fclose(fp);
  
	/* initialize the hmm model */
  fp = fopen(argv[3], "r");	
  if (fp == NULL) {
    fprintf(stderr, "Error: File %s not found \n", argv[3]);
    exit (1);
  }
  /*
  if (argc == 6) {
	  ReadHMM(fp, &hmm);
    fclose(fp);
    sscanf(argv[4], "%d", &seed);
    if (&seed == NULL) {
      fprintf(stderr, "Error: seed %s not valid \n", argv[4]);
      fprintf(stdout, "T: %d \n", T);
      //fprintf(stdout, "D: %d \n", hmm.D);
      //fprintf(stdout, "N: %d \n", hmm.N);
      //fprintf(stdout, "GC: %lf \n", GC[1]);
   	  exit (1);
	  }
    InitHMMwithInput(&hmm, seed, GC);
  }
  */
  if (argc == 5) {
    ReadOutHMM(fp, &hmm);
    fclose(fp);
    hmm.bg[2] = hmm.bg[3] = GC[2];
    hmm.bg[1] = hmm.bg[4] = GC[1];
  }

  if (hmm.bg == NULL){
    printf("PROB");
  }
  PrintHMM(stdout, &hmm);
  //fprintf(stdout, "N: %d \n", hmm.N);
  /* allocate memory */
  alpha = dmatrix(1, T, 1, hmm.N);
  beta = dmatrix(1, T, 1, hmm.N);
  gamma = dmatrix(1, T, 1, hmm.N);
  
  fprintf(stdout, "pass");
  /* call Baum Welch */
  BaumWelchBVN(&hmm, T, O1, O2, alpha, beta, gamma, &niter, P, peakPos);
  fprintf(stdout, "pass1");
  /* print the answer */
  PrintHMM(stdout, &hmm);
  fp = fopen(argv[argc-1], "w");
  if (fp == NULL) {
    fprintf(stderr, "Error: File %s not valid \n", argv[argc-1]);
     exit (1);
  }
  PrintHMM(fp, &hmm);
  
  /* free memory */
  free_ivector(O1, 1, T);
  free_dvector(O2, 1, T);
  free_dmatrix(alpha, 1, T, 1, hmm.N);
  free_dmatrix(beta, 1, T, 1, hmm.N);
  free_dmatrix(gamma, 1, T, 1, hmm.N);
  FreeHMM(&hmm);
}

void Usage(char *name)
{
	printf("Usage error. \n");
        
        printf("Usage1: ./esthmm <file.seq> <file.slop> <mod.hmm> <seed> <out.hmm>\n");
        printf("Usage2: ./esthmm <file.seq> <file.slop> <mod.hmm> <out.hmm>\n");
        
        printf("  seed - seed for random number genrator\n");
        printf("  mod.hmm - file with the initial model parameters\n");
        printf("  file.seq - file containing the obs. sequence\n");
        printf("  file.slop - file containing the obs. slop\n");
}
