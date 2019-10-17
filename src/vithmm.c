/*
 *  File: vithmm.c
 *
 *  The main function of viterbi step in TRACE.
 *
 *  The HMM structure and some codes are borrowed and modified from Kanungo's
 *  original HMM program.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "nrutil.h"
#include "hmm.h"
#include <sys/types.h>
#include <unistd.h>
#include <time.h>
#include <getopt.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

void Usage(char *name);

int main (int argc, char **argv)
{
  /* The following are read from input files */
  FILE *fp, *fp1, *fp2;
  HMM hmm; /* Trained HMM structure */
  int *O1; /* Observation sequence of PWM scores O[1..T] */
  double *O2, *O3; /* Observation sequence of counts and slopes O[1..T] */
  double *GC; /* GC content */
  int T; /* Total length*/
  int P; /* Total number of peaks */
  int *peakPos; /* Starting location of each peaks*/
  /* some numbers calculated during decoding */
  double **delta, *vprob;
  double proba, *logproba;
  double *logprobf;
  double *g; 
  int **pos;
  int *q;	/* State sequence q[1..T] */
  int **psi;
  gsl_vector * slop_vector, *counts_vector;
  gsl_matrix *pwm_matrix, *obs_matrix;
  int t, range, i, j, indexTF = 0;
  /* Set default numbers for max iteration and number of threads */
  THREAD_NUM = 40;
  /* Get all command-line options */
  extern char *optarg;
  extern int optind, opterr, optopt;
  /* Flags for all command line options */
  int errflg=0, vflg=0, fflg=0, tflg=0, pflg=0, nflg=0, dflg=0, aflg=0, bflg=0, eflg=0, xflg=0;
  char hmmfile[100], *slopefile, *countfile, *seqfile, *thresholdfile;
  char peakfile[100], motiffile[100], predfile[100], scorefile[100];
  int c;
  int ifMulti = 2; /* defalult model is multivaraint normal distribution
                      with problematic hidden state dropped */
  hmmfile[0] = '\0';
  peakfile[0] = '\0';
  scorefile[0] = '\0';
  motiffile[0] = '\0';
  predfile[0] = '\0';
  static struct option longopts[] = {
          {"peak-file", required_argument, NULL, 'P'},
          {"motif-file", required_argument, NULL, 'F'},
          {"scores-file", required_argument, NULL, 'A'},
          {"predictions-file", required_argument, NULL, 'B'},
          {"thread", required_argument, NULL, 'T'},
          {"threshold-file", required_argument, NULL, 'E'},
          {"model", required_argument, NULL, 'N'},
          {"TF-index", required_argument, NULL, 'X'},
          //{.name = "seq.file", .has_arg = no_argument, .val = 'Q'},
          //{.name = "count.file", .has_arg = no_argument, .val = 'C'},
          //{.name = "slope.file", .has_arg = no_argument, .val = 'S'},
          //{.name = "init.model", .has_arg = no_argument, .val = 'I'},
          {0,         0,                 0,  0 }
  };
  int option_index = 0;
  while ((c = getopt_long(argc, argv, "vhO:A:B:M:N:P:T:E:F:", longopts, &option_index)) != EOF){
    switch (c) {
	case 'v':
	  vflg++;
	  break;
	case 'h':
	  Usage(argv[0]);
	  exit(1);
	  break;
    case 'T':
      /* set number of threads */
      if (tflg)
        errflg++;
      else {
        tflg++;
        THREAD_NUM = atoi(optarg);
      }
        break;
    case 'X':
      /* set the index of TF */
      if (xflg)
        errflg++;
      else {
        xflg++;
        indexTF = atoi(optarg);
      }
      break;
    case 'E':
      /* get the threshold file name*/
      if (eflg)
        errflg++;
      else {
        eflg++;
        thresholdfile = optarg;
      }
      break;
    case 'A':
    /* get the output file name */
      if (aflg)
        errflg++;
      else {
        aflg++;
        strcat(scorefile, optarg);
      }
      break;
    case 'B':
    /* get the output file name */
      if (bflg)
        errflg++;
      else {
        bflg++;
        strcat(predfile, optarg);
      }
      break;
    case 'P':  
    /* get the intersect input file name*/
      if (pflg) 
        errflg++; 
      else { 
        pflg++;
        strcat(peakfile, optarg);
      } 
      break;
    case 'F':
      /* get the peak with motif sites file name */
      if (fflg)
        errflg++;
      else {
        fflg++;
        strcat(motiffile, optarg);
      }
      break;
    case 'N':  
    /* choose independ normal or multivariance, default is independent(0)*/
      if (nflg) 
        errflg++; 
      else { 
        nflg++;  
        ifMulti = atoi(optarg);
      } 
      break;
    case '?':
      errflg++;
    }
  }

  /* Check the required 4 input files were provided */
  if (argc - optind < 3){
    fprintf(stderr, "Error: required files not provided \n", seqfile);
    errflg++;
  }
  if (fflg + pflg < 1){
    fprintf(stderr, "Error: peak file or motif file is required \n", seqfile);
    errflg++;
  }
  if (errflg) {
    Usage(argv[0]);
    exit (1);
  }

  /* Get required input files */
  int index = optind;
  seqfile = argv[index++]; /* Sequence input file */
  countfile = argv[index++]; /* Counts file */
  slopefile = argv[index++]; /* Slopes file */
  strcat(hmmfile, argv[index++]); /* final model file */

  /* Read the observed sequence */
  checkFile(seqfile, "r");
  fp = fopen(seqfile, "r");
  GC = dvector(4);
  ReadSequence(fp, &T, GC, &O1, &P, &peakPos);
  fclose(fp);

  /* Read the slop file */
  checkFile(slopefile, "r");
  fp = fopen(slopefile, "r");
  slop_vector = gsl_vector_alloc(T);
  ReadTagFile(fp, T, slop_vector, 1.0);
  fclose(fp);

  /* Read the tag counts file */
  checkFile(countfile, "r");
  fp = fopen(countfile, "r");
  counts_vector = gsl_vector_alloc(T);
  ReadTagFile(fp, T, counts_vector, 1.0);
  fclose(fp);

  /* Get the learned HMM model */
  /* Read HMM input file */
  checkFile(hmmfile, "r");
  fp = fopen(hmmfile, "r");
  hmm.model = ifMulti;
  ReadM(fp, &hmm);
  ReadHMM(fp, &hmm);
  fclose(fp);
  if (GC){
    hmm.bg[0] = hmm.bg[3] = GC[0];
    hmm.bg[2] = hmm.bg[1] = GC[1];
  }
  else{
    hmm.bg[0]=hmm.bg[1]=hmm.bg[2]=hmm.bg[3]=0.25;
  }

  fprintf(stdout, "M: %d T:%d K: %d\n", hmm.M,T,hmm.K);

  /* Calculate PWM scores for each motif at each position */
  pwm_matrix = gsl_matrix_alloc(hmm.M, T);
  CalMotifScore_P(&hmm, pwm_matrix, O1, P, peakPos);
  free_ivector(O1, T);
  /* Put PWM scores, counts and/or slope into a matrix */
  obs_matrix = gsl_matrix_alloc(hmm.K, T);
  gsl_vector * tmp_vector = gsl_vector_alloc(T);
  for (i = 0; i < hmm.M; i++){
    gsl_matrix_get_row(tmp_vector, pwm_matrix, i);
    gsl_matrix_set_row(obs_matrix, i, tmp_vector);
  }
  gsl_matrix_free(pwm_matrix);
  gsl_vector_free(tmp_vector);

  gsl_matrix_set_row(obs_matrix, hmm.M, slop_vector);
  gsl_vector_free(slop_vector);

  gsl_matrix_set_row(obs_matrix, hmm.M+1, counts_vector);
  gsl_vector_free(counts_vector);
  
  hmm.thresholds = (double *) dvector(hmm.M);
  /* with a given threshold file, TRACE can limit the state labeling to the
     regions with PWM score higher than the provided value */
  if (eflg) {
    fp = fopen(thresholdfile, "r");  //TODO: thresholds function
    for (i = 0; i < hmm.M; i++) {
      if(fscanf(fp, "%lf\n", &(hmm.thresholds[i])) == EOF){
        fprintf(stderr, "Error: threshold file error \n");
        exit (1);
      }
    }
    fclose(fp);
  }
  else{
    for (i = 0; i < hmm.M; i++) {
      hmm.thresholds[i] = -INFINITY;
    }
  }

  /* Calculate emission probabilities */
  gsl_matrix * emission_matrix = gsl_matrix_alloc(hmm.N, T);
  if (hmm.model == 0) EmissionMatrix(&hmm, obs_matrix, P, peakPos, emission_matrix, T);
  if (hmm.model == 1) EmissionMatrix_mv(&hmm, obs_matrix, P, peakPos, emission_matrix, T);
  if (hmm.model == 2) EmissionMatrix_mv_reduce(&hmm, obs_matrix, P, peakPos, emission_matrix, T);

  double **alpha = dmatrix(hmm.N, T);
  double **beta = dmatrix(hmm.N, T);
  double **gamma = dmatrix(hmm.N, T);
  double  **posterior = dmatrix(T, hmm.N);
  logprobf = dvector(P);
  Forward_P(&hmm, T, alpha, logprobf, P, peakPos, emission_matrix);
  Backward_P(&hmm, T, beta, P, peakPos, emission_matrix);
  ComputeGamma_P(&hmm, T, alpha, beta, gamma);
  q = ivector(T);
  g = dvector(T);
  vprob = dvector(T);
  delta = dmatrix(T, hmm.N);
  psi = imatrix(T, hmm.N);
  logproba = dvector(P);
 
  /*                    */
  /* Start Viterbi step */
  /*                    */
  Viterbi(&hmm, T, g, alpha, beta, gamma, logprobf, delta, psi, q,
          vprob, logproba, posterior, P, peakPos, emission_matrix);
  gsl_matrix_free(emission_matrix);
  free_dmatrix(delta, T, hmm.N);
  free_dmatrix(gamma, hmm.N, T);
  free_dmatrix(alpha, hmm.N, T);
  free_dmatrix(beta, hmm.N, T);

  int TF_end;
  if (fflg){
    fp = fopen(motiffile, "r");
    if (fp == NULL) {
      fprintf(stderr, "Error: File path %s not valid \n", motiffile);
      exit (1);
    }
    if (aflg) fp1 = fopen(scorefile, "w");
    else {
      strcat(scorefile, motiffile);
      strcat(scorefile, "_with_probs.txt");
      fp1 = fopen(scorefile, "w");
    }
    if (fp1 == NULL) {
      fprintf(stderr, "Error: File path %s not valid \n", motiffile);
      exit (1);
    }
    printf("%s: %s\n", scorefile, motiffile);
    TF_end = getPosterior_motif(fp, fp1, T, peakPos, posterior, &hmm, q, vprob);
    fclose(fp);
    fclose(fp1);
  }
  if (pflg){
    fp = fopen(peakfile, "r");
    if (fp == NULL) {
      fprintf(stderr, "Error: File path %s not valid \n", peakfile);
      exit (1);
    }
    if (bflg) fp1 = fopen(predfile, "w");
    else{
      strcat(predfile, hmmfile);
      strcat(predfile, "_viterbi_results.txt");
      fp1 = fopen(predfile, "w");
    }
    getPosterior_all(fp, fp1, T, q, peakPos, posterior, &hmm);
    fclose(fp);
    fclose(fp1);
  }
  free_ivector(q, T);
  free_imatrix(psi, T, hmm.N);
  FreeHMM(&hmm);
}

void Usage(char *name)
{
  printf("Usage error. \n");
  printf("Usage: %s [-v] <seq.file> <slope.file> <counts.file> <final.model.file> "
         "-peak-file <peak_3.file> --motif-file <peak_7.file> "
         "-T <thread.num>\n", name);
  printf("  seq.file - file containing the obs. sequence\n");
  printf("  count.file - file containing the obs. tag counts\n");
  printf("  slope.file - file containing the obs. slope\n");
  printf("  final.model.file - output file containing the learned HMM\n");
  printf("  peak.file - file containing regions to detect TFBSs\n");
}
