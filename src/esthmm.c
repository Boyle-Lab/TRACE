/*
 *  File: esthmm.c
 *
 *  The main function of training step in TRACE.
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

void Usage(char *name);

int main (int argc, char **argv)
{
  /* The following are read from input files */
  FILE *fp, *fp1, *fp2;
  HMM hmm; /* Initialize the HMM structure */
  int P; /* Total number of peaks */
  int *peakPos; /* Starting location of each peaks*/
                /*there two are used to seperate values from each peak in concatenated data*/
  int T; /* Total length*/
  int *O1; /* Temporarily store the list of sequence, represented by numbers,
            * 1=A, 2=C, 3=G, 4=T */
  double *GC; /* GC content */
  gsl_matrix *pwm_matrix, *obs_matrix; /* Matrix of PWM scores, and observation data
             /* right now, there are only three observations: counts, slop and sequence
              might need to change the structure if more data are used in the future*/
  gsl_vector * slop_vector, *counts_vector;

  int niter; /* Numbers of iterations used in training */
  int i, j, k;
  int c, w = 0, indexTF;
  /* Set default numbers for max iteration and number of threads */
  MAXITERATION = 200;
  THREAD_NUM = 40;
  /* Get all command-line options */
  extern int optind, opterr, optopt;
  extern char *optarg;
  /* Flags for all command line options */
  int oflg=0, mflg=0, nflg=0, aflg =0, pflg =0, fflg=0, eflg=0, tflg=0, errflg=0, vflg=0;
  char *hmminitfile, *slopefile, *countfile, *seqfile;
  char *listfile, *motiffile, peakfile[50], outfile[50];
  char *thresholdfile;
  int ifMulti = 2, peakLength = 2;
  peakfile[0] = '\0';
  outfile[0] = '\0';
  static struct option longopts[] = {
    {"final-model", required_argument, NULL, 'O'},
    {"peak-file", required_argument, NULL, 'P'},
    {"motif-file", required_argument, NULL, 'F'},
    {"thread", required_argument, NULL, 'T'},
    {"threshold-file", required_argument, NULL, 'E'},
    {"max-inter", required_argument, NULL, 'M'},
    {"model", required_argument, NULL, 'N'},
          //{"predictions-file", required_argument, NULL, 'A'},
          //{.name = "seq.file", .has_arg = no_argument, .val = 'Q'},
          //{.name = "count.file", .has_arg = no_argument, .val = 'C'},
          //{.name = "slope.file", .has_arg = no_argument, .val = 'S'},
          //{.name = "init.model", .has_arg = no_argument, .val = 'I'},
    {0,         0,                 0,  0 }
  };
  int option_index = 0;
  while ((c = getopt_long(argc, argv, "vhO:A:M:N:P:T:E:F:", longopts, &option_index)) != EOF){
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
      case 'E':
      /* get the threshold file name*/
        if (eflg)
          errflg++;
        else {
          eflg++;
          thresholdfile = optarg;
        }
        break;
      case 'M':
      /* get the max number of iteration */
        if (mflg)
          errflg++;
        else {
          mflg++;
          MAXITERATION = atoi(optarg);
        }
        break;
      case 'O':
      /* get the HMM output file name*/
        if (oflg)
          errflg++;
        else {
          oflg++;
          strcat(outfile, optarg);
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
      case 'A':
      /* get the output file name */
        if (aflg)
          errflg++;
        else {
          aflg++;
          listfile = optarg;
        }
        break;
      case 'F':
        /* get the peak with motif sites file name */
        if (fflg)
          errflg++;
        else {
          fflg++;
          motiffile = optarg;
        }
        break;
      case 'P':
      /* get the peak file name*/
        if (pflg)
          errflg++;
        else {
          pflg++;
          strcat(peakfile, optarg);
        }
        break;
      case '?':
        errflg++;
    }
  }
  /* Check the required 4 input files were provided */
  if (argc - optind < 3){
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
  hmminitfile = argv[index++]; /* Initial model file */

  /* If trained model file name is not provided, use initial model file name
   * and end with "_final_model.txt" */
  if (!oflg){
    strcat(outfile,hmminitfile);
    strcat(outfile,"_final_model.txt");
  }

  /* Read the observed sequence */
  fp = fopen(seqfile, "r");
  if (fp == NULL) {
    fprintf(stderr, "Error: File %s not found \n", seqfile);
    exit (1);
  }
  GC = dvector(4);
  ReadSequence(fp, &T, GC, &O1, &P, &peakPos);
  fclose(fp);

  /* Read the slope file */
  fp = fopen(slopefile, "r");
  if (fp == NULL) {
    fprintf(stderr, "Error: File %s not found \n", slopefile);
    exit (1);
  }
  slop_vector = gsl_vector_alloc(T);
  ReadTagFile(fp, T, slop_vector, 1.0);
  fclose(fp);

  /* Read the tag counts file */
  fp = fopen(countfile, "r");
  if (fp == NULL) {
    fprintf(stderr, "Error: File %s not found \n", countfile);
    exit (1);
  }
  counts_vector = gsl_vector_alloc(T);
  ReadTagFile(fp, T, counts_vector, 2.0);
  fclose(fp);

  /* Initialize the HMM model */
  /* Read HMM input file */
  fp = fopen(hmminitfile, "r");
  if (fp == NULL) {
    fprintf(stderr, "Error: File %s not found \n", hmminitfile);
    exit (1);
  }
  hmm.model = ifMulti;
  ReadM(fp, &hmm);
  hmm.K = 2 + (hmm.inactive+1) * hmm.M; /* Number of data provided for each state:
                                         tag counts, slop, and PWM score for each TF*/
  //ReadInitHMM(fp, &hmm);
  ReadHMM(fp, &hmm);
  fclose(fp);
  if (GC){
      hmm.bg[0] = hmm.bg[3] = GC[0];
      hmm.bg[2] = hmm.bg[1] = GC[1];
  }
  else{
    hmm.bg[0]=hmm.bg[1]=hmm.bg[2]=hmm.bg[3]=0.25;
  }
  /* Check given file names are valid */
  fp = fopen(outfile, "w");
  if (fp == NULL) {
    fprintf(stderr, "Error: File %s not valid \n", outfile);
    exit (1);
  }
  fclose(fp);
  if (aflg) {
    fp = fopen(listfile, "w");
    if (fp == NULL) {
      fprintf(stderr, "Error: File %s not valid \n", listfile);
      exit(1);
    }
  }
  if (pflg){
    fp = fopen(peakfile, "r");
    if (fp == NULL) {
      fprintf(stderr, "Error: File %s not valid \n", peakfile);
      exit(1);
    }
    fclose(fp);
  }
  if (fflg){
    fp = fopen(motiffile, "r");
    if (fp == NULL) {
      fprintf(stderr, "Error: File %s not valid \n", motiffile);
      exit(1);
    }
    fclose(fp);
  }
  fprintf(stdout, "extra: %d M: %d N: %d T:%d K: %d\n", hmm.extraState,hmm.M,hmm.N,T,hmm.K);

  /* Calculate PWM scores for each motif at each position */
  if (hmm.M > 0){
    pwm_matrix = gsl_matrix_alloc(hmm.M, T);
    CalMotifScore_P(&hmm, pwm_matrix, O1, P, peakPos);
  }
  /* Set the initial mean parameters of PWM score feature based on actual calculation */
  gsl_vector * tmp_vector = gsl_vector_alloc(T);
  index = 0;
  for (i = 0; i < hmm.M; i++) {
    gsl_matrix_get_row(tmp_vector, pwm_matrix, i);
    for (j = 0; j < hmm.N; j++) {
      gsl_matrix_set(hmm.mean_matrix, i, j, gsl_vector_min(tmp_vector) / 2.0);
      gsl_matrix_set(hmm.var_matrix, i, j, 8.0);
    }
    for (j = index; j < index + hmm.D[i] * (hmm.inactive + 1); j++) {
      gsl_matrix_set(hmm.mean_matrix, i, j, gsl_vector_max(tmp_vector) / 2.0);
    }
  }
  /*
  int i, j, index;

    if (ifSet) {
      for (j = index; j < index + hmm.D[i]; j++) {
        gsl_matrix_set(hmm.mean_matrix, hmm.M, j, -0.5/2.0);
        gsl_matrix_set(hmm.mean_matrix, hmm.M+1, j, -0.3);
        //gsl_matrix_set(hmm.var_matrix, hmm.M, j, 2.0);
        gsl_matrix_set(hmm.var_matrix, hmm.M, j, 2.0);
      }
      for (j = index + hmm.D[i]; j < index + hmm.D[i] * 2; j++) {
        gsl_matrix_set(hmm.mean_matrix, hmm.M, j, 0.0);
        gsl_matrix_set(hmm.mean_matrix, hmm.M+1, j, 0.0);
        gsl_matrix_set(hmm.var_matrix, hmm.M, j, 2.0);
      }
    }
    index += hmm.D[i] * (hmm.inactive+1);
  }
  if (ifSet) {
    if (hmm.extraState < 4 * hmm.M) {
      for (i = hmm.N - hmm.extraState; i < hmm.N - hmm.extraState + 3; i++) {
        gsl_matrix_set(hmm.mean_matrix, hmm.M, i, 1.0);
      }
      for (i = hmm.N - hmm.extraState + 3;
           i < hmm.N - hmm.extraState + 6; i++) {
        gsl_matrix_set(hmm.mean_matrix, hmm.M, i, 0.0);
      }
    } else {
      j = 0;
      while (6 * (j+1) < hmm.extraState) {
        for (i = hmm.N - hmm.extraState + 3 * j;
             i < hmm.N - hmm.extraState + 3 * (j + 1); i++) {
          gsl_matrix_set(hmm.mean_matrix, hmm.M, i, 1.0/2.0);
          gsl_matrix_set(hmm.mean_matrix, hmm.M+1, i, 0.7);
        }
        gsl_matrix_set(hmm.mean_matrix, hmm.M, i-1, 1.0);
        gsl_matrix_set(hmm.mean_matrix, hmm.M+1, i-1, 1.0);
        for (i = hmm.N - hmm.extraState + 3 * (j + 1);
             i < hmm.N - hmm.extraState + 3 * (j + 2); i++) {
          gsl_matrix_set(hmm.mean_matrix, hmm.M, i, 0.0);
          gsl_matrix_set(hmm.mean_matrix, hmm.M+1, i, 0.2);
        }
        j++;
      }
      gsl_matrix_set(hmm.mean_matrix, hmm.M, hmm.N - 4, 1.0/2.0);
      gsl_matrix_set(hmm.mean_matrix, hmm.M, hmm.N - 3, -0.1);
      gsl_matrix_set(hmm.mean_matrix, hmm.M+1, hmm.N - 4, -0.3);
      gsl_matrix_set(hmm.mean_matrix, hmm.M+1, hmm.N - 3, -0.2);
    }
  }
  for (i = 0; i < hmm.N; i++){
    covarMatrix_GSL(&hmm, i, hmm.cov_matrix[i]);
  }
  free_ivector(O1, T);
*/
  /* Put PWM scores, counts and/or slope into a matrix */
  obs_matrix = gsl_matrix_alloc(hmm.K, T);
  for (i = 0; i < hmm.M; i++){
    gsl_matrix_get_row(tmp_vector, pwm_matrix, i);
    gsl_matrix_set_row(obs_matrix, i, tmp_vector);
  }
  gsl_matrix_set_row(obs_matrix, hmm.M, slop_vector);
  gsl_vector_free(slop_vector);
  gsl_matrix_set_row(obs_matrix, hmm.M+1, counts_vector);
  gsl_vector_free(counts_vector);
  if (hmm.M > 0){
    gsl_matrix_free(pwm_matrix);
    gsl_vector_free(tmp_vector);
  }

  hmm.thresholds = (double *) dvector(hmm.M);
  if (eflg) {
      fp = fopen(thresholdfile, "r");  //TODO: thresholds function
    printf("%s\n", thresholdfile);
    if (fp == NULL) {
      fprintf(stderr, "Error: File %s not valid \n", thresholdfile);
      exit(1);
    }
    for (i = 0; i < hmm.M; i++) {
      if(fscanf(fp, "%lf\n", &(hmm.thresholds[i])) == EOF){
        fprintf(stderr, "Error: threshold file error \n");
        exit (1);
      }
      printf("%d: %lf\n", i, hmm.thresholds[i]);
    }
    fclose(fp);
    printf("thresholds\n");
  }
  else{
    for (i = 0; i < hmm.M; i++) {
      hmm.thresholds[i] = -INFINITY;
    }
  }

  /* Start training step */
  double **alpha = dmatrix(hmm.N, T);
  double **beta = dmatrix(hmm.N, T);
  double **gamma = dmatrix(hmm.N, T);
  double *logprobf = dvector(P); /*vector containing log likelihood for each peak*/
  gsl_matrix * emission_matrix = gsl_matrix_alloc(hmm.N, T);
  BaumWelch(&hmm, T, obs_matrix, &niter, P, peakPos, logprobf, alpha, beta, gamma, emission_matrix);
  gsl_matrix_free(obs_matrix);
  fp = fopen(outfile, "w");
  /* Print the final model */
  PrintHMM(fp, &hmm);
  fclose(fp);

  if (pflg){
    /* Start decoding step */
    int *q = ivector(T); /* state sequence q[1..T] */
    int **psi = imatrix(T, hmm.N);
    double *g = dvector(T);
    double *vprob = dvector(T);
    double **delta = dmatrix(T, hmm.N);
    double *logproba = dvector(P);
    double  **posterior = dmatrix(T, hmm.N);
    Viterbi(&hmm, T, g, alpha, beta, gamma, logprobf, delta, psi, q,
            vprob, logproba, posterior, P, peakPos, emission_matrix);
    int TF_end;
    if (fflg){
      fp = fopen(motiffile, "r");
      strcat(motiffile,"_with_probs.txt");
      fp1 = fopen(motiffile, "w");
      if (fp1 == NULL) {
        fprintf(stderr, "Error: File %s not valid \n", motiffile);
        exit (1);
      }
      TF_end = getPosterior_all_P(fp, fp1, T, peakPos, posterior, &hmm, q, vprob);
      fclose(fp);
      fclose(fp1);
    }

    fp = fopen(peakfile, "r");
    if (fp == NULL) {
      fprintf(stderr, "Error: File %s not valid \n", peakfile);
      exit (1);
    }
    strcat(outfile,"_viterbi_results.txt");
    fp1 = fopen(outfile, "w");
    getPosterior_labels(fp, fp1, T, q, peakPos, posterior, &hmm);
    fclose(fp);
    fclose(fp1);
  }
}


void Usage(char *name)
{
  printf("Usage error. \n");
  printf("Usage1: %s [-v] ./esthmm <seq.file> <counts.file> <slope.file> "
         "<init.model.file> --final-model <final.model.file> -peak-file <peak_3.file> "
         "--thread <thread.num>\n", name);
  printf("Usage1: %s [-v] ./esthmm <sequence.file> <slope.file> <counts.file> "
         "<init.model.file> --final-model <final.model.file> -peak-file <peak_3.file> "
         "--motif-file <peak_7.file>  -T <thread.num>\n", name);
  printf("  seq.file - file containing the obs. sequence\n");
  printf("  count.file - file containing the obs. tag counts\n");
  printf("  slope.file - file containing the obs. slope\n");
  printf("  init.model.file - file with the initial model parameters\n");
  printf("  peak.file - file containing regions to detect TFBSs\n");
  printf("  final.model.file - output file containing the learned HMM\n");
}
