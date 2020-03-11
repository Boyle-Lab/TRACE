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
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort.h>

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
              might need to change the structure if more data are used in the future */
  gsl_vector *slop_vector, *counts_vector;

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
  int oflg=0, mflg=0, nflg=0, aflg =0, bflg =0, pflg =0, fflg=0, eflg=0, tflg=0, lflg=0, iflg=0;
  int errflg=0, vflg=0;
  char *slopefile, *countfile, *seqfile, *thresholdfile, *listfile;
  char hmminitfile[150], motiffile[150], peakfile[150], outfile[150], predfile[150], scorefile[150];
  int ifMulti = 2; /* defalult model is multivaraint normal distribution
                      with problematic hidden state dropped */
  int ifSkip = 0; /* if skip training step and start viterbi step
                      default is no */
  peakfile[0] = '\0';
  outfile[0] = '\0';
  predfile[0] = '\0';
  scorefile[0] = '\0';
  motiffile[0] = '\0';
  static struct option longopts[] = {
    {"final-model", required_argument, NULL, 'O'},
    {"initial-model", required_argument, NULL, 'I'},
    {"peak-file", required_argument, NULL, 'P'},
    {"motif-file", required_argument, NULL, 'F'},
    {"thread", required_argument, NULL, 'T'},
    {"threshold-file", required_argument, NULL, 'E'},
    {"max-inter", required_argument, NULL, 'M'},
    {"model", required_argument, NULL, 'N'},
    {"predictions-file", required_argument, NULL, 'B'},
    {"scores-file", required_argument, NULL, 'A'},
    {"viterbi", no_argument, NULL, 'V'},
          //{.name = "seq.file", .has_arg = no_argument, .val = 'Q'},
          //{.name = "count.file", .has_arg = no_argument, .val = 'C'},
          //{.name = "slope.file", .has_arg = no_argument, .val = 'S'},
          //{.name = "init.model", .has_arg = no_argument, .val = 'I'},
    {0,         0,                 0,  0 }
  };
  int option_index = 0;
  while ((c = getopt_long(argc, argv, "vhO:A:B:M:N:P:T:E:F:V:I:", longopts, &option_index)) != EOF){
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
      case 'I':
      /* get the initial HMM file name*/
        if (iflg)
          errflg++;
        else {
          iflg++;
          strcat(hmminitfile, optarg);
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
      case 'F':
        /* get the peak with motif sites file name */
        if (fflg)
          errflg++;
        else {
          fflg++;
          strcat(motiffile, optarg);
        }
        break;
      case 'P':
      /* get the peak file name */
        if (pflg)
          errflg++;
        else {
          pflg++;
          strcat(peakfile, optarg);
        }
        break;
      case 'V':
      /* if skip training step */
        if (lflg)
          errflg++;
        else {
          lflg++;
          ifSkip = 1;
        }
        break;
      case '?':
        errflg++;
    }
  }
  /* Check the required 3 input files were provided */
  if (argc - optind < 2){
    fprintf(stderr, "Error: required files not provided \n", seqfile);
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
  //hmminitfile = argv[index++]; /* Initial model file */

  /* Read the observed sequence */
  checkFile(seqfile, "r");
  fp = fopen(seqfile, "r");
  GC = dvector(4);
  ReadSequence(fp, &T, GC, &O1, &P, &peakPos);
  fclose(fp);

  /* Read the slope file */
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
  
  /* Check if user only wants to run decoding step*/
  if (ifSkip){
    /* Read HMM input file */
    if (!oflg){
      fprintf(stderr, "Error: final model file required \n");
      exit (1);
    }
    checkFile(outfile, "r");
    fp = fopen(outfile, "r");
    hmm.model = ifMulti;
    ReadM(fp, &hmm);
    ReadHMM(fp, &hmm);
    fclose(fp);
    /* Check region of interest file is provided for decoding */
    if (!pflg && !fflg){
      fprintf(stderr, "Error: peak file required \n");
      exit (1);
    }
  }
  else{
    /* Initialize the HMM model */
    /* Read HMM input file */
    if (!iflg){
      fprintf(stderr, "Error: initial model file required \n");
      exit (1);
    }
    checkFile(hmminitfile, "r");
    fp = fopen(hmminitfile, "r");
    hmm.model = ifMulti;
    ReadM(fp, &hmm);
    if (hmm.M == 0) ReadInitHMM(fp, &hmm);
    else ReadHMM(fp, &hmm);
    fclose(fp);
    /* Check given file names are valid */
    /* If trained model file name is not provided, use initial model file name
       with suffix "_final_model.txt" */
    if (!oflg){
      strcat(outfile,hmminitfile);
      strcat(outfile,"_final_model.txt");
    }
    checkFile(outfile, "w");
  }
  if (GC){
    hmm.bg[0] = hmm.bg[3] = GC[0];
    hmm.bg[2] = hmm.bg[1] = GC[1];
  }
  else{
    hmm.bg[0]=hmm.bg[1]=hmm.bg[2]=hmm.bg[3]=0.25;
  }
  /* Check given file names are valid */
  if (pflg){
    checkFile(peakfile, "r");
    /* If prediction file name is not provided, use initial model file name
       with suffix "_viterbi_results.txt" */
    if (!bflg) {
      strcat(predfile, outfile);
      strcat(predfile, "_viterbi_results.txt");
    }
    checkFile(predfile, "w");
  }
  if (fflg){
    checkFile(motiffile, "r");
    /* If prediction file name is not provided, use initial model file name
       with suffix "_with_probs.txt.txt" */
    if (!aflg) {
      strcat(scorefile, motiffile);
      strcat(scorefile, "_with_probs.txt");
    }
    checkFile(scorefile, "w");
  }
  
  /* Calculate PWM scores for each motif at each position */
  if (hmm.M > 0){
    pwm_matrix = gsl_matrix_alloc(hmm.M, T);
    CalMotifScore_P(&hmm, pwm_matrix, O1, P, peakPos);
  }
  gsl_vector * tmp_vector = gsl_vector_alloc(T);
  /* Set the initial mean parameters of PWM score feature based on max and min of actual calculation */
  index = 0;
  for (i = 0; i < hmm.M; i++) {
    gsl_matrix_get_row(tmp_vector, pwm_matrix, i);
    for (j = 0; j < hmm.N; j++) {
      gsl_matrix_set(hmm.mean_matrix, i, j, gsl_vector_min(tmp_vector) / 6.0);
      gsl_matrix_set(hmm.var_matrix, i, j, 8.0);
    }
    /* For each binding site state, the initial mean parameter for
     correspondong PWM score will be 1/2 of its greatest score */
    for (j = index; j < index + hmm.D[i] * (hmm.inactive + 1); j++) {
      gsl_matrix_set(hmm.mean_matrix, i, j, gsl_vector_max(tmp_vector) / 2.0);
    }
    index += hmm.D[i] * (hmm.inactive+1);
  }

  /* Set initial parameter values.
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
    
  /* Put PWM scores, counts and slopes into a matrix */
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
  }
  gsl_vector_free(tmp_vector);

  hmm.thresholds = (double *) dvector(hmm.M);
  /* with a given threshold file, TRACE can limit the state labeling to the
     regions with PWM score higher than the provided value */
  if (eflg) {
      fp = fopen(thresholdfile, "r");  //TODO: thresholds function
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
  
  /* matrix of alpha, beta and gamma in BW and viterbi algorithm */
  double **alpha = dmatrix(hmm.N, T);
  double **beta = dmatrix(hmm.N, T);
  double **gamma = dmatrix(hmm.N, T);
  double *logprobf = dvector(P); /* vector containing log likelihood for each peak */
  gsl_matrix * emission_matrix = gsl_matrix_alloc(hmm.N, T); /* matrix of emission probabilities */
  
  if (!ifSkip){
  /*                     */
  /* Start training step */
  /*                     */
    /* BW algorithm */
    BaumWelch(&hmm, T, obs_matrix, &niter, P, peakPos, logprobf, alpha, beta, gamma, emission_matrix);
    gsl_matrix_free(obs_matrix);
    fp = fopen(outfile, "w");
    /* Print the final model */
    PrintHMM(fp, &hmm);
    fclose(fp);
  }
  if (pflg || fflg){
  /*                     */
  /* Start decoding step */
  /*                     */
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
      fp1 = fopen(scorefile, "w");
      TF_end = getPosterior_motif(fp, fp1, T, peakPos, posterior, &hmm, q, vprob);
      fclose(fp);
      fclose(fp1);
    }
    if (pflg){
      fp = fopen(peakfile, "r");
      fp1 = fopen(predfile, "w");
      getPosterior_all(fp, fp1, T, q, peakPos, posterior, &hmm);
      fclose(fp);
      fclose(fp1);
    }
  }
}


void Usage(char *name)
{
  printf("Usage error. \n");
  printf("Usage1: %s [-v] ./TRACE <seq.file> <counts.file> <slope.file> "
         "--initial-model <init.model.file> --final-model <final.model.file> "
         "-peak-file <peak_3.file> --motif-file <peak_7.file>  --thread <thread.num>\n", name);
  printf("Usage2: %s [-v] ./TRACE --viterbi <seq.file> <counts.file> <slope.file> "
         "--final-model <final.model.file> -peak-file <peak_3.file> "
         "--motif-file <peak_7.file>  -T <thread.num>\n", name);
  printf("  seq.file - file containing the obs. sequence\n");
  printf("  count.file - file containing the obs. tag counts\n");
  printf("  slope.file - file containing the obs. slope\n");
  printf("  init.model.file - file with the initial model parameters\n");
  printf("  peak.file - file containing regions to detect TFBSs\n");
  printf("  final.model.file - output file containing the learned HMM\n");
}
