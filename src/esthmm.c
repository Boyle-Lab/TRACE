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

void Usage(char *name);

int main (int argc, char **argv)
{
  /* The following are read from input files */
  FILE *fp, *fp1, *fp2;
  HMM hmm;
  int P; /* total number of peaks */
  int *peakPos; /* Starting location of each peaks*/
                /*there two are used to seperate values from each peak in concatenated data*/
  int T; /* Total length*/
  int *O1; /* sequence, represented by number, 1=A, 2=C, 3=G, 4=T */
  double *GC; /* GC content */
  gsl_matrix *pwm_matrix, *obs_matrix; /* matrix of PWM scores, and observation data
             /* right now, there are only three observations: counts, slop and sequence
              might need to change the structure if more data are used in the future*/
  gsl_vector * slop_vector, *counts_vector;

  int niter; /*numbers of iterations, might need to set a max later*/
  int seed; /* seed for random number generator */
  int i, k;

  /* Get all command-line options */
  int c, w = 0, indexTF, readFile = 0;
  int errflg =0, iflg=0, sflg=0, cflg=0, qflg=0, lflg=0, oflg=0, vflg=0;
  int mflg=0, kflg=0, wflg=0, zflg=0, nflg=0, aflg =0, bflg =0, pflg =0, eflg=0;
  int dflg=0, tflg = 0, rflg = 0, hflg = 0, fflg=0;
  extern int optind, opterr, optopt;
  extern char *optarg;
  char *hmminitfile, *slopinitfile, *countsinitfile, *seqinitfile;
  char *listfile1, *listfile2, intsectfile[50], outfile[50], *TPfile;
  char *thresholdfile, *pValuefile;
  int ifMulti = 2, ifDbl=1, ifSet = 1, peakLength = 2;
  intsectfile[0] = '\0';
  outfile[0] = '\0';
  MAXITERATION = 200;
  THREAD_NUM = 40;
  while ((c= getopt(argc, argv, "vhI:S:C:Q:L:O:M:K:W:N:F:A:B:P:T:D:R:H:E:V:Z")) != EOF){
    switch (c) {
      case 'v':
        vflg++;
        break;
      case 'h':
        Usage(argv[0]);
        exit(1);
        break;
      case 'I':
      /* get the HMM input file name */
        if (iflg)
          errflg++;
        else {
          iflg++;
          hmminitfile = optarg;
        }
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
      case 'Z':
      /* set z score or not */
        if (zflg)
          errflg++;
        else {
          zflg++;
        }
        break;
      case 'S':
      /* if use the individually learned parameters */
        if (sflg)
          errflg++;
        else {
          sflg++;
          ifSet = atoi(optarg);
        }
        break;
      case 'W':
      /* set whether to use P(O) in BW */
        if (wflg)
          errflg++;
        else {
          wflg++;
          w = atoi(optarg);
        }
        break;
      case 'K':
      /* set random number generator seed */
        if (kflg)
          errflg++;
        else {
          kflg++;
          k = atoi(optarg);
        }
        break;
      case 'R':
        /* get the TP/FP input file name */
        if (rflg)
          errflg++;
        else {
          rflg++;
          TPfile = optarg;
        }
        break;
      case 'C':
      /* get the tag counts input file name */
        if (cflg)
          errflg++;
        else {
          cflg++;
          countsinitfile = optarg;
        }
        break;
      case 'Q':
      /* get the sequence input file name*/
        if (qflg)
          errflg++;
        else {
          qflg++;
          seqinitfile = optarg;
        }
        break;
      case 'L':
      /* get the slop input file name*/
        if (lflg)
          errflg++;
        else {
          lflg++;
          slopinitfile = optarg;
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
      /* get the pwm score input file name*/
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
      case 'D':
      /* if double the TFs */
        if (dflg)
          errflg++;
        else {
          dflg++;
          ifDbl = atoi(optarg);
        }
        break;
      case 'F':
      /* choose whether to read actrual number from file, default id not (0)*/
        if (fflg)
          errflg++;
        else {
          fflg++;
          readFile = atoi(optarg);
        }
        break;
      case 'H':
      /* length of states in peak*/
        if (hflg)
          errflg++;
        else {
          hflg++;
          peakLength = atoi(optarg);
        }
        break;
      case 'A':
      /* get the output file nam */
        if (aflg)
          errflg++;
        else {
          aflg++;
          listfile1 = optarg;
        }
        break;
      case 'B':
      /* get the output file name*/
        if (bflg)
          errflg++;
        else {
          bflg++;
          listfile2 = optarg;
        }
        break;
      case 'P':
      /* get the intersect input file name*/
        if (pflg)
          errflg++;
        else {
          pflg++;
          strcat(intsectfile, optarg);
        }
        break;
      case 'V':
      /* get the p-value file name*/
        if (vflg)
          errflg++;
        else {
          vflg++;
          pValuefile = optarg;
        }
        break;
      case '?':
        errflg++;
    }
  }

  /*al least one data file should be given */
  if ((!cflg && !qflg && !lflg) && (!mflg || !kflg) || !oflg) {
    errflg++; 
  }
  if (errflg) {
    Usage(argv[0]);
    exit (1);
  }
  if (qflg){
    /* read the observed sequence */
    fp = fopen(seqinitfile, "r");
    if (fp == NULL) {
      fprintf(stderr, "Error: File %s not found \n", seqinitfile);
      exit (1);
    }
    GC = dvector(4);
    ReadSequence(fp, &T, GC, &O1, &P, &peakPos); 
    fclose(fp);
  }
  if (lflg){
    /* read the slop file */
    fp = fopen(slopinitfile, "r");
    if (fp == NULL) {
      fprintf(stderr, "Error: File %s not found \n", slopinitfile);
      exit (1);
    }
    slop_vector = gsl_vector_alloc(T);
    ReadTagFile(fp, T, slop_vector, 1.0);
    fclose(fp);
  }
  if (cflg){
    /* read the tag counts file */
    fp = fopen(countsinitfile, "r");
    if (fp == NULL) {
      fprintf(stderr, "Error: File %s not found \n", countsinitfile);
      exit (1);
    }
    counts_vector = gsl_vector_alloc(T);
    ReadTagFile(fp, T, counts_vector, 2.0);
    fclose(fp);
  }
  /* initialize the hmm model */
  if (iflg) { 
    /*read HMM input file */
    fp = fopen(hmminitfile, "r");	
    if (fp == NULL) {
      fprintf(stderr, "Error: File %s not found \n", hmminitfile);
      exit (1);
    }
    hmm.model = ifMulti; 
    hmm.inactive = ifDbl; 
    ReadM(fp, &hmm);
    hmm.K = cflg + lflg + qflg * hmm.M; /*number of data that can include 
                                          tag counts, slop, PWM score for each TF*/
    if (readFile == 0) ReadInitHMM(fp, &hmm);
    if (readFile == 1) ReadHMM(fp, &hmm);
    fclose(fp);
    hmm.lPeak = peakLength;
    if (qflg){
      hmm.bg[0] = hmm.bg[3] = GC[0];
      hmm.bg[2] = hmm.bg[1] = GC[1];
    }
    else{
      hmm.bg[0]=hmm.bg[1]=hmm.bg[2]=hmm.bg[3]=0.25;
    }
  }

  if (aflg & bflg){
    fp1 = fopen(listfile1, "w");
    if (fp1 == NULL) {
      fprintf(stderr, "Error: File %s not valid \n", listfile1);
      exit (1);
    }
    fp2 = fopen(listfile2, "w");
    if (fp2 == NULL) {
      fprintf(stderr, "Error: File %s not valid \n", listfile2);
      exit (1);
    }
    fp = fopen(intsectfile, "r");
    if (fp == NULL) {
      fprintf(stderr, "Error: File %s not valid \n", intsectfile);
      exit(1);
    }
    fclose(fp);
    fclose(fp1);
    fclose(fp2);
  }
  fprintf(stdout, "extra: %d M: %d N: %d T:%d K: %d\n", hmm.extraState,hmm.M,hmm.N,T,hmm.K);

  /* Calculate PWM scores for each motif at each position */
  if (hmm.M > 0){
    pwm_matrix = gsl_matrix_alloc(hmm.M, T);
    CalMotifScore_P(&hmm, pwm_matrix, O1, P, peakPos);
  }
/*
  int i, j, index;
  gsl_vector * tmp_vector = gsl_vector_alloc(T);
  index = 0;
  for (i = 0; i < hmm.M; i++){
    gsl_matrix_get_row(tmp_vector, pwm_matrix, i);
    printf("max min %d: %f %f\n", i, gsl_vector_max(tmp_vector), gsl_vector_min(tmp_vector));
    for (j = 0; j < hmm.N; j++) {
      gsl_matrix_set(hmm.mean_matrix, i, j, gsl_vector_min(tmp_vector)/2.0);
      gsl_matrix_set(hmm.var_matrix, i, j, 8.0);
    }
    for (j = index; j < index+hmm.D[i] * (hmm.inactive+1); j++) {
      gsl_matrix_set(hmm.mean_matrix, i, j, gsl_vector_max(tmp_vector)/2.0);
    }
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
  gsl_vector * tmp_vector = gsl_vector_alloc(T);
  obs_matrix = gsl_matrix_alloc(hmm.K, T);
  for (i = 0; i < hmm.M; i++){
    gsl_matrix_get_row(tmp_vector, pwm_matrix, i);
    gsl_matrix_set_row(obs_matrix, i, tmp_vector);
  }
  if (lflg){
    gsl_matrix_set_row(obs_matrix, hmm.M, slop_vector);
    gsl_vector_free(slop_vector);
  }
  if (cflg){
    gsl_matrix_set_row(obs_matrix, hmm.M+1, counts_vector);
    gsl_vector_free(counts_vector);
  }
  if (hmm.M > 0){
    gsl_matrix_free(pwm_matrix);
    gsl_vector_free(tmp_vector);
  }

  /* Check required files are valid */
  fp = fopen(outfile, "w");
  if (fp == NULL) {
    fprintf(stderr, "Error: File %s not valid \n", outfile);
    exit (1);
  }
  fclose(fp);
  if (rflg){
    fp = fopen(TPfile, "r");
    getParameters_all_P(fp, &hmm, T, obs_matrix, P, peakPos);
    fclose(fp);
  }
  if (pflg){
    fp1 = fopen(intsectfile, "r");
    if (fp1 == NULL) {
      fprintf(stderr, "Error: File %s not valid \n", intsectfile);
      exit (1);
    }
    fclose(fp1);
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
      fscanf(fp, "%lf\n", &(hmm.thresholds[i]));
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
  double *logprobf = dvector(P); /*vector containing log likelihood 
                                        for each peak*/
  gsl_matrix * emission_matrix = gsl_matrix_alloc(hmm.N, T);
  BaumWelch(&hmm, T, obs_matrix, &niter, P, peakPos, logprobf, alpha, beta, gamma, emission_matrix);
  gsl_matrix_free(obs_matrix);
  fp = fopen(outfile, "w");
  /* Print the final model */
  PrintHMM(fp, &hmm);
  fclose(fp);
  
  /* Viterbi step */
  int	*q = ivector(T);	/* state sequence q[1..T] */
  int	**psi = imatrix(T, hmm.N);
  double *g = dvector(T);
  double *vprob = dvector(T);
  double **delta = dmatrix(T, hmm.N);
  double *logproba = dvector(P);
  double  **posterior = dmatrix(T, hmm.N);
  Viterbi(&hmm, T, g, alpha, beta, gamma, logprobf, delta, psi, q,
          vprob, logproba, posterior, P, peakPos, emission_matrix);
  int TF_end;
  if (pflg){
    if (aflg & bflg){
      fp1 = fopen(listfile1, "w");
      if (fp1 == NULL) {
        fprintf(stderr, "Error: File %s not valid \n", listfile1);
        exit (1);
      }
      fp2 = fopen(listfile2, "w");
      if (fp2 == NULL) {
        fprintf(stderr, "Error: File %s not valid \n", listfile2);
        exit (1);
      }
      fp = fopen(intsectfile, "r");
      if (fp == NULL) {
        fprintf(stderr, "Error: File %s not valid \n", intsectfile);
        exit (1);
      }
      TF_end = getPosterior_all_P(fp, fp1, fp2, T, peakPos, posterior, &hmm, q, vprob);
      fclose(fp);
      fclose(fp1);
      fclose(fp2);
      strcat(outfile,"_viterbi_result.txt");
      fp1 = fopen(outfile, "w");
      strcat(intsectfile,"_with.txt");
      fp = fopen(intsectfile, "r");
      if (fp == NULL) {
        fprintf(stderr, "Error: File %s not valid \n", intsectfile);
        exit (1);
      }
      getPosterior_labels(fp, fp1, T, q, peakPos, posterior, &hmm);
      fclose(fp);
      fclose(fp1);
    }
    else{
      fp = fopen(intsectfile, "r");
      if (fp == NULL) {
        fprintf(stderr, "Error: File %s not valid \n", intsectfile);
        exit (1);
      }
      strcat(outfile,"_viterbi_result.txt");
      fp1 = fopen(outfile, "w");
      getPosterior_labels(fp, fp1, T, q, peakPos, posterior, &hmm);
      fclose(fp);
      fclose(fp1);
    }

  }
}


void Usage(char *name)
{
  printf("Usage error. \n");
  printf("Usage1: %s [-v] ./esthmm -Q <sequence.file> -L <slope.file> -C <counts.file> "
         "-I <init.model.file> -O <final.model.file> -P <peak.file> -A <output1.file> "
         "-B <output2.file> -T <thread.num>\n", name);
  printf("Usage1: %s [-v] ./esthmm -Q <sequence.file> -L <slope.file> -C <counts.file> "
         "-I <init.model.file> -O <final.model.file> -P <peak.file> -T <thread.num>\n", name);
  printf("  init.model.file - file with the initial model parameters\n");
  printf("  file.counts - file containing the obs. tag counts\n");
  printf("  sequence.file - file containing the obs. sequence\n");
  printf("  slope.file - file containing the obs. slope\n");
  printf("  peak.file - file containing regions to detect TFBSs\n");
  printf("  final.model.file - output file containing the learned HMM\n");
}
