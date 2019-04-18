#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "nrutil.h"
#include "hmm.h"
#include <sys/types.h>
#include <unistd.h>
#include <time.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

static char rcsid[] = "$Id: testvit.c,v 1.3 1998/02/23 07:39:07 kanungo Exp kanungo $";

void Usage(char *name);

int main (int argc, char **argv)
{
  FILE	*fp,*fp1,*fp2;
  clock_t time;
  int 	t, T, range, i, j; 
  HMM  	hmm;
  int	*O1;	/* observation sequence O[1..T] */
  double *O2, *O3;
  double *GC;
  double **delta, *vprob;
  double 	proba, *logproba; 
  double	*logprobf;
  double *g; 
  int **pos;
  
  int P;
  int *peakPos;
  int	*q;	/* state sequence q[1..T] */
  int	**psi;
  
  int	errflg=0, hflg=0, cflg=0, qflg=0, lflg=0, oflg=0, rflg=0, sflg=0, vflg=0;
  int	tflg=0, pflg=0, zflg=0, nflg=0, dflg=0, aflg=0, bflg=0, eflg=0, xflg=0;
  extern char *optarg;
  extern int optind, opterr, optopt;
  int	indexTF = 0; 
  char	*hmmfile, *slopfile, *countsfile, *seqfile, *outfile, *intsectfile;
  char  *listfile1, *listfile2, *thresholdfile, *pValuefile;
  int c;
  int ifMulti = 0, ifDbl=0, ifSet = 0, peakLength = 2;

  MAXITERATION = 200;
  THREAD_NUM = 40;
  gsl_vector * slop_vector, *counts_vector; 
  gsl_matrix *pwm_matrix, *obs_matrix;

  while ((c= getopt(argc, argv, "vhH:C:Q:L:A:B:T:P:N:D:S:E:V:O:X:Z")) != EOF){
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
    case 'V':
      /* get the p-value file name*/
      if (vflg)
        errflg++;
      else {
        vflg++;
        pValuefile = optarg;
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
    case 'O':
    /* get the HMM input file name */
      if (oflg)
        errflg++; 
      else { 
        oflg++;
        hmmfile = optarg;
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
    case 'C':  
    /* get the tag counts input file name */
      if (cflg) 
        errflg++; 
      else { 
        cflg++;  
        countsfile = optarg;
      } 
      break;  
    case 'Q':  
    /* get the sequence input file name*/
      if (qflg) 
        errflg++; 
      else { 
        qflg++;  
        seqfile = optarg;
      } 
      break;    
    case 'L':  
    /* get the slop input file name*/
      if (lflg) 
        errflg++; 
      else { 
        lflg++;  
        slopfile = optarg;
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
        intsectfile = optarg;
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
    case 'N':  
    /* choose independ normal or multivariance, default is independent(0)*/
      if (nflg) 
        errflg++; 
      else { 
        nflg++;  
        ifMulti = atoi(optarg);
      } 
      break;
    case 'S':
      /* if use the individually learned parameters */
      if (sflg)
        errflg++;
      else {
        sflg++;
        ifSet=atoi(optarg);
      }
            break;
    case '?':
      errflg++;
    }
	}
 
  if (!hflg || (!cflg && !qflg && !lflg)) {
    errflg++; 
  }
	if (errflg) { 
		Usage(argv[0]);
		exit (1);
	}
  
  
  if (qflg){ 
    //int *O1; /* sequence, represented by number, 0=A, 1=C, 2=G, 3=T */
    /* read the observed sequence */
    fp = fopen(seqfile, "r");
    if (fp == NULL) {
      fprintf(stderr, "Error: File %s not found \n", seqfile);
      exit (1);
    }
    GC = dvector(4);
    ReadSequence(fp, &T, GC, &O1, &P, &peakPos);
    fclose(fp);
  }
  
  
  if (lflg){
    //double  *O2; /* slops */
    /* read the slop file */
    fp = fopen(slopfile, "r");
    if (fp == NULL) {
      fprintf(stderr, "Error: File %s not found \n", slopfile);
      exit (1);
    }
    slop_vector = gsl_vector_alloc(T);
    ReadTagFile(fp, T, slop_vector, 1.0);
    fclose(fp);
  }
  
  if (cflg){
    //double  *O3; /* tag counts */
    /* read the tag counts file */
    fp = fopen(countsfile, "r");
    if (fp == NULL) {
      fprintf(stderr, "Error: File %s not found \n", countsfile);
      exit (1);
    }
    counts_vector = gsl_vector_alloc(T);
    ReadTagFile(fp, T, counts_vector, 2.0);
    fclose(fp);
  }
  
  
                                      
  /* read the hmm model */
  if (oflg) {
    /*read HMM input file */
    fp = fopen(hmmfile, "r");	
    if (fp == NULL) {
      fprintf(stderr, "Error: File %s not found \n", hmmfile);
      exit (1);
    }
    hmm.model = ifMulti;
    hmm.inactive = ifDbl;
    ReadM(fp, &hmm);
    hmm.K = cflg + lflg + qflg * hmm.M; /*number of data that can include 
                                          tag counts, slop, PWM score for each TF*/
    ReadHMM(fp, &hmm);
    hmm.lPeak = peakLength;
    fclose(fp);
    if (qflg){
      hmm.bg[0] = hmm.bg[3] = GC[0];
      hmm.bg[2] = hmm.bg[1] = GC[1];
      
    }
    else{
      hmm.bg[0]=hmm.bg[1]=hmm.bg[2]=hmm.bg[3]=0.25;
    }
  }
  
  /*
    HMM  	subhmm;
    fp = fopen("2_all_out_MAX_all.txt", "r");	
    if (fp == NULL) {
      fprintf(stderr, "Error: File %s not found \n", subhmm);
      exit (1);
    }
    ReadM(fp, &subhmm);
    subhmm.K = cflg + lflg + qflg * subhmm.M;
    ReadHMM(fp, &subhmm);
    fclose(fp);
    for (i = 0; i < 20; i++){
      gsl_matrix_set(hmm.mean_matrix, 16, 412+i, gsl_matrix_get(subhmm.mean_matrix, 0, i));   
      gsl_matrix_set(hmm.mean_matrix, 80, 412+i, gsl_matrix_get(subhmm.mean_matrix, 1, i));   
      gsl_matrix_set(hmm.mean_matrix, 81, 412+i, gsl_matrix_get(subhmm.mean_matrix, 2, i)); 
      gsl_matrix_set(hmm.var_matrix, 16, 412+i, gsl_matrix_get(subhmm.var_matrix, 0, i)); 
      gsl_matrix_set(hmm.var_matrix, 80, 412+i, gsl_matrix_get(subhmm.var_matrix, 1, i));   
      gsl_matrix_set(hmm.var_matrix, 81, 412+i, gsl_matrix_get(subhmm.var_matrix, 2, i)); 
      //fprintf(stdout, "%d: %lf %lf %lf  %lf %lf %lf  %lf %lf %lf\n", i, gsl_matrix_get(subhmm.mean_matrix, 0, i),gsl_matrix_get(subhmm.mean_matrix, 1, i),gsl_matrix_get(subhmm.mean_matrix, 2, i),gsl_matrix_get(subhmm.var_matrix, 0, i),gsl_matrix_get(subhmm.var_matrix, 1, i),gsl_matrix_get(subhmm.var_matrix, 2, i),gsl_matrix_get(hmm.var_matrix, 16, 412+i),gsl_matrix_get(hmm.var_matrix, 80, 412+i),gsl_matrix_get(hmm.var_matrix, 81, 412+i));   
    }
    for (i = 0; i < 4; i++){
      gsl_matrix_set(hmm.mean_matrix, 16, 2114+i, gsl_matrix_get(subhmm.mean_matrix, 0, 20+i));  
      gsl_matrix_set(hmm.mean_matrix, 80, 2114+i, gsl_matrix_get(subhmm.mean_matrix, 1, 20+i)); 
      gsl_matrix_set(hmm.mean_matrix, 81, 2114+i, gsl_matrix_get(subhmm.mean_matrix, 2, 20+i));
      gsl_matrix_set(hmm.var_matrix, 16, 2114+i, gsl_matrix_get(subhmm.var_matrix, 0, 20+i));  
      gsl_matrix_set(hmm.var_matrix, 80, 2114+i, gsl_matrix_get(subhmm.var_matrix, 1, 20+i)); 
      gsl_matrix_set(hmm.var_matrix, 81, 2114+i, gsl_matrix_get(subhmm.var_matrix, 2, 20+i));
      //fprintf(stdout, "%d: %lf %lf %lf\n", 20+i, gsl_matrix_get(subhmm.mean_matrix, 0, 20+i),gsl_matrix_get(subhmm.mean_matrix, 1, 20+i),gsl_matrix_get(subhmm.mean_matrix, 2, 20+i));   
    }
    for (i = 0; i < hmm.N; i++){
      covarMatrix_GSL(&hmm, i, hmm.cov_matrix[i]);
    }
    //fprintf(stdout, "%lf %lf %lf %lf\n", gsl_matrix_get(hmm.log_A_matrix, 2115, 412),gsl_matrix_get(hmm.log_A_matrix, 2117, 422),gsl_matrix_get(subhmm.log_A_matrix, 21,0),gsl_matrix_get(subhmm.log_A_matrix, 23,10)); 
    
    FreeHMM(&subhmm);
 */
  //int ifSet = 0;
  if(ifSet) {
    int index = 0;
    int i = 0;
    int j;
    char file_name[1000];
    char motif[20];
    char motif_[20];
    HMM  	subhmm;
    fp = fopen("../data/36_initHMM.txt", "r");
    while (fscanf(fp, "%s\t%s", motif,motif_) != EOF) {

      sprintf(file_name, "9_cluster_out_d_%s_ind_all_add__check3.txt", motif);
      printf("file: %s\n", file_name);
      fp1 = fopen(file_name, "r");
      if (fp1 != NULL) {
        ReadM(fp1, &subhmm);
        subhmm.K = cflg + lflg + qflg * subhmm.M;
        ReadHMM(fp1, &subhmm);
        fclose(fp1);
        for (j = 0; j < hmm.D[i]*2; j++){
          gsl_matrix_set(hmm.mean_matrix, i, index+j, gsl_matrix_get(subhmm.mean_matrix, 0, j));
          gsl_matrix_set(hmm.mean_matrix, hmm.M, index+j, gsl_matrix_get(subhmm.mean_matrix, 9, j));
          gsl_matrix_set(hmm.mean_matrix, hmm.M+1, index+j, gsl_matrix_get(subhmm.mean_matrix, 10, j));
        }
        //for (j = 0; j < 4; j++){
        //  gsl_matrix_set(hmm.mean_matrix, i, 2050+i*4+j, gsl_matrix_get(subhmm.mean_matrix, 0, hmm.D[i]*2+j));
        //  gsl_matrix_set(hmm.mean_matrix, 80, 2050+i*4+j, gsl_matrix_get(subhmm.mean_matrix, 1, hmm.D[i]*2+j));
        //  gsl_matrix_set(hmm.mean_matrix, 81, 2050+i*4+j, gsl_matrix_get(subhmm.mean_matrix, 2, hmm.D[i]*2+j));
        //}
      }
      index+=hmm.D[i]*2;
      i++;

    }
    fclose(fp);
  }
  
     
  fp = fopen("test.txt", "w");	
  PrintHMM(fp, &hmm);
  fclose(fp);
  fprintf(stdout, "M: %d T:%d K: %d\n", hmm.M,T,hmm.K);

  pwm_matrix = gsl_matrix_alloc(hmm.M, T);
  CalMotifScore_P(&hmm, pwm_matrix, O1, P, peakPos);
  
  free_ivector(O1, T);
  time = clock();
  printf("time check0.2: %f \n", ((double)time) / CLOCKS_PER_SEC);

  obs_matrix = gsl_matrix_alloc(hmm.K, T);
  //int i;
  gsl_vector * tmp_vector = gsl_vector_alloc(T);
  
  for (i = 0; i < hmm.M; i++){
    gsl_matrix_get_row(tmp_vector, pwm_matrix, i);
    gsl_matrix_set_row(obs_matrix, i, tmp_vector);
  }
  gsl_matrix_free(pwm_matrix);
  gsl_vector_free(tmp_vector);
  if (lflg){
    gsl_matrix_set_row(obs_matrix, hmm.M, slop_vector);
    gsl_vector_free(slop_vector);
  }
  if (cflg){
    gsl_matrix_set_row(obs_matrix, hmm.M+1, counts_vector);
    gsl_vector_free(counts_vector);
  }
  
  time = clock();
  printf("time check0.3: %f \n", ((double)time) / CLOCKS_PER_SEC);

  hmm.thresholds = (double *) dvector(hmm.M);
  if (eflg) {
    fp = fopen(thresholdfile, "r");  //TODO: thresholds function

    //fscanf(fp, "pi:\n");
    for (i = 0; i < hmm.M; i++) {
      fscanf(fp, "%lf\t", &(hmm.thresholds[i]));
    }
    fclose(fp);
  }
  else{
    for (i = 0; i < hmm.M; i++) {
      hmm.thresholds[i] = -INFINITY;
    }
  }
  gsl_matrix * emission_matrix = gsl_matrix_alloc(hmm.N, T);
  
  time = clock();
  printf("time check0.4: %f \n", ((double)time) / CLOCKS_PER_SEC);
  
  if (hmm.model == 0) EmissionMatrix(&hmm, obs_matrix, P, peakPos, emission_matrix, T);
  if (hmm.model == 1) EmissionMatrix_mv(&hmm, obs_matrix, P, peakPos, emission_matrix, T);
  if (hmm.model == 2) EmissionMatrix_mv_reduce(&hmm, obs_matrix, P, peakPos, emission_matrix, T);
  
  double **alpha = dmatrix(hmm.N, T);
  double **beta = dmatrix(hmm.N, T);
  double **gamma = dmatrix(hmm.N, T);
  time = clock();
  printf("time check e: %f \n", ((double)time) / CLOCKS_PER_SEC);
  logprobf = dvector(P); 
  Forward_P(&hmm, T, alpha, logprobf, P, peakPos, emission_matrix);
  time = clock();
  printf("time check a: %f \n", ((double)time) / CLOCKS_PER_SEC);
  Backward_P(&hmm, T, beta, P, peakPos, emission_matrix);
  time = clock();
  printf("time check b: %f \n", ((double)time) / CLOCKS_PER_SEC);
  ComputeGamma_P(&hmm, T, alpha, beta, gamma);
  time = clock();
  printf("time checkg: %f \n", ((double)time) / CLOCKS_PER_SEC);
  
  q = ivector(T);
  g = dvector(T);
  vprob = dvector(T);
  delta = dmatrix(T, hmm.N);
  psi = imatrix(T, hmm.N);
  logproba = dvector(P);
  time = clock();
  fprintf(stdout, "time check2: %f \n", ((double)time) / CLOCKS_PER_SEC);
  double  **posterior;
  posterior = dmatrix(T, hmm.N);
  fprintf(stdout, "indexTF = %d\n", indexTF);
  Viterbi(&hmm, T, g, alpha, beta, gamma, logprobf, delta, psi, q, 
          vprob, logproba, posterior, P, peakPos, emission_matrix);
  gsl_matrix_free(emission_matrix);
  free_dmatrix(delta, T, hmm.N);
  free_dmatrix(gamma, hmm.N, T);
  free_dmatrix(alpha, hmm.N, T);
  free_dmatrix(beta, hmm.N, T);

  time = clock();
  fprintf(stdout, "time check3: %f \n", ((double)time) / CLOCKS_PER_SEC);
  fprintf(stdout, "flg %d %d %d\n", pflg, aflg, bflg);
  //fflush(stdout);
  int TF_end;
  if (pflg){
    if (aflg & bflg){
      fp1 = fopen(listfile1, "w");
      if (fp1 == NULL) {
        fprintf(stderr, "Error: File %s not valid \n", listfile1);
        exit (1);
      }
      //PrintHMM(fp, &hmm);

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
      if (xflg){
        getPosterior_one_P(fp, fp1, fp2, indexTF, T, peakPos, posterior, &hmm, q, vprob);
        fclose(fp);
        fclose(fp1);
        fclose(fp2);
      }
      else {
        TF_end = getPosterior_all_P(fp, fp1, fp2, T, peakPos, posterior, &hmm, q, vprob);
        fclose(fp);
        fclose(fp1);
        fclose(fp2);
        strcat(hmmfile, "_viterbi_result_v.txt");
        fp1 = fopen(hmmfile, "w");
        fp = fopen("../data/v_tmp_peak_with.bed",
                   "r");//TODO:change the peak file
        if (fp == NULL) {
          fprintf(stderr, "Error: File %s not valid \n", intsectfile);
          exit(1);
        }
        getPosterior_labels(fp, fp1, T, q, peakPos, posterior, &hmm);
        fclose(fp);
        fclose(fp1);
      }
    }

    else{
      fp = fopen(intsectfile, "r");
      if (fp == NULL) {
        fprintf(stderr, "Error: File %s not valid \n", intsectfile);
        exit (1);
      }
      strcat(hmmfile,"_viterbi_result_v.txt");
      fp1 = fopen(hmmfile, "w");
      getPosterior_labels(fp, fp1, T, q, peakPos, posterior, &hmm);
      fclose(fp);
      fclose(fp1);
    }

  }

  /*
  fprintf(fp, "Viterbi  MLE log prob = %E\n", logproba);
  fprintf(fp, "Optimal state sequence:\n");
  fprintf(fp, "T= %d\n", T);
  PrintSequence(fp, T, q);
  fprintf(fp,"------------------------------------\n");
  fprintf(fp,"The two log probabilites and optimal state sequences\n");
  fprintf(fp,"should identical (within numerical precision). \n");
  */
  //char tmp_str[1000]; 
  //sprintf(tmp_str, "%d", T);
  //strcat(tmp_str, hmmfile);

  //PrintProb(fp, P, posterior2); 
  free_ivector(q, T);
  //free_ivector(O1, T);
  //free_dvector(O2, T);

  free_imatrix(psi, T, hmm.N);

  FreeHMM(&hmm);
}

void Usage(char *name)
{
  printf("Usage error. \n");
  printf("Usage: %s [-v] -H <mod.hmm> -Q <file.seq> -C <file.counts> -L <file.slop> -O <file.out> -R <list.out>\n", name);
  printf("  mod.hmm - file with the model parameters\n");
  printf("  file.counts - file containing the obs. tag counts\n");
  printf("  file.seq - file containing the obs. sequence\n");
  printf("  file.slop - file containing the obs. slop\n");
  printf("  file.out - output file containing the learned HMM and states at each base\n");
  printf("  list.out - output file containing list of states at each base, probbability\n");
}
