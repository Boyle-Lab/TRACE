/*
 *  File: fileutils.c
 *
 *  functions involved in readign and writing files
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>

#include "nrutil.h"
#include "hmm.h"
#include "logmath.h"
#include <omp.h>

/* Read the sequence file to get length T, GC content, sequence O,
 * number of peaks P and peak start position peakPos */
void ReadSequence(FILE *fp, int *pT, double *GC, int **pO, int *pP, int **peakPos)
{
  int *O, unused_num, *peaks;
  int i;
  unused_num = fscanf(fp, "T= %d\n", pT);
  unused_num = fscanf(fp, "GC: ");
  for (i = 0; i < 4; i++) {
    if(fscanf(fp, "%lf\t", &GC[i]) == EOF){
      fprintf(stderr, "Error: sequence file error \n");
      exit (1);
    }
  }
  unused_num = fscanf(fp,"\n");
  O = ivector(*pT);
  for (i = 0; i < *pT; i++) {
    if(fscanf(fp,"%d", &O[i]) == EOF){
      fprintf(stderr, "Error: sequence file error \n");
      exit (1);
    }
  }
  unused_num = fscanf(fp,"\n");
  *pO = O;
  unused_num = fscanf(fp, "P= %d\n", pP);
  peaks = ivector(*pP + 1);
  for (i=0; i < *pP + 1; i++){
    if(fscanf(fp,"%d", &peaks[i]) == EOF){
      fprintf(stderr, "Error: sequence file error \n");
      exit (1);
    }
  }
  *peakPos = peaks;
}

/* Read count or slope file to store numbers in data_vector,
 * with a optional adjust, which is used to change the scale of original numbers*/
void ReadTagFile(FILE *fp, int T, gsl_vector * data_vector, double adjust)
{
  double tmp;
  int i;
  for (i=0; i < T; i++) {
    if(fscanf(fp,"%lf\t", &tmp) == EOF){
      fprintf(stderr, "Error: input file error \n");
      exit (1);
    }
    gsl_vector_set(data_vector, i, tmp*adjust);
  }
}


void PrintSequenceProb(FILE *fp, int T, int *O, double *vprob, double *g,
                       double **posterior, int indexTF)
{
  int i;
  fprintf(fp,"%d\t%lf\t%lf\t%lf", O[0], vprob[0], g[0], posterior[0][indexTF]);
  for (i=1; i < T; i++) {
    fprintf(fp,"\n%d\t%lf\t%lf\t%lf", O[i], vprob[i], g[i], posterior[i][indexTF]);
  }
}

/* check if the provided file path is valid */
void checkFile(char *filename, char *mode)
{
  FILE *fp = fopen(filename, mode);
  if (fp == NULL) {
    fprintf(stderr, "Error: File %s not valid \n", filename);
    exit(1);
  }
  fclose(fp);
}
