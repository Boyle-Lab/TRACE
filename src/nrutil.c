/*
**  File:   nrutil.c
**  codes were borrowed and modified from Kanungo's
**  original HMM program.
**  Tapas Kanungo, "UMDHMM: Hidden Markov Model Toolkit," in "Extended Finite State Models of Language," A. Kornai (editor), Cambridge University Press, 1999. http://www.kanungo.com/software/software.html.
**  Purpose: Memory allocation routines borrowed from the
**	book "Numerical Recipes" by Press, Flannery, Teukolsky,
**	and Vetterling.
**   state sequence and probablity of observing a sequence
**   given the model.
**
*/

#include <malloc.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include "hmm.h"
#include <math.h>

void nrerror(error_text)
char error_text[];
{
	void exit();

	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}



float *vector(nl,nh)
int nl,nh;
{
	float *v;

	v=(float *)calloc((unsigned) (nh-nl+1),sizeof(float));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl;
}

int *ivector(nh)
int nh;
{
	int *v;

	v=(int *)malloc((size_t)(sizeof(int)*nh));
	if (!v) nrerror("allocation failure in ivector()");
	return v;
}

/*
double *dvector(nl,nh)
int nl,nh;
{
	double *v;

	v=(double *)calloc((unsigned) (nh-nl+1),sizeof(double));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl;
}
*/

double *dvector(nh)
int nh;
{
	double *v;

	v=(double *)malloc((size_t)(sizeof(double)*nh));
	if (!v) nrerror("allocation failure in dvector()");
	return v;
}


float **matrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
{
	int i;
	float **m;

	m=(float **) calloc((unsigned) (nrh-nrl+1),sizeof(float*));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(float *) calloc((unsigned) (nch-ncl+1),sizeof(float));
		if (!m[i]) nrerror("allocation failure 2 in matrix()");
		m[i] -= ncl;
	}
	return m;
}

/*
double **dmatrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
{
	int i;
	double **m;

	m=(double **) calloc((unsigned) (nrh-nrl+1),sizeof(double*));
	if (!m) nrerror("allocation failure 1 in dmatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(double *) calloc((unsigned) (nch-ncl+1),sizeof(double));
		if (!m[i]) nrerror("allocation failure 2 in dmatrix()");
		m[i] -= ncl;
	}
	return m;
}
*/

double **dmatrix(nrh,nch)
int nrh,nch;
{
	int i;
	double **m;

	//m=(double **) calloc((unsigned) nrh,sizeof(double*));
  m=(double **) malloc((unsigned)(nrh*sizeof(double*)));
	if (!m) nrerror("allocation failure 1 in dmatrix()");

	for(i=0;i<nrh;i++) {
		m[i]=(double *) malloc((unsigned)(nch*sizeof(double)));
		if (!m[i]) {
      fprintf(stdout, "at %d:", i);
      nrerror("allocation failure 2 in dmatrix()");
    } 
	}
	return m;
}

int **imatrix(nrh,nch)
int nrh,nch;
{
	int i,**m;

	m=(int **)malloc((size_t)(sizeof(int*)*nrh));
	if (!m) nrerror("allocation failure 1 in imatrix()");

	for(i=0;i<nrh;i++) {
		m[i]=(int *)malloc((size_t)(sizeof(int)*nch));
		if (!m[i]) nrerror("allocation failure 2 in imatrix()");
	}
	return m;
}



float **submatrix(a,oldrl,oldrh,oldcl,oldch,newrl,newcl)
float **a;
int oldrl,oldrh,oldcl,oldch,newrl,newcl;
{
	int i,j;
	float **m;

	m=(float **) calloc((unsigned) (oldrh-oldrl+1),sizeof(float*));
	if (!m) nrerror("allocation failure in submatrix()");
	m -= newrl;

	for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+oldcl-newcl;

	return m;
}



void free_vector(v,nl,nh)
float *v;
int nl,nh;
{
	free((float*) (v+nl));
}

void free_ivector(v,nh)
int *v,nh;
{
	free((int*) (v));
  v = NULL;
}

void free_dvector(v,nh)
double *v;
int nh;
{
	free((double*) (v));
  v = NULL;
}



void free_matrix(m,nrl,nrh,ncl,nch)
float **m;
int nrl,nrh,ncl,nch;
{
	int i;

	for(i=nrh;i>=nrl;i--) free((float*) (m[i]+ncl));
	free((float*) (m+nrl));
}
/*
void free_dmatrix(m,nrl,nrh,ncl,nch)
double **m;
int nrl,nrh,ncl,nch;
{
	int i;

	for(i=nrh;i>=nrl;i--) free((m[i]+ncl));
	free((m+nrl));
}

void free_dmatrix(m,nrl,nrh,ncl,nch)
double **m;
int nrl,nrh,ncl,nch;
{
	int i;

	for(i=nrl;i<=nrh;i++) free(m[i]);
	//free((m));
}
*/
void free_dmatrix(m,nrh,nch)
double **m;
int nrh,nch;
{
	int i;

	//for(i=nrh;i>=nrl;i--) free((m[i]));
  for(i=0;i<nrh;i++) {
    free(m[i]);
  }
  free(m);
}

void free_imatrix(m,nrl,nrh,ncl,nch)
int **m;
int nrl,nrh,ncl,nch;
{
	int i;

	for(i=0;i<nrh;i++) {
    free((m[i]));
	}
  free((m));
  m = NULL;
}



void free_submatrix(b,nrl,nrh,ncl,nch)
float **b;
int nrl,nrh,ncl,nch;
{
	free((float*) (b+nrl));
}



float **convert_matrix(a,nrl,nrh,ncl,nch)
float *a;
int nrl,nrh,ncl,nch;
{
	int i,j,nrow,ncol;
	float **m;

	nrow=nrh-nrl+1;
	ncol=nch-ncl+1;
	m = (float **) calloc((unsigned) (nrow),sizeof(float*));
	if (!m) nrerror("allocation failure in convert_matrix()");
	m;
	for(i=0,j=0;i<=nrow-1;i++,j++) m[j]=a+ncol*i;
	return m;
}



void free_convert_matrix(b,nrl,nrh,ncl,nch)
float **b;
int nrl,nrh,ncl,nch;
{
	free((char*) (b+nrl));
}

void copyRow(rin,l,rout)
double *rin, *rout;
int l;
{
	int i;

	for(i=0;i<l;i++) {
    rout[i]=rin[i];
  }
}

void copyMatrix(min,nrl,ncl,mout)
double **min, **mout;
int nrl,ncl;
{
	int i,j;

	for(i=0;i<nrl;i++) {
    for(j=0;j<ncl;j++) {
      mout[i][j]=min[i][j];
    }
  }
}

/*print a matrix of (size1 x size2)*/
void printMatrix(fp, matrix, size1, size2)
FILE *fp;
double **matrix;
int size1, size2;

{
	int i,j;
	for (i = 0; i < size1; i++){
		for (j = 0; j < size2; j++){
			fprintf(fp, "%lf\t", matrix[i][j]);
		}
		fprintf(fp, "\n");
	}
}

void printVector(fp, matrix, size)
FILE *fp;
double *matrix;
int size;
{
	int i,j;
	for (i = 0; i < size; i++){
		fprintf(fp, "%lf\t", matrix[i]);
	}
	fprintf(fp, "\n");
}

/*print a matrix of (size1 x size2)*/
void printfMatrix(matrix, size1, size2)
float **matrix;
int size1, size2;
{
	int i,j;
	for (i = 0; i < size1; i++){
		for (j = 0; j < size2; j++){
			fprintf(stdout, "%lf ", matrix[i][j]);
		}
		fprintf(stdout, "\n");
	}
}

double listInsertnMax(list, insert, size)
double *list, insert;
int size;
{
	int i;
	double max = log(0.0);
	for (i = (size - 1); i > 0; i--){
		list[i] = list[i-1];
		max = MAX(max, list[i]);
	}
	list[0] = insert;
	max = MAX(max, list[0]);
	return max;
}
