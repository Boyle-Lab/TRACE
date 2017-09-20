#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "logmath.h"

typedef struct {
/* make sure to put TF motif at the first states */
  int N;		/* number of states;  Q={1,2,...,N} */
  int M; 		/* number of TFs*/
  int K; /*number of values in each state. 1 + N here (slop and pwm scrore for each TF)*/
         /*XXX: might need to change to take input later*/ 
  double	**A;	/* A[1..N][1..N]. a[i][j] is the transition prob
			   of going from state i at time t to state j
			   at time t+1 */
  //int P; /*number of the peaks*/
  //int peakPos; /*position for each peak */
  double **mu; /* mu[i] is mean of slop in state i */
  double **sigma; /* sigma[i] is variance in state i */
  double **rho; /*correlation*/
  double ***pwm; /* pwm[i][n][x] is probability of each base in each state from PWM for TF i*/
  int *D; /* D[i] is number of positions in motif */
  double *bg; /*background percentage for ACGT */
  double	*pi;	/* pi[1..N] pi[i] is the initial state distribution. */
} HMM;

//void ReadHMM(FILE *fp, HMM *phmm);
void ReadOutHMM(FILE *fp, HMM *phmm);
//void InitHMMwithInput(HMM *phmm, int seed, double *gc);
void PrintHMM(FILE *fp, HMM *phmm);
void InitHMM(HMM *phmm, int N, int M, int D, int seed);
//void CopyHMM(HMM *phmm1, HMM *phmm2);
void FreeHMM(HMM *phmm);
void printMatrix(double **matrix, int size1, int size2);
void ReadSequence(FILE *fp, int *pT, double **GC, int **pO, int *pP, int **peakPos);
void CalMotifScore(HMM *phmm, double **S, int *O1, int P, int *peakPos); 
void ReadSlop(FILE *fp, int *pT, double **pO);
void PrintSequence(FILE *fp, int T, int *O);
void PrintSequenceProb(FILE *fp, int T, int *O, double *vprob, double *g);
//void GenSequenceArray(HMM *phmm, int seed, int T, int *O, int *q);
//int GenInitalState(HMM *phmm);
//int GenNextState(HMM *phmm, int q_t);
//int GenSymbol(HMM *phmm, int q_t);
double ComputeEmission(HMM *phmm, int j, double *data); 
void covarMatrix(HMM *phmm, int state, double **matrix);

void ForwardMultiSeq(HMM *phmm, int T, double **S, double *O2, double **alpha, double *pprob, int P, int *peakPos);

void BackwardMultiSeq(HMM *phmm, int T, double **S, double *O2, double **beta, double *pprob, int P, int *peakPos);

void BaumWelchBVN(HMM *phmm, int T, int *O1, double *O2, double **alpha, 
  double **beta, double **gamma, int *pniter, int P, int *peakPos);

double *** AllocXi(int T, int N);
void FreeXi(double *** xi, int T, int N);

void ComputeGammaWithLog(HMM *phmm, int T, double **alpha, double **beta,
        double **gamma);
void ComputeXiWithLog(HMM* phmm, int T, double **S, double *O1, double **alpha, double **beta,
        double ***xi);

void ViterbiMultiSeq(HMM *phmm, int T, int *O1, double *O2, double *S, double *g, double **delta, int **psi,
        int *q, double *vprob, double *pprob, int P, int *peakPos);

/* random number generator related functions*/

int hmmgetseed(void);
void hmmsetseed(int seed);
double hmmgetrand(void);
 
#define MAX(x,y)        ((x) > (y) ? (x) : (y))
#define MIN(x,y)        ((x) < (y) ? (x) : (y))
 
