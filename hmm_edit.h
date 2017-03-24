

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "logmath.h"

typedef struct {
/* make sure to put TF motif at the first states */
	int N;		/* number of states;  Q={1,2,...,N} */
	int M; 		/* number of observation symbols; V={1,2,...,M}*/
	double	**A;	/* A[1..N][1..N]. a[i][j] is the transition prob
			   of going from state i at time t to state j
			   at time t+1 */
	//int P; /*number of the peaks*/
  //int peakPos; /*position for each peak */
  double *mu; /* mu[i] mean of the normal distribution */
  double *sigma; /* sigma[i] standard diviation of the normal distribution */
  double **pwm; /* pwm[n][x] probability of each base in each state from PWM */
  int D; /* number of positions in motif */
  double *CG; /*CG content */
	double	*pi;	/* pi[1..N] pi[i] is the initial state distribution. */
} HMM;
void ReadHMM(FILE *fp, HMM *phmm);
void ReadOutHMM(FILE *fp, HMM *phmm);
void InitHMMwithInput(HMM *phmm, int seed, double *gc);
void PrintHMM(FILE *fp, HMM *phmm);
void InitHMM(HMM *phmm, int N, int M, int D, int seed);
void CopyHMM(HMM *phmm1, HMM *phmm2);
void FreeHMM(HMM *phmm);

void ReadSequence(FILE *fp, int *pT, double **GC, int **pO, int *pP, int **peakPos);
void ReadSlop(FILE *fp, int *pT, double **pO);
void PrintSequence(FILE *fp, int T, int *O);
void PrintSequenceProb(FILE *fp, int T, int *O, double *vprob);
void GenSequenceArray(HMM *phmm, int seed, int T, int *O, int *q);
int GenInitalState(HMM *phmm);
int GenNextState(HMM *phmm, int q_t);
int GenSymbol(HMM *phmm, int q_t);

 
void ForwardMultiSeq(HMM *phmm, int T, int *O1, double *O2, double **alpha, double *pprob, int P, int *peakPos);
void BackwardMultiSeq(HMM *phmm, int T, int *O1, double *O2, double **beta, double *pprob, int P, int *peakPos);
double ComputeEmission(HMM *phmm, int j, int t, int *O1, double *O2); 
void BaumWelchMultiSeq(HMM *phmm, int T, int *O1, double *O2, double **alpha, 
  double **beta, double **gamma, int *pniter, int P, int *peakPos);

double *** AllocXi(int T, int N);
void FreeXi(double *** xi, int T, int N);
void ComputeGammaWithLog(HMM *phmm, int T, double **alpha, double **beta,
        double **gamma);
void ComputeXiWithLog(HMM* phmm, int T, int *O1, double *O2, double **alpha, double **beta,
        double ***xi);
void ViterbiMultiSeq(HMM *phmm, int T, int *O1, double *O2, double **delta, int **psi,
        int *q, double *vprob, double *pprob, int P, int *peakPos);

/* random number generator related functions*/

int hmmgetseed(void);
void hmmsetseed(int seed);
double hmmgetrand(void);
 
#define MAX(x,y)        ((x) > (y) ? (x) : (y))
#define MIN(x,y)        ((x) < (y) ? (x) : (y))
 

