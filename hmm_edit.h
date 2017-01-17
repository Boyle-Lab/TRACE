

#include <stdlib.h>
#include <stdio.h>
#include <math.h>


typedef struct {
/* make sure to put TF motif at the first states */
	int N;		/* number of states;  Q={1,2,...,N} */
	int M; 		/* number of observation symbols; V={1,2,...,M}*/
	double	**A;	/* A[1..N][1..N]. a[i][j] is the transition prob
			   of going from state i at time t to state j
			   at time t+1 */
	
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

void ReadSequence(FILE *fp, int *pT, double **GC, int **pO);
void ReadSlop(FILE *fp, int *pT, double **pO);
void PrintSequence(FILE *fp, int T, int *O);
void GenSequenceArray(HMM *phmm, int seed, int T, int *O, int *q);
int GenInitalState(HMM *phmm);
int GenNextState(HMM *phmm, int q_t);
int GenSymbol(HMM *phmm, int q_t);

 

void ForwardWithScale(HMM *phmm, int T, int *O1, double *O2, double **alpha,
        double *scale, double *pprob);
double ComputeEmission(HMM *phmm, int j, int t, int *O1, double *O2); 
void BackwardWithScale(HMM *phmm, int T, int *O1, double *O2, double **beta,
        double *scale, double *pprob);
void BaumWelch(HMM *phmm, int T, int *O1, double *O2, double **alpha, 
  double **beta, double **gamma, int *pniter, 
	double *plogprobinit, double *plogprobfinal);


double *** AllocXi(int T, int N);
void FreeXi(double *** xi, int T, int N);
void ComputeGamma(HMM *phmm, int T, double **alpha, double **beta,
        double **gamma);
void ComputeXi(HMM* phmm, int T, int *O1, double *O2, double **alpha, double **beta,
        double ***xi);
void Viterbi(HMM *phmm, int T, int *O1, double *O2, double **delta, int **psi,
        int *q, double *pprob);
void ViterbiLog(HMM *phmm, int T, int *O1, double *O2, double **delta, int **psi,
        int *q, double *pprob);

/* random number generator related functions*/

int hmmgetseed(void);
void hmmsetseed(int seed);
double hmmgetrand(void);
 
#define MAX(x,y)        ((x) > (y) ? (x) : (y))
#define MIN(x,y)        ((x) < (y) ? (x) : (y))
 

