#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "logmath.h"
#include <gsl/gsl_matrix.h>

typedef struct {
/* make sure to put TF motif at the first states */
  int N; /* number of states;  Q={1,2,...,N} */
  int M; /* number of TFs*/
  int K; /*number of values in each state. (slop and pwm scrore for each TF)*/
  int inactive;
  int model; /* 0: independent, 1: multivariance */
  int lPeak; /* length of states in peak surrounding FPs */
  int extraState; /* number of states that are not FPs */
  double **A; /* A[1..N][1..N]. a[i][j] is the transition prob
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
  double *thresholds; /* thresholds of pwm scores */
  gsl_matrix ** cov_matrix; /* N x N Covariance matrix */
  gsl_matrix * mean_matrix; /* K x N mean matrix, each row is a vector of means
                             * for each state for each feature*/
  gsl_matrix * var_matrix; /* K x N standard deviation matrix, each row is a vector
                            * of std for each state for each feature*/
  gsl_matrix * log_A_matrix; /* N x N transition prob matrix */
} HMM;

/* BaumWelch.c */
void BaumWelch(HMM *phmm, int T, gsl_matrix * obs_matrix, int *pniter, int P, 
               int *peakPos, double *logprobf, double **alpha, double **beta, 
               double **gamma, gsl_matrix * emission_matrix);
void UpdateVariance(HMM *phmm, gsl_matrix * obs_matrix, gsl_matrix * post_obs,
                    gsl_vector * prob_sum, gsl_matrix *prob_matrix, int T, 
                    int TF);
void UpdateCovariance(HMM *phmm, gsl_matrix * obs_matrix, 
                      gsl_matrix * post_obs, gsl_vector * prob_sum, 
                      gsl_matrix *prob_matrix, int T, int TF);
void UpdateVariance_2(HMM *phmm, gsl_matrix * obs_matrix, 
                      gsl_vector * prob_sum, gsl_matrix *prob_matrix, 
                      int T, int TF);
void UpdateCovariance_2(HMM *phmm, gsl_matrix * obs_matrix, 
                        gsl_vector * prob_sum, gsl_matrix *prob_matrix, 
                        int T, int TF);
void ComputeGamma(HMM *phmm, int T, gsl_matrix * alpha_matrix, 
                  gsl_matrix * beta_matrix, gsl_matrix * gamma_matrix);
void ComputeXi_sum(HMM* phmm, gsl_matrix * alpha_matrix, 
                   gsl_matrix * beta_matrix, gsl_vector * xi_sum_vector,
                   gsl_matrix * emission_matrix, int T);
void ComputeXi_sum_P(HMM* phmm, double **alpha, double **beta, double *xi_sum, 
                   gsl_matrix * emission_matrix, int T);
void ComputeGamma_P(HMM *phmm, int T, double **alpha, double **beta, 
                    double **gamma);

/*fwd_bwd.c*/
void Forward_P(HMM *phmm, int T, double **alpha, double *pprob, int P, 
               int *peakPos, gsl_matrix * emission_matrix);
void Backward_P(HMM *phmm, int T, double **beta, int P, int *peakPos, 
                gsl_matrix * emission_matrix);
                
/* sequence.c */
void ReadSequence(FILE *fp, int *pT, double *GC, int **pO, int *pP, 
                  int **peakPos);
void ReadTagFile(FILE *fp, int T, gsl_vector * data_vector, double adjust);
void CalMotifScore_P(HMM *phmm, gsl_matrix * S, int *O1, int P, int *peakPos);
void EmissionMatrix(HMM* phmm, gsl_matrix * obs_matrix, int P, int *peakPos, 
                    gsl_matrix * emission_matrix, int T);
void EmissionMatrix_mv(HMM* phmm, gsl_matrix * obs_matrix, int P, int *peakPos,
                       gsl_matrix * emission_matrix, int T);
void EmissionMatrix_mv_reduce(HMM* phmm, gsl_matrix * obs_matrix, int P, 
                              int *peakPos, gsl_matrix * emission_matrix, 
                              int T);
void covarMatrix_GSL(HMM *phmm, int state, gsl_matrix * cov_matrix);
void PrintSequenceProb(FILE *fp, int T, int *O, double *vprob, double *g, 
                       double **posterior, int indexTF);
                       
/* hmmutils.c */
void ReadM(FILE *fp, HMM *phmm);
void ReadInitHMM(FILE *fp, HMM *phmm);
void ReadHMM(FILE *fp, HMM *phmm);
void PrintHMM(FILE *fp, HMM *phmm);
void getRho(HMM *phmm);
void FreeHMM(HMM *phmm);
                  
/* viterbi.c */
void Viterbi(HMM *phmm, int T, double *g, double  **alpha, double	**beta, 
             double	**gamma, double  *logprobf, double **delta, 
             int **psi, int *q, double *vprob, double *pprob, 
             double **posterior, int P, int *peakPos,
             gsl_matrix * emission_matrix);
void getPosterior_P(FILE *fpIn, FILE *fpOut,int T, int *peakPos, 
                    double **posterior, int indexTF, HMM *phmm);
int getPosterior_all_P(FILE *fpIn, FILE *fpOut, int T, int *peakPos,
                       double **posterior, HMM *phmm, int *q, double *vprob);
void getPosterior_labels(FILE *fpIn, FILE *fpOut, int T, int *q,
                         int *peakPos, double **posterior, HMM *phmm);
                                                                        
int hmmgetseed(void);
void hmmsetseed(int seed); 
double hmmgetrand(void); 

#define MAX(x,y)        ((x) > (y) ? (x) : (y))
#define MIN(x,y)        ((x) < (y) ? (x) : (y))
#define SQRT_TWO_PI 2.5066282746310002
#define D_LOG2E 1.44269504088896340736
#define TINY 1.0e-20
int MAXITERATION;
int THREAD_NUM;