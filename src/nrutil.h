/*
**  File:   nrutil.c
**  codes were borrowed and modified from Kanungo's
**  original HMM program.
**  Tapas Kanungo, "UMDHMM: Hidden Markov Model Toolkit," in "Extended Finite State Models of Language," A. Kornai (editor), Cambridge University Press, 1999. http://www.kanungo.com/software/software.html.
**  Purpose: Memory allocation routines borrowed from the
**    book "Numerical Recipes" by Press, Flannery, Teukolsky,
**    and Vetterling.
**   state sequence and probablity of observing a sequence
**   given the model.
**
*/

float *vector();
float **matrix();
float **convert_matrix();
double *dvector();
double **dmatrix();
int *ivector();
int **imatrix();
float **submatrix();
void free_vector();
void free_dvector();
void free_ivector();
void free_matrix();
void free_dmatrix();
void free_imatrix();
void free_submatrix();
void free_convert_matrix();
void nrerror();
void copyRow();
void copyMatrix();
void printMatrix();
void printVector();
void printfMatrix();
double listInsertnMax();
