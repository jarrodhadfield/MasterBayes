#define _SPECRAND_H
#ifndef _GENERAL_H
#include "General.h"
#endif

int rmultinom_size1(double *prob, int mult_size);
int rmultinom_size1M(Matrix<double> prob, int mult_size);
void rdirichlet(double *count, int K, double *freq); 
Matrix<double> rmvnormM(Matrix<double> beta, Matrix<double> cholD, int K);
