#define _LIKFUNCTIONS_H
#ifndef _GENERAL_H
#include "General.h"
#endif

double LLE_G(int **Gobs, int **G, int nloci, int *id, int nsamp, int *categories, double **LE_mat, int mtype, double **A);

double LLP_B(int *offid, int noff, int nind, Matrix<double> X_design_betaDus [], Matrix<double> X_design_betaSus [], Matrix<double> X_design_betaDSus [], Matrix<double> X_design_betaDs [], Matrix<double> X_design_betaSs [], Matrix<double> X_design_betaDSs [], int *npar, int *DSuu, int *dam, int *sire, Matrix<double> beta, int *ntdam, int *ntsire, int *ndam, int *nsire,  std::map<int, int> Dams [], std::map<int, int> Sires [], int nusd, int *usdamcat, int nuss, int *ussirecat, Matrix<double> us,  Matrix<double> ratio[], int nmerge, int *mergeV, int *mergeUS, Matrix<double> mergeN [], int DSapprox);

double LLG_P(int *offid, int noff, Matrix<double> X_design_G [], int *dam, int *sire, int *nsire, std::map<int, int> Dams [], std::map<int, int> Sires []);

double LLN_P(int *offid, int noff, int nind, int *ndam, int *nsire, int *dam, int *sire, int nusd, int *usdamcat, int nuss, int *ussirecat, Matrix<double> us, Matrix<double> ratio[]);

double lmvnormM(Matrix<double> beta,  int nbeta, Matrix<double> mu, double log_sum_ev, Matrix<double> inv_sigma);
