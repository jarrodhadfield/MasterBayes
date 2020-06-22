#define _SAMPS_H
#ifndef _GENERAL_H
#include "General.h"
#endif

#ifndef _ERRORMAT_H
#include "ErrorMat.h"
#endif

#ifndef _SPECRAND_H
#include "SpecRand.h"
#endif

using namespace scythe;
using namespace std;

void sampS(int *offid, int noff, Matrix<double> X_design_GS [], int *npar, int *DSuu, Matrix<double> X_design_betaSus [], Matrix<double> X_design_betaDSus [], Matrix<double> X_design_betaSs [], Matrix<double> X_design_betaDSs [], int *dam, int *sire, Matrix<double> beta, Matrix<double> us, int *ussirecat, int nusd, int nuss,  int *ndam, int *nsire, int *ntdam, int *ntsire, Matrix<int> Sires_vec [], int nind, int nmerge, int *mergeV, int *mergeUS, Matrix<double> mergeN [], std::map<int, int> Dams [], std::map<int, int> Sires [], bool checkP);
