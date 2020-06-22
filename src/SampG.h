#define _SAMPG_H
#ifndef _GENERAL_H
#include "General.h"
#endif

#ifndef _ERRORMAT_H
#include "ErrorMat.h"
#endif

#ifndef _SPECRAND_H
#include "SpecRand.h"
#endif

void sampG(int nsamp, int **Gobs, int **G, int *nall, int nloci, int *id, double **A, int *categories, double **E_mat, int maxall, int maxrep, int *dam, int *sire, int nind, bool estA);
