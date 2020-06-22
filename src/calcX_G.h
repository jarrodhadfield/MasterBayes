#define _CALCX_G_H
#ifndef _GENERAL_H
#include "General.h"
#endif


void calcX_G(Matrix<double> X_design_G [], int *offid, int noff , int *ndam, int *nsire, int nind, Matrix<int> Dams_vec [], Matrix<int> Sires_vec [], int **G, int nloci, double **A, int mtype);

void calcX_Gcervus(Matrix<double> X_design_G [], int *offid, int noff , int *ndam, int *nsire, int nind, Matrix<int> Dams_vec [], Matrix<int> Sires_vec [], int **G, int nloci, double **A, double **E_mat, int mtype, int *nall);

void calcX_Gmm(Matrix<int> mmD [], Matrix<int> mmS [], int *offid, int noff , int *ndam, int *nsire, int nind, Matrix<int> Dams_vec [], Matrix<int> Sires_vec[],int **G, int nloci);

void calcX_GD(Matrix<double> X_design_GD [], int *offid, int noff , int *ndam, int nind, Matrix<int> Dams_vec [], int *sire, int **G, int nloci, double **A, int mtype);

void calcX_GS(Matrix<double> X_design_GD [], int *offid, int noff , int *nsire, int nind, Matrix<int> Sires_vec [], int *dam, int **G, int nloci, double **A, int mtype);

void calcX_GcervusS(Matrix<double> X_design_GS [], Matrix<double> X_design_G [], int *offid, int noff , int *ndam, int *nsire, std::map<int, int> Dams [], int *dam);

void calcX_GcervusD(Matrix<double> X_design_GD [], Matrix<double> X_design_G [], int *offid, int noff , int *ndam, int *nsire, std::map<int, int> Sires [], int *sire);

