#include "MCMCped.h" 

extern "C"{  

using namespace std;

void fillXG(
	int *nindP,		 // number of individuals sampled
        int *noffP,              // number of non-base (offspring) individuals
        int *ndamP,              // number of candidate dams per offspring
        int *nsireP,             // number of candidate sires per offspring
        int *nlociP,		 // number of loci
        int *nallP,              // number of alleles per locus
        int *maxallP,            // number of alleles at most polymorhic locus
        int *offidP,             // offspring id
        int *damidP,             // candidate dam id's for each offspring
        int *sireidP,		 // candidate sire id's for each offspring	
        double *X_design_GP,     // Mendelian transition probabilities dam and sire sampled
	double *AP,	         // starting allele frequencies
	double *E1P,	         // starting values of E1 and E2
	double *E2P,	         // starting values of E1 and E2
        int *GP,              // starting true genotypes    
        int *mtypeP
){         // if TRUE joint posterior distribution of P is written (default = marginal)

// pointers to single variables are redefined

int 	nind = nindP[0],
        noff = noffP[0],  	
	nloci = nlociP[0],      
        maxall = maxallP[0],
        *nall = nallP,
        i,
        index=0,
        mtype = mtypeP[0],
        ncat = 1;

// declare some variable sized arrays

int     *pG,
        **G;
double  *pA,
        **A,
        *pE_mat,
        **E_mat;
 
pG = new(nothrow) int [(1+int(mtype==1 || mtype==3))*nind*nloci];
if(pG==NULL)
{
 error("NO MEMORY for G\n");}
G = new(nothrow) int* [nind];
if(G==NULL)
{
 error("NO MEMORY for G\n");}
pA = new(nothrow) double [nloci*maxall];
if(pA==NULL)
{
 error("NO MEMORY for A\n");}
A = new(nothrow) double* [nloci];
if(A==NULL)
{
 error("NO MEMORY for A\n");
}
pE_mat = new(nothrow) double [ncat*(4+3*int(mtype==1 || mtype==3)+2*int(mtype==2))*nloci];
if(pE_mat==NULL){
 error("NO MEMORY for E_mat\n");
}
E_mat = new(nothrow) double* [nloci];
if(E_mat==NULL){
 error("NO MEMORY for E_mat\n");
}

for (i=0; i<nind; ++i){
G[i] = &pG[index];
index  += ((1+int(mtype==1 || mtype==3))*nloci);
}

index = 0;

for (i=0; i<nloci; ++i){
A[i] = &pA[index];
index  += maxall;
}

index = 0;
for (i=0; i<nloci; ++i){
E_mat[i] = &pE_mat[index];
index  += (ncat*(4+3*int(mtype==1 || mtype==3)+2*int(mtype==2)));
}
        
        Matrix<double> *X_design_G = new Matrix<double>[noff];

	map<int, int> *Dams = new map<int, int>[noff];     // two way indexing vectors
	map<int, int> *Sires = new map<int, int>[noff];    // map[i][dam_id] = n           dam_id is the n^th mother of the i^th individual
        Matrix<int> *Dams_vec = new Matrix<int>[noff];   // Matrix[i][n] = dam_id        the n^th mother of the i^th individul is dam.id 
	Matrix<int> *Sires_vec = new Matrix<int>[noff];

         read_G(GP, nind, nloci, G, mtype);
         read_A(AP, nloci, A, nall);
         read_stP(noff, ndamP, damidP, nsireP, sireidP,Dams,Sires,Dams_vec,Sires_vec); 
         Error_Mat(E1P[0], E2P[0], E_mat, ncat, nall, nloci, false, false, mtype);
         calcX_Gcervus(X_design_G, offidP, noff , ndamP, nsireP, nind, Dams_vec, Sires_vec, G, nloci, A, E_mat, mtype, nall);
         int cnt_ds=0;
         int p;

	 for(i = 0; i <noff ; i++){
           for(p = 0; p <(ndamP[i]*nsireP[i]); p++){
            X_design_GP[cnt_ds] = X_design_G[i][p];
            cnt_ds++;
           }
       }

   delete [] pG;
   delete [] G;
   delete [] pA;
   delete [] A;
   delete [] pE_mat;
   delete [] E_mat;
   delete [] Dams;
   delete [] Sires;
   delete [] Dams_vec;
   delete [] Sires_vec;
   delete [] X_design_G;

  }
}
