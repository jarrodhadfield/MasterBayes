#include "MCMCped.h" 

extern "C"{  

using namespace std;

void mismatches(
	int *nindP,		 // number of individuals sampled
        int *noffP,              // number of non-base (offspring) individuals
        int *ndamP,              // number of candidate dams per offspring
        int *nsireP,             // number of candidate sires per offspring
        int *nlociP,		 // number of loci
        int *offidP,             // offspring id
        int *damidP,             // candidate dam id's for each offspring
        int *sireidP,		 // candidate sire id's for each offspring	
        int *mmDP,             // number of misamtches per dam
        int *mmSP,             // number of misamtches per sire     
        int *st_GP              // starting true genotypes 
){         
// pointers to single variables are redefined

int 	nind = nindP[0],
        noff = noffP[0],  	
	nloci = nlociP[0],
        i,
        l=0;

// declare some variable sized arrays

int     *pG,
        **G;
 
pG = new(nothrow) int [2*nind*nloci];
if(pG==NULL)
{
 error("NO MEMORY for G\n");
}
G = new(nothrow) int* [nind];
if(G==NULL)
{
 error("NO MEMORY for G\n");
}

for (i=0; i<nind; ++i){
G[i] = &pG[l];
l  += (2*nloci);
}

     
        Matrix<int> *mmD = new Matrix<int>[noff];
        Matrix<int> *mmS = new Matrix<int>[noff];

	map<int, int> *Dams = new map<int, int>[noff];     // two way indexing vectors
	map<int, int> *Sires = new map<int, int>[noff];    // map[i][dam_id] = n           dam_id is the n^th mother of the i^th individual
        Matrix<int> *Dams_vec = new Matrix<int>[noff];   // Matrix[i][n] = dam_id        the n^th mother of the i^th individul is dam.id 
	Matrix<int> *Sires_vec = new Matrix<int>[noff];
	
        read_G(st_GP, nind, nloci, G, 1);

        read_stP(noff, ndamP, damidP, nsireP, sireidP,Dams,Sires,Dams_vec,Sires_vec);
         
        calcX_Gmm(mmD, mmS, offidP, noff , ndamP, nsireP, nind, Dams_vec, Sires_vec, G, nloci);
 
         int cnt=0;
  
	 for(i = 0; i <noff ; i++){
           for(int d = 0; d <ndamP[i]; d++){
            mmDP[cnt] = mmD[i][d];
            cnt++;
           }
         }
         
         cnt = 0;

 	 for(i = 0; i <noff ; i++){
           for(int s = 0; s <nsireP[i]; s++){
            mmSP[cnt] = mmS[i][s];
            cnt++;
           }
         }

       delete [] mmD;
       delete [] mmS;
       delete [] Dams;
       delete [] Sires;
       delete [] Dams_vec;
       delete [] Sires_vec;
       delete [] pG;
       delete [] G;


}
}
