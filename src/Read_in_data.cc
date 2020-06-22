#include "Read_in_data.h" 

void read_G(int *GP, int nsamp, int nloci, int **G, int mtype){

	int i;			// itterates through individuals
	int l;			// itterates through loci
	int records = 0;	// itterates through genotypes

        if(mtype==1 || mtype==3){
	  for(i = 0; i < nsamp; i++){
            for(l = 0; l < nloci; l++){	
	      G[i][l*2] = GP[records];
              records ++;        
              G[i][(l*2)+1] = GP[records];
              records ++;
            }
          }
        }else{		
	  for(i = 0; i < nsamp; i++){
            for(l = 0; l < nloci; l++){	
	      G[i][l] = GP[records];
              records ++;        
            }
          }
       }
}
 
void read_A(double *AP, int nloci, double **A, int *nall){

	int l;			// itterates through loci
        int a;
	int records = 0;	// itterates through genotypes
		
        records = 0;   
        for(l = 0; l < nloci; l++){
          for(a = 0; a < nall[l]; a++){                        
	    A[l][a] = AP[records];
            records ++;                                                
	  }
        }			
}

void read_X_beta(int noff, int *ndam, int *nsire, int *npar, double *X_design_betaDusP, double *X_design_betaSusP, double *X_design_betaDSusP,  double *X_design_betaDsP, double *X_design_betaSsP, double *X_design_betaDSsP, Matrix<double> X_design_betaDus [], Matrix<double> X_design_betaSus [], Matrix<double> X_design_betaDSus [], Matrix<double> X_design_betaDs [], Matrix<double> X_design_betaSs [], Matrix<double> X_design_betaDSs []){

	int i;
	int j;
	int cntdus = 0;
	int cntds = 0;
	int cntsus = 0;
	int cntss = 0;
	int cntdsus = 0;
	int cntdss = 0;


       for(i=0; i<noff; i++){                
         if(npar[0]>0){
	   X_design_betaDus [i]  =  ones<double>(ndam[i], npar[0]); 
           for(j=0; j<(npar[0]*ndam[i]); j++){
             X_design_betaDus[i][j] = X_design_betaDusP[cntdus];
             cntdus++;
           }
         }
         if(npar[1]>0){
	   X_design_betaDs [i]  =  ones<double>(ndam[i], npar[1]); 
           for(j=0; j<(npar[1]*ndam[i]); j++){
             X_design_betaDs[i][j] = X_design_betaDsP[cntds];
             cntds++;
           }
         }
         if(npar[2]>0){
	   X_design_betaSus [i]  =  ones<double>(nsire[i],npar[2]); 
           for(j=0; j<(npar[2]*nsire[i]); j++){
             X_design_betaSus[i][j] = X_design_betaSusP[cntsus];
             cntsus++;
           }
         }
         if(npar[3]>0){
	   X_design_betaSs [i]  =  ones<double>(nsire[i],npar[3]); 
           for(j=0; j<(npar[3]*nsire[i]); j++){
             X_design_betaSs[i][j] = X_design_betaSsP[cntss];
             cntss++;
           }
          }
         if(npar[4]>0){
	   X_design_betaDSus [i]  =  ones<double>(ndam[i]*nsire[i], npar[4]); 
           for(j=0; j<(npar[4]*ndam[i]*nsire[i]); j++){
             X_design_betaDSus[i][j] = X_design_betaDSusP[cntdsus];
             cntdsus++;
           }
         }
         if(npar[5]>0){
	   X_design_betaDSs [i]  =  ones<double>(ndam[i]*nsire[i], npar[5]); 
           for(j=0; j<(npar[5]*ndam[i]*nsire[i]); j++){
             X_design_betaDSs[i][j] = X_design_betaDSsP[cntdss];
             cntdss++;
           }
         }
      }
}
     
void read_stP(int noff, int *ndam, int *damid, int *nsire, int *sireid, map<int,int> Dams [], map<int,int> Sires [], Matrix<int> Dams_vec [], Matrix<int> Sires_vec []){

        int i;
        int e_i = 0;
	int d;
	int s;
	int d_cum_st = 0;
	int s_cum_st = 0;
	int d_cum_end;
	int s_cum_end;

	for(i=0; i<noff; i++){
	  d_cum_end = d_cum_st + ndam[i];
	  s_cum_end = s_cum_st + nsire[i];
	  e_i = 0;

		Dams_vec[i] = ones<int>(ndam[i],1);
		Sires_vec[i] = ones<int>(nsire[i],1);
		for(s=s_cum_st; s<s_cum_end; s++){
		Sires [i] [sireid[s]] = e_i;
		Sires_vec [i][e_i] = sireid[s];
                e_i ++;
		}
	  e_i = 0;
			for(d=d_cum_st; d<d_cum_end; d++){
			  Dams [i] [damid[d]] = e_i;
			  Dams_vec [i][e_i] = damid[d];
          e_i ++;
			}

	d_cum_st = d_cum_end;
	s_cum_st = s_cum_end;

	}

}


