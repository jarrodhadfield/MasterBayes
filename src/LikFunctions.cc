#include "LikFunctions.h"

double LLE_G(int **Gobs, int **G, int nloci, int *id, int nsamp, int *categories, double **LE_mat, int mtype, double **A){

	int l;			// itterates through loci
	int s;	                // itterates through duplicate samples
        int oa1;
        int oa2;
        int aa1;
        int aa2;
        int eterm;
        int cat;
        int ind;
        double LLB = 0.0;

        switch(mtype){

          case 1:
          for(s = 0; s < nsamp; s++){

            ind = id[s];
            cat = categories[s];
            for(l = 0; l < nloci; l++){
                          
              oa1 = Gobs[s][(l*2)];    
              oa2 = Gobs[s][(l*2)+1]; 
          
              aa1 = G[ind][(l*2)];
              aa2 = G[ind][(l*2)+1]; 

              if(oa1 != -999){
                if((oa1==aa1 &&  oa2==aa2) || (oa2==aa1 &&  oa1==aa2)){      
                  if(oa1==oa2){
                    LLB +=  log(A[l][oa1]*A[l][oa1]*exp(LE_mat[l][1+cat])+exp(LE_mat[l][cat]));    
                  }else{
                    LLB +=   log(2.0*A[l][oa1]*A[l][oa2]*exp(LE_mat[l][1+cat])+exp(LE_mat[l][cat]));    
                  }
                }else{
                  if(oa1==oa2){
                    LLB +=   log(A[l][oa1]*A[l][oa1]*exp(LE_mat[l][1+cat]));    
                  }else{
                    LLB +=   log(2.0*A[l][oa1]*A[l][oa2]*exp(LE_mat[l][1+cat]));    
                  }
                }               
              }     
            } 
          }
          break;


          case 2:
          for(s = 0; s < nsamp; s++){

            ind = id[s];
            cat = categories[s]*6;

            for(l = 0; l < nloci; l++){
                          
              oa1 = Gobs[s][l];    
              aa1 = G[ind][l];
              if(oa1 != -999){            
                LLB += LE_mat[l][oa1*3+aa1+cat];
              }  
            }      
          }
          break;

          case 3:

          for(s = 0; s < nsamp; s++){

            ind = id[s];
            cat = categories[s]*7;

            for(l = 0; l < nloci; l++){
                          
              oa1 = Gobs[s][(l*2)];    
              oa2 = Gobs[s][(l*2)+1];           
              aa1 = G[ind][(l*2)];
              aa2 = G[ind][(l*2)+1]; 
                          
              if(oa1 != -999){            
                if(oa1==oa2){                                      
                  eterm = 1;  
                  if(aa1==oa1 && aa2==oa1){
                    eterm = 0;
                  }
                  if(aa1!=oa1 && aa2!=oa1){
                    eterm = 2;
                  }                        
                }else{    
                  eterm = 5;
                  if(aa1==aa2){
                    eterm = 4;
                  }
                  if((aa1==oa1 && aa2==oa2) || (aa1==oa2 && aa2==oa1)){
                    eterm = 3;
                  }
                  if(aa1!=oa1 && aa2!=oa2 && aa1!=oa2 && aa2!=oa1){
                    eterm = 6;
                  }
                }
                LLB += LE_mat[l][eterm+cat];
              }  
            }      
          }
          break;
        }
                                                      
return LLB; 
}

double LLP_B(int *offid, int noff, int nind, Matrix<double> X_design_betaDus [], Matrix<double> X_design_betaSus [], Matrix<double> X_design_betaDSus [], Matrix<double> X_design_betaDs [], Matrix<double> X_design_betaSs [], Matrix<double> X_design_betaDSs [], int *npar, int *DSuu, int *dam, int *sire, Matrix<double> beta, int *ntdam, int *ntsire, int *ndam, int *nsire,  std::map<int, int> Dams [], std::map<int, int> Sires [], int nusd, int *usdamcat, int nuss, int *ussirecat, Matrix<double> us, Matrix<double> ratio[], int nmerge, int *mergeV, int *mergeUS, Matrix<double> mergeN [], int DSapprox){

 	int i;
        int dv;
        int sv;
        double dvar;
        int mvar;
        int cnt = 0;
	int d;
	int s;
        double n1;
        double n2;
	int pos_in_a;
	double ll_sum;
	double l_par;
	double ll_P_b = 0.0;
        bool damV = false;
        bool sireV = false;
        bool damsireV = false;

        Matrix<double> betaD (npar[0]+npar[1],1); 

        if((npar[0]+npar[1])>0){damV = true;}

           for(i = 0; i < (npar[0]+npar[1]); i++){   
              betaD [i] = beta[cnt];
              cnt++;
           }

        Matrix<double> betaS (npar[2]+npar[3],1); 

        if((npar[2]+npar[3])>0){sireV = true;}

           for(i = 0; i < (npar[2]+npar[3]); i++){   
              betaS [i] = beta[cnt];
              cnt++;
           }

        Matrix<double> betaDS (npar[4]+npar[5],1); 

        if((npar[4]+npar[5])>0){damsireV = true;}

           for(i = 0; i < (npar[4]+npar[5]); i++){   
              betaDS [i] = beta[cnt];
              cnt++;
           }

        Matrix<double> Dpred;
        Matrix<double> Spred;
        Matrix<double> DSpred;

	for(i = 0; i < noff; i++){   

/***************************/
/* beta reparameterisation */
/***************************/

           if(nmerge>0){
             for(dv = 0; dv < nmerge; dv++){  
               mvar = mergeV[dv];
               n1 = mergeN[dv](0,i);
               n2 = mergeN[dv](1,i);
               if(mvar<(npar[0]+npar[1])){
                 if(mergeUS[dv]==0){
                   n1 += us[usdamcat[i]];
                 }
                 if(mergeUS[dv]==1){
                   n2 += us[usdamcat[i]];
                 }
                 betaD[mvar] = beta[mvar]+log(n2/n1);
               }else{
                 if(mergeUS[dv]==0){
                   n1 += us[nusd+ussirecat[i]];
                 }
                 if(mergeUS[dv]==1){
                   n2 += us[nusd+ussirecat[i]];
                 }
                 betaS[mvar-(npar[0]+npar[1])] = beta[mvar]+log(n2/n1);
               }
             }    
           }


           if(damV){
             if(npar[0]>0 && npar[1]>0){
	         Dpred = cbind(X_design_betaDus[i], X_design_betaDs[i])*betaD;         
             }else{
               if(npar[0]>0){
                 Dpred = X_design_betaDus[i]*betaD;   
               }else{
                 Dpred = X_design_betaDs[i]*betaD;  
               } 
             }
           }

           if(sireV){
             if(npar[2]>0 && npar[3]>0){
	         Spred = cbind(X_design_betaSus[i], X_design_betaSs[i])*betaS;         
             }else{
               if(npar[2]>0){
                 Spred = X_design_betaSus[i]*betaS;   
               }else{
                 Spred = X_design_betaSs[i]*betaS;  
               } 
             }
           }

           d = dam[offid[i]];
           s = sire[offid[i]];          
           d = Dams[i][d]; 
           s = Sires[i][s];

           if(damsireV){
             if(DSapprox!=0){
               if(sireV==false){
                 Spred = Matrix<double>(ntsire[i],1);
               }
               for(sv=0; sv<ntsire[i]; sv++){ 
                 for(dv=0; dv<npar[4]; dv++){ 
                   Spred[sv] += X_design_betaDSus[i]((d*ntsire[i])+sv,dv)*betaDS[dv];
                 }
                 for(dv=0; dv<npar[5]; dv++){ 
                   Spred[sv] += X_design_betaDSs[i]((d*ntsire[i])+sv,dv)*betaDS[dv+npar[4]];
                 }
               }
               if(damV==false){
                 Dpred = Matrix<double>(ntdam[i],1);
               }
               for(dv=0; dv<ntdam[i]; dv++){ 
                 for(sv=0; sv<npar[4]; sv++){ 
                   Dpred[dv] += X_design_betaDSus[i]((dv*ntsire[i])+s,sv)*betaDS[sv];
                 }
                 for(sv=0; sv<npar[5]; sv++){ 
                   Dpred[dv] += X_design_betaDSs[i]((dv*ntsire[i])+s,sv)*betaDS[sv+npar[4]];
                 }
               }                  
             }else{
               if(npar[4]>0 && npar[5]>0){
	         DSpred = cbind(X_design_betaDSus[i], X_design_betaDSs[i])*betaDS;         
               }else{
                 if(npar[4]>0){
                    DSpred = X_design_betaDSus[i]*betaDS;   
                 }else{
                    DSpred = X_design_betaDSs[i]*betaDS;  
                 } 
               }
             }
           }  

           if(damsireV && DSapprox==0){    
             if(sireV==true || damV==true){
               cnt = 0;
               if(damV && sireV){              // D S and DS exist
                 for(dv = 0; dv < ntdam[i]; dv++){  
                   dvar = Dpred[dv];
                   for(sv = 0; sv < ntsire[i]; sv++){  
                     DSpred[cnt] += dvar+Spred[sv];
                     cnt ++;
                   }
                 }
               }else{    
                 if(damV){                         // D and DS exist
                   for(dv = 0; dv < ntdam[i]; dv++){  
                     dvar = Dpred[dv];
                     for(sv = 0; sv < ntsire[i]; sv++){  
                       DSpred[cnt] += dvar;
                       cnt ++;
                     }
                   }
                 }else{                                 // S and DS exist
                   for(dv = 0; dv < ntdam[i]; dv++){  
                     for(sv = 0; sv < ntsire[i]; sv++){  
                       DSpred[cnt] += Spred[sv];
                       cnt ++;
                     }
                   }
                 }
               }
             }

             DSpred = DSpred - (maxc(DSpred)[0] - 500.0);
             DSpred = exp(DSpred);

             pos_in_a = (d*ntsire[i])+s;
             l_par = log(DSpred[pos_in_a]);                
	     ll_sum = sumc(DSpred)[0];	     

             if(DSuu[0]==0 && DSuu[1]==0){ // dam and sire sampled

               ratio[0][i] = 0.0;
               ratio[1][i] = 0.0;

               if(nusd>0 || nuss>0){
                 if(nusd>0){
                   for(sv = 0; sv < ntsire[i]; sv++){  
                     ratio[0][i] += DSpred[((ndam[i]-1)*ntsire[i])+sv];
                   }
                 }
                 if(nuss>0){
                   for(dv = 0; dv < ntdam[i]; dv++){  
                     ratio[1][i] += DSpred[(dv*ntsire[i])+(nsire[i]-1)];
                   }
                 }
                 if(nusd>0 && nuss>0){
                   ratio[0][i] -=  DSpred[(ndam[i]-1)*ntsire[i]+nsire[i]-1];
                   ratio[1][i] -=  DSpred[(ndam[i]-1)*ntsire[i]+nsire[i]-1];
                   ll_sum -= (ratio[0][i]+ratio[1][i]+DSpred[(ndam[i]-1)*ntsire[i]+nsire[i]-1]);
                   ll_P_b += l_par - log(ll_sum + ratio[0][i]*us[usdamcat[i]]+ratio[1][i]*us[nusd+ussirecat[i]]+us[usdamcat[i]]*us[nusd+ussirecat[i]]*DSpred[(ndam[i]-1)*ntsire[i]+nsire[i]-1]);
                 }else{
                   ll_sum -= (ratio[0][i]+ratio[1][i]);
                   if(nusd>0){
                      ll_P_b += l_par - log(ll_sum + ratio[0][i]*us[usdamcat[i]]);
                   }else{
                      ll_P_b += l_par - log(ll_sum + ratio[1][i]*us[nusd+ussirecat[i]]);
                   }
                 }
               }else{
                 ll_P_b += l_par - log(ll_sum);
               }                
               if(nuss>0){
                ratio[1][i] = ll_sum/(ratio[1][i]*double(ntsire[i])-1.0); 
               }else{
                ratio[1][i] = 1.0;
               }
               if(nusd>0){
                 ratio[0][i] = ll_sum/(ratio[0][i]*double(ntdam[i])-1.0);
               }else{
                 ratio[0][i] = 1.0;
               }
             }else{                     
               if(DSuu[0]==1 && DSuu[1]==1){  // dam and sire unsampled
                 ratio[0][i] = 1.0;
                 ratio[1][i] = 1.0;
                 if(s!=(nsire[i]-1) && d!=(ndam[i]-1)){       
                   for(dv = 0; dv < ntdam[i]; dv++){  
                     ll_sum -= DSpred[(dv*ntsire[i])+(nsire[i]-1)];
                   }
                   for(sv = 0; sv < ntsire[i]; sv++){  
                     ll_sum -= DSpred[((ndam[i]-1)*ntsire[i])+sv];
                   }     
                   ll_P_b += l_par - log(ll_sum+DSpred[((ndam[i]-1)*ntsire[i])+(nsire[i]-1)]);
                 }
               }else{
                 if(DSuu[0]==1){              // dam unsampled only
                   ratio[0][i] = 1.0;
                   ratio[1][i] = 1.0;
                   if(d!=(ndam[i]-1)){
                     for(sv = 0; sv < ntsire[i]; sv++){  
                       ll_sum -= DSpred[((ndam[i]-1)*ntsire[i])+sv];
                     }      
                     ll_P_b += l_par - log(ll_sum);
                   }
                 }else{  
                   ratio[0][i] = 1.0;                    // sire unsampled only
                   ratio[1][i] = 1.0;
                   if(s!=(nsire[i]-1)){
                     for(dv = 0; dv < ntdam[i]; dv++){  
                       ll_sum -= DSpred[(dv*ntsire[i])+(nsire[i]-1)];
                     }
                     ll_P_b += l_par - log(ll_sum);
                   }
                 }
               }
             }  
           }else{
             if(sireV || (damsireV && DSapprox!=0)){                // S exists
               Spred = Spred - (maxc(Spred)[0] - 500.0);
               Spred = exp(Spred);
               l_par = log(Spred[s]);  
               ll_sum = sumc(Spred)[0];  
               if(DSuu[1]==0){  
                 if(nuss>0){
                   ll_P_b += l_par - log(ll_sum+(Spred[nsire[i]-1]*(us[nusd+ussirecat[i]]-1.0)));
                   ratio[1][i] = Spred[nsire[i]-1];
                   ratio[1][i] = (ll_sum-ratio[1][i])/((double(ntsire[i])-1.0)*ratio[1][i]);
                 }else{
                   ratio[1][i] = 1.0;
                   ll_P_b += l_par - log(ll_sum);
                 }
               }else{
                 ratio[1][i] = 1.0;
                 if(s!=(nsire[i]-1)){      
                   ll_P_b += l_par - log(ll_sum-Spred[nsire[i]-1]);
                 }
               }
             }
             if(damV || (damsireV && DSapprox!=0)){  
               Dpred = Dpred - (maxc(Dpred)[0] - 500.0);
               Dpred = exp(Dpred);
	       l_par = log(Dpred[d]);            // log(Pr(P|beta, X)) for the parent not normalised          
	       ll_sum = sumc(Dpred)[0];	    // sum(log(Pr(P|beta, X))) not normalised  
               if(DSuu[0]==0){ 
                 if(nusd>0){
                   ll_P_b += l_par - log(ll_sum+(Dpred[ndam[i]-1]*(us[usdamcat[i]]-1.0)));
                   ratio[0][i] = Dpred[ndam[i]-1];
                   ratio[0][i] = (ll_sum-ratio[0][i])/((double(ntdam[i])-1.0)*ratio[0][i]);
                 }else{
                   ratio[0][i] = 1.0;
                   ll_P_b += l_par - log(ll_sum);
                 }
               }else{
                 ratio[0][i] = 1.0;
                 if(d!=(ndam[i]-1)){                     
                   ll_P_b += l_par - log(ll_sum-Dpred[ndam[i]-1]);
                 } 
               }
             }
           }
         }	
       
return ll_P_b;
}

// CALCULATE LIKELIHOOD OF G GIVEN P
// **********************************


double LLG_P(int *offid, int noff, Matrix<double> X_design_G [], int *dam, int *sire, int *nsire, std::map<int, int> Dams [], std::map<int, int> Sires []){

	int i;
	int d;
	int s;
	int pos_in_a;
	double ll_G_P = 0.0;
	Matrix<double> a;
	for(i = 0; i < noff; i++){              
		a = X_design_G[i];      // log(Pr(P|beta, X)) not normalised     
		d = Dams[i][dam[offid[i]]];	
		s = Sires[i][sire[offid[i]]];	
		pos_in_a = (d*nsire[i])+s;
		ll_G_P += log(a[pos_in_a]);            // log(Pr(P|beta, X)) for the parent not normalised          
	}	
return ll_G_P;
}

double LLN_P(int *offid, int noff, int nind, int *ntdam, int *ntsire, int *dam, int *sire, int nusd, int *usdamcat, int nuss, int *ussirecat, Matrix<double> us, Matrix<double> ratio[]){

double llik = 0.0;


    for(int i = 0; i < noff; i++){
      if(nusd>0){
        if(dam[offid[i]]>=nind){
          llik += log(double(us[usdamcat[i]]));
          llik -= log(ratio[0][i]*double(ntdam[i]-1)+us[usdamcat[i]]);
        }else{
           llik -= log(ratio[0][i]*double(ntdam[i]-1)+us[usdamcat[i]]);
        }
      }

      if(nuss>0){
        if(sire[offid[i]]>=nind){
          llik += log(double(us[nusd+ussirecat[i]]));
          llik -= log(ratio[1][i]*double(ntsire[i]-1)+us[nusd+ussirecat[i]]);
        }else{
          llik -= log(ratio[1][i]*double(ntsire[i]-1)+us[nusd+ussirecat[i]]);
        }
      }
   }
return llik;
}                

#define LPIx2 1.837877066409345339082

double lmvnormM(Matrix<double> beta,  int nbeta, Matrix<double> mu, double log_sum_ev, Matrix<double> inv_sigma){

double llik;
 
   llik = (-(double(nbeta)*double(LPIx2)+log_sum_ev+t(beta-mu)*inv_sigma*(beta-mu))/2.0)[0];
   

return llik;
}                



