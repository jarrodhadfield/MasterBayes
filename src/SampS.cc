#include "SampS.h"

void sampS(int *offid, int noff, Matrix<double> X_design_GS [], int *npar, int *DSuu, Matrix<double> X_design_betaSus [], Matrix<double> X_design_betaDSus [], Matrix<double> X_design_betaSs [], Matrix<double> X_design_betaDSs [], int *dam, int *sire, Matrix<double> beta, Matrix<double> us, int *ussirecat, int nusd, int nuss, int *ndam, int *nsire, int *ntdam, int *ntsire, Matrix<int> Sires_vec [], int nind, int nmerge, int *mergeV, int *mergeUS, Matrix<double> mergeN [], std::map<int, int> Dams [], std::map<int, int> Sires [], bool checkP){

	int i;
	int d;
        int s;
        int j;
        int mvar;
        double n1;
        double n2;
        int cnt=0;
        int new_sire_pos;
        int dam_poss;
        double S_vec;
        double mean_vec;
        double n;
        double N;
        int *dsmatch = new int[noff];
         
         Matrix<double> betaSus (npar[2],1); 

         for(i = 0; i < npar[2]; i++){   
           betaSus [i] = beta[npar[0]+npar[1]+i];
         }

         Matrix<double> betaSs (npar[3],1); 

         for(i = 0; i < npar[3]; i++){   
           betaSs [i] = beta[npar[0]+npar[1]+npar[2]+i];
         }

         Matrix<double> betaDSus (npar[4],1);

         for(i = 0; i < npar[4]; i++){   
           betaDSus [i] = beta[npar[0]+npar[1]+npar[2]+npar[3]+i];
         }

         Matrix<double> betaDSs (npar[5],1); 

         for(i = 0; i < npar[5]; i++){   
           betaDSs [i] = beta[npar[0]+npar[1]+npar[2]+npar[3]+npar[4]+i];
         }

         Matrix<double> Spred;
	 Matrix<double> Spred_tmp;
	 Matrix<double> Spreds;
 
	for(i = 0; i < noff; i++){   
 
          dam_poss = Dams[i][dam[offid[i]]];

/***********************/
/* reparameterise beta */
/***********************/

 
           if(nmerge>0){
             for(d = 0; d < nmerge; d++){  
               mvar = mergeV[d];
               n1 = mergeN[d](0,i);
               n2 = mergeN[d](1,i);
               if(mvar<(npar[2]+npar[3])){
                 if(mergeUS[d]==0){
                   n1 += us[nusd+ussirecat[i]];
                 }
                 if(mergeUS[d]==1){
                   n2 += us[nusd+ussirecat[i]];
                 }
                 if(mvar<(npar[0]+npar[1]+npar[2])){
                   betaSus[mvar-(npar[0]+npar[1])] = beta[mvar]+log(n2/n1);
                 }else{
                   betaSs[mvar-(npar[0]+npar[1]+npar[2])] = beta[mvar]+log(n2/n1);
                 }
               }
             }    
           }

/***********************/
/* sample missing data */
/***********************/

              if(npar[2]!=0){   
                 Spred_tmp = X_design_betaSus[i]*betaSus;
              }else{
                if(npar[4]!=0 && dam[offid[i]]<nind){
                  Spred_tmp = Matrix<double>(ntsire[i],1);
                }
              }

              if(npar[4]!=0 && dam[offid[i]]<nind){ 
                for(d=0; d<ntsire[i]; d++){ 
                  for(s=0; s<npar[4]; s++){ 
                    Spred_tmp[d] += X_design_betaDSus[i]((dam_poss*ntsire[i])+d,s)*betaDSus[s];
                  }
                }
              }
              if(npar[2]!=0 || (npar[4]!=0 && dam[offid[i]]<nind && DSuu[1]==1)){ 
                Spred_tmp = Spred_tmp - (maxc(Spred_tmp)[0]-100.0);
                Spred = exp(Spred_tmp);
                mean_vec = meanc(Spred)[0];
                n = double(ntsire[i]-1);
                mean_vec -= Spred[nsire[i]-1]/(n+1.0);  
                if(n<2){
                  if(n==0){Spred[nsire[i]-1] = 1.0;}    
                  if(n==1){Spred[nsire[i]-1] = mean_vec*2.0;}
                }else{
                  mean_vec *= (n+1.0)/n;                  // mean linear predictor of sampled dams
                  S_vec = varc(Spred)[0];                 // variance of the vector = sample variance of linear -
                  N = n+us[ussirecat[i]+nusd];                  // predictor of sampled dams
                  S_vec *= N/(n*(N-n));
                  Spred[nsire[i]-1] = rnorm(mean_vec, sqrt(S_vec));
                  if(Spred[nsire[i]-1]<0.0){               // for those instances where the linear predictor goes negative
                    Spred[nsire[i]-1]=1e-100;
                  }
                }
                Spred_tmp[nsire[i]-1] = log(Spred[nsire[i]-1]);      
              }

/* combine complete and incomplete design matrices */

              Spreds = Matrix<double>(nsire[i],1);

              if(npar[3]!=0){   
                for(d=0; d<nsire[i]; d++){ 
                  for(s=0; s<npar[3]; s++){ 
                     Spreds[d] += X_design_betaSs[i](d,s)*betaSs[s];
                  }
                }
              }

              if(npar[5]!=0){ 
                for(d=0; d<nsire[i]; d++){ 
                  for(s=0; s<npar[5]; s++){ 
                    Spreds[d] += X_design_betaDSs[i]((dam_poss*ntsire[i])+d,s)*betaDSs[s];
                  }
                }
              } 
 
              if(npar[2]!=0 || (npar[4]!=0 && dam[offid[i]]<nind)){ 
                if(npar[3]!=0 || npar[5]!=0 || DSuu[1]==0){  
                  for(d=0; d<nsire[i]; d++){
                    Spreds[d] += Spred_tmp[d];
                  }
                }else{
                  for(d=0; d<nsire[i]; d++){
                    Spreds[d] = Spred[d];  // note that it is already exponentiated in this case
                  }
                }
              }

              if((npar[2]!=0 || (npar[4]!=0 && dam[offid[i]]<nind)) && npar[3]==0 && npar[5]==0 && DSuu[1]==1){ 
                for(d=0; d<nsire[i]; d++){
                  Spreds[d] *= X_design_GS[i][d]; 
                }
              }else{
                if((npar[2]+npar[3]+(npar[4]*int(dam[offid[i]]<nind))+npar[5])>0){ 
                  Spreds = Spreds - (maxc(Spreds)[0]-650.0);
                  Spreds = exp(Spreds);
                  for(d=0; d<nsire[i]; d++){
                    Spreds[d] *= X_design_GS[i][d]; 
                  }
                }else{
                   Spreds = X_design_GS[i];  
                }
              }
 
              if(nuss>0){
                Spreds[(nsire[i]-1)] *= us[ussirecat[i]+nusd];
              }
 
              if(checkP){
                dsmatch[0] = offid[i];
                cnt=1;
                for(d=0; d<cnt; d++){
                  for(s=0; s<noff; s++){
                    if(dam[offid[s]]==dsmatch[d] || sire[offid[s]]==dsmatch[d]){
                      dsmatch[cnt]=offid[s];
                      for(j=0; j<cnt; j++){
                        if(offid[s]==dsmatch[j]){
                          cnt--;
                          break;
                        }
                      }
                      cnt++;
                      for(j=0; j<nsire[i]; j++){
                        if(Sires_vec[i][j]==offid[s]){
                          Spreds[j]=0.0;
                          break;
                        }  
                      }
                    }
                  }  
                }
              }
             
              new_sire_pos = rmultinom_size1M(Spreds, nsire[i]);      // samples a parental combination
              sire[offid[i]] = Sires_vec[i][new_sire_pos];     // puts new dam.id into pedigree       
 
        }
       delete [] dsmatch;
}



