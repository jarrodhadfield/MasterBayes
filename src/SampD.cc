#include "SampD.h"

void sampD(int *offid, int noff, Matrix<double> X_design_GD [], int *npar, int *DSuu, Matrix<double> X_design_betaDus [], Matrix<double> X_design_betaDSus [], Matrix<double> X_design_betaDs [], Matrix<double> X_design_betaDSs [], int *dam, int *sire, Matrix<double> beta, Matrix<double> us, int *usdamcat, int nusd, int *ndam, int *nsire, int *ntdam, int *ntsire, Matrix<int> Dams_vec [], int nind, int nmerge, int *mergeV, int *mergeUS, Matrix<double> mergeN [], std::map<int, int> Dams [], std::map<int, int> Sires [], bool checkP){

	int i;
	int d;
        int s;
        int j;
        int mvar;
        double n1;
        double n2;
        int cnt=0;
        int new_dam_pos;
        int sire_poss;
        double S_vec;
        double mean_vec;
        double n;
        double N;
        int *dsmatch = new int[noff];
         
         Matrix<double> betaDus (npar[0],1); 

         for(i = 0; i < npar[0]; i++){   
           betaDus [i] = beta[i];
         }

         Matrix<double> betaDs (npar[1],1); 

         for(i = 0; i < npar[1]; i++){   
           betaDs [i] = beta[i+npar[0]];
         }

         Matrix<double> betaDSus (npar[4],1);

         for(i = 0; i < npar[4]; i++){   
           betaDSus [i] = beta[npar[0]+npar[1]+npar[2]+npar[3]+i];
         }

         Matrix<double> betaDSs (npar[5],1); 

         for(i = 0; i < npar[5]; i++){   
           betaDSs [i] = beta[npar[0]+npar[1]+npar[2]+npar[3]+npar[4]+i];
         }

         Matrix<double> Dpred;
	 Matrix<double> Dpred_tmp;
	 Matrix<double> Dpreds;
 
	for(i = 0; i < noff; i++){   

          sire_poss = Sires[i][sire[offid[i]]];

/***********************/
/* reparameterise beta */
/***********************/

           if(nmerge>0){
             for(d = 0; d < nmerge; d++){  
               mvar = mergeV[d];
               n1 = mergeN[d](0,i);
               n2 = mergeN[d](1,i);
                if(mvar<(npar[0]+npar[1])){
                 if(mergeUS[d]==0){
                   n1 += us[usdamcat[i]];
                 }
                 if(mergeUS[d]==1){
                   n2 += us[usdamcat[i]];
                 }
                 if(mvar<npar[0]){
                   betaDus[mvar] = beta[mvar]+log(n2/n1);
                 }else{
                   betaDs[mvar-npar[0]] = beta[mvar]+log(n2/n1);
                 }
               }
             }    
           }

/***********************/
/* sample missing data */
/***********************/

              if(npar[0]!=0){   
                 Dpred_tmp = X_design_betaDus[i]*betaDus;
              }else{
                if(npar[4]!=0 && sire[offid[i]]<nind){
                  Dpred_tmp = Matrix<double>(ntdam[i],1);
                }
              }
          
              if(npar[4]!=0 && sire[offid[i]]<nind){ 
                for(d=0; d<ntdam[i]; d++){ 
                  for(s=0; s<npar[4]; s++){ 
                    Dpred_tmp[d] += X_design_betaDSus[i]((d*ntsire[i])+sire_poss,s)*betaDSus[s];
                  }
                }
              }

              if(npar[0]!=0 || (npar[4]!=0 && sire[offid[i]]<nind && DSuu[0]==1)){ 
                Dpred_tmp = Dpred_tmp - (maxc(Dpred_tmp)[0]-100.0);
                Dpred = exp(Dpred_tmp);
                mean_vec = meanc(Dpred)[0];
                n = double(ntdam[i]-1);
                mean_vec -= Dpred[ndam[i]-1]/(n+1.0);  
                if(n<2){
                  if(n==0){Dpred[ndam[i]-1] = 1.0;}    
                  if(n==1){Dpred[ndam[i]-1] = mean_vec*2.0;}
                }else{
                  mean_vec *= (n+1.0)/n;                  // mean linear predictor of sampled dams
                  S_vec = varc(Dpred)[0];                 // variance of the vector = sample variance of linear -
                  N = n+us[usdamcat[i]];                  // predictor of sampled dams
                  S_vec *= N/(n*(N-n));
                  Dpred[ndam[i]-1] = rnorm(mean_vec, sqrt(S_vec));
                  if(Dpred[ndam[i]-1]<0.0){               // for those instances where the linear predictor goes negative
                    Dpred[ndam[i]-1]=1e-100;
                  }
                }
                Dpred_tmp[ndam[i]-1] = log(Dpred[ndam[i]-1]);      
              }
              
/* combine complete and incomplete design matrices */

              Dpreds = Matrix<double>(ndam[i],1);

              if(npar[1]!=0){   
                for(d=0; d<ndam[i]; d++){ 
                  for(s=0; s<npar[1]; s++){ 
                     Dpreds[d] += X_design_betaDs[i](d,s)*betaDs[s];
                  }
                }
              }

              if(npar[5]!=0){ 
                for(d=0; d<ndam[i]; d++){ 
                  for(s=0; s<npar[5]; s++){ 
                    Dpreds[d] += X_design_betaDSs[i]((d*ntsire[i])+sire_poss,s)*betaDSs[s];
                  }
                }
              } 

              if(npar[0]!=0 || (npar[4]!=0 && sire[offid[i]]<nind)){ 
                if(npar[1]!=0 || npar[5]!=0 || DSuu[0]==0){  
                  for(d=0; d<ndam[i]; d++){
                    Dpreds[d] += Dpred_tmp[d];
                  }
                }else{
                  for(d=0; d<ndam[i]; d++){
                    Dpreds[d] = Dpred[d]; // note that it is already exponentiated in this case
                  } 
                }
              }

              if((npar[0]!=0 || (npar[4]!=0 && sire[offid[i]]<nind)) && npar[1]==0 && npar[5]==0 && DSuu[0]==1){ 
                for(d=0; d<ndam[i]; d++){
                  Dpreds[d] *= X_design_GD[i][d]; 
                }
              }else{
                if((npar[0]+npar[1]+(npar[4]*int(sire[offid[i]]<nind))+npar[5])>0){ 
                   Dpreds = Dpreds - (maxc(Dpreds)[0]-100.0);
                   Dpreds = exp(Dpreds);
                   for(d=0; d<ndam[i]; d++){
                     Dpreds[d] *= X_design_GD[i][d]; 
                   }
                }else{
                   Dpreds = X_design_GD[i];  
                }
              }

              if(nusd>0){
                Dpreds[(ndam[i]-1)] *= us[usdamcat[i]];
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
                     for(j=0; j<ndam[i]; j++){
                       if(Dams_vec[i][j]==offid[s]){
                         Dpreds[j]=0.0;
                         break;
                       }
                     }
                   }
                 }
               }  
             }

            new_dam_pos = rmultinom_size1M(Dpreds, ndam[i]);      // samples a parental combination
		
            dam[offid[i]] = Dams_vec[i][new_dam_pos];     // puts new dam.id into pedigree       
 
        }

        delete [] dsmatch;
}



