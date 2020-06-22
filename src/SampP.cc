#include "SampP.h"

void sampP(int *offid, int noff, Matrix<double> X_design_G [], int *npar, int *DSuu, Matrix<double> X_design_betaDus [], Matrix<double> X_design_betaSus [], Matrix<double> X_design_betaDSus [], Matrix<double> X_design_betaDs [], Matrix<double> X_design_betaSs [], Matrix<double> X_design_betaDSs [], int *dam, int *sire, Matrix<double> beta, Matrix<double> us, int *usdamcat, int *ussirecat, int nusd, int nuss, int *ndam, int *nsire, int *ntdam, int *ntsire, Matrix<int> Dams_vec [], Matrix<int> Sires_vec [], int nind, int nmerge, int *mergeV, int *mergeUS, Matrix<double> mergeN [], bool checkP){

	int i;
	int d;
        int s;
        int b;
        int j;
        int mvar;
        double n1;
        double n2;
        int cnt=0;
        int newpar;
        int new_dam_pos;
        int new_sire_pos;
        double S_vec;
        double mean_vec;
        double n;
        double N;
        double tot_mean;
        bool damV;
        bool sireV;
        bool damsireV;
        int *dsmatch = new int[noff];
        bool mother_complete = TRUE;
        bool father_complete = TRUE;

        if(npar[0]!=0 || (npar[4]!=0 && DSuu[0]==1)){mother_complete=FALSE;}
        if(npar[2]!=0 || (npar[4]!=0 && DSuu[1]==1)){father_complete=FALSE;}

         Matrix<double> betaDus (npar[0],1); 

           for(i = 0; i < npar[0]; i++){   
              betaDus [i] = beta[cnt];
              cnt++;
           }
           
         Matrix<double> betaDs (npar[1],1); 

           for(i = 0; i < npar[1]; i++){   
              betaDs [i] = beta[cnt];
              cnt++;
           }

         if((npar[0]+npar[1])>0){damV=true;}else{damV=false;}

         Matrix<double> betaSus (npar[2],1); 
           for(i = 0; i < npar[2]; i++){   
              betaSus [i] = beta[cnt];
              cnt++;
           }

         Matrix<double> betaSs (npar[3],1); 
           for(i = 0; i < npar[3]; i++){   
              betaSs [i] = beta[cnt];
              cnt++;
           }

         if((npar[2]+npar[3])>0){sireV=true;}else{sireV=false;}

         Matrix<double> betaDSus (npar[4],1); 
           for(i = 0; i < npar[4]; i++){   
              betaDSus [i] = beta[cnt];
              cnt++;
           }

         Matrix<double> betaDSs (npar[5],1); 
           for(i = 0; i < npar[5]; i++){   
              betaDSs [i] = beta[cnt];
              cnt++;
           }

         if((npar[4]+npar[5])>0){damsireV=true;}else{damsireV=false;}

	Matrix<double> Dpred;
	Matrix<double> Spred;	
        Matrix<double> DSpred;
	Matrix<double> Dpred_tmp;
	Matrix<double> Spred_tmp;	
        Matrix<double> DSpred_tmp;
	Matrix<double> Dpreds;
	Matrix<double> Spreds;	
        Matrix<double> DSpreds;

	for(i = 0; i < noff; i++){   

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
               }else{
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

/* combine complete and incomplete design matrices */

           if(npar[4]!=0){
             DSpred_tmp = X_design_betaDSus[i]*betaDSus;
             if(npar[5]!=0){
               DSpred_tmp = DSpred_tmp + X_design_betaDSs[i]*betaDSs;
             }
           }else{
             if(npar[5]!=0){
               if(mother_complete && father_complete){  // Both parents have complete information
                 DSpred_tmp = Matrix<double>(nsire[i]*ndam[i],1);
                 for(d=0; d<ndam[i]; d++){
                   for(s=0; s<nsire[i]; s++){            
                     for(b=0; b<npar[5]; b++){
                       DSpred_tmp[(d*nsire[i])+s] += X_design_betaDSs[i]((d*ntsire[i])+s, b)*betaDSs[b]; 
                     }
                   }
                 }
               }else{   
                 DSpred_tmp = X_design_betaDSs[i]*betaDSs;
               }
             }
           }

           if(npar[0]!=0){
             Dpred_tmp = X_design_betaDus[i]*betaDus;
             if(npar[1]!=0){
               Dpred_tmp = Dpred_tmp + X_design_betaDs[i]*betaDs;
             }
           }else{
             if(npar[1]!=0){
               if(mother_complete){  // Mothers have complete information
                 Dpred_tmp = Matrix<double>(ndam[i],1);
                 for(d=0; d<ndam[i]; d++){
                   for(b=0; b<npar[1]; b++){
                     Dpred_tmp[d] += X_design_betaDs[i](d,b)*betaDs[b]; 
                   }
                 }  
               }else{
                 Dpred_tmp = X_design_betaDs[i]*betaDs;
               }
             }
           }

           
           if(npar[2]!=0){
             Spred_tmp = X_design_betaSus[i]*betaSus;
             if(npar[3]!=0){
               Spred_tmp =  Spred_tmp + X_design_betaSs[i]*betaSs;
             }
           }else{
             if(npar[3]!=0){
               if(father_complete){  // Fathers have complete information
                 Spred_tmp = Matrix<double>(nsire[i],1);
                 for(d=0; d<nsire[i]; d++){
                   for(b=0; b<npar[3]; b++){
                     Spred_tmp[d] += X_design_betaSs[i](d,b)*betaSs[b]; 
                   }
                 }  
               }else{
                 Spred_tmp = X_design_betaSs[i]*betaSs;
               }
             }
           }

/****************************/
/* combine with DS variable */
/****************************/
           
           if(damsireV){ 
             if(damV){  
               if(mother_complete && father_complete){  
                 for(d=0; d<ndam[i]; d++){
                   for(s=0; s<nsire[i]; s++){            
                     DSpred_tmp[(d*nsire[i])+s] += Dpred_tmp[d]; 
                   }
                 }
               }else{
                 for(d=0; d<ntdam[i]; d++){
                   for(s=0; s<ntsire[i]; s++){            
                     DSpred_tmp[(d*ntsire[i])+s] += Dpred_tmp[d]; 
                   }
                 }
               }
             }
             if(sireV){
               if(mother_complete && father_complete){   
                 for(d=0; d<ndam[i]; d++){
                   for(s=0; s<nsire[i]; s++){            
                     DSpred_tmp[(d*nsire[i])+s] += Spred_tmp[s]; 
                   }
                 }
               }else{
                 for(d=0; d<ntdam[i]; d++){
                   for(s=0; s<ntsire[i]; s++){            
                     DSpred_tmp[(d*ntsire[i])+s] += Spred_tmp[s]; 
                   }
                 }
               } 
             }
           }

/****************/
/* exponentiate */
/****************/

           if(damsireV){
             DSpred_tmp = DSpred_tmp - (maxc(DSpred_tmp)[0]-100.0);
             if(mother_complete==FALSE || father_complete==FALSE){
               DSpred = exp(DSpred_tmp);
             }
           }else{
             if(damV){
               Dpred_tmp = Dpred_tmp - (maxc(Dpred_tmp)[0]-100.0);
               if(mother_complete==FALSE){
                 Dpred = exp(Dpred_tmp);
               }
             }
             if(sireV){
               Spred_tmp = Spred_tmp - (maxc(Spred_tmp)[0]-100.0);
               if(father_complete==FALSE){
                 Spred = exp(Spred_tmp);
               }
             }
           }

/************************/
/* sample missing data  */
/************************/

           if(damsireV){
             tot_mean=0.0;    /*sample missing for DS variables */
             if(mother_complete==FALSE){ 
               for(s=0; s<nsire[i]; s++){
                 mean_vec = 0.0;
                 S_vec = 0.0;   
                 for(d=0; d<ntdam[i]; d++){                   
                   mean_vec += DSpred[(d*ntsire[i])+s];
                 }
                 mean_vec -= DSpred[((ndam[i]-1)*ntsire[i])+s];
                 if(ntdam[i]<3){
                   if(ntdam[i]==1){DSpred[((ndam[i]-1)*ntsire[i])+s] = 1.0;}    
                   if(ntdam[i]==2){DSpred[((ndam[i]-1)*ntsire[i])+s] = mean_vec;}
                 }else{
                   mean_vec /= (ntdam[i]-1);
                   for(d=0; d<ntdam[i]; d++){        
                     S_vec += pow(DSpred[(d*ntsire[i])+s]-mean_vec,2);
                   }
                   S_vec -= pow(DSpred[((ndam[i]-1)*ntsire[i])+s]-mean_vec,2);                        
                   S_vec  /= (ntdam[i]-1);
                   n = double(ntdam[i]-1);
                   N = n+us[usdamcat[i]];
                   S_vec *= N/(n*(N-n));
                   DSpred[((ndam[i]-1)*ntsire[i])+s] = rnorm(mean_vec, sqrt(S_vec));
                   if(DSpred[((ndam[i]-1)*ntsire[i])+s]<0.0){
                     DSpred[((ndam[i]-1)*ntsire[i])+s]=1e-100;
                   }
                 }
                 DSpred_tmp[((ndam[i]-1)*ntsire[i])+s] = log(DSpred[((ndam[i]-1)*ntsire[i])+s]);     
                 tot_mean += DSpred[((ndam[i]-1)*ntsire[i])+s];
               }
             }
             if(father_complete==FALSE){
               for(d=0; d<ndam[i]; d++){
                 mean_vec = 0.0;
                 S_vec = 0.0;   
                 for(s=0; s<ntsire[i]; s++){                   
                   mean_vec += DSpred[(d*ntsire[i])+s];
                 }
                 mean_vec -= DSpred[(d*ntsire[i])+(nsire[i]-1)];
                 if(ntsire[i]<3){
                   if(ntsire[i]==1){DSpred[(d*ntsire[i])+(nsire[i]-1)] = 1.0;}    
                   if(ntsire[i]==2){DSpred[(d*ntsire[i])+(nsire[i]-1)] = mean_vec;}
                 }else{
                   mean_vec /= (ntsire[i]-1);
                   for(s=0; s<ntsire[i]; s++){    
                     S_vec += pow(DSpred[(d*ntsire[i])+s]-mean_vec,2);
                   }
                   S_vec -= pow(DSpred[(d*ntsire[i])+(nsire[i]-1)]-mean_vec,2);          
                   S_vec  /= (ntsire[i]-1);
                   n = double(ntsire[i]-1);
                   N = n+us[ussirecat[i]+nusd];
                   S_vec *= N/(n*(N-n));
                   DSpred[(d*ntsire[i])+(nsire[i]-1)] = rnorm(mean_vec, sqrt(S_vec));
                   if(DSpred[(d*ntsire[i])+(nsire[i]-1)] <0.0){
                     DSpred[(d*ntsire[i])+(nsire[i]-1)]=1e-100;
                   }
                 }
                 DSpred_tmp[(d*ntsire[i])+(nsire[i]-1)] = log(DSpred[(d*ntsire[i])+(nsire[i]-1)]);         
                 tot_mean += DSpred[(d*ntsire[i])+(nsire[i]-1)];
               }
             }  
             if(mother_complete==FALSE && father_complete==FALSE){
               DSpred[((ndam[i]-1)*ntsire[i])+(nsire[i]-1)] = tot_mean/(ndam[i]+nsire[i]);
               DSpred_tmp[((ndam[i]-1)*ntsire[i])+(nsire[i]-1)] = log(tot_mean/(ndam[i]+nsire[i]));
             }
           }else{
             if(sireV && father_complete==FALSE){                  
               mean_vec = meanc(Spred)[0];
               n = double(ntsire[i]-1);
               mean_vec -= Spred[nsire[i]-1]/(n+1.0);
               if(n<2){
                 if(n==0){Spred[nsire[i]-1] = 1.0;}    
                 if(n==1){Spred[nsire[i]-1] = mean_vec*2.0;}
               }else{
                 mean_vec *= (n+1.0)/n;
                 S_vec = varc(Spred)[0];
                 N = n+us[ussirecat[i]+nusd];
                 S_vec *= N/(n*(N-n));
                 Spred[nsire[i]-1] = rnorm(mean_vec, sqrt(S_vec));
                 if(Spred[nsire[i]-1]<0.0){
                   Spred[nsire[i]-1]=1e-100;
                 }
               }
               Spred_tmp[nsire[i]-1] = log(Spred[nsire[i]-1]);   
             }
             if(damV && mother_complete==FALSE){
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
           }

/**************************************/
/* combine matrices and genetic data  */
/**************************************/

           if(damsireV==false && sireV==false && damV==false){   // no variables
             DSpreds = X_design_G[i];
           }else{
             DSpreds = log(X_design_G[i]);
           }


           if(damsireV){
             if(mother_complete && father_complete){
               DSpreds =  DSpreds + DSpred_tmp;
             }else{
               for(d=0; d<ndam[i]; d++){
                 for(s=0; s<nsire[i]; s++){            
                   DSpreds[(d*nsire[i])+s] += DSpred_tmp[(d*ntsire[i])+s]; 
                 }
               }
             }
           }else{
             cnt=0;
             if(sireV && damV){ 
               for(d=0; d<ndam[i]; d++){
                 for(s=0; s<nsire[i]; s++){
                   DSpreds[cnt] += Spred_tmp[s]+Dpred_tmp[d];
                   cnt++;
                 }
               }
             }else{
               if(sireV){ 
                 for(d=0; d<ndam[i]; d++){
                   for(s=0; s<nsire[i]; s++){
                     DSpreds[cnt] += Spred_tmp[s];
                     cnt++;
                   }
                 }
               }
               if(damV){ 
                 for(d=0; d<ndam[i]; d++){
                   for(s=0; s<nsire[i]; s++){
                     DSpreds[cnt] += Dpred_tmp[d];
                     cnt++;
                   }
                 }
               }
             }
           }  

           if(damsireV || sireV ||  damV){ 
             DSpreds = DSpreds - (maxc(DSpreds)[0]-650.0);
             DSpreds = exp(DSpreds);
           }

           if(nuss>0){
             for(d=0; d<ndam[i]; d++){
                DSpreds[(nsire[i]*(d+1))-1] *= us[ussirecat[i]+nusd];
             }
           }

           if(nusd>0){
             for(s=0; s<nsire[i]; s++){
               DSpreds[(nsire[i]*(ndam[i]-1))+s] *= us[usdamcat[i]];
             }
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
                         new_dam_pos = nsire[i]*j;
                         for(b=0; b<nsire[i]; b++){
                           DSpreds[new_dam_pos+b]=0.0;
                         }
                         break;
                       }
                     }
                     for(j=0; j<nsire[i]; j++){
                       if(Sires_vec[i][j]==offid[s]){
                         new_sire_pos = j;
                         for(b=0; b<ndam[i]; b++){
                           DSpreds[new_sire_pos+nsire[i]*b]=0.0;
                         }
                         break;
                       }
                     }
                   }
                }
              }  
            }

            newpar = rmultinom_size1M(DSpreds, nsire[i]*ndam[i]);      // samples a parental combination

            new_dam_pos = newpar/nsire[i];
            new_sire_pos = newpar-(new_dam_pos*nsire[i]);
		
            dam[offid[i]] = Dams_vec[i][new_dam_pos];     // puts new dam.id into pedigree
            sire[offid[i]] = Sires_vec[i][new_sire_pos];  // puts new sire.id into pedigre        
 
        }
     delete [] dsmatch;
}



