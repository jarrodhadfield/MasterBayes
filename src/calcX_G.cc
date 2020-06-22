#include "calcX_G.h"

void calcX_G(Matrix<double> X_design_G [], int *offid, int noff , int *ndam, int *nsire, int nind, Matrix<int> Dams_vec [], Matrix<int> Sires_vec [], int **G, int nloci, double **A, int mtype){

        int i;
        int l;
        int d;
        int s;
        int records_off;
        int off_poss;
        int sire_poss;
        int dam_poss;
        int o1=0;
        int o2=0;
        int d1=0;
        int d2=0;
        int s1=0;
        int s2=0;
        double o1inp;
        double o2inp;
        double sl;
        double fl;
        double TacP;
        double p;
        double q;


    for(i = 0; i < noff ; i++){         // iterate through offspring
      X_design_G[i] = ones<double>(ndam[i]*nsire[i], 1); 
    }

     switch(mtype){
 
     case 1:
     case 3:

       for(l = 0; l<nloci ; l++){             // iterates through loci  

         for(i = 0; i < noff ; i++){         // iterate through offspring
           off_poss = offid[i];   
           o1 = G[off_poss][l*2];
           o2 = G[off_poss][(l*2)+1];

            for(d = 0; d < ndam[i]; d++){    // gets position of offspring genotype in G

              dam_poss = Dams_vec[i][d];       // gets position of dam genotype in G

              if(dam_poss<nind){                 // if the dam is sampled (and not a zero_prob) 
               d1 = G[dam_poss][l*2];
               d2 = G[dam_poss][(l*2)+1]; 
              }

              for(s = 0; s < nsire[i]; s++){      // iterate through sires

               records_off = (nsire[i]*d)+s;
               sire_poss = Sires_vec[i][s];           // gets position of sire genotype in G   

               if(sire_poss<nind){     
                 s1 = G[sire_poss][l*2];
                 s2 = G[sire_poss][(l*2)+1]; 
               }
         
// both sexes sampled 
                 
               if(sire_poss<nind && dam_poss<nind){        
   
                 TacP = 0.0;
       
                 if((d1==o1 && s1==o2) || (d1==o2 && s1==o1)){TacP +=0.25;}
                 if((d1==o1 && s2==o2) || (d1==o2 && s2==o1)){TacP +=0.25;}
                 if((d2==o1 && s1==o2) || (d2==o2 && s1==o1)){TacP +=0.25;}
                 if((d2==o1 && s2==o2) || (d2==o2 && s2==o1)){TacP +=0.25;}

                 X_design_G[i][records_off] *=  TacP;
               }                              

// female sampled, male unsampled 

               if(sire_poss>=nind && dam_poss<nind){                                  
                 sl = 0.0;
                 fl = 0.0;
                 if(d1==o1 || d1==o2){sl += 0.5;}
                 if(d2==o1 || d2==o2){sl += 0.5;}
                 if(sl>0.25){
                   if(o1!=o2){
                     if(o1==d1 || o1==d2){o1inp = 1.0;}else{o1inp = 0.0;}
                     if(o2==d1 || o2==d2){o2inp = 1.0;}else{o2inp = 0.0;}
                     if(o1inp>0.5){
                       fl += A[l][o2];
                     }
                     if(o2inp>0.5){
                       fl += A[l][o1];
                     }
                     fl /= o1inp+o2inp;
                   }else{
                     fl = A[l][o1];
                   }
                 }  
                 X_design_G[i][records_off] *=  sl*fl;
               } 

// male sampled, female unsampled 

               if(sire_poss<nind && dam_poss>=nind){                                   
                 sl = 0.0;
                 fl = 0.0;
                 if(s1==o1 || s1==o2){sl += 0.5;}
                 if(s2==o1 || s2==o2){sl += 0.5;}
                 if(sl>0.25){
                   if(o1!=o2){
                     if(o1==s1 || o1==s2){o1inp = 1.0;}else{o1inp = 0.0;}
                     if(o2==s1 || o2==s2){o2inp = 1.0;}else{o2inp = 0.0;}
                     if(o1inp>0.5){
                       fl += A[l][o2];
                     }
                     if(o2inp>0.5){
                       fl += A[l][o1];
                     }
                     fl /= o1inp+o2inp;
                   }else{
                     fl = A[l][o1];
                   }
                 }
                 X_design_G[i][records_off] *=  sl*fl;
               }
              // both parents unsampled 

               if(sire_poss>=nind && dam_poss>=nind){
                 X_design_G[i][records_off] *=  A[l][o1]*A[l][o2];
                 if(o1!=o2){
                   X_design_G[i][records_off] *=  2.0;
                 }
               }
             }
           }  
         }
       }
     break;

     case 2:
        for(l = 0; l<nloci ; l++){             // iterates through loci  

          p = A[l][0];                        // sl frequency of NULL allle
          q = A[l][1];

          for(i = 0; i < noff ; i++){         // iterate through offspring
            off_poss = offid[i];   
            o1 = G[off_poss][l];

            for(d = 0; d < ndam[i]; d++){    // gets position of offspring genotype in G

              dam_poss = Dams_vec[i][d];       // gets position of dam genotype in G

              if(dam_poss<nind){                 // if the dam is sampled (and not a zero_prob) 
               d1 = G[dam_poss][l];
              }

              for(s = 0; s < nsire[i]; s++){      // iterate through sires
                records_off = (nsire[i]*d)+s;
                sire_poss = Sires_vec[i][s];           // gets position of sire genotype in G   

                if(sire_poss<nind){     
                  s1 = G[sire_poss][l];
                }
         
// both sexes sampled 
                 
                if(sire_poss<nind && dam_poss<nind){ 
       
                  TacP = 0.0;

                  if(d1+s1==2){
                    if(d1==1 && s1==1){
                      TacP = 0.25;
                      if(o1==1){
                        TacP += 0.25;
                      }
                    }else{
                      if(o1==1){
                        TacP = 1.0;
                      }
                    }
                  }else{
                    if(o1==d1){
                      TacP +=0.5;
                    }
                    if(o1==s1){
                      TacP +=0.5;
                    }
                  }
                  X_design_G[i][records_off] *=  TacP;
                }                              

// female sampled, male unsampled 

                if(sire_poss>=nind && dam_poss<nind){    

                  TacP = 0.0;
                  if(o1==1){
                    if(d1<2){
                      TacP = q;
                    }
                    if(d1>0){
                      TacP += p;
                    }
                  }else{
                    if(o1==0 && d1<2){
                      TacP=p;
                    }
                    if(o1==2 && d1>0){
                      TacP=q;
                    }
                  }
                  if(d1==1){
                    TacP /= 2.0;
                  }
                  X_design_G[i][records_off] *=  TacP;
                } 

// male sampled, female unsampled 

                if(sire_poss<nind && dam_poss>=nind){                                   
                  TacP = 0.0;
                  if(o1==1){
                    if(s1<2){
                      TacP = q;
                    }
                    if(s1>0){
                      TacP += p;
                    }
                  }else{
                    if(o1==0 && s1<2){
                      TacP=p;
                    }
                    if(o1==2 && s1>0){
                      TacP=q;
                    }
                  }
                  if(s1==1){
                    TacP /= 2.0;
                  }
                  X_design_G[i][records_off] *=  TacP;
                } 

                if(sire_poss>=nind && dam_poss>=nind){
                  if(o1==0){
                    X_design_G[i][records_off] *=  p*p;
                  }
                  if(o1==1){
                    X_design_G[i][records_off] *=  2*p*q;
                  }
                  if(o1==2){
                    X_design_G[i][records_off] *=  q*q;
                  }
                }
              }
            }  
          }
        }
        break;
        }
}

