#include "calcX_G.h"

void calcX_Gmm(Matrix<int> mmD [], Matrix<int> mmS [], int *offid, int noff , int *ndam, int *nsire, int nind, Matrix<int> Dams_vec [], Matrix<int> Sires_vec[],int **G, int nloci){

        int i;
        int l;
        int d;
        int s;
        int off_poss;
        int sire_poss;
        int dam_poss;
        int o1;
        int o2;
        int d1;
        int d2;
        int s1;
        int s2;

       for(i = 0; i < noff ; i++){         // iterate through offspring
         mmD[i] = ones<int>(ndam[i], 1)-1; 
         mmS[i] = ones<int>(nsire[i], 1)-1;
       }

          for(l = 0; l<nloci ; l++){             // iterates through loci  

           for(i = 0; i < noff ; i++){         // iterate through offspring
             off_poss = offid[i];   
             o1 = G[off_poss][l*2];
             o2 = G[off_poss][(l*2)+1];
 
             if(o1!=-999){
               for(d = 0; d < ndam[i]; d++){    // gets position of offspring genotype in G
 
                 dam_poss = Dams_vec[i][d];       // gets position of dam genotype in G
  
                 if(dam_poss<nind){                 // if the dam is sampled (and not a zero_prob) 
                   d1 = G[dam_poss][l*2];
                   d2 = G[dam_poss][(l*2)+1]; 
                   if(d1!=-999){
                     if(o1!=d1 && o2!=d2 && o1!=d2 && o2!=d1){
                       mmD[i][d]++;
                     }
                   } 
                 }
               }

               for(s = 0; s < nsire[i]; s++){      // iterate through sires

                 sire_poss = Sires_vec[i][s];           // gets position of sire genotype in G   
  
                 if(sire_poss<nind){     
                   s1 = G[sire_poss][l*2];
                   s2 = G[sire_poss][(l*2)+1]; 
                   if(s1!=-999){
                     if(o1!=s1 && o2!=s2 && o1!=s2 && o2!=s1){
                       mmS[i][s]++;
                     }
                   }
                 }
               }
             } 
           }
         }
}

