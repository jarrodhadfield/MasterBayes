#include "ErrorMat.h"

void Error_Mat(Matrix<double> E1, Matrix<double> E2, double **E_mat, int ncat, int *nall, int nloci, bool LogL, bool perlocus, int mtype){

/* if codom = FALSE then E1 is the probability 
of a 1 allele being scored as 0, and E2 is the 
probability of a 0 allele being scored as 1 */

   int l;
   int cat;
   double e1;
   double e2;
   double ve2;
   int loop;

   switch(mtype){

     case 1:
     for(l= 0; l < nloci; l++){                  // iterates through loci
       for(cat= 0; cat < ncat; cat++){          // iterates through loci
         if(perlocus==FALSE){
           e2 = E2[cat]*(2.0-E2[cat]);
           loop = cat*7;
         }else{
           e2 = E2[(l*ncat)+cat]*(2.0-E2[(l*ncat)+cat]);
           loop = cat*7;
         }
         if(LogL==true){
           E_mat[l][loop] = log(1-e2);
           E_mat[l][1+loop] = log(e2);
           E_mat[l][2+loop] = log(E2[cat]);
         }else{
           E_mat[l][loop] = 1-e2;
           E_mat[l][1+loop] = e2;
           E_mat[l][2+loop] = E2[cat];
         }
       }
     }
     break;

     case 2:
     for(l= 0; l < nloci; l++){                  // iterates through loci
       for(cat= 0; cat < ncat; cat++){          // iterates through loci
         if(perlocus==FALSE){
           e1 =  E1[cat];
           e2 =  E2[cat];
           loop = cat*6;
         }else{
           e1 = E1[(l*ncat)+cat];
           e2 =  E2[(l*ncat)+cat];
           loop = cat*6;
         }
         if(LogL==true){
           E_mat[l][loop] = log(pow(1.0-e2,2.0));
           E_mat[l][1+loop] = log(e1*(1.0-e2));
           E_mat[l][2+loop] = log(pow(e1,2.0));
           E_mat[l][3+loop] = log(e2*(2.0-e2));
           E_mat[l][4+loop] = log(1.0+e1*(e2-1.0));
           E_mat[l][5+loop] = log(1.0-pow(e1,2.0));
         }else{
           E_mat[l][loop] = pow(1.0-e2,2.0);
           E_mat[l][1+loop] = e1*(1.0-e2);
           E_mat[l][2+loop] = pow(e1,2.0);
           E_mat[l][3+loop] = e2*(2.0-e2);
           E_mat[l][4+loop] = 1.0+e1*(e2-1.0);
           E_mat[l][5+loop] = 1.0-pow(e1,2.0);
         }
       }
     }
     break;

     case 3:

     for(l= 0; l < nloci; l++){                  // iterates through loci
       for(cat= 0; cat < ncat; cat++){          // iterates through loci
         if(perlocus==FALSE){
           ve2 = E2[cat];
           e1 =  E1[cat]/(1.0+E1[cat]);
           e2 = ve2/(double(nall[l])-1.0);
           loop = cat*7;
         }else{
           ve2 = E2[(l*ncat)+cat];
           e1 =  E1[(l*ncat)+cat]/(1.0+E1[(l*ncat)+cat]);
           e2 = ve2/(double(nall[l])-1.0);
           loop = cat*7;
         }
         if(LogL==true){
           E_mat[l][loop] = log(pow(1.0-ve2,2.0));
           E_mat[l][1+loop] = log(e2*(1.0-ve2)+(e1*pow(1.0-ve2-e2,2.0)));
           E_mat[l][2+loop] = log(pow(e2,2.0));
           E_mat[l][3+loop] = log(pow(1.0-ve2,2.0)+pow(e2,2.0)-(2.0*e1*pow(1.0-ve2-e2,2.0)));
           E_mat[l][4+loop] = log(2.0*e2*(1.0-ve2));
           E_mat[l][5+loop] = log(e2*(1.0-ve2+e2));
           E_mat[l][6+loop] = log(2.0*pow(e2,2.0));
         }else{
           E_mat[l][loop] = pow(1.0-ve2,2.0);
           E_mat[l][1+loop] = (e2*(1.0-ve2))+(e1*pow(1.0-ve2-e2,2.0));
           E_mat[l][2+loop] = pow(e2,2.0);
           E_mat[l][3+loop] = pow(1.0-ve2,2.0)+pow(e2,2.0)-(2.0*e1*pow(1.0-ve2-e2,2.0));
           E_mat[l][4+loop] = 2.0*e2*(1.0-ve2);
           E_mat[l][5+loop] = e2*(1.0-ve2+e2);
           E_mat[l][6+loop] = 2.0*pow(e2,2.0);
         }
       }
     }
     break;

   }
 }

 
