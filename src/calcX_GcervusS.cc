#include "calcX_G.h"

void calcX_GcervusS(Matrix<double> X_design_GS [], Matrix<double> X_design_G [], int *offid, int noff , int *ndam, int *nsire, std::map<int, int> Dams [], int *dam){
 
         int i, d, dam_poss;

         for(i = 0; i < noff ; i++){          
           X_design_GS[i] = Matrix<double>(nsire[i], 1);              
           dam_poss = Dams[i][dam[offid[i]]];
           for(d = 0; d < nsire[i]; d++){   
              X_design_GS[i][d] = X_design_G[i][(dam_poss*nsire[i])+d];
           }
         }

}

