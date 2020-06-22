#include "calcX_G.h"

void calcX_GcervusD(Matrix<double> X_design_GD [], Matrix<double> X_design_G [], int *offid, int noff , int *ndam, int *nsire, std::map<int, int> Sires [], int *sire){
 
         int i, d, sire_poss;

         for(i = 0; i < noff ; i++){          
           X_design_GD[i] = Matrix<double>(ndam[i], 1);              
           sire_poss = Sires[i][sire[offid[i]]];
           for(d = 0; d < ndam[i]; d++){   
              X_design_GD[i][d] = X_design_G[i][(d*nsire[i])+sire_poss];
           }
         }
}

