#include "calcX_G.h"

void calcX_Gcervus(Matrix<double> X_design_G [], int *offid, int noff , int *ndam, int *nsire, int nind, Matrix<int> Dams_vec [], Matrix<int> Sires_vec [], int **G, int nloci, double **A, double **E_mat, int mtype, int *nall){

        int i;
        int l;
        int d;
        int s;
        int na;
        int dj;
        int dk;
        int sj;
        int sk;
        int oj;
        int ok;
        int hit;
        int records_off;
        int off_poss;
        int sire_poss;
        int dam_poss;
        int o1 = -999;
        int o2 = -999;
        int d1 = -999;
        int d2 = -999;
        int s1 = -999;
        int s2 = -999;
        double o1inp;
        double o2inp;
        double sl;
        double fl;
        double p;
        double q;
        double T[6];
        int Aobs[7];
        double P[7];
        double TacP;
        double pG;
        double sslfl=0.0;
        double dslfl=0.0;

       for(i = 0; i < noff ; i++){         // iterate through offspring
         X_design_G[i] = ones<double>(ndam[i]*nsire[i], 1); 
       }

       switch(mtype){

       case 1:

       // Microsatellite with CERVUS errors

         for(l = 0; l<nloci ; l++){             // iterates through loci  

           for(i = 0; i < noff ; i++){         // iterate through offspring

             off_poss = offid[i];              // gets position of offspring genotype in G
             o1 = G[off_poss][l*2];
             o2 = G[off_poss][(l*2)+1];

             for(d = 0; d < ndam[i]; d++){    
 
               dam_poss = Dams_vec[i][d];       // gets position of dam genotype in G

               if(dam_poss<nind){
                 d1 = G[dam_poss][l*2];
                 d2 = G[dam_poss][(l*2)+1]; 
               }else{
                 d1 = -999;
               }
               for(s = 0; s < nsire[i]; s++){      // iterate through sires

                 records_off = (nsire[i]*d)+s;
                 sire_poss = Sires_vec[i][s];           // gets position of sire genotype in G   

                 if(sire_poss<nind){ 
                   s1 = G[sire_poss][l*2];
                   s2 = G[sire_poss][(l*2)+1]; 
                 }else{
                   s1 = -999;
                 }

                 if(o1!=-999){
                   pG = A[l][o1]*A[l][o2];
                   if(o1!=o2){
                     pG *= 2.0;
                   }
                   if(d1!=-999){ // if dam has genotype
                     sl = 0.0;
                     fl = 0.0;
                     if(d1==o1 || d1==o2){sl += 0.5;}
                     if(d2==o1 || d2==o2){sl += 0.5;}
                     if(sl>0.25){
                       if(o1!=o2){
                         if(o1==d1 || o1==d2){o1inp = 1.0;}else{o1inp = 0.0;}
                         if(o2==d1 || o2==d2){o2inp = 1.0;}else{o2inp = 0.0;}
                         if(o1inp>0.5){fl += A[l][o2];}  
                         if(o2inp>0.5){fl += A[l][o1];} 
                         fl /= o1inp+o2inp;
                       }else{
                         fl = A[l][o1];
                       }
                     }  
                     dslfl = sl*fl;
                   }
            
                   if(s1!=-999){ // if sire has genotype
                     sl = 0.0;
                     fl = 0.0;
                     if(s1==o1 || s1==o2){sl += 0.5;}
                     if(s2==o1 || s2==o2){sl += 0.5;}
                     if(sl>0.25){
                       if(o1!=o2){
                         if(o1==s1 || o1==s2){o1inp = 1.0;}else{o1inp = 0.0;}
                         if(o2==s1 || o2==s2){o2inp = 1.0;}else{o2inp = 0.0;}
                         if(o1inp>0.5){fl += A[l][o2];}
                         if(o2inp>0.5){fl += A[l][o1];}
                         fl /= o1inp+o2inp;
                       }else{
                         fl = A[l][o1];
                       }
                     }  
                     sslfl = sl*fl;
                   }

                   if(s1==-999 && d1==-999){       // if nobody sampled
                     X_design_G[i][records_off] *= pG; 
                   }
                   if(d1!=-999 && s1==-999){      // if dam sampled only
                     X_design_G[i][records_off] *=  pow(E_mat[l][0],2.0)*dslfl+pG*E_mat[l][1]*(2.0-E_mat[l][1]);
                   }
                   if(s1!=-999 &&  d1==-999){      // if sire sampled only
                     X_design_G[i][records_off] *=  pow(E_mat[l][0],2.0)*sslfl+pG*E_mat[l][1]*(2.0-E_mat[l][1]);
                   }
                   if(s1!=-999 && d1!=-999){      // if both sampled only

                     TacP = 0.0;
          
                     if((d1==o1 && s1==o2) || (d1==o2 && s1==o1)){TacP +=0.25;}
                     if((d1==o1 && s2==o2) || (d1==o2 && s2==o1)){TacP +=0.25;}
                     if((d2==o1 && s1==o2) || (d2==o2 && s1==o1)){TacP +=0.25;}
                     if((d2==o1 && s2==o2) || (d2==o2 && s2==o1)){TacP +=0.25;}

                     if(dam_poss!=sire_poss){  // non-selfers
                       X_design_G[i][records_off] *=  pow(E_mat[l][0],3.0)*TacP+E_mat[l][1]*pow(E_mat[l][0],2.0)*(sslfl+dslfl)+pG*(E_mat[l][1]+pow(E_mat[l][1],2.0)-pow(E_mat[l][1],3.0));
                     }else{                   // selfers
                       X_design_G[i][records_off] *=  pow(E_mat[l][0],2.0)*TacP+E_mat[l][1]*pG+0.5*E_mat[l][1]*(1.0-E_mat[l][1])*(pG+double(o1==o2)*pow(pG,0.5));
                     }
                   } 
                 }
               }
             }
           }
         }
       break;

  // AFLP

       case 2:     
         for(l = 0; l<nloci ; l++){             // iterates through loci  

             p=A[l][0]; // frequency of 0 allele
             q=A[l][1]; // frequency of 1 allele

             T[0] = p/(1.0+p);                           // Pr(0|1,0) 
             T[1] = 1.0/(1.0+p);                         // Pr(1|1,0)  
             T[2] = pow(T[0],2.0);                       // Pr(0|1,1) 
             T[3] = (1.0+(2.0*p))/pow(1.0+p,2.0);        // Pr(1|1,1)  
                                                         // Pr(0|0,0) = 1
             T[4] = 0.5*T[0];                            // Pr(0|1)  when selfing
             T[5] = 1.5*T[0]+(1-p)/(1+p);                // Pr(1|1)  when selfing P(?|1,0) does not exist and P(0|0)=1 and P(1|0)=0

             P[0] = pow(p,2.0)*E_mat[l][0];                                    // Pr(act = 0 | obs = 0 )
             P[1] = pow(p,2.0)*E_mat[l][3];                                    // Pr(act = 0 | obs = 1 )
             P[2] = (pow(1.0-p,2.0)*E_mat[l][2])+(2.0*p*(1.0-p)*E_mat[l][1]);  // Pr(act = 1 | obs = 0 )
             P[3] = (pow(1.0-p,2.0)*E_mat[l][5])+(2.0*p*(1.0-p)*E_mat[l][4]);  // Pr(act = 1 | obs = 1 )

             o1inp = P[0]+P[2];
             P[0] /= o1inp;                              
             P[2] /= o1inp;
             o1inp = P[1]+P[3];
             P[1] /= o1inp;              
             P[3] /= o1inp;


           for(i = 0; i < noff ; i++){         // iterate through offspring
             off_poss = offid[i];   
             o1 = G[off_poss][l];

             for(d = 0; d < ndam[i]; d++){    // gets position of offspring genotype in G
 
               dam_poss = Dams_vec[i][d];       // gets position of dam genotype in G
 
               if(dam_poss<nind){
                 d1 = G[dam_poss][l];
               }else{
                 d1 = -999;
               }

               for(s = 0; s < nsire[i]; s++){      // iterate through sires

                 records_off = (nsire[i]*d)+s;
                 sire_poss = Sires_vec[i][s];           // gets position of sire genotype in G   

                 if(sire_poss<nind){ 
                   s1 = G[sire_poss][l];
                 }else{
                   s1 = -999;
                 }

// all sampled and all genotypes known

                 if(o1!=-999){
                   if(d1!=-999 && s1!=-999){
                    if(dam_poss==sire_poss){  // selfing - so don't double count Pr(act|obs) for parent
                      TacP = P[o1]*P[d1];
                      TacP += ((T[4]*P[o1])+(T[5]*P[o1+2]))*P[d1+2];
                    }else{
                      TacP = P[o1]*P[d1]*P[s1];
                      TacP += ((T[0]*P[o1])+(T[1]*P[o1+2]))*((P[d1+2]*P[s1])+(P[d1]*P[s1+2]));
                      TacP += ((T[2]*P[o1])+(T[3]*P[o1+2]))*P[d1+2]*P[s1+2];
                    }
                   }else{
                     if(d1!=-999 || s1!=-999){
                       if(d1==-999){
                         d1=s1;
                       }
                       TacP = P[o1]*P[d1]*p;
                       TacP += P[o1]*P[d1+2]*p*(p/(1.0+p));
                       TacP += P[o1+2]*P[d1]*(1-p);
                       TacP += P[o1+2]*P[d1+2]*((1+p*(1-p))/(1+p));
                     }else{
                       TacP = P[o1]+P[o1+2];
                     }
                   }
                   X_design_G[i][records_off] *= TacP;
                 }
               }
             }
           }
         } 
        break;
 // Microsatellite with Wang errors

        case 3:

        for(l = 0; l<nloci ; l++){             // iterates through loci  

           for(i = 0; i < noff ; i++){         // iterate through offspring

             off_poss = offid[i];   
             o1 = G[off_poss][l*2];
             o2 = G[off_poss][(l*2)+1];

             for(d = 0; d < ndam[i]; d++){    // gets position of offspring genotype in G
 
               dam_poss = Dams_vec[i][d];       // gets position of dam genotype in G

               if(dam_poss<nind){
                 d1 = G[dam_poss][l*2];
                 d2 = G[dam_poss][(l*2)+1]; 
               }else{
                 d1 = -999;
               }

               for(s = 0; s < nsire[i]; s++){      // iterate through sires

                 records_off = (nsire[i]*d)+s;
                 sire_poss = Sires_vec[i][s];           // gets position of sire genotype in G   

                 if(sire_poss<nind){ 
                   s1 = G[sire_poss][l*2];
                   s2 = G[sire_poss][(l*2)+1]; 
                 }else{
                   s1 = -999;
                 }

                 if(o1!=-999){  // probability equal if offspring genotype unobserved.

// get the set of observed alleles and their frequencies

                  Aobs[0] = o1;
                  P[0] = A[l][o1];
                  A[l][o1] = 0.0;                  
                  na = 1;

                  if(o1!=o2){
                   Aobs[1] = o2;
                   P[1] = A[l][o2];
                   A[l][o2] = 0.0;  
                   na = 2;
                  }

                  if(d1!=-999){
                    hit=0;
                    for(dk = 0; dk < na; dk++){   
                      if(d1==Aobs[dk]){
                        hit=1;
                      }
                    }

                    if(hit==0){
                      Aobs[na]=d1;
                      P[na] = A[l][d1];
                      A[l][d1] = 0.0;   
                      na++;
                    }

                    hit=0;
                    for(dk = 0; dk < na; dk++){   
                      if(d2==Aobs[dk]){
                        hit=1;
                      }
                    }

                    if(hit==0){
                      Aobs[na]=d2;
                      P[na] = A[l][d2];
                      A[l][d2] = 0.0;   
                      na++;
                    }
                  }

                  if(s1!=-999){
                    hit=0;
                    for(sk = 0; sk < na; sk++){   
                      if(s1==Aobs[sk]){
                        hit=1;
                      }
                    }

                    if(hit==0){
                      Aobs[na]=s1;
                      P[na] = A[l][s1];
                      A[l][s1] = 0.0;   
                      na++;
                    }

                    hit=0;
                    for(sk = 0; sk < na; sk++){   
                      if(s2==Aobs[sk]){
                        hit=1;
                      }
                    }

                    if(hit==0){
                      Aobs[na]=s2;
                      P[na] = A[l][s2];
                      A[l][s2] = 0.0;   
                      na++;
                    }
                  }

// collapse unobserved alleles if they exist

                  if(na!=nall[l]){
                    P[na] = 0.0;
                    for(dj = 0; dj<nall[l]; dj++){
                      P[na] += A[l][dj];
                      if(A[l][dj]>(1e-10)){
                        Aobs[na]=dj;
                      }
                    }
                    na++;
                  }

// iterate over potential genotypes

                  o1inp = 0.0; 
                  o2inp = 0.0;

                  for(dj = 0; dj<na; dj++){   
                    for(dk = dj; dk<na; dk++){  
                      for(sj = 0; sj<na; sj++){   
                        for(sk = sj; sk<na; sk++){
                          for(oj = 0; oj<na; oj++){   
                            for(ok = oj; ok<na; ok++){

                              pG = 1.0;

// P(O_o | O)
                              if(o1==o2){
                                if(o1==Aobs[oj] && o1==Aobs[ok]){
                                  pG *= E_mat[l][0];
                                }else{
                                  if(o1!=Aobs[oj] && o1!=Aobs[ok]){
                                    pG *= E_mat[l][2];
                                  }else{
                                    pG *= E_mat[l][1];
                                  }
                                }    
                              }else{ 
                                if((o1==Aobs[oj] && o2==Aobs[ok]) || (o2==Aobs[oj] && o1==Aobs[ok])){                                      
                                  pG *= E_mat[l][3];
                                }else{
                                  if((o1!=Aobs[oj] && o2!=Aobs[ok]) || (o2!=Aobs[oj] && o1!=Aobs[ok])){   
                                    pG *= E_mat[l][6];  
                                  }else{
                                    if(oj==ok){
                                      pG *= E_mat[l][4];
                                    }else{
                                      pG *= E_mat[l][5];
                                    }  
                                  }  
                                } 
                              }
// P(M_o | M)
                              if(d1!=-999){
                                if(d1==d2){
                                  if(d1==Aobs[dj] && d1==Aobs[dk]){
                                    pG *= E_mat[l][0];
                                  }else{
                                    if(d1!=Aobs[dj] && d1!=Aobs[dk]){
                                      pG *= E_mat[l][2];
                                    }else{
                                      pG *= E_mat[l][1];
                                    }
                                  }    
                                }else{ 
                                  if((d1==Aobs[dj] && d2==Aobs[dk]) || (d2==Aobs[dj] && d1==Aobs[dk])){                                      
                                    pG *= E_mat[l][3];
                                  }else{
                                    if((d1!=Aobs[dj] && d2!=Aobs[dk]) || (d2!=Aobs[dj] && d1!=Aobs[dk])){   
                                      pG *= E_mat[l][6];  
                                    }else{
                                      if(dj==dk){
                                        pG *= E_mat[l][4];
                                      }else{
                                        pG *= E_mat[l][5];
                                      }  
                                    }  
                                  } 
                                }
                              }
// P(S_o | S)

                              if(s1!=-999 && dam_poss!=sire_poss){
                                if(s1==s2){
                                  if(s1==Aobs[sj] && s1==Aobs[sk]){
                                    pG *= E_mat[l][0];
                                  }else{
                                    if(s1!=Aobs[sj] && s1!=Aobs[sk]){
                                      pG *= E_mat[l][2];
                                    }else{
                                      pG *= E_mat[l][1];
                                    }
                                  }    
                                }else{ 
                                  if((s1==Aobs[sj] && s2==Aobs[sk]) || (s2==Aobs[sj] && s1==Aobs[sk])){                                      
                                    pG *= E_mat[l][3];
                                  }else{
                                    if((s1!=Aobs[sj] && s2!=Aobs[sk]) || (s2!=Aobs[sj] && s1!=Aobs[sk])){   
                                      pG *= E_mat[l][6];  
                                    }else{
                                      if(sj==sk){
                                        pG *= E_mat[l][4];
                                      }else{
                                        pG *= E_mat[l][5];
                                      }  
                                    }  
                                  } 
                                }
                              }  

                              pG *= P[dj]*P[dk];  
                              if(dj!=dk){pG *= 2.0;}

                              TacP = 0.0;

                              if(dam_poss!=sire_poss){  // If not selfing 

                                pG *= P[sj]*P[sk]; 
                                if(sj!=sk){pG *= 2.0;}
       
                                if((dj==oj && sj==ok) || (dj==ok && sj==oj)){TacP += 0.25;}
                                if((dj==oj && sk==ok) || (dj==ok && sk==oj)){TacP += 0.25;}
                                if((dk==oj && sj==ok) || (dk==ok && sj==oj)){TacP += 0.25;}
                                if((dk==oj && sk==ok) || (dk==ok && sk==oj)){TacP += 0.25;}

                              }else{
                                
				if(dj==sj && dk==sk){
                                  if((dj==oj && sj==ok) || (dj==ok && sj==oj)){TacP += 0.25;}
                                  if((dj==oj && sk==ok) || (dj==ok && sk==oj)){TacP += 0.25;}
                                  if((dk==oj && sj==ok) || (dk==ok && sj==oj)){TacP += 0.25;}
                                  if((dk==oj && sk==ok) || (dk==ok && sk==oj)){TacP += 0.25;}
                                }else{
                                  pG = 0.0;
                                }
                              }
                              o1inp +=  pG;
                              o2inp +=  pG*TacP;
                            }
                          }
                        }
                      }
                    }
                  }
                  X_design_G[i][records_off] *=  o2inp/o1inp;
// read back in used allele frequencies
                  for(dj = 0; dj<(na-1); dj++){ 
                    A[l][Aobs[dj]] = P[dj];
                  }

                  if(A[l][Aobs[(na-1)]]<(1e-10)){
                    A[l][Aobs[(na-1)]] = P[(na-1)];
                  }
                }
              }
            }
          }
        }
        break;


       }
}

