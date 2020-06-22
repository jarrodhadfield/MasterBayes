#include "SampDomG.h"

void sampDomG(int nsamp, int **Gobs, int **G, int *nall, int nloci, int *id, double **A, int *categories, double **E_mat, int maxall, int maxrep, int *dam, int *sire, int nind, bool estA){

        int l;                    // locus
        int i;                    // sample
        int o; 
        int ind = 0;              // individual
        int n = 0;                // number of samples within indiviual
        int ns;                   // sample within indiviual
        int obs_G[50];           // observed genotypes of an individual
        int da;                   // genotypes of individual's dam
        int sa;                   // genotypes of individual's sire
        int off_1;
        int sp_1;
        int cat_tmp[50];      // categories to which those genotypes belong
        double Ppart[3];  // Probability vector for true genotypes
        int sample_cat;           // category for observed genotype
        double p;                 // allele frequnecy of allele 1
        double q;                 // allele frequnecy of allele 2
        int end = 1;                // a bit of fuckwittery
        int itt_set;              // nth sampled allele per indivdiual
        double newA[2];      // allele counts from true genotypes
        int rel_status; 
        int samp_A;
        int *no_off = new int[nind];

        int** par = new int*[nind];
        for(i = 0; i < nind; ++i){
           par[i] = new int[50];
        }

        for(o=0; o < nind; o++){
           no_off[o] = 0;
        }

        for(o=0; o < nind; o++){                         
           if(sire[o]<nind){
              l = sire[o];
              par[l][no_off[l]] = o;
              no_off[l]++;
           } 
           if(dam[o]<nind && dam[o]!=sire[o]){
              l = dam[o];
              par[l][no_off[l]] = o;
              no_off[l]++;
           }
        }

        for(l= 0; l < nloci; l++){                  // iterates through loci

          p = A[l][0];
          q = A[l][1];
          newA[0] = 1.0;
          newA[1] = 1.0;

          for(i=0; i <nsamp; i++){                 // iterates through samples

               if(Gobs[i][l]!=-999){               
                 obs_G[n] = Gobs[i][l];       // gets observed genotypes 
                 cat_tmp[n] = categories[i]*6;
                 n++;      
               } 

                if(i == (nsamp-1)){end=-i;}         // if this is the last sample       
           
/* last record for that individual so start sampling */

if(id[i+end]!=ind){                 // the last record for that individual


// get dam genotype //     
                           
                if(dam[ind]<nind){                 
                 da = G[dam[ind]][l];
                }else{
                 da = -999;
                }
// get sire genotype //
                if(sire[ind]<nind){               
                 sa = G[sire[ind]][l];
                }else{
                 sa = -999;
                }
                              
                rel_status=0;    // no parents //
                if(da!=-999){ 
                  rel_status++;  // has a mother //
                }
                if(sa!=-999){ 
                  rel_status++;  // has a father //
                }
// rel status: 0: no par 1: one par 2: two par  

        switch(rel_status){

          case 0:  
          Ppart[0] = p*p;
          Ppart[1] = 2.0*p*q;
          Ppart[2] = q*q;
          break;

          case 1:
          if(da==-999){
            da=sa;
          }
          if(da==0){
            Ppart[0]=p;
            Ppart[1]=q;
            Ppart[2]=0.0;
          }
          if(da==1){
            Ppart[0]=p/2.0;
            Ppart[1]=(p+q)/2.0;
            Ppart[2]=q/2.0;
          }
          if(da==2){
            Ppart[0]=0.0;
            Ppart[1]=p;
            Ppart[2]=q;
          }
          break;

          case 2:
          if(da==1 || sa==1){
            if(da==1 && sa==1){
                Ppart[0]=0.25;
                Ppart[1]=0.5;
                Ppart[2]=0.25;
            }else{
              if(da==0 || sa==0){
                Ppart[0]=0.5;
                Ppart[1]=0.5;
                Ppart[2]=0.0;
              }else{               
                Ppart[0]=0.0;
                Ppart[1]=0.5;
                Ppart[2]=0.5;
              }
            }
          }else{
            if(da==sa){
              Ppart[1]=0.0;
              Ppart[da]=1.0;
              Ppart[2-da]=0.0;
            }else{
              Ppart[0]=0.0;
              Ppart[1]=1.0;
              Ppart[2]=0.0;
            }
          }
          break;
        }

        for(itt_set=0; itt_set<no_off[ind]; itt_set++){              
 
          o = par[ind][itt_set];

          off_1 = G[o][l];       // offspring and spouse genotypes //

          sp_1 = -999;

          if(sire[o]==dam[o]){
            sp_1=-998;
          }else{
            if(sire[o]!=ind){
              if(sire[o]<nind){
                sp_1 = G[sire[o]][l];
              }
            }else{
              if(dam[o]<nind){
                sp_1 = G[dam[o]][l];
              }
            } 
          }

          if(sp_1>-998){
            if(off_1==1){
              if(sp_1==0){
                Ppart[0] *= 0.0;
                Ppart[1] *= 0.5;
                Ppart[2] *= 1.0;
              }
              if(sp_1==2){
                Ppart[0] *= 1.0;
                Ppart[1] *= 0.5;
                Ppart[2] *= 0.0;
              }
            }else{
              if(off_1==0){
                Ppart[0] *= 1.0;
                Ppart[1] *= 0.5;
                Ppart[2] *= 0.0;
              }else{
                Ppart[0] *= 0.0;
                Ppart[1] *= 0.5;
                Ppart[2] *= 1.0;
              }
            }
          }else{
            if(sp_1==-999){
              if(off_1==1){
                Ppart[0] *= q;
                Ppart[1] *= (p+q)/2.0;
                Ppart[2] *= p;
              }else{
                if(off_1==0){
                  Ppart[0] *= p;
                  Ppart[1] *= p/2.0;
                  Ppart[2] *= 0.0;
                }else{
                  Ppart[0] *= 0.0;
                  Ppart[1] *= q/2.0;
                  Ppart[2] *= q;
                }
              }
            }else{
              if(off_1==0){
                Ppart[1] *= 0.25;
                Ppart[2] = 0.0;
              }
              if(off_1==1){
                Ppart[0] = 0.0;
                Ppart[2] = 0.0;
              }
              if(off_1==2){
                Ppart[0] = 0.0;
                Ppart[1] *= 0.25;
              }
            }
          }
        }

        for(ns=0; ns < n; ns++){    // iterates through samples because cat may be different
          sample_cat = cat_tmp[ns];
          if(obs_G[ns]==0){
            Ppart[0] *= E_mat[l][sample_cat];
            Ppart[1] *= E_mat[l][1+sample_cat]; 
            Ppart[2] *= E_mat[l][2+sample_cat]; 
          }else{
            Ppart[0] *= E_mat[l][3+sample_cat];
            Ppart[1] *= E_mat[l][4+sample_cat]; 
            Ppart[2] *= E_mat[l][5+sample_cat]; 
          }
        }

        samp_A = rmultinom_size1(Ppart, 3);
        if(rel_status==0){
          newA[0] += 2.0-double(samp_A);
          newA[1] += double(samp_A);
        }
        if(rel_status==1){
          if(samp_A==1){
            if(da==0){
              newA[1]++;
            }
            if(da==1){
              newA[0]+=0.5;
              newA[1]+=0.5;
            }
            if(da==2){
              newA[0]++;
            }
          }else{
            if(samp_A==0){
              newA[0]++;
            }else{
              newA[1]++;
            }
          }
        }
        G[ind][l] = samp_A;
        ind++;         
        n = 0; 
      }                    // end this is the last record if statement
   }   
                           // end sample loop
     if(estA==TRUE){
       A[l][0] = rbeta(newA[0], newA[1]);
       A[l][1] = 1.0-A[l][0];
     } 

     end = 1;
     ind = 0;
  }                         // end loci loop
         for(i = 0; i < nind; i++){   
           delete [] par [i];
         }
         delete [] no_off;
         delete [] par;
}
 
