#include "SampGC.h"

void sampGC(int nsamp, int **Gobs, int **G, int *nall, int nloci, int *id, double **A, int *categories, double **E_mat, int maxall, int maxrep, int *dam, int *sire, int nind, bool estA){
         
        int l;                    // locus
        int i;                    // sample
        int o; 
        int ind = 0;              // individual
        int n = 0;                // number of samples within indiviual
        int ns;                   // sample within indiviual
        int da[2];                // genotypes of individual's dam
        int sa[2];                // genotypes of individual's sire
        int S_1;
        int S_2;
        int off_1;
        int off_2;
        int sp_1;
        int sp_2;
        double tac;
        int which_vec = 0;
        int offhom_sizeS;
        int *cat_tmp = new int[maxrep];      // categories to which those genotypes belong
        double *Ppart = new double[maxall*maxall];        // Probability vector for true genotypes
        int a1;                   // observed allele 1
        int a2;                   // observed allele 2
        int sample_cat;           // category for observed genotype
        int na1;                  // temp allele 1
        int na2;                  // temp allele 2
        double u;                 // allele frequnecy of allele 1
        double v;                 // allele frequnecy of allele 2
        double z;                 // allele frequnecy of allele 1+2
        int nl;                   // number of alleles at locus l
        int cnt;                  // temporary counter
        int q = 1;                // a bit of fuckwittery
        int misstype = 0;         // are there dupliacate genotypes that differ
        int *all_set = new int[maxall];      // temporary vector for counting sampled alleles
        int *unique_set = new int[maxall];   // set of sampled alleles
        int set_size;             // number of sampled alleles per individual 
        int itt_set;              // nth sampled allele per indivdiual
        double *PA = new double[maxall];        // temporary vector of allele frequencies
        double *newA = new double[maxall];      // allele counts from true genotypes
        int gen_status; 
        int rel_status; 
        int samp_A;
        int no_base;
        int *no_off = new int[nind];
        int no_all_in_S;
        int IBL [2];             // variable for recording the number of alleles in offspring produced by selfing 
                                 // if there is only a single allele in the selfed offspring IBL[0] is the allele
                                 // and IBL[1] is the number of IBL[0] alleles in the selfed offspring.  If there 
                                 // are 2 alleles then IBL[1] is the second allele and the genotype of the indiviual is known.

        int** obs_G = new int*[maxrep];
        for(i = 0; i < maxrep; ++i){
           obs_G[i] = new int[2];
        }

        double** vec = new double*[2];
        for(i = 0; i < 2; ++i){
           vec[i] = new double[maxall];
        }

        int** par = new int*[nind];
        for(i = 0; i < nind; ++i){
           par[i] = new int[100];
        }


/*
        int oldG1;         // some variables that are used for bug checking
        int oldG2;
        bool problem;
        int j;
*/
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
        
        IBL[0] = -999;
        IBL[1] = -999;

        for(l= 0; l < nloci; l++){                  // iterates through loci

          nl = nall[l];                               // number of alleles at that locus

          for(i= 0; i < nl; i++){                     // create temporary vector of locus allele frequencies
            PA[i] = A[l][i];
            newA[i] = 1.0;
          }      

          for(i=0; i <nsamp; i++){                 // iterates through samples
               if(Gobs[i][l*2]!=-999){               
                 obs_G[n][0] = Gobs[i][l*2];       // gets observed genotypes 
                 obs_G[n][1] = Gobs[i][(l*2)+1];    
                 cat_tmp[n] = categories[i]*7;
             
                 if(n>0 && (obs_G[n][0]!=obs_G[0][0] || obs_G[n][1]!=obs_G[0][1])){
                   misstype=1;  // if sampled once, or sampled > 1 but identically then no APPARENT misstypes
                 }
                 n++;      
               } 

                if(i == (nsamp-1)){q=-i;}         // if this is the last sample       
           
/* last record for that individual so start sampling */

if(id[i+q]!=ind){                 // the last record for that individual

//               oldG1 = G[ind][l*2];          // some variables that are used for bug checking
//               oldG2 = G[ind][(l*2)+1];

// get dam genotype //                
                if(dam[ind]<nind){                 
                 da[0]=G[dam[ind]][l*2];
                 da[1]=G[dam[ind]][(l*2)+1];
                }else{
                 da[0]=-999;
                }
// get sire genotype //
                if(sire[ind]<nind){               
                 sa[0]=G[sire[ind]][l*2];
                 sa[1]=G[sire[ind]][(l*2)+1];
                }else{
                 sa[0]=-999;
                }
                              
// get offspring/spouse genotypes //

                if(n==0){
                 gen_status=0;  // no observed genotype //
                }else{
                  if(misstype==0){
                    gen_status=1;   // observed but consistent genotype (HET) //
                    if(obs_G[0][0]==obs_G[0][1]){
                     gen_status=2;  // observed but consistent genotype (HOM) //
                    }
                  }else{
                   gen_status=3;    // observed and inconsistent genotype //
                  }
                }
                
// gen status: 0: no Gobs 1: Gobs HET consistent 2: Gobs HOM consistent 3: Gobs inconsistent

                if(no_off[ind]==0){
                  rel_status=0;    // no relatives //
                  if(da[0]!=-999){ 
                    rel_status++;  // has a mother and no offspring //
                  }
                  if(sa[0]!=-999){ 
                    rel_status++;  // has a father and no offspring (1 if mother not recorded 2 if it is)/
                  }
                }else{
                  rel_status = 4;  // has offspring and one or no parents //
                  if(da[0]!=-999 && sa[0]!=-999){ 
                    rel_status= 3;  // has offspring and both parents mother //
                  }
                }
                if(rel_status==1 && gen_status!=0){rel_status=4;}

// rel status: 0: no rel 1: one par but no observed genotype!! 2: both par  3: offspring + both par 4: everything else 

switch(rel_status){
    
 case 0: // no contextual information provided by relatives //

      switch(gen_status){

      case 0: // no observed genotype //

                 samp_A = rmultinom_size1(PA, nl);
                 newA[samp_A] ++;
                 G[ind][l*2] = samp_A;
                 samp_A = rmultinom_size1(PA, nl);
                 newA[samp_A] ++;
                 G[ind][(l*2)+1] = samp_A;

      break;

      case 1: // consistent heterozygous genotypes
             
                a1 = obs_G[0][0];             // observed first allele    
                a2 = obs_G[0][1];             // observed second allele    
                u = PA[a1];                 // frequency of first allele 
                v = PA[a2];                 // frequency of second allele 

                 z = u+v;
                 PA[a1]=0.0;   
                 PA[a2]=0.0;
                 Ppart[0] = 2.0*u*v; 
                 Ppart[1] = (u*u)+(v*v);
                 Ppart[2] = 2.0*z*(1.0-z); 
                 Ppart[3] = (1.0-z)*(1.0-z); 
                 for(ns= 0; ns < n; ns++){    // iterates through samples because cat may be different
                  sample_cat = cat_tmp[ns];
                  Ppart[0] *= 1.0+(E_mat[l][sample_cat]/(2.0*u*v*E_mat[l][1+sample_cat]));
                 }
                 switch(rmultinom_size1(Ppart, 4)){
                     case 0:
                        G[ind][l*2] = a1;
                        G[ind][(l*2)+1] = a2;
                        newA[a1] ++;
                        newA[a2] ++;
                        break;
                     case 1:
                        if((int)rbinom(1.0, (u*u)/((u*u)+(v*v)))==1){ 
                         G[ind][l*2] = a1;
                         G[ind][(l*2)+1] = a1;
                         newA[a1] +=2.0;
                        }else{
                         G[ind][l*2] = a2;
                         G[ind][(l*2)+1] = a2;
                         newA[a2] +=2.0;
                        }
                        break;
                     case 2:
                       samp_A = rmultinom_size1(PA, nl);
                       G[ind][(l*2)] = samp_A;
                       newA[samp_A] ++;
                       if((int)rbinom(1.0, u/z)==1){ 
                        G[ind][(l*2)+1] = a1;
                        newA[a1] ++;
                       }else{
                        G[ind][(l*2)+1] = a2;
                        newA[a2] ++;
                       }
                       break;
                     case 3:
                       samp_A = rmultinom_size1(PA, nl);
                       newA[samp_A]++;
                       G[ind][(l*2)] = samp_A;
                       samp_A = rmultinom_size1(PA, nl);
                       newA[samp_A]++;
                       G[ind][(l*2)+1] = samp_A;
                     break;
                 }   // end switch statement
                 PA[a1] = u;
                 PA[a2] = v;
      break;

      case 2: // consistent homozygous genotype(s)

                
                 a1 = obs_G[0][0];             // observed genotype               
                 u = PA[a1];                 // frequency of first allele 
                  
                 PA[a1] = 0.0;
                 Ppart[0] = u*u; 
                 Ppart[1] = 2.0*u*(1.0-u); 
                 Ppart[2] = (1.0-u)*(1.0-u); 
                 for(ns= 0; ns < n; ns++){    // iterates through samples because cat may be different
                  sample_cat = cat_tmp[ns];
                  Ppart[0] *= 1.0+(E_mat[l][sample_cat]/(u*u*E_mat[l][1+sample_cat]));
                 }
                  
                   switch(rmultinom_size1(Ppart, 3)){
                     case 0:  
                        G[ind][l*2] = a1;
                        G[ind][(l*2)+1] = a1;
                        newA[a1]+=2.0;
                        break;
                     case 1:
                        G[ind][l*2] = a1;
                        newA[a1]++;
                        samp_A = rmultinom_size1(PA, nl);
                        G[ind][(l*2)+1] = samp_A;
                        newA[samp_A]++;
                        break;
                     case 2:
                        samp_A = rmultinom_size1(PA, nl);
                        G[ind][(l*2)] = samp_A;
                        newA[samp_A]++;
                        samp_A = rmultinom_size1(PA, nl);
                        G[ind][(l*2)+1] = samp_A;
                        newA[samp_A]++;
                        break;
                   }
                PA[a1] = u;          
      break;


      case 3: // inconsistent genotypes //
           
              for(itt_set= 0; itt_set < nl; itt_set++){
               all_set[itt_set]=0;      // initialise sets
              }

              for(ns= 0; ns < n; ns++){
               all_set[obs_G[ns][0]]=1;
               all_set[obs_G[ns][1]]=1;
              }

              set_size = 0;

              for(itt_set= 0; itt_set < nl; itt_set++){
                 if(all_set[itt_set]==1){
                  unique_set[set_size]=itt_set;
                  set_size ++;
                 }
              }  // gets the unique set of alleles and there length

               // if the set of ind's observed alleles is 2 //

               if(set_size==2){        
                a1 = unique_set[0];
                a2 = unique_set[1];
                u = PA[a1];
                v = PA[a2];
                z = u+v;
                Ppart[0] = 2.0*u*v;
                Ppart[1] = u*u;
                Ppart[2] = v*v;
                Ppart[3] = 2.0*u*(1.0-z);
                Ppart[4] = 2.0*v*(1.0-z);
                Ppart[5] = (1.0-z)*(1.0-z); 

                for(ns=0; ns < n; ns++){               // itterate through samples 
                   sample_cat = cat_tmp[ns];
                   if(obs_G[ns][0]==a1 && obs_G[ns][1]==a1){       // if obs is homozygote a1/a1 
                     Ppart[1] *=1.0+(E_mat[l][sample_cat]/(u*u*E_mat[l][1+sample_cat])); 
                   }
                   if(obs_G[ns][0]==a2 && obs_G[ns][1]==a2){       // if obs is homozygote a2/a2
                     Ppart[2] *=1.0+(E_mat[l][sample_cat]/(v*v*E_mat[l][1+sample_cat]));  
                   }        
                   if((obs_G[ns][0]==a1 && obs_G[ns][1]==a2) || (obs_G[ns][0]==a2 && obs_G[ns][1]==a1)){   // if obs is heterozygote
                     Ppart[0] *=1.0+(E_mat[l][sample_cat]/(2.0*u*v*E_mat[l][1+sample_cat]));            
                   }                                     // end het/hom ifelse statement
                }                                       // end sample loop

              PA[a1]=0.0;  
              PA[a2]=0.0;
              switch(rmultinom_size1(Ppart, 6)){
               case 0:
                G[ind][l*2] = a1;
                G[ind][(l*2)+1] = a2;
                newA[a1]++;
                newA[a2]++;
               break;
               case 1:
                G[ind][l*2] = a1;
                G[ind][(l*2)+1] = a1;
                newA[a1]+=2.0;
               break;
               case 2:
                G[ind][l*2] = a2;
                G[ind][(l*2)+1] = a2;
                newA[a2]+=2.0;
               break;
               case 3:
               G[ind][l*2] = a1;
               newA[a1]++;
               G[ind][(l*2)+1] = rmultinom_size1(PA, nl);
               samp_A = rmultinom_size1(PA, nl);
               G[ind][(l*2)+1] = samp_A;
               newA[samp_A]++;
               break;
               case 4:
               G[ind][l*2] = a2;
               newA[a2]++;
               G[ind][(l*2)+1] = rmultinom_size1(PA, nl);
               samp_A = rmultinom_size1(PA, nl);
               G[ind][(l*2)+1] = samp_A;
               newA[samp_A]++;
               break;
               case 5:
               samp_A = rmultinom_size1(PA, nl);
               G[ind][(l*2)] = samp_A;
               newA[samp_A]++;
               samp_A = rmultinom_size1(PA, nl);
               G[ind][(l*2)+1] = samp_A;
               newA[samp_A]++;
              } //end switch statement
              PA[a1] = u;
              PA[a2] = v;
              } // end allele set size==2 elseif statement 

              // if the set of ind's observed alleles is greater than 2 //

              if(set_size>2){   
               cnt = 0;
               z=0.0;
                 for(itt_set= 0; itt_set < set_size; itt_set++){
                  a1 = unique_set[itt_set];
                  u = PA[a1];
                  PA[a1] = 0.0;
                    for(na1= 0; na1 < nl; na1++){
                      Ppart[cnt] = u*A[l][na1];
                      if(a1!=na1){Ppart[cnt] *=2.0;}             // if heterozygous multiply by 2
                      cnt++;                          
                    }
                  z += u;                            // cumulative allele frequencies
                 }
                 for(itt_set= 0; itt_set < set_size; itt_set++){
                   a1 = unique_set[itt_set];
                   for(na1= 0; na1 < set_size; na1++){
                     if(itt_set!=na1){
                       Ppart[(na1*nl)+a1] *=0.5;
                     }
                   }        // since heterozygotes from the unique set appear twice 
                 }

                 Ppart[set_size*nl] = (1-z)*(1-z);  // genotypes not containing sampled alleles

                 for(ns=0; ns < n; ns++){
                    a1 = obs_G[ns][0];
                    a2 = obs_G[ns][1];
                    sample_cat = cat_tmp[ns];
                    cnt = 0;
                    for(itt_set= 0; itt_set < set_size; itt_set++){
                     na2 = unique_set[itt_set];
                       for(na1= 0; na1 < nl; na1++){ 
                           if((a1==na1 && a2==na2) || (a2==na1 && a1==na2)){  
                              if(a1==a2){
                                Ppart[cnt] *=1.0+(E_mat[l][sample_cat]/(A[l][a1]*A[l][a1]*E_mat[l][1+sample_cat]));     
                              }else{
                                Ppart[cnt] *= 1.0+(E_mat[l][sample_cat]/(2.0*A[l][a1]*A[l][a2]*E_mat[l][1+sample_cat]));  
                              }
                           }                       
                           cnt ++;
                        }
                    }
                }  //  end of sample within individual loop

                na2 = rmultinom_size1(Ppart, (set_size*nl)+1);
                if(na2!=set_size*nl){   // if true genotype does not contain sampled alleles
                  samp_A = unique_set[na2/nl];
                  newA[samp_A]++;
                  G[ind][l*2]=samp_A;
                  samp_A = na2-((na2/nl)*nl);
                  newA[samp_A]++;
                  G[ind][(l*2)+1] = samp_A;
                }else{
                 samp_A =  rmultinom_size1(PA, nl);
                 G[ind][l*2] = samp_A;
                 newA[samp_A]++;
                 samp_A =  rmultinom_size1(PA, nl);
                 G[ind][(l*2)+1] = samp_A;           
                 newA[samp_A]++;
                }

                for(itt_set= 0; itt_set < set_size; itt_set++){
                 a1 = unique_set[itt_set];
                 PA[a1] = A[l][a1];  // reset temporary allele frequency
                }
              }    
        // end set_size > 2 ifesle statement  
              break;
          }
break;
               
case 1: // contextual information provided by a single parent//

         if(da[0]==-999){
          da[0]=sa[0];
          da[1]=sa[1];
         }

         G[ind][l*2] = da[(int)rbinom(1.0, 0.5)];
         samp_A = rmultinom_size1(PA, nl); 
         G[ind][(l*2)+1] = samp_A;       
         newA[samp_A] ++;

//      gen_status is always zero.  May need to extend to gen_status>0 for rel_status=1 to make algorithm more efficient.
//      see older versions (<2.1) at this place for some old unfinished code.

break;

case 2: // contextual information provided by both parents only /

            Ppart[0]=0.25; 
            Ppart[1]=0.25; 
            Ppart[2]=0.25;
            Ppart[3]=0.25;
            cnt = 0;

          if(gen_status!=0){ // observed genotypes//

            for(S_1=0; S_1<2; S_1++){
              for(S_2=0; S_2<2; S_2++){  
                na1=da[S_1];
                na2=sa[S_2];
                for(ns=0; ns < n; ns++){               // itterate through samples 
                  sample_cat = cat_tmp[ns];
                  a1 = obs_G[ns][0];
                  a2 = obs_G[ns][1];
                  if((a1==na1 && a2==na2) || (a2==na1 && a1==na2)){  
                    if(a1==a2){
                      Ppart[cnt] *= 1.0+(E_mat[l][sample_cat]/(A[l][a1]*A[l][a1]*E_mat[l][1+sample_cat]));     
                    }else{
                      Ppart[cnt] *= 1.0+(E_mat[l][sample_cat]/(2.0*A[l][a1]*A[l][a2]*E_mat[l][1+sample_cat]));  
                    }
                  }
                }
                cnt++;
              }
            }
          }

          cnt = rmultinom_size1(Ppart, 4);
          G[ind][l*2] = da[cnt/2];
          if(cnt==0 || cnt==2){
            G[ind][(l*2)+1] = sa[0];
          }else{
            G[ind][(l*2)+1] = sa[1];
          }

break;     

case 3:

// both parents recorded and offspring - only 4 possible (ordered) genotypes with uniform prior//

             Ppart[0] = 1.0;
             Ppart[1] = 1.0;
             Ppart[2] = 1.0;
             Ppart[3] = 1.0;

             for(itt_set=0; itt_set < no_off[ind]; itt_set++){               

               o = par[ind][itt_set];

               off_1 = G[o][l*2];       // offspring and spouse genotypes //
               off_2 = G[o][(l*2)+1];

               sp_1 = -999;
               sp_2 = -999;

               if(sire[o]==dam[o]){      // new inbreeding code - 28/3/08 //
                 sp_1 = -998;
                 sp_2 = -998;
               }else{
                 if(sire[o]!=ind){
                   if(sire[o]<nind){
                     sp_1 = G[sire[o]][l*2];
                     sp_2 = G[sire[o]][(l*2)+1];
                   }
                 }else{
                   if(dam[o]<nind){
                     sp_1 = G[dam[o]][l*2];
                     sp_2 = G[dam[o]][(l*2)+1];
                   }
                 }
               }

               u=PA[off_1];
               v=PA[off_2];
               cnt=0;
               for(S_1=0; S_1<2; S_1++){
                 for(S_2=0; S_2<2; S_2++){  
                   if(sp_1>(-998)){
                     tac = 0.0;
                     if((da[S_1]==off_1 && off_2==sp_1) || (da[S_1]==off_2 && off_1==sp_1)){tac+=0.25;}
                     if((da[S_1]==off_1 && off_2==sp_2) || (da[S_1]==off_2 && off_1==sp_2)){tac+=0.25;}
                     if((sa[S_2]==off_1 && off_2==sp_1) || (sa[S_2]==off_2 && off_1==sp_1)){tac+=0.25;}
                     if((sa[S_2]==off_1 && off_2==sp_2) || (sa[S_2]==off_2 && off_1==sp_2)){tac+=0.25;}
                     Ppart[cnt] *= tac;
                   }else{
                     if(sp_1==-999){
                       if(off_1==off_2){
                         Ppart[cnt] *= (int(da[S_1]==off_1)+int(sa[S_2]==off_1))*u;
                       }else{
                         z=0.0;
                         if(da[S_1]==off_1 || sa[S_2]==off_1){
                           z+=v;
                         }
                         if(da[S_1]==off_2 || sa[S_2]==off_2){
                           z+=u;
                           if(da[S_1]==off_1 || sa[S_2]==off_1){
                             z/=2.0;
                           }
                         }                                    
                         Ppart[cnt] *= z; 
                       }
                     }else{
                       if(off_1!=IBL[0]){
                         if(IBL[0]==-999){
                           IBL[0]=off_1;
                           IBL[1]=-1;
                         }else{
                           IBL[1]=off_1;
                         }
                       }else{
                         if(IBL[1]<0){
                           IBL[1]--;
                         }
                       }
                     }
                   }
                   cnt++;
                 }
               }
             }          

             if(IBL[0]!=-999){
               if(IBL[1]<0){
                 Ppart[0] *= pow(double((da[0]==IBL[0])+(sa[0]==IBL[0])), -2.0*IBL[1]);
                 Ppart[1] *= pow(double((da[1]==IBL[0])+(sa[0]==IBL[0])), -2.0*IBL[1]);
                 Ppart[2] *= pow(double((da[0]==IBL[0])+(sa[1]==IBL[0])), -2.0*IBL[1]);
                 Ppart[3] *= pow(double((da[1]==IBL[0])+(sa[1]==IBL[0])), -2.0*IBL[1]);
               }else{
                 if(da[0]!=IBL[0] && da[0]!=IBL[1]){
                   Ppart[0] = 0.0;
                   Ppart[2] = 0.0;
                 }
                 if(da[1]!=IBL[0] && da[1]!=IBL[1]){
                   Ppart[1] = 0.0;
                   Ppart[3] = 0.0;
                 }
                 if(sa[0]!=IBL[0] && sa[0]!=IBL[1]){
                   Ppart[0] = 0.0;
                   Ppart[1] = 0.0;
                 }
                 if(sa[1]!=IBL[0] && sa[1]!=IBL[1]){
                   Ppart[2] = 0.0;
                   Ppart[3] = 0.0;
                 }
               }
             }

             IBL[0]=-999;
             IBL[1]=-999;

              cnt = 0;

              for(S_1=0; S_1<2; S_1++){
                for(S_2=0; S_2<2; S_2++){ 
                  na1=da[S_1];
                  na2=sa[S_2];                   
                  for(ns=0; ns < n; ns++){               // itterate through samples 
                    sample_cat = cat_tmp[ns];
                    a1 = obs_G[ns][0];
                    a2 = obs_G[ns][1];
                    if((a1==na1 && a2==na2) || (a2==na1 && a1==na2)){  
                      if(a1==a2){
                        Ppart[cnt] *= 1.0+(E_mat[l][sample_cat]/(A[l][a1]*A[l][a1]*E_mat[l][1+sample_cat]));     
                      }else{
                        Ppart[cnt] *= 1.0+(E_mat[l][sample_cat]/(2.0*A[l][a1]*A[l][a2]*E_mat[l][1+sample_cat]));  
                      }
                    }
                  }
                  cnt++;
                }
              }
              cnt = rmultinom_size1(Ppart, 4);
              G[ind][l*2] = da[cnt/2];
              if(cnt==0 || cnt==2){
                G[ind][(l*2)+1] = sa[0];
              }else{
                G[ind][(l*2)+1] = sa[1];
              }

break;
 
case 4: // has offspring but no parents, or one parent //

      if(da[0]!=-999 || sa[0]!=-999){
        no_base = 0;
        if(da[0]==-999){
         da[0]=sa[0];
         da[1]=sa[1];
        }
        S_1 = da[0];
        if(da[0]==da[1]){           // if parent is homozygous must be da[0]/something with prior = A
          S_2 = -999;               // n(S)=1
          for(o=0; o<nl; o++){               
            vec[0][o] = PA[o];
          }
        }else{                     // if parent is heterozygous must be da[0]/something or da[1]/something with priors = A/2
          S_2 = da[1];             // n(S)=2
          for(o=0; o<nl; o++){               
            vec[0][o] = PA[o];
            vec[1][o] = PA[o];
           }
           vec[0][S_2] = PA[S_1]+PA[S_2];
           vec[1][S_1] = PA[S_1]+PA[S_2];
         }
       }else{       // neither parent known
          no_base = -999;
          o = par[ind][0];
          S_1 = G[o][l*2]; 
          S_2 = G[o][(l*2)+1]; 
          if(S_1==S_2){  // if offspring is HOM must be off[0]/something with prior = A*u
            S_2 = -999;
            u = PA[S_1];
            for(o=0; o<nl; o++){               
             vec[0][o] = PA[o]*u;
            }
            vec[0][S_1] *= 0.5;   // divide HOM frequency by 2 (= multiplying all HET by 2)
          }else{             // if offspring is HET must be off[0]/something or off[1]/something with prior = A*u/A*v
            u = PA[S_1];
            v = PA[S_2];
            for(o=0; o<nl; o++){              
             vec[0][o] = PA[o]*u;
             vec[1][o] = PA[o]*v;
            }
            vec[0][S_1] *= 0.5;   // divide homozygous genotype frequency by 2 (= multiplying all HET by 2)
            vec[1][S_2] *= 0.5;   
          }                        // don't need to worry about the fact vec[0][S_2] = vec[1][S_1]
       }

        which_vec=0;

// itterate through offspring/spouses if they exist //

       for(itt_set=0; itt_set< no_off[ind]; itt_set++){          

         if(S_1==-999){break;}             // if n(S) = 0 G already known so finish

         o = par[ind][itt_set];

         off_1 = G[o][l*2];       // offspring and spouse genotypes //
         off_2 = G[o][(l*2)+1];

         sp_1 = -999;
         sp_2 = -999;

         if(sire[o]==dam[o]){      // new inbreeding code - 28/3/08 //
           sp_1 = -998;
           sp_2 = -998;
         }else{
           if(sire[o]!=ind){
             if(sire[o]<nind){
               sp_1 = G[sire[o]][l*2];
               sp_2 = G[sire[o]][(l*2)+1];
             }
           }else{
             if(dam[o]<nind){
               sp_1 = G[dam[o]][l*2];
               sp_2 = G[dam[o]][(l*2)+1];
             }
           }
         }

         offhom_sizeS = (2*(off_1!=off_2))+int(S_2!=-999); // 0: off=HOM and n(s)=1, 1: off=HOM and n(s)=2
                                                           // 2: off=HET and n(s)=1, 3: off=HET and n(s)=2

         if(sp_1==-998 && (offhom_sizeS==2 || offhom_sizeS==3)){   // selfers that produces HET offspring must also be HET.
            G[ind][(l*2)] = off_1;
            G[ind][(l*2)+1] = off_2;
            if(no_base==-999){         
              newA[off_1]++;
              newA[off_2]++; 
            }
            S_1=-999;
            break;
         }

         switch(offhom_sizeS){

// off HOM and n(s)=1//
         case 0:         
           if(off_1!=S_1){       // if off[0] nin S then genotype is known
             if(no_base==-999){         
              newA[S_1]++;
              newA[off_1]++; 
             }else{
               if(da[0]==off_1 || da[1]==off_1){
                 if(da[0]==S_1 || da[1]==S_1){
                   newA[off_1]+=0.5;
                   newA[S_1]+=0.5;
                 }else{
                   newA[S_1]++;
                 }
               }else{
                 newA[off_1]++;
               }
             }
             G[ind][(l*2)] = S_1;
             G[ind][(l*2)+1] = off_1;
             S_1 = -999;
           }else{                                             // if off[0] in S then multiply homozygote pairing by 2 because:
             vec[which_vec][off_1] *=2.0*(1.0+(sp_1==-998));  // if sp is hom off_1/off_2 then HOM Tac=1 and HET Tac=0.5
           }                                                  // if sp is het off_1/sp_2 then HOM Tac=0.5 and HET Tac=0.25
           break;                                             // for hermaphroditic systems multiply homozygotes by 4 

// off HOM and n(s)=2 we can reduce the set size //

         case 1:
          if(off_1==S_1 || off_1==S_2){
            if(off_1==S_1){  
             which_vec=0;                     
             S_2=-999;
            }else{
             which_vec=1;           
             S_1=S_2;
             S_2=-999;
            } 
            vec[which_vec][S_1] *=2.0*(1.0+(sp_1==-998));   // see above 
          }else{
            which_vec=0;
            u = vec[0][off_1];
            v = vec[1][off_1];
            for(o=0; o < nl; o++){              
             vec[0][o] = 0.0;
            }
            vec[0][S_1]=u;
            vec[0][S_2]=v;       
            S_1 = off_1;
            S_2 = -999;
          }  
        break;


// off HET and n(s)=1//  

         case 2:  

          if(off_1 == S_1 || off_2 ==S_1){  // if either off alleles in S //
            if(sp_1!=-999){                 // if spouse is sampled //
                if((sp_1== off_1 && sp_2== off_2) || (sp_1== off_2 && sp_2== off_1)){   
                   vec[which_vec][off_1] *=2.0; 
                   vec[which_vec][off_2] *=2.0;
                }else{   // if sp differes from off then we know an ind allele for sure // 
                  if(off_1!=sp_1 && off_1!=sp_2){ // if off_1
                    if(off_1!=S_1){               // and off_1 is not in S 
                      G[ind][(l*2)] = S_1;
                      G[ind][(l*2)+1] = off_1;
                      if(no_base==-999){         
                        newA[S_1]++;
                        newA[off_1]++; 
                      }else{
                        if(da[0]==off_1 || da[1]==off_1){
                          if(da[0]==S_1 || da[1]==S_1){
                            newA[off_1]+=0.5;
                            newA[S_1]+=0.5;
                          }else{
                            newA[S_1]++;
                          }
                        }else{
                          newA[off_1]++;
                        }
                      }
                      S_1 = -999;
                    }else{                        // if off_1 is in S then sp is HET off_2/something
                      vec[which_vec][S_1] *=2.0;   // Tac(S_1/something | S_1/S_1, something/something_else) =0.5
                    }                             // Tac(S_1/something | S_1/something, something/something_else) =0.25
                  }else{                     // Tac(S_1/something | S_1/something_else, something/something_else) =0.25
                    if(off_2!=S_1){              
                      G[ind][(l*2)] = S_1;
                      G[ind][(l*2)+1] = off_2;
                      if(no_base==-999){         
                        newA[S_1]++;
                        newA[off_2]++; 
                      }else{
                        if(da[0]==off_2 || da[1]==off_2){
                          if(da[0]==S_1 || da[1]==S_1){
                            newA[off_2]+=0.5;
                            newA[S_1]+=0.5;
                          }else{
                            newA[S_1]++;
                          }
                        }else{
                          newA[off_2]++;
                        }
                      }
                      S_1 = -999;
                    }else{
                      vec[which_vec][S_1] *=2.0;
                    }
                  }
                }
            }else{                       // spouse is not sampled //
              u = PA[off_1];
              v = PA[off_2];
              vec[which_vec][S_1] *= 2.0;
              if(off_1!=S_1){ // if off_1
                vec[which_vec][off_1] *= (u+v)/u;
              }else{
                vec[which_vec][off_2] *= (u+v)/v;
              }
            }
          }else{                         // if off alleles not in S //
            if(sp_1!=-999){                 // if spouse is sampled //
              if((off_1==sp_1 && off_2== sp_2) || (off_1== sp_2 && off_2== sp_1)){   
                u = vec[which_vec][off_1];
                v = vec[which_vec][off_2];
                for(o=0; o < nl; o++){              
                 vec[which_vec][o] = 0.0;
                }
                 vec[which_vec][off_1]=u;
                 vec[which_vec][off_2]=v;  
              }else{
                if(sp_1!=off_1 && sp_2!=off_1){
                  G[ind][(l*2)] = S_1;
                  G[ind][(l*2)+1] = off_1;
             if(no_base==-999){         
              newA[S_1]++;
              newA[off_1]++; 
             }else{
               if(da[0]==off_1 || da[1]==off_1){
                 if(da[0]==S_1 || da[1]==S_1){
                   newA[off_1]+=0.5;
                   newA[S_1]+=0.5;
                 }else{
                   newA[S_1]++;
                 }
               }else{
                 newA[off_1]++;
               }
             }  
                }else{
                  G[ind][(l*2)] = S_1;
                  G[ind][(l*2)+1] = off_2;
             if(no_base==-999){         
              newA[S_1]++;
              newA[off_2]++; 
             }else{
               if(da[0]==off_2 || da[1]==off_2){
                 if(da[0]==S_1 || da[1]==S_1){
                   newA[off_2]+=0.5;
                   newA[S_1]+=0.5;
                 }else{
                   newA[S_1]++;
                 }
               }else{
                 newA[off_2]++;
               }
             } 
                }  
                 S_1 = -999;                      
              }
            }else{                        // spouse is not sampled //
              u = vec[which_vec][off_1];
              v = vec[which_vec][off_2]; 
              for(o=0; o < nl; o++){              
               vec[which_vec][o] = 0.0;
              }
              vec[which_vec][off_1] = u*PA[off_2];
              vec[which_vec][off_2] = v*PA[off_1]; 
            }                    
          } 
          break;   

          // off HET and n(s)=2//  
          case 3:   

            // IF spouse IS sampled //

            no_all_in_S = int(S_1==off_1)+int(S_1==off_2)+int(S_2==off_1)+int(S_2==off_2);

            switch(no_all_in_S){

            case 0:

            if(sp_1!=-999){  // spouse is sampled //
              if((sp_1== off_1 && sp_2== off_2) || (sp_1== off_2 && sp_2== off_1)){   
                 u = vec[0][off_1];
                 v = vec[0][off_2];
                 for(o=0; o < nl; o++){               
                   vec[0][o] = 0.0;                        // all multiplied by u//
                 }
                 vec[0][off_1]=u;
                 vec[0][off_2]=v;
                 u = vec[1][off_1];
                 v = vec[1][off_2];
                 for(o=0; o < nl; o++){               
                   vec[1][o] = 0.0;                        // all multiplied by u//
                 }
                 vec[1][off_1]=u;
                 vec[1][off_2]=v;
               }else{   // if sp differes from off then we know an ind allele for sure // 
                 if(off_1!=sp_1 && off_1!=sp_2){ // if off_1
                   which_vec=0;
                   u = vec[0][off_1];
                   v = vec[1][off_1];
                   for(o=0; o < nl; o++){               
                     vec[0][o] = 0.0;                        // all multiplied by u//
                   }
                   vec[0][S_1]=u;
                   vec[0][S_2]=v;
                   S_1=off_1;
                   S_2=-999;
                 }else{                        // if off_1 is in S then sp is HET off_2/something
                   which_vec=0;
                   u = vec[0][off_2];
                   v = vec[1][off_2];
                   for(o=0; o < nl; o++){               
                     vec[0][o] = 0.0;                        // all multiplied by u//
                   }
                   vec[0][S_1]=u;
                   vec[0][S_2]=v;
                   S_1=off_2;
                   S_2=-999;
                 }
               }
             }else{     // spouse is unsampled //
                   u = vec[0][off_1];
                   v = vec[0][off_2];
                   for(o=0; o < nl; o++){               
                     vec[0][o] = 0.0;                        // all multiplied by u//
                   }
                   vec[0][off_1] = u*PA[off_2]; 
                   vec[0][off_2] = v*PA[off_1]; 

                   u = vec[1][off_1];
                   v = vec[1][off_2];
                   for(o=0; o < nl; o++){               
                     vec[1][o] = 0.0;                        // all multiplied by u//
                   }
                   vec[1][off_1] = u*PA[off_2]; 
                   vec[1][off_2] = v*PA[off_1]; 
             }
            break;

            case 1:
            if(sp_1!=-999){  // spouse is sampled //
              if((sp_1==off_1 && sp_2== off_2) || (sp_1== off_2 && sp_2== off_1)){   
                if(S_1!=off_1 && S_1!=off_2){
                  u = vec[0][off_1];
                  v = vec[0][off_2];
                  for(o=0; o < nl; o++){               
                   vec[0][o] = 0.0;                        // all multiplied by u//
                  }
                  vec[0][off_1]=u;
                  vec[0][off_2]=v; 
                  vec[1][off_1]*=2.0;
                  vec[1][off_2]*=2.0;
                }else{
                  u = vec[1][off_1];
                  v = vec[1][off_2];
                  for(o=0; o < nl; o++){               
                   vec[1][o] = 0.0;                        // all multiplied by u//
                  }
                  vec[1][off_1]=u;
                  vec[1][off_2]=v;       
                  vec[0][off_1]*=2.0;
                  vec[0][off_2]*=2.0;
                }
              }else{
                if(off_1!=sp_1 && off_1!=sp_2){ // if off_1
                   if(off_1==S_1 || off_1==S_2){
                     if(off_1==S_1){  
                      which_vec=0;                           
                      S_2=-999;
                     }else{
                      which_vec=1;
                      S_1=S_2;
                      S_2=-999;  
                     } 
                     vec[which_vec][S_1]*=2.0;
                   }else{
                     which_vec=0;
                     u = vec[0][off_1];
                     v = vec[1][off_1];
                     for(o=0; o < nl; o++){               
                      vec[0][o] = 0.0;                        // all multiplied by u//
                     }
                     vec[0][S_1]=u;
                     vec[0][S_2]=v;       
                     S_1=off_1;
                     S_2=-999;
                   }
                }else{                     // Tac(S_1/something | S_1/something_else, something/something_else) =0.25
                   if(off_2==S_1 || off_2==S_2){
                     if(off_2==S_1){  
                      which_vec=0;                           
                      S_2=-999;
                     }else{
                      which_vec=1;
                      S_1=S_2;
                      S_2=-999;
                     } 
                     vec[which_vec][S_1]*=2.0;
                   }else{
                     which_vec=0;
                     u = vec[0][off_2];
                     v = vec[1][off_2];
                     for(o=0; o < nl; o++){               
                      vec[0][o] = 0.0;                        // all multiplied by u//
                     }
                     vec[0][S_1]=u;
                     vec[0][S_2]=v;       
                     S_1=off_2;
                     S_2=-999;
                   }
                }
              }
            }else{ // spouse is unsampled //
              if(off_1==S_1 || off_1==S_2){
              u = PA[off_1];
              v = PA[off_2];
                if(off_1==S_1){  
                   vec[0][off_1] *=2.0;
                   vec[0][off_2] *= (u+v)/v;
                   v =  vec[1][off_2]*(u/v);
                   u =  vec[1][off_1];
                   for(o=0; o < nl; o++){               
                      vec[1][o] = 0.0;                        // all multiplied by u//
                   }
                   vec[1][off_1] = u;
                   vec[1][off_2] = v;
                }else{
                   vec[1][off_1] *=2.0;
                   vec[1][off_2] *= (u+v)/v;
                   v =  vec[0][off_2]*(u/v);
                   u =  vec[0][off_1];
                   for(o=0; o < nl; o++){               
                      vec[0][o] = 0.0;                        // all multiplied by u//
                   }
                   vec[0][off_1] = u;
                   vec[0][off_2] = v;
                }
              }else{
              v = PA[off_1];
              u = PA[off_2];
                if(off_2==S_1){  
                   vec[0][off_2] *=2.0;
                   vec[0][off_1] *= (u+v)/v;
                   v =  vec[1][off_1]*(u/v);
                   u =  vec[1][off_2];
                   for(o=0; o < nl; o++){               
                      vec[1][o] = 0.0;                        // all multiplied by u//
                   }
                   vec[1][off_2] = u;
                   vec[1][off_1] = v;
                }else{
                   vec[1][off_2] *=2.0;
                   vec[1][off_1] *= (u+v)/v;
                   v =  vec[0][off_1]*(u/v);
                   u =  vec[0][off_2];
                   for(o=0; o < nl; o++){               
                      vec[0][o] = 0.0;                        // all multiplied by u//
                   }
                   vec[0][off_2] = u;
                   vec[0][off_1] = v;
                }              }
            }

            break;

            case 2:
            if(sp_1!=-999){  // spouse is sampled //
              if((sp_1== off_1 && sp_2== off_2) || (sp_1== off_2 && sp_2== off_1)){   // if spouse = offspring  //  
                 vec[0][S_1] *=2.0;                                 // Tac(S_1/S_2|S_1/S_2|S_1/S_2) =0.5
                 vec[0][S_2] *=2.0;                                 // Tac(S_1/S_2|S_1/S_1|S_1/S_2) =0.5
                 vec[1][S_1] *=2.0;                                 // Tac(S_1/S_2|S_2/S_something|S_1/S_2) =0.25
                 vec[1][S_2] *=2.0; 
              }else{                                                                  // if spouse neq offspring  //
                if(off_1!=sp_1 && off_1!=sp_2){ // off_1 is from ind 
                  if(off_1==S_1){
                   which_vec=0;
                   S_2=-999;  
                  }else{
                   which_vec=1;
                   S_1=S_2;
                   S_2=-999;
                  } 
                }else{                       // off_2 is from ind 
                  if(off_2==S_1){
                    which_vec=0;
                    S_2=-999;  
                  }else{
                    which_vec=1;
                    S_1=S_2;
                    S_2=-999;   
                  }
                }
                vec[which_vec][S_1] *=2.0; 
              }
            }else{
             u = PA[S_1];
             v = PA[S_2];
                 for(o=0; o < nl; o++){   
                  vec[0][o] *= v/2.0; 
                  vec[1][o] *= u/2.0;
                 }
                  vec[0][S_1] *= 2.0; 
                  vec[0][S_2] *= (u+v)/v; 
                  vec[1][S_2] *= 2.0; 
                  vec[1][S_1] *= (u+v)/u; 
            }
            break;
            }
          break;
          } // end switch off HET/HOM n(S)1/2 
        } // end spouse/offspring loop

        // arbitrary priors set now add info from genotype data //

        if(S_1!=-999){     // if S_1 is -999 genotype already known
          if(S_2==-999){       // n(S)=1
            for(o=0; o<nl; o++){
              Ppart[o] = vec[which_vec][o];
            }
            for(ns=0; ns < n; ns++){               
              sample_cat = cat_tmp[ns];  
              a1 =  obs_G[ns][0];
              a2 =  obs_G[ns][1];  
              if(a1==a2){
                if(a1==S_1){   
                  Ppart[S_1] *= 1.0+(E_mat[l][sample_cat]/(A[l][a1]*A[l][a1]*E_mat[l][1+sample_cat]));     
                }
              }else{
                if(a1==S_1){   
                  Ppart[a2] *= 1.0+(E_mat[l][sample_cat]/(2.0*A[l][a1]*A[l][a2]*E_mat[l][1+sample_cat]));    
                }
                if(a2==S_1){  
                  Ppart[a1] *= 1.0+(E_mat[l][sample_cat]/(2.0*A[l][a1]*A[l][a2]*E_mat[l][1+sample_cat])); 
                }
              }
            }

            G[ind][(l*2)] = S_1;
            samp_A = rmultinom_size1(Ppart, nl);
            G[ind][(l*2)+1] = samp_A;
            if(no_base==-999){         
              newA[S_1]++;
              newA[samp_A]++; 
            }else{
              if(da[0]==samp_A || da[1]==samp_A){
                if(da[0]==S_1 || da[1]==S_1){
                  newA[samp_A]+=0.5;
                  newA[S_1]+=0.5;
                }else{
                  newA[S_1]++;
                }
              }else{
                newA[samp_A]++;
              }
            }
          }else{         // n(S)=2
            vec[1][S_1]=0.0;  
            for(o=0; o<nl; o++){
               Ppart[o] = vec[0][o];
               Ppart[nl+o] = vec[1][o];
             }                                              
             for(ns=0; ns < n; ns++){             
               sample_cat = cat_tmp[ns];
               a1 =  obs_G[ns][0];
               a2 =  obs_G[ns][1]; 
               if(a1==a2){
                 if(a1==S_1){
                   Ppart[a1] *= 1.0+(E_mat[l][sample_cat]/(A[l][a1]*A[l][a1]*E_mat[l][1+sample_cat]));    
                 }
                 if(a1==S_2){
                   Ppart[a1+nl] *= 1.0+(E_mat[l][sample_cat]/(A[l][a1]*A[l][a1]*E_mat[l][1+sample_cat]));    
                 }
               }else{
                 if(a1==S_1){
                   Ppart[a2] *= 1.0+(E_mat[l][sample_cat]/(2.0*A[l][a1]*A[l][a2]*E_mat[l][1+sample_cat]));    
                 }
                 if(a2==S_1){
                   Ppart[a1] *= 1.0+(E_mat[l][sample_cat]/(2.0*A[l][a1]*A[l][a2]*E_mat[l][1+sample_cat]));    
                 }
                 if(a1==S_2){
                   Ppart[a2+nl] *= 1.0+(E_mat[l][sample_cat]/(2.0*A[l][a1]*A[l][a2]*E_mat[l][1+sample_cat]));    
                 }
                 if(a2==S_2){
                   Ppart[a1+nl] *= 1.0+(E_mat[l][sample_cat]/(2.0*A[l][a1]*A[l][a2]*E_mat[l][1+sample_cat]));    
                 }
               }
             }
             cnt = rmultinom_size1(Ppart, nl*2);
             if(cnt/nl==0){
               G[ind][(l*2)] = S_1;
               G[ind][(l*2)+1] = cnt;
               if(no_base==-999){         
                newA[S_1]++;
                newA[cnt]++; 
               }else{
                 if(da[0]==cnt || da[1]==cnt){
                   if(da[0]==S_1 || da[1]==S_1){
                     newA[cnt]+=0.5;
                     newA[S_1]+=0.5;
                   }else{
                     newA[S_1]++;
                   }
                 }else{
                   newA[cnt]++;
                 }
               }
             }else{
               G[ind][(l*2)] = S_2;
               G[ind][(l*2)+1] = cnt-nl;
               if(no_base==-999){         
                newA[S_2]++;
                newA[cnt-nl]++; 
               }else{
                 if(da[0]==S_2 || da[1]==S_2){
                   if(da[0]==(cnt-nl) || da[1]==(cnt-nl)){
                     newA[S_2]+=0.5;
                     newA[cnt-nl]+=0.5;
                   }else{
                     newA[cnt-nl]++;
                   }
                 }else{
                   newA[S_2]++;
                 }
               }
             }
           }
         }  
         break;
       }

/****************************************************************/
/*  Some bug-checking code for detecting problems with sampG    */
/*  remember to uncomment the declared variables at the start   */
/*  and also the old geneotypes before sampling on lines 110 -  */
/*  111 and their declarations on line 58-61                    */
/****************************************************************/

/*
       problem= FALSE;

       if(no_off[ind]>0){
         for(j=0; j < no_off[ind]; j++){
           o = par[ind][j];
           off_1 = G[o][l*2];       // offspring and spouse genotypes //
           off_2 = G[o][(l*2)+1];
           if(G[ind][(l*2)]!=off_1 & G[ind][(l*2)]!=off_2 & G[ind][(l*2)+1]!=off_1 & G[ind][(l*2)+1]!=off_2){
             problem = TRUE;
           }
         }
       }
       if(G[ind][(l*2)]>=nl | G[ind][(l*2)+1]>=nl){
         problem = TRUE;
       }
       if(da[0]!=-999){
         if(G[ind][(l*2)]!=da[0] & G[ind][(l*2)]!=da[1] & G[ind][(l*2)+1]!=da[0] & G[ind][(l*2)+1]!=da[1]){
           problem = TRUE;
         }
       }
       if(problem==TRUE){
         cout << "something is wrong at" << endl;
         cout << "locus" << endl;
         cout << l << endl;
         cout << "for individual" << endl;
         cout << ind << endl;
         cout << "The old genotype was:" << endl;
         cout << oldG1 << endl;
         cout << oldG2 << endl;
         cout << "The new genotype is:" << endl;
         cout << G[ind][(l*2)] << endl;
         cout << G[ind][(l*2)+1] << endl;
         cout << "The observed genotypes are:" << endl;
         for(ns=0; ns < n; ns++){               
           cout <<  obs_G[ns][0] << endl;
           cout <<  obs_G[ns][1] << endl;
         }
        cout << "gen_status is:" << endl;
         cout << gen_status<< endl;
         cout << "rel status is:" << endl;
         cout << rel_status << endl;


         cout << "S_1 is recorded as:" << endl;
         cout << S_1 << endl;
         cout << "S_2 is recorded as:" << endl;
         cout << S_2 << endl;
         if(da[0]!=-999){
           cout << "A parent exists" << endl;
           if(G[ind][(l*2)]!=da[0] & G[ind][(l*2)]!=da[1] & G[ind][(l*2)+1]!=da[0] & G[ind][(l*2)+1]!=da[1]){
             cout << "but the sampled genotype is not comaptible with the parent sampled" << endl;
             cout << "The parental genotype are" << endl;
             cout << da[0] << endl;
             cout << da[1] << endl;
             cout << sa[0] << endl;
             cout << sa[1] << endl;
           }else{
             cout << "and the parental genotype is compatible"  << endl;
             cout << "The parental genotypes are" << endl;
             cout << da[0] << endl;
             cout << da[1] << endl;
             cout << sa[0] << endl;
             cout << sa[1] << endl;
           }
         }else{
           cout << "no parent exists"  << endl;
         }
         if(no_off[ind]>0){
           cout << "Offspring exist with genotype" << endl;
           bool problem = FALSE;
           for(j=0; j < no_off[ind]; j++){
             cout << "with genotype" << endl; 
             o = par[ind][j];
             off_1 = G[o][l*2];       // offspring and spouse genotypes //
             off_2 = G[o][(l*2)+1];
             cout << off_1 << endl;
             cout << off_2 << endl;
             sp_1 = -999;
             sp_2 = -999;
             cout << "and spouse" << endl; 
               if(sire[o]!=ind){
                 if(sire[o]<nind){
                   sp_1 = G[sire[o]][l*2];
                   sp_2 = G[sire[o]][(l*2)+1];
                 }
               }else{
                 if(dam[o]<nind){
                   sp_1 = G[dam[o]][l*2];
                   sp_2 = G[dam[o]][(l*2)+1];
                 }
               }
             cout << sp_1 << endl;
             cout << sp_2 << endl;
             if(G[ind][(l*2)]!=off_1 & G[ind][(l*2)]!=off_2 & G[ind][(l*2)+1]!=off_1 & G[ind][(l*2)+1]!=off_2){
               problem = TRUE;
             }
           }
         }else{
           cout << "No Offspring exist" << endl;
         }
       }
*/ 
          ind++;         
          n = 0; 
          misstype = 0;
        }                    // end this is the last record if statement
     }   
                    // end sample loop
      if(estA==TRUE){
        rdirichlet(newA, nl, PA);
        for(i=0; i < nl; i++){  
          A[l][i] = PA[i];
        }
     } 

     q = 1;
     ind = 0;
  }                         // end loci loop
  
        for(i = 0; i < maxrep; ++i){
           delete [] obs_G[i];
        }
        for(i = 0; i < 2; ++i){
           delete [] vec[i];
        }
        for(i = 0; i < nind; ++i){
           delete [] par[i];
        }
  delete [] obs_G;
  delete [] vec;
  delete [] par;


  delete [] cat_tmp;
  delete [] Ppart;
  delete [] all_set;
  delete [] unique_set;
  delete [] PA;       
  delete [] newA;    
  delete [] no_off; 
}
 
