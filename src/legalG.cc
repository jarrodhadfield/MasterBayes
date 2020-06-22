#include "MCMCped.h" 

extern "C"{  

using namespace std;

void legalG(
	int *nindP,		 // number of individuals sampled
        int *damP,             // number of candidate dams per offspring
        int *sireP,            // number of candidate sires per offspring
        int *nlociP,		 // number of loci
        int *nallP,              // number of alleles per locus
        int *maxallP,
	double *st_AP,	         // starting allele frequencies
        int *st_GP,              // starting true genotypes   
        bool *legalP, 
        int *mtypeP
){        

// pointers to single variables are redefined

int 	nind = nindP[0], 	
	nloci = nlociP[0],      
        *nall = nallP,
        maxall = maxallP[0],
        i,
        l,
        o,
        o1,
        o2,
        d1,
        d2,
        s1,
        s2, 
        index=0;

double  phom,
        qhom,
        pqhet;

double  *Pall = new double[maxall];

int mtype = mtypeP[0];

// declare some variable sized arrays

int     *pG,
        **G;
double  *pA,
        **A;
 
pG = new(nothrow) int [(1+int(mtype==1 || mtype==3))*nind*nloci];
if(pG==NULL)
{
 error("NO MEMORY for G\n");
}
G = new(nothrow) int* [nind];
if(G==NULL)
{
 error("NO MEMORY for G\n");
}
pA = new(nothrow) double [nloci*maxall];
if(pA==NULL)
{
 error("NO MEMORY for A\n");
}
A = new(nothrow) double* [nloci];
if(A==NULL)
{
 error("NO MEMORY for A\n");
}

for (i=0; i<nind; ++i){
G[i] = &pG[index];
index  += ((1+int(mtype==1 || mtype==3))*nloci);
}

index = 0;

for (i=0; i<nloci; ++i){
A[i] = &pA[index];
index  += maxall;
}
        
    bool TacP = true;
    bool o1ind;
    bool o1ins;
    bool o2ind;
    bool o2ins;

    read_G(st_GP, nind, nloci, G, mtype);
    read_A(st_AP, nloci, A, nall);

    int** par = new int*[nind];
    for(i = 0; i < nind; ++i){
      par[i] = new int[50];
    }

     int *no_off = new int[nind];
     int itt_set;

     for(i=0; i < nind; i++){
           no_off[i] = 0;
     }

     for(i=0; i < nind; i++){                         
       if(sireP[i]<nind){
         l = sireP[i];
         par[l][no_off[l]] = i;
         no_off[l]++;
       } 
       if(damP[i]<nind && damP[i]!=sireP[i]){
         l = damP[i];
         par[l][no_off[l]] = i;
         no_off[l]++;
       }
     }


  switch(mtype){

    case 1:
    case 3:

    for(i=0; i<nind; ++i){
      for(l=0; l<nloci; ++l){
        if(damP[i]<nind || sireP[i]<nind){

          o1 = G[i][l*2];
          o2 = G[i][(l*2)+1];     

          if(damP[i]<nind && sireP[i]<nind){

            d1 = G[damP[i]][l*2];
            d2 = G[damP[i]][(l*2)+1];     

            s1 = G[sireP[i]][l*2];
            s2 = G[sireP[i]][(l*2)+1];     

            TacP=false;

            if((d1==o1 && s1==o2) || (d1==o2 && s1==o1)){TacP=true;}
            if((d1==o1 && s2==o2) || (d1==o2 && s2==o1)){TacP=true;}
            if((d2==o1 && s1==o2) || (d2==o2 && s1==o1)){TacP=true;}
            if((d2==o1 && s2==o2) || (d2==o2 && s2==o1)){TacP=true;}

            if(TacP==false){

              legalP[0]=false;

              if(o1==d1 || o1==d2){
                o1ind = true;
              }else{
                o1ind = false;
              }

              if(o2==d1 || o2==d2){
                o2ind = true;
              }else{
                o2ind = false;
              }

              if(o1==s1 || o1==s2){
                o1ins = true;
              }else{
                o1ins = false;
              }

              if(o2==s1 || o2==s2){
                o2ins = true;
              }else{
                o2ins = false;
              }

              if(o1ind==true || o2ind==true || o1ins==true || o2ins==true){ 
                if(o1ind==true || o1ins==true){
                  if(o1ind==true){
                    if((int)rbinom(1.0, 0.5)==1){
                      G[i][(l*2)+1] = s1;
                    }else{
                      G[i][(l*2)+1] = s2;
                    }
                  }else{
                    if((int)rbinom(1.0, 0.5)==1){
                      G[i][(l*2)+1] = d1;
                    }else{
                      G[i][(l*2)+1] = d2;
                    }
                  }
                }else{
                  if(o2ind==true){
                    if((int)rbinom(1.0, 0.5)==1){
                      G[i][(l*2)] = s1;
                    }else{
                      G[i][(l*2)] = s2;
                    }
                  }else{
                    if((int)rbinom(1.0, 0.5)==1){
                      G[i][(l*2)] = d1;
                    }else{
                      G[i][(l*2)] = d2;
                    }
                  }
                }
              }else{
                if((int)rbinom(1.0, 0.5)==1){
                  G[i][(l*2)] = d1;
                }else{
                  G[i][(l*2)] = d2;
                }
                if((int)rbinom(1.0, 0.5)==1){
                  G[i][(l*2)+1] = s1;
                }else{
                  G[i][(l*2)+1] = s2;
                }
              }
            }
          }else{            
            if(damP[i]<nind){
              d1 = G[damP[i]][l*2];
              d2 = G[damP[i]][(l*2)+1];     
              if(o1!=d1 &&  o1!=d2 && o2!=d1 && o2!=d2){
                legalP[0]=false;
                if((int)rbinom(1.0, 0.5)==1){
                  G[i][(l*2)+1] = d1;
                }else{
                  G[i][(l*2)+1] = d2;
                }
              }
            }else{
              s1 = G[sireP[i]][l*2];
              s2 = G[sireP[i]][(l*2)+1];     
              if(o1!=s1 &&  o1!=s2 && o2!=s1 && o2!=s2){
                legalP[0]=false;
                if((int)rbinom(1.0, 0.5)==1){
                  G[i][(l*2)+1] = s1;
                }else{
                  G[i][(l*2)+1] = s2;
                }
              }
            } 
            if(G[i][(l*2)]==-999){
              G[i][(l*2)] = rmultinom_size1(A[l], nall[l]);
            }           
          }
        }else{  
          if(no_off[i]==0){      // if founders have no offspring and missing genotypes
            if(G[i][l*2]==-999){
              G[i][l*2] = rmultinom_size1(A[l], nall[l]);
              G[i][(l*2)+1] = rmultinom_size1(A[l], nall[l]);
            }
          }else{                // if founders have offspring and missing genotypes
            if(G[i][l*2]==-999){
              for(itt_set=0; itt_set < nall[l]; itt_set++){ 
                Pall[itt_set] = 0.001;
              }
              index=0;
              for(itt_set=0; itt_set < no_off[i]; itt_set++){ 
                o = par[i][itt_set];
                o1 = G[o][l*2];       // offspring and spouse genotypes //
                o2 = G[o][(l*2)+1];
                s1 = -999;
                s2 = -999;
                if(sireP[o]!=i){
                  if(sireP[o]<nind){
                    s1 = G[sireP[o]][l*2];
                    s2 = G[sireP[o]][(l*2)+1];
                  }
                }else{
                  if(damP[o]<nind){
                    s1 = G[damP[o]][l*2];
                    s2 = G[damP[o]][(l*2)+1];
                  }
                }
                if(o1!=-999){
                  if(index!=2 && o1!=s1 && o1!=s2 && o1!=G[i][(l*2)]){
                     G[i][(l*2)+index]=o1;       // offspring and spouse genotypes //
                     index++;
                  }else{
                     Pall[o1] ++;
                  }
                  if(index!=2 && o2!=s1 && o2!=s2 && o2!=G[i][(l*2)]){
                     G[i][(l*2)+index]=o2;       // offspring and spouse genotypes //
                     index++;
                  }else{
                     Pall[o2] ++;
                  }
                }
              }
              if(index==0){
                G[i][l*2] = rmultinom_size1(Pall, nall[l]);
                if(Pall[G[i][l*2]]>1.0){Pall[G[i][l*2]]--;}
                G[i][(l*2)+1] = rmultinom_size1(Pall, nall[l]);
              }
              if(index==1){
                G[i][(l*2)+1] = rmultinom_size1(Pall, nall[l]);
              }
            }
          }
        }
      }
    }     
    break;

    case 2:
     for(i=0; i<nind; ++i){
      for(l=0; l<nloci; ++l){

        phom = pow(A[l][0], 2.0);
        qhom = pow(A[l][1], 2.0);
        pqhet = 2.0*A[l][0]*A[l][1];
 
        if(damP[i]<nind || sireP[i]<nind){
 
          o1 = G[i][l];
 
          if(damP[i]<nind && sireP[i]<nind){

            d1 = G[damP[i]][l];

            s1 = G[sireP[i]][l];

            TacP=true;

            if(d1==0 || s1==0){
              if(s1==0){s1=d1;}
              if(s1==0){
                if(o1!=0){
                  TacP=false;
                  G[i][l]  = 0;
                }
              }
              if(s1==1){
                if(o1==2 || o1==-999){
                  TacP=false;
                  if((int)rbinom(1.0, qhom/(qhom+pqhet))==1){
                    G[i][l] =0;
                  }else{
                    G[i][l] =1;
                  }
                }
              }
              if(s1==2){
                if(o1!=1){
                  TacP=false;
                  G[i][l]  = 1;
                }
              }
            }else{
              if(d1==2 || s1==2){
                if(s1==2){s1=d1;}
                if(s1==2){
                  if(o1!=2){
                    TacP=false;
                    G[i][l]  = 2;
                  }
                }
                if(s1==1){
                  if(o1==0 || o1==-999){
                    TacP=false;
                    if((int)rbinom(1.0, phom/(phom+pqhet))==1){
                      G[i][l] =2;
                    }else{
                      G[i][l] =1;
                    }
                  }
                }
              }else{
                if(o1==-999){
                  G[i][l] =  int(rbinom(1.0, 0.5));     
                  G[i][l] += int(rbinom(1.0, 0.5));   
                }
              }
            }
          }else{  // only one parent
            if(damP[i]<nind){
              d1=G[damP[i]][l];
            }else{
              d1=G[sireP[i]][l];
            }
            if(d1==0){
              if(o1==2 || o1==-999){
                if((int)rbinom(1.0, phom/(phom+0.5*pqhet))==1){
                  G[i][l]=0;
                }else{
                  G[i][l]=1;
                }
              }
            }
            if(d1==1 && o1==-999){
              if((int)rbinom(1.0, 0.5)==1){
                G[i][l]=1;
              }else{
                G[i][l] = 2*(int)rbinom(1.0, pow(qhom,0.5)/(pow(phom,0.5)+pow(qhom,0.5)));
              }
            }
            if(d1==2){
              if(o1==0 || o1==-999){
                if((int)rbinom(1.0, qhom/(qhom+0.5*pqhet))==1){
                  G[i][l]=2;
                }else{
                  G[i][l]=1;
                }
              }
            }
          }
        }else{   // neither parent
          if(G[i][l]==-999){
            G[i][l] = (int)rbinom(1.0, pow(qhom,0.5));     
            G[i][l] += (int)rbinom(1.0, pow(qhom,0.5));     
          } 
        }       
      }
    }
   
    if(TacP==false){legalP[0]=false;}
    break;
  }

    int records = 0;           

    if(mtype==1 || mtype==3){
      for(l = 0; l < nloci; l++){  
        for(i = 0; i < nind; i++){	                        
          st_GP[records] = G[i][(l*2)];  
          records ++;
          st_GP[records] = G[i][(l*2)+1];  
          records ++;
        }
      }
    }else{
      for(l = 0; l < nloci; l++){  
        for(i = 0; i < nind; i++){	                        
          st_GP[records] = G[i][l];  
          records ++;
        }
      }
    }
        for(i = 0; i < nind; ++i){
           delete [] par[i];
        }
        delete [] Pall;
        delete [] par;
        delete [] no_off;
        delete [] pG;
        delete [] G;
        delete [] pA;
        delete [] A;


 }
 }




















