#include "MCMCped.h" 

extern "C"{  

using namespace std;

void MCMCped(
	int *nindoffP,		 // number of individuals sampled and offspring sampled
        int *ndamP,              // number of candidate dams per offspring
        int *nsireP,             // number of candidate sires per offspring
        int *ntdamP,              // number of candidate dams per offspring
        int *ntsireP,             // number of candidate sires per offspring
        int *nsampP,             // number of samples 
        int *nlociP,		 // number of loci
        int *nallP,              // number of alleles per locus
        int *maxallP,            // number of alleles at most polymorhic locus
        int *maxrepP,            // maximum number of repeat samples per individual
        int *ncatP,              // number of categories over which E1 and E2 varies
        int *nusdP,              // number of categories over which number of unsampled dams vary
        int *nussP,              // number of categories over which number of unsampled sires vary
        int *nparP,              // number of  parameters to estimate
        int *beta_mapP,           // number of linked parameters
        int *nbetaP,
        double *mergeNP,
        int *mergeVP,
        int *mergeUSP,
        int *nmergeP,
 	int *nittP,		    // number of itterations
	int *thinP,		    // thinning interval
	int *burninP,     	    // burn in
        int *idP,                   // numeric id relating samples to individuals
	int *GobsP,                 // observed genotypes	
        int *offidP,                // offspring id
        int *damidP,                // candidate dam id's for each offspring
        int *sireidP,		    // candidate sire id's for each offspring
        double *X_design_betaDusP,  // design matrices for dam variables for each offspring	
        double *X_design_betaSusP,  // design matrices for sire variables for each offspring	
        double *X_design_betaDSusP, // design matrices for dam:sire interactions/relational variables 
        double *X_design_betaDsP,  // design matrices for dam variables for each offspring	
        double *X_design_betaSsP,  // design matrices for sire variables for each offspring	
        double *X_design_betaDSsP, // design matrices for dam:sire interactions/relational variables 
	double *st_AP,	         // starting allele frequencies
	double *st_E1P,	         // starting values of E1 and E2
	double *st_E2P,	         // starting values of E1 and E2
	double *st_betaP,        // starting vector of parameters (beta)
        double *st_usP,
        int *st_GP,              // starting true genotypes    
        int *st_damP,            // starting vector of dams
        int *st_sireP,           // starting vector of sires
	double *post_AP,	 // posterior distribution of allele frequencies
	double *post_E1P,	 // posterior distribution of E1
	double *post_E2P,	 // posterior distribution of E2
	double *post_betaP,	 // posterior distribution of beta
	double *post_usP,	 // posterior distribution of beta
        int *post_GP,            // posterior distribution of true genotypes
        int *post_PP,            // posterior distribution of pedigrees
	double *prior_E1P,       // prior distribution of E1 beta(a,b)
	double *prior_E2P,       // prior distribution of E2 beta(a,b)
	double *prior_beta_muP,  // prior distribution of beta mean(MVN) = mu, var(MVN) = sigma
	double *prior_beta_invsigmaP,
	double *log_detP,
	double *prior_us_muP,  // prior distribution of USdam mean(lnorm) = mu, var(lnorm) = sigma
	double *prior_us_sigmaP,
        double *int_E1P,         // standard deviation of normal candidate generating function for E1
        double *int_E2P,         // standard deviation of normal candidate generating function for E2
        double *int_betaP,       // standard deviation of normal candidate generating function for beta
        double *int_usP,       // standard deviation of normal candidate generating function for beta
        int *categoriesP,        // categoreies over which E varies 
        int *usdamcatP,
        int *ussirecatP,
        int *estimatingP,            // logicals inidcating whether parameters should be estimated or fixed
        int *store_postP){

int 	nind = nindoffP[0],
        noff = nindoffP[1],
        nsamp = nsampP[0],        	
	nloci = nlociP[0],
        nbeta = nbetaP[0],
        maxall = maxallP[0],
        maxrep = maxrepP[0],
        ncat = ncatP[0],
        nusd = nusdP[0],
        nuss = nussP[0],
        *npar = nparP,
        nitt = nittP[0], 
        thin = thinP[0],
        burnin = burninP[0],
        *nall = nallP,
        *id = idP,
        *categories = categoriesP,
        *usdamcat = usdamcatP,
        *ussirecat = ussirecatP,
        acceptB = 1000,
        acceptUS = 1000,
        acceptE1 = 1000,
        acceptE2 = 1000,
        tall=0,
        i,
        l,
        cnt,
        d,
        s;

bool    estP = bool(estimatingP[0]),           
        estG = bool(estimatingP[1]), 
        estA = bool(estimatingP[2]), 
        estE1 = bool(estimatingP[3]), 
        estE2 = bool(estimatingP[4]), 
        estbeta = bool(estimatingP[5]),      
        estUSdam = bool(estimatingP[6]), 
        estUSsire = bool(estimatingP[7]), 
        estUS = (estUSdam || estUSsire),
        perlocus = bool(estimatingP[8]),
        USdamsire = bool(estimatingP[10]),
        checkP = bool(estimatingP[11]),
        jointP = bool(estimatingP[12]),

        est_pE1 = FALSE,           
        est_pE2 = FALSE,
        est_pbeta = FALSE,           
        est_pus = FALSE,
        writeG = bool(store_postP[0]),   
        writeA = bool(store_postP[1]),     
        writeJP = bool(store_postP[2]),  
        verbose = bool(store_postP[3]);

int mtype = estimatingP[9],
    DSapprox = estimatingP[13];


/*************************************
* declare some variable sized arrays *
*************************************/

int     *pG,
        **G,
        *pGobs,
        **Gobs,       
        *ppost_G,
        **post_G;
double  *pA,
        **A,
        *pE_mat,
        **E_mat,
        *pLE_mat,
        **LE_mat;

  pG = new(nothrow) int [(1+int(mtype==1 || mtype==3))*nind*nloci];
  if(pG==NULL){
    error("NO MEMORY for G\n");
  }

  G = new(nothrow) int* [nind];
  if(G==NULL){
   error("NO MEMORY for G\n");
  }

  pGobs = new(nothrow) int [(1+int(mtype==1 || mtype==3))*nsamp*nloci];
  if(pGobs==NULL){
    error("NO MEMORY for Gobs\n");
  }

  Gobs = new(nothrow) int* [nsamp];
  if(Gobs==NULL){
   error("NO MEMORY for Gobs\n");
  }

  pA = new(nothrow) double [nloci*maxall];
  if(pA==NULL){
   error("NO MEMORY for A\n");
  }

  A = new(nothrow) double* [nloci];
  if(A==NULL){
   error("NO MEMORY for A\n");
  }

  pE_mat = new(nothrow) double [ncat*(4+3*int(mtype==1 || mtype==3)+2*int(mtype==2))*nloci];
  if(pE_mat==NULL){
   error("NO MEMORY for E_mat\n");
  }

  E_mat = new(nothrow) double* [nloci];
  if(E_mat==NULL){
   error("NO MEMORY for E_mat\n");
  }

  pLE_mat = new(nothrow) double [ncat*(4+3*int(mtype==1 || mtype==3)+2*int(mtype==2))*nloci];
  if(pLE_mat==NULL){
   error("NO MEMORY for LE_mat\n");
  }

  LE_mat = new(nothrow) double* [nloci];
  if(LE_mat==NULL){
   error("NO MEMORY for LE_mat\n");
  }

  cnt = 0;
  for (i=0; i<nind; ++i){
    G[i] = &pG[cnt];
    cnt  += ((1+int(mtype==1 || mtype==3))*nloci);
  }

  cnt = 0;
  for (i=0; i<nsamp; ++i){
    Gobs[i] = &pGobs[cnt];
    cnt  += ((1+int(mtype==1 || mtype==3))*nloci);
  }

  cnt = 0;
  for (i=0; i<nloci; ++i){
    A[i] = &pA[cnt];
    cnt  += maxall;
  }

  cnt = 0;
  for (i=0; i<nloci; ++i){
    E_mat[i] = &pE_mat[cnt];
    cnt  += (ncat*(4+3*int(mtype==1 || mtype==3)+2*int(mtype==2)));
  }

  cnt = 0;
  for (i=0; i<nloci; ++i){
    LE_mat[i] = &pLE_mat[cnt];
    cnt  += (ncat*(4+3*int(mtype==1 || mtype==3)+2*int(mtype==2)));
  }

  if(writeG==true){

    tall = 0;

    for(i = 0; i < nloci; i++){
      tall += int(0.5*nallP[i]*(nallP[i]+1.0));
    }

    ppost_G = new(nothrow) int [nind*tall];
    if(ppost_G==NULL){
      error("NO MEMORY for posterior_G\n");
    }

    post_G = new(nothrow) int* [nind];
    if(post_G==NULL){
      error("NO MEMORY for posterior_G\n");
    }

    cnt = 0;
    for (i=0; i<nind; ++i){
      post_G[i] = &ppost_G[cnt];
      cnt  += tall;
    }

    for (i=0; i<nind; ++i){
      for (l=0; l<tall; ++l){
        post_G[i][l] = 0;
      }
    }
  }

  Matrix<int> *post_P = new Matrix<int>[noff]; 

  if(writeJP==false){

    for (i=0; i<noff; ++i){
      post_P[i] = Matrix<int>(ndamP[i]*nsireP[i], 1); 
    }
  }


/*************************************
* declare some variable sized arrays *
*************************************/
          
           int tot_par = npar[0]+npar[1]+npar[2]+npar[3]+npar[4]+npar[5],
               nus = nusd+nuss,
               DSuu[2];  // when DS exists and is formed by mate relational npar[4] may be greater than 0
                         // but ambiguity exists as to whether both sexes are unsampled, or one, and if one
                         // - which one  DSuu[0] is 1 if missing Dam data exists, DSuu[1] is 1 if missing Sire data
                         // exists hack hack hack
               if(npar[0]>0 || (npar[4]>0 && nusd>0)){
                 DSuu[0] =1;
               }else{
                 DSuu[0] =0;
               }
               if(npar[2]>0 || (npar[4]>0 && nuss>0)){
                 DSuu[1] =1;
               }else{
                 DSuu[1] =0;
               }

               int ncatnloci = ncat*nloci*int(perlocus) + ncat*(1-int(perlocus));

Matrix<double> E1_0 (ncatnloci,1,st_E1P), 	        // starting vector of allelic dropout rate (new)
               E1_1 (ncatnloci,1,st_E1P),	        // starting vector of allelic dropout rate (old)
	       E2_0 (ncatnloci,1,st_E2P), 	        // starting vector of stochastic error rate (new)
  	       E2_1 (ncatnloci,1,st_E2P),	        // starting vector of stochastic error rate (old)
               us_0 (nus,1,st_usP),             // starting vector of unsampled population size (new)
               us_1 (nus,1,st_usP),             // starting vector of unsampled population size (old)
               beta_0 (nbeta,1,st_betaP), 	// starting vector of beta (new)
               beta_1 (nbeta,1,st_betaP), 	// starting vector of beta (old)
               beta_mapped (tot_par, 1),
               int_E1,	        
               int_E2,	        
               int_beta,	        
               int_us,	        
               prior_E1,
               prior_E2,
               prior_beta_mu,
               prior_beta_invsigma,
               prior_us_mu,
               prior_us_sigma;

               Matrix<double> *X_design_G = new Matrix<double>[noff];
               Matrix<double> *X_design_GD = new Matrix<double>[noff];
               Matrix<double> *X_design_GS = new Matrix<double>[noff];
               Matrix<double> *X_design_betaDus = new Matrix<double>[noff];
               Matrix<double> *X_design_betaSus = new Matrix<double>[noff];
               Matrix<double> *X_design_betaDSus = new Matrix<double>[noff];
               Matrix<double> *X_design_betaDs = new Matrix<double>[noff];
               Matrix<double> *X_design_betaSs = new Matrix<double>[noff];
               Matrix<double> *X_design_betaDSs = new Matrix<double>[noff];
               Matrix<double> *ratio_0 = new Matrix<double>[2];
               Matrix<double> *ratio_1 = new Matrix<double>[2];

        double llB_0 = 0.0,
               llB_1 = 0.0,
               llE_0 = 0.0,
               llE_1 = 0.0,
               llUS_0 = 0.0,
               llUS_1 = 0.0,
               log_det = 0.0;  

/************************/    
/* MH tuning parameters */
/************************/

               if(estE1){    
                 int_E1 = Matrix<double>(ncatnloci,ncatnloci, int_E1P);
               }
               if(estE2){    
                 int_E2 = Matrix<double>(ncatnloci,ncatnloci, int_E2P);
               }
               if(estbeta){    
                 int_beta = Matrix<double>(nbeta,nbeta, int_betaP);
               }
               if(estUS){    
                 int_us = Matrix<double>(nus,nus, int_usP);
               }

/***********************/    
/* Prior specification */
/***********************/

	       if(int(prior_E1P[0])!=999){       
          	 prior_E1 = Matrix<double> (ncatnloci,2, prior_E1P);
                 est_pE1 = TRUE;
               }
               if(int(prior_E2P[0])!=999){       
          	 prior_E2 = Matrix<double>(ncatnloci,2, prior_E2P);
                 est_pE2 = TRUE;
               }
               if(int(prior_beta_muP[0])!=999){      
          	 prior_beta_mu = Matrix<double> (nbeta,1, prior_beta_muP);
          	 prior_beta_invsigma = Matrix<double> (nbeta,nbeta, prior_beta_invsigmaP);
                 log_det = log_detP[0];
                 est_pbeta = TRUE;
               }

               if(int(prior_us_muP[0])!=999){         
          	 prior_us_mu = Matrix<double>  (nus,1, prior_us_muP);
          	 prior_us_sigma = Matrix<double> (nus,1, prior_us_sigmaP);
                 est_pus = TRUE;
               }

/*************************/    
/* Merging specification */
/*************************/
               int nmerge = nmergeP[0];
               int *mergeV = mergeVP;
               int *mergeUS = mergeUSP;
               Matrix<double> *mergeN = new Matrix<double>[nmerge];	       

               if(nmerge>0){    
                 for(i=0; i<nmerge; i++){   
                   cnt=0;  
                   mergeN [i] = Matrix<double>(2,noff);           
                   for(l=0; l<(2*noff); l++){
                     mergeN [i][l] = mergeNP[cnt];
                     cnt++;   
                   }
                 }
               }
               if(noff>0){
                 ratio_0[0] = ones<double>(noff,1);
                 ratio_0[1] = ones<double>(noff,1);
                 ratio_1[0] = ones<double>(noff,1);
                 ratio_1[1] = ones<double>(noff,1);
               }

	map<int, int> *Dams = new map<int, int>[noff];     // two way indexing vectors
	map<int, int> *Sires = new map<int, int> [noff];    // map[i][dam_id] = n   dam_id is the n^th mother of the i^th individual
        int *dam = st_damP;            // starting vector of dams
	int *sire = st_sireP;		
        Matrix<int> *Dams_vec = new Matrix<int>[noff];   // Matrix[i][n] = dam_id    the n^th mother of the i^th individul is dam.id 
	Matrix<int> *Sires_vec = new Matrix<int>[noff];
	

	GetRNGstate();                                 // get seed for random number generation
 
        if(nloci!=0){  
          read_G(st_GP, nind, nloci, G, mtype);
          read_A(st_AP, nloci, A, nall);
        }

        if(nsamp!=0 && (estG==TRUE || estE1==TRUE || estE2==TRUE || estA==TRUE)){
          read_G(GobsP, nsamp, nloci, Gobs, mtype);
        }   
 
        read_stP(noff, ndamP, damidP, nsireP, sireidP,Dams,Sires,Dams_vec,Sires_vec);
        		        
        if(nbeta>0){
           read_X_beta(noff,ntdamP,ntsireP,npar,X_design_betaDusP, X_design_betaSusP,X_design_betaDSusP, X_design_betaDsP, X_design_betaSsP,X_design_betaDSsP, X_design_betaDus,X_design_betaSus,X_design_betaDSus, X_design_betaDs, X_design_betaSs,X_design_betaDSs);

           for(i = 0; i < tot_par; i++){
             beta_mapped[i] = beta_0[beta_mapP[i]];
           }

           
           llB_0 = LLP_B(offidP,noff,nind,X_design_betaDus,X_design_betaSus,X_design_betaDSus,X_design_betaDs,
X_design_betaSs,X_design_betaDSs,npar, DSuu, dam,sire,beta_mapped,ntdamP,ntsireP,ndamP, nsireP, Dams,Sires, nusd,  usdamcat, nuss, ussirecat, us_0, ratio_0, nmerge, mergeV, mergeUS, mergeN, DSapprox);

           if(est_pbeta){
               llB_0 += lmvnormM(beta_0,  nbeta, prior_beta_mu, log_det, prior_beta_invsigma);
             }
           llB_1 = llB_0;
           ratio_1[0] = ratio_0[0];
           ratio_1[1] = ratio_0[1];
         }
             
        if(estUS==TRUE){
          llUS_0 = LLN_P(offidP, noff, nind, ntdamP, ntsireP, dam, sire, nusd, usdamcat, nuss, ussirecat, us_0, ratio_1);
          if(est_pus){
            for(i = 0; i < (nusd+nuss*(1-int(USdamsire))); i++){
              llUS_0 += dlnorm(us_0[i], prior_us_mu[i], prior_us_sigma[i],1);
            }
          }
          llUS_1 = llUS_0;
        }

        if(estG==TRUE){      
          Error_Mat(E1_0, E2_0, E_mat, ncat, nall, nloci, false, perlocus, mtype);
          if(jointP){
            calcX_G(X_design_G, offidP, noff , ndamP, nsireP, nind, Dams_vec, Sires_vec, G, nloci, A, mtype);
          }else{
            calcX_GD(X_design_GD, offidP, noff , ndamP, nind, Dams_vec, sire, G, nloci, A, mtype);
            calcX_GS(X_design_GS, offidP, noff , nsireP, nind, Sires_vec, dam, G, nloci, A, mtype);
          }
        }else{
          if(estP==TRUE){
           Error_Mat(E1_0, E2_0, E_mat, ncat, nall, nloci, false, perlocus, mtype);
           calcX_Gcervus(X_design_G, offidP, noff , ndamP, nsireP, nind, Dams_vec, Sires_vec, G, nloci, A,E_mat, mtype, nall);         
           if(jointP==FALSE){
               calcX_GcervusS(X_design_GS, X_design_G, offidP, noff, ndamP, nsireP, Dams, dam);
               calcX_GcervusD(X_design_GD, X_design_G, offidP, noff, ndamP, nsireP, Sires, sire);
           }
         }
       }

        if(estE1==TRUE || estE2==TRUE){    
            Error_Mat(E1_0, E2_0, LE_mat, ncat, nall, nloci, true, perlocus,mtype); 
            llE_0 = LLE_G(Gobs, G, nloci,id,nsamp,categories,LE_mat, mtype, A);
          if(est_pE1){		
            for(i = 0; i < ncatnloci; i++){  
              llE_0 += dbeta(E1_0[i], prior_E1P[i], prior_E1P[i+ncatnloci],1);
            }
          }
          if(est_pE2){		
            for(i = 0; i < ncatnloci; i++){  
              llE_0 += dbeta(E2_0[i], prior_E2P[i], prior_E2P[i+ncatnloci],1);
            }
          }
          llE_1 = llE_0;
        }
      
	int itt;
	int write_postE = 0;
        int write_postA = 0; 
        int write_postB = 0; 
        int write_postG = 0; 
        int write_postP = 0;
        int write_postUS = 0;
        double m_ll;
        int tall_tmp = 0;
        int nl;
        int a1;
        int a2;
        int na1;
        int na2;
          
        if(verbose==TRUE){
          Rprintf("\n Starting parameterisation\n");
          if(estbeta==TRUE){
	   Rprintf("\n                      beta = \n");      
	    for (i=0; i<nbeta; ++i){
	     Rprintf("                            %10.5f\n", beta_0[i]);
            }  
          }
          if(estUS==TRUE){
	    Rprintf("\n      unsampled population = \n");
	    for (i=0; i<nus; ++i){
	      Rprintf("                            %10.5f\n", us_0[i]);
            }
          }
          if(estE1==TRUE){
	   Rprintf("\n                        E1 = \n");
	    for (i=0; i<ncatnloci; ++i){
	      Rprintf("                            %10.5f\n", E1_0[i]);
            }  
          }
         if(estE2==TRUE){
           Rprintf("\n                        E2 = \n");
	    for (i=0; i<ncatnloci; ++i){
	      Rprintf("                            %10.5f\n", E2_0[i]);
            }
          }
        }

//***************************************************************************************************************
//************************************ MCMC MCMC MCMC MCMC MCMC MCMC MCMC ***************************************
//***************************************************************************************************************

 	for(itt = 0; itt < nitt; itt++){       // start the fucker running     

/**************
* sample P    *
***************/
	    if(estP==TRUE){ 
              for(i = 0; i < tot_par; i++){
                beta_mapped[i] = beta_1[beta_mapP[i]];
              }
              if(jointP){
                sampP(offidP,noff,X_design_G,npar, DSuu, X_design_betaDus, X_design_betaSus,X_design_betaDSus, X_design_betaDs, X_design_betaSs,X_design_betaDSs, dam,sire,beta_mapped,us_1, usdamcat,ussirecat, nusd, nuss, ndamP,nsireP,ntdamP,ntsireP, Dams_vec,Sires_vec, nind, nmerge, mergeV, mergeUS, mergeN, checkP);
              }else{
                sampS(offidP,noff,X_design_GS,npar, DSuu, X_design_betaSus, X_design_betaDSus, X_design_betaSs, X_design_betaDSs, dam,sire,beta_mapped,us_1, ussirecat, nusd, nuss, ndamP,nsireP,ntdamP,ntsireP, Sires_vec, nind, nmerge, mergeV, mergeUS, mergeN, Dams, Sires, checkP);
                if(estG){
                   calcX_GD(X_design_GD, offidP, noff , ndamP, nind, Dams_vec, sire, G, nloci, A, mtype);
                }else{
                   calcX_GcervusD(X_design_GD, X_design_G, offidP, noff, ndamP, nsireP, Sires, sire);
                }
                sampD(offidP,noff,X_design_GD,npar, DSuu, X_design_betaDus, X_design_betaDSus, X_design_betaDs, X_design_betaDSs, dam,sire,beta_mapped,us_1, usdamcat, nusd, ndamP,nsireP,ntdamP,ntsireP, Dams_vec, nind, nmerge, mergeV, mergeUS, mergeN, Dams, Sires, checkP);
                if(estG==FALSE){
                  calcX_GcervusS(X_design_GS, X_design_G, offidP, noff, ndamP, nsireP, Dams, dam);                
                }
              }
              if(estbeta==TRUE){
                llB_1 = LLP_B(offidP,noff,nind,X_design_betaDus,X_design_betaSus,X_design_betaDSus,X_design_betaDs,
X_design_betaSs,X_design_betaDSs,npar, DSuu, dam,sire,beta_mapped,ntdamP,ntsireP,ndamP,nsireP,Dams,Sires, nusd,  usdamcat, nuss, ussirecat, us_1, ratio_1, nmerge, mergeV, mergeUS, mergeN, DSapprox);
                if(est_pbeta){
                  llB_1 += lmvnormM(beta_1,  nbeta, prior_beta_mu, log_det, prior_beta_invsigma);
                }
               }
              if(estUS==TRUE){
                llUS_1 = LLN_P(offidP, noff, nind, ntdamP, ntsireP, dam, sire, nusd, usdamcat, nuss, ussirecat, us_1, ratio_1);
                if(est_pus){
                  for(i = 0; i < (nusd+nuss*(1-int(USdamsire))); i++){
                    llUS_1 += dlnorm(us_1[i], prior_us_mu[i], prior_us_sigma[i],1);
                  }
                }
              }
            }

/**************
* sample G & A*
***************/


	    if(estG==TRUE){  

              switch(mtype){ 
                 case 1:  
                 sampGC(nsamp,Gobs,G,nall,nloci,id,A,categories,E_mat,maxall,maxrep,dam,sire, nind, estA);
                 break;
                 case 2:
                 sampDomG(nsamp,Gobs,G,nall,nloci,id,A,categories,E_mat,maxall,maxrep,dam,sire, nind, estA);
                 break;
                 case 3:
                 sampG(nsamp,Gobs,G,nall,nloci,id,A,categories,E_mat,maxall,maxrep,dam,sire, nind, estA);
                 break;
              }

              if(estP==TRUE){
                if(jointP){
                  calcX_G(X_design_G, offidP, noff , ndamP, nsireP, nind, Dams_vec, Sires_vec, G, nloci,A, mtype);
                }else{
                  calcX_GD(X_design_GD, offidP, noff , ndamP, nind, Dams_vec, sire, G, nloci, A, mtype);
                  calcX_GS(X_design_GS, offidP, noff , nsireP, nind, Sires_vec, dam, G, nloci, A, mtype);
                }
              }

              if(estE1==TRUE || estE2==TRUE){
                llE_1 = LLE_G(Gobs, G, nloci,id,nsamp,categories,LE_mat, mtype, A);
                if(est_pE1){		
                  for(i = 0; i < ncatnloci; i++){  
                   llE_1 += dbeta(E1_1[i], prior_E1P[i], prior_E1P[i+ncatnloci],1);
                  }
                }
                if(est_pE2){		
                  for(i = 0; i < ncatnloci; i++){  
                   llE_1 += dbeta(E2_1[i], prior_E2P[i], prior_E2P[i+ncatnloci],1);
                  }
                }    
              }
            }

/***********
* sample E1 *
************/

              if(estE1==TRUE){ 

                E1_0 = fabs(rmvnormM(E1_1, int_E1, ncatnloci));  

                Error_Mat(E1_0, E2_0, LE_mat, ncat, nall, nloci, true, perlocus, mtype);

                llE_0 = LLE_G(Gobs, G, nloci,id,nsamp,categories,LE_mat, mtype,A);

                if(est_pE1){		
                  for(i = 0; i < ncatnloci; i++){  
                    llE_0 += dbeta(E1_0[i], prior_E1P[i], prior_E1P[i+ncatnloci],1);
                  }
                }

                if(est_pE2){		
                  for(i = 0; i < ncatnloci; i++){  
                    llE_0 += dbeta(E2_0[i], prior_E2P[i], prior_E2P[i+ncatnloci],1);
                  }
                }

	        m_ll = std::min(1.0, llE_0-llE_1); 
               
	        if(m_ll<log(runif(0.0,1.0))){       
                  llE_0 = llE_1;  
                  E1_0 = E1_1; 
                  Error_Mat(E1_0, E2_0, LE_mat, ncat, nall, nloci, true, perlocus, mtype);
                  acceptE1 --;
                }else{                               
                  Error_Mat(E1_0, E2_0, E_mat, ncat, nall, nloci, false, perlocus, mtype);
                }
                llE_1 = llE_0;
                E1_1 = E1_0;                  
              }

/***********
* sample E2 *
************/
       
              if(estE2==TRUE){ 
                E2_0 = fabs(rmvnormM(E2_1, int_E2, ncatnloci));  
   
                Error_Mat(E1_0, E2_0, LE_mat, ncat, nall, nloci, true, perlocus, mtype);

                llE_0 = LLE_G(Gobs, G, nloci,id,nsamp,categories,LE_mat, mtype,A);
 
                if(est_pE1){		
                  for(i = 0; i < ncatnloci; i++){  
                    llE_0 += dbeta(E1_0[i], prior_E1P[i], prior_E1P[i+ncatnloci],1);
                  }
                }

                if(est_pE2){		
                  for(i = 0; i < ncatnloci; i++){  
                    llE_0 += dbeta(E2_0[i], prior_E2P[i], prior_E2P[i+ncatnloci],1);
                  }
                }

	        m_ll = std::min(1.0, llE_0-llE_1); 
               
	        if(m_ll<log(runif(0.0,1.0))){    
	          llE_0 = llE_1;  
                  E2_0 = E2_1; 
                  Error_Mat(E1_0, E2_0, LE_mat, ncat, nall, nloci, true, perlocus, mtype);
                  acceptE2 --;
                }else{                             
                  Error_Mat(E1_0, E2_0, E_mat, ncat, nall, nloci, false, perlocus, mtype);
                }
                llE_1 = llE_0;
                E2_1 = E2_0;
              }
 
/**************
* sample beta *
***************/

            if(estbeta==TRUE){ 

              beta_0 = rmvnormM(beta_1, int_beta, nbeta);

              for(i = 0; i < tot_par; i++){
                beta_mapped[i] = beta_0[beta_mapP[i]];
              }

              llB_0 = LLP_B(offidP,noff,nind,X_design_betaDus,X_design_betaSus,X_design_betaDSus,X_design_betaDs,
X_design_betaSs,X_design_betaDSs,npar,DSuu,dam,sire,beta_mapped,ntdamP,ntsireP,ndamP,nsireP,Dams,Sires, nusd,  usdamcat, nuss, ussirecat, us_0, ratio_0, nmerge, mergeV, mergeUS, mergeN, DSapprox);

              if(est_pbeta){
                llB_0 += lmvnormM(beta_0,  nbeta, prior_beta_mu, log_det, prior_beta_invsigma);
              }

 	      m_ll = std::min(1.0, llB_0-llB_1); 
	      if(m_ll<log(runif(0.0,1.0))){
               llB_0 = llB_1;
	       beta_0 = beta_1;		
               ratio_0[0] = ratio_1[0];
               ratio_0[1] = ratio_1[1];
               acceptB --;
              }                
              llB_1 = llB_0;
              beta_1 = beta_0;	
              ratio_1[0] = ratio_0[0];
              ratio_1[1] = ratio_0[1];
              if(estUS==TRUE){
                llUS_1 = LLN_P(offidP, noff, nind, ntdamP, ntsireP, dam, sire, nusd, usdamcat, nuss, ussirecat, us_1, ratio_1);
                if(est_pus){
                  for(i = 0; i < (nusd+nuss*(1-int(USdamsire))); i++){
                     llUS_1 += dlnorm(us_1[i], prior_us_mu[i], prior_us_sigma[i],1);
                  }
                }
              }
            }

/*************************
* sample population size *
**************************/

           if(estUS==TRUE){ 
             us_0 = fabs(rmvnormM(us_1, int_us, nus));
             if(estUSdam==FALSE){
                for(i = 0; i < nusd; i++){
                   us_0[i]= us_1[i];
                }  
             }
             if(estUSsire==FALSE){
                for(i = 0; i < nuss; i++){
                   us_0[nusd+i]= us_1[nusd+i];
                }
             }
             if(USdamsire==TRUE){
                for(i = 0; i < nusd; i++){
                   us_0[nusd+i]= us_0[i];
                }
             }

             llUS_0 = LLN_P(offidP, noff, nind, ntdamP, ntsireP, dam, sire, nusd, usdamcat, nuss, ussirecat, us_0, ratio_1);
             
             if(est_pus){
               for(i = 0; i < (nusd+nuss*(1-int(USdamsire))); i++){
                  llUS_0 += dlnorm(us_0[i], prior_us_mu[i], prior_us_sigma[i],1);
                }
             }

	     m_ll = std::min(1.0, llUS_0-llUS_1); 

              if(m_ll<log(runif(0.0,1.0))){
                llUS_0 = llUS_1;
	        us_0 = us_1;	
                acceptUS--;
              }
              llUS_1 = llUS_0;
              us_1 = us_0;
            }

/************************** 
* write posterior samples *
***************************/
  
	   if(itt>=burnin && itt%thin == 0){

             for(i = 0; i < ncatnloci; i++){
               if(estE1==TRUE){	
	        post_E1P[write_postE]= E1_0[i];
               }
               if(estE2==TRUE){	
                post_E2P[write_postE]= E2_0[i];
               }
               write_postE ++;
	     }
             
             if(writeA==TRUE && estA==TRUE){ 
                 for(l=0; l<nloci; l++){  
                    for(i=0; i<nall[l]; i++){
                     post_AP[write_postA] = A[l][i];                 // record    
                     write_postA ++;
                    } 
                 }
              }
              if(estbeta==TRUE){ 
                 for (i=0; i<nbeta; ++i){
                  post_betaP[write_postB] = beta_0[i];                 // record   
                  write_postB ++; 
                 }
              }              
              if(estUSdam){
                for(i = 0; i < nusd; i++){
                   post_usP[write_postUS]= us_0[i];
                   write_postUS ++; 
                }  
              }
              if(estUSsire){
                for(i = 0; i < nuss; i++){
                  post_usP[write_postUS]= us_0[nusd+i];
                  write_postUS ++; 
                }
              }
              if(writeG==true){
                if(mtype==1 || mtype==3){
                  tall_tmp = 0;
                  for(l=0; l<nloci; l++){ 
                    nl = nall[l];
                    for(i=0; i<nind; i++){
                     na1 = G[i][(l*2)];
                     na2 = G[i][(l*2)+1];
                     if(na1<na2){
                       a1=na1;
                       a2=na2;
                     }else{
                       a1=na2;
                       a2=na1;
                     }  
                     post_G[i][int(a1*(nl-(0.5*(a1+1)))+tall_tmp+a2)] ++;
                    }
                  tall_tmp += int(0.5*nl*(nl+1.0));  
                  }
                }else{
                  for(l=0; l<nloci; l++){ 
                    for(i=0; i<nind; i++){
                      na1 = G[i][l];
                      post_G[i][(l*3)+na1] ++;
                    }
                  }
                }
              }
            if(estP==true){      
              if(writeJP==true){
                 for(i=0; i<noff; i++){
                    post_PP[(write_postP*noff*2)+i] = dam[offidP[i]];
                    post_PP[(write_postP*noff*2)+i+noff] = sire[offidP[i]];
                 }
                 write_postP ++;
              }else{
                 for(i=0; i<noff; i++){
                   d = dam[offidP[i]];
                   d = Dams[i][d];
                   s = sire[offidP[i]];
                   s = Sires[i][s];
                   write_postP = (nsireP[i]*d)+s;
                   post_P[i][write_postP] ++;
                 }
              }
           }
        }

// some summary information: don't panic!

        if(itt%1000==0 && itt!=0 && verbose==TRUE){  
          Rprintf("\n\nMCMC iteration %d of %d \n", itt+1, nitt);
          if(estbeta==TRUE){
	   Rprintf("\n                      beta = \n");      
	    for (i=0; i<nbeta; ++i){
	     Rprintf("                            %10.5f\n", beta_0[i]);
            }  
             Rprintf("\nMetropolis acceptance rate =%10.5f\n", long(acceptB)/1000.0);
          }
          if(estUS==TRUE){
	    Rprintf("\n      unsampled population = \n");
	    for (i=0; i<nus; ++i){
	      Rprintf("                            %10.5f\n", us_0[i]);
            }
            Rprintf("\nMetropolis acceptance rate =%10.5f\n", long(acceptUS)/1000.0);
          }
          if(estE1==TRUE){
	   Rprintf("\n                        E1 = \n");
	    for (i=0; i<ncatnloci; ++i){
	      Rprintf("                            %10.5f\n", E1_0[i]);
            }  
            Rprintf("\nMetropolis acceptance rate =%10.5f\n", long(acceptE1)/1000.0);
          }
         if(estE2==TRUE){
           Rprintf("\n                        E2 = \n");
	    for (i=0; i<ncatnloci; ++i){
	      Rprintf("                            %10.5f\n", E2_0[i]);
            }
           Rprintf("\nMetropolis acceptance rate =%10.5f\n", long(acceptE2)/1000.0);
          }
          acceptB=1000;
          acceptE1=1000;
          acceptE2=1000;
          acceptUS=1000;
          R_FlushConsole();
          R_ProcessEvents();
        }
     }
	
/************************************** END MCMC END MCMC END MCMC ********************************************/

        if(writeG==true){
	 write_postG = 0;
            for(l = 0; l < tall ; l++){
               for(i = 0; i < nind ; i++){	
	         post_GP[write_postG] = post_G[i][l];
                 write_postG++;
	       }
	   }
        }
     if(estP==TRUE){
        if(writeJP==false){
	 write_postP = 0;
               for(i = 0; i < noff; i++){	
                  for(l = 0; l < ndamP[i]*nsireP[i]; l++){   
	             post_PP[write_postP] = post_P[i][l];
                      write_postP++;
	       }
	   }
        }
     }
    
     
// read back posterior distribution of genotypes and allele frequencies to R

	PutRNGstate();


delete [] pG;
delete [] G;
delete [] pGobs;
delete [] Gobs;
delete [] pA;
delete [] A;
delete [] pE_mat;
delete [] E_mat;
delete [] pLE_mat;
delete [] LE_mat;
if(writeG==true){
delete [] ppost_G;
delete [] post_G;
}

        delete [] mergeN;
        delete [] post_P;
        delete [] X_design_G;
        delete [] X_design_GD;
        delete [] X_design_GS;
        delete [] X_design_betaDus;
        delete [] X_design_betaSus;
        delete [] X_design_betaDSus;
        delete [] X_design_betaDs;
        delete [] X_design_betaSs;
        delete [] X_design_betaDSs;
        delete [] Dams;
        delete [] Sires;
        delete [] Dams_vec;
        delete [] Sires_vec;
        delete [] ratio_0;
        delete [] ratio_1;
}
} // extern "C" 


 
