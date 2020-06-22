#include "SpecRand.h"

int rmultinom_size1(double *prob, int mult_size){

  double p_tot = 0.0;
  double pp;
  int k;
  int nk = 0;

  for(k= 0; k < mult_size; k++){
    p_tot += prob[k];
  }

  for(k= 0; k < mult_size; k++){	
    pp = prob[k]/p_tot;
    if(pp>=1.0){pp = 1.0;}
    nk = (int)rbinom(1.0, pp);
    if(nk==1){
      break;
    }
    p_tot -= prob[k];
  }	
  return k;	
}


int rmultinom_size1M(Matrix<double> prob, int mult_size){

double p_tot = 0.0;
double pp;
int k;
int nk = 0;

for(k= 0; k < mult_size; k++){
p_tot += prob[k];}

	for(k= 0; k < mult_size; k++){	
	  pp = prob[k]/p_tot;
			if(pp>=1.0){pp = 1.0;}
		  	  nk = (int)rbinom(1.0, pp);
				if(nk==1){
			  	  break;
				}
		  	  p_tot -= prob[k];
	}	
        return k;	
}


// GENERATE RANDOM DIRICHLET
// *************************

void rdirichlet(double *count, int K, double *freq){

double sumgamma = 0.0;
int i;

for(i=0; i<K; i++){
freq[i] = rgamma(count[i], 1.0);
sumgamma += freq[i];
}

for(i=0; i<K; i++){
freq[i]  = freq[i]/sumgamma;
}

}

// GENERATE RANDOM MULTIVARIATE NORMAL DEVIATES
// ********************************************

Matrix<double> rmvnormM(Matrix<double> beta, Matrix<double> cholD, int K){

Matrix<double> X (1,K); 
int i;	

for(i=0; i<K; i++){
X[i] = rnorm(0.0,1.0);
}

return beta+t(X*cholD);
}


