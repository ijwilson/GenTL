// Time-stamp: <2006-04-12 23:16:30 uijw>
#ifndef SUBDIVPAC_H__
#define SUBDIVPAC_H__
#include "pac.h"    
class subdivpac  : public pac {
public:
  /** the constructor                                                  */
  subdivpac(TNT::Array2D<int> &hps, const std::vector<double> &dis 	    
	    ,const TNT::Array1D<int> &loc, int nprm=10	    
	    ,const unsigned int rsd=1, bool plm=true,bool nlr=false,bool alL=false
	    ,double theta=-1.0):    
    pac(hps,dis,nprm,rsd,alL,theta),pop(loc),path_(n,sites),global_(n,sites),
    perfectlocalmatch(plm),nolocalrecomb(nlr)
#ifndef USE_R
    ,r(rsd) 
#endif
  {
    current_pop=pop.copy();    
    npops_=0;    
    for (int i=0;i<loc.dim();i++) 
      if (loc[i]>npops_) npops_=loc[i];    
    npops_++;  
    //  if (nolocalrecomb) std::cout << "No local Recomb \n";
    //  if (perfectlocalmatch) std::cout << "Perfect Local Copy\n";
  }  
  /** operator for a single value of rho                                */  
  double operator()(double rho, const TNT::Array1D<double> &F) {    
    return likelihood(rho,F);
  }  
  /** operator for a single value of rho                                */  
  double operator()(double rho, double F, double theta) {
    this->thetahat=theta;
    TNT::Array1D<double> Fv(npops(),F);

    //   std::cout << rho << " " << Fv << " " << theta << std::endl;

    return likelihood(rho,Fv);
  }  
  /** operator for an mcmc analysis                                     */  
  double operator()(double rho, const TNT::Array1D<double>  &F,
		    const std::pair<int,int> &changeorder);
  /** interface for optimisation algorithms                             */  
  double operator()(std::vector<double> &rho_f) {    
    TNT::Array1D<double> F(rho_f.size()-1,&rho_f[0]+1);    
    return likelihood(rho_f[0],F);  
  }  
  /** Number of populations                                             */
  int npops() const {    
    return npops_;  
  }  
  /** Accessor functions to return t the current path                   */
  TNT::Array2D<int> path() const {
    return path_;
  }
  /** Accessor function to return  the current "globalness"             */
  TNT::Array2D<bool> global() const {
    return global_;
  }
  /** Accessor function to return  the current "globalness"             */
  bool global(int i, int j) const {
    // what is the index of the starting order that corresponds to row i?
    int original = (h[i]-st[0])/sites;
    // std::cout << i << " " << original << std::endl;
    return global_[original][j];
  }
  /** get the array of "average" localness                              */  
  TNT::Array2D<double> localise(double rho, const TNT::Array1D<double> &F
				, double *recombinations=0,bool scale=true,double scalefac=1.0); 
  /** get the best path and localness                                   */
  TNT::Array2D<int> bestpath(double rho , const TNT::Array1D<double> &F 
			     ,std::vector<int> &perm);
  double single(double rho, const TNT::Array1D<double> &F, int seed=0);
  /** get an array of "samples" at the different positions              */  
  void sample(std::vector< TNT::Array3D<int> > &samp,const TNT::Array1D<double> &samppos
	      , int samplespercalc, double rho, const TNT::Array1D<double> &F,bool scale=true); 
 
  /** the calculation function - others are interfaces   */  
  double likelihood(double rho, const TNT::Array1D<double>  &F) {    
    calls_ +=1;    
    // use a rng to get the same set of random numbers for each run of the     
    // method - slightly slower but the time is not significant compared     
    // to the rest of the algorithm    
#ifdef USE_R
    rng rp;
#else
    rng rp(seed);
#endif       
    // and reset the ordering to its inital value    h=st;    
    current_pop=pop.copy();        
    if (avlogL) { 
      // do we average the likelihood or log-likelihood?      
      double sm=0.0;      
      for (int i=0;i<nperm;i++) {	
        std::vector<int> prms=rp.integer_permutations(n,n);	
        for (int j=0;j<n;j++) {	  
        h[j]=st[prms[j]];	  
        current_pop[j]=pop[prms[j]];	
        }	
        sm += calc(rho,F,path_,global_);      
      }      
      return sm-log(nperm);     
    } else {        
      std::vector<double> a(nperm);      
      for (int i=0;i<nperm;i++) {	
	  std::vector<int> prms=rp.integer_permutations(n,n);	
	  for (int j=0;j<n;j++) {	  
	    h[j]=st[prms[j]];	  
	    current_pop[j]=pop[prms[j]];	
	  }	
	  a[i] = calc(rho,F,path_,global_);      
	}      
      double scale=*(std::max_element(a.begin(),a.end()));      
      double sm=0.0;      
      for (int i=0;i<nperm;i++) sm += exp(a[i]-scale);      
      return log(sm) + scale - log(double(nperm));    
    }  
  }  
  /** sample cases and controls */
  TNT::Array3D<int> SampleCaseControlSNP(TNT::Array1D<double> &weightings
					 ,TNT::Array1D<int> &pos
					 ,TNT::Array1D<double> &d1, double rho 
					 ,const TNT::Array1D<double> &F
					 ,TNT::Array1D<int> &CaseControl,int seed);
protected:  
  const TNT::Array1D<int> &pop;  
  TNT::Array1D<int> current_pop;  
  int npops_;  
  TNT::Array2D<int> path_;										  
  TNT::Array2D<bool> global_;  	
  bool perfectlocalmatch,nolocalrecomb;
  rng r;  
  subdivpac(); // default constructor -- don't want it used  
  /** calculate the likelihood for a single rho for one   * permutation of the sample  */  
  double calc(double rho, const TNT::Array1D<double> &F, TNT::Array2D<int> &path	      
	    , TNT::Array2D<bool> &global,bool resamplepath=true); 
  TNT::Array2D<int> sampleSNP(double &lprob,double rho, const TNT::Array1D<double> & F
				       , TNT::Array1D<int> &pos
				       , TNT::Array1D<double> &d1, TNT::Array2D<int> &path
				       , TNT::Array2D<bool> &global,int minfreq);

};
#endif
