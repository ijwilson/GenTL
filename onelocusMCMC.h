#ifndef MCMC_OL_H__
#define MCMC_OL_H__

#include  "gsl_rand.h"
#include "gsl_distributions.h"
#include "tnt/tnt.h"

class onelocusMCMC {
public:
  onelocusMCMC(const std::vector<std::string> &a, const TNT::Array1D<int> &loc
       , int seed,ctsdistribution *theta,ctsdistribution *f
       , int main, bool oneF=false);

  onelocusMCMC(int *x, int *loc, int len, ctsdistribution *theta,ctsdistribution *f);
  double likelihood() {
    return calc_likelihood(F,theta,g_);
  }
  double sampleg(TNT::Array1D<bool> &gg);
  double lpsampleg(TNT::Array1D<bool> &gg);
  double lpg(double *F, bool *gg);
  double likelihooda(double thet, bool *gg);
  bool check();
  //
  bool updateF(double tune);
  bool updateg();
  bool updatetheta(double tune);
  bool updateg2();
  bool updateg3();
  // priors
  ctsdistribution *thetaprior,*Fprior;
  TNT::Array1D<int> x;           
  TNT::Array1D<int> location; 
  int mainland;
  bool onef;
  int npops,n,k;   
  // auxilary parameters
  TNT::Array1D<bool> g_;          // global?
  TNT::Array1D<int> order;       // the ordering of data
  // parameters
  TNT::Array1D<double> F;
  double theta;
  // details
  rng r;
  // cache
  double llike, lprobg;

private:
  double calc_likelihood(double *F,double thet, bool *gg);
  void kickstart();
};



#endif
