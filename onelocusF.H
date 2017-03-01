/// @file
#ifndef ONELOCUSF_H__
#define ONELOCUSF_H__

#include "tnt/tnt.h"
#include "gsl_rand.H"

/** One locus version of the PAC model.
 * 
 * One locus version of the pac model for subdivided
 * population, as a proof of concept that this works 
 * as well as Fst say */
class onelocusF {
public:
  onelocusF(const std::vector<std::string> &a
	,const TNT::Array1D<int> &loc);

  onelocusF(int *a, int *loc, int len);
  /** For a rewrite use the equations from the paper 
   * This is w (from Liu) which is the sequential probability 
   * of the data given the missing values - so independent of F 
   */  
  double w(double theta);
  double likelihood(double thet, bool *gg);
  /// likelihood when F = 0. 
  double likelihoodF0(double theta);
  double likelihoodresampleg(double *F, double theta,bool *g,double &lpsamp); 
 
  // double resampleg(double *F,rng &r);
  //double lpsampleg(double *F);

  double resampleg_b(double *F,double theta,rng &r);
  double lpsampleg_b(double *F,double theta);

  double lpg(double *F);
  
  double likelihoodb(double theta, double *F);

  /** Now for the straight calculation of the likelihood we need to first
   *  sample nperm missing data values
   */
  double likelihood( TNT::Array1D<double> &F, double theta, rng &r, int nperms);

  std::vector<std::vector<int> > Tabulate()  {
    int *st1 = location;
    int *st2 = x;
    return tabulate(st1,st1+n,st2,st2+n);

  }
  /** permute the order of sampling.   */
  void perm(rng &r) {
    permute(&order[0],n,r);
  }

  double lpg(double *F, bool *gg);
  // 
  /// return the number of populations in the sample
  int npops() const {
    return npops_;
  }
private:
  TNT::Array1D<int>  x;           
  TNT::Array1D<int>  location; 
  int n,k,npops_;
  // auxilary parameters
  TNT::Array1D<bool> g_;          /// is this individual global?
  TNT::Array1D<int> order;        /// the ordering of data
};

TNT::Array1D<int> convdata1(const std::vector<std::string> &a);

#endif
