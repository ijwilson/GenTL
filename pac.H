// @file   Time-stamp: <05/12/15 18:59:10 ijw>
#ifndef PAC_H__
#define PAC_H__
#include "tnt/tnt.h"
#include "gsl_rand.H"
#include <vector>
#include <iostream>
#include <iterator>  // for ostream_iterator
#include <algorithm> // for max_element
#include <numeric>   // for accumulate
#include <set>
/** A class to incorporate the calculation of PAC 
 * likelihoods as given in Li & Stephens Genetics paper  */
class pac {
public:
  /** constructor for the pac object                                      */
  pac( TNT::Array2D<int> &hps, const std::vector<double> &dis, int nprm=10
      , const unsigned int rsd=1, bool alL=false,double theta=-1.0):
    n(hps.dim1()),sites(hps.dim2()),seed(rsd),nperm(nprm),calls_(0)
    ,avlogL(alL),haps(hps),h(hps.dim1()),d(dis) {
    if (theta>0.0) thetahat=theta;
    else {
      theta=0.0;
      for (int i=1;i<n;i++) theta += 1./(double(i));
      // The value of thetabar given in Li and Stephens
      thetahat=1./theta;
    }
    for (int i=0;i<n;i++) h[i]=haps[i];
    st=h;  // the starting order
  };
  /** simple destructor                                                   */
  ~pac(){};
  /** operator for an mcmc analysis                                     */
  double operator()(double rho, std::pair<int,int> &changeorder) {
    int *tmp = h[changeorder.first];
    h[changeorder.first]=h[changeorder.second];
    h[changeorder.second]=tmp;
    return calc(rho);
  }
  /** if a failure then put the order of h back                         */
  void reorder(std::pair<int,int> changeorder) {
    int *tmp = h[changeorder.first];
    h[changeorder.first]=h[changeorder.second];
    h[changeorder.second]=tmp;
  }
  /** operator for a single value of rho                                */
  double operator()(double rho,double theta=-1.0) {
    if (theta>0) thetahat=theta;
    calls_ +=1;
    // use a rng to get the same set of random numbers for each run of the 
    // method - slightly slower but the time is not significant compared 
    // to the rest of the algorithm
#ifdef USE_R
    rng r;
#else
    rng r(seed);
#endif   
    // and reset the ordering to its inital value
    reset();
    if (avlogL) { // do we average the likelihood or log-likelihood?
      double sm=0.0;
      for (int i=0;i<nperm;i++) {
	nextperm(r);
	sm += calc(rho);
      }
      return sm-log(nperm); 
    } else {  
      std::vector<double> a(nperm);
      for (int i=0;i<nperm;i++) {
	nextperm(r);
	a[i] = calc(rho);
      }
      double scale=*(std::max_element(a.begin(),a.end()));
      double sm=0.0;
      for (int i=0;i<nperm;i++) sm += exp(a[i]-scale);
      return log(sm) + scale - log(double(nperm));
    }
  }
  /** sample from the pac model with parameter rho                        */
  TNT::Array2D<int> sample(const double &rho, int seed=1) {
#ifdef USE_R
    rng r;
#else
    rng r(seed);
#endif
    TNT::Array2D<int> a(n,sites);
    SamplePath(a,rho,r);
    return a;
  }
  /** get an estimate of the numbers and positions of recombinations                   */
  TNT::Array1D<double> recombinations(const double &rho, int seed=1, bool scale=true);
  /** print out the haplotypes and distances                              */
  void print(std::ostream &o) {
    o << haps;
    o << std::endl;
    copy(d.begin(),d.end(),std::ostream_iterator<double>(o," "));
    o << std::endl;
  }
  /** how many calls have been made ? */
  int calls() {
    return calls_;
  }
  /** Number of distinct haplotypes                                      */
  int distincthaps() {
    std::set<std::string> hps;
    for (int i=0;i<n;i++) {
      std::ostringstream oss;
      std::copy(haps[i],haps[i]+sites,std::ostream_iterator<int>(oss,""));
      hps.insert(oss.str());
    }
    return hps.size();
  }
  /** sampleSNP for case control data (in genotypes)                    */
  TNT::Array3D<int> SampleCaseControlSNP(TNT::Array1D<double> &weightings
					 ,int nsamp, TNT::Array1D<int> &pos, TNT::Array1D<double> &d1
					 , double rho, TNT::Array1D<int> &CaseControl, int seed);
protected:
  int n,sites,seed,nperm; 
  int calls_;                                //< how many calls have been made
  bool avlogL;                               //< average logL or L?
  double thetahat;                           //< what value of thetahat do we use?
  TNT::Array2D<int> &haps;                   //< Reference to the haplotype matrix
  std::vector<int *> h;                      //< current order of haplotypes
  std::vector<int *> st;                     //< the starting order of haplotypes
  const std::vector<double> &d;              //< reference to the distances between sites
  pac();                                     //< default constructor -- don't want it used
  /** calculate the likelihood for a single rho for one
   * permutation of the sample                          */
  double calc(double rho);  
  double SamplePath(TNT::Array2D<int> &a,const double &rho, rng &r);
  TNT::Array3D<int> SampleSNP(double &lprob,int nsamp, TNT::Array1D<int> &pos
				   , TNT::Array1D<double> &d1, const double &rho
				   , rng &r, int minfreq);
  /** rest the order to the original                                                    */
  void reset() {
    h=st;
  }
  /** Get the next set of orderings for the haplotypes                                  */
  void nextperm(rng &r) {
    std::vector<int> prms=r.integer_permutations(n,n);
    for (int j=0;j<n;j++) {
      h[j]=st[prms[j]];
    }
  }
};

/** a pac class that allows missing data                                                */
class pacmissing: public pac {
  pacmissing(TNT::Array2D<int> &hps, const std::vector<double> &dis
	     , int nprm=10, const unsigned int rsd=1, bool alL=false
	     ,double theta=-1.0)
    :pac(hps,dis,nprm,rsd,alL,theta),path(n,sites)
#ifndef USE_R
    ,localr(1)
#endif
 {
    for (int i=0;i<n;i++)
      for (int j=0;j<sites;j++) 
	if (h[i][j]<0) miss_.push_back(std::pair<int,int>(i,j));
	  
  }
private:
  std::vector<std::pair<int,int> > miss_;  // which data are missing - pair (individual,site)
  double SamplePath(TNT::Array2D<int> &samp,const double &rho, rng &r);
  TNT::Array3D<int> SampleSNP(double &lprob,int nsamp, TNT::Array1D<int> &pos
			      , TNT::Array1D<double> &d1, const double &rho, rng &r, int minfreq);
  double calc(double rho) {
    return SamplePath(path,rho,localr);
  }
  TNT::Array2D<int> path;
  rng localr;
  // functions
  void reset() {
    pac::reset();  // reset the ordering 
    for (size_t i=0;i<miss_.size();i++) 
      h[miss_[i].first][miss_[i].second]=-1;
  }
};

double rescale(double *a, int count);
void getequalspacedpositionsanddistances(const TNT::Array1D<double> &distances
					 , TNT::Array1D<int> &pos,TNT::Array1D<double> &gap,TNT::Array1D<double> &x);
#endif
