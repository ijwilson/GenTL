/**               @file                                                  */
/** Templates for simple pact object                                     */
//   Time-stamp: <2012-01-27 17:57:21 ijw>                 
#ifndef PACT_H__
#define PACT_H__
#include "tnt/tnt.h"
#include "gsl_rand.h"
#include "GenEnums.h"

#include <deque>

/** Calculate product of approximate conditionals likelihood.
 * A class to incorporate the calculation of PAC 
 * likelihoods as given in Li & Stephens Genetics paper  
 * a templatised version that takes haplotypes as a type T
 * all that is required from this class is that it has: 
 *   - operator(ind,site)
 *   - sites()               -- number of sites
 *   - n()                   -- number of individuals
 *   - ismissing(ind,site)   -- is this SNP missing?
 *   - reset()               -- reset all missing SNPs
 *   - tostring(ind)         -- return a string of the haplotype
 */
template<typename T>
class pact {
public:  
  /** constructor for the pac object                                    */
  pact(T &hps, const TNT::Array1D<double> &dis, int nprm=10
       , const unsigned int rsd=1,double theta=-1.0):
    n(hps.n()),sites(hps.sites()),seed(rsd),nperm(nprm),calls_(0)
    ,haps(hps),d(dis),order(hps.n()) {
    for (int i=0;i<n;i++) order[i]=i;
    if (theta>0.0) thetahat=theta;
    else {
      theta=0.0;
      for (int i=1;i<n;i++) theta += 1./(double(i));
      // The value of thetabar given in Li and Stephens
      thetahat=1./theta;
    }
  };
  /** simple destructor                                                 */
  ~pact(){};
  /** operator for an mcmc analysis                                     */
  double operator()(double rho, std::pair<int,int> &changeorder) {
    int tmp = order[changeorder.first];
    order[changeorder.first] = order[changeorder.second];
    order[changeorder.second]=tmp;
    return calc(rho);
  }
  /** if a failure then put the order of h back                         */
  void reorder(std::pair<int,int> changeorder) {
    int tmp = order[changeorder.first];
    order[changeorder.first]=order[changeorder.second];
    order[changeorder.second]=tmp;
  }
  /** operator for a single value of rho                                */
  double operator()(double rho,double theta=-1.0) {
    if (theta>0) thetahat=theta;
    calls_ +=1;
    // use a rng to get the same set of random numbers for each run of the 
    // method - slightly slower but the time is not significant compared 
    // to the rest of the algorithm
    rng r(seed);   
    rng r2(seed);
    // and reset the ordering to its inital value
    for (int i=0;i<n;i++) order[i]=i;
    std::vector<double> a(nperm);
    for (int i=0;i<nperm;i++) {
      permute(&order[0],n,r);
      haps.reset();
      a[i] = calc(rho,r2);
    }
    double scale=*(std::max_element(a.begin(),a.end()));
    double sm=0.0;
    for (int i=0;i<nperm;i++) sm += exp(a[i]-scale);
    return log(sm) + scale - log(double(nperm));
  }
  /** sample from the pac model with parameter rho                        */
  TNT::Array2D<int> sample(const double &rho, int seed=1) {
    rng r(seed);
    TNT::Array2D<int> a(n,sites);
    SamplePath(a,rho,r);
    return a;
  }
  /** get an estimate of the numbers and positions of recombinations      */
  TNT::Array1D<double> recombinations(const double &rho, int seed=1
                                      , bool scale=true);
  /** print out the haplotypes and distances                              */
  void print(std::ostream &o) {
    o << haps;
    o << std::endl;
    copy(&d[0],&d[0]+d.dim(),std::ostream_iterator<double>(o," "));
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
      hps.insert(haps.tostring(i));
    }
    return hps.size();
  }
  /** sampleSNP for case control data (in genotypes)                    */
  TNT::Array3D<int> SampleCaseControlSNP(TNT::Array1D<double> &weightings
                                         ,int nsamp, TNT::Array1D<int> &pos
                                         , TNT::Array1D<double> &d1, double rho
                                         , TNT::Array1D<int> &CaseControl
                                         , int seed);
  /** utility functions - so as to avoid explicit calls for children   */
  int nsites() const {
    return sites;
  }
  int ss() const {
    return n;
  }
public:
  int n,sites,seed,nperm; 
  int calls_;                           //< how many calls have been made
  double thetahat;                      //< what value of thetahat do we use?
  T &haps;                              //< Reference to the haplotype 
  const TNT::Array1D<double>  &d;       //< ref to distances between sites
public:
  TNT::Array1D<int> order;              //< current ordering of the sites
protected:
  //
  pact();                               //< default ctor, don't want it used
  /** calculate the likelihood for a single rho for one
   * permutation of the sample                          */
  double getalpha(double rho,TNT::Array2D<double> &alpha, double &gammamatch, int k);
  double getalpha(double rho,TNT::Array2D<double> &alpha, double &gammamatch, int k
                  ,std::deque<int> &missingvals);
  double calc(double rho); 
  double calc(double rho,rng &r);
  double SamplePath(TNT::Array2D<int> &a,const double &rho, rng &r);
  double SamplePathMissing(TNT::Array2D<int> &a,const double &rho, rng &r);
  TNT::Array3D<int> SampleSNP(double &lprob,int nsamp, TNT::Array1D<int> &pos
                              , TNT::Array1D<double> &d1
                              , const double &rho, rng &r, int minfreq);  
  TNT::Array3D<int> SampleSNPMissing(double &lprob,int nsamp
                                     , TNT::Array1D<int> &pos
                                     , TNT::Array1D<double> &d1
                                     , const double &rho, rng &r
                                     , int minfreq);
};

void getequalspacedpositionsanddistances(const TNT::Array1D<double> &distances
                                         , TNT::Array1D<int> &pos
                                         ,TNT::Array1D<double> &gap
                                         ,TNT::Array1D<double> &x);

// utility function to rescale probabilities so that the total probability is 1
// returns the logarithm of the rescaling factor for rescaling

/// The function object multiplies an element by a Factor
template <class Type>
struct DivValue
{
public:
  // Constructor initializes the value to multiply by
  DivValue ( const Type& _Val ) : Factor ( 1./_Val ) {
  }
  // The function call for the element to be multiplied
  void operator ( ) ( Type& elem ) const
  {
    elem *= Factor;
  }
  Type Factor;   // The value to multiply by
};

template <typename T>
T rescalesm(T *a, int count) 
{   
  T scalefac=std::accumulate(a,a+count,0.0);
  assert(scalefac>0.0);

  std::for_each(a,a+count,DivValue<T>(scalefac));
    
  return log(scalefac);
}
/****************************************************************************/
/*  calculate the likelihood for a single arrangement of the samples          
 * should only be called when there is no missing data                      */
/****************************************************************************/
template<typename T>
double pact<T>::calc(double rho)
{
  // alpha[x,s] defined to be the probability of sampling the first s sites of
  // your haplotype, and of sampling haplotype x
  if (haps.nmissing()>0) 
    throw std::domain_error("This calc version only works"
                            " for non-missing data"); 
  TNT::Array2D<double> alpha(sites,n-1);
  double lprob=0.0,gammamatch;
  // now consider each of the subsequent samples in turn.
  // Because my arrays are 0 offset haplotype k is the k+1 st sample
  // so use the same k as in Ni and Stephens' appendix 
  for (int k=1;k<n;k++) {
    double lscalefactor=getalpha(alpha,gammamatch,k);
    double prob=std::accumulate(alpha[sites-1],alpha[sites-1]+k,0.0);
    lprob+=log(prob)+lscalefactor;
  }
  return lprob;
}
/****************************************************************************/
/** get an estimate of the numbers and positions of recombinations          */
/****************************************************************************/
template<typename T>
TNT::Array1D<double> pact<T>::recombinations(const double &rho, 
                                             int seed, bool scale) 
{
  rng r(seed);  
  for (int i=0;i<n;i++) order[i]=0;
  TNT::Array1D<double> recombs(sites-1,0.0);
  TNT::Array2D<int> recs(nperm,sites-1,0);
  TNT::Array2D<int> a(n,sites);
  std::vector<double> ll(nperm);
    
  for (int i=0;i<nperm;i++) {
    ll[i]=SamplePath(a,rho,r);
    for (int j=0;j<n;j++)
      for (int k=0;k<sites-1;k++)
        recs[i][k]+=int(a[j][k]!=a[j][k+1]);
  }
  if (!scale ) for (int i=0;i<nperm;i++) ll[i]=0.0;	  
  double sc=*(std::max_element(ll.begin(),ll.end()));
  for (int i=0;i<nperm;i++) ll[i] = exp(ll[i]-sc);  
  double divby=std::accumulate(ll.begin(),ll.end(),0.0);
    
  for (int i=0;i<nperm;i++) { 
    for (int s=0;s<sites-1;s++) 
      recombs[s]+=ll[i]*recs[i][s];
  }
  std::for_each(&recombs[0],&recombs[0]+sites-1,DivValue<double>(divby));
  return recombs;
}
/****************************************************************************/
/** calculate the likelihood and sample the paths of those that are copied  
 * again, only called when there is no missing data                         */
/****************************************************************************/
template<typename T>
double pact<T>::SamplePath(TNT::Array2D<int> &samp,const double &rho, rng &r)
{
  if (haps.nmissing()>0) 
    throw std::domain_error("This SamplePath only for use without missing data");
  // alpha[x,s] defined to be the probability of sampling the first s sites of
  // your haplotype, and of sampling haplotype x 
  TNT::Array2D<double> alpha(sites,n-1);
  // the first sample has probability 1
  double lprob=0.0;
  double gammamatch;
  // now consider each of the subsequent samples in turn.
  // Because my arrays are 0 offset haplotype k is the k+1 st sample
  // so use the same k as in Li and Stephens' appendix 
  for (int k=1;k<n;k++) {
    double lscalefac=getalpha(rho,alpha,gammamatch,k);
    // calculate the prob for haplotype k by summing over all the 
    // paths to the last site
    double prob=std::accumulate(alpha[sites-1],alpha[sites-1]+k,0.0);
    lprob+=log(prob)+lscalefac;
    // and now for the "backwards" part of the algorithm, 
    // as in Boys et al 2000 - the last site is "site-1"
    std::vector<double> lp(k);
    samp[k][sites-1] = gen_from_p(alpha[sites-1],alpha[sites-1]+k,r);
    for (int s=sites-2;s>=0;s--) {
      for (int x=0;x<k;x++) {
        double p_stay = exp(-rho*d[s]/k);
        lp[x]=alpha[s][x];
        // is it the same as the "next" site
        if (samp[k][s+1]==x) lp[x] *= (p_stay+(1.-p_stay)/k); 
        else lp[x] *= (1.-p_stay)/k; 
      }
      samp[k][s]=gen_from_pr(lp,k,r);
    }
  }
  return lprob;
}
/** routine to calculate alpha for a particular k                     */
template<typename T> 
double pact<T>::getalpha(double rho, TNT::Array2D<double> &alpha,double &gammamatch, int k)
{
  // some quantities we only need to calculate once for each k
  gammamatch=double(k)/(k+thetahat)+0.5*thetahat/(k+thetahat);
  double gammamiss=1.-gammamatch;   
  // we can calculate the probabilities for the first site directly 
  // loop over the haplotypes from which the site could be copied
  for (int x=0;x<k;x++) {
    if (haps.value(order[x],0)==haps.value(order[k],0)) 
      alpha[0][x]=gammamatch/double(k);
    else alpha[0][x]=gammamiss/double(k);
  }
  // fix as in MacDonald and Zucchini so that the mean alpha is 1
  double logscalefactor=rescalesm(alpha[0],k);    
  // calculate alpha for the rest of the sites
  for (int s=1;s<sites;s++) {
    // the "second term" from Li and Stephens 
    double sm=std::accumulate(alpha[s-1],alpha[s-1]+k,0.0);
    sm /= (double)(k);       
    double p=exp(-rho*d[s-1]/(double(k)));
    
    for (int x=0;x<k;x++) {
      alpha[s][x]=p*alpha[s-1][x]+(1-p)*sm;
      if (haps.value(order[x],s)==haps.value(order[k],s)) 
        alpha[s][x]*=gammamatch;
      else alpha[s][x]*=gammamiss;
    }
    // the scale factor of MacDonald and Zucchini
    logscalefactor += rescalesm(alpha[s],k);
  }
  return logscalefactor;
}
/** routine to calculate alpha for a particular k keeping track of  missing values   */
template<typename T> 
double pact<T>::getalpha(double rho,TNT::Array2D<double> &alpha,double &gammamatch
                         , int k,std::deque<int> &missingvals)
{
  assert(missingvals.size()==0);
  // some quantities we only need to calculate once for each k
  gammamatch=double(k)/(k+thetahat)+0.5*thetahat/(k+thetahat);
  double gammamiss=1.-gammamatch;     // = 0.5*thetahat/(k+thetahat);
  // we can calculate the probabilities for the first site directly 
  // loop over the haplotypes from which the site could be copied
  if (haps.ismissing(order[k],0)) { 
    missingvals.push_back(0);
    for (int x=0;x<k;x++) {
      alpha[0][x]=1./double(k);
    }
  } else {
    for (int x=0;x<k;x++) {
      if (haps.value(order[x],0)==haps.value(order[k],0)) 
        alpha[0][x]=gammamatch/double(k);
      else alpha[0][x]=gammamiss/double(k);
    }
  }
  // fix as in MacDonald and Zucchini so that the mean alpha is 1
  double logscalefactor=rescalesm(alpha[0],k);    
  // calculate alpha for the rest of the sites
  for (int s=1;s<sites;s++) {
    // the "second term" from Li and Stephens 
    double sm=std::accumulate(alpha[s-1],alpha[s-1]+k,0.0)/(double(k));
    double p=exp(-rho*d[s-1]/(double(k)));
    if (haps.ismissing(order[k],s)) {
      missingvals.push_back(s);
      for (int x=0;x<k;x++) {
        alpha[s][x]=p*alpha[s-1][x]+(1-p)*sm;
      }
    } else {
      for (int x=0;x<k;x++) {
        alpha[s][x]=p*alpha[s-1][x]+(1-p)*sm;
        if (haps.value(order[x],s)==haps.value(order[k],s)) 
          alpha[s][x]*=gammamatch;
        else alpha[s][x]*=gammamiss;
      }
    }
    // the scale factor of MacDonald and Zucchini
    logscalefactor+=rescalesm(alpha[s],k);
  }
  return logscalefactor;
}
/*****************************************************************************/
/*  calculate the likelihood for a single arrangement of the samples          
 *  should only be called when there is missing data                         */
/*****************************************************************************/
template <typename T>
double pact<T>::calc(double rho,rng &r)
{
  // alpha[x,s] defined to be the probability of sampling the first s sites of
  // your haplotype, and of sampling haplotype x 
  TNT::Array2D<double> alpha(sites,n-1);
  std::vector<double> lprob(n,0.0);
  // the first sample has probability 1
  lprob[0]=0.0;
  double gammamatch;
  for (int s=0;s<sites;s++) {
    if (haps.ismissing(order[0],s)) {
      if (r()<0.5) 
        haps.change(order[0],s,zero);
      else haps.change(order[0],s,one);
    }
  }
  // now consider each of the subsequent samples in turn.
  // Because my arrays are 0 offset haplotype k is the k+1 st sample
  // so use the same k as in Ni and Stephens' appendix 
  for (int k=1;k<n;k++) {
    std::deque<int> missingvals;
    double lscalefactor = getalpha(rho,alpha,gammamatch,k,missingvals);
    // calculate the prob for haplotype k by summing over all the 
    // paths to the last site
    double prob=std::accumulate(alpha[sites-1],alpha[sites-1]+k,0.0);
    lprob[k]=log(prob)+lscalefactor;
    // now sample the missing sites
    if (missingvals.size()>0) {
      std::vector<double> lp(k);
      // first sample the last missing site
      int lst=missingvals.back();
      int samp=gen_from_p(alpha[lst],alpha[lst]+k,r);
      if (r() < gammamatch) 
        haps.change(order[k],lst,haps.value(order[samp],lst));
      else  
        haps.change(order[k],lst,
                    flip(haps.value(order[samp],lst)));
      if (missingvals.size()>=1) {
        for (int s=lst-1;s>=missingvals[0];s--) {
          for (int x=0;x<k;x++) {
            double p_stay = exp(-rho*d[s]/k);
            lp[x]=alpha[s][x];
            // is it the same as the "next" site
            if (samp==x) lp[x] *= (p_stay+(1.-p_stay)/k); 
            else lp[x] *= (1.-p_stay)/k; 
          }
          samp=gen_from_pr(lp,k,r);
          if (haps.ismissing(order[k],s)) {
            if (r() < gammamatch) 
              haps.change(order[k],s,haps.value(order[samp],s));
            else  
              haps.change(order[k],s,flip(haps.value(order[samp],s)));
          }
        }
      }
    }
  }
  assert(haps.check());

  return std::accumulate(lprob.begin(),lprob.end(),0.0);
}
/****************************************************************************/
/** calculate the likelihood and sample the paths of those that are copied  
 * again, only called when there is no missing data                         */
/****************************************************************************/
template<typename T>
double pact<T>::SamplePathMissing(TNT::Array2D<int> &samp,const double &rho, rng &r)
{
  // alpha[x,s] defined to be the probability of sampling the first s sites of
  // your haplotype, and of sampling haplotype x 
  TNT::Array2D<double> alpha(sites,n-1);
  // the first sample has probability 1
  double lprob=0.0,gammamatch;
  for (int s=0;s<sites;s++) {
    if (haps.ismissing(order[0],s)) {
      if (r()<0.5) 
        haps.change(order[0],s,zero);
      else haps.change(order[0],s,one);
    }
  }
  // now consider each of the subsequent samples in turn.
  // Because my arrays are 0 offset haplotype k is the k+1 st sample
  // so use the same k as in Li and Stephens' appendix 
  for (int k=1;k<n;k++) {
    std::deque<int> missingvals;    
    double logscalefactor=getalpha(rho,alpha,gammamatch,k,missingvals);
    // calculate the prob for haplotype k by summing over all the 
    // paths to the last site
    double prob=std::accumulate(alpha[sites-1],alpha[sites-1]+k,0.0);
    lprob+=log(prob)+logscalefactor;
    // and now for the "backwards" part of the algorithm, 
    // as in Boys et al 2000 - the last site is "site-1"
    std::vector<double> lp(k);
    samp[k][sites-1] = gen_from_p(alpha[sites-1],alpha[sites-1]+k,r);
    for (int s=sites-2;s>=0;s--) {
      for (int x=0;x<k;x++) {
        double p_stay = exp(-rho*d[s]/k);
        lp[x]=alpha[s][x];
        // is it the same as the "next" site
        if (samp[k][s+1]==x) lp[x] *= (p_stay+(1.-p_stay)/k); 
        else lp[x] *= (1.-p_stay)/k; 
      }
      samp[k][s]=gen_from_pr(lp,k,r);
      if (haps.ismissing(order[k],s)) {
        if (r() < gammamatch) 
          haps.change(order[k],s,haps.value(order[samp],s));
        else  
          haps.change(order[k],s,flip(haps.value(order[samp],s)));
      } 
    }
  }
  return lprob;
}
/*****************************************************************************/
/** Get samples of an unseen SNP between the pos and pos+1 SNP with 
 *   relative distances d1                                                   */
/*****************************************************************************/
template<typename T>
TNT::Array3D<int> pact<T>::SampleSNPMissing(double &lprob,int nsamp
                                            , TNT::Array1D<int> &pos
                                            , TNT::Array1D<double> &d1
                                            , const double &rho, rng &r
                                            , int minfreq)
{
  int nSNP=pos.dim();
  TNT::Array3D<int> path(nSNP,nsamp,n);
  TNT::Array1D<double> gammamatch(n);
  // alpha[x,s] defined to be the probability of sampling the first s sites of
  // your haplotype, and of sampling haplotype x 
  TNT::Array2D<double> alpha(sites,n-1);
  // the first sample has probability 1
  lprob=0.0;
  for (int s=0;s<sites;s++) {
    if (haps.ismissing(order[0],s)) {
      if (r()<0.5) 
        haps.change(order[0],s,zero);
      else haps.change(order[0],s,one);
    }
  }
  for (int SNP=0;SNP<nSNP;SNP++) 
    for (int jj=0;jj<nsamp;jj++) 
      path[SNP][jj][0]=0;
  // now consider each of the subsequent samples in turn.
  // Because my arrays are 0 offset haplotype k is the k+1 st sample
  // so use the same k as in Li and Stephens' appendix 
  for (int k=1;k<n;k++) {
    std::deque<int> missingvals;    
    double logscalefactor=getalpha(rho,alpha,gammamatch[k],k,missingvals);

    // calculate the prob for haplotype k by summing over all the 
    // paths to the last site
    double prob=std::accumulate(alpha[sites-1],alpha[sites-1]+k,0.0);
    lprob+=log(prob)+logscalefactor;
    std::vector<double> lp(k);
    TNT::Array1D<int> samp(sites);
    samp[sites-1] = gen_from_p(alpha[sites-1],alpha[sites-1]+k,r);
    for (int s=sites-2;s>=0;s--) {
      for (int x=0;x<k;x++) {
        double p_stay = exp(-rho*d[s]/k);
        lp[x]=alpha[s][x];
        // is it the same as the "next" site
        if (samp[s+1]==x) lp[x] *= (p_stay+(1.-p_stay)/k); 
        else lp[x] *= (1.-p_stay)/k; 
      }
      samp[s]=gen_from_pr(lp,k,r);
      if (haps.ismissing(order[k],s)) {
        if (r() < gammamatch[k]) 
          haps.change(order[k],s,haps.value(order[samp[s]],s));
        else  
          haps.change(order[k],s,flip(haps.value(order[samp[s]],s)));
      } 
    }
    for (int SNP=0;SNP<nSNP;SNP++) {
      int position=pos[SNP];// SNP lies between columns position and position+1
      double d2=d[position]-d1[SNP]; 
      assert(d1[SNP]<d[position]); 
      assert(fabs(d1[SNP]+d2-d[position])<1E-8);
      TNT::Array1D<double> newalpha(k);
      // start by recalculating for the different values
      // the "second term" from Li and Stephens
      double sm=std::accumulate(alpha[position]
                                ,alpha[position]+k,0.0)/(double(k));
      double p=exp(-rho*d1[SNP]/(double(k)));
      for (int x=0;x<k;x++) 
        newalpha[x] = p*alpha[position][x]+(1-p)*sm;
      
      // and now for the "backwards" part of the algorithm   
      std::vector<double> lp(k);
      for (int jj=0;jj<nsamp;jj++) {
        // sample a path "above"
        for (int x=0;x<k;x++) {
          double p_stay = exp(-rho*d2/k);
          lp[x]=newalpha[x];
          if (samp[position+1]==x) lp[x] *= (p_stay+(1.-p_stay)/k); 
          else lp[x] *= (1.-p_stay)/k; 
        }
        path[SNP][jj][k] = gen_from_pr(lp,k,r);
      }
    }
  }
  TNT::Array3D<int> sampSNP(nSNP,nsamp,n);
  // alpha[x,s] defined to be 
  for (int SNP=0;SNP<nSNP;SNP++) {
    for (int jj=0;jj<nsamp;) {
      bool keep=false;
      if (r()<0.5) sampSNP[SNP][jj][0]=1;
      else sampSNP[SNP][jj][0]=0;
      for (int k=1;k<n;k++) {
        if (r()<gammamatch[k]) 
          sampSNP[SNP][jj][k]=sampSNP[SNP][jj][path[SNP][jj][k]];
        else {
          sampSNP[SNP][jj][k]=1-sampSNP[SNP][jj][path[SNP][jj][k]]; //mutation 
          keep=true;  // we can try to keep this set!
        }
      }
      if (keep) {
        int freq=0;
        for (int i=0;i<n;i++) freq+=sampSNP[SNP][jj][i]==1;
        if (freq>=minfreq&&n-freq> minfreq) jj++;// keep this and get the next
      }
    }
  }
  return sampSNP;
}
/*****************************************************************************/
/** Get samples of an unseen SNP between the pos and pos+1 SNP with 
 *   relative distances d1                                                   */
/*****************************************************************************/
template<typename T>
TNT::Array3D<int> pact<T>::SampleSNP(double &lprob,int nsamp
                                     , TNT::Array1D<int> &pos
                                     , TNT::Array1D<double> &d1
                                     , const double &rho, rng &r, int minfreq)
{
  lprob=0.0;
  assert(pos.dim()==d1.dim()); 
  int nSNP=pos.dim();
  TNT::Array3D<int> path(nSNP,nsamp,n);
  TNT::Array1D<double> gammamatch(n);
  // alpha[x,s] defined to be the probability of sampling the first s sites of
  // your haplotype, and of sampling haplotype x 
  TNT::Array2D<double> alpha(sites,n-1);
  // the first sample has probability 1
  // double lprob=0.0;
  for (int SNP=0;SNP<nSNP;SNP++) 
    for (int jj=0;jj<nsamp;jj++) 
      path[SNP][jj][0]=0;
  // now consider each of the subsequent samples in turn.
  // Because my arrays are 0 offset haplotype k is the k+1 st sample
  // so use the same k as in Li and Stephens' appendix 
  for (int k=1;k<n;k++) {
    double logscalefactor=getalpha(rho,alpha,gammamatch[k],k);
    // some quantities we only need to calculate once for each k
    // calculate the prob for haplotype k by summing over all the 
    // paths to the last site
    double prob=std::accumulate(alpha[sites-1],alpha[sites-1]+k,0.0);
    lprob+=log(prob)+logscalefactor;
    // std::cerr << lprob << std::endl;
    // and now for the "backwards" part of the algorithm, 
    // as in Boys et al 2000 - the last site is "site-1"
    std::vector<double> lp(k);
    TNT::Array1D<int> samp(sites);
    samp[sites-1] = gen_from_p(alpha[sites-1],alpha[sites-1]+k,r);
    for (int s=sites-2;s>=0;s--) {
      for (int x=0;x<k;x++) {
        double p_stay = exp(-rho*d[s]/k);
        lp[x]=alpha[s][x];
        // is it the same as the "next" site
        if (samp[s+1]==x) lp[x] *= (p_stay+(1.-p_stay)/k); 
        else lp[x] *= (1.-p_stay)/k; 
      }
      samp[s]=gen_from_pr(lp,k,r);
    }
    for (int SNP=0;SNP<nSNP;SNP++) {
      int position=pos[SNP]; // SNP lies between columns position and position+1
      double d2=d[position]-d1[SNP]; 
      assert(d1[SNP]<d[position]); 
      assert(fabs(d1[SNP]+d2-d[position])<1E-8);
      TNT::Array1D<double> newalpha(k);
      // start by recalculating for the different values
      // the "second term" from Li and Stephens
      double sm=std::accumulate(alpha[position]
                                ,alpha[position]+k,0.0)/(double(k));
      double p=exp(-rho*d1[SNP]/(double(k)));
      
      for (int x=0;x<k;x++) 
        newalpha[x] = p*alpha[position][x]+(1-p)*sm;
      
      // and now for the "backwards" part of the algorithm   
      std::vector<double> lp(k);
      for (int jj=0;jj<nsamp;jj++) {
        // sample a path "above"
        for (int x=0;x<k;x++) {
          double p_stay = exp(-rho*d2/k);
          lp[x]=newalpha[x];
          if (samp[position+1]) lp[x] *= (p_stay+(1.-p_stay)/k); 
          else lp[x] *= (1.-p_stay)/k; 
        }
        path[SNP][jj][k] = gen_from_pr(lp,k,r);
      }
    }
  }
  TNT::Array3D<int> sampSNP(nSNP,nsamp,n);
  // alpha[x,s] defined to be 
  for (int SNP=0;SNP<nSNP;SNP++) {
    for (int jj=0;jj<nsamp;) {
      bool keep=false;
      if (r()<0.5) sampSNP[SNP][jj][0]=1;
      else sampSNP[SNP][jj][0]=0;
      for (int k=1;k<n;k++) {
        if (r()<gammamatch[k]) 
          sampSNP[SNP][jj][k]=sampSNP[SNP][jj][path[SNP][jj][k]];
        else {
          sampSNP[SNP][jj][k]=1-sampSNP[SNP][jj][path[SNP][jj][k]]; //mutation 
          keep=true;  // we can try to keep this set!
        }
      }
      if (keep) {
        int freq=0;
        for (int i=0;i<n;i++) freq+=sampSNP[SNP][jj][i]==1;
        if (freq>=minfreq&&n-freq> minfreq) jj++;// keep this and get the next
      }
    }
  }
  return sampSNP;
}
/*****************************************************************************/
/** sampleSNP for case control data (in genotypes)                           */
/* The ordering of genotypes is (if controls are labelled 0 and cases 1)
 *   
 *      Control_00 Control_01 Control_11 Case_00 Case_01 Case_11
 *  
 *****************************************************************************/
template <typename T>
TNT::Array3D<int> 
pact<T>::SampleCaseControlSNP(TNT::Array1D<double> &weightings
                              ,int nsamp, TNT::Array1D<int> &pos
                              ,TNT::Array1D<double> &d1, double rho
                              , TNT::Array1D<int> &CaseControl,int seed)
{
  rng r(seed);
  int nSNP=pos.dim();
  TNT::Array3D<int> tables(pos.dim(),nsamp*nperm,6,0);
  std::vector<int> rev(n);    
  haps.reset();
  for (int i=0;i<nperm;i++) {
    for (int j=0;j<n;j++) order[j]=j;
    permute(&order[0],n,r);
    for (int j=0;j<n;j++) rev[order[j]]=j;
    haps.reset();
    TNT::Array3D<int> res =  SampleSNPMissing(weightings[i],nsamp,pos,d1,rho,r,n/10);

    for (int SNP=0;SNP<nSNP;SNP++) {
      for (int j=0;j<nsamp;j++) {
        for (int k=0;k<n;k+=2) {
          if (res[SNP][j][rev[k]]==0&&res[SNP][j][rev[k+1]]==0) 
            tables[SNP][i*nsamp+j][CaseControl[k]*3]+=1;
	    
          else if (res[SNP][j][rev[k]]==1&&res[SNP][j][rev[k+1]]==1) 
            tables[SNP][i*nsamp+j][2+CaseControl[k]*3]+=1;
	    
          else  tables[SNP][i*nsamp+j][1+CaseControl[k]*3]+=1;
        }
      }
    }
  }    
  double scale=*(std::max_element(&weightings[0],&weightings[0]+nperm));
  for (int i=0;i<nperm;i++) 
    weightings[i] = exp(weightings[i]-scale);
    
  return tables;
}

#endif
