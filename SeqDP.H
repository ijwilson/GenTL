#ifndef SEQDP_H__
#define SEQDP_H__

#include <vector>

#include "gsl_rand.H"
#include "gsl_distributions.H"

#include "flowSeq.H"
#include "flowGram.H"

/** An approximate sequential Dirichlet process for 
 * alignment                                          */

class SeqDP {
public:
  SeqDP(const std::vector<flowGram> &FG, ctsdistribution *Alphaprior, mvctsdistribution *zero
	, const std::vector< ctsdistribution *> &pp ,rng &myr, int maxcopies,const std::vector<int> &inittag)
    :tauprior(pp),zeroprior(zero),alphaprior(Alphaprior),r(myr),InitialTag(inittag),maxCopies(maxcopies)
    ,initialised(false) {
    for (size_t i=0;i<FG.size();i++) 
      fg.push_back(FG[i].intensity);
  }

  //  void add(const flowGram &f);
  void gibbs(int index);
  /** Call the genotypes using new parameters parameters                    */
  void initialise();
  
  double alpha() const {
    return alpha_;
  }
  /** return the number of different "alleles" */
  int nalls() const {
    return static_cast<int>(seq_.size());
  }
  /** returns the number of samples            */
  size_t n() const {
    return fg.size();
  }
  const std::vector<int> &copies() const {
    return copies_;
  }
  /** recall the genotypes based on the current clusters                   */
  void callGenotypes();
  /** Set the parameters                                                   */
  void SetParameters(double Alpha, const  std::vector<double> &Zeropars, const  std::vector<double> &Sigmapars)  {
    alpha_=Alpha;
    zeropars=Zeropars;
    assert(Sigmapars.size()==tauprior.size());
    sigmapars=Sigmapars;
  }


  /** update the sequencing parameters                                      */
  void gibbsUpdateSeqPars() {
    CalcPostPar();
    SampleParameters(); 
  };
  /** update a single sequence proportional to its probability             */

  void gibbsCallSequence(int which);

  /** The log-posterior for the current set of parameters     */
  double lposterior() const {
    return llikelihood()+lprior();
  }
  /**  Prior density for the parameters                      */
  double lprior() const;
  double llikelihood() const;


  /** printing functions                                                        */
  void printCluster(std::ofstream &oc) {
    printvector(oc,lab_,'\n');
  }
  void printSeq(std::ofstream &os,size_t minsize=1) {
    for (int ii=0;ii<nalls();ii++) {
      std::vector<int> cl = getCluster(ii);
	if (cl.size()>=minsize) {
	  os << "> ";
	  printvector(os,cl,false);
	  os << "\n";
	  os << seq_[ii].sequence() << std::endl;
	}
      }
  }  
  /** utility to print the parameters                        */
  std::string printParameters(const char *sep=" ") const {
    std::ostringstream oss;
    oss << alpha_ << sep;
    std::copy(zeropars.begin(),zeropars.end(),std::ostream_iterator<double>(oss,sep));
    std::copy(sigmapars.begin(),sigmapars.end(),std::ostream_iterator<double>(oss,sep));
    return oss.str();
  }
  void initialiseParameters();
  void sequentialStart(bool initParameters=true);

private:
  std::vector<flowSeq> seq_;
  std::vector<std::vector<double> > fg;
  std::vector<int> copies_;
  std::vector<int> lab_;

  double alpha_;
  //  std::vector<double> normgammpar,sigma;

  std::vector<ctsdistribution *> tauprior;    // priors
  mvctsdistribution *zeroprior;
  ctsdistribution *alphaprior;
  // the current sampled parameters
public:
  std::vector<double> zeropars;
  std::vector<double> sigmapars;              // model parameters
private:
  // the current posterior hyperparameters
  std::vector<double> zeropost;
  std::vector<std::vector<double> > taupost;
  rng &r;
  std::vector<int> InitialTag;                            // do we use an initial tag
                                                          // the initial tag insisted on,
  int maxCopies;                                          // but not printed 
  bool initialised;
 
private:
  // private functions
  void SampleParameters();
  void CalcPostPar();
 std::vector<int> getCluster(int wh) {
    assert(static_cast<size_t>(wh)<seq_.size());
    std::vector<int> s;
    for (size_t ii=0;ii<n();ii++) 
      if (lab_[ii]==wh) s.push_back(ii);
    return s;
 }
};



#endif
