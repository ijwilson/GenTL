#ifndef STEPSTONE_H__
#define STEPSTONE_H__

#include <iosfwd>
#include <vector>

#include "gsl_rand.h"



/** The example mutation class  -stepping stone */
class StepStone {
public:
  StepStone(int nloci):d(nloci) {
  }
  StepStone(){};
  std::ostream &print(std::ostream &o) const;
    
  ~StepStone(){};
    std::vector<int> d;

  
  StepStone mutate(const std::vector<double> &theta,double time, rng &r) const {
    assert(nloci()==theta.size());
    StepStone tmp(nloci());
    for (size_t i=0;i<nloci();i++) {
      int muts=r.rpoisson(theta[i]*time*0.5);
      if (muts>0) tmp.d[i] = d[i]+2*r.rBinom(0.5,muts)-muts;
      else tmp.d[i]=d[i];
    }
    return tmp;
  }

  void initial(int loci,rng &r) {
    d.resize(loci);
    for (size_t i=0;i<nloci();i++) d[i]=10;
  }
    
  /** utility function for the number of loci   */
  size_t nloci() const{
    return d.size();
  }

};

std::ostream &operator<<(std::ostream &o, const StepStone &a);



#endif
