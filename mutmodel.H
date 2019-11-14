#ifndef MUTMODEL_H__
#define MUTMODEL_H__

#include <vector>

#include "gsl_distributions.h"
#include "gsl_random.h"

namespace GenTL {

  enum bases  {
    A=1, G=2, C=4, T=5
  };


  class sequenceMutationModel {
    sequenceMutationModel(ctsdistribution &mu, size_t sequence_length,rng &r)
      :overall_mu(mu),local_mu(length) {
      for (size_t i=0;i<sequence_length;i++) local_mu[i]=mu.sample(r);
    
    };
    void mutate(std::vector<bases> &below, const std::vector<bases> above, 
		double time)=0;
  private:
    ctsdistribution &overall_mu;
    std::vector<double> local_mu;
    rng &lr;
  };


  class JCmutation:public sequenceMutationModel {
    FCmutation(ctsdistribution &rate, size_t sequence_length,rng &r )
      :sequenceMutationModel(rate,sequence_length,r) {};
  };

  class F84mutation:public sequenceMutationModel {
  };


  class F81mutation:public sequenceMutationModel {
  };
}


#endif
