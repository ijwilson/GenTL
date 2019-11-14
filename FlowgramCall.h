#ifndef FLOWGTRAMCALL_H__
#define  FLOWGTRAMCALL_H__
#include "flowSeq.h"
#include "gsl_rand.h"



flowSeq randomCall(const std::vector<double> &flowgram,
                   const std::vector<double> &normgammpar,
                   const std::vector<double> &sigma,
                   rng &r,
		   int maxCopies,
		   const std::vector<int> &InitialTag,
                   const std::string &FlowOrder="TACG"
                   );

flowSeq maxCall(const std::vector<double> &flowgram,
                   const std::vector<double> &normgammpar,
                   const std::vector<double> &sigma, 
                   const std::string &FlowOrder="TACG"
                );

/** call a flowseq form a vector of flowgrams indexed by wh                        */
flowSeq randomCallv(const std::vector< std::vector<double> > &flowgramv,
		    const std::vector<int> &wh,
		    const std::vector<double> &normgammpar,
		    const std::vector<double> &sigma, 
		    rng &r,
		    int maxCopies,
		   const std::vector<int> &InitialTag,
		    const std::string &FlowOrder="TACG"
		    );
#endif
