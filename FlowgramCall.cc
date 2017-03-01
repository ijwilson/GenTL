#include "FlowgramCall.H"
#include "gsl/gsl_randist.h"
#include "gsl_rand.H"
#include "util.H"

#include "tnt/tnt.h"

#include <cmath>

flowSeq maxCall(const std::vector<double> &flowgram,
		const double normgammpar[2],
		const std::vector<double> &sigma, 
		const std::string &FlowOrder
		)
{
  flowSeq fs(flowgram.size(),FlowOrder);
  size_t maxsig = sigma.size(); // but these are from 0 to mx-1
  for (size_t ii=0;ii<flowgram.size();ii++) {
    size_t mx = static_cast<int>(flowgram[ii]+0.5)*2+2;
    std::vector<double> p(mx+1);    // vector of probabilities

    p[0] = gsl_ran_lognormal_pdf(flowgram[ii]+0.005,normgammpar[0],normgammpar[1]);
    if (mx<maxsig) {
      for (size_t jj=1;jj<=mx;jj++) 
        p[jj] =  gsl_ran_gaussian_pdf (flowgram[ii]-static_cast<double>(jj),sigma[jj-1]);
    } else {
      for (size_t jj=1;jj<=mx;jj++) {
        double s=(jj>=maxsig)?sigma[maxsig-1]:sigma[jj-1];
        p[jj] =  gsl_ran_gaussian_pdf (flowgram[ii]-static_cast<double>(jj),s);
      }
    }
  
    fs[ii] = std::distance(p.begin(),std::max_element(p.begin(),p.end()));    
    if (flowgram[ii]==0.0&&fs[ii]==1) {
      printvector(std::cerr,p);
      std::cerr << fs[ii] << std::endl;
    }
  }
  // now I need to check that this is viable - check that there are no sequences of three zeros
  // and if so recall one of them as a one
  
  int count=0;
  for (size_t i=0;i<flowgram.size();i++) {
    if (fs[i]==0) ++count;
    else count=0;
    if (count==3) {
      std::vector<double> p(3);
      assert(fs[i-2]==0&&fs[i-1]==0&&fs[i]==0);
      p[0] = gsl_ran_gaussian_pdf (flowgram[i-2]-1.0,sigma[0]);
      p[1] = gsl_ran_gaussian_pdf (flowgram[i-1]-1.0,sigma[0]);
      p[2] = gsl_ran_gaussian_pdf (flowgram[i]-1.0,sigma[0]);
      int wh = std::distance(p.begin(),std::max_element(p.begin(),p.end()));  
      //    std::cerr <<   flowgram[i-2] << " " << flowgram[i-1] << " " << flowgram[i] << " " << wh << std::endl;
      //printvector(p);
      assert(wh>=0&&wh<=2);
      fs[i-2+wh] = 1;
      if (flowgram[i-2+wh]==0.0) {
        std::cerr <<  "=========================================================\n"
                  << flowgram[i-2] << " " << flowgram[i-1] << " " << flowgram[i] << " " << wh << std::endl;
        printvector(std::cerr,p);
      }
      i+=wh-2;  // and restart with zero count
      count=0;
    }
  } 
  return fs;
}

flowSeq randomCall(const std::vector<double> &flowgram,
                   const std::vector<double> &normgammpar,
                   const std::vector<double> &sigma, 
                   rng &r,
		   int maxCopy,
                   const std::vector<int> &InitialTag,
		   const std::string &FlowOrder
                   )
{
  size_t tagsize = InitialTag.size();
  flowSeq fs(flowgram.size()-tagsize,FlowOrder);

  size_t maxsig = sigma.size(); // but these are from 0 to mx-1
  for (size_t ii=tagsize;ii<flowgram.size();ii++) {
    size_t mx = static_cast<int>(flowgram[ii]+0.5)*2+2;
    std::vector<double> p(mx+1);    // vector of probabilities
    p[0] = gsl_ran_lognormal_pdf(flowgram[ii]+0.005,normgammpar[0],normgammpar[1]);
    if (mx<maxsig) {
      for (size_t jj=1;jj<=mx;jj++) 
        p[jj] =  gsl_ran_gaussian_pdf (flowgram[ii]-static_cast<double>(jj),sigma[jj-1]);
    } else {
      for (size_t jj=1;jj<=mx;jj++) {
        double s=(jj>=maxsig)?sigma[maxsig-1]:sigma[jj-1];
        p[jj] =  gsl_ran_gaussian_pdf (flowgram[ii]-static_cast<double>(jj),s);
      }
    }
    fs[ii-tagsize] = gen_from_p(p,r);
    if (fs[ii-tagsize]>maxCopy) fs[ii-tagsize]=maxCopy;
  }
  // now I need to check that this is viable - check that there are no sequences of three zeros
  // and if so recall one of them as a one
  // there has been a filter so this should be fine!!!
  
  int count=0;
  for (size_t i=tagsize;i<flowgram.size();i++) {
    if (fs[i-tagsize]==0) {
      ++count;
      if (count==3) {
	std::vector<double> p(3);
	p[0] = gsl_ran_gaussian_pdf (flowgram[i-2]-1.0,sigma[0]);
	p[1] = gsl_ran_gaussian_pdf (flowgram[i-1]-1.0,sigma[0]);
	p[2] = gsl_ran_gaussian_pdf (flowgram[i]-1.0,sigma[0]);
	int wh = gen_from_p(p,r);
	fs[i-2+wh-tagsize] = 1;
	i+=wh-2;  // and restart with zero count
	count=0;
      }
    } else count=0;
  }
#ifndef NDEBUG
  fs.checkseq("After random Call",maxCopy);
#endif
  return fs;
}


/** call a flowseq from a vector of flowgrams indexed by wh                        */
flowSeq randomCallv(const std::vector< std::vector<double> > &flowgramv,
		    const std::vector<int> &wh,
		    const std::vector<double> &normgammpar,
		    const std::vector<double> &sigma, 
		    rng &r,
		    int maxCopy,
		    const std::vector<int> &InitialTag,
       
		     const std::string &FlowOrder
		    )
{
  size_t len = flowgramv[wh[0]].size();
  size_t tagsize=InitialTag.size();
#ifndef NDEBUG
  for (size_t jj=1;jj<wh.size();jj++)  assert( flowgramv[wh[jj]].size()==len);
#endif
  flowSeq fs(len-tagsize,FlowOrder);
  size_t maxsig = sigma.size(); // but these are from 0 to mx-1
  for (size_t ii=tagsize;ii<len;ii++) {
    size_t mx = static_cast<int>(flowgramv[wh[0]][ii]+0.5)*2+2;
    std::vector<double> p(mx+1,0.0);    // vector of log probabilities
    TNT::Array2D<double> tmpp(wh.size(),mx+1);

    for (size_t jj=0;jj<wh.size();jj++) {
      double xx = gsl_ran_lognormal_pdf(flowgramv[wh[jj]][ii]+0.005,normgammpar[0],normgammpar[1]);
      tmpp[jj][0] = xx;
      if (xx>0) {
	p[0] += log(xx);
      } else {
	p[0] = -1E99;
      }
      if (mx<maxsig) {
	for (size_t kk=1;kk<=mx;kk++) {
	  double xx =  gsl_ran_gaussian_pdf(flowgramv[wh[jj]][ii]-static_cast<double>(kk),sigma[kk-1]);
	  tmpp[jj][kk] = xx;
	  if (xx>0) p[kk] +=  log(xx);
	  else p[kk] = -1E99;
	}
      } else {
	for (size_t kk=1;kk<=mx;kk++) {
	  double s=(kk>=maxsig)?sigma[maxsig-1]:sigma[kk-1];
	  double xx = gsl_ran_gaussian_pdf(flowgramv[wh[jj]][ii]-static_cast<double>(kk),s);
	  tmpp[jj][kk] =xx ;
	   if (xx>0) p[kk] +=  log(xx);
	  else p[kk] = -1E99;
	}
      }
    }
    double mxp=*std::max_element(p.begin(),p.end());
    for (size_t iii=0;iii<p.size();iii++) p[iii]=exp(p[iii]-mxp);
    //  printvector(p);
    if (std::accumulate(p.begin(),p.end(),0.0)<=0.0) {
      std::cerr << "problem here\n";
      std::ofstream outerr("out.error");
      outerr  << tmpp;
    }
    int ff =  gen_from_p(p,r);
    if (ff>maxCopy) ff=maxCopy;
    fs[ii-tagsize] = ff;
 
  }
  // now I need to check that this is viable - check that there are no sequences of three zeros
  // and if so recall one of them as a one
  // there has been a filter so this should be fine!!!
  
  int count=0;
  for (size_t i=tagsize;i<len;i++) {
    if (fs[i-tagsize]==0) {
      ++count;
      if (count==3) {
	std::vector<double> p(3,0.0);
	for (size_t jj=0;jj<wh.size();jj++) {
	  double x =  gsl_ran_gaussian_pdf (flowgramv[wh[jj]][i-2]-1.0,sigma[0]);
	  if (x>0.0) p[0] += log(x);
	  else p[0]=-1E99;
	  x = gsl_ran_gaussian_pdf (flowgramv[wh[jj]][i-1]-1.0,sigma[0]);
	  if (x>0.0) p[1] += log(x);
	  else p[1]=-1E99;
	  x   = gsl_ran_gaussian_pdf (flowgramv[wh[jj]][i]-1.0,sigma[0]);
	  if (x>0.0) p[2] += log(x);
	  else p[2]=-1E99;
	}
	double mxp=*std::max_element(p.begin(),p.end());
	for (size_t iii=0;iii<p.size();iii++) p[iii]=exp(p[iii]-mxp);
	int wh = gen_from_p(p,r);
	fs[i-2+wh-tagsize] = 1;
	i+=wh-2;  // and restart with zero count
	count=0;
      }
    } else count=0;
  }
#ifndef NDEBUG
  fs.checkseq("After random Call vector",maxCopy);
#endif
  return fs;
}
