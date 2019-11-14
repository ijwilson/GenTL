#include "flowGram.h"
#include "SffInfo.h"
#include "newio.h"
#include "flowSeq.h"
#include <sstream>
#include <cstdlib>  // for system
#include <cstdio>  // tmpnam  
#include <cmath> 
#include <numeric>  // for accumulate
#include <cassert>

flowGram::flowGram(const std::vector<std::string> &lines)
{
  std::vector<std::string> tok;
  Tokenize(lines[0],tok," >=");
  name=tok[0];
  region=atoi(tok[4].c_str());
  run=tok[6];

  std::vector<std::string> tokb;
  Tokenize(tok[2],tokb,"_");
  x=atoi(tokb[0].c_str());
  y=atoi(tokb[1].c_str());
  
  for (size_t ii=1;ii<lines.size();ii++) {
    tok.clear();
    Tokenize(lines[ii],tok,", ");
    //   printvector(tok);
    for (size_t jj=0;jj<tok.size();jj+=2) {
      Bases.push_back(tok[jj][0]);
      intensity.push_back(atof(tok[jj+1].c_str()));
    }
  }
}

/**************************************************************************************************/
std::vector<std::vector<flowGram> > 
SplitbyPrimers(const std::vector<flowGram> xx, const std::vector<std::string> &primers
	       , const std::vector<double> &zp,
	       const std::vector<double> &pp, double minBayes) 
{
  size_t np=primers.size();
  std::vector<std::vector<flowGram> > aa(np+1); // the results
  std::vector<flowSeq> fsprimers;
  for (size_t ii=0;ii<np;ii++)  fsprimers.push_back(flowSeq(primers[ii]));

  for (size_t ii=0;ii<xx.size();ii++) {  // loop through the flowgrams
    std::vector<double> lp(np);          
    for (size_t jj=0;jj<np;jj++)  lp[jj] = fsprimers[jj].lprobsimple(xx[ii].intensity,zp,pp);
    // find maximum and rescale 
    size_t max_index= std::distance(lp.begin(),std::max_element(lp.begin(),lp.end()));
    double max_value=lp[max_index];
    for (size_t jj=0;jj<np;jj++) lp[jj] = exp(lp[jj]-max_value);
    double sum=std::accumulate(lp.begin(),lp.end(),0.0);
    double denominator = sum-1.0;
    if ((fabs(denominator)<1E-40)||(1.0/denominator>minBayes)) {
      aa[max_index].push_back(xx[ii]);
    } else {
      aa[np].push_back(xx[ii]);
    }
  }
  return aa;
}
/**************************************************************************************************
std::ostream &flowGram::print(std::ostream &o,const std::string &gap) const {
  for (size_t ii=0;ii<Bases.size()-1;ii++) {
    o << Bases[ii] << ":" << intensity[ii] << gap;
  }
  o << Bases.back() << ":" << intensity.back();
  return o;
}
*/
