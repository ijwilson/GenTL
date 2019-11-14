
#include <iostream>
#include <sstream>
#include <cmath>
#include "gsl/gsl_randist.h"
#include "flowSeq.h"
#include "util.h" // for error
/**                                                                              */
flowSeq::flowSeq(size_t len,const std::string &Floworder)
 :copies(len),order(Floworder) {
}
/**                                                                              */
flowSeq::flowSeq(const std::string &seq,bool TrimFloatingStart, bool TrimFloatingEnd,const std::string &Name,const std::string &FlowOrder)
:order(FlowOrder),name(Name)
{
  std::string::const_iterator cc=seq.begin();
 
  if (TrimFloatingStart) {
	//  std::cerr << seq << std::endl;
	  
	  // step though until you are sure you have been through at least one
	  // start at the first G, or the second A, or the second T or the 
	char cnt[4]={1,0,0,0};   // which have had the chance to be seen before
//	std::cerr << " Dropping ";
	while(true) {
	  assert(*cc=='A'||*cc=='T'||*cc=='G'||*cc=='C');   
	  int nxt = order.find(*cc);  // what is the position of the next
	  
	  if (cnt[nxt]>0) break;      // had the change to accept this as start
	  else {
	    //std::cerr << *cc;
	    for (int jj=1;jj<=nxt;jj++) cnt[jj]++;    // or throw and wait for another		
	    // get rid of the rest of these
	    char curr=*cc;
	    while (*cc==curr) cc++;
	  }
	}
	//	std::cerr << "\nDropped :";
	//	std::cerr << std::string(seq.begin(),cc) << std::endl;
//	std::cerr << "  Starting with " << *cc << " now getting flowgram" << std::endl;
  }
  
  while (true) {
    assert(*cc=='A'||*cc=='T'||*cc=='G'||*cc=='C');   
    for (int i=0;i<4;i++) {
      int count=0;
      while (*cc==order[i]) {
	count+=1;
	if (++cc==seq.end()) break;
      }
      copies.push_back(count);
    }
    if (cc == seq.end()) break;
  }
  
  if (!TrimFloatingEnd) {
    if (length()%4!=0) 
      copies.erase(copies.begin(),copies.begin()+4*length()/4);
  }
}

flowSeq::flowSeq(const std::string &seq,const string &Name,const string &FlowOrder)
 :order(FlowOrder),name(Name)
{
   std::string::const_iterator cc=seq.begin();
    
  while (true) {
    assert(*cc=='A'||*cc=='T'||*cc=='G'||*cc=='C');   
    for (int i=0;i<4;i++) {
      int count=0;
      while (*cc==order[i]) {
	count+=1;
	if (++cc==seq.end()) break;
      }
      copies.push_back(count);
    }
    if (cc == seq.end()) break;
  }
}

/**                                                                              */
std::ostream &flowSeq::print(std::ostream &o,const std::string &sep) const
{ 
//   for (size_t ii=0;ii<tag.size();ii++) {
//     if (tag[ii]==0) {  // no copies so use lognormal
//       double xx =  gsl_ran_lognormal_pdf(flowgram[ii]+0.005,zp[0],zp[1]);
//       if (xx>0) temp += log(xx);
//       else return -1E99;
//     } else {
//       double xx;
//       if (static_cast<size_t>(tag[ii])>pp.size()) 
// 	xx = gsl_ran_gaussian_pdf(flowgram[ii]-tag[ii],pp.back());
//       else
//       	xx =  gsl_ran_gaussian_pdf(flowgram[ii]-tag[ii],pp[tag[ii]]-1);
//       if (xx>0) temp += log(xx);
//       else return -1E99;
//     }
//   }
  size_t ii;
  for (ii=0;ii<copies.size()-1;ii++) {
    o << order[ii%4] << ":" << copies[ii] << sep;
  }
  o << order[ii%4] << ":" << copies[ii];

  return o;
}
/**                                                                              */
std::ostream &operator<<(std::ostream &out, const flowSeq &a)
{
  return a.print(out,",");
}

/** convert a flowseq back into a sequence    */
std::string flowSeq::sequence() const
{
  std::ostringstream oss;
  for (size_t ii=0;ii<copies.size();ii++) {
    for (int jj=0;jj<copies[ii];jj++) oss << order[ii%4];
  }
  return oss.str();
}
/**                                                                              */
double flowSeq::lprob(const std::vector<double> &flowgram, const std::vector<double> &zp
		      , const std::vector<double> &pp,const std::vector<int> &tag) const
/** We shall assume that the alignment is correct here.  Other probabilities that sum over
 * other alignments shall have to wait               
 * The first parameter is the mean of the lognormal distribution for 0 copies
 * the rest are the precisions                                                               */
{
  double temp=0.0;
  //  std::cerr << flowgram.size() << " " << copies.size() << " " << tag.size() << std::endl;
  if (flowgram.size() != copies.size()+tag.size()) {
    std::cerr << " stop here";
  }
  assert(flowgram.size() == copies.size()+tag.size());

  std::vector<int>::const_iterator cc=tag.begin();
  std::vector<double>::const_iterator ff=flowgram.begin();

  while (ff!=flowgram.end()) {
    if (cc==tag.end()) cc=copies.begin();

   if (*cc==0) {  // no copies so use lognormal
      double xx =  gsl_ran_lognormal_pdf(*ff+0.005,zp[0],zp[1]);
      if (xx>0) temp += log(xx);
      else return -1E99;
    } else {
      double xx;
      if (static_cast<size_t>(*cc)>pp.size()) 
	xx = gsl_ran_gaussian_pdf(*ff-*cc,pp.back());
      else
      	xx =  gsl_ran_gaussian_pdf(*ff-*cc,pp[*cc-1]);
      if (xx>0) temp += log(xx);
      else return -1E99;
    }
   ff++;
   cc++;
  }
  return temp;
}
/** proability from the start and just covering the length of the flowsequence           */
double flowSeq::lprobsimple(const std::vector<double> &flowgram, const std::vector<double> &zp
		      , const std::vector<double> &pp) const
/** We shall assume that the alignment is correct here.  Other probabilities that sum over
 * other alignments shall have to wait               
 * The first parameter is the mean of the lognormal distribution for 0 copies
 * the rest are the precisions                                                               */
{
  double temp=0.0;
 
  for (size_t ii=0;ii<length();ii++) {
     if (copies[ii]==0) {  // no copies so use lognormal
       double xx =  gsl_ran_lognormal_pdf(flowgram[ii]+0.005,zp[0],zp[1]);
       if (xx>0) temp += log(xx);
       else return -1E99;
     } else {
       double xx;
       if (static_cast<size_t>(copies[ii])>pp.size()) 
	 xx = gsl_ran_gaussian_pdf(flowgram[ii]-copies[ii],pp.back());
       else
	 xx =  gsl_ran_gaussian_pdf(flowgram[ii]-copies[ii],pp[copies[ii]-1]);
       if (xx>0) temp += log(xx);
       else return -1E99;
     }
   }
   return temp;
}


//extern int maxCopies;

void flowSeq::checkseq(std::string message,int maxCopies) const
{
  assert(maxCopies>0);
  int count=0;
  for  (size_t ii=0;ii<length();ii++) {
    if (copies[ii]<0) 
      error("problem with flowseq, should not have below zero entries, have",copies[ii],"\n");
    if (copies[ii]>maxCopies) 
      warning("maximum number of copies greater than ",maxCopies,", have ",copies[ii],"\n");
    if (copies[ii]==0) {
      count++;
      if (count==3) 
	error("should not have three consecutive zeros ",message ,"\n") ;
    } else count=0;
  }	
}
