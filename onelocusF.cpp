#include "onelocusF.h"
#include <map>
#include <stdexcept>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#ifdef USE_R
#include "R.h"
#endif

/** Started a major averhaul on 3/4/05
 * I'm pretty sure that I got this wrong by incorporating
 * the prior for g.  Get it right !!!!!!1                   */


/** Constructor                                                                    */
onelocusF::onelocusF(const std::vector<std::string> &a, const TNT::Array1D<int> &loc)
  :n(a.size()),g_(loc.dim()),order(loc.dim())
{
  location=loc.copy();
  x=convdata1(a);
  npops_=*(std::max_element(&location[0],&location[0]+location.dim()))+1;
  k=*(std::max_element(&x[0],&x[0]+x.dim()))+1;
  
  for (int i=0;i<n;i++) order[i]=i;
}
/** Constructor for use from within R                                              */
onelocusF::onelocusF( int *a,  int *loc, int len)
  :x(len,a),location(len,loc),n(len),g_(len),order(len)
{
  npops_=*(std::max_element(&location[0],&location[0]+location.dim()))+1;
  k=*(std::max_element(&x[0],&x[0]+x.dim()))+1;
  for (int i=0;i<n;i++) order[i]=i;
}
/** calculate the "prior" probability of the global vector cc conditional
 * on the given order and F                                                        */
double onelocusF::lpg(double *F)
{
  TNT::Array1D<int> local(npops_,0);
  local[location[order[0]]] =1;
  double lp=0.0;
  for (int i=1;i<n;i++) {
    int pop=location[order[i]];
    if (fabs(F[pop])>1E-9) {
      lp -= log(1.+(local[pop]-1.)*F[pop]);
		if (g_[order[i]]) {
		lp += log(1.-F[pop]);
		}  else {
			assert(local[pop]>0);
			lp += log(local[pop]*F[pop]); 
		} 
    } else assert(g_[order[i]]);//else lp+=0.0;  
    local[pop]+=1;
  } 
  return lp;
}
/*  calculate the likelihood of the data conditional on theta and g                */
double onelocusF::likelihood(double thet, bool *gg)
{
  TNT::Array2D<int> local(npops_,k,0);
  TNT::Array1D<int> global(k,0);
  TNT::Array1D<int> ss(npops_,0);
  // the first sample
  int totg=1;
  global[x[order[0]]]=1;
  local[location[order[0]]][x[order[0]]]=1;
  ss[location[order[0]]]=1;
  double lp=0.0;
  // now loop through the rest
  for (int j=1;j<n;j++) {
    int pop=location[order[j]];
    int al=x[order[j]]; 
    if (gg[order[j]]) { // global copy
      if (global[al]==0) {             // new allele  
       
	lp += log(thet)- log(totg+thet); 
      } else {
	lp += log(global[al]) - log(totg+thet);
      }
      global[al]+=1;
      totg+=1;
    } else {           // local copy
      lp += log(double(local[pop][al])/(double(ss[pop])));
    }
    local[pop][al]+=1;
    ss[pop]+=1;
  }
  return lp;
}
/*  calculate the likelihood of the data conditional on theta and 
 * globalness up to (but not including) the previous value                 
 */
double onelocusF::likelihoodb(double thet, double *F)
{
  TNT::Array2D<int> local(npops_,k,0);
  TNT::Array1D<int> global(k,0);
  TNT::Array1D<int> ss(npops_,0);
  // the first sample
  int totg=1;
  global[x[order[0]]]=1;
  local[location[order[0]]][x[order[0]]]=1;
  ss[location[order[0]]]=1;
  double lp=0.0;
  // now loop through the rest
  for (int j=1;j<n;j++) {
    int pop=location[order[j]];
    int al=x[order[j]]; 
    // we don't know if this is global or not
    double pglobal=(1.-F[pop])/(1.+(ss[pop]-1.)*F[pop]);
    if (global[al]==0) { // must be global
      lp += log(pglobal)+log(thet)- log(totg+thet);
    } else {
      double plocal=local[pop][al]*F[pop]/(1.+(ss[pop]-1.)*F[pop]);
      lp+=log(plocal+pglobal*global[al]/(thet+totg));
    } 
    if (global[order[j]]) { // global copy
      global[al]+=1;
      totg+=1;
    } 
    local[pop][al]+=1;
    ss[pop]+=1;
  }
  return lp;
}

/******************************************************************/
/******************************************************************/
/******************************************************************/
/**
   Convert the data from a vector of vector strings to a 
   TNT::Array2D<int> with each row a sample            
*/
TNT::Array1D<int> convdata1(const std::vector<std::string> &a)
{
  TNT::Array1D<int> x(a.size());
  
  std::map<std::string,int> lmap;
  int nv=0;
  for (size_t j=0;j<a.size();j++) {
    if (lmap.find(a[j])==lmap.end()) {
      lmap[a[j]]=nv++;
      x[j]=nv-1;
    } else {
      x[j]=lmap[a[j]];
    }
  }
  return x;
}
/** what is the likelihood of the data given the 
 * global status 
 * which does not depend on F?!
 */
double onelocusF::w(double theta)
{
  TNT::Array2D<int> local(npops_,k,0);
  TNT::Array1D<int> global(k,0);
  TNT::Array1D<int> ss(npops_,0);
  // the first sample
  int totg=1;
  global[x[order[0]]]=1;
  local[location[order[0]]][x[order[0]]]=1;
  ss[location[order[0]]]=1;
  double lp=0.0;
  // now loop through the rest
  for (int j=1;j<n;j++) {
    int pop=location[order[j]];
    int al=x[order[j]];
    if (g_[order[j]]) { // global copy
      if (global[al]==0) {             // new allele  
	lp += log(theta)- log(totg+theta); 
      } else {
	lp += log(global[al]) - log(totg+theta);
      }
      global[al]+=1;
      totg+=1;
    } else {           // local copy
      lp += log(double(local[pop][al])/(double(ss[pop])));
    }
    local[pop][al]+=1;
    ss[pop]+=1;
  }
  return lp;
}
double onelocusF::likelihoodF0(double theta)
{
  TNT::Array1D<int> global(k,0);
  // the first sample
  global[x[order[0]]]=1;
  double lp=0.0;
  int totg=1;
  // now loop through the rest
  for (int j=1;j<n;j++) {
    int al=x[order[j]];
    if (global[al]==0) {             // new allele  
      lp += log(theta)- log(totg+theta); 
    } else {
      lp += log(global[al]) - log(totg+theta);
    }
    global[al]+=1;
    totg+=1;
  }
  return lp;
}
double onelocusF::resampleg_b(double *F,double theta,rng &r)
{
  TNT::Array2D<int> local(npops(),k,0);
  TNT::Array1D<int> global(k,0);
  TNT::Array1D<int> ss(npops(),0);
  // the first sample
  int totg=1;
  global[x[order[0]]]=1;
  local[location[order[0]]][x[order[0]]]=1;
  ss[location[order[0]]]=1;
  g_[order[0]]=true;
  double lp=0.0;
  // now loop through the rest
  for (int j=1;j<n;j++) {
    int pop=location[order[j]];
    int allele=x[order[j]];
    if (fabs(F[pop])>1E-9) {
      if (local[pop][allele]==0) {
	g_[order[j]]=true;  // global copy
	global[allele]+=1;
	totg+=1;
      } else {
	double pglobal= (1.-F[pop])*double(global[allele])/(double(totg)+theta);
	double plocal=double(local[pop][allele])*F[pop];
	double p=plocal/(plocal+pglobal);
	if (r()<p) {
	  g_[order[j]]=false;
	  lp += log(p);
	} else {
	  g_[order[j]]=true;
	  lp += log(1.-p);
	  global[allele]+=1;
	  totg+=1;
	}
      }
    } else {
      g_[order[j]]=true;
      global[allele]+=1;
      totg+=1;
    }
    local[pop][allele]+=1;
    ss[pop]+=1;
  }
  return lp;
}
double onelocusF::lpsampleg_b(double *F,double theta)
{
 TNT::Array2D<int> local(npops_,k,0);
  TNT::Array1D<int> global(k,0);
  TNT::Array1D<int> ss(npops_,0);
  // the first sample
  int totg=1;
  global[x[order[0]]]=1;
  local[location[order[0]]][x[order[0]]]=1;
  ss[location[order[0]]]=1;
  double lp=0.0;
  // now loop through the rest
  for (int j=1;j<n;j++) {
    int pop=location[order[j]];
    int al=x[order[j]];
   if (fabs(F[pop])>1E-9) {
     if (local[pop][al]==0) {
       global[al]+=1;
       totg+=1;
     } else {
       double pglobal= (1.-F[pop])*double(global[al])/(double(totg)+theta);
       double plocal=double(local[pop][al])*F[pop];
       double p=plocal/(plocal+pglobal);
       if (!g_[order[j]]) {
	 lp += log(p);
       } else {
	 lp += log(1.-p);
	 global[al]+=1;
	 totg+=1;
       }
     }
   } else {
     if (!g_[order[j]]) {
       throw std::runtime_error("Error - local with F=0.0\n");
     }
     global[al]+=1;
     totg+=1;
   }
   local[pop][al]+=1;
   ss[pop]+=1;
  }
  
  return lp;
}

  /** Now for the straight calculation of the likelihood we need to first
   *  sample nperm missing data values
   */
double onelocusF::likelihood( TNT::Array1D<double> &F, double theta, rng &r, int nperms) 
{
  assert(F.dim()==npops());
  double ret=0.0;

  for (int i=0;i< nperms;i++) {
    perm(r);
    resampleg_b(F,theta,r);
    ret += exp(w(theta));
  }
  return log(ret)-log(static_cast<double>(nperms));
}
/*******************************************************************************8
 ******************  Old versions *********************************************8
 *****************************************************************************/
/** 
    sample g conditional on F, allowing only those arrangements that are 
    allowable                 

double onelocusF::resampleg(double *F, rng &r)
{
   TNT::Array2D<int> local(npops_,k,0);
   TNT::Array1D<int> ss(npops_,0);
   // the first sample
   g_[order[0]]=true;
   local[location[order[0]]][x[order[0]]]=1;
   ss[location[order[0]]]=1;
   double lp=0.0;
   // now loop through the rest
   for (int j=1;j<n;j++) {
    int pop=location[order[j]];
    int al=x[order[j]];
    if (fabs(F[pop])>1E-9) {
      if (local[pop][al]==0) {  // got to be global
	g_[order[j]]=true;
      } else {
	double p1=F[pop]*ss[pop]/(1.+(ss[pop]-1.)*F[pop]); // probability local
	if (r()<p1) {
	  g_[order[j]]=false;
	  lp += log(p1);
	} else {
	  g_[order[j]]=true;
	  lp += log(1.-p1);
	}
      }
    } else {
      g_[order[j]]=true;
      // lp += 0.0;
    }
    local[pop][al]+=1;
    ss[pop]+=1;
   }
   return lp;
}

*/
/** 
    what is the probability of sampling global vector ?   
 
double onelocusF::lpsampleg(double *F) 
{
  TNT::Array2D<int> local(npops_,k,0);
  TNT::Array1D<int> ss(npops_,0);
  // the first sample
  assert(g_[order[0]]==true);
  local[location[order[0]]][x[order[0]]]=1;
  ss[location[order[0]]]=1;
  double lp=0.0;
  // now loop through the rest
  for (int j=1;j<n;j++) {
    int pop=location[order[j]];
    int al=x[order[j]];
    if (fabs(F[pop])>1E-9) {
      if (local[pop][al]>0) {
	double p1=F[pop]*ss[pop]/(1.+(ss[pop]-1.)*F[pop]); // probability local
	if (g_[order[j]]) { // global
	  lp += log(1.-p1);
	} else {
	  lp += log(p1);
	}
      } // else lp += 0.0
    } else assert(g_[order[j]]);
    local[pop][al]+=1;
    ss[pop]+=1;
  }
  return lp;
}
*/
