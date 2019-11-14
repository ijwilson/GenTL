#include "onelocusMCMC.h"
#include "onelocusF.h"
#include <map>
       
#ifndef USE_R
onelocusMCMC::onelocusMCMC(const std::vector<std::string> &a, const TNT::Array1D<int> &loc
	   , int seed,ctsdistribution *th,ctsdistribution *f,int main, bool oneF)
  :thetaprior(th),Fprior(f),mainland(main),onef(oneF),n(a.size())
  ,g_(loc.dim()),order(loc.dim()),r(seed)
{
  location=loc.copy();
  x=convdata1(a);
  kickstart();
}
#endif

onelocusMCMC::onelocusMCMC(int *xx, int *loc, int len, ctsdistribution *th
			   ,ctsdistribution *f)
  :thetaprior(th),Fprior(f),x(len,xx),location(len,loc)
  ,mainland(-1),onef(true),n(len),g_(len),order(len)
#ifndef USE_R
  ,r(1) 
#endif
{
  kickstart();
}


/** Put the startup code for both constructors into a single function */
void onelocusMCMC::kickstart() {

  npops=*(std::max_element(&location[0],&location[0]+location.dim()))+1;
  k=*(std::max_element(&x[0],&x[0]+x.dim()))+1;
  for (int i=0;i<n;i++) order[i]=i;
  permute(&order[0],n,r);
  F=TNT::Array1D<double>(npops,0.1);
  if (Fprior->proper()) {
    for (int i=0;i<npops;i++) {
      F[i]=Fprior->sample(r);
      if (F[i]<=0.0||F[i]>=1.0) F[i]=r()*0.1;
    }
  }

  if (onef) for (int i=1;i<npops;i++) F[i]=F[0];
 if (thetaprior->proper()) { 
   theta=thetaprior->sample(r);
   if (theta<0.0) theta=-theta;
 } else theta=2.0;
  sampleg(g_);
  llike=likelihood();
  lprobg=lpg(F,g_);
}


/**
   P(g|order,F)
*/
double onelocusMCMC::lpg(double *F, bool *gg)
{
  TNT::Array1D<int> local(npops,0);
  local[location[order[0]]] =1;
  double lp=0.0;
  for (int i=1;i<n;i++) {
    int l=location[order[i]];
    if (l!=mainland) {
      lp -= log(1.+(local[l]-1.)*F[l]);
      if (gg[order[i]]) {
	lp += log(1.-F[l]);
      }  else {
	assert(local[l]>0);
	lp += log(local[l]*F[l]); 
      } 
      local[l]+=1;
    } 
  }
  return lp;
}

bool onelocusMCMC::updateF(double tune)
{
  int wh=r.rint(npops);
  double oldF=F[wh]; 
  F[wh]+=r.normal()*tune;
  while (F[wh]<0.0||F[wh]>1.0) {
    if (F[wh]<0.0) F[wh]=-F[wh];
    if (F[wh]>1.0) F[wh]=2.-F[wh];
  }  
  if (onef) {
    for (int i=0;i<npops;i++) F[i]=F[wh];
  }
  // double newll=likelihood();
  double newlpg=lpg(F,g_);
  double x=exp(newlpg-lprobg
	       +Fprior->log_pdf(F[wh])-Fprior->log_pdf(oldF));
  if (x>1.||r()<x) {
    // llike=newll;
    lprobg=newlpg;
    return true;
  } else {
    F[wh]=oldF;
    if (onef) for (int i=0;i<npops;i++) F[i]=oldF; 
    return false;
  }
}
/** Try to update g by swapping a pair                */
bool onelocusMCMC::updateg()
{
  double lpback=lpsampleg(g_);
  
  std::pair<int,int>wh =r.sample2int(n);
  std::swap(order[wh.first],order[wh.second]);
  TNT::Array1D<bool> gg(n);
  double lpsamp=sampleg(gg);
  assert(fabs(lpsamp-lpsampleg(gg))<1E-10);
  double newlpg=lpg(F,gg);
  double newll=likelihooda(theta,gg);
  double x=exp(newll-llike
	       +newlpg-lprobg
	       -lpsamp+lpback);
  if (x>1.||r()<x) {
    llike=newll;
    lprobg=newlpg;
    g_.inject(gg);
    return true;
  } else {
    std::swap(order[wh.first],order[wh.second]);
    return false;
  } 
}
bool onelocusMCMC::updateg2()
{
  double lpback=lpsampleg(g_);
 
  TNT::Array1D<bool> gg(n);
  double lpsamp=sampleg(gg);
  assert(fabs(lpsamp-lpsampleg(gg))<1E-10);
  double newlpg=lpg(F,gg);
  double newll=likelihooda(theta,gg);

  double x=exp(newll-llike
	       +newlpg-lprobg
	       -lpsamp+lpback);
  //   std::cout << x << " probabilities newll "<<newll <<" llike "<< llike
  //	    << "    newlpg " << newlpg << " lprobg " << lprobg 
  //    << " lpsamp " << lpsamp << " lpback " << lpback << std::endl;
  if (x>1.||r()<x) {
    llike=newll;
    lprobg=newlpg;
    g_.inject(gg);
    return true;
  } else {
    return false;
  } 
}

bool onelocusMCMC::updateg3()   // independence sampler
{
  double lpback=lpsampleg(g_);
  TNT::Array1D<int> oldorder=order.copy();
  permute(&order[0],order.dim(),r);
  TNT::Array1D<bool> gg(n);
  double lpsamp=sampleg(gg);
  //  assert(fabs(lpsamp-lpsampleg(gg))<1E-10);
  double newlpg=lpg(F,gg);
  double newll=likelihooda(theta,gg);

  double x=exp(newll-llike
	       +newlpg-lprobg
 	       -lpsamp+lpback);
  //      std::cout << x << " probabilities newll "<<newll <<" llike "<< llike
  //	    << "    newlpg " << newlpg << " lprobg " << lprobg 
  //    << " lpsamp " << lpsamp << " lpback " << lpback << std::endl;
  if (x>1.||r()<x) {
    llike=newll;
    lprobg=newlpg;
    g_.inject(gg);
    return true;
  } else {
    order.inject(oldorder);
    return false;
  } 
}
/** Attempt an update of theta */
bool onelocusMCMC::updatetheta(double tune)
{
  double newtheta=theta+tune*r.normal();
  if (newtheta<0.0) newtheta=-newtheta;
  double newll=likelihooda(newtheta,g_);
  double x=exp(newll-llike
	       +thetaprior->log_pdf(newtheta)-thetaprior->log_pdf(theta));
  if (x>1||r()<x) {
    theta=newtheta;
    llike=newll;
    return true;
  } else {
    return false;
  } 
}
/** Check that everything matches still ...  */
bool onelocusMCMC::check()
{
  if (fabs(llike-likelihooda(theta,g_))>1E-8)
    std::cerr << "error with llike" << std::endl;
  if (fabs(lprobg-lpg(F,g_))>1E-8) 
    std::cerr << "error with lprobg" << std::endl;
  return true;
}
/** Calculate the likelihood of a particular arrangement              */
double onelocusMCMC::calc_likelihood(double *FF,double thet, bool *gg)
{
  TNT::Array2D<int> local(npops,k,0);
  TNT::Array1D<int> global(k,0);
  TNT::Array1D<int> ss(npops,0);
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
    if (pop!=mainland) {
      if (gg[order[j]]) { // global copy
	lp +=  log(1.-FF[pop]) - log(1.+(ss[pop]-1.)*FF[pop]);
	if (global[al]==0) {             // new allele  
	  lp += log(thet)- log(totg+thet); 
	} else {
	  lp += log(global[al]) - log(totg+thet);
	}
	//  std::cerr << global << "and totg = " << totg << std::endl;
	assert(std::accumulate(&global[0],&global[0]+global.dim(),0)==totg);
	global[al]+=1;
	totg+=1;
      } else {           // local copy
	lp += log(FF[pop]) + log(local[pop][al]) - log(1.+(ss[pop]-1.)*FF[pop]);
      }
      local[pop][al]+=1;
      ss[pop]+=1;
    } else {
      if (global[al]==0) {             // new allele  
	lp += log(thet)- log(totg+thet); 
      } else {
	lp += log(global[al]) - log(totg+thet);
      }
      global[al]+=1;
      totg+=1;
    }
  }
  return lp;
}

double onelocusMCMC::likelihooda(double thet, bool *gg)
{
  TNT::Array2D<int> local(npops,k,0);
  TNT::Array1D<int> global(k,0);
  TNT::Array1D<int> ss(npops,0);
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
    if (pop!=mainland) {
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
    } else {
      if (global[al]==0) {             // new allele  
	lp += log(thet)- log(totg+thet); 
      } else {
	lp += log(global[al]) - log(totg+thet);
      }
      global[al]+=1;
      totg+=1;
    }
  }
  return lp;
}
/** Sample a new set of g's proportional to P(g|y,order,F) */
double onelocusMCMC::sampleg(TNT::Array1D<bool> &gg) 
{
  TNT::Array2D<int> local(npops,k,0);
  TNT::Array1D<int> ss(npops,0);
  // the first sample
  gg[order[0]]=true;
  local[location[order[0]]][x[order[0]]]=1;
  ss[location[order[0]]]=1;
  double lp=0.0;
  // now loop through the rest
  for (int j=1;j<n;j++) {
    int pop=location[order[j]];
    int al=x[order[j]];
    if (pop!=mainland) {
      if (local[pop][al]==0) {  // got to be global
	  gg[order[j]]=true;
      } else {
	double p1=F[pop]*ss[pop]/(1.+(ss[pop]-1.)*F[pop]); // probability local
	if (r()<p1) {
	  gg[order[j]]=false;
	  lp += log(p1);
	} else {
	  lp += log(1.-p1);
	  gg[order[j]]=true;
	}
      }
      local[pop][al]+=1;
      ss[pop]+=1;
    } else {
      gg[order[j]]=true;
    }
  }
  return lp;
}
/** 
    what is the probability of sampling global vector gg?   
 */
double onelocusMCMC::lpsampleg(TNT::Array1D<bool> &gg) 
{
  TNT::Array2D<int> local(npops,k,0);
  TNT::Array1D<int> ss(npops,0);
  // the first sample
  assert(gg[order[0]]==true);
  local[location[order[0]]][x[order[0]]]=1;
  ss[location[order[0]]]=1;
  double lp=0.0;
  // now loop through the rest
  for (int j=1;j<n;j++) {
    int pop=location[order[j]];
    int al=x[order[j]];
    if (pop!=mainland) {
      if (local[pop][al]>0) {
	double p1=F[pop]*ss[pop]/(1.+(ss[pop]-1.)*F[pop]); // probability local
	if (gg[order[j]]) { // global
	  lp += log(1.-p1);
	} else {
	  lp += log(p1);
	}
      } // else lp += 0.0
      local[pop][al]+=1;
      ss[pop]+=1;
    } 
  }
  return lp;
}
