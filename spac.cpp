#include "spac.H"
#include "utilcpp/newio.h"
#include <iomanip>
/** operator for an mcmc analysis                                     */
double subdivpac::operator()(double rho, const TNT::Array1D<double> &F, const std::pair<int,int> &changeorder) 
{
  if (changeorder.first!=changeorder.second) {
    int *tmp = h[changeorder.first];
    h[changeorder.first]=h[changeorder.second];
    h[changeorder.second]=tmp;
    int l=current_pop[changeorder.first];
    current_pop[changeorder.first]=current_pop[changeorder.second];
    current_pop[changeorder.second]=l;
  }
  return calc(rho,F,path_,global_);      
}
/** Calculate the likelihood for the data conditional on the global/local array.  
    Theta here is generally thetahat.
 */
double likelihood(double theta)
{

}




/** calculate the likelihood for one rearrangement of the haplotypes 
 * which is assumed to have already happened.    We also have to sample 
 * (as we go) whether or not each position is in the "local" or global population */
double subdivpac::calc(double rho, const TNT::Array1D<double> & F, TNT::Array2D<int> &path
		       , TNT::Array2D<bool> &global,bool resamplepath)
{
  // path holds the path through the haplotypes, we might as well keep it 
  // as it costs little.  Global holds the copy type-local or from the global pool 
  // alpha[x,s] defined as the probability of sampling the first s sites of
  // your haplotype, and of sampling haplotype x. The first k hold the 
  // probabilites of a local copy of the site, the last k copying from the global pool.
  // Some (or most) of these will be zero (zero locally and globally if they are 
  // local to  another population, zero globally if they are local to the current
  // population, only non zero for both if they are both local and in the global pool)
  TNT::Array2D<double> alpha(sites,2*n-2);
  std::vector<int> popcount(npops_,0);     // the count in each of the populations
  TNT::Array1D<int> globalcount(sites,1); // the count in the "global" pool
  
  popcount[current_pop[0]] +=1;           // add one to the local count
  for (int s=0;s<sites;s++) {
    global[0][s]=true;
    path[0][s]=0;
  }

  double lprob=0.0;                       // the first sample has probability 1
  // now consider each of the subsequent samples in turn.
  // Because my arrays are 0 offset haplotype k is the k+1 st sample
  // so use the same k as in Li and Stephens' appendix
  for (int k=1;k<n;k++) {
    int location = current_pop[k];  // the location of the haplotype under consideration
    // what is the probability that a recombination is copied from 
    // a "local" chromosome?
    double plocalcopy=0.0;  // value if this is the first haplotype from this subpopulation
    if (popcount[location]>0) 
      plocalcopy=F[location]*double(popcount[location])/(1.+(double)(popcount[location]-1.)*F[location]);
#ifdef LOUD
    std::cout << "plocalcopy = " << plocalcopy << " local popsize " 
	      << popcount[location] <<  " location " << location << std::endl;
#endif

    // what are the match probabilities for the global population
    double gammamatchglobal=double(globalcount[0])/(double(globalcount[0])+thetahat)
      +0.5*thetahat/(double(globalcount[0])+thetahat);
    double gammamissglobal=1.- gammamatchglobal;//   0.5*thetahat/(double(globalcount[0])+thetahat);
    double gammamatchlocal;    
    if (perfectlocalmatch)
      gammamatchlocal=1.0;       
    else
      gammamatchlocal=k/(k+thetahat) + 0.5*thetahat/(k+thetahat);;

    double gammamisslocal= 1.-gammamatchlocal;         

    // we can calculate the probabilities for the first site directly
    // there are no transitions and we assume exchangeability of haplotypes 
    // loop over the haplotypes from which the site could be copied
    for (int x=0;x<k;x++) {
      if (h[x][0]==h[k][0]) {
	if (current_pop[x]==location) 
	  alpha[0][x] = plocalcopy*gammamatchlocal/double(popcount[location]);  
	else alpha[0][x]=0.0;
	if (global[x][0]) 
	  alpha[0][k+x] = (1.-plocalcopy)*gammamatchglobal/double(globalcount[0]);
	else alpha[0][k+x]=0.0;
      } else {
	if (current_pop[x]==location) 
	  alpha[0][x] = plocalcopy*gammamisslocal/double(popcount[location]);
	else alpha[0][x]=0.0;
	if (global[x][0]) 
	  alpha[0][k+x] = (1.-plocalcopy)*gammamissglobal/double(globalcount[0]);
	else alpha[0][k+x]=0.0;
      }
    }
    double logscalefactor=rescale(alpha[0],2*k); //MacDonald & Zucchini fix
    // calculate alpha for the rest of the sites
    for (int s=1;s<sites;s++) {
      // first the probability of recombination.  I am not sure 
      // about the use of localcopies here but it seems as good as anything
      // justification is that you can recombine with anything in the same
      // population, and all the globals (but these may include ones from the 
      // same population twice) which is a problem      
      double p_global_norec=exp(-rho*d[s-1]/(double(globalcount[s-1])));
      double p_local_norec=nolocalrecomb?1.0:exp(-rho*d[s-1]/double(popcount[location]));
  
      //std::cout << k << " " << p_local_norec << " " << p_global_norec << std::endl;

      gammamatchglobal=double(globalcount[s])/(double(globalcount[s])+thetahat)
	+0.5*thetahat/(double(globalcount[s])+thetahat);
      gammamissglobal=0.5*thetahat/(double(globalcount[s])+thetahat);

      assert(p_local_norec>=0.0&&p_local_norec<=1.0);
      assert(p_global_norec>=0.0&&p_global_norec<=1.0);
      // instead of 1 term in li & Stephens we need 4
      double smlocal=0.0;
      double smglobal=0.0;
      for (int ii=0;ii<k;ii++) {
	// sum up the probabilites for locals
	if (current_pop[ii]==location) smlocal+=alpha[s-1][ii];
	if (global[ii][s-1]) smglobal+=alpha[s-1][k+ii];
      }
      double smlocallocal=0.0;
      double smgloballocal=0.0;
      if (popcount[location]>0) 
	smlocallocal =smlocal*(1.-p_local_norec)*plocalcopy/double(popcount[location]);
     
      double smlocalglobal = smlocal*(1.-p_local_norec)*(1.-plocalcopy)/double(globalcount[s]);
      double smglobalglobal = smglobal*(1.-p_global_norec)*(1.-plocalcopy)/double(globalcount[s]); 
      if (popcount[location]>0) 
	smgloballocal = smglobal*(1.-p_global_norec)*plocalcopy/double(popcount[location]);

      for (int x=0;x<k;x++) {
	if (current_pop[x]==location)
	  alpha[s][x]=p_local_norec*alpha[s-1][x]+smlocallocal+smgloballocal;
	else alpha[s][x]=0.0;
	
	if (global[x][s])
	  alpha[s][k+x]=p_global_norec*alpha[s-1][k+x]+smlocalglobal+smglobalglobal; 
	else alpha[s][k+x]=0.0;

	if (h[x][s]==h[k][s]) {
	  alpha[s][x]*=gammamatchlocal;
	  alpha[s][x+k] *= gammamatchglobal;
	} else {
	  alpha[s][x]*=gammamisslocal;
	  alpha[s][x+k] *= gammamissglobal;
	}
      }
      logscalefactor+=rescale(alpha[s],2*k); // The scale factor of MacDonald and Zucchini
    }

    // calculate the prob for haplotype k by summing over all the 
    // paths to the last site
    double prob=std::accumulate(alpha[sites-1],alpha[sites-1]+2*k,0.0);
    lprob+=log(prob)+logscalefactor;
 
    if (log(prob)+logscalefactor>0.0) {
      std::cerr << "problem, added a probability > 1 " << lprob 
		<< " log(prob) = " << log(prob) << " logscalefactor = " << logscalefactor 
		<< " added " << log(prob)+logscalefactor << std::endl;
      exit(EXIT_FAILURE);
    }

    // Method as in Boys et al 2000 - the last site is "site-1".  sample from the "equilibrium"
    path[k][sites-1] = gen_from_p(alpha[sites-1],alpha[sites-1]+2*k,r);
    // need to calculate the transition probabilities here
    std::vector<double> lp(2*k);  
    for (int s=sites-2;s>=0;s--) {
      double scaled_rec_rate=globalcount[s];
      double p_global_norec=exp(-rho*d[s]/scaled_rec_rate);
      double p_local_norec;
      if (nolocalrecomb) p_local_norec=1.0;
      else p_local_norec=exp(-rho*d[s]/double(popcount[location]));
  
      int xprime=path[k][s+1];
      for (int x=0;x<k;x++) { // loop throught all the possible haplotypes
	if (xprime<k) { // current path (ahead) local
	  if (popcount[location]==0) {
	    std::cerr << "site " << s+1 << " xprime = " << xprime << std::endl;
	  }
	  if (xprime==x) {       // staying on the same path
	    lp[x] = alpha[s][x]*
	      (p_local_norec+(1.-p_local_norec)*plocalcopy/double(popcount[location]));
	  } else if (popcount[location]>0) {               // move
	    lp[x] = alpha[s][x]*
	      (1.-p_local_norec)*plocalcopy/double(popcount[location]);
	  } else lp[x]=0.0;
	  lp[k+x] = alpha[s][k+x]*
	    (1.-p_local_norec)*(1.-plocalcopy)/double(globalcount[s+1]); 
	} else {  // current path global
	  if (!global[xprime-k][s+1]) {
	    std::cout << "site " << s+1 << " xprime = " << xprime 
		      << " alpha[s+1][xprime] "<< alpha[s+1][xprime] 
		      << " alpha[s][xprime] "<< alpha[s][xprime] << std::endl;
	    exit(EXIT_FAILURE);
	  }
	  if (xprime==k+x) // copy from x as global
	    lp[k+x] = alpha[s][k+x]*
	      (p_global_norec+(1.-p_global_norec)*(1.-plocalcopy)/double(globalcount[s+1]));
	  else 
	    lp[k+x] = alpha[s][k+x]*
	      (1.-p_global_norec)*(1.-plocalcopy)/double(globalcount[s+1]);
	  if (popcount[location]>0)
	    lp[x] = alpha[s][x]*
	      (1.-p_global_norec)*plocalcopy/double(popcount[location]);
	  else lp[x]=0.0;
	}
      }
      path[k][s]=gen_from_pr(lp,2*k,r);
    }
    for (int ii=0;ii<sites;ii++) {
      if (path[k][ii]>=k) {
	global[k][ii]=true;
	globalcount[ii]+=1;
      } else global[k][ii]=false;
    }
    popcount[location]+=1;
  }
  
  
  return lprob;
}
/**********************************************************************/
/** Funciton to get the relative localness of each chromome           */
/**********************************************************************/
TNT::Array2D<double> subdivpac::localise(double rho
					 , const TNT::Array1D<double> &F
					 ,double *recombinations
					 ,bool scale, double scalefac)
{
  TNT::Array2D<double> retl(global_.dim1(),global_.dim2(),0.0);
  TNT::Array3D<double> lcl(nperm,global_.dim1(),global_.dim2(),0.0);
  TNT::Array2D<double> recs(nperm,global_.dim2()-1);
  std::vector<double> a(nperm);

  calls_ +=1;
  // use a rng to get the same set of random numbers for each run of the 
  // method - slightly slower but the time is not significant compared 
  // to the rest of the algorithm
  rng rp(seed);   
  current_pop=pop.copy();
  
  for (int i=0;i<nperm;i++) {
    std::vector<int> prms=rp.integer_permutations(n,n);
    for (int j=0;j<n;j++) {
      h[j]=st[prms[j]];
      current_pop[j]=pop[prms[j]];
    }
    a[i]=calc(rho,F,path_,global_)/scalefac; 
    for (int j=0;j<n;j++) {
      lcl[i][prms[j]][0] += global_[j][0]==false;
      for (int s=1;s<global_.dim2();s++) {
	lcl[i][prms[j]][s] += global_[j][s]==false;
	recs[i][s-1] += double(path_[j][s-1]!=path_[j][s]);
      }
    }
  }

  if (recombinations) 
    for (int i=0;i<<global_.dim2()-1;i++)
      recombinations[i]=0.0;
  
  if (!scale ) for (int i=0;i<nperm;i++) a[i]=0.0;
  double sc=*(std::max_element(a.begin(),a.end()));
  for (int i=0;i<nperm;i++) a[i] = exp(a[i]-sc);  
  double divby=std::accumulate(a.begin(),a.end(),0.0);

  for (int i=0;i<nperm;i++) { 
    for (int j=0;j<n;j++) {
      for (int s=0;s<global_.dim2();s++) {
	retl[j][s]+=a[i]*lcl[i][j][s];
      }
    }
    if (recombinations) 
      for (int s=0;s<global_.dim2()-1;s++) 
	recombinations[s]+=a[i]*recs[i][s];
  }
  for (int j=0;j<n;j++) {
    for (int s=0;s<global_.dim2();s++) {
      retl[j][s]/=divby;
    }
  }
  if (recombinations) 
    for (int s=0;s<global_.dim2()-1;s++) recombinations[s]/=divby;
  
  return retl;
}

/** do a single calculation in the original order                      */
double subdivpac::single(double rho, const TNT::Array1D<double> &F,int seed)
{
  h=st;
  current_pop=pop.copy();
  if (seed>0) r.reset(seed);
  return calc(rho,F,path_,global_);
}
/**********************************************************************/
/** Function to return a sampled path for the best ordering of 
 * haplotypes.  It returns an array with the number of the arrangement 
 * and perm contains the order (so gives the meaning of the numbers)  */
/**********************************************************************/
 TNT::Array2D<int>   subdivpac::bestpath(double rho
					 , const TNT::Array1D<double> &F
					 ,std::vector<int> &perm)
{
  TNT::Array2D<int> bpath(n,sites);

  calls_ +=1;
  // use a rng to get the same set of random numbers for each run of the 
  // method - slightly slower but the time is not significant compared 
  // to the rest of the algorithm
  rng rp(seed);   
  // and reset the ordering to its inital value
  // h=st;
  current_pop=pop.copy();

  double maxa=-1;
 
  for (int i=0;i<nperm;i++) {
    std::vector<int> rev(n);
    std::vector<int> prms=r.integer_permutations(n,n);
    for (int j=0;j<n;j++) {
      h[j]=st[prms[j]];
      current_pop[j]=pop[prms[j]];
      rev[prms[j]]=j;
    }
    double a=calc(rho,F,path_,global_);
    if (exp(a)>maxa) {
      perm=prms;
      for (int j=0;j<n;j++)
	for (int k=0;k<sites;k++) {
	  bpath[j][k]=path_[j][k];
	}
    }
  }
  return bpath;
}
/*********************************************************************************************/
/** calculate the likelihood for one rearrangement of the haplotypes 
 * which is assumed to have already happened.    We also have to sample 
 * (as we go) whether or not each position is in the "local" or global population */
TNT::Array2D<int> subdivpac::sampleSNP(double &lprob,double rho, const TNT::Array1D<double> & F
				       , TNT::Array1D<int> &pos
				       , TNT::Array1D<double> &d1, TNT::Array2D<int> &path
				       , TNT::Array2D<bool> &global,int minfreq)
{
  /** path holds the path through the haplotypes, we might as well keep it 
   * as it costs little.  Global holds the copy type-local or from the global pool 
   * alpha[x,s] defined as the probability of sampling the first s sites of
   * your haplotype, and of sampling haplotype x. The first k hold the 
   * probabilites of a local copy of the site, the last k copying from the global pool.
   * Some (or most) of these will be zero (zero locally and globally if they are 
   * local to  another population, zero globally if they are local to the current
   * population, only non zero for both if they are both local and in the global pool)  */

  TNT::Array2D<double> alpha(sites,2*n-2);
  std::vector<int> popcount(npops_,0);    // the count in each of the populations
  TNT::Array1D<int> globalcount(sites,1); // the count in the "global" pool

  // variables for the sampleSNP part of this function
  int nSNP=pos.dim();
  TNT::Array1D<int> sampleglobalcount(nSNP,1); // all initially 1
  TNT::Array2D<bool> sampleglobal(nSNP,n);
  TNT::Array2D<int> samplepath(nSNP,n);

  // The first sample
  for (int s=0;s<sites;s++) {
    global[0][s]=true;                    // the first sample is global
    path[0][s]=0;                         // mark the path
  }
  for (int SNP=0;SNP<nSNP;SNP++) {
    sampleglobal[SNP][0]=true;
    samplepath[SNP][0]=0;
  }
  popcount[current_pop[0]] +=1;           // add one to the local count
  lprob=0.0;                              // the first sample has probability 1
  // now consider each of the subsequent samples in turn.
  // Because my arrays are 0 offset haplotype k is the k+1 st sample
  // so use the same k as in Li and Stephens' appendix

  for (int k=1;k<n;k++) {
    int location = current_pop[k];  // the location of the haplotype under consideration
    // what is the probability that a recombination is copied from a "local" chromosome?
    double plocalcopy=0.0;  // value if this is the first haplotype from this subpopulation
    if (popcount[location]>0) 
      plocalcopy=F[location]*double(popcount[location])
	/(1.+(double)(popcount[location]-1.)*F[location]);

    // what are the match probabilities for the global population
    double gammamatchglobal=double(globalcount[0])/(double(globalcount[0])+thetahat)
      +0.5*thetahat/(double(globalcount[0])+thetahat);
    double gammamissglobal=1.-gammamatchglobal;
    double gammamatchlocal;       // constant values througout
  if (perfectlocalmatch)
    gammamatchlocal=1.0;
  else 
    gammamatchlocal=k/(k+thetahat) + 0.5*thetahat/(k+thetahat);;
  double gammamisslocal= 1.- gammamatchlocal;

    // we can calculate the probabilities for the first site directly
    // there are no transitions and we assume exchangeability of haplotypes 
    // loop over the haplotypes from which the site could be copied
    for (int x=0;x<k;x++) {
      if (h[x][0]==h[k][0]) {
	if (current_pop[x]==location) 
	  alpha[0][x] = plocalcopy*gammamatchlocal/double(popcount[location]);  
	else alpha[0][x]=0.0;
	if (global[x][0]) 
	  alpha[0][k+x] = (1.-plocalcopy)*gammamatchglobal/double(globalcount[0]);
	else alpha[0][k+x]=0.0;
      } else {
	if (current_pop[x]==location) 
	  alpha[0][x] = plocalcopy*gammamisslocal/double(popcount[location]);
	else alpha[0][x]=0.0;
	if (global[x][0]) 
	  alpha[0][k+x] = (1.-plocalcopy)*gammamissglobal/double(globalcount[0]);
	else alpha[0][k+x]=0.0;
      }
    }
    double logscalefactor=rescale(alpha[0],2*k);//MacDonald & Zucchini fix
    // calculate alpha for the rest of the sites
    for (int s=1;s<sites;s++) {
      // first the probability of recombination.  I am not sure 
      // about the use of localcopies here but it seems as good as anything
      // justification is that you can recombine with anything in the same
      // population, and all the globals (but these may include ones from the 
      // same population twice) which is a problem

      double scaled_rec_rate=double(globalcount[s-1]);
      double p_global_norec=exp(-rho*d[s-1]/scaled_rec_rate);
      double p_local_norec;
      if (nolocalrecomb) p_local_norec=1.0;
      else p_local_norec=exp(-rho*d[s-1]/double(popcount[location]));

      gammamatchglobal=double(globalcount[s])/(double(globalcount[s])+thetahat)
	+0.5*thetahat/(double(globalcount[s])+thetahat);
      gammamissglobal=0.5*thetahat/(double(globalcount[s])+thetahat);

      assert(p_local_norec>=0.0&&p_local_norec<=1.0);
      assert(p_global_norec>=0.0&&p_global_norec<=1.0);

      // instead of one term in li & Stephens we need 4
      double smlocal=0.0;
      double smglobal=0.0;
      for (int ii=0;ii<k;ii++) {
	// sum up the probabilites for locals
	if (current_pop[ii]==location) smlocal+=alpha[s-1][ii];
	if (global[ii][s-1]) smglobal+=alpha[s-1][k+ii];
      }
      double smlocallocal=0.0;
      double smgloballocal=0.0;
      if (popcount[location]>0) 
	smlocallocal =smlocal*(1.-p_local_norec)*plocalcopy/double(popcount[location]);
     
      double smlocalglobal = smlocal*(1.-p_local_norec)*(1.-plocalcopy)/double(globalcount[s]);
      double smglobalglobal = smglobal*(1.-p_global_norec)*(1.-plocalcopy)/double(globalcount[s]); 
      if (popcount[location]>0) 
	smgloballocal = smglobal*(1.-p_global_norec)*plocalcopy/double(popcount[location]);

      for (int x=0;x<k;x++) {
	if (current_pop[x]==location)
	  alpha[s][x]=p_local_norec*alpha[s-1][x]+smlocallocal+smgloballocal;
	else alpha[s][x]=0.0;
	
	if (global[x][s])
	  alpha[s][k+x]=p_global_norec*alpha[s-1][k+x]+smlocalglobal+smglobalglobal; 
	else alpha[s][k+x]=0.0;

	if (h[x][s]==h[k][s]) {
	  alpha[s][x]*=gammamatchlocal;
	  alpha[s][x+k] *= gammamatchglobal;
	} else {
	  alpha[s][x]*=gammamisslocal;
	  alpha[s][x+k] *= gammamissglobal;
	}
      }
      logscalefactor+=rescale(alpha[s],2*k); // The scale factor of MacDonald and Zucchini
    }

    // calculate the prob for haplotype k by summing over all the 
    // paths to the last site
    double prob=std::accumulate(alpha[sites-1],alpha[sites-1]+2*k,0.0);
    lprob+=log(prob)+logscalefactor;

    // Now for the "backwards" part of the algorithm, if required.  Method as in 
    // Boys et al 2000 - the last site is "site-1".  sample from the "equilibrium"
    path[k][sites-1] = gen_from_p(alpha[sites-1],alpha[sites-1]+2*k,r);
    // need to calculate the transition probabilities here
    std::vector<double> lp(2*k);  
    for (int s=sites-2;s>=0;s--) {
      double scaled_rec_rate=globalcount[s];
      double p_global_norec=exp(-rho*d[s]/scaled_rec_rate);
      
      double p_local_norec;
      if (nolocalrecomb) p_local_norec=1.0;
      else p_local_norec=exp(-rho*d[s]/double(popcount[location]));
      
      int xprime=path[k][s+1];
      for (int x=0;x<k;x++) { // loop throught all the possible haplotypes
	if (xprime<k) { // current path (ahead) local
	  if (popcount[location]==0) {
	    std::cerr << "site " << s+1 << " xprime = " << xprime << std::endl;
	  }
	  if (xprime==x) {       // staying on the same path
	    lp[x] = alpha[s][x]*
	      (p_local_norec+(1.-p_local_norec)*plocalcopy/double(popcount[location]));
	  } else if (popcount[location]>0) {               // move
	    lp[x] = alpha[s][x]*
	      (1.-p_local_norec)*plocalcopy/double(popcount[location]);
	  } else lp[x]=0.0;
	  lp[k+x] = alpha[s][k+x]*
	    (1.-p_local_norec)*(1.-plocalcopy)/double(globalcount[s+1]); 
	} else {  // current path global
	  if (xprime==k+x) // copy from x as global
	    lp[k+x] = alpha[s][k+x]*
	      (p_global_norec+(1.-p_global_norec)*(1.-plocalcopy)/double(globalcount[s+1]));
	  else 
	    lp[k+x] = alpha[s][k+x]*
	      (1.-p_global_norec)*(1.-plocalcopy)/double(globalcount[s+1]);
	  if (popcount[location]>0)
	    lp[x] = alpha[s][x]*
	      (1.-p_global_norec)*plocalcopy/double(popcount[location]);
	  else lp[x]=0.0;
	}
      }
      path[k][s]=gen_from_pr(lp,2*k,r);
    }
    
    for (int SNP=0;SNP<nSNP;SNP++) {
      int position=pos[SNP];
      assert(d1[SNP]<d[position]); 
      double d2=d[position]-d1[SNP];
      TNT::Array1D<double> newalpha(2*k);
  
      double scaled_rec_rate=double(globalcount[position]);
      double p_global_norec=exp(-rho*d[position]/scaled_rec_rate); 
      double p_local_norec;
      if (nolocalrecomb) p_local_norec=1.0;
      else p_local_norec=exp(-rho*d[position]/double(popcount[location]));
      //double p_local_norec=exp(-rho*d[position]/double(k));         //?!
      double plocalcopy=0.0;  // value if this is the first haplotype from this subpopulation
      if (popcount[location]>0) plocalcopy=F[location]*double(popcount[location])
	/(1.+(double)(popcount[location]-1.)*F[location]);
      // start by recalculating for the different values
      // the "second term" from Li and Stephens four terms now rather than one
      double smlocal=0.0,smglobal=0.0;
      for (int ii=0;ii<k;ii++) {// sum up the probabilites for locals
	if (current_pop[ii]==location) smlocal+=alpha[position][ii];
	if (global[ii][position]) smglobal+=alpha[position][k+ii];
      }
      double smlocallocal=0.0,smgloballocal=0.0;
      if (popcount[location]>0) 
	smlocallocal =smlocal*(1.-p_local_norec)*plocalcopy/double(popcount[location]);
    
      double smlocalglobal = smlocal*(1.-p_local_norec)*(1.-plocalcopy)/double(globalcount[position]);
      double smglobalglobal = smglobal*(1.-p_global_norec)*(1.-plocalcopy)/double(globalcount[position]); 
      if (popcount[location]>0) 
	smgloballocal = smglobal*(1.-p_global_norec)*plocalcopy/double(popcount[location]);
      for (int x=0;x<k;x++) {
	if (current_pop[x]==location)
	  newalpha[x]=p_local_norec*alpha[position][x]+smlocallocal+smgloballocal;
	else newalpha[x]=0.0;	
	if (global[x][position])
	  newalpha[k+x]=p_global_norec*alpha[position][k+x]+smlocalglobal+smglobalglobal; 
	else newalpha[k+x]=0.0;
      }

      // and now for the "backwards" part of the algorithm   
      std::vector<double> lp(2*k);
      // get the path above
      int pathabove= path[k][position+1];
      scaled_rec_rate=sampleglobalcount[SNP];
      p_global_norec=exp(-rho*d2/scaled_rec_rate);
      //p_local_norec=exp(-rho*d2/(double(k)));//??!!
      if (nolocalrecomb) p_local_norec=1.0;
      else p_local_norec=exp(-rho*d2/double(popcount[location]));
      for (int x=0;x<k;x++) {
	if (pathabove<k) { // current path (ahead) local
	  if (popcount[location]==0) {
	    std::cerr << "site " << position << " xprime = " << pathabove << std::endl;
	  }
	  if (pathabove==x) {       // staying on the same path
	    lp[x] = newalpha[x]*
	      (p_local_norec+(1.-p_local_norec)*plocalcopy/double(popcount[location]));
	  } else if (popcount[location]>0) {               // move
	    lp[x] = newalpha[x]*
	      (1.-p_local_norec)*plocalcopy/double(popcount[location]);
	  } else lp[x]=0.0;
	  lp[k+x] = newalpha[k+x]*
	    (1.-p_local_norec)*(1.-plocalcopy)/double(globalcount[position+1]); 
	} else {  // current path global
	  if (!global[pathabove-k][position+1]) {
	    std::cout << "site " << position+1 << " pathabove = " << pathabove 
		      << " alpha[s+1][pathabove] "<< alpha[position+1][pathabove] 
		      << " newalpha[pathabove] "<< newalpha[pathabove] << std::endl;
	    exit(EXIT_FAILURE);
	  }
	  if (pathabove==k+x) // copy from x as global
	    lp[k+x] = newalpha[k+x]*
	      (p_global_norec+(1.-p_global_norec)*(1.-plocalcopy)/double(globalcount[position+1]));
	  else 
	    lp[k+x] = newalpha[k+x]*
	      (1.-p_global_norec)*(1.-plocalcopy)/double(globalcount[position+1]);
	  if (popcount[location]>0)
	    lp[x] = newalpha[x]*
	      (1.-p_global_norec)*plocalcopy/double(popcount[location]);
	  else lp[x]=0.0;
	}
      }
      samplepath[SNP][k]=gen_from_pr(lp,2*k,r);
    } // end of SNP loop
    // now make all necessary adjustments !?!
    popcount[location]+=1;
    for (int ii=0;ii<sites;ii++) {
      if (path[k][ii]>=k) {
	global[k][ii]=true;
	globalcount[ii]+=1;
      } else global[k][ii]=false;
    }
    for (int SNP=0;SNP<nSNP;SNP++) {
      if (samplepath[SNP][k]>=k) {
	sampleglobal[SNP][k]=true;
	sampleglobalcount[SNP]+=1;
	samplepath[SNP][k] -=k;  // I think this is correct!?! 
      } else  sampleglobal[SNP][k]=true;
    }
  } 	
  
  TNT::Array2D<int> sampSNP(nSNP,n);
  TNT::Array1D<double> gammamatchlocal(n);
  for (int k=1;k<n;k++) {
  if (perfectlocalmatch)
    gammamatchlocal[k]=1.0;       // constant values througout
  else
    gammamatchlocal[k]=k/(k+thetahat) + 0.5*thetahat/(k+thetahat);;
  }

  for (int SNP=0;SNP<nSNP;SNP++) {
    for (;;) {
      int sampglobalcount=1;
      bool keep=false;
      if (r()<0.5) sampSNP[SNP][0]=1;
      else sampSNP[SNP][0]=0;
      for (int k=1;k<n;k++) {
	if (sampleglobal[SNP][k]) {
	  double gammamatch=double(sampglobalcount)/(double(sampglobalcount)+thetahat)
	    +0.5*thetahat/(double(sampglobalcount)+thetahat);
	  if (r()<gammamatch) sampSNP[SNP][k]=sampSNP[SNP][samplepath[SNP][k]];
	  else {
	    sampSNP[SNP][k]=1-sampSNP[SNP][samplepath[SNP][k]]; // a mutation 
	    keep=true;  // we can try to keep this set!
	  }
	  sampglobalcount++;
	}
	else {
	  if (r()<gammamatchlocal[k]) sampSNP[SNP][k]=sampSNP[SNP][samplepath[SNP][k]];
	  else {
	    sampSNP[SNP][k]=1-sampSNP[SNP][samplepath[SNP][k]]; // a mutation 
	    keep=true;  // we can try to keep this set!
	  }
	}
      }
      if (keep) {  // at least one mutation
	int freq=0;
	for (int i=0;i<n;i++) freq+=sampSNP[SNP][i]==1;
	if (freq>=minfreq&&n-freq> minfreq) break;  // keep this and get the next
      }
    }
  }

  return sampSNP;
}
/******************************************************************************/
/** sampleSNP for case control data (in genotypes)                            */
/* The ordering of genotypes is (if controls are labelled 0 and cases 1)
 *   
 *      Control_00 Control_01 Control_11 Case_00 Case_01 Case_11
 *  
******************************************************************************/
TNT::Array3D<int> subdivpac::SampleCaseControlSNP(TNT::Array1D<double> &weightings
						  , TNT::Array1D<int> &pos
						  ,TNT::Array1D<double> &d1
						  , double rho
						  ,  const TNT::Array1D<double> &F,
						  TNT::Array1D<int> &CaseControl
						  ,int seed)
  {
    rng r(seed);
    int nSNP=pos.dim();
    TNT::Array3D<int> tables(pos.dim(),nperm,6,0);
    current_pop=pop.copy();  
    for (int i=0;i<nperm;i++) {
      std::vector<int> rev(n);
      std::vector<int> prms=r.integer_permutations(n,n);
      for (int j=0;j<n;j++) {
	h[j]=st[prms[j]];
	current_pop[j]= pop[prms[j]];
	rev[prms[j]]=j;
      }
     
      TNT::Array2D<int> res =  sampleSNP(weightings[i],rho,F
					 ,pos,d1,path_,global_,n/5);
      for (int SNP=0;SNP<nSNP;SNP++) {
	for (int k=0;k<n;k+=2) {
	  if (res[SNP][rev[k]]==0&&res[SNP][rev[k+1]]==0) 
	    tables[SNP][i][CaseControl[k]*3]+=1;
	  
	  else if (res[SNP][rev[k]]==1&&res[SNP][rev[k+1]]==1) 
	    tables[SNP][i][2+CaseControl[k]*3]+=1;
	  
	  else  tables[SNP][i][1+CaseControl[k]*3]+=1;
	}
      }
    }    
    double scale=*(std::max_element(&weightings[0],&weightings[0]+nperm));
    for (int i=0;i<nperm;i++) 
      weightings[i] = exp(weightings[i]-scale);
    
    return tables;
  }
