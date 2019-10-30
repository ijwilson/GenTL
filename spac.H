// @file
// Time-stamp: <2012-01-27 11:34:53 ijw>
#ifndef SUBDIVPACT_H__
#define SUBDIVPACT_H__
#include "GenTL/pact.H"

/** templatised subdivided pac class.
 */
template <typename T>
class spact  : public pact<T> {
public:
  /** the constructor                                                           */
  spact(T &hps, const TNT::Array1D<double> &dis 	    
	    ,const TNT::Array1D<int> &loc	 
	,rng &r,double theta=-1.0,bool simpp=true):    
    pact<T>(hps,dis,1,1,theta)
    ,pop(loc),path_(this->n,this->sites)
    ,global_(this->n,this->sites),simpleprior(simpp) {   
    npops_=0;    
    for (int i=0;i<loc.dim();i++) if (loc[i]>npops_) npops_=loc[i];    
    permute(&this->order[0],this->order.dim(),r);
    npops_++;
  }  
  /** produce a random permutation of the ordering                              */
  void perm(rng &r) {
    permute(&this->order[0],this->n,r);
  }
  /** calculate the log of the pre-data probability of getting the 
   * observed path and globalness                                               */  
  double logpriorpath(double rho, const double *F) {
  //  double lpsimple= logpriorpathsimple(rho,F,path_,global_);
  //  double lp= logpriorpath(rho,F,path_,global_);
  //  std::cout << lpsimple << " " << lp << " " << lpsimple-lp << std::endl;
    if (simpleprior) return logpriorpathsimple(rho,F,path_,global_);
    else return logpriorpath(rho,F,path_,global_);
  }
  /** what is the prior probability of this configuration of globals and paths  */
  double logpriorpath(double rho, const double *F
		  ,const  TNT::Array2D<int> &path
		  , const TNT::Array2D<bool> &global) const;
  double logpriorpathsimple(double rho, const double *F
		  ,const  TNT::Array2D<int> &path
		  , const TNT::Array2D<bool> &global) const;
  /** Sample a new path
   * returns the probability of this particular sample                          */
  double samplepath(double rho, const double *F,rng &r) {
    return samplepath(this->thetahat,rho,F,path_,global_,r);
  }
  /**  How many populations ?                                                  */
  int npops() const {
    return npops_;
  }
  /** The likelihood of the path P(D|path,theta) - independend of F and rho   */
  double likelihoodpath(double theta=-1)  {
    if (theta<0.0) theta=this->thetahat;
    return likelihoodpath(theta,path_,global_);
  }

  double likelihoodpath(double theta,const  TNT::Array2D<int> &path
			, const TNT::Array2D<bool> &global) const ;
  
  double samplepath(double theta, double rho,const double *F
		    ,TNT::Array2D<int> &path,TNT::Array2D<bool> &global,rng &r);

  /** The probability of sampling this path conditional on F and rho        */
  double lpsamplepath(double rho, const double *F) const {
    return lpsamplepath(this->thetahat,rho,F,path_,global_);
  }
  /** The probability of sampling this path conditional on F and rho        */
  double lpsamplepath(double theta, double rho, const double *F
		       ,const TNT::Array2D<int> &path
		      ,const TNT::Array2D<bool> &global) const ;



  std::pair<double,double> sample(double theta,double rho, const double *F
		,TNT::Array2D<bool> &global,rng &r);

  std::pair<double,double> sample(double rho, const double *F,rng &r) {
    return sample(this->thetahat,rho,F,global_,r);
  }

  std::pair<double,double> lpsample(double theta, double rho, const double *F
		      ,const TNT::Array2D<bool> &global) const ;
  std::pair<double,double> lpsample(double rho, const double *F) const {
    return lpsample(this->thetahat,rho,F,global_);
  }

  double logprior(double rho, const double *F
		  , const TNT::Array2D<bool> &global) const;

  double logprior(double rho, const double *F) const {
    return logprior(rho,F,global_);
  }

  /** Accessor functions to return t the current path                   */
  TNT::Array2D<int> path() const {
    return path_;
  }
  /** Accessor function to return  the current "globalness"             */
  TNT::Array2D<bool> global() const {
    return global_;
  }
  /** change the path to these new values                               */
  void changepath(const TNT::Array2D<int> &path,const TNT::Array2D<bool> &global) {
    global_.inject(global);
    path_.inject(path);
  }
  /** change the path to these new values                               */
  void changeglobal(const TNT::Array2D<bool> &global) {
    global_.inject(global);
  }
protected:  
  const TNT::Array1D<int> &pop;  
  int npops_;  
  TNT::Array2D<int> path_;										  
  TNT::Array2D<bool> global_;  
  bool simpleprior;	
  // rng r;  
  spact(); // default constructor -- don't want it used  
  /** calculate the likelihood for a single rho for one   * permutation of the sample  */  
//   double calc(double rho, const TNT::Array1D<double> &F, TNT::Array2D<int> &path	      
// 	    , TNT::Array2D<bool> &global,bool resamplepath=true); 

  double getalpha(double theta,double rho, const double *F,
		  TNT::Array2D<double> &alpha,
		  const TNT::Array2D<bool> &global, 
		  const TNT::Array1D<int> &popcount,
		  const TNT::Array1D<int> &globalcount, int k) const ;

  double getalpha_prior(double rho, const double *F,
			  TNT::Array2D<double> &alpha,
			  const TNT::Array2D<bool> & global, 
			  const TNT::Array1D<int> &popcount,
			const TNT::Array1D<int> &globalcount, int k) const;
};

/** calculate the likelihood for one rearrangement of the haplotypes 
 * which is assumed to have already happened.    We also have to sample 
 * (as we go) whether or not each position is in the "local" or global population */
template <typename T>
std::pair<double,double> spact<T>::sample(double theta,double rho, const double *F
			    ,TNT::Array2D<bool> &global,rng &r)
{
  // path holds the path through the haplotypes, we might as well keep it 
  // as it costs little.  Global holds the copy type-local or from the global pool 
  // alpha[x,s] defined as the probability of sampling the first s sites of
  // your haplotype, and of sampling haplotype x. The first k hold the 
  // probabilites of a local copy of the site, the last k copying from the global pool.
  // Some (or most) of these will be zero (zero locally and globally if they are 
  // local to  another population, zero globally if they are local to the current
  // population, only non zero for both if they are both local and in the global pool)

  TNT::Array2D<double> alpha(pact<T>::sites,2*pact<T>::n-2);
  TNT::Array1D<int> popcount(npops_,0);     // the count in each of the populations
  TNT::Array1D<int> globalcount(this->sites,1); // the count in the "global" pool
 
  double lpsample=0.0,lprob=0.0;
  popcount[pop[this->order[0]]] +=1;           // add one to the local count
  for (int s=0;s<pact<T>::sites;s++) global[0][s]=true;

  // now consider each of the subsequent samples in turn.
  // Because my arrays are 0 offset haplotype k is the k+1 st sample
  // so use the same k as in Li and Stephens' appendix
  for (int k=1;k<this->n;k++) {
    int location=pop[this->order[k]];
    double plocalcopy= 
      (popcount[location]>0)?F[location]*double(popcount[location])
    /(1.+(double)(popcount[location]-1.)*F[location]):0.0;
    // get the alpha's for the forward part of the algorithm
    double logscalefactor=getalpha(theta,rho,F,alpha,global,popcount,globalcount,k);  
    double prob=std::accumulate(alpha[this->sites-1],alpha[this->sites-1]+2*k,0.0);
    lprob+=log(prob)+logscalefactor;
    double lp[2];
    lp[0]=std::accumulate(&alpha[this->sites-1][0],&alpha[this->sites-1][k],0.0);
    lp[1]=std::accumulate(&alpha[this->sites-1][0]+k,&alpha[this->sites-1][0]+2*k,0.0);

    //std::cout << "lp " << lp[0] << " " << lp[1] ;

    TNT::Array1D<double> lpabove(k,0.0);
    if (r()<lp[0]/(lp[0]+lp[1])) {
      lpsample += log(lp[0]) - log(lp[0]+lp[1]);
      global[k][this->sites-1]=false;
      for (int jj=0;jj<k;jj++) lpabove[jj]=alpha[this->sites-1][jj];      
    } else {
      lpsample += log(lp[1]) - log(lp[0]+lp[1]);
      global[k][this->sites-1]=true;
      for (int jj=0;jj<k;jj++) lpabove[jj]=alpha[this->sites-1][jj+k];      
    }
    //std::cout << " lpsample " << lpsample << std::endl;
    //   std::cout << "##################lpabove###############################\n";
    //   std::cout << lpabove;
    
  
    TNT::Array1D<double> currlp(2*k);
    for (int s=this->sites-2;s>=0;s--) {
      for (int jj=0;jj<2*k;jj++) currlp[jj]=0.0;
      double p_local_norec=exp(-rho*this->d[s]/(double(popcount[location])));
      double p_global_norec=exp(-rho*this->d[s]/(double(globalcount[s])));
      for (int x=0;x<k;x++) { // loop throught all the possible haplotypes
	                      // this haplotype could be copied from
	if (global[x][s]) {                          // copying from a global x
	  //	  std::cout << "copying from global\n";
	  for (int xprime=0;xprime<k;xprime++) {     
	    // loop through the paths above that it could be copied from
	    if (global[k][s+1]) {                    //  global x --> xprime global 
	      if (x==xprime) { 
		//		std::cout << "here\n";    
		currlp[x+k]+=lpabove[xprime]*alpha[s][x+k]
		  *(p_global_norec+(1.0-p_global_norec)*(1.-plocalcopy)/double(globalcount[s+1]));
		//		std::cout << "newcurrlp---------------------------------" 
		//			  << std::endl << currlp << std::cout;
	      } else { 
		currlp[x+k]+=lpabove[xprime]*alpha[s][x+k]
		  *(1.0-p_global_norec)*(1.-plocalcopy)/double(globalcount[s+1]);
	      }
	      if (pop[this->order[x]]==location) {  // only local for same population
		currlp[x]+=lpabove[xprime]*alpha[s][x]
		  *(1.0-p_local_norec)*(1.-plocalcopy)/double(globalcount[s+1]);
	      } 
	    } else {                                // global x --> xprime local
	      currlp[x+k]+=lpabove[xprime]*alpha[s][x+k]
		*(1.0-p_global_norec)*plocalcopy/double(popcount[location]);
	      if (pop[this->order[x]]==location) {
		currlp[x]+=lpabove[xprime]*alpha[s][x]
		  *(1.0-p_global_norec)*plocalcopy/double(popcount[location]);
	      } 
	    }
	  }
	} else { // copying from a local chromosome
	  for (int xprime=0;xprime<k;xprime++) {
	    // loop through the paths above that it could be copied from
	    if (global[k][s+1]) {                  // local x --> xprime global 
	      currlp[x+k]+=0.0; // must be copied globally
	      currlp[x]+=lpabove[xprime]*alpha[s][x]
		*(1.-p_local_norec)*(1.-plocalcopy)/double(globalcount[s+1]);
	    } else {                               // local x --> local xprime 
	      if (x==xprime) {
		currlp[x]+=lpabove[xprime]*alpha[s][x]
		  *(p_local_norec+(1.0-p_local_norec)*(1.-plocalcopy)/double(popcount[location]));
	      } else {
		if (pop[this->order[x]]==location)
		  currlp[x]+=lpabove[xprime]*alpha[s][x]
		    *(1.0-p_local_norec)*plocalcopy/double(popcount[location]);
	      } 
	    }
	  }
	}
      }
      //std::cout << lpabove;
      
      lp[0]=std::accumulate(&currlp[0],&currlp[0]+k,0.0);
      lp[1]=std::accumulate(&currlp[0]+k,&currlp[0]+2*k,0.0);
      // std::cout << k << " " <<  lp[0] << " " << lp[1] << std::endl;
      if (r()<lp[0]/(lp[0]+lp[1])) {
	lpsample += log(lp[0]) - log(lp[0]+lp[1]);
	global[k][s]=false;
	for (int jj=0;jj<k;jj++) lpabove[jj]=currlp[jj];
      } else {
	lpsample += log(lp[1]) - log(lp[0]+lp[1]);
	global[k][s]=true;
	for (int jj=0;jj<k;jj++) lpabove[jj]=currlp[jj+k];
      }
      //std::cout << k << " " << s << " "  << lpsample << std::endl ;
 
    }
        //  lpabove.inject(currlp);
    for (int jj=0;jj<this->sites;jj++) if (global[k][jj]) globalcount[jj]+=1;

    popcount[location]+=1;
  }
  return std::pair<double,double>(lprob,lpsample);
}
/** calculate the likelihood for one rearrangement of the haplotypes 
 * which is assumed to have already happened.    We also have to sample 
 * (as we go) whether or not each position is in the "local" or global population 
 * returns the likelihood for the data and the probability of sampling the globals
 * under this model 
*/
template <typename T>
std::pair<double,double> spact<T>::lpsample(double theta,double rho, const double *F
			    ,const TNT::Array2D<bool> &global) const 
{
  TNT::Array2D<double> alpha(this->sites,2*this->n-2);
  TNT::Array1D<int> popcount(npops_,0);     // the count in each of the populations
  TNT::Array1D<int> globalcount(this->sites,1); // the count in the "global" pool
 
  double lpsample=0.0,lprob=0.0;
  popcount[pop[this->order[0]]] +=1;           // add one to the local count
  for (int k=1;k<this->n;k++) {
    int location=pop[this->order[k]];
    double plocalcopy= 
      (popcount[location]>0)?F[location]*double(popcount[location])
      /(1.+(double)(popcount[location]-1.)*F[location]):0.0;
    // get the alpha's for the forward part of the algorithm
    double logscalefactor= getalpha(theta,rho,F,alpha,global,popcount,globalcount,k);  // don't need the probability
    double prob=std::accumulate(alpha[this->sites-1],alpha[this->sites-1]+2*k,0.0);
    lprob+=log(prob)+logscalefactor;
 
    double lp[2];
    lp[0]=std::accumulate(&alpha[this->sites-1][0],&alpha[this->sites-1][0]+k,0.0);
    lp[1]=std::accumulate(&alpha[this->sites-1][0]+k,&alpha[this->sites-1][0]+2*k,0.0);

 
    TNT::Array1D<double> lpabove(k,0.0);
    if (global[k][this->sites-1]) {
      lpsample += log(lp[1]) - log(lp[0]+lp[1]);
      for (int jj=0;jj<k;jj++) lpabove[jj]=alpha[this->sites-1][jj+k];
    } else {
      lpsample += log(lp[0]) - log(lp[0]+lp[1]);
      for (int jj=0;jj<k;jj++) lpabove[jj]=alpha[this->sites-1][jj];
    }

 
    TNT::Array1D<double> currlp(2*k);
    for (int s=this->sites-2;s>=0;s--) {
      for (int jj=0;jj<2*k;jj++) currlp[jj]=0.0;
      double p_local_norec=exp(-rho*this->d[s]/(double(popcount[location])));
      double p_global_norec=exp(-rho*this->d[s]/(double(globalcount[s])));
      for (int x=0;x<k;x++) { // loop throught all the possible haplotypes
	                      // this haplotype could be copied from
	if (global[x][s]) {                          // copying from a global x
	  for (int xprime=0;xprime<k;xprime++) {     
	    // loop through the paths above that it could be copied from
	    if (global[k][s+1]) {                    //  global x --> xprime global 
	      if (x==xprime) {     
		currlp[x+k]+=lpabove[xprime]*alpha[s][x+k]
		  *(p_global_norec+(1.0-p_global_norec)*(1.-plocalcopy)/double(globalcount[s+1]));
	      } else { 
		currlp[x+k]+=lpabove[xprime]*alpha[s][x+k]
		  *(1.0-p_global_norec)*(1.-plocalcopy)/double(globalcount[s+1]);
	      }
	      if (pop[this->order[x]]==location) {  // only local for same population
		currlp[x]+=lpabove[xprime]*alpha[s][x]
		  *(1.0-p_local_norec)*(1.-plocalcopy)/double(globalcount[s+1]);
	      } 
	    } else {                                // global x --> xprime local
	      currlp[x+k]+=lpabove[xprime]*alpha[s][x+k]
		*(1.0-p_global_norec)*plocalcopy/double(popcount[location]);
	      if (pop[this->order[x]]==location) {
		currlp[x]+=lpabove[xprime]*alpha[s][x]
		  *(1.0-p_global_norec)*plocalcopy/double(popcount[location]);
	      } 
	    }
	  }
	} else { // copying from a local chromosome
	  for (int xprime=0;xprime<k;xprime++) {
	    // loop through the paths above that it could be copied from
	    if (global[k][s+1]) {                  // local x --> xprime global 
	      currlp[x+k]+=0.0; // must be copied globally
	      currlp[x]+=lpabove[xprime]*alpha[s][x]
		*(1.-p_local_norec)*(1.-plocalcopy)/double(globalcount[s+1]);
	    } else {                               // local x --> local xprime 
	      if (x==xprime) {
		currlp[x]+=lpabove[xprime]*alpha[s][x]
		  *(p_local_norec+(1.0-p_local_norec)*(1.-plocalcopy)/double(popcount[location]));
	      } else {
		if (pop[this->order[x]]==location)
		  currlp[x]+=lpabove[xprime]*alpha[s][x]
		    *(1.0-p_local_norec)*plocalcopy/double(popcount[location]);
	      } 
	    }
	  }
	}
      }
      lp[0]=std::accumulate(&currlp[0],&currlp[k],0.0);
      lp[1]=std::accumulate(&currlp[k],&currlp[0]+2*k,0.0);
      // std::cout << k << " " <<  lp[0] << " " << lp[1] << std::endl;
      if (global[k][s]) {
	lpsample += log(lp[1]) - log(lp[0]+lp[1]);
	for (int jj=0;jj<k;jj++) lpabove[jj]=currlp[jj+k];
      } else {
	lpsample += log(lp[0]) - log(lp[0]+lp[1]);
	for (int jj=0;jj<k;jj++) lpabove[jj]=currlp[jj];
      }
     }

    for (int ii=0;ii<this->sites;ii++) if (global[k][ii]) globalcount[ii]+=1; 
    popcount[location]+=1;

  }
  return std::pair<double,double> (lprob,lpsample);
}
/** return the log of probability of sampling global from the  prior 
    -- without a path                                                 */
template <typename T>
double spact<T>::logprior(double rho, const double *F
			  , const TNT::Array2D<bool> &global) const
{
  TNT::Array2D<double> alpha(this->sites,2*this->n-2);
  TNT::Array1D<int> popcount(npops_,0);     // the count in each of the populations
  TNT::Array1D<int> globalcount(this->sites,1); // the count in the "global" pool
 
  double lpsample=0.0;
  popcount[pop[this->order[0]]] +=1;           // add one to the local count
  for (int k=1;k<this->n;k++) {
    int location=pop[this->order[k]];
    double plocalcopy= 
      (popcount[location]>0)?F[location]*double(popcount[location])
      /(1.+(double)(popcount[location]-1.)*F[location]):0.0;   
    // get the alpha's for the forward part of the algorithm
    getalpha_prior(rho,F,alpha,global,popcount,globalcount,k);  // don't need the probability
    
    double lp[2];
    lp[0]=std::accumulate(&alpha[this->sites-1][0],&alpha[this->sites-1][0]+k,0.0);
    lp[1]=std::accumulate(&alpha[this->sites-1][k],&alpha[this->sites-1][0]+2*k,0.0);
    
    TNT::Array1D<double> lpabove(k,0.0);
    if (global[k][this->sites-1]) {
      lpsample += log(lp[1]) - log(lp[0]+lp[1]);
      for (int jj=0;jj<k;jj++) lpabove[jj]=alpha[this->sites-1][jj+k];
    } else {
      lpsample += log(lp[0]) - log(lp[0]+lp[1]);
      for (int jj=0;jj<k;jj++) lpabove[jj]=alpha[this->sites-1][jj];
    }
    TNT::Array1D<double> currlp(2*k);
    for (int s=this->sites-2;s>=0;s--) {
      //  if (global[k][s+1]) std::cout << "global\n";
      // else std::cout << "local\n";  
      // std::cout << lpabove;
      for (int jj=0;jj<2*k;jj++) currlp[jj]=0.0;
      double p_local_norec=exp(-rho*this->d[s]/(double(popcount[location])));
      double p_global_norec=exp(-rho*this->d[s]/(double(globalcount[s])));
      for (int x=0;x<k;x++) { // loop throught all the possible haplotypes
	                      // this haplotype could be copied from
	if (global[x][s]) {                          // copying from a global x
	  for (int xprime=0;xprime<k;xprime++) {     
	    // loop through the paths above that it could be copied from
	    if (global[k][s+1]) {                    //  global x --> xprime global 
	      if (x==xprime) {     
		currlp[x+k]+=lpabove[xprime]*alpha[s][x+k]
		  *(p_global_norec+(1.0-p_global_norec)*(1.-plocalcopy)/double(globalcount[s+1]));
	      } else { 
		currlp[x+k]+=lpabove[xprime]*alpha[s][x+k]
		  *(1.0-p_global_norec)*(1.-plocalcopy)/double(globalcount[s+1]);
	      }
	      if (pop[this->order[x]]==location) {  // only local for same population
		currlp[x]+=lpabove[xprime]*alpha[s][x]
		  *(1.0-p_local_norec)*(1.-plocalcopy)/double(globalcount[s+1]);
	      } 
	    } else {                                // global x --> xprime local
	      currlp[x+k]+=lpabove[xprime]*alpha[s][x+k]
		*(1.0-p_global_norec)*plocalcopy/double(popcount[location]);
	      if (pop[this->order[x]]==location) {
		currlp[x]+=lpabove[xprime]*alpha[s][x]
		  *(1.0-p_global_norec)*plocalcopy/double(popcount[location]);
	      } 
	    }
	  }
	} else { // copying from a local chromosome
	  for (int xprime=0;xprime<k;xprime++) {
	    // loop through the paths above that it could be copied from
	    if (global[k][s+1]) {                  // local x --> xprime global 
	      currlp[x+k]+=0.0; // must be copied globally
	      currlp[x]+=lpabove[xprime]*alpha[s][x]
		*(1.-p_local_norec)*(1.-plocalcopy)/double(globalcount[s+1]);
	    } else {                               // local x --> local xprime 
	      if (x==xprime) {
		currlp[x]+=lpabove[xprime]*alpha[s][x]
		  *(p_local_norec+(1.0-p_local_norec)*(1.-plocalcopy)/double(popcount[location]));
	      } else {
		if (pop[this->order[x]]==location)
		  currlp[x]+=lpabove[xprime]*alpha[s][x]
		    *(1.0-p_local_norec)*plocalcopy/double(popcount[location]);
	      } 
	    }
	  }
	}
      }
      lp[0]=std::accumulate(&currlp[0],&currlp[k],0.0);
      lp[1]=std::accumulate(&currlp[k],&currlp[2*k],0.0);
      if (global[k][s]) {
	lpsample += log(lp[1]) - log(lp[0]+lp[1]);
	globalcount[s]+=1;
	for (int jj=0;jj<k;jj++) lpabove[jj]=currlp[jj+k];
      } else {
	lpsample += log(lp[0]) - log(lp[0]+lp[1]);
	for (int jj=0;jj<k;jj++) lpabove[jj]=currlp[jj];
      }
    }

    for (int ii=0;ii<this->sites;ii++) if (global[k][ii]) globalcount[ii]+=1; 
    popcount[location]+=1;
  }
  return lpsample;
}
/** calculate the likelihood for one rearrangement of the haplotypes 
 * which is assumed to have already happened.    We also have to sample 
 * (as we go) whether or not each position is in the "local" or global population */
template <typename T>
double spact<T>::samplepath(double theta,double rho, const double *F
			    ,TNT::Array2D<int> &path,TNT::Array2D<bool> &global,rng &r)
{
  // path holds the path through the haplotypes, we might as well keep it 
  // as it costs little.  Global holds the copy type-local or from the global pool 
  // alpha[x,s] defined as the probability of sampling the first s sites of
  // your haplotype, and of sampling haplotype x. The first k hold the 
  // probabilites of a local copy of the site, the last k copying from the global pool.
  // Some (or most) of these will be zero (zero locally and globally if they are 
  // local to  another population, zero globally if they are local to the current
  // population, only non zero for both if they are both local and in the global pool)

  TNT::Array2D<double> alpha(this->sites,2*this->n-2);
  TNT::Array1D<int> popcount(npops_,0);     // the count in each of the populations
  TNT::Array1D<int> globalcount(this->sites,1); // the count in the "global" pool
 
  double lpsample=0.0;
  popcount[pop[this->order[0]]] +=1;           // add one to the local count
  for (int s=0;s<pact<T>::sites;s++)  {
    global[0][s]=true;
    path[0][s]=0;
  }

   // now consider each of the subsequent samples in turn.
  // Because my arrays are 0 offset haplotype k is the k+1 st sample
  // so use the same k as in Li and Stephens' appendix
  for (int k=1;k<this->n;k++) {
    int location=pop[this->order[k]];
    // get the alpha's for the forward part of the algorithm
    getalpha(theta,rho,F,alpha,global,popcount,globalcount,k);  
    double plocalcopy= 
      (popcount[location]>0)?F[location]*double(popcount[location])
      /(1.+(double)(popcount[location]-1.)*F[location]):0.0;
    // Forwards backwards - first sample the last site 
    double psample;
    path[k][this->sites-1] 
      = gen_from_p(alpha[this->sites-1],alpha[this->sites-1]+2*k,r,psample);
    assert(psample>0.0);
    lpsample+=log(psample);
    
    //   std::cout << "k = " << k << " " << psample << " " << lpsample << std::endl;
    //   for (int ii=0;ii<2*k;ii++) std::cout << alpha[this->sites-1][ii]<< " ";
    //   std::cout << std::endl;
    // need to calculate the transition probabilities here
    std::vector<double> lp(2*k);  
    for (int s=this->sites-2;s>=0;s--) {
      double p_local_norec=exp(-rho*this->d[s]/(double(popcount[location])));
      double p_global_norec=exp(-rho*this->d[s]/(double(globalcount[s])));
      int xprime=path[k][s+1];
      for (int x=0;x<k;x++) { // loop throught all the possible haplotypes
	if (xprime<k) { // current path (ahead) local
	  assert(popcount[location]>0);
	  if (xprime==x) {       // staying on the same path
	    lp[x] = alpha[s][x]*
	      (p_local_norec+(1.-p_local_norec)*plocalcopy/double(popcount[location]));
	  } else  {              // move
	    lp[x] = alpha[s][x]*
	      (1.-p_local_norec)*plocalcopy/double(popcount[location]);
	  } 
	  lp[k+x] = alpha[s][k+x]*
	    (1.-p_global_norec)*plocalcopy/double(popcount[location]);
	  
	} else {  //  path ahead global
	  assert(global[xprime-k][s+1]);  // must be a copy of a global
	  if (xprime==k+x) { // copy from x as global
	    lp[k+x] = alpha[s][k+x]*
	      (p_global_norec+(1.0-p_global_norec)*(1.-plocalcopy)/double(globalcount[s+1]));
	  } else {
	    lp[k+x] = alpha[s][k+x]*
	      (1.-p_global_norec)*(1.-plocalcopy)/double(globalcount[s+1]);
	  }
	  lp[x] = alpha[s][x]*
	    (1.-p_local_norec)*(1.-plocalcopy)/double(globalcount[s+1]);
	}
      }
      path[k][s]=gen_from_pr(lp,2*k,r,psample);
      assert(psample>0.0);
      lpsample += log(psample);
      //  std::cout << psample << " " << lpsample << std::endl;
    }
    for (int ii=0;ii<this->sites;ii++) {
      if (path[k][ii]>=k) {
	global[k][ii]=true;
	globalcount[ii]+=1;
      } else global[k][ii]=false;
    }
    popcount[location]+=1;
  }
  return lpsample;
}
/** calculate the likelihood for one rearrangement of the haplotypes 
 * which is assumed to have already happened.    We also have to sample 
 * (as we go) whether or not each position is in the "local" or global population */
template <typename T>
double spact<T>::lpsamplepath(double theta,double rho, const double *F
			      ,const TNT::Array2D<int> &path,
			      const TNT::Array2D<bool> &global) const 
{
  TNT::Array2D<double> alpha(this->sites,2*this->n-2);
  TNT::Array1D<int> popcount(npops_,0);     // the count in each of the populations
  TNT::Array1D<int> globalcount(this->sites,1); // the count in the "global" pool
 
  double lpsample=0.0,totp;
  popcount[pop[this->order[0]]] +=1;           // add one to the local count
  for (int k=1;k<this->n;k++) {
    int location=pop[this->order[k]];
    // get the alpha's for the forward part of the algorithm
    getalpha(theta,rho,F,alpha,global,popcount,globalcount,k);  // don't need the probability
   
    double plocalcopy= 
      (popcount[location]>0)?F[location]*double(popcount[location])
      /(1.+(double)(popcount[location]-1.)*F[location]):0.0;
    // Method as in Boys et al 2000 - the last site is "site-1".  sample from the "equilibrium"

    totp = std::accumulate(alpha[this->sites-1],alpha[this->sites-1]+2*k,0.0);
    lpsample += log(alpha[this->sites-1][path[k][this->sites-1]])-log(totp);
 
    // need to calculate the transition probabilities here
    std::vector<double> lp(2*k);  
    for (int s=this->sites-2;s>=0;s--) {
      double p_local_norec=exp(-rho*this->d[s]/(double(popcount[location])));
      double p_global_norec=exp(-rho*this->d[s]/(double(globalcount[s])));
      int xprime=path[k][s+1];   // xprime is the path (ahead) - where it has come from
      for (int x=0;x<k;x++) {    // loop throught all the possible haplotypes that is could be copied from
	if (xprime<k) {          // current path (ahead) local
	  assert(popcount[location]>0);
	  if (xprime==x) {       // staying on the same path
	    lp[x] = alpha[s][x]*
	      (p_local_norec+(1.-p_local_norec)*plocalcopy/double(popcount[location]));
	  } else  {              // move
	    lp[x] = alpha[s][x]*
	      (1.-p_local_norec)*plocalcopy/double(popcount[location]);
	  } 
	  lp[k+x] = alpha[s][k+x]*
	    (1.-p_global_norec)*plocalcopy/double(popcount[location]); 
	} else {  // current path global
	  assert(global[xprime-k][s+1]);
	  if (xprime==k+x) { // copy from x as global
	    lp[k+x] = alpha[s][k+x]*
	      (p_global_norec+(1.0-p_global_norec)*(1.-plocalcopy)/double(globalcount[s+1]));
	  } else {
	    lp[k+x] = alpha[s][k+x]*
	      (1.-p_global_norec)*(1.-plocalcopy)/double(globalcount[s+1]);
	  }
	  //	  if (popcount[location]>0)
	    lp[x] = alpha[s][x]*
	      (1.-p_local_norec)*(1.-plocalcopy)/double(globalcount[s+1]);
	    //	  else lp[x]=0.0;
	}
      }
      
      totp = std::accumulate(&lp[0],&lp[0]+2*k,0.0);
      lpsample += log(lp[path[k][s]])-log(totp);
    }
    for (int ii=0;ii<this->sites;ii++) 
      if (global[k][ii]) globalcount[ii]+=1;
 
    popcount[location]+=1;
  }
  return lpsample;
}
/** calculate the likelihood for the haplotypes conditional
 * on the path - no dependency on F 
 */
template <typename T>
double spact<T>::likelihoodpath(double theta,const  TNT::Array2D<int> &path
				,const TNT::Array2D<bool> &global) const 
{
  double lp=0.0;
  TNT::Array1D<int> globalcount(this->sites,1);

  for (int k=1;k<this->n;k++) {
    for (int s=0;s<this->sites;s++) {
      int wh=path[k][s];
      if (global[k][s]) {
	assert(path[k][s]>=k);
	if (this->haps[this->order[k]][s]==this->haps[this->order[wh-k]][s]) { //matches
	  lp += log(double(globalcount[s])/(double(globalcount[s])+theta)
		    +0.5*theta/(double(globalcount[s])+theta)); 
	} else {
	  lp += log(0.5*theta/(double(globalcount[s])+theta));
	}
	globalcount[s]+=1;
      } else {
	assert(wh<k);
	assert(this->haps[this->order[k]][s]==this->haps[this->order[wh]][s]);
      }
    }
  }
  return lp;
}

/** A simpler, and much quicker, calculation of the prior for a path */
template <typename T>
double spact<T>::logpriorpathsimple(double rho, const double *F
		  ,const  TNT::Array2D<int> &path
		  , const TNT::Array2D<bool> &global) const
{
  TNT::Array1D<int> popcount(npops_,0);         // count in each population
  TNT::Array1D<int> globalcount(this->sites,1); // count in  "global" pool
 
  double lpsample=0.0;                         // first path has probability 1
  popcount[pop[this->order[0]]] =1;            // one in local count
  // std::cerr << "Problem -- this is incorrect......" << std::endl;
  // But why is it incorrect ?  Because sometimes a global is forced?
  // but I can't remember why - if the number of globals changes 
  // I think that this needs some more work.
 
  for (int k=1;k<this->n;k++) {
    int location=pop[this->order[k]];
    double plocalcopy= 
      (popcount[location]>0)?F[location]*double(popcount[location])
      /(1.+(double)(popcount[location]-1.)*F[location]):0.0;
   
    if (global[k][0]) 
      lpsample += log((1.-plocalcopy)/(double(globalcount[0])));
    else lpsample += log(plocalcopy/(double(popcount[location])));

    for (int s=1;s<this->sites;s++) {
      int xprime=path[k][s-1]; // current path (behind)
      if (global[k][s-1]) {    // path (behind) global
	double p_global_norec=exp(-rho*this->d[s-1]/(double(globalcount[s-1])));
	if (xprime==path[k][s]) // no change in a global path
	  lpsample += log(p_global_norec 
			  +(1.-p_global_norec)*(1.-plocalcopy)/double(globalcount[s])
			  );
	else {
	  if (global[k][s]) lpsample +=  log(1.-p_global_norec) 
	    + log(1.-plocalcopy) -log(double(globalcount[s]));
	  else lpsample += log(1.-p_global_norec) 
	    +log(plocalcopy) - log(double(popcount[location]));
	}
      } else {                 // path behind local
	double p_local_norec=exp(-rho*this->d[s-1]/(double(popcount[location])));
	if (xprime==path[k][s])
	  lpsample += log(p_local_norec 
			  + (1.-p_local_norec)*plocalcopy/double(popcount[location])
			  );
	else {
	   lpsample += log(1.-p_local_norec);
	   if (global[k][s]) lpsample += log(1.-plocalcopy) - log(double(globalcount[s]));
	   else lpsample += log(plocalcopy) - log(double(popcount[location]));
	}
      }
    }
    for (int ii=0;ii<this->sites;ii++) 
      if (global[k][ii]) globalcount[ii]+=1;
    popcount[location]+=1;
  }
  return lpsample;
}
/** Calculate the probability of the prior path using the same machinery as for 
    the likelihood calculations                                                      */
template <typename T>
double spact<T>::logpriorpath(double rho, const double *F
		  ,const  TNT::Array2D<int> &path
		  , const TNT::Array2D<bool> &global) const
{
  TNT::Array2D<double> alpha(this->sites,2*this->n-2);
  TNT::Array1D<int> popcount(npops_,0);     // the count in each of the populations
  TNT::Array1D<int> globalcount(this->sites,1); // the count in the "global" pool
 
  double lpsample=0.0,totp;
  popcount[pop[this->order[0]]] +=1;           // add one to the local count
  for (int k=1;k<this->n;k++) {
    int location=pop[this->order[k]];
    // get the alpha's for the forward part of the algorithm
    getalpha_prior(rho,F,alpha,global,popcount,globalcount,k);  // don't need the probability
    double plocalcopy= 
      (popcount[location]>0)?F[location]*double(popcount[location])
      /(1.+(double)(popcount[location]-1.)*F[location]):0.0;
    // Method as in Boys et al 2000 - the last site is "site-1".  sample from the "equilibrium"
    totp = std::accumulate(alpha[this->sites-1],alpha[this->sites-1]+2*k,0.0);
    lpsample += log(alpha[this->sites-1][path[k][this->sites-1]])-log(totp);
    // need to calculate the transition probabilities here
    std::vector<double> lp(2*k);  
    for (int s=this->sites-2;s>=0;s--) {
      double p_local_norec=exp(-rho*this->d[s]/(double(popcount[location])));
      double p_global_norec=exp(-rho*this->d[s]/(double(globalcount[s])));
      bool glob=global[k][s+1];
      for (int x=0;x<k;x++) { // loop throught all the possible haplotypes
	if (!glob) { // path ahead local
	  int xprime=path[k][s+1];
	  assert(popcount[location]>0);
	  if (xprime==x) {       // copying to the same path
	    assert(pop[this->order[x]]==location);
	    lp[x] = alpha[s][x]*
	      (p_local_norec
	       +(1.-p_local_norec)*plocalcopy/double(popcount[location]));
	  } else  {              // calculations for moving from path x to a local path
	    lp[x] = alpha[s][x]*
	      (1.-p_local_norec)*plocalcopy/double(popcount[location]);
	  }
	  lp[x+k]=alpha[s][x+k]*
	    (1.-p_global_norec)*plocalcopy/double(popcount[location]);
	} else {  // path ahead global
	  int xprime=path[k][s+1]-k;
	  assert(global[xprime][s+1]);
	  if (xprime==x) { // copy from x as global
	    assert(global[x][s+1]);
	    lp[k+x] = alpha[s][k+x]*(p_global_norec
		       +(1.0-p_global_norec)*(1.-plocalcopy)/double(globalcount[s+1])
				     );
	  } else {
	    lp[k+x] = alpha[s][k+x]*
	      (1.-p_global_norec)*(1.-plocalcopy)/double(globalcount[s+1]);
	  }
	  lp[x] = alpha[s][x]*
	      (1.-p_local_norec)*(1.-plocalcopy)/double(globalcount[s+1]);
	}
      }
      totp = std::accumulate(&lp[0],&lp[0]+2*k,0.0);
      lpsample += log(lp[path[k][s]])-log(totp);
    }
    for (int ii=0;ii<this->sites;ii++) 
      if (global[k][ii]) globalcount[ii]+=1;
    
    popcount[location]+=1;
  }
  return lpsample;
}

/**  A function to calculate the alpha matrix for sample k for (k>0)
 *   This assumes that we shall sample the global/local status of sample
 *   k rather than assuming it to be known            
 */
template <typename T>
double spact<T>::getalpha(double theta,double rho, const double *F,
			  TNT::Array2D<double> &alpha,
			  const TNT::Array2D<bool> &global, 
			  const TNT::Array1D<int> &popcount,
			  const TNT::Array1D<int> &globalcount, int k) const
{
  int location = pop[this->order[k]];  // the location of haplotype k
  // the probability that a recombination is copied from a "local" chromosome
  double plocalcopy= 
    (popcount[location]>0)?F[location]*double(popcount[location])
    /(1.+(double)(popcount[location]-1.)*F[location]):0.0;

  // what are the match probabilities for the global population
  double gammamatchglobal=(double(globalcount[0])+0.5*theta)/(double(globalcount[0])+theta);
  // and the local population
  double gammamatchlocal=1.0;      
  // we can calculate the probabilities for the first site directly
  // there are no transitions and we assume exchangeability of haplotypes 
  // loop over the haplotypes from which the site could be copied
  for (int x=0;x<k;x++) {
    if (pop[this->order[x]]==location) 
      alpha[0][x] = plocalcopy/double(popcount[location]);  
    else alpha[0][x]=0.0;
    
    if (global[x][0])  // can't be copied globally from a local copy
      alpha[0][k+x]= (1.-plocalcopy)/double(globalcount[0]); 
    else alpha[0][k+x]=0.0;

    if (this->haps[this->order[x]][0]==this->haps[this->order[k]][0]) {
      alpha[0][x] *= gammamatchlocal;
      alpha[0][k+x] *= gammamatchglobal;
    } else {
      alpha[0][x] *= (1.0-gammamatchlocal);
      alpha[0][k+x] *= (1.0-gammamatchglobal);
    }
  }
   double logscalefactor=rescalesm(alpha[0],2*k); //MacDonald & Zucchini fix
  // calculate alpha for the rest of the sites
  for (int s=1;s<this->sites;s++) {
    // first the probability of recombination between this site and the last
    double p_global_norec=exp(-rho*this->d[s-1]/(double(globalcount[s-1])));
    double p_local_norec=exp(-rho*this->d[s-1]/(double(popcount[location]))); 

    // std::cout << this->d.size() << " d: " << this->d[s-1] << "  pnr:  " <<   p_global_norec << " " <<  p_local_norec << std::endl;
    
    gammamatchglobal=(double(globalcount[s])+0.5*theta)/(double(globalcount[s])+theta);
    // instead of 1 term in li & Stephens we need 4
    double smlocal=0.0;
    double smglobal=0.0;
    for (int ii=0;ii<k;ii++) {
      // sum up the probabilites for locals
      if (pop[this->order[ii]]==location) smlocal+=alpha[s-1][ii];
      if (global[ii][s-1]) smglobal+=alpha[s-1][k+ii];
    }
 
    double smlocallocal=
      (popcount[location]>0)
      ?smlocal*(1.-p_local_norec)*plocalcopy/double(popcount[location])
      :0.0;

    double smgloballocal=
      (popcount[location]>0)
      ?smglobal*(1.-p_global_norec)*plocalcopy/double(popcount[location])
      :0.0;
    
    double smlocalglobal = smlocal*(1.-p_local_norec)*(1.-plocalcopy)/double(globalcount[s]);
    double smglobalglobal = smglobal*(1.-p_global_norec)*(1.-plocalcopy)/double(globalcount[s]); 
    
    for (int x=0;x<k;x++) {
      if (pop[this->order[x]]==location)
	alpha[s][x]=p_local_norec*alpha[s-1][x]+smlocallocal+smgloballocal;
      else alpha[s][x]=0.0;
      if (global[x][s])  // can't be copied globally from a local copy
	alpha[s][k+x]=p_global_norec*alpha[s-1][k+x]+smlocalglobal+smglobalglobal; 
      else alpha[s][k+x]=0.0;
      
      if (this->haps[this->order[x]][s]==this->haps[this->order[k]][s]) {
	alpha[s][x]*=gammamatchlocal;
	alpha[s][x+k] *= gammamatchglobal;
      } else {
	alpha[s][x]*=(1.-gammamatchlocal);
	alpha[s][x+k] *= (1.-gammamatchglobal);
      }
    }
    logscalefactor+=rescalesm(alpha[s],2*k); 
    //  std::cout << "lsf " << "s = " << s << " " << logscalefactor << std::endl;
  }
  return logscalefactor;
}
/**  A function to calculate the alpha matrix for sample k for (k>0)
 *   This assumes that we shall sample the global/local status of sample
 *   k rather than assuming it to be known            
 */
template <typename T>
double spact<T>::getalpha_prior(double rho, const double *F,
			  TNT::Array2D<double> &alpha,
			  const TNT::Array2D<bool> &global, 
			  const TNT::Array1D<int> &popcount,
			  const TNT::Array1D<int> &globalcount, int k) const
{
  int location = pop[this->order[k]];  // the location of haplotype k
  // the probability that a recombination is copied from a "local" chromosome
  double plocalcopy= 
    (popcount[location]>0)?F[location]*double(popcount[location])
    /(1.+(double)(popcount[location]-1.)*F[location]):0.0;
  // we can calculate the probabilities for the first site directly
  // there are no transitions and we assume exchangeability of haplotypes 
  // loop over the haplotypes from which the site could be copied
  for (int x=0;x<k;x++) {
    if (pop[this->order[x]]==location) 
      alpha[0][x] = plocalcopy/double(popcount[location]);  
    else alpha[0][x]=0.0;
    
    if (global[x][0])  // can't be copied globally from a local copy
      alpha[0][k+x]= (1.-plocalcopy)/double(globalcount[0]); 
    else alpha[0][k+x]=0.0;
  }
    double logscalefactor=rescalesm(alpha[0],2*k); //MacDonald & Zucchini fix
  // calculate alpha for the rest of the sites
  for (int s=1;s<this->sites;s++) {
    // first the probability of recombination between the last site and this
    double p_global_norec=exp(-rho*this->d[s-1]/(double(globalcount[s-1])));
    double p_local_norec=exp(-rho*this->d[s-1]/(double(popcount[location])));         
    
    // instead of 1 term in li & Stephens we need 4
    double smlocal=0.0;     // sum of all local probabilities
    double smglobal=0.0;    // sum of global probabilities

    for (int ii=0;ii<k;ii++) {
      // sum up the probabilites for locals
      if (pop[this->order[ii]]==location) {
	smlocal+=alpha[s-1][ii];
      } else 	assert(alpha[s-1][ii]<=0.0);
      if (global[ii][s-1]) {
	smglobal+=alpha[s-1][k+ii];
      } else assert(alpha[s-1][k+ii]<=0.0);
    }

    double smlocallocal=
      (popcount[location]>0)
      ?smlocal*(1.-p_local_norec)*plocalcopy/double(popcount[location])
      :0.0;
  
    double smgloballocal=
      (popcount[location]>0)
      ?smglobal*(1.-p_global_norec)*plocalcopy/double(popcount[location])
      :0.0;
    
    double smlocalglobal = 
      smlocal*(1.-p_local_norec)*(1.-plocalcopy)/double(globalcount[s]);
    double smglobalglobal = 
      smglobal*(1.-p_global_norec)*(1.-plocalcopy)/double(globalcount[s]); 

    for (int x=0;x<k;x++) {
      if (pop[this->order[x]]==location)
	alpha[s][x]=p_local_norec*alpha[s-1][x]+smlocallocal+smgloballocal;
      else alpha[s][x]=0.0;
      if (global[x][s])  // can't be copied globally from a local copy
	alpha[s][k+x]=p_global_norec*alpha[s-1][k+x]+smlocalglobal+smglobalglobal; 
      else alpha[s][k+x]=0.0;
    }
    logscalefactor+=rescalesm(alpha[s],2*k); 
  }
  return logscalefactor;
}
/**  A function to calculate the alpha matrix for sample k for (k>0)
 *   This assumes that we shall sample the global/local status of sample
 *   k rather than assuming it to be known            
 */
// template <typename T>
// double spact<T>::getalpha_priorb(double rho, const double *F,
// 			  TNT::Array2D<double> &alpha,
// 			  const TNT::Array2D<bool> &global, 
// 			  const TNT::Array1D<int> &popcount,
// 			  const TNT::Array1D<int> &globalcount, int k) const
// {
//   int location = pop[this->order[k]];  // the location of haplotype k
//   // the probability that a recombination is copied from a "local" chromosome
//   double plocalcopy= 
//     (popcount[location]>0)?F[location]*double(popcount[location])
//     /(1.+(double)(popcount[location]-1.)*F[location]):0.0;
//   // we can calculate the probabilities for the first site directly
//   // there are no transitions and we assume exchangeability of haplotypes 
//   // loop over the haplotypes from which the site could be copied
  
//   alpha[0][0]=plocalcopy/double(popcount[location]);
//   alpha[0][1]=plocalcopy - alpha[0][0];
//   alpha[0][2]=(1.-plocalcopy)/double(globalcount[0]);
//   alpha[0][3]=(1.-plocalcopy)-alpha[0][2];
  
//   // calculate alpha for the rest of the sites
//   for (int s=1;s<this->sites;s++) {
//     // first the probability of recombination between the last site and this
//     double p_global_norec=exp(-rho*this->d[s-1]/(double(globalcount[s-1])));
//     double p_local_norec=exp(-rho*this->d[s-1]/(double(popcount[location])));         

//     double smlocallocal=
//       (popcount[location]>0)
//       ?(alpha[s-1][0]+alpha[s-1][1])*(1.-p_local_norec)*plocalcopy/double(popcount[location])
//       :0.0;
  
//     double smgloballocal=
//       (popcount[location]>0)
//       ?(alpha[s-1][2]+alpha[s-1][3])*(1.-p_global_norec)*plocalcopy/double(popcount[location])
//       :0.0;
    
//     double smlocalglobal = 
//       (alpha[s-1][0]+alpha[s-1][1])*(1.-p_local_norec)*(1.-plocalcopy)/double(globalcount[s]);
//     double smglobalglobal = 
//       (alpha[s-1][2]+alpha[s-1][3])*(1.-p_global_norec)*(1.-plocalcopy)/double(globalcount[s]); 

//     for (int x=0;x<k;x++) {
//       if (pop[this->order[x]]==location)
// 	alpha[s][x]=p_local_norec*alpha[s-1][x]+smlocallocal+smgloballocal;
//       else alpha[s][x]=0.0;
//       if (global[x][s])  // can't be copied globally from a local copy
// 	alpha[s][k+x]=p_global_norec*alpha[s-1][k+x]+smlocalglobal+smglobalglobal; 
//       else alpha[s][k+x]=0.0;
//     }
//     logscalefactor+=rescalesm(alpha[s],2*k); 
//   }
//   return logscalefactor;
// }
#endif
/**************************************************************************************************
*
*
*
*
*
*
*
*
*
*
*
*
*
*
*
*
****************************************************************************************************/
/** calculate the likelihood for one rearrangement of the haplotypes 
 * which is assumed to have already happened.    We also have to sample 
 * (as we go) whether or not each position is in the "local" or global population */
// template <typename T>
// double spact<T>::calc(double rho, const TNT::Array1D<double> & F, TNT::Array2D<int> &path
// 		       , TNT::Array2D<bool> &global,bool resamplepath)
// {
//   // path holds the path through the haplotypes, we might as well keep it 
//   // as it costs little.  Global holds the copy type-local or from the global pool 
//   // alpha[x,s] defined as the probability of sampling the first s sites of
//   // your haplotype, and of sampling haplotype x. The first k hold the 
//   // probabilites of a local copy of the site, the last k copying from the global pool.
//   // Some (or most) of these will be zero (zero locally and globally if they are 
//   // local to  another population, zero globally if they are local to the current
//   // population, only non zero for both if they are both local and in the global pool)
//   TNT::Array2D<double> alpha(sites,2*n-2);
//   std::vector<int> popcount(npops_,0);     // the count in each of the populations
//   TNT::Array1D<int> globalcount(sites,1); // the count in the "global" pool
  
//   popcount[pop[order[0]]] +=1;           // add one to the local count
//   for (int s=0;s<sites;s++) {
//     global[0][s]=true;
//     path[0][s]=0;
//   }

//   double lprob=0.0;                       // the first sample has probability 1
//   // now consider each of the subsequent samples in turn.
//   // Because my arrays are 0 offset haplotype k is the k+1 st sample
//   // so use the same k as in Li and Stephens' appendix
//   for (int k=1;k<n;k++) {
//     int location = pop[order[k]];  // the location of the haplotype under consideration
//     // what is the probability that a recombination is copied from 
//     // a "local" chromosome?
//     double plocalcopy=0.0;  // value if this is the first haplotype from this subpopulation
//     if (popcount[location]>0) 
//       plocalcopy=F[location]*double(popcount[location])/(1.+(double)(popcount[location]-1.)*F[location]);
// #ifdef LOUD
//     std::cout << "plocalcopy = " << plocalcopy << " local popsize " 
// 	      << popcount[location] <<  " location " << location << std::endl;
// #endif

//     // what are the match probabilities for the global population
//     double gammamatchglobal=double(globalcount[0])/(double(globalcount[0])+thetahat)
//       +0.5*thetahat/(double(globalcount[0])+thetahat);
//     double gammamissglobal=1.- gammamatchglobal;//   0.5*thetahat/(double(globalcount[0])+thetahat);

// #ifdef PERFECTLOCALMATCH 
//     double gammamatchlocal=1.0;       // constant values througout
//     double gammamisslocal=0.0;
// #else
//     // get some values for local matches and misses - use k here (v. dodgy) 
//     double gammamatchlocal=k/(k+thetahat) + 0.5*thetahat/(k+thetahat);;
//     double gammamisslocal= 1.-gammamatchlocal;          // = 0.5*thetahat/(k+thetahat);
// #endif

//     // we can calculate the probabilities for the first site directly
//     // there are no transitions and we assume exchangeability of haplotypes 
//     // loop over the haplotypes from which the site could be copied
//     for (int x=0;x<k;x++) {
//       if (haps[order[x]][0]==haps[order[k]][0]) {
// 	if (pop[order[x]]==location) 
// 	  alpha[0][x] = plocalcopy*gammamatchlocal/double(popcount[location]);  
// 	else alpha[0][x]=0.0;
// 	if (global[x][0]) 
// 	  alpha[0][k+x] = (1.-plocalcopy)*gammamatchglobal/double(globalcount[0]);
// 	else alpha[0][k+x]=0.0;
//       } else {
// 	if (pop[order[x]]==location) 
// 	  alpha[0][x] = plocalcopy*gammamisslocal/double(popcount[location]);
// 	else alpha[0][x]=0.0;
// 	if (global[x][0]) 
// 	  alpha[0][k+x] = (1.-plocalcopy)*gammamissglobal/double(globalcount[0]);
// 	else alpha[0][k+x]=0.0;
//       }
//     }
//     double logscalefactor=rescale(alpha[0],2*k); //MacDonald & Zucchini fix
//     // calculate alpha for the rest of the sites
//     for (int s=1;s<sites;s++) {
//       // first the probability of recombination.  I am not sure 
//       // about the use of localcopies here but it seems as good as anything
//       // justification is that you can recombine with anything in the same
//       // population, and all the globals (but these may include ones from the 
//       // same population twice) which is a problem
      
//       double scaled_rec_rate=double(globalcount[s-1]);//+double(popcount[location]);
//       double p_local_norec=exp(-rho*d[s-1]/scaled_rec_rate);
//       double p_global_norec=exp(-rho*d[s-1]/double(k));         //?!
   
//       gammamatchglobal=double(globalcount[s])/(double(globalcount[s])+thetahat)
// 	+0.5*thetahat/(double(globalcount[s])+thetahat);
//       gammamissglobal=0.5*thetahat/(double(globalcount[s])+thetahat);

//       assert(p_local_norec>=0.0&&p_local_norec<=1.0);
//       assert(p_global_norec>=0.0&&p_global_norec<=1.0);
//       // instead of 1 term in li & Stephens we need 4
//       double smlocal=0.0;
//       double smglobal=0.0;
//       for (int ii=0;ii<k;ii++) {
// 	// sum up the probabilites for locals
// 	if (pop[order[ii]]==location) smlocal+=alpha[s-1][ii];
// 	if (global[ii][s-1]) smglobal+=alpha[s-1][k+ii];
//       }
//       double smlocallocal=0.0;
//       double smgloballocal=0.0;
//       if (popcount[location]>0) 
// 	smlocallocal =smlocal*(1.-p_local_norec)*plocalcopy/double(popcount[location]);
     
//       double smlocalglobal = smlocal*(1.-p_local_norec)*(1.-plocalcopy)/double(globalcount[s]);
//       double smglobalglobal = smglobal*(1.-p_global_norec)*(1.-plocalcopy)/double(globalcount[s]); 
//       if (popcount[location]>0) 
// 	smgloballocal = smglobal*(1.-p_global_norec)*plocalcopy/double(popcount[location]);

//       for (int x=0;x<k;x++) {
// 	if (pop[order[x]]==location)
// 	  alpha[s][x]=p_local_norec*alpha[s-1][x]+smlocallocal+smgloballocal;
// 	else alpha[s][x]=0.0;
	
// 	if (global[x][s])
// 	  alpha[s][k+x]=p_global_norec*alpha[s-1][k+x]+smlocalglobal+smglobalglobal; 
// 	else alpha[s][k+x]=0.0;

// 	if (haps[order[x]][s]==haps[order[k]][s]) {
// 	  alpha[s][x]*=gammamatchlocal;
// 	  alpha[s][x+k] *= gammamatchglobal;
// 	} else {
// 	  alpha[s][x]*=gammamisslocal;
// 	  alpha[s][x+k] *= gammamissglobal;
// 	}
//       }
//       logscalefactor+=rescale(alpha[s],2*k); // The scale factor of MacDonald and Zucchini
//     }

//     // calculate the prob for haplotype k by summing over all the 
//     // paths to the last site
//     double prob=std::accumulate(alpha[sites-1],alpha[sites-1]+2*k,0.0);
//     lprob+=log(prob)+logscalefactor;
 
//     if (log(prob)+logscalefactor>0.0) {
//       std::cerr << "problem, added a probability > 1 " << lprob 
// 		<< " log(prob) = " << log(prob) << " logscalefactor = " << logscalefactor 
// 		<< " added " << log(prob)+logscalefactor << std::endl;
//       exit(EXIT_FAILURE);
//     }

//     // Method as in Boys et al 2000 - the last site is "site-1".  sample from the "equilibrium"
//     path[k][sites-1] = gen_from_p(alpha[sites-1],alpha[sites-1]+2*k,r);
//     // need to calculate the transition probabilities here
//     std::vector<double> lp(2*k);  
//     for (int s=sites-2;s>=0;s--) {
//       double scaled_rec_rate=globalcount[s];
//       double p_local_norec=exp(-rho*d[s]/scaled_rec_rate);
//       double p_global_norec=exp(-rho*d[s]/(double(k)));//??!!
//       int xprime=path[k][s+1];
//       for (int x=0;x<k;x++) { // loop throught all the possible haplotypes
// 	if (xprime<k) { // current path (ahead) local
// 	  if (popcount[location]==0) {
// 	    std::cerr << "site " << s+1 << " xprime = " << xprime << std::endl;
// 	  }
// 	  if (xprime==x) {       // staying on the same path
// 	    lp[x] = alpha[s][x]*
// 	      (p_local_norec+(1.-p_local_norec)*plocalcopy/double(popcount[location]));
// 	  } else if (popcount[location]>0) {               // move
// 	    lp[x] = alpha[s][x]*
// 	      (1.-p_local_norec)*plocalcopy/double(popcount[location]);
// 	  } else lp[x]=0.0;
// 	  lp[k+x] = alpha[s][k+x]*
// 	    (1.-p_local_norec)*(1.-plocalcopy)/double(globalcount[s+1]); 
// 	} else {  // current path global
// 	  if (!global[xprime-k][s+1]) {
// 	    std::cout << "site " << s+1 << " xprime = " << xprime 
// 		      << " alpha[s+1][xprime] "<< alpha[s+1][xprime] 
// 		      << " alpha[s][xprime] "<< alpha[s][xprime] << std::endl;
// 	    exit(EXIT_FAILURE);
// 	  }
// 	  if (xprime==k+x) // copy from x as global
// 	    lp[k+x] = alpha[s][k+x]*
// 	      (p_global_norec+(1.-p_global_norec)*(1.-plocalcopy)/double(globalcount[s+1]));
// 	  else 
// 	    lp[k+x] = alpha[s][k+x]*
// 	      (1.-p_global_norec)*(1.-plocalcopy)/double(globalcount[s+1]);
// 	  if (popcount[location]>0)
// 	    lp[x] = alpha[s][x]*
// 	      (1.-p_global_norec)*plocalcopy/double(popcount[location]);
// 	  else lp[x]=0.0;
// 	}
//       }
//       path[k][s]=gen_from_pr(lp,2*k,r);
//     }
//     for (int ii=0;ii<sites;ii++) {
//       if (path[k][ii]>=k) {
// 	global[k][ii]=true;
// 	globalcount[ii]+=1;
//       } else global[k][ii]=false;
//     }
//     popcount[location]+=1;
//   }
  
  
//   return lprob;
// }
// /** calculate the likelihood for one rearrangement of the haplotypes 
//  * which is assumed to have already happened.    We also have to sample 
//  * (as we go) whether or not each position is in the "local" or global population */
// template <typename T>
// double spact<T>::sample(double rho, const TNT::Array1D<double> & F, TNT::Array2D<bool> &global)
// {
//   // path holds the path through the haplotypes, we might as well keep it 
//   // as it costs little.  Global holds the copy type-local or from the global pool 
//   // alpha[x,s] defined as the probability of sampling the first s sites of
//   // your haplotype, and of sampling haplotype x. The first k hold the 
//   // probabilites of a local copy of the site, the last k copying from the global pool.
//   // Some (or most) of these will be zero (zero locally and globally if they are 
//   // local to  another population, zero globally if they are local to the current
//   // population, only non zero for both if they are both local and in the global pool)
//   TNT::Array2D<double> alpha(sites,2*n-2);
//   std::vector<int> popcount(npops_,0);     // the count in each of the populations
//   TNT::Array1D<int> globalcount(sites,1); // the count in the "global" pool
//   TNT::Array2D<int> path(n,sites,0);
//   double lpsample=0.0;
//   double lprob=0.0;
  
//   popcount[pop[order[0]]] +=1;           // add one to the local count
//   for (int s=0;s<sites;s++)  global[0][s]=true;

//   double lprob=0.0;      // the first sample has probability 1
//   // now consider each of the subsequent samples in turn.
//   // Because my arrays are 0 offset haplotype k is the k+1 st sample
//   // so use the same k as in Li and Stephens' appendix
//   for (int k=1;k<n;k++) {
//     double lscalefactor = getalpha(rho,alpha,gammamatch,k,missingvals);
//     // calculate the prob for haplotype k by summing over all the 
//     // paths to the last site
//     double prob=std::accumulate(alpha[sites-1],alpha[sites-1]+2*k,0.0);
//     lprob+=log(prob)+lscalefactor;
    
//     // Method as in Boys et al 2000 - the last site is "site-1".  sample from the "equilibrium"
//     double psample;
//     path[k][sites-1] = gen_from_p(alpha[sites-1],alpha[sites-1]+2*k,r,psample);
//     lpsample+=log(psample);
//     // need to calculate the transition probabilities here
//     std::vector<double> lp(2*k);  
//     for (int s=sites-2;s>=0;s--) {
//       double scaled_rec_rate=globalcount[s];
//       double p_local_norec=exp(-rho*d[s]/(double(popcount[location])));
//       double p_global_norec=exp(-rho*d[s]/(double(globalcount[s])));
//       int xprime=path[k][s+1];
//       for (int x=0;x<k;x++) { // loop throught all the possible haplotypes
// 	if (xprime<k) { // current path (ahead) local
// 	  if (popcount[location]==0) {
// 	    std::cerr << "site " << s+1 << " xprime = " << xprime << std::endl;
// 	  }
// 	  if (xprime==x) {       // staying on the same path
// 	    lp[x] = alpha[s][x]*
// 	      (p_local_norec+(1.-p_local_norec)*plocalcopy/double(popcount[location]));
// 	  } else if (popcount[location]>0) {               // move
// 	    lp[x] = alpha[s][x]*
// 	      (1.-p_local_norec)*plocalcopy/double(popcount[location]);
// 	  } else lp[x]=0.0;
// 	  lp[k+x] = alpha[s][k+x]*
// 	    (1.-p_local_norec)*(1.-plocalcopy)/double(globalcount[s+1]); 
// 	} else {  // current path global
// 	  if (!global[xprime-k][s+1]) {
// 	    std::cout << "site " << s+1 << " xprime = " << xprime 
// 		      << " alpha[s+1][xprime] "<< alpha[s+1][xprime] 
// 		      << " alpha[s][xprime] "<< alpha[s][xprime] << std::endl;
// 	    exit(EXIT_FAILURE);
// 	  }
// 	  if (xprime==k+x) // copy from x as global
// 	    lp[k+x] = alpha[s][k+x]*
// 	      (p_global_norec+(1.-p_global_norec)*(1.-plocalcopy)/double(globalcount[s+1]));
// 	  else 
// 	    lp[k+x] = alpha[s][k+x]*
// 	      (1.-p_global_norec)*(1.-plocalcopy)/double(globalcount[s+1]);
// 	  if (popcount[location]>0)
// 	    lp[x] = alpha[s][x]*
// 	      (1.-p_global_norec)*plocalcopy/double(popcount[location]);
// 	  else lp[x]=0.0;
// 	}
//       }
//       path[k][s]=gen_from_pr(lp,2*k,r,psample);
//       lpsample += psample;
//     }
//     for (int ii=0;ii<sites;ii++) {
//       if (path[k][ii]>=k) {
// 	global[k][ii]=true;
// 	globalcount[ii]+=1;
//       } else global[k][ii]=false;
//     }
//     popcount[location]+=1;
//   }
//   return lprob;
// }
/** calculate the likelihood for one rearrangement of the haplotypes 
 * which is assumed to have already happened.    We also have to sample 
   * (as we go) whether or not each position is in the "local" or global population */
// template <typename T>
//   double spact<T>::lpsample(double rho, const TNT::Array1D<double> & F
// 			 , const TNT::Array2D<bool> &global)
// {
//   double lpsamp=0.0;
//   int location = pop[order[k]];  // the location of the haplotype under consideration
//     // what is the probability that a recombination is copied from 
//     // a "local" chromosome?
//     double plocalcopy=0.0;  // value if this is the first haplotype from this subpopulation
//     if (popcount[location]>0) 
//       plocalcopy=F[location]*double(popcount[location])/(1.+(double)(popcount[location]-1.)*F[location]);
//     // what are the match probabilities for the global population
//     double gammamatchglobal=double(globalcount[0])/(double(globalcount[0])+thetahat)
//       +0.5*thetahat/(double(globalcount[0])+thetahat);
//     double gammamissglobal=1.- gammamatchglobal;//   0.5*thetahat/(double(globalcount[0])+thetahat);

// #ifdef PERFECTLOCALMATCH 
//     double gammamatchlocal=1.0;       // constant values througout
//     double gammamisslocal=0.0;
// #else
//     // get some values for local matches and misses - use k here (v. dodgy) 
//     double gammamatchlocal=k/(k+thetahat) + 0.5*thetahat/(k+thetahat);;
//     double gammamisslocal= 1.-gammamatchlocal;          // = 0.5*thetahat/(k+thetahat);
// #endif

//     // we can calculate the probabilities for the first site directly
//     // there are no transitions and we assume exchangeability of haplotypes 
//     // loop over the haplotypes from which the site could be copied
//     for (int x=0;x<k;x++) {
//       if (haps[order[x]][0]==haps[order[k]][0]) {
// 	if (pop[order[x]]==location) 
// 	  alpha[0][x] = plocalcopy*gammamatchlocal/double(popcount[location]);  
// 	else alpha[0][x]=0.0;
// 	if (global[x][0]) 
// 	  alpha[0][k+x] = (1.-plocalcopy)*gammamatchglobal/double(globalcount[0]);
// 	else alpha[0][k+x]=0.0;
//       } else {
// 	if (pop[order[x]]==location) 
// 	  alpha[0][x] = plocalcopy*(1.-gammamatchlocal)/double(popcount[location]);
// 	else alpha[0][x]=0.0;
// 	if (global[x][0]) 
// 	  alpha[0][k+x] = (1.-plocalcopy)*(1.-gammamatchglobal)/double(globalcount[0]);
// 	else alpha[0][k+x]=0.0;
//       }
//     }
//     double logscalefactor=rescale(alpha[0],2*k); //MacDonald & Zucchini fix
//     // calculate alpha for the rest of the sites
//     for (int s=1;s<sites;s++) {
//       // first the probability of recombination.  I am not sure 
//       // about the use of localcopies here but it seems as good as anything
//       // justification is that you can recombine with anything in the same
//       // population, and all the globals (but these may include ones from the 
//       // same population twice) which is a problem
      
//       double scaled_rec_rate=double(globalcount[s-1]);//+double(popcount[location]);
//       double p_local_norec=exp(-rho*d[s-1]/scaled_rec_rate);
//       double p_global_norec=exp(-rho*d[s-1]/double(k));         //?!
   
//       gammamatchglobal=double(globalcount[s])/(double(globalcount[s])+thetahat)
// 	+0.5*thetahat/(double(globalcount[s])+thetahat);
//       gammamissglobal=0.5*thetahat/(double(globalcount[s])+thetahat);

//       assert(p_local_norec>=0.0&&p_local_norec<=1.0);
//       assert(p_global_norec>=0.0&&p_global_norec<=1.0);
//       // instead of 1 term in li & Stephens we need 4
//       double smlocal=0.0;
//       double smglobal=0.0;
//       for (int ii=0;ii<k;ii++) {
// 	// sum up the probabilites for locals
// 	if (pop[order[ii]]==location) smlocal+=alpha[s-1][ii];
// 	if (global[ii][s-1]) smglobal+=alpha[s-1][k+ii];
//       }
//       double smlocallocal=0.0;
//       double smgloballocal=0.0;
//       if (popcount[location]>0) 
// 	smlocallocal =smlocal*(1.-p_local_norec)*plocalcopy/double(popcount[location]);
     
//       double smlocalglobal = smlocal*(1.-p_local_norec)*(1.-plocalcopy)/double(globalcount[s]);
//       double smglobalglobal = smglobal*(1.-p_global_norec)*(1.-plocalcopy)/double(globalcount[s]); 
//       if (popcount[location]>0) 
// 	smgloballocal = smglobal*(1.-p_global_norec)*plocalcopy/double(popcount[location]);

//       for (int x=0;x<k;x++) {
// 	if (pop[order[x]]==location)
// 	  alpha[s][x]=p_local_norec*alpha[s-1][x]+smlocallocal+smgloballocal;
// 	else alpha[s][x]=0.0;
	
// 	if (global[x][s])
// 	  alpha[s][k+x]=p_global_norec*alpha[s-1][k+x]+smlocalglobal+smglobalglobal; 
// 	else alpha[s][k+x]=0.0;

// 	if (haps[order[x]][s]==haps[order[k]][s]) {
// 	  alpha[s][x]*=gammamatchlocal;
// 	  alpha[s][x+k] *= gammamatchglobal;
// 	} else {
// 	  alpha[s][x]*=gammamisslocal;
// 	  alpha[s][x+k] *= gammamissglobal;
// 	}
//       }
//       logscalefactor+=rescale(alpha[s],2*k); // The scale factor of MacDonald and Zucchini
//     }

        

//     // Method as in Boys et al 2000 - the last site is "site-1".  sample from the "equilibrium"
//     double psample;
//     path[k][sites-1] = gen_from_p(alpha[sites-1],alpha[sites-1]+2*k,r,psample);
//     lpsample+=log(psample);
//     // need to calculate the transition probabilities here
//     std::vector<double> lp(2*k);  
//     for (int s=sites-2;s>=0;s--) {
//       double scaled_rec_rate=globalcount[s];
//       double p_local_norec=exp(-rho*d[s]/scaled_rec_rate);
//       double p_global_norec=exp(-rho*d[s]/(double(k)));//??!!
//       int xprime=path[k][s+1];
//       for (int x=0;x<k;x++) { // loop throught all the possible haplotypes
// 	if (xprime<k) { // current path (ahead) local
// 	  if (popcount[location]==0) {
// 	    std::cerr << "site " << s+1 << " xprime = " << xprime << std::endl;
// 	  }
// 	  if (xprime==x) {       // staying on the same path
// 	    lp[x] = alpha[s][x]*
// 	      (p_local_norec+(1.-p_local_norec)*plocalcopy/double(popcount[location]));
// 	  } else if (popcount[location]>0) {               // move
// 	    lp[x] = alpha[s][x]*
// 	      (1.-p_local_norec)*plocalcopy/double(popcount[location]);
// 	  } else lp[x]=0.0;
// 	  lp[k+x] = alpha[s][k+x]*
// 	    (1.-p_local_norec)*(1.-plocalcopy)/double(globalcount[s+1]); 
// 	} else {  // current path global
// 	  if (!global[xprime-k][s+1]) {
// 	    std::cout << "site " << s+1 << " xprime = " << xprime 
// 		      << " alpha[s+1][xprime] "<< alpha[s+1][xprime] 
// 		      << " alpha[s][xprime] "<< alpha[s][xprime] << std::endl;
// 	    exit(EXIT_FAILURE);
// 	  }
// 	  if (xprime==k+x) // copy from x as global
// 	    lp[k+x] = alpha[s][k+x]*
// 	      (p_global_norec+(1.-p_global_norec)*(1.-plocalcopy)/double(globalcount[s+1]));
// 	  else 
// 	    lp[k+x] = alpha[s][k+x]*
// 	      (1.-p_global_norec)*(1.-plocalcopy)/double(globalcount[s+1]);
// 	  if (popcount[location]>0)
// 	    lp[x] = alpha[s][x]*
// 	      (1.-p_global_norec)*plocalcopy/double(popcount[location]);
// 	  else lp[x]=0.0;
// 	}
//       }
//       path[k][s]=gen_from_pr(lp,2*k,r,psample);
//       lpsample += psample;
//     }
//     for (int ii=0;ii<sites;ii++) {
//       if (path[k][ii]>=k) {
// 	global[k][ii]=true;
// 	globalcount[ii]+=1;
//       } else global[k][ii]=false;
//     }
//     popcount[location]+=1;
//   }
//   return lprob;
// }

/** calculate the likelihood for one rearrangement of the haplotypes conditional
 * on the globalness - should be no dependency on F 
 */
// template <typename T>
//   double spact<T>::likelihood(double theta, TNT::Array2D<bool> &global)
// {
//   // alpha[x,s] defined as the probability of sampling the first s sites of
//   // your haplotype, and of sampling haplotype x.
//   TNT::Array2D<double> alpha(sites,n-1);
//   std::vector<int> popcount(npops_,0);     // the count in each of the populations
//   TNT::Array1D<int> globalcount(sites,1);  // the count in the "global" pool
  
//   popcount[pop[order[0]]] +=1;             // add one to the local count
  
//   double lprob=0.0;                        // the first sample has probability 1
//   // now consider each of the subsequent samples in turn.
//   // Because my arrays are 0 offset haplotype k is the k+1 st sample
//   // so use the same k as in Li and Stephens' appendix
//   for (int k=1;k<n;k++) {
//     int location = pop[order[k]];  // the location of the haplotype under consideration
//     // what is the probability that a recombination is copied from a "local" chromosome?
//     double plocalcopy=
//       F[location]*double(popcount[location])/(1.+(double)(popcount[location]-1.)*F[location]);
//     // what are the match probabilities for the global population
//     double gammamatchglobal=double(globalcount[0])/(double(globalcount[0])+theta)
//       +0.5*theta/(double(globalcount[0])+theta);
//     double gammamatchlocal = 1.0;       
//     // we can calculate the probabilities for the first site directly
//     // there are no transitions and we assume exchangeability of haplotypes 
//     // loop over the haplotypes from which the site could be copied
//     if (global[k][0]) {
//       for (int x=0;x<k;x++) {
// 	if (global[x]) {
// 	  if (haps[hapsder[x]][0]==haps[order[k]][0])   // matches 
// 	    alpha[0][x] = gammamatchglobal/(double(globalcount[0]));
// 	  else  
// 	    alpha[0][x] = (1.-gammamatchglobal)/(double(globalcount[0]));
// 	} else alpha[0][x]=0.0;
//       }
//     } else {  // local copy
//       for (int x=0;x<k;x++) {
// 	if (pop[order[x]]==location) {
// 	  if (haps[order[x]][0]==haps[order[k]][0]) 
// 	    alpha[0][x]=1.0;
// 	  else alpha[0][x]=0.0;
// 	} else alpha[0][x]=0.0;
//       }
//     }
//     double logscalefactor=rescale(alpha[0],2*k); //MacDonald & Zucchini fix
//     // calculate alpha for the rest of the sites
//     for (int s=1;s<sites;s++) {
//       // first the probability of recombination. 
//       // These values can be justified by reference to the ESD
//       double p_local_norec=exp(-rho*d[s-1]/double(globalcount[s-1]));
//       double p_global_norec=exp(-rho*d[s-1]/double(popcount[location]));        
//       gammamatchglobal=double(globalcount[s])/(double(globalcount[s])+theta)
// 	+0.5*theta/(double(globalcount[s])+theta);
//       // now calculate the second term from Li & Stephens
//       // not quite as simple as we have to work out every combination 
//       // copy from local and global 
//       double smfromlocal=0.0,smfromglobal=0.0;
//       if (global[k][s]) {         // copying to a global 
// 	for (int ii=0;ii<k;ii++) {
// 	  if (global[x][s-1]) {     // both global
// 	    smfromglobal += alpha[s-1][ii]*(1.- p_global_norec);
// 	  } else {                       // target global copied from local
// 	    if (location==pop[order[x]]) // must be the same location
// 	      smfromlocal += alpha[s-1][ii]*(1.- p_local_norec);
// 	  }
// 	}
// 	// now calculate the alpha values from each of the difference chromosomes
// 	for (int x=0;x<k;x++) {
// 	  if (global[x][s]) {    
// 	    alpha[s][x] = p_global_norec*alpha[s-1][x] 
// 	      + (1.-p_global_norec)*smfromglobal*(1.-plocalcopy)/double(globalcount[s]);;

// 	    if (haps[order[x]][0]==haps[order[k]][0])   // matches 
// 	      alpha[s][x]*=gammamatchglobal;
// 	    else alpha[s][x]*=(1.-gammamatchglobal);
// 	  } else {
// 	    alpha[s][x] = p_local_norec*alpha[s-1][x] 
// 	      + (1.-p_local_norec)*smfromlocal*(1.-plocalcopy)/double(globalcount[s]);;
// 	  }
// 	    else 
// 	      alpha[s][x]= ;
// 	    if (pop[order[x]]==location)
// 	    alpha[s][x]=p_local_norec*alpha[s-1][x]+sm;
// 	  else alpha[s][x]=0.0;
	  
// 	  if (global[x][s])
// 	  alpha[s][k+x]=p_global_norec*alpha[s-1][k+x]+smlocalglobal+smglobalglobal; 
// 	else alpha[s][k+x]=0.0;

// 	if (haps[order[x]][s]==h[haps[order[k]][s]) {
// 	  alpha[s][x]*=gammamatchlocal;
// 	  alpha[s][x+k] *= gammamatchglobal;
// 	} else {
// 	  alpha[s][x]*=gammamisslocal;
// 	  alpha[s][x+k] *= gammamissglobal;
// 	}
//       } else {
// 	for (int ii=0;ii<k;ii++) {
// 	  if (global[x][s]) {  // target local copied from global
// 	    smfromglobal += alpha[s-1][ii]*(1.- p_global_norec);
// 	  } else {                    // both local
// 	    smfromlocal += alpha[s-1][ii]*(1.- p_local_norec);
// 	  }
// 	}
//       }
   
//       for (int x=0;x<k;x++) {
// 	if (pop[order[x]]==location)
// 	  alpha[s][x]=p_local_norec*alpha[s-1][x]+sm;
// 	else alpha[s][x]=0.0;
	
// 	if (global[x][s])
// 	  alpha[s][k+x]=p_global_norec*alpha[s-1][k+x]+smlocalglobal+smglobalglobal; 
// 	else alpha[s][k+x]=0.0;

// 	if (haps[order[x]][s]==h[haps[order[k]][s]) {
// 	  alpha[s][x]*=gammamatchlocal;
// 	  alpha[s][x+k] *= gammamatchglobal;
// 	} else {
// 	  alpha[s][x]*=gammamisslocal;
// 	  alpha[s][x+k] *= gammamissglobal;
// 	}
//       }
//       logscalefactor+=rescale(alpha[s],2*k); // The scale factor of MacDonald and Zucchini
//     }

//     // calculate the prob for haplotype k by summing over all the 
//     // paths to the last site
//     double prob=std::accumulate(alpha[sites-1],alpha[sites-1]+2*k,0.0);
//     lprob+=log(prob)+logscalefactor;
 
//     if (log(prob)+logscalefactor>0.0) {
//       std::cerr << "problem, added a probability > 1 " << lprob 
// 		<< " log(prob) = " << log(prob) << " logscalefactor = " << logscalefactor 
// 		<< " added " << log(prob)+logscalefactor << std::endl;
//       exit(EXIT_FAILURE);
//     }
    
//     return lprob;
// 
/** calculate the likelihood for one rearrangement of the haplotypes 
 * which is assumed to have already happened.    We also have to sample 
 * (as we go) whether or not each position is in the "local" or global population */
// template <typename T>
//   double spact<T>::lpsample(double rho, const TNT::Array1D<double> & F
// 			 , const TNT::Array2D<bool> &global)
// {
//   double lpsamp=0.0;
//   int location = pop[order[k]];  // the location of the haplotype under consideration
//     // what is the probability that a recombination is copied from 
//     // a "local" chromosome?
//     double plocalcopy=0.0;  // value if this is the first haplotype from this subpopulation
//     if (popcount[location]>0) 
//       plocalcopy=F[location]*double(popcount[location])/(1.+(double)(popcount[location]-1.)*F[location]);
//     // what are the match probabilities for the global population
//     double gammamatchglobal=double(globalcount[0])/(double(globalcount[0])+thetahat)
//       +0.5*thetahat/(double(globalcount[0])+thetahat);
//     double gammamissglobal=1.- gammamatchglobal;//   0.5*thetahat/(double(globalcount[0])+thetahat);

// #ifdef PERFECTLOCALMATCH 
//     double gammamatchlocal=1.0;       // constant values througout
//     double gammamisslocal=0.0;
// #else
//     // get some values for local matches and misses - use k here (v. dodgy) 
//     double gammamatchlocal=k/(k+thetahat) + 0.5*thetahat/(k+thetahat);;
//     double gammamisslocal= 1.-gammamatchlocal;          // = 0.5*thetahat/(k+thetahat);
// #endif

//     // we can calculate the probabilities for the first site directly
//     // there are no transitions and we assume exchangeability of haplotypes 
//     // loop over the haplotypes from which the site could be copied
//     for (int x=0;x<k;x++) {
//       if (haps[order[x]][0]==haps[order[k]][0]) {
// 	if (pop[order[x]]==location) 
// 	  alpha[0][x] = plocalcopy*gammamatchlocal/double(popcount[location]);  
// 	else alpha[0][x]=0.0;
// 	if (global[x][0]) 
// 	  alpha[0][k+x] = (1.-plocalcopy)*gammamatchglobal/double(globalcount[0]);
// 	else alpha[0][k+x]=0.0;
//       } else {
// 	if (pop[order[x]]==location) 
// 	  alpha[0][x] = plocalcopy*(1.-gammamatchlocal)/double(popcount[location]);
// 	else alpha[0][x]=0.0;
// 	if (global[x][0]) 
// 	  alpha[0][k+x] = (1.-plocalcopy)*(1.-gammamatchglobal)/double(globalcount[0]);
// 	else alpha[0][k+x]=0.0;
//       }
//     }
//     double logscalefactor=rescale(alpha[0],2*k); //MacDonald & Zucchini fix
//     // calculate alpha for the rest of the sites
//     for (int s=1;s<sites;s++) {
//       // first the probability of recombination.  I am not sure 
//       // about the use of localcopies here but it seems as good as anything
//       // justification is that you can recombine with anything in the same
//       // population, and all the globals (but these may include ones from the 
//       // same population twice) which is a problem
      
//       double scaled_rec_rate=double(globalcount[s-1]);//+double(popcount[location]);
//       double p_local_norec=exp(-rho*d[s-1]/scaled_rec_rate);
//       double p_global_norec=exp(-rho*d[s-1]/double(k));         //?!
   
//       gammamatchglobal=double(globalcount[s])/(double(globalcount[s])+thetahat)
// 	+0.5*thetahat/(double(globalcount[s])+thetahat);
//       gammamissglobal=0.5*thetahat/(double(globalcount[s])+thetahat);

//       assert(p_local_norec>=0.0&&p_local_norec<=1.0);
//       assert(p_global_norec>=0.0&&p_global_norec<=1.0);
//       // instead of 1 term in li & Stephens we need 4
//       double smlocal=0.0;
//       double smglobal=0.0;
//       for (int ii=0;ii<k;ii++) {
// 	// sum up the probabilites for locals
// 	if (pop[order[ii]]==location) smlocal+=alpha[s-1][ii];
// 	if (global[ii][s-1]) smglobal+=alpha[s-1][k+ii];
//       }
//       double smlocallocal=0.0;
//       double smgloballocal=0.0;
//       if (popcount[location]>0) 
// 	smlocallocal =smlocal*(1.-p_local_norec)*plocalcopy/double(popcount[location]);
     
//       double smlocalglobal = smlocal*(1.-p_local_norec)*(1.-plocalcopy)/double(globalcount[s]);
//       double smglobalglobal = smglobal*(1.-p_global_norec)*(1.-plocalcopy)/double(globalcount[s]); 
//       if (popcount[location]>0) 
// 	smgloballocal = smglobal*(1.-p_global_norec)*plocalcopy/double(popcount[location]);

//       for (int x=0;x<k;x++) {
// 	if (pop[order[x]]==location)
// 	  alpha[s][x]=p_local_norec*alpha[s-1][x]+smlocallocal+smgloballocal;
// 	else alpha[s][x]=0.0;
	
// 	if (global[x][s])
// 	  alpha[s][k+x]=p_global_norec*alpha[s-1][k+x]+smlocalglobal+smglobalglobal; 
// 	else alpha[s][k+x]=0.0;

// 	if (haps[order[x]][s]==haps[order[k]][s]) {
// 	  alpha[s][x]*=gammamatchlocal;
// 	  alpha[s][x+k] *= gammamatchglobal;
// 	} else {
// 	  alpha[s][x]*=gammamisslocal;
// 	  alpha[s][x+k] *= gammamissglobal;
// 	}
//       }
//       logscalefactor+=rescale(alpha[s],2*k); // The scale factor of MacDonald and Zucchini
//     }

        

//     // Method as in Boys et al 2000 - the last site is "site-1".  sample from the "equilibrium"
//     double psample;
//     path[k][sites-1] = gen_from_p(alpha[sites-1],alpha[sites-1]+2*k,r,psample);
//     lpsample+=log(psample);
//     // need to calculate the transition probabilities here
//     std::vector<double> lp(2*k);  
//     for (int s=sites-2;s>=0;s--) {
//       double scaled_rec_rate=globalcount[s];
//       double p_local_norec=exp(-rho*d[s]/scaled_rec_rate);
//       double p_global_norec=exp(-rho*d[s]/(double(k)));//??!!
//       int xprime=path[k][s+1];
//       for (int x=0;x<k;x++) { // loop throught all the possible haplotypes
// 	if (xprime<k) { // current path (ahead) local
// 	  if (popcount[location]==0) {
// 	    std::cerr << "site " << s+1 << " xprime = " << xprime << std::endl;
// 	  }
// 	  if (xprime==x) {       // staying on the same path
// 	    lp[x] = alpha[s][x]*
// 	      (p_local_norec+(1.-p_local_norec)*plocalcopy/double(popcount[location]));
// 	  } else if (popcount[location]>0) {               // move
// 	    lp[x] = alpha[s][x]*
// 	      (1.-p_local_norec)*plocalcopy/double(popcount[location]);
// 	  } else lp[x]=0.0;
// 	  lp[k+x] = alpha[s][k+x]*
// 	    (1.-p_local_norec)*(1.-plocalcopy)/double(globalcount[s+1]); 
// 	} else {  // current path global
// 	  if (!global[xprime-k][s+1]) {
// 	    std::cout << "site " << s+1 << " xprime = " << xprime 
// 		      << " alpha[s+1][xprime] "<< alpha[s+1][xprime] 
// 		      << " alpha[s][xprime] "<< alpha[s][xprime] << std::endl;
// 	    exit(EXIT_FAILURE);
// 	  }
// 	  if (xprime==k+x) // copy from x as global
// 	    lp[k+x] = alpha[s][k+x]*
// 	      (p_global_norec+(1.-p_global_norec)*(1.-plocalcopy)/double(globalcount[s+1]));
// 	  else 
// 	    lp[k+x] = alpha[s][k+x]*
// 	      (1.-p_global_norec)*(1.-plocalcopy)/double(globalcount[s+1]);
// 	  if (popcount[location]>0)
// 	    lp[x] = alpha[s][x]*
// 	      (1.-p_global_norec)*plocalcopy/double(popcount[location]);
// 	  else lp[x]=0.0;
// 	}
//       }
//       path[k][s]=gen_from_pr(lp,2*k,r,psample);
//       lpsample += psample;
//     }
//     for (int ii=0;ii<sites;ii++) {
//       if (path[k][ii]>=k) {
// 	global[k][ii]=true;
// 	globalcount[ii]+=1;
//       } else global[k][ii]=false;
//     }
//     popcount[location]+=1;
//   }
//   return lprob;
// }

/** calculate the likelihood for one rearrangement of the haplotypes conditional
 * on the globalness - should be no dependency on F 
 */
// template <typename T>
//   double spact<T>::likelihood(double theta, TNT::Array2D<bool> &global)
// {
//   // alpha[x,s] defined as the probability of sampling the first s sites of
//   // your haplotype, and of sampling haplotype x.
//   TNT::Array2D<double> alpha(sites,n-1);
//   std::vector<int> popcount(npops_,0);     // the count in each of the populations
//   TNT::Array1D<int> globalcount(sites,1);  // the count in the "global" pool
  
//   popcount[pop[order[0]]] +=1;             // add one to the local count
  
//   double lprob=0.0;                        // the first sample has probability 1
//   // now consider each of the subsequent samples in turn.
//   // Because my arrays are 0 offset haplotype k is the k+1 st sample
//   // so use the same k as in Li and Stephens' appendix
//   for (int k=1;k<n;k++) {
//     int location = pop[order[k]];  // the location of the haplotype under consideration
//     // what is the probability that a recombination is copied from a "local" chromosome?
//     double plocalcopy=
//       F[location]*double(popcount[location])/(1.+(double)(popcount[location]-1.)*F[location]);
//     // what are the match probabilities for the global population
//     double gammamatchglobal=double(globalcount[0])/(double(globalcount[0])+theta)
//       +0.5*theta/(double(globalcount[0])+theta);
//     double gammamatchlocal = 1.0;       
//     // we can calculate the probabilities for the first site directly
//     // there are no transitions and we assume exchangeability of haplotypes 
//     // loop over the haplotypes from which the site could be copied
//     if (global[k][0]) {
//       for (int x=0;x<k;x++) {
// 	if (global[x]) {
// 	  if (haps[hapsder[x]][0]==haps[order[k]][0])   // matches 
// 	    alpha[0][x] = gammamatchglobal/(double(globalcount[0]));
// 	  else  
// 	    alpha[0][x] = (1.-gammamatchglobal)/(double(globalcount[0]));
// 	} else alpha[0][x]=0.0;
//       }
//     } else {  // local copy
//       for (int x=0;x<k;x++) {
// 	if (pop[order[x]]==location) {
// 	  if (haps[order[x]][0]==haps[order[k]][0]) 
// 	    alpha[0][x]=1.0;
// 	  else alpha[0][x]=0.0;
// 	} else alpha[0][x]=0.0;
//       }
//     }
//     double logscalefactor=rescale(alpha[0],2*k); //MacDonald & Zucchini fix
//     // calculate alpha for the rest of the sites
//     for (int s=1;s<sites;s++) {
//       // first the probability of recombination. 
//       // These values can be justified by reference to the ESD
//       double p_local_norec=exp(-rho*d[s-1]/double(globalcount[s-1]));
//       double p_global_norec=exp(-rho*d[s-1]/double(popcount[location]));        
//       gammamatchglobal=double(globalcount[s])/(double(globalcount[s])+theta)
// 	+0.5*theta/(double(globalcount[s])+theta);
//       // now calculate the second term from Li & Stephens
//       // not quite as simple as we have to work out every combination 
//       // copy from local and global 
//       double smfromlocal=0.0,smfromglobal=0.0;
//       if (global[k][s]) {         // copying to a global 
// 	for (int ii=0;ii<k;ii++) {
// 	  if (global[x][s-1]) {     // both global
// 	    smfromglobal += alpha[s-1][ii]*(1.- p_global_norec);
// 	  } else {                       // target global copied from local
// 	    if (location==pop[order[x]]) // must be the same location
// 	      smfromlocal += alpha[s-1][ii]*(1.- p_local_norec);
// 	  }
// 	}
// 	// now calculate the alpha values from each of the difference chromosomes
// 	for (int x=0;x<k;x++) {
// 	  if (global[x][s]) {    
// 	    alpha[s][x] = p_global_norec*alpha[s-1][x] 
// 	      + (1.-p_global_norec)*smfromglobal*(1.-plocalcopy)/double(globalcount[s]);;

// 	    if (haps[order[x]][0]==haps[order[k]][0])   // matches 
// 	      alpha[s][x]*=gammamatchglobal;
// 	    else alpha[s][x]*=(1.-gammamatchglobal);
// 	  } else {
// 	    alpha[s][x] = p_local_norec*alpha[s-1][x] 
// 	      + (1.-p_local_norec)*smfromlocal*(1.-plocalcopy)/double(globalcount[s]);;
// 	  }
// 	    else 
// 	      alpha[s][x]= ;
// 	    if (pop[order[x]]==location)
// 	    alpha[s][x]=p_local_norec*alpha[s-1][x]+sm;
// 	  else alpha[s][x]=0.0;
	  
// 	  if (global[x][s])
// 	  alpha[s][k+x]=p_global_norec*alpha[s-1][k+x]+smlocalglobal+smglobalglobal; 
// 	else alpha[s][k+x]=0.0;

// 	if (haps[order[x]][s]==h[haps[order[k]][s]) {
// 	  alpha[s][x]*=gammamatchlocal;
// 	  alpha[s][x+k] *= gammamatchglobal;
// 	} else {
// 	  alpha[s][x]*=gammamisslocal;
// 	  alpha[s][x+k] *= gammamissglobal;
// 	}
//       } else {
// 	for (int ii=0;ii<k;ii++) {
// 	  if (global[x][s]) {  // target local copied from global
// 	    smfromglobal += alpha[s-1][ii]*(1.- p_global_norec);
// 	  } else {                    // both local
// 	    smfromlocal += alpha[s-1][ii]*(1.- p_local_norec);
// 	  }
// 	}
//       }
   
//       for (int x=0;x<k;x++) {
// 	if (pop[order[x]]==location)
// 	  alpha[s][x]=p_local_norec*alpha[s-1][x]+sm;
// 	else alpha[s][x]=0.0;
	
// 	if (global[x][s])
// 	  alpha[s][k+x]=p_global_norec*alpha[s-1][k+x]+smlocalglobal+smglobalglobal; 
// 	else alpha[s][k+x]=0.0;

// 	if (haps[order[x]][s]==h[haps[order[k]][s]) {
// 	  alpha[s][x]*=gammamatchlocal;
// 	  alpha[s][x+k] *= gammamatchglobal;
// 	} else {
// 	  alpha[s][x]*=gammamisslocal;
// 	  alpha[s][x+k] *= gammamissglobal;
// 	}
//       }
//       logscalefactor+=rescale(alpha[s],2*k); // The scale factor of MacDonald and Zucchini
//     }

//     // calculate the prob for haplotype k by summing over all the 
//     // paths to the last site
//     double prob=std::accumulate(alpha[sites-1],alpha[sites-1]+2*k,0.0);
//     lprob+=log(prob)+logscalefactor;
 
//     if (log(prob)+logscalefactor>0.0) {
//       std::cerr << "problem, added a probability > 1 " << lprob 
// 		<< " log(prob) = " << log(prob) << " logscalefactor = " << logscalefactor 
// 		<< " added " << log(prob)+logscalefactor << std::endl;
//       exit(EXIT_FAILURE);
//     }
    
//     return lprob;
// }
// }}}

