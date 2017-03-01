#include "pac.H"
#include "utilityfunctionals.H"
#include <cassert>
#include <vector>
#include <iterator>
// utility function to rescale probabilities so that the total probability is 1
// returns the logarithm of the rescaling factor for rescaling
double rescale(double *a, int count) 
{   
  double scalefac=std::accumulate(a,a+count,0.0);
#ifndef NDEBUG  
   if (scalefac>1.0001) {
     std::cout << "warning -- sum of probabilities = " << scalefac << "  > 1 +++++++++++++" 
	       << std::endl;
     for (int i=0;i<count;i++) std::cout << a[i] << " ";
     std::cout << std::endl;
   }
#endif
   for (int i=0;i<count;i++) a[i]/=scalefac;
   return log(scalefac);
}
/******************************************************************************/
/*  calculate the likelihood for a single arrangement of the samples          */
/******************************************************************************/
double pac::calc(double rho)
{
  // alpha[x,s] defined to be the probability of sampling the first s sites of
  // your haplotype, and of sampling haplotype x 
  TNT::Array2D<double> alpha(sites,n-1);
  std::vector<double> lprob(n,0.0);
  // the first sample has probability 1
  lprob[0]=0.0;
  // now consider each of the subsequent samples in turn.
  // Because my arrays are 0 offset haplotype k is the k+1 st sample
  // so use the same k as in Ni and Stephens' appendix 
  for (int k=1;k<n;k++) {
    // some quantities we only need to calculate once for each k
    double gammamatch=double(k)/(k+thetahat)+0.5*thetahat/(k+thetahat);
    double gammamiss=1.-gammamatch;     // = 0.5*thetahat/(k+thetahat);
    // we can calculate the probabilities for the first site directly 
    // loop over the haplotypes from which the site could be copied
    for (int x=0;x<k;x++) {
      if (h[x][0]==h[k][0]) alpha[0][x]=gammamatch/double(k);
      else alpha[0][x]=gammamiss/double(k);
    }
    // fix as in MacDonald and Zucchini so that the mean alpha is 1
    double logscalefactor=rescale(alpha[0],k);    
    // calculate alpha for the rest of the sites
    for (int s=1;s<sites;s++) {
      // the "second term" from Li and Stephens 
      double sm=std::accumulate(alpha[s-1],alpha[s-1]+k,0.0);
      sm /= (double)(k);       
      double p=exp(-rho*d[s-1]/(double(k)));
     
      for (int x=0;x<k;x++) {
        alpha[s][x]=p*alpha[s-1][x]+(1-p)*sm;
	if (h[x][s]==h[k][s]) alpha[s][x]*=gammamatch;
	else alpha[s][x]*=gammamiss;
      }
      // the scale factor of MacDonald and Zucchini
      logscalefactor+=rescale(alpha[s],k);
    }
    // calculate the prob for haplotype k by summing over all the 
    // paths to the last site
    double prob=std::accumulate(alpha[sites-1],alpha[sites-1]+k,0.0);
    lprob[k]=log(prob)+logscalefactor;
  }
  return std::accumulate(lprob.begin(),lprob.end(),0.0);
}
/****************************************************************************/
/** calculate the likelihood and sample the paths of those that are copied  */
/****************************************************************************/
double pac::SamplePath(TNT::Array2D<int> &samp,const double &rho, rng &r)
{
  // alpha[x,s] defined to be the probability of sampling the first s sites of
  // your haplotype, and of sampling haplotype x 
  TNT::Array2D<double> alpha(sites,n-1);
  // the first sample has probability 1
  double lprob=0.0;
  // now consider each of the subsequent samples in turn.
  // Because my arrays are 0 offset haplotype k is the k+1 st sample
  // so use the same k as in Li and Stephens' appendix 
  for (int k=1;k<n;k++) {
     // some quantities we only need to calculate once for each k
    double gammamatch=double(k)/(k+thetahat)+0.5*thetahat/(k+thetahat);
    double gammamiss=0.5*thetahat/(k+thetahat);
    // we can calculate the probabilities for the first site directly 
    // loop over the haplotypes from which the site could be copied
    for (int x=0;x<k;x++) {
      if (h[x][0]==h[k][0]) alpha[0][x]=gammamatch/double(k);
      else alpha[0][x]=gammamiss/double(k);
    }
    // fix as in MacDonald and Zucchini so that the mean alpha is 1
    double logscalefactor=rescale(alpha[0],k);    
    // calculate alpha for the rest of the sites
    for (int s=1;s<sites;s++) {
      // the "second term" from Li and Stephens 
      double sm=std::accumulate(alpha[s-1],alpha[s-1]+k,0.0);
      sm /= (double)(k);       
      double p=exp(-rho*d[s-1]/(double(k)));
      for (int x=0;x<k;x++) {
	alpha[s][x]=p*alpha[s-1][x]+(1-p)*sm;
	if (h[x][s]==h[k][s]) alpha[s][x]*=gammamatch;
	else alpha[s][x]*=gammamiss;
      }
      // the scale factor of MacDonald and Zucchini
      logscalefactor+=rescale(alpha[s],k);
    }
    // calculate the prob for haplotype k by summing over all the 
    // paths to the last site
    double prob=std::accumulate(alpha[sites-1],alpha[sites-1]+k,0.0);
    lprob+=log(prob)+logscalefactor;
    // and now for the "backwards" part of the algorithm, 
    // as in Boys et al 2000 - the last site is "site-1"
    std::vector<double> lp(k);
    samp[k][sites-1] = gen_from_p(alpha[sites-1],alpha[sites-1]+k,r);
    for (int s=sites-2;s>=0;s--) {
	for (int x=0;x<k;x++) {
	  double p_stay = exp(-rho*d[s]/k);
	  lp[x]=alpha[s][x];
	  // is it the same as the "next" site
	  if (samp[k][s+1]==x) lp[x] *= (p_stay+(1.-p_stay)/k); 
	  else lp[x] *= (1.-p_stay)/k; 
	}
	samp[k][s]=gen_from_pr(lp,k,r);
    }
  }
  return lprob;
}
/******************************************************************************/
/** get an estimate of the numbers and positions of recombinations            */
/******************************************************************************/
TNT::Array1D<double> pac::recombinations(const double &rho, 
    int seed, bool scale) 
{
#ifndef USE_R
    rng r(seed);
#else
    rng r;
#endif  
    h=st;
    TNT::Array1D<double> recombs(sites-1,0.0);
    TNT::Array2D<int> recs(nperm,sites-1,0);
    TNT::Array2D<int> a(n,sites);
    std::vector<double> ll(nperm);
    
    for (int i=0;i<nperm;i++) {
      ll[i]=SamplePath(a,rho,r);
      for (int j=0;j<n;j++)
	for (int k=0;k<sites-1;k++)
	  recs[i][k]+=int(a[j][k]!=a[j][k+1]);
    }
    if (!scale ) for (int i=0;i<nperm;i++) ll[i]=0.0;	  
    double sc=*(std::max_element(ll.begin(),ll.end()));
    for (int i=0;i<nperm;i++) ll[i] = exp(ll[i]-sc);  
    double divby=std::accumulate(ll.begin(),ll.end(),0.0);
    
    for (int i=0;i<nperm;i++) { 
      for (int s=0;s<sites-1;s++) 
	recombs[s]+=ll[i]*recs[i][s];
    }
    for (int s=0;s<sites-1;s++) recombs[s]/=divby;
    return recombs;
}
/******************************************************************************/
/** Get samples of an unseen SNP between the pos and pos+1 SNP with 
 *   relative distances d1                                                    */
/******************************************************************************/
TNT::Array3D<int> pac::SampleSNP(double &lprob,int nsamp, TNT::Array1D<int> &pos
  , TNT::Array1D<double> &d1, const double &rho, rng &r, int minfreq)
{
  lprob=0.0;
  assert(pos.dim()==d1.dim()); 
  int nSNP=pos.dim();
  TNT::Array3D<int> path(nSNP,nsamp,n);
  TNT::Array1D<double> gammamatch(n);
  // alpha[x,s] defined to be the probability of sampling the first s sites of
  // your haplotype, and of sampling haplotype x 
  TNT::Array2D<double> alpha(sites,n-1);
  // the first sample has probability 1
  // double lprob=0.0;
  for (int SNP=0;SNP<nSNP;SNP++) 
    for (int jj=0;jj<nsamp;jj++) 
      path[SNP][jj][0]=0;
  // now consider each of the subsequent samples in turn.
  // Because my arrays are 0 offset haplotype k is the k+1 st sample
  // so use the same k as in Li and Stephens' appendix 
  for (int k=1;k<n;k++) {
     // some quantities we only need to calculate once for each k
    gammamatch[k]=double(k)/(k+thetahat)+0.5*thetahat/(k+thetahat);
    // we can calculate the probabilities for the first site directly 
    // loop over the haplotypes from which the site could be copied
    for (int x=0;x<k;x++) {
      if (h[x][0]==h[k][0]) alpha[0][x]=gammamatch[k]/double(k);
      else alpha[0][x]=(1.-gammamatch[k])/double(k);
    }
    // fix as in MacDonald and Zucchini so that the mean alpha is 1
    double logscalefactor=rescale(alpha[0],k);    
    // calculate alpha for the rest of the sites
    for (int s=1;s<sites;s++) {
      // the "second term" from Li and Stephens 
      double sm=std::accumulate(alpha[s-1],alpha[s-1]+k,0.0);
      sm /= (double)(k);       
      double p=exp(-rho*d[s-1]/(double(k)));
      for (int x=0;x<k;x++) {
	alpha[s][x]=p*alpha[s-1][x]+(1-p)*sm;
	if (h[x][s]==h[k][s]) alpha[s][x]*=gammamatch[k];
	else alpha[s][x]*=(1.-gammamatch[k]);
      }
      // the scale factor of MacDonald and Zucchini
      logscalefactor+=rescale(alpha[s],k);
    }
    // calculate the prob for haplotype k by summing over all the 
    // paths to the last site
    double prob=std::accumulate(alpha[sites-1],alpha[sites-1]+k,0.0);
    lprob+=log(prob)+logscalefactor;
    //    std::cerr << lprob << std::endl;
    for (int SNP=0;SNP<nSNP;SNP++) {
      int position=pos[SNP];// this SNP lies between columns position and position+1
      double d2=d[position]-d1[SNP]; 
      assert(d1[SNP]<d[position]); 
      assert(fabs(d1[SNP]+d2-d[position])<1E-8);
      TNT::Array1D<double> newalpha(k);
      // start by recalculating for the different values
      // the "second term" from Li and Stephens
      double sm=std::accumulate(alpha[position]
                                ,alpha[position]+k,0.0)/(double(k));
      double p=exp(-rho*d1[SNP]/(double(k)));
      for (int x=0;x<k;x++) 
        newalpha[x] = p*alpha[position][x]+(1-p)*sm;
      
      // and now for the "backwards" part of the algorithm   
      std::vector<double> lp(k);
      for (int jj=0;jj<nsamp;jj++) {
	// sample a path "above"
	int pathabove= gen_from_p(alpha[position+1],alpha[position+1]+k,r);
	for (int x=0;x<k;x++) {
	  double p_stay = exp(-rho*d2/k);
	  lp[x]=newalpha[x];
	  if (pathabove==x) lp[x] *= (p_stay+(1.-p_stay)/k); 
	  else lp[x] *= (1.-p_stay)/k; 
	}
	path[SNP][jj][k] = gen_from_pr(lp,k,r);
      }
    }
  }
  TNT::Array3D<int> sampSNP(nSNP,nsamp,n);
  // alpha[x,s] defined to be 
  for (int SNP=0;SNP<nSNP;SNP++) {
    for (int jj=0;jj<nsamp;) {
      bool keep=false;
      if (r()<0.5) sampSNP[SNP][jj][0]=1;
      else sampSNP[SNP][jj][0]=0;
      for (int k=1;k<n;k++) {
	if (r()<gammamatch[k]) sampSNP[SNP][jj][k]=sampSNP[SNP][jj][path[SNP][jj][k]];
	else {
	  sampSNP[SNP][jj][k]=1-sampSNP[SNP][jj][path[SNP][jj][k]]; // a mutation 
	  keep=true;  // we can try to keep this set!
	}
      }
      if (keep) {
	int freq=0;
	for (int i=0;i<n;i++) freq+=sampSNP[SNP][jj][i]==1;
	if (freq>=minfreq&&n-freq> minfreq) jj++;  // keep this and get the next
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
TNT::Array3D<int> pac::SampleCaseControlSNP(TNT::Array1D<double> &weightings
					    ,int nsamp, TNT::Array1D<int> &pos
 ,TNT::Array1D<double> &d1, double rho, TNT::Array1D<int> &CaseControl,int seed)
  {
#ifdef USE_R
    rng r;
#else
    rng r(seed);
#endif
    int nSNP=pos.dim();
    TNT::Array3D<int> tables(pos.dim(),nsamp*nperm,6,0);
    
    for (int i=0;i<nperm;i++) {
      h=st; // reset the permutations 
     
      std::vector<int> rev(n);
      {
	std::vector<int> prms=r.integer_permutations(n,n);
	for (int j=0;j<n;j++) {
	  h[j]=st[prms[j]];
	  rev[prms[j]]=j;
	}
      }
     
      TNT::Array3D<int> res =  SampleSNP(weightings[i],nsamp,pos,d1,rho,r,n/5);

      for (int SNP=0;SNP<nSNP;SNP++) {
	for (int j=0;j<nsamp;j++) {
	  for (int k=0;k<n;k+=2) {
	    if (res[SNP][j][rev[k]]==0&&res[SNP][j][rev[k+1]]==0) 
	      tables[SNP][i*nsamp+j][CaseControl[k]*3]+=1;
	    
	    else if (res[SNP][j][rev[k]]==1&&res[SNP][j][rev[k+1]]==1) 
	      tables[SNP][i*nsamp+j][2+CaseControl[k]*3]+=1;
	    
	    else  tables[SNP][i*nsamp+j][1+CaseControl[k]*3]+=1;
	  }
	}
      }
    }    
    double scale=*(std::max_element(&weightings[0],&weightings[0]+nperm));
    for (int i=0;i<nperm;i++) 
      weightings[i] = exp(weightings[i]-scale);
    
    return tables;
  }
/****************************************************************************
 *
 *  utility function to get data in the correct format for SampleSNP
 *
 ******************************************************************************/
void getequalspacedpositionsanddistances(const TNT::Array1D<double> &distances
					 , TNT::Array1D<int> &pos,TNT::Array1D<double> &gap
					 , TNT::Array1D<double> &x)
{
  int n=pos.dim();
  double pmin=distances[0];
  double pmax=distances[distances.dim()-1]; 
  double trygap=(pmax-pmin)/(n+1.);
  int index=0;
  for (int i=0;i<n;i++) {
    double position=pmin+(i+1.)*trygap;
    for (;;) {
      assert(distances[index]<position);
      if (distances[index+1] >= position) break;
      index++;
      assert(index<distances.dim()-1);
    }
    if (fabs(distances[index]-position)<1E-8) {
      position -= 1E-6; // so not equal
    }
    pos[i]=index;
    gap[i]=position-distances[index];
    assert(gap[i]>0.0);
    x[i]=position;
  }
}
/****************************************************************************/
/** calculate the likelihood and sample the paths of those that are copied  
    for data with missing values   -- not written yet                       */
/****************************************************************************/
double pacmissing::SamplePath(TNT::Array2D<int> &samp,const double &rho, rng &r)
{
  // alpha[x,s] defined to be the probability of sampling the first s sites of
  // your haplotype, and of sampling haplotype x 
  TNT::Array2D<double> alpha(sites,n-1);
  // the first sample has probability 1
  double lprob=0.0;
 // now guess at the first missing values
  for (int s=0;s<sites;s++) {
    if (h[0][s]<0) { // missing
      if (r()<0.5) h[0][s]=1;
      else h[0][s]=0;
    }
  }
  // now consider each of the subsequent samples in turn.
  // Because my arrays are 0 offset haplotype k is the k+1 st sample
  // so use the same k as in Li and Stephens' appendix 
 
  for (int k=1;k<n;k++) {
     // some quantities we only need to calculate once for each k
    double gammamatch=double(k)/(k+thetahat)+0.5*thetahat/(k+thetahat);
    double gammamiss=0.5*thetahat/(k+thetahat);
    // we can calculate the probabilities for the first site directly 
    // loop over the haplotypes from which the site could be copied
    if (h[k][0]<0) { //this site missing
      for (int x=0;x<k;x++)  alpha[0][x]=1./double(k);
    } else {
      for (int x=0;x<k;x++) {
	if (h[x][0]==h[k][0]) alpha[0][x]=gammamatch/double(k);
	else alpha[0][x]=gammamiss/double(k);
      }
    }
    // fix as in MacDonald and Zucchini so that the mean alpha is 1
    double logscalefactor=rescale(alpha[0],k);    
    // calculate alpha for the rest of the sites
    for (int s=1;s<sites;s++) {
      // the "second term" from Li and Stephens 
      double sm=std::accumulate(alpha[s-1],alpha[s-1]+k,0.0);
      sm /= (double)(k);       
      double p=exp(-rho*d[s-1]/(double(k)));
      for (int x=0;x<k;x++) {
	alpha[s][x]=p*alpha[s-1][x]+(1-p)*sm;
	if (h[k][s]>0) { // not missing 
	  if (h[x][s]==h[k][s]) alpha[s][x]*=gammamatch;
	  else alpha[s][x]*=gammamiss;
	}
      }
      // the scale factor of MacDonald and Zucchini
      logscalefactor+=rescale(alpha[s],k);
    }
    // calculate the prob for haplotype k by summing over all the 
    // paths to the last site
    double prob=std::accumulate(alpha[sites-1],alpha[sites-1]+k,0.0);
    lprob+=log(prob)+logscalefactor;
    // and now for the "backwards" part of the algorithm, 
    // as in Boys et al 2000 - the last site is "site-1"
    std::vector<double> lp(k);
    samp[k][sites-1] = gen_from_p(alpha[sites-1],alpha[sites-1]+k,r);
    for (int s=sites-2;s>=0;s--) {
	for (int x=0;x<k;x++) {
	  double p_stay = exp(-rho*d[s]/k);
	  lp[x]=alpha[s][x];
	  // is it the same as the "next" site
	  if (samp[k][s+1]==x) lp[x] *= (p_stay+(1.-p_stay)/k); 
	  else lp[x] *= (1.-p_stay)/k; 
	}
	samp[k][s]=gen_from_pr(lp,k,r);
	if (h[k][s]<0) {
	  if (r()<gammamatch)
	    h[k][s]=h[samp[k][s]][s];
	  else 
	    h[k][s]=1-h[samp[k][s]][s];
      }
    }
  }
  return lprob;
}


/******************************************************************************/
/** Get samples of an unseen SNP between the pos and pos+1 SNP with 
 *   relative distances d1 
 */
/******************************************************************************/
TNT::Array3D<int> pacmissing::SampleSNP(double &lprob,int nsamp, TNT::Array1D<int> &pos
  , TNT::Array1D<double> &d1, const double &rho, rng &r, int minfreq)
{
  lprob=0.0;
  assert(pos.dim()==d1.dim()); 
  int nSNP=pos.dim();
  TNT::Array3D<int> path(nSNP,nsamp,n);
  TNT::Array1D<double> gammamatch(n);
  // alpha[x,s] defined to be the probability of sampling the first s sites of
  // your haplotype, and of sampling haplotype x 
  TNT::Array2D<double> alpha(sites,n-1);
  // the first sample has probability 1
  // double lprob=0.0;
  for (int SNP=0;SNP<nSNP;SNP++) 
    for (int jj=0;jj<nsamp;jj++) 
      path[SNP][jj][0]=0;
  // now consider each of the subsequent samples in turn.
  // Because my arrays are 0 offset haplotype k is the k+1 st sample
  // so use the same k as in Li and Stephens' appendix 
  for (int k=1;k<n;k++) {
     // some quantities we only need to calculate once for each k
    gammamatch[k]=double(k)/(k+thetahat)+0.5*thetahat/(k+thetahat);
    // we can calculate the probabilities for the first site directly 
    // loop over the haplotypes from which the site could be copied
    for (int x=0;x<k;x++) {
      if (h[x][0]==h[k][0]) alpha[0][x]=gammamatch[k]/double(k);
      else alpha[0][x]=(1.-gammamatch[k])/double(k);
    }
    // fix as in MacDonald and Zucchini so that the mean alpha is 1
    double logscalefactor=rescale(alpha[0],k);    
    // calculate alpha for the rest of the sites
    for (int s=1;s<sites;s++) {
      // the "second term" from Li and Stephens 
      double sm=std::accumulate(alpha[s-1],alpha[s-1]+k,0.0);
      sm /= (double)(k);       
      double p=exp(-rho*d[s-1]/(double(k)));
      for (int x=0;x<k;x++) {
	alpha[s][x]=p*alpha[s-1][x]+(1-p)*sm;
	if (h[x][s]==h[k][s]) alpha[s][x]*=gammamatch[k];
	else alpha[s][x]*=(1.-gammamatch[k]);
      }
      // the scale factor of MacDonald and Zucchini
      logscalefactor+=rescale(alpha[s],k);
    }
    // calculate the prob for haplotype k by summing over all the 
    // paths to the last site
    double prob=std::accumulate(alpha[sites-1],alpha[sites-1]+k,0.0);
    lprob+=log(prob)+logscalefactor;
    //    std::cerr << lprob << std::endl;
    for (int SNP=0;SNP<nSNP;SNP++) {
      int position=pos[SNP];// this SNP lies between columns position and position+1
      double d2=d[position]-d1[SNP]; 
      assert(d1[SNP]<d[position]); 
      assert(fabs(d1[SNP]+d2-d[position])<1E-8);
      TNT::Array1D<double> newalpha(k);
      // start by recalculating for the different values
      // the "second term" from Li and Stephens
      double sm=std::accumulate(alpha[position]
                                ,alpha[position]+k,0.0)/(double(k));
      double p=exp(-rho*d1[SNP]/(double(k)));
      for (int x=0;x<k;x++) 
        newalpha[x] = p*alpha[position][x]+(1-p)*sm;
      
      // and now for the "backwards" part of the algorithm   
      std::vector<double> lp(k);
      for (int jj=0;jj<nsamp;jj++) {
	// sample a path "above"
	int pathabove= gen_from_p(alpha[position+1],alpha[position+1]+k,r);
	for (int x=0;x<k;x++) {
	  double p_stay = exp(-rho*d2/k);
	  lp[x]=newalpha[x];
	  if (pathabove==x) lp[x] *= (p_stay+(1.-p_stay)/k); 
	  else lp[x] *= (1.-p_stay)/k; 
	}
	path[SNP][jj][k] = gen_from_pr(lp,k,r);
      }
    }
  }
  TNT::Array3D<int> sampSNP(nSNP,nsamp,n);
  // alpha[x,s] defined to be 
  for (int SNP=0;SNP<nSNP;SNP++) {
    for (int jj=0;jj<nsamp;) {
      bool keep=false;
      if (r()<0.5) sampSNP[SNP][jj][0]=1;
      else sampSNP[SNP][jj][0]=0;
      for (int k=1;k<n;k++) {
	if (r()<gammamatch[k]) sampSNP[SNP][jj][k]=sampSNP[SNP][jj][path[SNP][jj][k]];
	else {
	  sampSNP[SNP][jj][k]=1-sampSNP[SNP][jj][path[SNP][jj][k]]; // a mutation 
	  keep=true;  // we can try to keep this set!
	}
      }
      if (keep) {
	int freq=0;
	for (int i=0;i<n;i++) freq+=sampSNP[SNP][jj][i]==1;
	if (freq>=minfreq&&n-freq> minfreq) jj++;  // keep this and get the next
      }
    }
  }
  return sampSNP;
}
