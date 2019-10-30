#include "SeqDP.H"
#include "FlowgramCall.H"
#include "progressBar.H"

double llDirichletProcess(double a, const std::vector<int> &y,int tot) 
{
  assert(std::accumulate(y.begin(),y.end(),0)==tot);
  size_t k=y.size();
  double lp = k*log(a) + lgamma(a)  -lgamma(static_cast<double>(tot+a));
  for (size_t i=0;i<k;i++) lp +=lgamma(static_cast<double>(y[i]));
  return lp;
}


/** Take a sequence away and add it at random                  */
void SeqDP::gibbs(int index) {
  //  std::cerr << "seq " << index << " alpha = " << alpha() << " have " << nalls() << " alleles: ";

  // std::cerr << "now have " << nalls() << " alleles\n";
  int prev=lab_[index];
  if (copies_[prev]==1) {
    //    std::cerr << "removing a singleton\n";
    for (size_t jj=0;jj<lab_.size();jj++) if (lab_[jj]==nalls()-1) lab_[jj]=prev;
    copies_[prev] = copies_.back(); 
    copies_.pop_back();
    seq_[prev] = seq_.back();
    seq_.pop_back(); 
    //  std::cerr << "removed it\n";
    //assert(nalls()==na-1);  
  } else {
    copies_[prev]--;
  }
  
  std::vector<double> lp(nalls()+1);
  for (int ii=0;ii<nalls();ii++) {
    lp[ii] = seq_[ii].lprob(fg[index],zeropars,sigmapars,InitialTag);
  }
  
  flowSeq fs = randomCall(fg[index],zeropars,sigmapars,r,maxCopies,InitialTag);

  lp[nalls()] = fs.lprob(fg[index],zeropars,sigmapars,InitialTag);
  //   printvector(lp);
  double lpmax = *(std::max_element(lp.begin(),lp.end()));
  for (int ii=0;ii<nalls();ii++) {
    lp[ii] = exp(lp[ii]-lpmax)*copies_[ii]/(alpha()+n()-1);
  } 
  lp[nalls()] = exp(lp[nalls()]-lpmax)*alpha()/(alpha()+n()-1);
  //   printvector(lp);
  int wh = gen_from_p(lp,r);
  if (wh==nalls()) { // make sure it is new
    std::vector<flowSeq>::iterator ff=std::find(seq_.begin(),seq_.end(),fs);
    if (ff != seq_.end()) {
      size_t wh = std::distance(seq_.begin(),ff);
      lab_[index] = wh;
      copies_[wh]+=1;
      return;
    }
    //    std::cerr << "new sequence \n";
    copies_.push_back(1);
    seq_.push_back(fs);
  } else {
    //     std::cerr << "add to sequence " << wh << " with " << copies_[wh]<< " copies\n";
    copies_[wh]++;
  }
  lab_[index] = wh;
   
}


void SeqDP::sequentialStart(bool initParameters) {
  progressBar p(std::cerr,"sequential start",n());
  if (initParameters) initialiseParameters();
  seq_.push_back(randomCall(fg[0],zeropars,sigmapars,r,maxCopies,InitialTag));
  copies_.push_back(1);
  lab_.push_back(0);
  for (size_t jj=1;jj<n();jj++) {
    p.update(jj);
     std::vector<double> lp(nalls()+1);
     for (int ii=0;ii<nalls();ii++) {
       lp[ii] = seq_[ii].lprob(fg[jj],zeropars,sigmapars,InitialTag);
     }
     // assert(maxCopies>0);
     flowSeq fs = randomCall(fg[jj],zeropars,sigmapars,r,maxCopies,InitialTag);
     lp[nalls()] = fs.lprob(fg[jj],zeropars,sigmapars,InitialTag);
     double lpmax = *(std::max_element(lp.begin(),lp.end()));
     for (int ii=0;ii<nalls();ii++) {
       lp[ii] = exp(lp[ii]-lpmax)*copies_[ii]/(alpha()+n());
     } 
     lp[nalls()] = exp(lp[nalls()]-lpmax)*alpha()/(alpha()+n());
     
     int wh = gen_from_p(lp,r);
     if (wh==nalls()) {
       copies_.push_back(1);
       seq_.push_back(fs);
     } else {
       copies_[wh]++;
     }
     lab_.push_back(wh);
  }
  p.finish();
}


void SeqDP::initialiseParameters()
{
  if (not initialised) {
    initialised=true;
    alpha_=alphaprior->sample(r);
    SampleParameters();
  }
}


/** Call the genotypes using new parameters parameters     
 * call all samples independently but clustering samples that are the same               
 */
void SeqDP::initialise()
{
  initialiseParameters();
  //  printvector(zeropars);
  // printvector(sigmapars);
  for (size_t ii=0;ii<n();ii++) {
    // call the flowsequence
    flowSeq current = randomCall(fg[ii],zeropars,sigmapars,r,maxCopies,InitialTag);
    std::vector<flowSeq>::iterator ff=std::find(seq_.begin(),seq_.end(),current);
    if (ff==seq_.end()) {                 // this flowseq not seen before
      lab_.push_back(seq_.size());       // the next label
      seq_.push_back(current);
      copies_.push_back(1);
    } else  {
      lab_.push_back(std::distance(seq_.begin(),ff));
      copies_[lab_[ii]]+=1;
    }		
  }
}
/** gather all posterior calculations in one place to make it all happen 
 *  quickly and easily   */
void SeqDP::CalcPostPar()
{
  // we are assume that zeroprior is a gamma distribution that is conjugate 
  // to the lognormal distribution for intensities when we have a flowSeq with count 0
  int maxd = tauprior.size();

  std::vector<double> xx(3,0.0);
  std::vector<double> yy(maxd,0.0);
  std::vector<int> nn(maxd,0);

  assert(fg[0].size()==seq_[0]().size()+InitialTag.size());  
  size_t tagsize=InitialTag.size();
  for (size_t ii=0;ii<n();ii++) {
    size_t len = fg[ii].size();
    assert(len>tagsize);
    for (size_t jj=0;jj<tagsize;jj++) {
      if (InitialTag[jj]==0) {
	double lx=log(fg[ii][jj]+0.001);
	xx[0] += 1.0;                              // count n
	xx[1] += lx;                               // sum_x
	xx[2] += lx*lx;
      } else {
	size_t wh=std::min(InitialTag[jj],maxd)-1;
	nn[wh] += 1;
	yy[wh] += pow((fg[ii][jj]-InitialTag[jj]),2.0);
      }
    }
    flowSeq &lfs=seq_[lab_[ii]];
    for (size_t jj=tagsize;jj<len;jj++) {
      if (lfs[jj-tagsize]==0) {
	double lx=log(fg[ii][jj]+0.001);
   	xx[0] += 1.0;                              // count n
   	xx[1] += lx;                               // sum_x
   	xx[2] += lx*lx;
      } else {
   	size_t wh=std::min(lfs[jj-tagsize],maxd)-1;
   	nn[wh] += 1;
   	yy[wh] += pow((fg[ii][jj]-lfs[jj-tagsize]),2.0);
      }
    }
  }
  //     std::vector<int>::const_iterator cc=InitialTag.begin();
  //     std::vector<double>::const_iterator ff=fg[ii].begin();

  //     while (ff!=fg[ii].end()) {
  //       if (*cc==0) {
  // 	double lx=log(*ff+0.001);
  // 	xx[0] += 1.0;                              // count n
  // 	xx[1] += lx;                               // sum_x
  // 	xx[2] += lx*lx;
  //       } else {
  // 	size_t wh=std::min(*cc,maxd)-1;
  // 	nn[wh] += 1;
  // 	yy[wh] += pow((*ff-*cc),2.0);
  //       }
  //       cc++;
  //       if (cc==InitialTag.end()) cc=seq_[lab_[ii]]().begin();
  //       ff++;
  //     }
  //   }
  //     size_t len = fg[ii].size();
  //     for (size_t jj=0;jj<len;jj++) {
  //       if (lfs[jj]==0) {
  // 	double lx=log(fg[ii][jj]+0.001);
  // 	xx[0] += 1.0;                              // count n
  // 	xx[1] += lx;                               // sum_x
  // 	xx[2] += lx*lx;
  //       } else {
  // 	size_t wh=std::min(lfs[jj],maxd)-1;
  // 	nn[wh] += 1;
  // 	yy[wh] += pow((fg[ii][jj]-lfs[jj]),2.0);
  //       }
  //     }

  assert((xx[2]-xx[1]*xx[1]/xx[0])>0.0);
  std::vector<double> pars =zeroprior->parameters();

  zeropost.resize(4); 
  zeropost[0]  = (pars[0]*pars[1]+xx[1])/(pars[1] + xx[0]);
  zeropost[1]  = pars[1]+xx[0];
  zeropost[2] =  pars[2]+xx[0]/2.0;
  zeropost[3] =  pars[3]+(xx[2]-xx[1]*xx[1]/xx[0])/2.0 + (xx[0]*pars[1]*pow((xx[1]/xx[0]-pars[0]),2)/(2.0*(xx[0]+pars[1])));
  
  taupost.resize(0);
  for (int i=0;i<maxd;i++) {
    std::vector<double> pr = tauprior[i]->parameters();
    pr[0] += static_cast<double>(nn[i])/2.0;
    pr[1] += pr[1] + yy[i]/2.0;
    taupost.push_back(pr);	
  }
}
/** Sample parameters using either the priors (if the posteriors have not 
 * been calculated, or the posteriors                        */
void SeqDP::SampleParameters()
{
  zeropars.resize(2);
  sigmapars.resize(tauprior.size());    
  if (zeropost.size()==0) { // sample from priors
    //  std::cerr << "sampling from priors" << std::endl;
    zeropars = zeroprior->sample(r);
    zeropars[1] = 1./sqrt(zeropars[1]);
    for (size_t ii=0;ii<sigmapars.size();ii++) sigmapars[ii] =  1./sqrt(tauprior[ii]->sample(r));
  } else {  // sample from posterior
    //   std::cerr << " sampling from posterior with parameters: ";
    //    printvector(zeropost);
    for (size_t ii=0;ii<sigmapars.size();ii++) 
      sigmapars[ii] =  1./sqrt(r.rgamma(taupost[ii][0],1./taupost[ii][1]));
    
    double ztau = r.rgamma(zeropost[2],1./zeropost[3]);
    //   std::cerr << "ztau = " << ztau << std::endl;
    zeropars[0] =  zeropost[0]+r.normal()*sqrt(1.0/(zeropost[1]*ztau));
    zeropars[1] = 1./sqrt(ztau);
    // zeropars = zeroprior->sample(r);
    //   zeropars[1] = 1./sqrt(zeropars[1]);

  }
}


/** Call the genotypes using new parameters parameters                   
    use the current clustering of genotypes to get the new ones     */
void SeqDP::callGenotypes()
{
  std::vector<std::vector<int> > wh(nalls());
  for (size_t jj=0;jj<n();jj++)  {   // cluster flowGrams with the same sequence
    assert(lab_[jj]<nalls()&&lab_[jj]>=0);
    wh[lab_[jj]].push_back(jj);
  }
  for (int  ii=0;ii<nalls();ii++) {
    assert(wh[ii].size()>0);
    // call the flowsequence
    seq_[ii] = randomCallv(fg,wh[ii],zeropars,sigmapars,r,maxCopies,InitialTag);	
  }	
}

/** Call the genotypes using new parameters parameters                   
    use the current clustering of genotypes to get the new ones     */
void SeqDP::gibbsCallSequence(int which)
{
  assert(which>=0&&which<nalls());
  std::vector<int>  samp;
  for (size_t jj=0;jj<n();jj++)  {   // cluster flowGrams with the same sequence
    if (lab_[jj]==which) samp.push_back(jj);
  }
  assert(samp.size()>0);
  seq_[which] = randomCallv(fg,samp,zeropars,sigmapars,r,maxCopies,InitialTag);	
}


double SeqDP::llikelihood() const {
  double temp=0.0;
  for (size_t ii=0;ii<n();ii++) 
    temp += seq_[lab_[ii]].lprob(fg[ii],zeropars,sigmapars,InitialTag);
  return temp;
}
/** What is the prior probability of the model parameters    */
double SeqDP::lprior() const {
  // need to transform
  std::vector<double> tr=zeropars;
  tr[1] = 1.0/pow(zeropars[1],2.0);
  double tmp=zeroprior->log_pdf(tr);
  for (size_t ii=0;ii< tauprior.size() ;ii++) 
    tmp += tauprior[ii]->log_pdf(1./pow(sigmapars[ii],2.0));

  tmp += alphaprior->log_pdf(alpha_);
  tmp += llDirichletProcess(alpha_,copies_,n());
  return tmp;
}
