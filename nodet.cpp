#include "coalescent.H"
#include "nodet.H"
#include <algorithm>  // for max_element
#include <map>

namespace GenTL {
  template <>                 
  void node<std::vector<char> >::write_leaves(int pos)
  {
    if (isleaf()) { // at a leaf
      data()[pos]=1;
    } else {
      left->write_leaves(pos);
      right->write_leaves(pos);
    } 
    return;
  }

 /** Recursively find the position of the mutation              */
  template <> 
  double node<std::vector<char> >::recursemutate(double mtime, const int pos)
  {
    assert(!isleaf());
    
    mtime -= time-left->time;

    if (mtime<0.0) { // mutation happens here
      left->write_leaves(pos);
      return mtime;
    } else { // no mutation here so try to the left
      if (!left->isleaf()) {
	mtime = left->recursemutate(mtime,pos);
	if (mtime<0.0) return mtime;
      }
    }
    // no mutation so try to the right
    mtime -= time-right->time;
    
    if (mtime<0.0) {
      right->write_leaves(pos);
      return mtime;
    } else {
      if (!right->isleaf()) 
	mtime = right->recursemutate(mtime,pos);
    }
    return mtime;
  }
 /** First mutation function - infinite site     */
  template  <>
  void tree<std::vector<char> >::mutate(int nloci, rng &r)
  {
    // first allocate memory for the node data
    for (size_t i=0;i<sample_.size();i++) 
      sample_[i].data().assign(nloci,char(0));

    double len=root->length();

    for (int i=0;i<nloci;i++) {
      double mt=len*r();
#ifndef NDEBUG
      double mutated=root->recursemutate(mt,i);
      assert(mutated<0.0);
#else 
      root->recursemutate(mt,i);
#endif
    }
  }
  
 /** second mutation function - using infinite sites
  * to simulate binary sites by xor - ing the results     
  * This assumes that all the number of mutations are at least 
 - note that this returns an extra columns that should be ignored*/
  template  <>
  void tree<std::vector<char> >::mutate(const std::vector<int> &nmuts, rng &r)
  {
    // first allocate memory for the node data
    int nloci=nmuts.size();
    for (size_t i=0;i<sample_.size();i++) 
      sample_[i].data().assign(nloci+1,char(0));
    // we may need some extra working space

    double len=root->length();
    for (int i=0;i<nloci;i++) {
      double mt=len*r();
      // make sure it is all set to zero
      // as it may have been used as a scratchpad for the last locus
      for (size_t k=0;k<sample_.size();k++)  sample_[k].data()[i]=0;
      if (root->recursemutate(mt,i)>=0.0) {
        throw std::runtime_error("should return negative number in tree::mutate\n");
      }
      
      for (int j=1;j<nmuts[i];j++) {
	double mt2=len*r();
	for (size_t k=0;k<sample_.size();k++)  sample_[k].data()[i+1]=0;
	root->recursemutate(mt2,i+1);
	for (size_t k=0;k<sample_.size();k++) 
	  sample_[k].data()[i] =  sample_[k].data()[i] xor sample_[k].data()[i+1] ;
      } 
    }
  }
 /** Recursively find the position of the mutation              */
  template <> 
  void node<std::vector<char> >::recursemutate(const std::vector<double> &theta, rng &r)
  {
    assert(!isleaf());
    
    double ltime=time-left->time;
    double rtime=time-right->time;
    for (size_t i=0;i<theta.size();i++) {
      char d=data()[i];
      int muts=r.rpoisson(ltime*theta[i]/2.);
      left->data()[i]=(d+muts)%2;
      muts=r.rpoisson(rtime*theta[i]/2.);
      right->data()[i]=(d+muts)%2;
    }
    if (!left->isleaf()) {
      left->recursemutate(theta,r);
    }
    if (!right->isleaf()) {
      right->recursemutate(theta,r);
    }  
  }
   /** the other recurse mutation function               */
//   template <> 
//   void node<std::vector<char> >::recursemutate(size_t start, const std::vector<double> &theta, rng &r)
//   {
//     assert(!isleaf());
    
//     double ltime=time-left->time;
//     double rtime=time-right->time;
//     size_t end=start+theta.size();
//     for (size_t i=start;i<end;i++) {
//       thet=theta[i-start];
//       char d=data()[i];
//       int muts=r.rpoisson(ltime*thet/2.);
//       left->data()[i]=(d+muts)%2;
//       muts=r.rpoisson(rtime*thet/2.);
//       right->data()[i]=(d+muts)%2;
//     }
//     if (!left->isleaf()) {
//     //       left->recursemutate(start,theta,r);
//     }
//     if (!right->isleaf()) {
//     //      right->recursemutate(start,theta,r);
//     }  
//   }
 /** third mutation function - for high mutation rates*/
  template  <>
  void tree<std::vector<char> >::mutate(const std::vector<double> &theta, rng &r)
  {
    //ijw
      // first allocate memory for the node data
    int nloci=theta.size();
    for (size_t i=0;i<sample_.size();i++) 
      sample_[i].data().assign(nloci,char(0));

    for (size_t i=0;i<ancestors_.size();i++) 
      ancestors_[i].data().assign(nloci,char(0));

    root->recursemutate(theta,r);
  }
 /** third mutation function - for high mutation rates*/
//   template  <>
//   void tree<std::vector<char> >::mutate(size_t start, const std::vector<double> &theta, rng &r)
//   {
//     //ijw assumes that the memory is already assigned
//     root->recursemutate(r,theta,r);
//   }

 /** A utility function to do the simulations for a bp                        */
  template <>
  void tree<double>::bpsimulate(int n, ctsdistribution &nexteventprior, rng &r)
  {
    double t=0.0;
    std::vector<node<double> *> lines(n);
    for (int i=0;i<n;i++) {
      lines[i]=&sample_[i];
      sample_[i].label=i;
      sample_[i].data()=1.0;
    }
    for (int nleft=n;nleft>1;nleft--) {
      int pos=n-nleft;
      double plustime=nexteventprior.sample(r);
      if (plustime<1E-8) plustime =1E-7;
      t += plustime;
  
      std::pair<int,int> chromo = r.sample2intsorted(0,nleft-1);
      ancestors_[pos].left=lines[chromo.first];
      ancestors_[pos].right=lines[chromo.second];
      ancestors_[pos].time=t;
      ancestors_[pos].data()=1.0;
      if (r()<0.5) 
        ancestors_[pos].label=lines[chromo.second]->label;
      else 
        ancestors_[pos].label=lines[chromo.first]->label;
      lines[chromo.second]=&ancestors_[pos];
      lines[chromo.first]=lines[nleft-1];
    }
    return;
  }
 template <>
  void node<std::vector<char> >::change_below(int pos,char value)
  {
    if (isleaf()) { // at a leaf
      data()[pos]=value;
    } else {
      left->change_below(pos,value);
      right->change_below(pos,value);
    }
    return;
  }

  /** Obtain a set of mutually exclusive haplogroups from the tree              
   * We shall ignore mutations that set them up for now - and assume that 
   * we have enough sites to determine the groups - otherwise we would have to 
   * do them from data rather than from the true tree                          

   Do this by recursively finding the split that separates the maximum to the left and right
  */

 template <>
 std::vector<int> tree<std::vector<char> >::haplogroups(int pos,double minsplit)
 {
   getnodefreq();

   int minsize = static_cast<int>(static_cast<double>(sample_.size())*minsplit);
      
   int hapgroup=0;
   std::vector<node<std::vector<char> > >::reverse_iterator i=ancestors_.rbegin();
   (*i).left->change_below(pos,hapgroup++);
   (*i).right->change_below(pos,hapgroup++);
   i++;
   while (i!=ancestors_.rend()) {
     if ((*i).countleft[0]>minsize)
       (*i).left->change_below(pos,hapgroup++);
     if ((*i).countright[0]>minsize)
       (*i).right->change_below(pos,hapgroup++); 
     i++;
   }

   std::vector<int> ret(sample_.size());

   std::map<char,int> hapg;
   int mx=0;
   for (size_t ii=0;ii<sample_.size();ii++) {
     if (hapg.find(sample_[ii].data()[pos])==hapg.end()) {
       hapg[sample_[ii].data()[pos]]=mx++;
     } 
   }

   for (size_t ii=0;ii<sample_.size();ii++) 
     ret[ii] =  hapg[sample_[ii].data()[pos]];

   return ret;
 }

}
