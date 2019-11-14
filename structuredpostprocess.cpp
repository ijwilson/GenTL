#include "structuredpostprocess.h"
#include "rtree.h"
#include <cassert>
#include <list>
#include <algorithm>
#include <vector>
#include <iterator>

void descendents(GenTL::rnode<std::vector<int> > *any, int mxpos)
{
  // allocate enough memory
  any->data().resize(mxpos);
  position mn=std::min<int>(any->right_->active.first(),any->left_->active.first());
  position mx=std::max<int>(any->right_->active.last(),any->left_->active.last());
 
  for (int i=0;i<mn;i++) any->data()[i]=0;
  for (int i=mn;i<=mx;i++) {
    if (any->right_->active.present(i)) {
      if  (any->left_->active.present(i))
	any->data()[i] = any->left_->coal->data()[i]|any->right_->coal->data()[i];
      else any->data()[i] = any->right(i)->data()[i];
    } else if  (any->left_->active.present(i)) {
      any->data()[i] 
	= any->left(i)->data()[i];
    } else  any->data()[i] =0;
  }
  for (int i=mx+1;i<mxpos;i++) 	any->data()[i]=0;
}

void PostprocessStructTree(GenTL::rtree<std::vector<int> > &t
			     ,const std::vector<int> &samp,int mxpos) {
    
  // Allocate the locations of the leaves
  // There should be far better ways of doing this !!
  for (int i=0;i<samp.size();i++) {
    t.sample[i].data().resize(mxpos);
    for (int j=0;j<mxpos;j++) t.sample[i].data()[j] = (int)(1<<samp[i]);
    std::copy(t.sample[i].data().begin(),t.sample[i].data().end(),std::ostream_iterator<int>(std::cout," "));
  }
  // now go through the ancestors and extract the data
 
  list<GenTL::rnode<std::vector<int> > *>::iterator aI=t.ancestors.begin();
  int count=0;
  while (aI!=t.ancestors.end()) {
    descendents(*aI,mxpos);
    aI++;
  }
} 
/** A debugging function -- print out the ancestors                             */
void   printancestors(std::ostream &o,rtree<std::vector<int> > &t)
{
    list<GenTL::rnode<std::vector<int> > *>::iterator aI=t.ancestors.begin();
    while (aI!=t.ancestors.end()) {
      o << (*aI)->time << " ";
      std::copy((*aI)->data().begin(),(*aI)->data().end(),std::ostream_iterator<int>(o," "));
      o << std::endl;
      aI++;
    }
}
