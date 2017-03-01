#ifndef STRUCTUREDPOSTPROCESS_H__
#define STRUCTUREDPOSTPROCESS_H__

#include <map>
#include <vector>
#include "gentldeclr.h"


/** note, a maximum of 16 populations                       */
void PostprocessStructTree(GenTL::rtree<std::vector<int> > &t
			     ,const std::vector<int> &samp, int maxpos);
			 
void     printancestors(std::ostream &o,GenTL::rtree<std::vector<int> > &t);

#endif
